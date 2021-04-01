MODULE python_wrapper

  USE mo_kind          , ONLY: dp
  USE mo_decomposition , ONLY: grid_init
  USE mo_cloud         , ONLY: sucloud
  USE mo_memory_g3b    , ONLY: construct_g3b, geosp
  USE mo_constants     , ONLY: amd, amo3
  USE mo_hyb           , ONLY: inihyb,vct_a,vct_b
  USE mo_radiation_parameters, ONLY: nb_sw, &
       ih2o, ico2, ich4, io3, io2, in2o, icfc,       &
       ighg, iaero, iza, co2vmr, ch4vmr, o2vmr,      &
       n2ovmr, cfcvmr, co2mmr, ch4mmr, o2mmr,        &
       n2ommr,  cemiss, diff
  USE mo_radiation     , ONLY: setup_radiation_new, rrtm_interface, pre_radiation

  IMPLICIT NONE

CONTAINS

  SUBROUTINE echam_radiation_offline(lat, lon, nlev, ntime, lolandh_in, loglach_in, xl_in, xi_in, &
       & aclc_in, p0_in, rhumidity_in, cdnc_in, t_surf_in, albedo_in, t_in, q_in, ao3_in, geosp_in, &
       & mu0_in, pf, ph, cdnc_cal, sups, supt, suptc, supsc, tdws, flt, fls, fltc, flsc)

    ! Dimensions of input data
    INTEGER, INTENT(in)   :: lat
    INTEGER, INTENT(in)   :: lon
    INTEGER, INTENT(in)   :: nlev
    INTEGER, INTENT(in)   :: ntime

    ! Logical input
    LOGICAL, INTENT(in)   :: cdnc_cal
    
    ! 2D input
    REAL(dp), INTENT(in)  :: lolandh_in(lon,lat)             ! Land-Sea fraction
    REAL(dp), INTENT(in)  :: loglach_in(lon,lat)             ! Glacier fraction

    ! 3D input
    REAL(dp), INTENT(in) :: p0_in(lon,lat,ntime)             ! pressure at the surface (Pa)
    REAL(dp), INTENT(in) :: t_surf_in(lon,lat,ntime)         ! temperature at the surface (K)
    REAL(dp), INTENT(in) :: albedo_in(lon,lat,ntime)         ! SW surface albedo
    REAL(dp), INTENT(in) :: geosp_in(lon,lat,ntime)          ! Surface geopotential (presently not used)
    REAL(dp), INTENT(in) :: mu0_in(lon,lat,ntime)               ! Cosine of solar zenith angle

    ! 4D input
    REAL(dp), INTENT(in) :: xl_in(lon,lat,nlev,ntime)        ! liquid water mixing ratio (kg/kg)
    REAL(dp), INTENT(in) :: xi_in(lon,lat,nlev,ntime)        ! ice mixing ratio (kg/kg)
    REAL(dp), INTENT(in) :: aclc_in(lon,lat,nlev,ntime)      ! cloud cover fraction
    REAL(dp), INTENT(in) :: rhumidity_in(lon,lat,nlev,ntime) ! relativ humidity
    REAL(dp), INTENT(in), OPTIONAL :: cdnc_in(lon,lat,nlev,ntime)      ! cloud cond. nuclei (m^-3)
    REAL(dp), INTENT(in) :: t_in(lon,lat,nlev,ntime)         ! temperature (K)
    REAL(dp), INTENT(in) :: q_in(lon,lat,nlev,ntime)         ! water vapor mixing ratio (kg/kg)
    REAL(dp), INTENT(in) :: ao3_in(lon,lat,nlev,ntime)       ! O3 mass mixing ratio (kg/kg)
    REAL(dp), INTENT(in) :: pf(lon,lat,nlev,ntime)           ! Pressure on full levels
    REAL(dp), INTENT(in) :: ph(lon,lat,nlev+1,ntime)         ! Pressure on half levels

    ! ----------------
    ! output variables from radiation. Set to intent(out) for output
    ! ----------------

    ! 3D output
    REAL(dp), INTENT(out) :: supt(lon,lat,ntime)             ! W/m2, surface upward LW flux
    REAL(dp), INTENT(out) :: sups(lon,lat,ntime)             ! W/m2, surface upward SW flux
    REAL(dp), INTENT(out) :: suptc(lon,lat,ntime)            ! W/m2, surface upward LW flux, clear sky
    REAL(dp), INTENT(out) :: supsc(lon,lat,ntime)            ! W/m2, surface upward SW flux, clear sky
    REAL(dp), INTENT(out) :: tdws(lon,lat,ntime)             ! W/m2, TOA SW irradiation
    REAL(dp)              :: semit(lon,lat,ntime)            ! surface emissivity

    ! 4D output
    REAL(dp), INTENT(out) :: flt(lon,lat,nlev+1,ntime)       ! W/m2, LW net flux
    REAL(dp), INTENT(out) :: fls(lon,lat,nlev+1,ntime)       ! W/m2, SW net flux
    REAL(dp), INTENT(out) :: fltc(lon,lat,nlev+1,ntime)      ! W/m2, LW net flux, clear sky
    REAL(dp), INTENT(out) :: flsc(lon,lat,nlev+1,ntime)      ! W/m2, SW net flux, clear sky
    REAL(dp) :: pswnir(lon,lat,ntime)                        ! net surface near infrared flux
    REAL(dp) :: pswdifnir(lon,lat,nlev,ntime)                ! fraction of diffuse near infrared
    REAL(dp) :: pswvis(lon,lat,ntime)                        ! net surface visible (250 - 680 nm)flux
    REAL(dp) :: pswdifvis(lon,lat,nlev,ntime)                ! fraction of diffuse visible
    REAL(dp) :: srfvisdn(lon,lat,nlev,ntime)                 ! surf. visible downward
    REAL(dp) :: srfpardn(lon,lat,ntime)                      ! surf. PAR downw.
    REAL(dp) :: frcpardif(lon,lat,nlev,ntime)                ! fraction of diffuse PAR


    ! ----------------
    ! Internally used variables
    ! ----------------

    ! More dimensions
    INTEGER, PARAMETER :: naer=5                           
    INTEGER, PARAMETER :: ncfc=2
    INTEGER, PARAMETER :: ntrac=1

    ! Loop variables
    INTEGER :: jk     ! level
    INTEGER :: j      ! lon
    INTEGER :: i      ! time
    INTEGER :: krow   ! lat

    ! Switches for greenhouse gases and aerosols. This needs some more work
    INTEGER :: ih2o_in
    INTEGER :: ico2_in
    INTEGER :: ich4_in
    INTEGER :: io3_in
    INTEGER :: io2_in
    INTEGER :: in2o_in
    INTEGER :: icfc_in
    INTEGER :: iaero_in

    ! Pi
    REAL(dp) :: pi

    ! 1D
    INTEGER  :: ktype(lon)                                   ! type of convection, not used
    REAL(dp) :: zprat, zn1, zn2
    
    ! 2D / Land/sea and glacier masks
    LOGICAL  :: loland(lon,lat)
    LOGICAL  :: loglac(lon,lat)
    REAL(dp) :: pgeom1(lon,nlev)                             ! geopotential above ground, not used

    ! 3D
    REAL(dp) :: pxtm1(lon,nlev,ntrac)                        ! tracer mass mixing ratios, not used
    
    ! 4D
    REAL(dp) :: qs(lon,lat,nlev,ntime)
    REAL(dp) :: ao3(lon,lat,nlev,ntime)
    REAL(dp) :: co2(lon,lat,nlev,ntime)
    REAL(dp) :: ch4(lon,lat,nlev,ntime)
    REAL(dp) :: o2(lon,lat,nlev,ntime)
    REAL(dp) :: n2o(lon,lat,nlev,ntime)
    REAL(dp) :: cfc(lon,lat,nlev,ncfc,ntime)
    REAL(dp) :: aer(lon,lat,nlev,naer,ntime)
    REAL(dp) :: aclcov(lon,lat,ntime)
    REAL(dp) :: aclcov1(lon,lat,ntime)
    REAL(dp) :: dpr(lon,lat,nlev, ntime)
    REAL(dp) :: th(lon,lat,nlev+1, ntime)
    REAL(dp) :: xl(lon,lat,nlev,ntime)        
    REAL(dp) :: xi(lon,lat,nlev,ntime)        
    REAL(dp) :: aclcac(lon,lat,nlev,ntime)    
    REAL(dp) :: mu0(lon,lat,ntime)            
    REAL(dp) :: rhumidity(lon,lat,nlev,ntime) 
    REAL(dp) :: cdnc(lon,lat,nlev,ntime)      
    REAL(dp) :: q(lon,lat,nlev,ntime)
    
    ! Initialize some variables
    ktype(:) = 0._dp
    pgeom1(:,:) = 0._dp
    pxtm1(:,:,:) = 0._dp
    pi = 4._dp * ATAN(1._dp)

    ! set negativ values to 0
    xl(:,:,:,:) = MAX(xl_in(:,:,:,:),0._dp)
    xi(:,:,:,:) = MAX(xi_in(:,:,:,:),0._dp)
    aclcac(:,:,:,:) = MAX(aclc_in(:,:,:,:),EPSILON(0._dp))
    mu0 = MAX(mu0_in(:,:,:),EPSILON(0._dp))
    rhumidity(:,:,:,:) = MAX(rhumidity_in(:,:,:,:),EPSILON(0._dp))
    q(:,:,:,:)=MAX(q_in(:,:,:,:),0._dp)

    ! Calculat sat. specific humidity
    qs(:,:,:,:) = q(:,:,:,:)/rhumidity(:,:,:,:)
    
    ! Initialize grid
    CALL grid_init(lon,lat,nlev,ntime)

    ! initialize vertical grid
    CALL inihyb

    ! allocate memory for geosp
    CALL construct_g3b
    geosp = geosp_in

    ! initialize some cloud parameters
    CALL sucloud

    ! initialize radiation and cloud optics
    CALL pre_radiation

    ih2o_in=1
    ico2_in=2
    ich4_in=2
    io3_in=1
    io2_in=2
    in2o_in=2
    iaero_in=0
    CALL setup_radiation_new(ih2o_in,ico2_in,ich4_in,io3_in,io2_in,in2o_in,iaero_in)
    iza = 1
    icfc_in=0


    ! Create land-sea and glacier mask
    DO krow=1,lat
       DO j=1,lon
          loland(j,krow) = .FALSE.
          loglac(j,krow) = .FALSE.
          IF ( lolandh_in(j,krow) > 0.5_dp ) THEN
             loland(j,krow) = .TRUE.
          END IF
          IF ( loglach_in(j,krow) > 0.5_dp ) THEN
             loglac(j,krow) = .TRUE.
          END IF
       END DO
    END DO

    IF (.NOT. cdnc_cal) THEN
       cdnc(:,:,:,:)=MAX(cdnc_in(:,:,:,:),EPSILON(0._dp))
    ELSE
       DO i=1,ntime
          DO krow=1,lat
             DO j=1,lon
                IF (loland(j,krow).AND.(.NOT.loglac(j,krow))) THEN
                   zn1 = 20._dp
                   zn2 = 180._dp
                ELSE
                   zn1 = 20._dp
                   zn2 = 80._dp
                ENDIF
                DO jk = 1, nlev
                   zprat=(MIN(8._dp,80000._dp/pf(j,krow,jk,i)))**2._dp
                   IF (pf(j,krow,jk,i) .LT. 80000._dp) THEN
                      cdnc(j,krow,jk,i)=1.E6_dp*(zn1+(zn2-zn1)*(EXP(1._dp-zprat)))
                   ELSE
                      cdnc(j,krow,jk,i)=zn2*1.e6_dp
                   END IF
                END DO
             END DO
          END DO
       END DO
    END IF

    ! Initialize greenhouse gases and aerosols
    co2(:,:,:,:) = co2mmr
    ch4(:,:,:,:) = ch4mmr
    o2(:,:,:,:) = o2mmr
    ao3(:,:,:,:) = ao3_in(:,:,:,:)*amo3/amd
    n2o(:,:,:,:) = n2ommr
    aer(:,:,:,:,:) = 0._dp

    SELECT CASE (icfc)
    CASE (0)
       cfc(:,:,:,:,:)=EPSILON(1._dp)
    CASE (2)
       cfc(:,:,:,1,:)=cfcvmr(1)
       cfc(:,:,:,2,:)=cfcvmr(2)
    CASE default
       STOP 'cfc: this "icfc" is not supported'
    END SELECT

    ! Calculate total cloud cover from cloud fraction per level
    DO j=1,lon
       aclcov1(j,1:lat,:) = 1._dp - aclcac(j,1:lat,1,:)
       DO jk = 2, nlev
          aclcov1(j,1:lat,:) = aclcov1(j,1:lat,:)                           &
               &     *(1._dp-MAX(aclcac(j,1:lat,jk,:),aclcac(j,1:lat,jk-1,:))) &
               &     /(1._dp-MIN(aclcac(j,1:lat,jk-1,:),1.-EPSILON(1._dp)))
       END DO
       aclcov(j,1:lat,:) = 1.-aclcov1(j,1:lat,:)
    END DO


    ! Thickness of layers
    dpr(:,:,1,:)=ph(:,:,2,:)
    DO jk=2,nlev
       dpr(:,:,jk,:)=ph(:,:,jk+1,:)-ph(:,:,jk,:)
    END DO


    ! Temperature at half levels
    ! --------------------------
    ! interpolate between full levels
    DO jk=2,nlev
       th(:,:,jk,:)= (t_in(:,:,jk-1,:)*pf(:,:,jk-1,:)*(pf(:,:,jk,:)-ph(:,:,jk,:)) &
            &     +t_in(:,:,jk,:)  *pf(:,:,jk,:)  *(ph(:,:,jk,:)-pf(:,:,jk-1,:)) ) &
            &   /(            ph(:,:,jk,:)  *(pf(:,:,jk,:)-pf(:,:,jk-1,:)) )
    END DO
    ! and extrapolate to TOA and surface
    th(:,:,1,:) = t_in(:,:,1,:)-pf(:,:,1,:)*(t_in(:,:,1,:)-th(:,:,2,:))/(pf(:,:,1,:)-ph(:,:,2,:))
    th(:,:,nlev+1,:) = t_surf_in(:,:,:)
    
    
    ! Run radiation code
    DO i=1,ntime
       DO krow=1,lat
          CALL rrtm_interface( &
               & iaero                     ,lon                       ,lon                   ,nlev                   ,& 
               & krow                      ,ntrac                     ,ktype                 ,nb_sw                  ,&
               & loland                    ,loglac                    ,cemiss                ,diff                   ,&
               & mu0(:,krow,i)             ,pgeom1                    ,albedo_in(:,krow,i)   ,albedo_in(:,krow,i)    ,&
               & albedo_in(:,krow,i)       ,albedo_in(:,krow,i)       ,albedo_in(:,krow,i)   ,albedo_in(:,krow,i)    ,&
               & pf(:,krow,:,i)            ,ph(:,krow,:,i)            ,p0_in(:,krow,i)       ,t_in(:,krow,:,i)       ,&
               & th(:,krow,:,i)            ,t_surf_in(:,krow,i)       ,q(:,krow,:,i)         ,qs(:,krow,:,i)         ,&
               & xl(:,krow,:,i)            ,xi(:,krow,:,i)            ,cdnc(:,krow,:,i)      ,aclcac(:,krow,:,i)     ,&
               & aclcov(:,krow,i)          ,ao3(:,krow,:,i)           ,co2(:,krow,:,i)       ,ch4(:,krow,:,i)        ,&
               & n2o(:,krow,:,i)           ,cfc(:,krow,:,:,i)         ,o2(:,krow,:,i)        ,pxtm1                  ,&
               & flt(:,krow,:,i)           ,fls(:,krow,:,i)           ,fltc(:,krow,:,i)      ,flsc(:,krow,:,i)       ,&
               & supt(:,krow,i)            ,sups(:,krow,i)            ,suptc(:,krow,i)       ,supsc(:,krow,i)        ,&
               & semit(:,krow,i)           ,tdws(:,krow,i)            ,pswnir(:,krow,i)      ,pswdifnir(:,krow,:,i)  ,&
               & pswvis(:,krow,i)          ,pswdifvis(:,krow,:,i)     ,srfvisdn(:,krow,:,i)  ,srfpardn(:,krow,i)     ,&
               & frcpardif(:,krow,:,i)                                                                                )
       END DO
    END DO


  END SUBROUTINE echam_radiation_offline

END MODULE python_wrapper
