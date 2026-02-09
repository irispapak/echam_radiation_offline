PROGRAM radiation
  USE mo_kind          , ONLY: sp,dp
  USE mo_constants     , ONLY: api,amd,amco2,amch4,amo3,amn2o,     &
       &                       g,cpd,vtmpc2
  USE mo_cloud         , ONLY: sucloud
  USE mo_tropopause    , ONLY: WMO_tropopause
  USE mo_newcld_optics,  ONLY: setup_newcld_optics, newcld_optics
  USE mo_hyb	       , ONLY: inihyb,vct_a,vct_b
  USE mo_decomposition,  ONLY: grid_init
  USE mo_memory_g3b,     ONLY: construct_g3b, geosp
  USE mo_radiation
  USE mo_radiation_parameters

  IMPLICIT NONE

  INTEGER, PARAMETER 	:: lat=96
  INTEGER, PARAMETER 	:: lon=192
  INTEGER, PARAMETER 	:: nlev=47
  INTEGER, PARAMETER 	:: naer=5
  INTEGER, PARAMETER 	:: nnewaer=0
  INTEGER, PARAMETER 	:: ncfc=2
  INTEGER, PARAMETER 	:: ntime=4
  INTEGER, PARAMETER    :: ntrac=1                  !LT
  INTEGER :: ktype(lon)=0                           !LT  < type of convection, not used
  REAL(dp) :: pgeom1(lon,nlev)=0._dp                !LT  < geopotential above ground, not used
  REAL(dp) :: pxtm1(lon,nlev,ntrac)=0._dp           !LT  < tracer mass mixing ratios, not used
  REAL(dp), PARAMETER 	:: pi=3.1415926		    !Pi

  character(len=20) :: filename
   
  character(len=160):: infilelsm=("T63GR15_lsm.nc")   !land-sea mask
  character(len=160):: infileglc=("T63GR15_glc.nc")   !glacier mask


  REAL(dp), PARAMETER :: nan=-999.0


  LOGICAL  :: loland(lon,lat)
  LOGICAL  :: loglac(lon,lat)
  REAL(dp)  :: lolandh(lon,lat)
  REAL(dp)  :: loglach(lon,lat)

  REAL(dp) :: mu0(lon,lat,ntime)    		! solar zenith angle
  REAL(dp) :: albedo(lon,lat,ntime)     	! SW albedo of surface
  REAL(dp) :: albedoo(lon,lat,1)

  REAL(dp) :: tsurf1(lon,lat,ntime)     	! surface temperature
  REAL(dp) :: tsurf2(lon,lat,ntime)     	! surface temperature

  REAL(dp) :: h0=7.            	  		! km, scale height
  REAL(dp) :: p0(lon,lat,ntime)
  REAL(dp) :: zf(lon,lat,nlev,ntime)    	! km, height
  REAL(dp) :: zh(lon,lat,nlev+1,ntime)  	! km, height at half levels
  REAL(dp) :: dz(lon,lat,nlev,ntime)    	! km, height thickness
  REAL(dp) :: zsrf(lon,lat,ntime)       	! Pa, height at the surface
  REAL(dp) :: pf(lon,lat,nlev,ntime)    	! Pa, pressure
  REAL(dp) :: ph(lon,lat,nlev+1,ntime)  	! Pa, pressure at half levels
  REAL(dp) :: dpr(lon,lat,nlev,ntime)   	! Pa, pressure thickness
  REAL(dp) :: psrf(lon,lat,ntime)       	! Pa, pressure at the surface

  REAL(dp) :: t(lon,lat,nlev,ntime)    		! K, temperature from file
  REAL(dp) :: tf(lon,lat,nlev,ntime)   		! K, temperature 
  REAL(dp) :: tf1(lon,lat,nlev,ntime)   	! K, temperature ! for LR feedback calculation
  REAL(dp) :: tf2(lon,lat,nlev,ntime)   	! K, temperature ! for LR feedback calculation

  REAL(dp) :: th(lon,lat,nlev+1,ntime)  	! K, temperature at half levels
  REAL(dp) :: th1(lon,lat,nlev+1,ntime) 	! K, temperature at half levels ! for LR feedback calculation
  REAL(dp) :: th2(lon,lat,nlev+1,ntime) 	! K, temperature at half levels ! for LR feedback calculation
  REAL(dp) :: tsrf(lon,lat,ntime)       	! K, temperature at the surface
  REAL(dp) :: dtsrf(lon,lat,ntime)      	! K, temperature at the surface ! for LR feedback calculation
  REAL(dp) :: tsrf2(lon,lat,ntime)      	! K, temperature at the surface ! for LR feedback calculation
  REAL(dp) :: tsrf1(lon,lat,ntime)      	! K, temperature at the surface ! for LR feedback calculation

  REAL(dp) :: q(lon,lat,nlev,ntime)     	! g/g, water vapor mixing ratio
  REAL(dp) :: qs(lon,lat,nlev,ntime)    	! g/g, water vapor saturation mixing ratio
  REAL(dp) :: xl(lon,lat,nlev,ntime)    	! g/g, liquid water mixing ratio
  REAL(dp) :: xi(lon,lat,nlev,ntime)    	! g/g, ice mixing ratio
  REAL(dp) :: rhumidity(lon,lat,nlev,ntime)	! relativ humidity

  REAL(dp) :: cdnc(lon,lat,nlev,ntime)  	! cloud cond. nuclei
  REAL(dp) :: icnc(lon,lat,nlev,ntime)  	! ice crystal number concentration
  REAL(dp) :: aclcac(lon,lat,nlev,ntime)	! m2/m2, cloud cover fraction
  REAL(dp) :: aclc(lon,lat,nlev,ntime)		! m2/m2, cloud cover fraction
  REAL(dp) :: aclcov(lon,lat,ntime)     	! m2/m2, total cloud cover from file
  REAL(dp) :: aclcov1(lon,lat,ntime)    	! m2/m2, total cloud cover internaly calcualted

  INTEGER  :: iaerh(lon,lat,nlev,ntime)
  REAL(dp) :: aer(lon,lat,nlev,naer,ntime)

  REAL(dp) :: ao3(lon,lat,nlev,ntime)   	! g/g, O3 mass mixing ratio
  REAL(dp) :: co2(lon,lat,nlev,ntime)   	! g/g, CO2 mass mixing ratio
  REAL(dp) :: ch4(lon,lat,nlev,ntime)   	! g/g, CH4 mass mixing ratio
  REAL(dp) :: o2(lon,lat,nlev,ntime)   	! g/g, O2 mass mixing ratio
  REAL(dp) :: n2o(lon,lat,nlev,ntime)   	! g/g, N2O mass mixing ratio
  REAL(dp) :: cfc(lon,lat,nlev,ncfc,ntime)	! m3/m3, CFC11 and CFC12 vol. mixing ratio


  REAL(dp), DIMENSION(lat,nlev,6) :: zswext ! extinction (sw)
  REAL(dp), DIMENSION(lat,nlev,6) :: zswssa ! single scattering albedo (sw)
  REAL(dp), DIMENSION(lat,nlev,6) :: zswasy ! asymmetry factor (sw)
  REAL(dp), DIMENSION(lat,nlev,16) :: zext  ! extinction (lw) - 16 RRTM Bands
  REAL(dp), DIMENSION(lon,lat,nlev,6) :: zswext1 ! extinction (sw)
  REAL(dp), DIMENSION(lon,lat,nlev,6) :: zswssa1 ! single scattering albedo (sw)
  REAL(dp), DIMENSION(lon,lat,nlev,6) :: zswasy1 ! asymmetry factor (sw)
  REAL(dp), DIMENSION(lon,lat,nlev,16) :: zext1  ! extinction (lw) - 16 RRTM Bands

  ! output variables
  ! ----------------
  REAL(dp) :: flt(lon,lat,nlev+1,ntime) 	! W/m2, LW net flux
  REAL(dp) :: fls(lon,lat,nlev+1,ntime) 	! W/m2, SW net flux
  REAL(dp) :: fltc(lon,lat,nlev+1,ntime)	! W/m2, LW net flux, clear sky
  REAL(dp) :: flsc(lon,lat,nlev+1,ntime)	! W/m2, SW net flux, clear sky

  REAL(dp) :: supt(lon,lat,ntime)       	! W/m2, surface upward LW flux
  REAL(dp) :: sups(lon,lat,ntime)       	! W/m2, surface upward SW flux
  REAL(dp) :: suptc(lon,lat,ntime)      	! W/m2, surface upward LW flux, clear sky
  REAL(dp) :: supsc(lon,lat,ntime)      	! W/m2, surface upward SW flux, clear sky

  REAL(dp) :: semit(lon,lat,ntime)      	! surface emissivity
  REAL(dp) :: tdws(lon,lat,ntime)       	! W/m2, TOA SW irradiation
  REAL(dp) :: msemit(lon,lat,ntime)     	! surface emissivity
  REAL(dp) :: mtdws(lon,lat,ntime)      	! W/m2, TOA SW irradiation

  REAL(dp) :: dflt(lon,lat,nlev,ntime)  	! W/m2, divergence of LW net flux
  REAL(dp) :: dfls(lon,lat,nlev,ntime)  	! W/m2, divergence of SW net flux
  REAL(dp) :: dm(lon,lat,nlev,ntime)    	! kg/m2, density in layer
  REAL(dp) :: cp(lon,lat,nlev,ntime)    	! W/K/kg, cp(q)
  REAL(dp) :: qradt(lon,lat,nlev,ntime) 	! K/s, LW heating rate
  REAL(dp) :: qrads(lon,lat,nlev,ntime) 	! K/s, SW heating rate

  REAL(dp) :: xlon(lon), xlat(lat)		! Longitude/latitude arrays
  REAL(dp) :: xtime(ntime)			! time axis
  INTEGER  :: ilev(nlev), i, j, k, jk

  ! shortwave radiation : visible, NIR and fraction of diffuse visible, NIR
  REAL(dp) :: pswnir(lon,lat,ntime)            	! net surface near infrared flux
  REAL(dp) :: pswdifnir(lon,lat,nlev,ntime)    	! fraction of diffuse near infrared
  REAL(dp) :: pswvis(lon,lat,ntime)            	! net surface visible (250 - 680 nm)flux
  REAL(dp) :: pswdifvis(lon,lat,nlev,ntime)    	! fraction of diffuse visible
!LT
  REAL(dp) :: srfvisdn(lon,lat,nlev,ntime)    	! surf. visible downward
  REAL(dp) :: srfpardn(lon,lat,ntime)            	! surf. PAR downw.
  REAL(dp) :: frcpardif(lon,lat,nlev,ntime)    	! fraction of diffuse PAR

  REAL(dp) :: EqT, Tsol		   		! Equation of time, solar time
  REAL(dp) :: hang, dec, lati, loni  		! hour angle, declination, actual latitude, longitude
  REAL(dp) :: zang(lon,lat,ntime)	     	! solar zenith angle
  REAL(dp) :: tloc, tloc_c	             	! local time, local time correction for longitude

  REAL(dp) :: tropo(lon,lat,ntime)

  INTEGER                    :: ierr, myid      ! for message passing
  INTEGER                    :: krow
  REAL(dp) :: ralphr(lon,lat,nlev), delpr(lon,lat,nlev),          &
                             rdelpr(lon,lat,nlev),  zrlnpr(lon,lat,nlev)

  ! initialize grid

  CALL grid_init(lon,lat,nlev,ntime)

  ! allocate memory

  CALL construct_g3b

  ! initialize vertical grid
  
  CALL inihyb

  ! initialize some cloud parameters
 
  CALL sucloud
 
  
  ! initialize radiation and cloud optics
 
  CALL pre_radiation
  CALL setup_radiation


  !some namelist input
 

  in: SELECT CASE (input)
  CASE (0)
     infile=infile
  CASE (1)
     infile=infile
  CASE default
     STOP 'not supported'
  END SELECT in

  dy: SELECT CASE (day)
  CASE (0)
     day=1
  CASE (1)
     nday=nday
  CASE default
     STOP 'not supported'
  END SELECT dy


  !read land-sea and glacier mask

  CALL read_nc3 (infilelsm, lon, lat, 1, xlon, xlat, xtime,  'SLM',lolandh)
  CALL read_nc3 (infileglc, lon, lat, 1, xlon, xlat, xtime,  'GLAC',loglach)
  do jk=1,lat
      do j=1,lon
       loland(j,jk) = .FALSE.
       loglac(j,jk) = .FALSE.
       IF ( lolandh(j,jk) > 0.5 ) THEN
       loland(j,jk) = .TRUE.
       ENDIF
       IF ( loglach(j,jk) > 0.5 ) THEN
       loglac(j,jk) = .TRUE.
       ENDIF       
      enddo
  enddo


  ! ------------------------------------
  ! read cloud water/ice from file

  CALL read_nc (infile, lon, lat, nlev, ntime, xlon, xlat, ilev, xtime, 'xl',xl)
  CALL read_nc (infile, lon, lat, nlev, ntime, xlon, xlat, ilev, xtime, 'xi',xi)
  xl(:,:,:,:)=max(xl(:,:,:,:),0._dp) ! set negativ values to 0
  xi(:,:,:,:)=max(xi(:,:,:,:),0._dp) ! set negativ values to 0

  ! read cloud cover fraction

  CALL read_nc (infile, lon, lat, nlev, ntime, xlon, xlat, ilev, xtime, 'aclc',aclc)
  aclcac(:,:,:,:)=max(aclc(:,:,:,:),EPSILON(0._dp))

  !read surface pressure

  CALL read_nc3 (infile, lon, lat, ntime, xlon, xlat, xtime,  'aps',p0)

  ! read rel. humidity to calculate saturation specific humidity

  CALL read_nc (infile, lon, lat, nlev, ntime, xlon, xlat, ilev, xtime, 'rhumidity',rhumidity)
  rhumidity(:,:,:,:)=max(rhumidity(:,:,:,:),EPSILON(0._dp)) ! set negative values to EPSILON(1.)

  ! read cloud droplet number concentration

  CALL read_nc (infile, lon, lat, nlev, ntime, xlon, xlat, ilev, xtime, 'acdnc',cdnc)
  cdnc(:,:,:,:)=max(cdnc(:,:,:,:),EPSILON(0._dp)) ! set negative values to EPSILON(1.)

  ! read ice crystal number concentration

  CALL read_nc (infile, lon, lat, nlev, ntime, xlon, xlat, ilev, xtime, 'aicnc',icnc)
  icnc(:,:,:,:)=max(icnc(:,:,:,:),EPSILON(0._dp)) ! set negative values to EPSILON(1.)


  ! read albedo

  CALL read_nc3 (infile, lon, lat, 1, xlon, xlat, xtime,  'albedo',albedoo)
  DO i=1,ntime
   albedo(:,:,i)=albedoo(:,:,1)
  ENDDO
  
  ! read surface temperature

  CALL read_nc3 (infile, lon, lat, ntime, xlon, xlat, xtime,  'tsurf',tsurf1)
  CALL read_nc3 (infile, lon, lat, ntime, xlon, xlat, xtime,  'tsurf',tsurf2)

  !
  ! vertical coordinate system of input files
  ! -----------------------------------------
!
!  ! in ECHAM the pressure at the full levels is defined in the middle
!  ! between the upper and lower layer interfaces: 
!  ! pf(jk)=(ph(jk)+ph(jk+1))/2 with ph(1)=0. and ph(nlev+1)=psrf
! Pressure Levels are now constructed from vertical coordinat table A/B
! and surface pressure, like in ECHAM

  do i=1,ntime
    do jk=1,lat
      do j=1,lon
       ph (j,jk,:,i) = vct_a(:)+ vct_b(:)* p0(j,jk,i)
      enddo
    enddo
  enddo
  DO jk=1,nlev
     pf(:,:,jk,:)=0.5_dp*(ph(:,:,jk,:)+ph(:,:,jk+1,:))
  END DO

  ! calculate total cloud cover from cloud fraction per level

  do i=1,lon
     aclcov1(i,1:lat,:) = 1. - aclcac(i,1:lat,1,:)
     DO jk = 2, nlev
        aclcov1(i,1:lat,:) = aclcov1(i,1:lat,:)                           &
             &     *(1.-MAX(aclcac(i,1:lat,jk,:),aclcac(i,1:lat,jk-1,:))) &
             &     /(1.-MIN(aclcac(i,1:lat,jk-1,:),1.-EPSILON(1._dp)))
     END DO
     aclcov(i,1:lat,:) = 1.-aclcov1(i,1:lat,:)
  end do

  zenith: SELECT CASE (iza)
  CASE (1)
  ! calculate the solar zenith angle. Following:
  ! http://solardat.uoregon.edu/SolarRadiationBasics.html
  ! daniel.klocke@zmaw.de
  ! calcluate equation of time for solar time difference
  ! calculate difference of solar time to local time (gmt) other wise longitut correction is needed
  Tloc=-6._dp				! Tloc is set to -6, to start on 0 in the first loop
      print *,'calculating zenith angles for Julian day:',nday
  DO j=1,ntime
    Tloc=Tloc+6._dp
    IF (Tloc > 23. ) THEN
      Tloc=0._dp
      nday=nday+1.
      print *,'calculating zenith angles for Julian day:',nday
    ENDIF
      print *,Tloc,'oclock'
    IF (nday > 0. .and. nday .le. 106.) then
      EqT = -14.2 * sin( pi * (nday + 7.)/111.)
    ELSEIF (nday .gt. 106. .and. nday .le. 166.) then
      EqT = 4. * sin( pi * (nday - 106.)/59.)
    ELSEIF (nday .gt. 166. .and. nday .le. 246.) then
      EqT = -6.5 * sin( pi * (nday - 166.)/80.)
    ELSEIF (nday .gt. 246. .and. nday .le. 366.) then
      EqT = 16.4 * sin( pi * (nday - 247.)/113.)
    ENDIF
			
    DO k=1,lon     				! zenith angle calculation for longitude
      loni = (+180. - 360./lon * k )		! longitude
      Tloc_c = Tloc - (loni - 180)/15.		! time correction for longitude
      Tsol = Tloc_c + EqT/60.			! solar time
      ! hour angle
      hang = pi * (12.-Tsol)/12.
      ! declination
      dec = (23.45*pi/180.) * (cos(2.*pi/365. * (nday-172.)))
      ! solar zenith angle 
      DO i=1,lat
        lati = (90. - 180./lat * i )* (pi/180.)
        zang(k,i,j) =  acos(sin(lati)*sin(dec)+(cos(lati)*cos(dec)*cos(hang)))*180./pi
      ENDDO
    ENDDO
  ENDDO
  mu0(:,:,:)=zang(:,:,:)
  mu0(:,:,:)=min(mu0(:,:,:),90._dp)
  CASE (2)
     mu0(:,:,:)=za0
  CASE default
     STOP 'zenith: this "iza" is not supported'
  END SELECT zenith
  mu0(:,:,:)=COS(mu0(:,:,:)/180.*api)

  ! profiles of temperature and gases 
  ! at full levels, depending on switches 
  ! ------------------------------------------
!LT
!  temp: SELECT CASE (itemp)
!  CASE (1)
     ! read temperature
     CALL read_nc(infile, lon, lat, nlev, ntime, xlon, xlat, ilev, xtime, 't',t)
     tf2=t
     CALL read_nc(infile, lon, lat, nlev, ntime, xlon, xlat, ilev, xtime, 't',t)
     tf1=t
!  CASE (2)
!     tf(:,:,:,:)=t0
!  CASE default
!     STOP 'temp: this "itemp" is not supported'
!  END SELECT temp

  h2o: SELECT CASE (ih2o)
  CASE (0)
     q(:,:,:,:)=EPSILON(1._dp)
     xl(:,:,:,:)=0._dp
     xi(:,:,:,:)=0._dp
  CASE (1)
     ! read specific humidity
     CALL read_nc(infile, lon, lat, nlev, ntime, xlon, xlat, ilev, xtime, 'q',q)
  CASE default
     STOP 'h2o: this "ih2o" is not supported'
  END SELECT h2o
  q(:,:,:,:)=max(q(:,:,:,:),0._dp) ! set negativ values to 0

  ! calculate saturation specific humidity
  qs(:,:,:,:) = q(:,:,:,:)/rhumidity(:,:,:,:)


  cd: SELECT CASE (ico2)
  CASE (0)
     co2(:,:,:,:)=EPSILON(1._dp)
  CASE (2)
     co2(:,:,:,:)=co2mmr
  CASE default
     STOP 'co2: this "ico2" is not supported'
  END SELECT cd


  methane: SELECT CASE (ich4)
  CASE (0)
     ch4(:,:,:,:)=EPSILON(1._dp)
  CASE (1)
     ! read ch4 in volume mixing ratio
!     CALL read_nc(infile, lon, lat, nlev, ntime, xlon, xlat, ilev, xtime, 'CH4',ch4)
     ! convert to mass mixing ratio
!     ch4=ch4*amch4/amd
  CASE (2)
     ch4(:,:,:,:)=ch4mmr
  CASE default
     STOP 'ch4: this "ich4" is not supported'
  END SELECT methane

!LT
  oxygen: SELECT CASE (io2)
  CASE (0)
     o2(:,:,:,:)=EPSILON(1._dp)
  CASE (1)
     ! read o2 in volume mixing ratio
!     CALL read_nc(infilex, lon, lat, nlev, ntime, xlon, xlat, ilev, xtime, 'O2',o2)
     ! convert to mass mixing ratio
!     o2=o2*amo2/amd
  CASE (2)
     o2(:,:,:,:)=o2mmr
  CASE default
     STOP 'o2: this "io2" is not supported'
  END SELECT oxygen

  ozone: SELECT CASE (io3)
  CASE (0)
     ao3(:,:,:,:)=EPSILON(1._dp)
  CASE (1)
     ! read o3 volume mixing ratio in ppmv
     CALL read_nc(infile, lon, lat, nlev, ntime, xlon, xlat, ilev, xtime, 'ao3',ao3)
     ! convert from ppmv to volume mixing ratio
     ao3(:,:,:,:)=ao3(:,:,:,:)
     ! convert to mass mixing ratio
     ao3(:,:,:,:)=ao3(:,:,:,:)*amo3/amd
  CASE default
     STOP 'o3: this "io3" is not supported'
  END SELECT ozone


  cn2o: SELECT CASE (in2o)
  CASE (0)
     n2o(:,:,:,:)=EPSILON(1._dp)
  CASE (1)
     CALL read_nc(infilex, lon, lat, nlev, ntime, xlon, xlat, ilev, xtime, 'N2O',n2o)
  CASE (2)
     n2o(:,:,:,:)=n2ommr
  CASE default
     STOP 'n2o: this "in2o" is not supported'
  END SELECT cn2o


  cfcs: SELECT CASE (icfc)
  CASE (0)
     cfc(:,:,:,:,:)=EPSILON(1._dp)
  CASE (2)
     cfc(:,:,:,1,:)=cfcvmr(1)
     cfc(:,:,:,2,:)=cfcvmr(2)
  CASE default
     STOP 'cfc: this "icfc" is not supported'
  END SELECT cfcs


  aerosol: SELECT CASE (iaero)
  CASE (0)
     aer(:,:,:,:,:)=0._dp
  CASE (3)
     aer(:,:,:,:,:)=0._dp
  CASE default
     STOP 'aerosols: this "iaero" is not supported'
  END SELECT aerosol
  ! aerosol values are at least epsilon
  aer(:,:,:,:,:)=MAX(aer(:,:,:,:,:),EPSILON(1._dp))

  !define surface pressure
  do i=1,ntime
  psrf(:,:,i)=p0(:,:,i)
  enddo


  ! thickness of layers
  ! -------------------

  dpr(:,:,1,:)=ph(:,:,2,:)
  DO jk=2,nlev
     dpr(:,:,jk,:)=ph(:,:,jk+1,:)-ph(:,:,jk,:)
  END DO

  ! temperature at half levels
  ! --------------------------
  ! interpolate between full levels

  DO jk=2,nlev
     th1(:,:,jk,:)= ( tf1(:,:,jk-1,:)*pf(:,:,jk-1,:)*(pf(:,:,jk,:)-ph(:,:,jk,:)  ) &
          &     +tf1(:,:,jk,:)  *pf(:,:,jk,:)  *(ph(:,:,jk,:)-pf(:,:,jk-1,:)) ) &
          &   /(            ph(:,:,jk,:)  *(pf(:,:,jk,:)-pf(:,:,jk-1,:)) )
  END DO

  ! and extrapolate to TOA and surface


   th1(:,:,1,:)     =tf1(:,:,1,:)-pf(:,:,1,:)*(tf1(:,:,1,:)-th1(:,:,2,:))/(pf(:,:,1,:)-ph(:,:,2,:))



  th1(:,:,nlev+1,:)=tsurf1(:,:,:)

  !surface
  
  tsrf1(:,:,:)=tsurf1(:,:,:)


 !!! calcualte perturbed temperature profiles for lapse rate feedback

  DO jk=2,nlev
     th2(:,:,jk,:)= ( tf2(:,:,jk-1,:)*pf(:,:,jk-1,:)*(pf(:,:,jk,:)-ph(:,:,jk,:)  ) &
          &     +tf2(:,:,jk,:)  *pf(:,:,jk,:)  *(ph(:,:,jk,:)-pf(:,:,jk-1,:)) ) &
          &   /(            ph(:,:,jk,:)  *(pf(:,:,jk,:)-pf(:,:,jk-1,:)) )
  END DO

 !!! extrapolate perturbed temperature profile to TOA and surface

  th2(:,:,1,:)     =tf2(:,:,1,:)-pf(:,:,1,:)*(tf2(:,:,1,:)-th2(:,:,2,:))/(pf(:,:,1,:)-ph(:,:,2,:))


  th2(:,:,nlev+1,:)=tsurf2(:,:,:)

  ! surface
 
  tsrf2(:,:,:)=tsurf2(:,:,:)



  ! read tropopause height

  CALL read_nc3 (infile, lon, lat, ntime, xlon, xlat, xtime, 'tropo',tropo)


  !read geopotential height

  CALL read_nc3 (infile, lon, lat, ntime, xlon, xlat, xtime,  'geosp',geosp)



! case of NO lapse rate feedback calculation:

   DO jk=1,nlev
  	tf(:,:,jk,:)=tf1(:,:,jk,:)
   END DO
   DO jk=1,nlev+1
 	th(:,:,jk,:)=th1(:,:,jk,:)
   END DO

   tsrf(:,:,:)=tsurf1(:,:,:)

!end case of NO lapse rate feedback calculation



! lapse rate feedback
!!!!!!! set temperature below tropopause to peturbed temperature profile
!
!  do jk=1,nlev
!  	where (pf(:,:,jk,:) > tropo(:,:,:))
!		tf(:,:,jk,:)=tsrf1(:,:,:)+(tf2(:,:,jk,:)-tsrf2(:,:,:))
!	else where
!  		tf(:,:,jk,:)=tf1(:,:,jk,:)
!	end where
!  end do
!
!  do jk=1,nlev+1
!	where(ph(:,:,jk,:) > tropo(:,:,:))
!		th(:,:,jk,:)=tsrf1+(th2(:,:,jk,:)-tsrf2(:,:,:))
!	else where
!		th(:,:,jk,:)=th1(:,:,jk,:)
!	end where
!  end do
!
!
!  tsrf(:,:,:)=tsurf1(:,:,:)
!
!
!
!
!!!!!!!!!!!!!! end extra calculations for lapse rate feedback





print *,'-------------------------------------------'
print *,'------------setting done!------------------'
print *,'-------------------------------------------'


  DO i=1,ntime


  print *,'calculating radiation for timestep: ',i,' from: ',ntime

     DO krow=1,lat
!
!LT
     CALL rrtm_interface( &
         & iaero                     ,lon                       ,lon                   ,nlev                   ,& 
         & krow                      ,ntrac                     ,ktype                 ,nb_sw                  ,&
         & loland                    ,loglac                    ,cemiss                ,diff                   ,&
         & mu0(:,krow,i)             ,pgeom1                    ,albedo(:,krow,i)      ,albedo(:,krow,i)       ,&
         & albedo(:,krow,i)          ,albedo(:,krow,i)          ,albedo(:,krow,i)      ,albedo(:,krow,i)       ,&
         & pf(:,krow,:,i)            ,ph(:,krow,:,i)            ,psrf(:,krow,i)        ,tf(:,krow,:,i)         ,&
         & th(:,krow,:,i)            ,tsrf(:,krow,i)            ,q(:,krow,:,i)         ,qs(:,krow,:,i)         ,&
         & xl(:,krow,:,i)            ,xi(:,krow,:,i)            ,cdnc(:,krow,:,i)      ,icnc(:,krow,:,i)       ,&
         & aclcac(:,krow,:,i)     ,&
         & aclcov(:,krow,i)          ,ao3(:,krow,:,i)           ,co2(:,krow,:,i)       ,ch4(:,krow,:,i)        ,&
         & n2o(:,krow,:,i)           ,cfc(:,krow,:,:,i)         ,o2(:,krow,:,i)        ,pxtm1                  ,&
         & flt(:,krow,:,i)           ,fls(:,krow,:,i)           ,fltc(:,krow,:,i)      ,flsc(:,krow,:,i)       ,&
         & supt(:,krow,i)            ,sups(:,krow,i)            ,suptc(:,krow,i)       ,supsc(:,krow,i)        ,&
         & semit(:,krow,i)           ,tdws(:,krow,i)            ,pswnir(:,krow,i)      ,pswdifnir(:,krow,:,i)  ,&
         & pswvis(:,krow,i)          ,pswdifvis(:,krow,:,i)     ,srfvisdn(:,krow,:,i)  ,srfpardn(:,krow,i)     ,&
         & frcpardif(:,krow,:,i)                                                                                )

	
    ENDDO
  ENDDO		! loop time


  !
  ! write profiles
  116 format("flth_",F0.0,"nc")
  Write (filename,116)nday
  CALL write_nc(filename,lon, lat, nlev, ntime, xlon, xlat, ilev, xtime,'flth',flt(:,:,1:nlev,:))
  117 format("flthc_",F0.0,"nc")
  Write (filename,117)nday
  CALL write_nc(filename,lon, lat, nlev, ntime, xlon, xlat, ilev, xtime,'flthc',fltc(:,:,1:nlev,:))
  118 format("flsh_",F0.0,"nc")
  Write (filename,118)nday
  CALL write_nc(filename,lon, lat, nlev, ntime, xlon, xlat, ilev, xtime,'flsh',fls(:,:,1:nlev,:))
  119 format("flshc_",F0.0,"nc")
  Write (filename,119)nday
  CALL write_nc(filename,lon,lat,nlev,ntime,xlon,xlat,ilev,xtime,'flshc',flsc(:,:,1:nlev,:))
  !
  ! write surface fields
!  120 format("supt_",F0.0,"nc")
!  Write (filename,120)nday
!  CALL write_nc3(filename,lon, lat, ntime, xlon, xlat, xtime, 'supt',supt)
!  121 format("sups_",F0.0,"nc")
!  Write (filename,121)nday
!  CALL write_nc3(filename,lon, lat, ntime, xlon, xlat, xtime, 'sups',sups)
!  122 format("semit_",F0.0,"nc")
!  Write (filename,122)nday
!  CALL write_nc3(filename,lon, lat, ntime, xlon, xlat, xtime, 'semit',semit)
!  123 format("tdws_",F0.0,"nc")
!  Write (filename,123)nday
!  CALL write_nc3(filename,lon, lat, ntime, xlon, xlat, xtime, 'tdws',tdws)
!  124 format("suptc_",F0.0,"nc")
!  Write (filename,124)nday
!  CALL write_nc3(filename,lon, lat, ntime, xlon, xlat, xtime, 'suptc',suptc)
!  125 format("supsc_",F0.0,"nc")
!  Write (filename,125)nday
!  CALL write_nc3(filename,lon, lat, ntime, xlon, xlat, xtime, 'supsc',supsc)


  ! compute radiative heating rates
  ! -------------------------------
!  DO jk=1,nlev
!     dflt(:,:,jk,:) =-(flt(:,:,jk+1,:)-flt(:,:,jk,:))
!     dfls(:,:,jk,:) =-(fls(:,:,jk+1,:)-fls(:,:,jk,:))
!     dm(:,:,jk,:)   = dpr(:,:,jk,:)/g
!     cp(:,:,jk,:)   = cpd/(1.+vtmpc2*q(:,:,jk,:))
!     qradt(:,:,jk,:)=dflt(:,:,jk,:)/dm(:,:,jk,:)/cp(:,:,jk,:)
!     qrads(:,:,jk,:)=dfls(:,:,jk,:)/dm(:,:,jk,:)/cp(:,:,jk,:)
!  END DO
  !
  ! write profiles
!  126 format("dflt_",F0.0,"nc")
!  Write (filename,126)nday
!  CALL write_nc(filename,lon, lat, nlev, ntime, xlon, xlat, ilev, xtime, 'dfltf',dflt(:,:,1:nlev,:))
!  127 format("dfls_",F0.0,"nc")
!  Write (filename,127)nday
!  CALL write_nc(filename,lon, lat, nlev, ntime, xlon, xlat, ilev, xtime, 'dflsf',dfls(:,:,1:nlev,:))
!  128 format("dm_",F0.0,"nc")
!  Write (filename,128)nday
!  CALL write_nc(filename,lon, lat, nlev, ntime, xlon, xlat, ilev, xtime, 'dmf',dm)
!  129 format("cp_",F0.0,"nc")
!  Write (filename,129)nday
!  CALL write_nc(filename,lon, lat, nlev, ntime, xlon, xlat, ilev, xtime, 'cpf',cp)
!  130 format("ql_",F0.0,"nc")
!  Write (filename,130)nday
!  CALL write_nc(filename,lon, lat, nlev, ntime, xlon, xlat, ilev, xtime, 'qlf',qradt*86400)
!  131 format("qs_",F0.0,"nc")
!  Write (filename,131)nday
!  CALL write_nc(filename,lon, lat, nlev, ntime, xlon, xlat, ilev, xtime, 'qsf',qrads*86400)
 
contains

!Deallocate variable 4 dimension
subroutine dealloc_var4d(varname, var)

  implicit none

  real(dp), allocatable,  intent(inout), dimension(:,:,:) ::  var
  
  character (len = *), intent(in) :: varname

  integer :: dealloc_status
 
  deallocate(var, stat=dealloc_status)

  if(dealloc_status .ne. 0) then
      write(*,'("Error while deallocating memory: ",A)') varname
      stop
  endif
  
  return
 
end subroutine dealloc_var4d


subroutine read_nc (infile, nlon, nlat, nmlev, ntime, xlon, xlat, ilev, xtime,  varname, var)
  

   implicit none

   include 'netcdf.inc'

   integer                              :: nlon, nlat, nmlev, ntime, ntime2
   real(dp), dimension(nlon)        :: xlon
   real(dp), dimension(nlat)        :: xlat
   real(dp), dimension(ntime)       :: xtime
   integer, dimension(nmlev)            :: ilev
   real(dp), dimension(nlon,nlat,nmlev,ntime) :: var
   character(len=*), intent(in) 	:: infile, varname

   character(len=256) 			:: nctitle, ncechamver, ncdate 
   character(len=12) 			:: wdate, wtime
   character(len=25) 			:: timestr
   
   integer 				:: vtime(8)

   
   real(dp), dimension(ntime) ::  time

   integer, dimension(4) :: start, count

   integer :: ncstatus, alloc_status, dealloc_status, ncomode
   integer :: ncid, latdimid, londimid, timedimid, timestep, mlevdimid
   integer :: lonid, latid, timeid, mlevid, varid
   integer :: dim1, dim2, dim3

   !Open file
   call nc_check(nf_open(infile, nf_nowrite, ncid))
     
   !Read grid
   call nc_check(nf_inq_dimid(ncid, "time", timedimid))
   !call nc_check(nf_get_dimlen(ncid, timedimid, ntime2))
   

  ! write(*,*) ncnlon, ncnlat, ncnmlev, ncnilev, ncntime


   !Read longitude
   call nc_check(nf_inq_varid(ncid, "lon", lonid))
   call nc_check(nf_get_vara_double(ncid, lonid, 1, nlon, xlon))

   !Read latitude
   call nc_check(nf_inq_varid(ncid, "lat", latid))
   call nc_check(nf_get_vara_double(ncid, latid, 1, nlat, xlat))

   !Read level
   call nc_check(nf_inq_varid(ncid, "lev", mlevid))
   call nc_check(nf_get_vara_int(ncid, mlevid, 1, nmlev, ilev))

   !Read time
   call nc_check(nf_inq_varid(ncid, "time", timeid))
   call nc_check(nf_get_vara_double(ncid, timeid, 1, ntime, xtime))

   !Read global attributes
!   call nc_check(nf_get_att_text(ncid, nf_global, "title", nctitle))
!   call nc_check(nf_get_att_text(ncid, nf_global, "echam_version", ncechamver))
!   call nc_check(nf_get_att_text(ncid, nf_global, "date_time", ncdate))


   call nc_check(nf_inq_varid(ncid, varname, varid))
   count(1) = nlon
   count(2) = nlat
   count(3) = nmlev
   count(4) = ntime
   start(1) = 1
   start(2) = 1
   start(3) = 1
   start(4) = 1

!   do timestep = 1, ntime
!      start(4) = timestep
      call nc_check(nf_get_vara_double(ncid, varid, start, count, var))
!   end do
  
   !Close file
   call nc_check(nf_close(ncid))

end subroutine read_nc

!-----------------------------------------------------------------------
!Modules
!-----------------------------------------------------------------------

!Netcdf error check
subroutine nc_check(status)

  implicit none

  include 'netcdf.inc'

  integer, intent(in) :: status

  if(status .ne. nf_noerr) then
      write(*,'(A)') nf_strerror(status)
      stop
   endif

  return

end subroutine nc_check

!-----------------------------------------------------------------------
!Error check
subroutine err_check(status, err_name)

  implicit none

  integer, intent(in) :: status
  character (len = *), intent(in) :: err_name

  if(status .ne. 0) then
      write(*,'("Error: ",A)') err_name
      stop
   endif

  return

end subroutine err_check

!-----------------------------------------------------------------------

!Read dimension
subroutine nc_getdim(fileid, dimname, dimvar)

  implicit none

  include 'netcdf.inc'
  
  integer ,intent(in) :: fileid, dimvar
  character (len = *), intent(in) :: dimname

  integer :: ncstatus, dimid
 
          
  !Fetch ID for dim
  call nc_check(nf_inq_dimid(fileid, dimname, dimid))

  !Read dim
  call nc_check(nf_inq_dimlen(fileid, dimid, dimvar))
      
  return
 
end subroutine nc_getdim


subroutine read_nc3 (infile, nlon, nlat, ntime, xlon, xlat, xtime,  varname, var3)


   implicit none

   include 'netcdf.inc'

   integer                              :: nlon, nlat, ntime
   real(dp), dimension(nlon)        :: xlon
   real(dp), dimension(nlat)        :: xlat
   real(dp), dimension(ntime)       :: xtime
   real(dp), dimension(nlon,nlat,ntime) :: var3
   character(len=*), intent(in) 	:: infile, varname

   character(len=256) 			:: nctitle, ncechamver, ncdate 
   character(len=12) 			:: wdate, wtime
   character(len=25) 			:: timestr
   
   integer 				:: vtime(8)

   real(dp), dimension(ntime) ::  time

   integer, dimension(3) :: start, count

   integer :: ncstatus, alloc_status, dealloc_status, ncomode
   integer :: ncid, latdimid, londimid, timedimid, timestep
   integer :: lonid, latid, timeid, varid
   integer :: dim1, dim2

   !Open file
   call nc_check(nf_open(infile, nf_nowrite, ncid))
     
   !Read longitude
   call nc_check(nf_inq_varid(ncid, "lon", lonid))
   call nc_check(nf_get_vara_double(ncid, lonid, 1, nlon, xlon))

   !Read latitude
   call nc_check(nf_inq_varid(ncid, "lat", latid))
   call nc_check(nf_get_vara_double(ncid, latid, 1, nlat, xlat))

   !Read time

   call nc_check(nf_inq_varid(ncid, "time", timeid))
   call nc_check(nf_get_vara_double(ncid, timeid, 1, ntime, xtime))

   !Read global attributes
!   call nc_check(nf_get_att_text(ncid, nf_global, "title", nctitle))
!   call nc_check(nf_get_att_text(ncid, nf_global, "echam_version", ncechamver))
!   call nc_check(nf_get_att_text(ncid, nf_global, "date_time", ncdate))


   call nc_check(nf_inq_varid(ncid, varname, varid))

   count(1) = nlon
   count(2) = nlat
   count(3) = ntime
   start(1) = 1
   start(2) = 1
   start(3) = 1

!   do timestep = 1, ntime
!      start(3) = timestep
      call nc_check(nf_get_vara_double(ncid, varid, start, count, var3))
!   end do
    
   !Close file
   call nc_check(nf_close(ncid))

end subroutine read_nc3




subroutine write_nc (outfile, nlon, nlat, nmlev, ntime, xlon, xlat, ilev, xtime,  varname, var)
   implicit none

   include 'netcdf.inc'

   integer                              :: nlon, nlat, nmlev, ntime
   real(dp), dimension(nlon)        :: xlon
   real(dp), dimension(nlat)        :: xlat 
   real(dp), dimension(ntime)       :: xtime
   integer, dimension(nmlev)            :: ilev
   real(dp), dimension(nlon,nlat,nlev,ntime) :: var
   character(len=*), intent(in) 	:: outfile, varname

   character(len=256) 			:: nctitle, ncechamver, ncdate 
   character(len=12) 			:: wdate, wtime
   character(len=25) 			:: timestr
   
   integer 				:: vtime(8)

   
!   real(dp), allocatable, dimension(:) ::  lon
!   real(dp), allocatable, dimension(:) ::  lat
!   real(dp), allocatable, dimension(:) ::  mlev
!   real(dp), allocatable, dimension(:) ::  time

   integer, dimension(4) :: start, count, vardims

   integer :: ncstatus, alloc_status, dealloc_status, ncomode
   integer :: ncid, latdimid, londimid, timedimid, timestep, mlevdimid
   integer :: lonid, latid, timeid, mlevid, varid
!   integer :: dim1, dim2, dim3

   ncid = 0
   !Fetch time
   call get_time(timestr)

   !Create netcdf file
   call nc_check(nf_create(outfile, nf_clobber, ncid))

   !Define global attributes
   call nc_check(nf_put_att_text(ncid, nf_global, "title", len(nctitle), nctitle))
   call nc_check(nf_put_att_text(ncid, nf_global, "echam_version", len(ncechamver), ncechamver))
   call nc_check(nf_put_att_text(ncid, nf_global, "date_time", len(trim(timestr)), trim(timestr)))

   !Define grid
   call nc_check(nf_def_dim(ncid, 'lon', nlon, londimid))
   call nc_check(nf_def_dim(ncid, 'lat', nlat, latdimid))
   call nc_check(nf_def_dim(ncid, 'lev', nmlev, mlevdimid))
   call nc_check(nf_def_dim(ncid, 'time', nf_unlimited, timedimid))

   vardims(1) = londimid
   vardims(2) = latdimid
   vardims(3) = mlevdimid
   vardims(4) = timedimid

   !Define variables
   call nc_check(nf_def_var(ncid, 'lon', nf_double, 1, londimid, lonid))
   call nc_check(nf_put_att_text(ncid, lonid, "long_name", &
        &len("longitude"), "longitude"))
   call nc_check(nf_put_att_text(ncid, lonid, "units", &
        &len("degrees_east"), "degrees_east"))

   call nc_check(nf_def_var(ncid, 'lat', nf_double, 1, latdimid, latid))
   call nc_check(nf_put_att_text(ncid, latid, "long_name", &
        &len("latitude"), "latitude"))
   call nc_check(nf_put_att_text(ncid, latid, "units", &
        &len("degrees_north"), "degrees_north"))

   call nc_check(nf_def_var(ncid, 'lev', nf_int, 1, mlevdimid, mlevid))
   call nc_check(nf_put_att_text(ncid, mlevid, "long_name", &
        &len("hybrid level at layer midpoints"), "hybrid level at layer midpoints"))
   call nc_check(nf_put_att_text(ncid, mlevid, "units", &
        &len("level"), "level"))

   call nc_check(nf_def_var(ncid, 'time', nf_double, 1, timedimid, timeid))
   call nc_check(nf_put_att_text(ncid, timeid, "units", &
        &len("1997-02-01 00:00"), "1997-02-01 00:00"))
   call nc_check(nf_put_att_text(ncid, timeid, "calendar", &
        &len("gregorian"), "gregorian"))
   
   call nc_check(nf_def_var(ncid, varname, nf_double, 4, vardims, varid))
   call nc_check(nf_put_att_text(ncid, varid, "long_name", len(varname), varname))
   call nc_check(nf_put_att_text(ncid, varid, "unit", len("hPa"), "hPa"))
   call nc_check(nf_put_att_text(ncid, varid, "grid_type", len("gaussian"), "gaussian"))

  

   call nc_check(nf_enddef(ncid))

   call nc_check(nf_set_fill(ncid, nf_nofill, ncomode))

   !Write grid
   call nc_check(nf_put_var_double(ncid, lonid, xlon))
   call nc_check(nf_put_var_double(ncid, latid, xlat))
   call nc_check(nf_put_var_int(ncid, mlevid, ilev))
   call nc_check(nf_put_vara_double(ncid, timeid, 1, ntime, xtime))

!These settings tell netcdf to write one timestep of data. (The
!setting of start(3) inside the loop below tells netCDF which
!timestep to write.)

   count(1) = nlon
   count(2) = nlat
   count(3) = nmlev
   count(4) = ntime
   start(1) = 1
   start(2) = 1
   start(3) = 1
   start(4) = 1

!Write the pretend data. This will write our surface pressure data.
!The arrays only hold one timestep worth
!of data. We will just rewrite the same data for each timestep. In
!a real(sp) application, the data would change between timesteps.

!   do timestep = 1, ntime
!      start(4) = timestep
!      call nc_check(nf_put_vara_double(ncid, timeid, timestep, 1, xtime))
      call nc_check(nf_put_vara_double(ncid, varid, start, count, var))
!   end do

   !Close file
   call nc_check(nf_close(ncid))

end subroutine write_nc



!Create time string
subroutine get_time(string)

  implicit none

  character (len = *), intent(out) :: string

  integer :: vtime(8)

  call date_and_time(values = vtime)

  write(string,100) vtime(1),vtime(2),vtime(3),vtime(5),vtime(6),vtime(7)

  100 format(I4.4,I2.2,I2.2,' ',I2.2,I2.2,I2.2)

  return

end subroutine get_time

subroutine write_nc3 (outfile, nlon, nlat, ntime, xlon, xlat, xtime, varname, var3)

   implicit none

   include 'netcdf.inc'
   
   integer                              :: nlon, nlat, ntime
   real(dp), dimension(nlon)        :: xlon
   real(dp), dimension(nlat)        :: xlat
   real(dp), dimension(ntime)       :: xtime
   real(dp), dimension(nlon,nlat,ntime) :: var3
   character(len=*), intent(in) 	:: outfile, varname

   character(len=256) 			:: nctitle, ncechamver, ncdate 
   character(len=12) 			:: wdate, wtime
   character(len=25) 			:: timestr
   
   integer 				:: vtime(8)

   
!   real(dp), allocatable, dimension(:) ::  lon
!   real(dp), allocatable, dimension(:) ::  lat
!   real(dp), allocatable, dimension(:) ::  mlev
!   real(dp), allocatable, dimension(:) ::  time

   integer, dimension(3) :: start, count, vardims

   integer :: ncstatus, alloc_status, dealloc_status, ncomode
   integer :: ncid, latdimid, londimid, timedimid, timestep
   integer :: lonid, latid, timeid, varid
   integer :: dim1, dim2

   ncid = 0
   !Fetch time
   call get_time(timestr)

   !Create netcdf file
   call nc_check(nf_create(outfile, nf_clobber, ncid))

   !Define global attributes
   call nc_check(nf_put_att_text(ncid, nf_global, "title", len(nctitle), nctitle))
   call nc_check(nf_put_att_text(ncid, nf_global, "echam_version", len(ncechamver), ncechamver))
   call nc_check(nf_put_att_text(ncid, nf_global, "date_time", len(trim(timestr)), trim(timestr)))

   !Define grid
   call nc_check(nf_def_dim(ncid, 'lon', nlon, londimid))
   call nc_check(nf_def_dim(ncid, 'lat', nlat, latdimid))
   call nc_check(nf_def_dim(ncid, 'time', nf_unlimited, timedimid))

   vardims(1) = londimid
   vardims(2) = latdimid
   vardims(3) = timedimid

   !Define variables
   call nc_check(nf_def_var(ncid, 'lon', nf_double, 1, londimid, lonid))
   call nc_check(nf_put_att_text(ncid, lonid, "long_name", &
        &len("longitude"), "longitude"))
   call nc_check(nf_put_att_text(ncid, lonid, "units", &
        &len("degrees_east"), "degrees_east"))

   call nc_check(nf_def_var(ncid, 'lat', nf_double, 1, latdimid, latid))
   call nc_check(nf_put_att_text(ncid, latid, "long_name", &
        &len("latitude"), "latitude"))
   call nc_check(nf_put_att_text(ncid, latid, "units", &
        &len("degrees_north"), "degrees_north"))

   call nc_check(nf_def_var(ncid, 'time', nf_double, 1, timedimid, timeid))
   call nc_check(nf_put_att_text(ncid, timeid, "units", &
        &len("1997-02-01 00:00"), "1997-02-01 00:00"))
   call nc_check(nf_put_att_text(ncid, timeid, "calendar", &
        &len("gregorian"), "gregorian"))

   call nc_check(nf_def_var(ncid, varname, nf_double, 3, vardims, varid))
   call nc_check(nf_put_att_text(ncid, varid, "long_name", len(varname), varname))
   call nc_check(nf_put_att_text(ncid, varid, "unit", len("hPa"), "hPa"))
   call nc_check(nf_put_att_text(ncid, varid, "grid_type", len("gaussian"), "gaussian"))

  

   call nc_check(nf_enddef(ncid))

   call nc_check(nf_set_fill(ncid, nf_nofill, ncomode))


   !Write grid
   call nc_check(nf_put_var_double(ncid, lonid, xlon))
   call nc_check(nf_put_var_double(ncid, latid, xlat))
   call nc_check(nf_put_vara_double(ncid, timeid, 1, ntime, xtime))

!These settings tell netcdf to write one timestep of data. (The
!setting of start(3) inside the loop below tells netCDF which
!timestep to write.)

   count(1) = nlon
   count(2) = nlat
   count(3) = ntime
   start(1) = 1
   start(2) = 1
   start(3) = 1


!Write the pretend data. This will write our surface pressure data.
!The arrays only hold one timestep worth
!of data. We will just rewrite the same data for each timestep. In
!a real(sp) application, the data would change between timesteps.

!print *,'ntime',ntime
!   do timestep = 1, ntime
!      start(3) = timestep
      call nc_check(nf_put_vara_double(ncid, timeid, 1, ntime, xtime))
      call nc_check(nf_put_vara_double(ncid, varid, start, count, var3))
!   end do

   !Close file
   call nc_check(nf_close(ncid))



end subroutine write_nc3

!Deallocate variable 3 dimension
subroutine dealloc_var3d(varname, var)

  implicit none

  real(sp), allocatable,  intent(inout), dimension(:,:) ::  var
  
  character (len = *), intent(in) :: varname

  integer :: dealloc_status
 
  deallocate(var, stat=dealloc_status)

  if(dealloc_status .ne. 0) then
      write(*,'("Error while deallocating memory: ",A)') varname
      stop
  endif
  
  return
 
end subroutine dealloc_var3d



END PROGRAM radiation


