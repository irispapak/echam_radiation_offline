!>
!! @brief Prepares and provides cloud optical properties
!!
!! @remarks
!!   This code makes use of an outdated (ECHAM4) band structures which are then
!!   artificially mapped to the band structure of the given radiation scheme.
!!
!! @author Bjorn Stevens, MPI-M, Hamburg (2009-09-19)
!!
!! $ID: n/a$
!!
!! @par Origin
!!   Rewrite and synthesis of ECHAM5 code, particularly rad_int.f90 whose
!!   contributers included:  M.A. Giorgetta, MPI-M (2002-05); U. Schulzweida,
!!   MPI-M (2002-05); P. Stier MPI-M \& Caltech (2004-04, 2006-07), M. Thomas
!!   MPI-M (2007-06); U. Schlese, MPI-M (2007-06); M. Esch, MPI-M (2007-06); 
!!   S.J. Lorenz, MPI-M (2007-11).
!! 
!! @par Copyright
!!   2002-2009 by the Deutsche Wetterdienst (DWD) and the Max-Planck-Institut
!!   for Meteorology (MPI-M).  This software is provided for non-commerical 
!!   use only.  See the LICENSE and the WARRANTY conditions
!! 
!! @par License
!!   The use of ICON is hereby granted free of charge for an unlimited time,
!!   provided:
!!   <ol>
!!    <li> Its use is limited to own non-commercial and non-violent purposes;
!!    <li> The code is not re-distributed without the consent of DWD and MPI-M;
!!    <li> This header appears in all copies of the code;
!!    <li> You accept the warranty conditions (see WARRANTY).
!!   </ol>
!!   Commericial use of the code is allowed subject to a separate licensing 
!!   agreement with the DWD and MPI-M
!!
!! @par Warranty
!!   This code is distributed in the home that it will be useful, but WITHOUT
!!   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
!!   FITNESS FOR A PARTICULAR PURPOSE.  
!
MODULE mo_cloud_optics

  USE mo_kind,           ONLY: wp
  USE mo_constants,      ONLY: api,rhoh2o
!LT
!  USE mo_exception,      ONLY: finish

!  USE mo_cosp_simulator, ONLY: cosp_reffl, cosp_reffi, locosp

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: setup_cloud_optics, cloud_optics

  INTEGER, PARAMETER :: jpband=16

  REAL (wp), PARAMETER ::    &
       ccwmin = 1.e-7_wp,    & !< min condensate for lw cloud opacity
       reimn  = 10.0_wp,     & !< min ice effective radius (microns)
       reimx  = 150.0_wp,    & !< max ice effective radius (microns)
       relmn  = 4.0_wp,      & !< min liquid effective radius (microns)
       relmx  = 24.0_wp,     & !< max liquid effective radius (microns)
       zkap_cont = 1.143_wp, & !< continental (Martin et al. ) breadth parameter
       zkap_mrtm = 1.077_wp    !< maritime (Martin et al.) breadth parameter

  LOGICAL, SAVE   :: l_variable_inhoml = .TRUE.

  REAL (wp), SAVE :: &
       zasic,        & !< scaling constant for ice-cloud asymmetry factor 
       zinpar,       & !< exponent for variable lquid-cloud inohomgeneity
       zinhomi,      & !< ice-cloud inohomogeneity factor
       cnst_zinhoml    !< constant value for liquid-cloud inhomogeneity factor

CONTAINS
  !-----------------------------------------------------------------------------
  !>
  !! @brief sets resolution dependent parameters for cloud optics
  !
  SUBROUTINE setup_cloud_optics(lcouple, lipcc, ntrunc, etah)

    LOGICAL, INTENT (IN)   :: & 
         lcouple,             & !< logical for coupled model
         lipcc                  !< logical for ipcc version of model
    INTEGER, INTENT (IN)   :: &
         ntrunc                 !< spectral truncation (triangular)
    REAL (wp), INTENT (IN) :: &
         etah(:)                !< normalized vertical coordinate

    LOGICAL :: low_vertical_resolution
    INTEGER :: kk

    zasic   = 0.85_wp
    zinpar  = 0.12_wp
    zinhomi = 0.90_wp
    !
    ! -- define low vertical resolution as tenth model level from the surface
    !    being at a pressure less than half that of the surface pressure
    !
    kk = max(1,size(etah) - 10)
!    write(*,*) 'kk:'
!    write(*,*) kk
!    write(*,*) 'etah(kk):'
!    write(*,*) etah(kk)

    low_vertical_resolution = etah(kk) < 0.5
    
    IF (low_vertical_resolution) THEN
      IF(ntrunc <= 31) THEN
        zinpar  = 0.06_wp
        zinhomi = 0.80_wp
      ELSE IF(ntrunc == 42) THEN
        zinpar  = 0.07_wp
        zinhomi = 0.85_wp
      ELSE 
        zinpar  = 0.08_wp
        zinhomi = 0.85_wp
      ENDIF
    ELSE 
      IF(ntrunc <= 42) THEN
        zinpar  = 0.08_wp
        zinhomi = 0.85_wp
      ELSE IF(ntrunc >= 319) THEN
        zinhomi = 1.00_wp
        zasic   = 0.91_wp
      ELSE
        zinpar  = 0.10_wp
        zinhomi = 0.85_wp
      ENDIF
    ENDIF
    !
    ! --- slightly different tuning constants for coupled model
    !
    IF(lcouple .OR. lipcc) THEN
      zasic        = 0.85_wp
      cnst_zinhoml = 0.70_wp
      IF (low_vertical_resolution) THEN
        zinhomi    = 0.70_wp 
      ELSE
        zinhomi    = 1.00_wp 
      END IF
      l_variable_inhoml = .FALSE.
    ELSE
      l_variable_inhoml = .TRUE.
    END IF

  END SUBROUTINE setup_cloud_optics
  !-----------------------------------------------------------------------------
  !>
  !! @brief Calculates cloud optical from cloud physical properties
  !!
  !! @remarks  
  !!   Currently this model assumes four bands in the SW and maps these to the
  !!   six-band model of ECHAM5.
  !
  SUBROUTINE cloud_optics( &
          & laglac       ,laland          ,kproma         ,kbdim         ,&
          & klev         ,ktype           ,nb_lw          ,nb_sw         ,&
          & diff         ,icldlyr         ,zlwp           ,ziwp          ,&
          & zlwc         ,ziwc            ,zcdnc          ,tau_lw        ,&
          & tau_sw       ,omg             ,asy                           ,&
          & zradip       ,zradlp )

    INTEGER, INTENT (IN)    :: &
         kproma,               & !< actual length of grid point block
         kbdim,                & !< maximum length of grid point block
         klev,                 &
         ktype(kbdim),         & !< type of convection
         nb_sw,                & 
         nb_lw,                &
         icldlyr(kbdim,klev)
    LOGICAL, INTENT (IN)    :: & 
         laglac(kbdim),        & !< logical for glacier points
         laland(kbdim)           !< logical for land points
    REAL (wp), INTENT (IN)  :: &
         diff,                 &
         zlwp(kbdim,klev),     & !< liquid water path
         ziwp(kbdim,klev),     & !< ice water path
         zcdnc(kbdim,klev),    & !< cloud drop number concentration
         zlwc(kbdim,klev),     & !< liquid water content
         ziwc(kbdim,klev)        !< ice water content

    REAL (wp), INTENT (OUT) ::        &
         tau_lw(kbdim,klev,nb_lw),   & !< LW optical depth
         tau_sw(kbdim,nb_sw,klev),   & !< SW optical depth
         omg(kbdim,nb_sw,klev),      & !< cloud single scattering albedo
         asy(kbdim,nb_sw,klev),      & !< cloud asymmetry factor
         zradlp(kbdim,klev),         & !< effective radius of liquid
         zradip(kbdim,klev)            !< effective radius of ice

    LOGICAL, PARAMETER   :: l_ice_optics = .TRUE.
    INTEGER, PARAMETER   :: n_mdl_bnds = 4,                   &
         & i14band_map(14) = (/4,4,3,3,3,3,3,2,2,1,1,1,1,4/), &
         & i06band_map(6)  = (/1,1,1,2,3,4/)

    REAL (wp), PARAMETER :: &
         a4l0(4) = (/  1.8362_wp,  2.0731_wp,  1.8672_wp,  1.07870_wp/), &
         a4l1(4) = (/ -1.0665_wp, -1.1079_wp, -1.0420_wp, -0.79772_wp/), &
         a4i0(4) = (/  1.9787_wp,  2.1818_wp,  1.9608_wp,  1.25580_wp/), &
         a4i1(4) = (/ -1.0365_wp, -1.0611_wp, -1.0212_wp, -0.88622_wp/)

    REAL (wp), PARAMETER :: &
         b4l0(4) = (/ 1.00000_wp,   1.00000_wp,   0.99936_wp,   0.90032_wp /), &
         b4l1(4) = (/-2.2217e-7_wp,-1.6712e-5_wp,-0.0013632_wp,-0.091955_wp/), &
         b4i0(4) = (/ 1.00000_wp,   0.99999_wp,   0.99975_wp,   0.89779_wp /), &
         b4i1(4) = (/-1.1430e-7_wp,-7.9238e-6_wp,-0.0016620_wp,-0.080200_wp/), &
         b4i2(4) = (/ 0.00000_wp,   0.00000_wp,   6.9726e-6_wp, 0.00000_wp /)

    REAL (wp), PARAMETER :: &
         c4l0(4) = (/ 0.780630_wp, 0.741020_wp, 0.707300_wp, 0.705540_wp/), &
         c4l1(4) = (/ 0.126000_wp, 0.163150_wp, 0.182990_wp, 0.887980_wp/), &
         c4l2(4) = (/-0.042412_wp,-0.050268_wp,-0.045693_wp,-1.821400_wp/), &
         c4l3(4) = (/ 0.000000_wp, 0.000000_wp, 0.000000_wp, 1.577500_wp/), &
         c4l4(4) = (/ 0.000000_wp, 0.000000_wp, 0.000000_wp,-0.462930_wp/), &
         c4i0(4) = (/ 0.796020_wp, 0.771760_wp, 0.746910_wp, 0.776140_wp/), &
         c4i1(4) = (/ 0.101830_wp, 0.119950_wp, 0.135140_wp, 0.151860_wp/), &
         c4i2(4) = (/-0.028648_wp,-0.030557_wp,-0.027140_wp,-0.031452_wp/)

    REAL (wp), PARAMETER, DIMENSION(jpband) ::rebcug = (/ &
         &  0.718_wp, 0.726_wp, 1.136_wp, 1.320_wp,       &
         &  1.505_wp, 1.290_wp, 0.911_wp, 0.949_wp,       &
         &  1.021_wp, 1.193_wp, 1.279_wp, 0.626_wp,       &
         &  0.647_wp, 0.668_wp, 0.690_wp, 0.690_wp     /)

    REAL (wp), PARAMETER, DIMENSION(jpband) ::rebcuh = (/ &
         &  0.0069_wp, 0.0060_wp, 0.0024_wp, 0.0004_wp,   &
         & -0.0016_wp, 0.0003_wp, 0.0043_wp, 0.0038_wp,   &
         &  0.0030_wp, 0.0013_wp, 0.0005_wp, 0.0054_wp,   &
         &  0.0052_wp, 0.0050_wp, 0.0048_wp, 0.0048_wp /)

    INTEGER  :: iband, jk, jl, kk, k
    REAL(wp) :: zip, zlp, ztol, ztoi, zol, zoi, zgl, zgi

    REAL(wp) :: &
         zmsald,             & !< for LW cloud liquid optical depth
         zmsaid,             & !< for LW cloud ice optical depth
         zmacl(kbdim,klev),  & !< for LW cloud liquid optical depth
         zmaci(kbdim,klev),  & !< for LW cloud ice optical depth
         zkap(kbdim),        & !< breath parameter for scaling effective radius
         zlwpt(kbdim),       & !< liquid water path
         zinhoml(kbdim),     & !< cloud inhomogeneity factor (liquid)
         zscratch,           & ! 
         ztau_sw(kbdim,n_mdl_bnds,klev), & !< SW optical depth
         zomg(kbdim,n_mdl_bnds,klev),    & !< cloud single scattering albedo
         zasy(kbdim,n_mdl_bnds,klev)       !< cloud asymmetry factor
    !
    ! 1.0 Basic cloud properties
    ! --------------------------------
    IF (l_variable_inhoml) THEN
      zlwpt(1:kproma) = 0.0_wp
      DO jk = 1, klev
        DO jl = 1, kproma
          zlwpt(jl) = zlwpt(jl)+zlwp(jl,jk)
        END DO
      END DO
      WHERE (zlwpt(1:kproma) > 1.0_wp) 
        zinhoml(1:kproma) = zlwpt(1:kproma)**(-zinpar)
      ELSEWHERE
        zinhoml(1:kproma) = 1.0_wp
      END WHERE
    ELSE
      DO jl = 1, kproma
        IF(ktype(jl) .EQ. 0) THEN
!LT: changed from 0.8 to 0.7 
          zinhoml(jl) = 0.70_wp
        ELSE
          zinhoml(jl) = 0.60_wp
        END IF
      END DO
    END IF

    WHERE (laland(1:kproma).AND.(.NOT.laglac(1:kproma))) 
      zkap(1:kproma)=zkap_cont ! continental breadth factor
    ELSEWHERE
      zkap(1:kproma)=zkap_mrtm ! maritime breadth factor 
    END WHERE

    zscratch = 1.0e6_wp*(3.0e-9_wp/(4.0_wp*api*rhoh2o))**(1.0_wp/3.0_wp)
    DO jk=1,klev
      DO jl=1,kproma
        zradip(jl,jk) = MAX(reimn,MIN(reimx,83.8_wp*ziwc(jl,jk)**0.216_wp))
        zradlp(jl,jk) = MAX(relmn,MIN(relmx,zscratch*zkap(jl)*(zlwc(jl,jk)    &
             /zcdnc(jl,jk))**(1.0_wp/3.0_wp)))
        zmacl(jl,jk)=0.025520637_wp+0.2854650784_wp*EXP(-0.088968393014_wp    &
             *zradlp(jl,jk))
        zmaci(jl,jk)=0.020219423_wp+0.2058619832_wp*EXP(-0.067631070625_wp    &
             *zradip(jl,jk))
      END DO
    END DO
    !
    ! 2.0 Cloud LW optical depth in different bands
    !     Ice cloud emissivity after Ebert and Curry (1992) with diffusivity
    !     factor  diff for ice_optics otherwise ice cloud emissivity after
    !     Rockel et al. (1991)
    ! --------------------------------
    DO iband=1,nb_lw
      DO jk=1,klev
        DO jl=1,kproma
          IF (zlwp(jl,jk)+ziwp(jl,jk)>ccwmin) THEN
            zmsald=zmacl(jl,jk)
            IF (l_ice_optics) THEN
              zmsaid=(rebcuh(iband)+rebcug(iband)/zradip(jl,jk))*diff
            ELSE
              zmsaid=zmaci(jl,jk)
            END IF
            tau_lw(jl,jk,iband) = zmsald*zlwp(jl,jk)*zinhoml(jl)              &
                 &                 + zmsaid*ziwp(jl,jk)*zinhomi
          ELSE
            tau_lw(jl,jk,iband)=0.0_wp
          END IF
!        write(*,*) 'ccwmin:'
!        write(*,*) ccwmin
!        write(*,*) 'zinhoml:'
!        write(*,*) zinhoml(jl)
!        write(*,*) 'zinhomi:'
!        write(*,*) zinhomi
!        write(*,*) 'cld_tau_lw:'
!        write(*,*) tau_lw(jl,jk,iband)
        END DO
      END DO
    END DO
    !
    ! 3.0 Cloud Optical Properties (optical depth, single scattering albedo 
    !     (omega) and Asymmetry Parameter (sometimes called gamma) for SW calcs
    ! --------------------------------
    ztau_sw(1:kproma,1:n_mdl_bnds,1:klev)  = 0.0_wp
    zomg(1:kproma,1:n_mdl_bnds,1:klev)     = 1.0_wp
    zasy(1:kproma,1:n_mdl_bnds,1:klev)     = 0.0_wp

    DO iband = 1,n_mdl_bnds
      DO jk=1,klev
        DO jl=1,kproma
          IF (icldlyr(jl,jk)==1 .AND. (zlwp(jl,jk)+ziwp(jl,jk))>ccwmin) THEN

            zlp = LOG10(zradlp(jl,jk))
            zip = LOG10(zradip(jl,jk))

            ztol = zlwp(jl,jk)*a4l0(iband)*zradlp(jl,jk)**a4l1(iband)
            ztoi = ziwp(jl,jk)*a4i0(iband)*zradip(jl,jk)**a4i1(iband)
            ztau_sw(jl,iband,jk) = ztol*zinhoml(jl) + ztoi*zinhomi

            IF (iband < 4) THEN
              zol    = b4l0(iband) + zradlp(jl,jk)*b4l1(iband)
              zoi    = b4i0(iband) + zradip(jl,jk)*(b4i1(iband)+zradip(jl,jk)  &
                   &                               *b4i2(iband))
            ELSE
              zol    = b4l0(iband) * zradlp(jl,jk)**b4l1(iband)
              zoi    = b4i0(iband) * zradip(jl,jk)**b4i1(iband)
            END IF

            zscratch = (ztol*zol+ztoi*zoi)
            zomg(jl,iband,jk) = zscratch/(ztol+ztoi)
            zgl = c4l0(iband) + zlp*(c4l1(iband) + zlp*(c4l2(iband)            &
                 &                 + zlp*(c4l3(iband) + zlp*c4l4(iband))))
            zgi = zasic*(c4i0(iband) + zip*(c4i1(iband) + zip*c4i2(iband)))
            zasy(jl,iband,jk) = (ztol*zol*zgl+ztoi*zoi*zgi)/zscratch
          END IF
        END DO
      END DO
    END DO
!LT
!    IF (nb_lw /= jpband) CALL finish('cloud_optics lw','band mismatch')
    SELECT CASE (nb_sw)
    CASE(6)
      DO k=1,nb_sw
        kk = i06band_map(k)
        tau_sw(1:kproma,k,1:klev) = ztau_sw(1:kproma,kk,1:klev)
        omg(1:kproma,k,1:klev)    = zomg(1:kproma,kk,1:klev)
        asy(1:kproma,k,1:klev)    = zasy(1:kproma,kk,1:klev)
      end do
    CASE(14)
      DO k=1,nb_sw
        kk = i14band_map(k)
        tau_sw(1:kproma,k,:) = ztau_sw(1:kproma,kk,1:klev)
        omg(1:kproma,k,:)    = zomg(1:kproma,kk,1:klev)
        asy(1:kproma,k,:)    = zasy(1:kproma,kk,1:klev)
      END DO
    CASE DEFAULT
!      CALL finish('cloud_optics','no map to bands')
    END SELECT

  END SUBROUTINE cloud_optics

END MODULE mo_cloud_optics
