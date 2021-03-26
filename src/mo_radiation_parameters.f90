!++mgs: new module 14.03.2010
!++mgs: added decl_sun_cur (for MOZ photolysis) 02.06.2010
!>
!! @brief Module to provide parameters to radiation routines and avoid circular dependencies.
!!
!! @remarks
!!   This module contains the public parameters provided by the radiation module
!!   mo_radiation.
!!
!! @author Bjorn Stevens, MPI-M, Hamburg (2009-09-19):
!!
!!         Martin Schultz, FZJ, Juelich (2010-04-13):
!!              extracted parameters from mo_radiation
!!
!! $ID: n/a$
!!
!! @par Origin
!!   see mo_radiation.f90
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
MODULE mo_radiation_parameters

  USE mo_kind,            ONLY: wp, dp
  USE mo_constants,       ONLY: g, api, cpd

IMPLICIT NONE

  PRIVATE

  PUBLIC :: l_srtm, l_lrtm, l_newoptics, ldiur, lradforcing
  PUBLIC :: lyr_perp, yr_perp, nmonth, isolrad, nb_sw
  PUBLIC :: ih2o, ico2, ich4, io3, io2, in2o, icfc, ighg, iaero
  PUBLIC :: co2vmr, ch4vmr, o2vmr, n2ovmr, cfcvmr
  PUBLIC :: co2mmr, ch4mmr, o2mmr, n2ommr            
  !  PUBLIC :: twopi, deg2rad, deg_day, ch4_v, n2o_v
  PUBLIC :: cemiss, diff
  PUBLIC :: iza, za0, input, infile, inper, infilex, day, nday, year, nyear      
  PUBLIC :: psct, psctm, ssi_factor, flx_ratio_cur, flx_ratio_rad, decl_sun_cur, radctl 
  ! 1.0 NAMELIST global variables and parameters
  ! --------------------------------
! 
  LOGICAL :: l_srtm      = .TRUE.  !< USE AER RRTM Shortwave Model (else ECHAM5)
  LOGICAL :: l_lrtm      = .TRUE.  !< USE New (V4) LRTM Model (else ECHAM5 RRTM)
  LOGICAL :: l_newoptics = .TRUE.  !< USE New (V4) LRTM Model (else ECHAM5 RRTM)
  LOGICAL :: ldiur       = .TRUE.  !< diurnal cycle
  LOGICAL :: lradforcing(2) = (/.FALSE.,.FALSE./) !< &! switch on/off diagnostic 
             !of instantaneous aerosol solar (lradforcing(1)) and 
             !thermal (lradforcing(2)) radiation forcing
  LOGICAL :: lyr_perp    = .FALSE. !< switch to specify perpetual vsop87 year
  INTEGER :: yr_perp     = -99999  !< year if (lyr_perp == .TRUE.)
  INTEGER :: nmonth      =  0      !< 0=annual cycle; 1-12 for perpetual month
  ! nmonth currently works for zonal mean ozone and the orbit (year 1987) only
!LT
!  INTEGER :: isolrad     =  3      !< mode of solar constant calculation
  INTEGER :: isolrad     =  0      !< mode of solar constant calculation
                                   !< default is rrtm solar constant
  INTEGER :: nb_sw      !< number of shortwave bands, set in setup
  !
  ! --- Switches for radiative agents
  !
  INTEGER :: ih2o  = 1  !< water vapor, clouds and ice for radiation
  INTEGER :: ico2  = 2  !< carbon dioxide
  INTEGER :: ich4  = 3  !< methane
  INTEGER :: io3   = 3  !< ozone
  INTEGER :: io2   = 2  !< molecular oxygen
  INTEGER :: in2o  = 3  !< nitrous oxide
  INTEGER :: icfc  = 2  !< cfc11 and cfc12
  INTEGER :: ighg  = 0  !< greenhouse gase scenario
  INTEGER :: iaero = 2  !< aerosol model
  INTEGER :: input     =   1				!LT
  character(len=160) :: infile 				!LT
  INTEGER :: inper     =   1				!LT
  character(len=160) :: infilex 			!LT
  INTEGER :: day      =   1				!LT
  REAL(dp) :: nday    =   1._dp				!LT
  INTEGER :: year     =   1				!LT
  INTEGER :: nyear    =   1990	   			!LT
  ! Zenith angle 				        ! LT
  INTEGER :: iza       =   1				! LT
  REAL(dp):: za0       =  45._dp			! LT
  !
  ! --- Default gas volume mixing ratios - 1990 values (CMIP5)
  !
!  REAL(wp) :: co2vmr    =  353.9e-06_wp !< CO2
!  REAL(wp) :: ch4vmr    = 1693.6e-09_wp !< CH4
!  REAL(wp) :: o2vmr     =    0.20946_wp !< O2
!  REAL(wp) :: n2ovmr    =  309.5e-09_wp !< N20
   REAL(wp) :: cfcvmr(2) = (/252.8e-12_wp,466.2e-12_wp/)  !< CFC 11 and CFC 12
!take 1850 values:
   REAL(wp) :: co2vmr    =  284.7e-06_wp !< 1850 concentration
   REAL(wp) :: ch4vmr    = 791.0e-09_wp !< 1850 concentration
   REAL(wp) :: o2vmr     =    0.20946_wp !< O2
   REAL(wp) :: n2ovmr    =  275.4e-09_wp !< 1850 concentration
!
  !
  ! 2.0 Non NAMELIST global variables and parameters
  ! --------------------------------
  REAL(wp), PARAMETER :: twopi    = 2.0_wp*api
  REAL(wp), PARAMETER :: deg2rad  = api/180.0_wp
!LT
!  REAL(wp), PARAMETER :: deg_day  = g*ndaylen/cpd
  REAL(wp), PARAMETER :: ch4_v(3) = (/1.25e-01_wp,  683.0_wp, -1.43_wp  /)
  REAL(wp), PARAMETER :: n2o_v(3) = (/1.20e-02_wp, 1395.0_wp, -1.43_wp/)
  !
  ! --- radiative transfer parameters
  !
  REAL(wp), PARAMETER :: cemiss = 0.996_wp  !< LW Emissivity Factor
  REAL(wp), PARAMETER :: diff   = 1.66_wp   !< LW Diffusivity Factor
  !
  !
  !++hs
  REAL(wp) :: psct                          !< local (orbit relative and possibly
  !                                            time dependent) solar constant
  REAL(wp) :: psctm                         !< orbit and time dependent solar constant for radiation time step
  REAL(wp) :: ssi_factor(14)                !< fraction of TSI in the 14 RRTM SW bands
  !--hs
  REAL(wp) :: flx_ratio_cur, flx_ratio_rad
  REAL(wp) :: decl_sun_cur                  !< solar declination at current time step
  !
  ! 3.0 Variables computed by routines in mo_radiation (export to submodels)
  ! --------------------------------
  !
  REAL(wp) :: co2mmr, ch4mmr, o2mmr, n2ommr                ! setup_radiation

!  public solar_parameters

  !
!LT 
!     ------------------------------------------------------------------
NAMELIST /radctl/ &
  nmonth,         &! index for annual cycle or perpetual month experiments
                   !  0      : annual cycle on
                   !  1 - 12 : perpetual month January - December
                   !  !!!!!!!! only with PCMDI-Orbit !!!!!!!!!!!!
  ldiur,          &! true for diurnal cycle on
  isolrad,        &! isolrad = 0 : rrtm solar constant
                   ! isolrad = 1 : dependent spectrally resolved solar constant read from file
                   ! isolrad = 2 : solar constant
                   ! isolrad = 3 : constant for amip runs
  ih2o,           &! ih2o = 0 : no H2O in radiation computation, i.e.
                   !            specific humidity = cloud water = cloud ice = 0
                   ! ih2o = 1 : use prognostic specific humidity, cloud water and cloud ice
  ico2,           &! ico2 = 0 : no CO2 in radiation computation
                   ! ico2 = 1 : use prognostic CO2 mass mixing ratio of tracer co2
                   ! ico2 = 2 : uniform volume mixing ratio co2vmr
                   ! ico2 = 4 : uniform volume mixing ratio in scenario run (ighg)
  ich4,           &! ich4 = 0 : no CH4 in radiation computation
                   ! ich4 = 2 : uniform volume mixing ratio ch4vmr
                   ! ich4 = 3 : troposphere: ch4vmr; decay with elevation above
                   ! ich4 = 4 : uniform volume mixing ratio in scenario run (ighg)
  io3,            &! io3  = 0 : no O3 in radiation computation
                   ! io3  = 2 : spectral climatology, as in  ECHAM4
                   ! io3  = 3 : gridpoint climatology from NetCDF file
                   ! io3  = 4 : gridpoint climatology from IPCC-NetCDF file
  io2,            &! io2  = 0 : no O2 in radiation computation
                   ! io2  = 2 : O2    volume mixing ratio o2vmr
  in2o,           &! in2o = 0 : no N2O in radiation computation
                   ! in2o = 2 : uniform volume mixing ratio n2ovmr
                   ! in2o = 3 : troposphere: n2ovmr; decay with elevation above
                   ! in2o = 4 : uniform volume mixing ratio in scenario run (ighg)
  icfc,           &! icfc = 0 : no CFCs in radiation computation
                   ! icfc = 2 : uniform volume mixing ratios cfcvmr(1:2) for :
                   !               CFC11     CFC12
                   ! icfc = 4 : uniform volume mixing ratios in scenario run (ighg)
  ighg,           &! ighg = 0 : no scenario
                   ! ighg = 1 : scenario A1B
                   ! ighg = 2 : scenario B1
                   ! ighg = 3 : scenario A2
  iaero,          &! iaero= 0 : no aerosols in radiation computation
                   ! iaero= 3 : Kinne aerosol climatology
  co2vmr,         &! CO2 volume mixing ratio for ico2=2
  ch4vmr,         &! CH4 volume mixing ratio for ich4=2,3
  o2vmr,          &! O2  volume mixing ratio for io2=2
  n2ovmr,         &! N2O volume mixing ratio for in2o=2,3
  cfcvmr,         &! CFC volume mixing ratios for icfc=2
  l_lrtm,         &! true: USE New (V4) LRTM Model (else ECHAM5 RRTM)
  l_srtm,         &! true: USE AER RRTM Shortwave Model (else ECHAM5)
  l_newoptics,    &! true: USE new cloud optics
  iza,            &!LT
  za0,            &!LT
  input,          &!LT
  infile,         &!LT
  inper,          &!LT
  infilex,        &!LT
  day,            &!LT
  nday,           &!LT
  year,           &!LT
  nyear           !LT
!     ------------------------------------------------------------------


!contains

  !---------------------------------------------------------------------------
  !>
  !! @brief Scans a block and fills with solar parameters
  !! 
  !! @remarks: This routine calculates the solar zenith angle for each
  !! point in a block of data.  For simulations with no dirunal cycle 
  !! the cosine of the zenith angle is set to its average value (assuming 
  !! negatives to be zero and for a day divided into nds intervals).  
  !! Additionally a field is set indicating the fraction of the day over 
  !! which the solar zenith angle is greater than zero.  Otherwise the field 
  !! is set to 1 or 0 depending on whether the zenith angle is greater or 
  !! less than 1. 
  !
!  SUBROUTINE solar_parameters(decl_sun, dist_sun, time_of_day                  &
!       &                     ,sinlon, sinlat, coslon, coslat                   &
!       &                     ,flx_ratio, cos_mu0, daylght_frc, l_refraction)
!
!    LOGICAL, INTENT(in), OPTIONAL :: l_refraction !< add atmospheric refraction
!    REAL(wp), INTENT(in)  :: &
!         decl_sun,           & !< delination of the sun
!         dist_sun,           & !< distance from the sun in astronomical units
!         time_of_day,        & !< time_of_day (in radians)
!         sinlon(:,:),        & !< sines of longitudes
!         sinlat(:,:),        & !< and latitudes
!         coslon(:,:),        & !< cosines of longitudes
!         coslat(:,:)           !< and latitudes
!    REAL(wp), INTENT(out) :: &
!         flx_ratio,          & !< ratio of actual to average solar constant
!         cos_mu0(:,:),       & !< cos_mu_0, cosine of the solar zenith angle
!         daylght_frc(:,:)      !< daylight fraction (0 or 1) with diurnal cycle
!
!    INTEGER     :: i, j
!    REAL(wp)    :: zen1, zen2, zen3, z1, z2, z3, xx

!    LOGICAL :: ll_refraction = .TRUE. !< default atmospheric refractino
!    LOGICAL, SAVE      :: initialized = .FALSE.
!    INTEGER, PARAMETER :: nds = 128 !< number of diurnal samples
!    REAL (wp), SAVE :: cosrad(nds), sinrad(nds)
!    REAL (wp)       :: xsmpl(nds), xnmbr(nds)

!    IF (.NOT.initialized) THEN
!      DO i = 1, nds
!        xx = twopi*(i-1.0_wp)/nds
!        sinrad(i) = SIN(xx)
!        cosrad(i) = COS(xx)
!      END DO
!      initialized = .TRUE.
!    END IF
!    IF (PRESENT(l_refraction)) ll_refraction = l_refraction

!    flx_ratio = 1.0_wp/dist_sun**2
!    zen1 = SIN(decl_sun)
!    zen2 = COS(decl_sun)*COS(time_of_day)
!    zen3 = COS(decl_sun)*SIN(time_of_day)

!    IF (ldiur) THEN
!      cos_mu0(:,:)     =  zen1*sinlat(:,:)                &
!           &             -zen2*coslat(:,:)*coslon(:,:)    &
!           &             +zen3*coslat(:,:)*sinlon(:,:) 
!      daylght_frc(:,:) = 1.0_wp
!      WHERE (cos_mu0(:,:) < 0.0_wp)
!        cos_mu0(:,:)     = 0.0_wp
!        daylght_frc(:,:) = 0.0_wp
!      END WHERE
!    ELSE
!      DO j = 1, SIZE(cos_mu0,2)
!        DO i = 1, SIZE(cos_mu0,1)

!          z1 =  zen1*sinlat(i,j)
!          z2 = -zen2*coslat(i,j)
!          z3 =  zen3*coslat(i,j)

!          xsmpl(:) = z1 + z2*cosrad(:) + z3*sinrad(:)
!          xnmbr(:) = 1.0_wp
!          WHERE (xsmpl(:) < EPSILON(1.0_wp))
!            xsmpl(:) = 0.0_wp
!            xnmbr(:) = 0.0_wp
!          END WHERE

!          cos_mu0(i,j)     = SUM(xsmpl(:))
!          daylght_frc(i,j) = SUM(xnmbr(:))
!        END DO
!      END DO

!      WHERE (daylght_frc(:,:) > EPSILON(1.0_wp))
!        cos_mu0(:,:)     = cos_mu0(:,:)/daylght_frc(:,:)
!        daylght_frc(:,:) = daylght_frc(:,:)/nds
!      END WHERE
!    END IF

!  END SUBROUTINE solar_parameters


END MODULE mo_radiation_parameters
