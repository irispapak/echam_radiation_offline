!>
!! @brief Module to provide interface to radiation routines. 
!!
!! @remarks
!!   This module contains routines that provide the interface between ECHAM
!!   and the radiation code.  Mostly it organizes and calculates the 
!!   information necessary to call the radiative transfer solvers.
!!
!! @author Bjorn Stevens, MPI-M, Hamburg (2009-09-19): 
!!
!!         Hauke Schmidt, MPI-M, Hamburg (2009-12-18): Few modifications to
!!              allow specific solar irradiance for AMIP-type and preindustrial 
!!              simulations.
!!         Luis Kornblueh, MPI-M, Hamburg (2010-04-06): Never ever use write 
!!              directly 
!!         Martin Schultz, FZJ, Juelich (2010-04-13):
!!              Extracted public parameters into new module mo_radiation_parameters
!!              to avoid circular dependencies in submodels
!!                                      (2010-06-03):
!!              Added submodel calls, decl_sun_cur
!!
!! $ID: n/a$
!!
!! @par Origin
!!   Major segments of this code combines and rewrites (for the ICON standard) 
!!   code previously contained in the ECHAM5 routines rad_int.f90, 
!!   radiation.f90 and prerad.f90.  Modifications were also made to provide
!!   a cleaner interface to the aerosol and cloud properties. Contributors to
!!   the code from which the present routines were derived include:  M. Jarraud,
!!   ECMWF (1983-06); M.A. Giorgetta, MPI-M (2002-05); U. Schulzweida,  MPI-M
!!   (2002-05); P. Stier MPI-M \& Caltech (2004-04, 2006-07), M. Thomas MPI-M 
!!   (2007-06); U. Schlese, MPI-M (2007-06); M. Esch, MPI-M (2007-06); S.J. 
!!   Lorenz, MPI-M (2007-11); T. Raddatz, MPI-M (2006-05); I. Kirchner.
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
MODULE mo_radiation

  USE mo_kind,            ONLY: wp
  USE mo_constants,       ONLY: g, api, cpd, rd, avo, rae, vtmpc1,             &
       &                        amco2, amch4, amn2o, amo3, amo2, amd, amw
  USE mo_control,         ONLY: lcouple, lipcc, lmidatm, nn
!  USE mo_time_base,       ONLY: get_calendar_type, JULIAN
  USE mo_io_units,        ONLY: nnml,nin
  USE mo_param_switches,  ONLY: lrad
!
!  USE mo_time_control,    ONLY: l_trigrad, trigrad,                            &
!       &                        ndaylen, l_orbvsop87, get_orbit_times,         &
!       &                        p_bcast_event, current_date, next_date,        &
!       &                        previous_date, radiation_date,                 &
!       &                        prev_radiation_date,get_date_components,       &
!       &                        lresume, lstart
!  USE mo_convect_tables,  ONLY: prepare_ua_index,lookup_ua
!  USE mo_geoloc,          ONLY: amu0_x, rdayl_x, amu0m_x, rdaylm_x, coslon_2d, &
!       &                        sinlon_2d, sinlat_2d, coslat_2d
!  USE mo_orbit,           ONLY: cecc, cobld, clonp, orbit_kepler, orbit_vsop87
!  USE mo_solar_irradiance,ONLY: get_solar_irradiance, set_solar_irradiance, &
!                                get_solar_irradiance_m, set_solar_irradiance_m
  USE mo_cloud_optics,    ONLY: setup_cloud_optics, cloud_optics
  USE mo_newcld_optics,   ONLY: setup_newcld_optics, newcld_optics
!  USE mo_aero_tanre,      ONLY: aerosol_optics_tanre
  USE mo_aero_kinne,      ONLY: set_aop_kinne  
!  USE mo_aero_volc,       ONLY: add_aop_volc
  USE mo_hyb,             ONLY: cetah
  USE mo_parrrtm,         ONLY: jpband, jpinpx, jpxsec 
  USE mo_echam5_sw,       ONLY: nsw, su_sw4, sw
  USE mo_srtm_config,     ONLY: jpsw, setup_srtm, ssi_default, ssi_amip,       &
       &                        ssi_preind
  USE mo_srtm,            ONLY: srtm_srtm_224gp
  USE mo_lrtm,            ONLY: lrtm
  USE mo_lrtm_setup,      ONLY: lrtm_setup 
  USE mo_radiation_parameters, ONLY: l_srtm, l_lrtm, l_newoptics, ldiur,&
                                     nmonth, isolrad, nb_sw,   &
                                     ih2o, ico2, ich4, io3, io2, in2o, icfc,      &
                                     ighg, iaero, nmonth, iaero,                  &
                                     co2vmr, ch4vmr, o2vmr, n2ovmr, cfcvmr,       &
                                     co2mmr, ch4mmr, o2mmr, n2ommr,               &
                                     cemiss, diff,                                &
                                     psct, psctm, ssi_factor,                     &
                                     flx_ratio_cur, flx_ratio_rad, radctl
  USE mo_exception,       ONLY: message, message_text, finish

  IMPLICIT NONE
  
  PRIVATE
 
  PUBLIC :: setup_radiation, rrtm_interface, pre_radiation, setup_radiation_new
  !
  !
CONTAINS
  !-----------------------------------------------------------------------------
  !>
  !! @brief Prepares information for radiation call
  !
!LT
   SUBROUTINE pre_radiation
!
!    LOGICAL  :: l_rad_call, l_write_solar
!    INTEGER  :: icurrentyear, icurrentmonth, iprevmonth, i
!    REAL(wp) :: rasc_sun, decl_sun, dist_sun, orbit_date, time_of_day, zrae
     REAL(wp) :: solcm
!    !
!    ! 1.0 Compute orbital parameters for current time step
!    ! --------------------------------
!    l_rad_call = .FALSE.
!    CALL get_orbit_times(l_rad_call, lyr_perp, nmonth, yr_perp, time_of_day, &
!         &               orbit_date)
!
!    IF (l_orbvsop87) THEN 
!      CALL orbit_vsop87 (orbit_date, rasc_sun, decl_sun, dist_sun)
!    ELSE
!      CALL orbit_kepler (orbit_date, rasc_sun, decl_sun, dist_sun)
!    END IF
!    decl_sun_cur = decl_sun       ! save for aerosol and chemistry submodels
!    CALL solar_parameters(decl_sun, dist_sun, time_of_day, &
!         &                sinlon_2d, sinlat_2d, coslon_2d, coslat_2d, &
!         &                flx_ratio_cur, amu0_x, rdayl_x)
!     SELECT CASE (isolrad)
!     CASE (0)
!       solc = SUM(ssi_default)
!    CASE (1)
!      CALL get_solar_irradiance(current_date, next_date)
!      CALL set_solar_irradiance(solc)
!    CASE (2)
!      solc = SUM(ssi_preind)
!    CASE (3)
!      solc = SUM(ssi_amip)
!     CASE default
!      WRITE (message_text, '(a,i2,a)') &
!           'isolrad = ', isolrad, ' in radctl namelist is not supported'
!      CALL message('pre_radiation', message_text)
!     END SELECT
!
!
!    !
!    ! 2.0 Prepare time dependent quantities for rad (on radiation timestep)
!    ! --------------------------------
!    IF (lrad .AND. l_trigrad) THEN
!      l_rad_call = .TRUE.
!      CALL get_orbit_times(l_rad_call, lyr_perp, nmonth, yr_perp, time_of_day , &
!           &               orbit_date)
!
!      IF ( l_orbvsop87 ) THEN 
!        CALL orbit_vsop87 (orbit_date, rasc_sun, decl_sun, dist_sun)
!      ELSE
!        CALL orbit_kepler (orbit_date, rasc_sun, decl_sun, dist_sun)
!      END IF
!      CALL solar_parameters(decl_sun, dist_sun, time_of_day, &
!           &                sinlon_2d, sinlat_2d, coslon_2d, coslat_2d, &
!           &                flx_ratio_rad ,amu0m_x, rdaylm_x)
!      ! consider curvature of the atmosphere for high zenith angles
!      zrae = rae*(rae+2.0_wp)
!      amu0m_x(:,:)  = rae/(SQRT(amu0m_x(:,:)**2+zrae)-amu0m_x(:,:))
!      !jsr&hs: For the calculation of radiative transfer, a maximum zenith angle
!      !        of about 84 degrees is applied in order to avoid to much overshooting
!      !        when the extrapolation of the radiative fluxes from night time
!      !        regions to daytime regions is done for time steps at which no
!      !        radiation calculation is performed. This translates into cosines
!      !        of the zenith angle > 0.1.
!      !        This approach limits the calculation of the curvature effect
!      !        above, and may have to be reconsidered when radiation is cleaned
!      !        up.
!      amu0m_x(:,:) = MAX(amu0m_x(:,:),0.1_wp)
!      !
!      ! --- Prepare Ozone climatology
!      !
!      SELECT CASE (io3)
!      CASE (3) 
!        CALL pre_o3clim_3(nmonth)
!      CASE (4) 
!        CALL pre_o3clim_4
!      END SELECT
!      !
!
!    !++jsr&hs
!    ! 3.0 Prepare possibly time dependent total solar and spectral irradiance
!    ! --------------------------------
!    ! ATTENTION: 
!    ! This part requires some further work. Currently, a solar constant of
!    ! 1361.371 is used as default. This is the TSI averaged over the
!    ! years 1979 to 1988, and should be used for AMIP type runs. If lcouple is
!    ! true, a solar constant of 1360.875 is used, the average for the years 1844
!    ! to 1856. This should be used for a preindutrial control run.
!    ! The spectral distribution of this TSI is currently also prescribed for
!    ! these two cases depending on the lcouple switch.
!    ! For transient CMIP5 simulations, the available time
!    ! varying TSI and SSI has to be read in and used here.
!
    SELECT CASE (isolrad)
     CASE (0)
!       solcm = SUM(ssi_default)
!       IF (l_srtm) ssi_factor = ssi_default
!take preindustrial values:
        solcm = SUM(ssi_preind)
        IF (l_srtm) ssi_factor = ssi_preind
!    CASE (1)
!      CALL get_solar_irradiance_m(prev_radiation_date, radiation_date, nb_sw, l_srtm)
!      CALL set_solar_irradiance_m(solcm, ssi_factor, nb_sw, l_srtm)
!    CASE (2)
!      solcm = SUM(ssi_preind)
!      IF (l_srtm) ssi_factor = ssi_preind
!    CASE (3)
!      solcm = SUM(ssi_amip)
!      IF (l_srtm) ssi_factor = ssi_amip
     CASE default
!      WRITE (message_text, '(a,i2,a)') &
!           'isolrad = ', isolrad, ' in radctl namelist is not supported'
!      CALL message('pre_radiation', message_text)
     END SELECT
!LT    psctm = flx_ratio_rad*solcm
     psctm = solcm
     IF (l_srtm) THEN
       ssi_factor(:) = ssi_factor(:)/solcm
     END IF
!
!    ! output of solar constant every month
!
!    CALL get_date_components(current_date, month=icurrentmonth, &
!         year=icurrentyear)
!    CALL get_date_components(previous_date, month=iprevmonth)
!    l_write_solar = icurrentmonth/=iprevmonth
!    IF (l_write_solar .OR. lresume .OR. lstart) THEN
!      CALL message('','')
!      WRITE (message_text,'(a,i0,a,i2.2,a,f6.1)') &
!           'Total solar constant [W/m^2] for ',      &
!           icurrentyear, '-',                        & 
!           icurrentmonth, ' = ', solc
!      CALL message('',message_text)
!      CALL message('','')
!      IF (l_srtm) THEN
!        DO i = 1, nb_sw
!          WRITE (message_text,'(a,i2,a,f7.5)') &
!               '   solar constant fraction: band ', i, &
!               ' = ', ssi_factor(i)
!          CALL message('',message_text)
!        END DO
!      END IF
!    END IF
!    !--jsr&hs
!    END IF
!
   END SUBROUTINE pre_radiation
  !---------------------------------------------------------------------------
  !>
  !! @brief Organizes the calls to the ratiation solver
  !! 
  !! @remarks This routine organises the input/output for the radiation
  !! computation.  The state of radiatively active constituents is set as the
  !! input. Output are flux transmissivities and emissivities at all the half
  !! levels of the grid (respectively ratio solar flux/solar input and ratio 
  !! thermal flux/local black-body flux). This output will be used in radheat
  !! at all time steps until the next full radiation time step.
  !
!LT
!  SUBROUTINE radiation( &
!       &  kproma            ,kbdim           ,klev             ,klevp1        &
!       & ,krow              ,ktrac           ,ktype            ,loland        &
!       & ,loglac            ,alb_vis         ,alb_nir          ,alb_vis_dir   &
!       & ,alb_nir_dir       ,alb_vis_dif     ,alb_nir_dif      ,tk_sfc        &
!       & ,pp_hl             ,pp_fl           ,tk_fl            ,qm_vap        &
!       & ,qm_liq            ,qm_ice          ,pgeom1           ,co2           &
!       & ,cdnc              ,cld_frc          ,pxtm1           ,cld_cvr       &
!       & ,nir_sfc           ,nir_dff_sfc     ,vis_sfc          ,vis_dff_sfc   &
!       & ,dpar_sfc          ,par_dff_sfc     ,emter_toa        ,trsol_toa     &
!       & ,emter_sfc         ,trsol_sfc       ,emter            ,trsol         &
!       & ,ozone                                                               )
!
!    INTEGER, INTENT(IN)  :: kproma, kbdim, klev, klevp1, krow, ktrac
!    INTEGER, INTENT(IN)  :: ktype(kbdim)  !< type of convection
!
!    LOGICAL, INTENT(IN)  :: &
!         loland(:),         & !< land mask
!         loglac(:)            !< glacier mask
!
!    REAL(wp), INTENT(IN) :: &
!         alb_vis(:),        & !< surface albedo for visible range
!                              !< (direct and diffuse combined)
!         alb_nir(:),        & !< surface albedo for NIR range
!                              !< (direct and diffuse combined)
!         alb_vis_dir(:),    & !< surface albedo for visible range and direct light
!         alb_nir_dir(:),    & !< surface albedo for NIR range and direct light
!         alb_vis_dif(:),    & !< surface albedo for visible range and diffuse light
!         alb_nir_dif(:),    & !< surface albedo for NIR range and diffuse light
!         tk_sfc(:),         & !< Surface temperature
!         pp_hl(:,:),        & !< pressure at half levels [Pa]
!         pp_fl(:,:),        & !< Pressure at full levels [Pa]
!         tk_fl(:,:),        & !< Temperature on full levels [K]
!         qm_vap(:,:),       & !< Water vapor mixing ratio 
!         qm_liq(:,:),       & !< Liquid water mixing ratio
!         qm_ice(:,:),       & !< Ice water mixing ratio
!         pgeom1(kbdim,klev),& !< geopotential above ground
!         co2(:,:),          & !< Carbon Dioxide mixing ratio
!         cdnc(:,:),         & !< Cloud drop number concentration
!         cld_frc(:,:),      & !< Cloud fraction
!         pxtm1(kbdim,klev,ktrac) !< tracer mass mixing ratios 
!
!    REAL(wp), INTENT(OUT) ::      &
!         cld_cvr(:),              & !< Cloud cover in a column
!         nir_sfc(:),              & !< SW Near Ifrared 
!         nir_dff_sfc(:),          & !< Near-Infrared diffuse fraction of SW NIR
!         vis_sfc(:),              & !< SW Visible (250-680 nm) 
!         vis_dff_sfc(:),          & !< Near-Infrared diffuse fraction of SW vis
!         dpar_sfc(:),             & !< surf. PAR downw.
!         par_dff_sfc(:),          & !< fraction of diffuse PAR
!         emter_toa(kbdim,2),      & !< TOA terrestrial emissivity
!         trsol_toa(kbdim,2),      & !< TOA solar transmissivity
!         emter_sfc(kbdim,klevp1), & !< Surface terrestrial emissivity
!         trsol_sfc(kbdim,klevp1), & !< Surface solar transmissivity
!         emter(kbdim,klevp1),     & !< Terrestrial emissivity
!         trsol(kbdim,klevp1),     & !< Solar transmissivity
!         ozone(kbdim,klev)          !< Ozone 
!
!    INTEGER  :: jk, jl, idx(kbdim), iaero_f
!
!    REAL(wp) ::                     &
!         cos_mu0(kbdim),            & !< Cos of local zenith angle
!         ppd_hl(kbdim,klev),        & !< pressure diff between half levels [Pa]
!         pp_sfc(kbdim),             & !< surface pressure [Pa}
!         tk_hl(kbdim,klevp1),       & !< Tempeature at half levels [Pa]
!         xq_sat(kbdim,klev),        & !< Saturation mixing ratio 
!         xq_vap(kbdim,klev),        & !< Water vapor mixing ratio
!         xq_liq(kbdim,klev),        & !< Liquid water mixing ratio
!         xq_ice(kbdim,klev),        & !< Ice mixing ratio
!         xm_co2(kbdim,klev),        & !< CO2 mixing ratio
!         xm_o3(kbdim,klev),         & !< Ozone mixing ratio
!         xm_o2(kbdim,klev),         & !< O2 mixing ratio
!         xm_ch4(kbdim,klev),        & !< Methane mixing ratio
!         xm_n2o(kbdim,klev),        & !< Nitrous Oxide mixing ratio
!         xm_cfc(kbdim,klev,2),      & !< CFC mixing ratio
!         emiss_sfc(kbdim),          & !< Sfc emissivity
!         toa_solar_irr(kbdim),      & !< TOA solar irradiance
!         flx_uplw_sfc(kbdim),       & !< Srfc upward lw flux  [Wm2]
!         flx_upsw_sfc(kbdim),       & !< Srfc upward sw flux  [Wm2]
!         flx_uplw_clr_sfc(kbdim),   & !< Srfc upward lw flux (clear sky) [Wm2]
!         flx_upsw_clr_sfc(kbdim),   & !< Srfc upward sw flux (clear sky) [Wm2]
!         flx_dnlw(kbdim,klevp1),    & !< Net dwnwrd LW flux [Wm2]
!         flx_dnsw(kbdim,klevp1),    & !< Net dwnwrd SW flux [Wm2]
!         flx_dnlw_clr(kbdim,klevp1),& !< Net dn LW flux (clear sky) [Wm2]
!         flx_dnsw_clr(kbdim,klevp1),& !< Net dn SW flux (clear sky) [Wm2]
!         dvis_sfc(kbdim),           & !< surface downward visible flux
!         ua(kbdim),                 & !<
!         y1(kbdim)                    !< 1D Scratch Array for diagnostics
! for forcing calculation
!    REAL(wp) ::                       &
!         flx_dnlw_f(kbdim,klevp1),    & !< Net dwnwrd LW flux [Wm2] for forcing
!         flx_dnsw_f(kbdim,klevp1),    & !< Net dwnwrd SW flux [Wm2] for forcing
!         flx_dnlw_clr_f(kbdim,klevp1),& !< Net dn LW flux (clear sky) [Wm2] for forcing
!         flx_dnsw_clr_f(kbdim,klevp1)   !< Net dn SW flux (clear sky) [Wm2] for forcing
!
!    !
!    ! 1.0 calculate variable input parameters (location and state variables)
!    ! --------------------------------
!    ! 
!    ! --- solar zenith angle
!    !
!    !hs lcl_solar_cnst    = flx_ratio_rad*solc     
!    cos_mu0(1:kproma) = amu0m_x(1:kproma,krow)
!    !
!    ! --- Pressure (surface and distance between half levels)
!    !
!    pp_sfc(1:kproma)   = pp_hl(1:kproma,klevp1)
!    ppd_hl(1:kproma,:) = pp_hl(1:kproma,2:klev+1)-pp_hl(1:kproma,1:klev)
!    !
!    ! --- temperature at half levels
!    !
!    DO jk=2,klev
!      DO jl = 1, kproma
!        tk_hl(jl,jk) = (tk_fl(jl,jk-1)*pp_fl(jl,jk-1)*( pp_fl(jl,jk)          &
!             & - pp_hl(jl,jk) ) + tk_fl(jl,jk)*pp_fl(jl,jk)*( pp_hl(jl,jk)    &
!             & - pp_fl(jl,jk-1))) /(pp_hl(jl,jk)*(pp_fl(jl,jk) -pp_fl(jl,jk-1)))
!      END DO
!    END DO
!    DO jl = 1, kproma
!      tk_hl(jl,klevp1) = tk_sfc(jl)
!      tk_hl(jl,1)      = tk_fl(jl,1)-pp_fl(jl,1)*(tk_fl(jl,1) - tk_hl(jl,2))  &
!           &             / (pp_fl(jl,1)-pp_hl(jl,2))
!    END DO
!    !
!    ! --- phases of water substance
!    !
!    xq_vap(1:kproma,:) = MAX(qm_vap(1:kproma,:),EPSILON(1.0_wp))
!
!    DO jk = 1, klev
!      CALL prepare_ua_index('radiation', kproma, tk_fl(:,jk), idx(:))
!      CALL lookup_ua(kproma, idx(:), ua(:))
!      xq_sat(1:kproma,jk) = ua(1:kproma)/pp_fl(1:kproma,jk)
!      xq_sat(1:kproma,jk) = MIN(xq_sat(1:kproma,jk),0.5_wp)
!      xq_sat(1:kproma,jk) = xq_sat(1:kproma,jk)/(1.0_wp-vtmpc1                &
!           &                * xq_sat(1:kproma,jk))
!      xq_sat(1:kproma,jk) = MAX(2.0_wp*EPSILON(1.0_wp),xq_sat(1:kproma,jk))
!    END DO
!    xq_liq(1:kproma,:) = MAX(qm_liq(1:kproma,:),0.0_wp)       ! cloud liquid
!    xq_ice(1:kproma,:) = MAX(qm_ice(1:kproma,:),0.0_wp)       ! cloud ice
!    !
!    ! --- cloud cover
!    ! 
!    cld_cvr(1:kproma) = 1.0_wp - cld_frc(1:kproma,1)
!    DO jk = 2, klev
!      cld_cvr(1:kproma) = cld_cvr(1:kproma)                                    &
!           &        *(1.0_wp-MAX(cld_frc(1:kproma,jk),cld_frc(1:kproma,jk-1))) &
!           &        /(1.0_wp-MIN(cld_frc(1:kproma,jk-1),1.0_wp-EPSILON(1.0_wp)))
!    END DO
!    cld_cvr(1:kproma) = 1.0_wp-cld_cvr(1:kproma)   
!    !
!    ! --- gases
!    ! 
!    xm_co2(1:kproma,:)   = gas_profile(kproma, klev, ico2, gas_mmr = co2mmr,    &
!         &  gas_scenario = ghg_co2mmr, gas_val = co2)
!    xm_ch4(1:kproma,:)   = gas_profile(kproma, klev, ich4, gas_mmr = ch4mmr,    &
!         &  gas_scenario = ghg_ch4mmr, pressure = pp_fl, xp = ch4_v)
!    xm_n2o(1:kproma,:)   = gas_profile(kproma, klev, in2o, gas_mmr = n2ommr,    &
!         &  gas_scenario = ghg_n2ommr, pressure = pp_fl, xp = n2o_v)
!    xm_cfc(1:kproma,:,1) =  gas_profile(kproma, klev, icfc, gas_mmr=cfcvmr(1),  &
!         &  gas_scenario = ghg_cfcvmr(1))
!    xm_cfc(1:kproma,:,2) =  gas_profile(kproma, klev, icfc, gas_mmr=cfcvmr(2),  &
!         &  gas_scenario = ghg_cfcvmr(2))
!    xm_o2(1:kproma,:)    = gas_profile(kproma, klev, io2,  gas_mmr = o2mmr)     
!
!    ozon: SELECT CASE (io3)
!    CASE (0)
!      xm_o3(1:kproma,:) = EPSILON(1.0_wp)
!    CASE (2)
!      xm_o3(1:kproma,:) = o3_lwb(krow,ppd_hl,pp_hl)
!    CASE (3)
!      xm_o3(1:kproma,:) = o3clim(krow,kproma,kbdim,klev,pp_hl,pp_fl)
!    CASE (4)
!      xm_o3(1:kproma,:) = o3clim(krow,kproma,kbdim,klev,pp_hl,pp_fl)
!    CASE default
!      CALL finish('radiation','o3: this "io3" is not supported')
!    END SELECT ozon
!    ozone(1:kproma,:) = xm_o3(1:kproma,:)
!
!    ! 2.0 Call interface to radiation solver
!    ! --------------------------------
!    IF (lradforcing(1).OR.lradforcing(2)) THEN
!    ! call first without aerosols and keep relevant quantities only, then standard call
!       iaero_f=0
!       CALL rrtm_interface( &
!            & iaero_f         ,kproma          ,kbdim           ,klev            ,& 
!            & krow            ,ktrac           ,ktype           ,nb_sw           ,&
!            & loland          ,loglac          ,cemiss          ,diff            ,&
!            & cos_mu0         ,pgeom1          ,alb_vis         ,alb_nir         ,&
!            & alb_vis_dir     ,alb_nir_dir     ,alb_vis_dif     ,alb_nir_dif     ,&
!            & pp_fl           ,pp_hl           ,pp_sfc          ,tk_fl           ,&
!            & tk_hl           ,tk_sfc          ,xq_vap          ,xq_sat          ,&
!            & xq_liq          ,xq_ice          ,cdnc            ,cld_frc         ,&
!            & cld_cvr         ,xm_o3           ,xm_co2          ,xm_ch4          ,&
!            & xm_n2o          ,xm_cfc          ,xm_o2           ,pxtm1           ,&
!            & flx_dnlw_f      ,flx_dnsw_f      ,flx_dnlw_clr_f  ,flx_dnsw_clr_f  ,&
!            & flx_uplw_sfc    ,flx_upsw_sfc    ,flx_uplw_clr_sfc,flx_upsw_clr_sfc,&
!            & emiss_sfc       ,toa_solar_irr   ,nir_sfc         ,nir_dff_sfc     ,&
!            & vis_sfc         ,vis_dff_sfc     ,dvis_sfc        ,dpar_sfc        ,&
!            & par_dff_sfc                                                        )
!       CALL prepare_forcing(kproma           ,kbdim             ,klevp1        &
!                           ,krow             ,flx_dnlw_f        ,flx_dnsw_f    &
!                           ,flx_dnlw_clr_f   ,flx_dnsw_clr_f    ,cos_mu0       )
!    END IF
!    CALL rrtm_interface( &
!         & iaero           ,kproma          ,kbdim           ,klev            ,& 
!         & krow            ,ktrac           ,ktype           ,nb_sw           ,&
!         & loland          ,loglac          ,cemiss          ,diff            ,&
!         & cos_mu0         ,pgeom1          ,alb_vis         ,alb_nir         ,&
!         & alb_vis_dir     ,alb_nir_dir     ,alb_vis_dif     ,alb_nir_dif     ,&
!         & pp_fl           ,pp_hl           ,pp_sfc          ,tk_fl           ,&
!         & tk_hl           ,tk_sfc          ,xq_vap          ,xq_sat          ,&
!         & xq_liq          ,xq_ice          ,cdnc            ,cld_frc         ,&
!         & cld_cvr         ,xm_o3           ,xm_co2          ,xm_ch4          ,&
!         & xm_n2o          ,xm_cfc          ,xm_o2           ,pxtm1           ,&
!         & flx_dnlw        ,flx_dnsw        ,flx_dnlw_clr    ,flx_dnsw_clr    ,&
!         & flx_uplw_sfc    ,flx_upsw_sfc    ,flx_uplw_clr_sfc,flx_upsw_clr_sfc,&
!         & emiss_sfc       ,toa_solar_irr   ,nir_sfc         ,nir_dff_sfc     ,&
!         & vis_sfc         ,vis_dff_sfc     ,dvis_sfc        ,dpar_sfc        ,&
!         & par_dff_sfc                                                         )
!    !
!    ! 3.0 Diagnostics
!    ! --------------------------------
!    !
!    ! --- Total fluxes
!    y1(1:kproma) = (psctm*cos_mu0(1:kproma))
!    trsol(1:kproma,1:klevp1) = flx_dnsw(1:kproma,1:klevp1) / &
!         &                     SPREAD(y1(1:kproma),2,klevp1)
!    emter(1:kproma,:) = flx_dnlw(1:kproma,:)
!    !
!    ! --- fluxes for JSBACH
!    nir_sfc(1:kproma) = nir_sfc(1:kproma) / y1(1:kproma)
!    vis_sfc(1:kproma) = vis_sfc(1:kproma) / y1(1:kproma)
!    IF (l_srtm) THEN
!      dpar_sfc(1:kproma) = dpar_sfc(1:kproma) / y1(1:kproma)
!    ELSE
!      dpar_sfc(1:kproma) = dvis_sfc(1:kproma) / y1(1:kproma)
!      par_dff_sfc(1:kproma) = vis_dff_sfc(1:kproma)
!    ENDIF
!    !
!    ! --- Clear sky fluxes
!    emter_toa(1:kproma,1) = flx_dnlw_clr(1:kproma,1)
!    emter_toa(1:kproma,2) = flx_dnlw_clr(1:kproma,klevp1)
!    emter_sfc(1:kproma,:) = flx_dnlw_clr(1:kproma,:)
!    trsol_toa(1:kproma,1) = flx_dnsw_clr(1:kproma,1) / y1(1:kproma)
!    trsol_toa(1:kproma,2) = flx_dnsw_clr(1:kproma,klevp1) / y1(1:kproma)
!    trsol_sfc(1:kproma,1:klevp1) = flx_dnsw_clr(1:kproma,1:klevp1) / &
!         &                         SPREAD(y1(1:kproma),2,klevp1)
!
!  END SUBROUTINE radiation
!LT
!  !---------------------------------------------------------------------------
!  !>
!  !! GAS_PROFILE:  Determines Gas distributions based on case specification
!  !! 
!  !! @par Revsision History 
!  !! B. Stevens (2009-08). 
!  !! H. Schmidt (2010-08): profile calculation added for scenario case.
!  !!
!  !! Description: This routine calculates the gas distributions for one of
!  !! five cases:  (0) no gas present; (1) prognostic gas; (2) specified 
!  !! mixing ratio; (3) mixing ratio decaying with height given profile;
!  !! (4) scenario run with different mixing ratio, if profile parameters are
!  !! given a vertical profile is calculated as in (3).
!  !
!  FUNCTION gas_profile (kproma, klev, igas, gas_mmr, gas_scenario, gas_mmr_v, &
!       &                gas_scenario_v, gas_val, xp, pressure)
!
!    INTEGER, INTENT (IN) :: kproma, klev, igas
!    REAL (wp), OPTIONAL, INTENT (IN) :: gas_mmr, gas_scenario
!    REAL (wp), OPTIONAL, INTENT (IN) :: pressure(:,:), xp(3)
!    REAL (wp), OPTIONAL, INTENT (IN) :: gas_mmr_v(:,:)
!    REAL (wp), OPTIONAL, INTENT (IN) :: gas_scenario_v(:,:)
!    REAL (wp), OPTIONAL, INTENT (IN) :: gas_val(:,:)
!
!    REAL (wp) :: gas_profile(kproma,klev), zx_d, zx_m
!    LOGICAL :: gas_initialized
!
!    gas_initialized = .FALSE.
!    SELECT CASE (igas)
!    CASE (0)
!      gas_profile(1:kproma,:) = EPSILON(1.0_wp)
!      gas_initialized = .TRUE.
!    CASE (1)
!      IF (PRESENT(gas_val)) THEN
!        gas_profile(1:kproma,:) = MAX(gas_val(1:kproma,:), EPSILON(1.0_wp))
!        gas_initialized = .TRUE.
!      END IF
!    CASE (2)
!      IF (PRESENT(gas_mmr)) THEN
!        gas_profile(1:kproma,:) = gas_mmr
!        gas_initialized = .TRUE.
!      ELSE IF (PRESENT(gas_mmr_v)) THEN
!        gas_profile(1:kproma,:) = gas_mmr_v(1:kproma,:)
!        gas_initialized = .TRUE.
!      END IF
!    CASE (3)
!      IF (PRESENT(gas_mmr) .AND. PRESENT(xp) .AND. PRESENT(pressure)) THEN
!        zx_m = (gas_mmr+xp(1)*gas_mmr)*0.5_wp
!        zx_d = (gas_mmr-xp(1)*gas_mmr)*0.5_wp
!        gas_profile(1:kproma,:)=(1-(zx_d/zx_m)*TANH(LOG(pressure(1:kproma,:)   &
!             &                  /xp(2)) /xp(3))) * zx_m
!        gas_initialized = .TRUE.
!      END IF
!    CASE (4)
!      IF (PRESENT(gas_scenario)) THEN
!        IF (PRESENT(xp) .AND. PRESENT(pressure)) THEN
!          ! comment H. Schmidt: If the respective parameters are present, a vertical 
!          ! profile is calculated as in option (3). This allows a seamless
!          ! continuation of preindustrial control with scenarios. The treatment here is
!          ! inconsistent with having two different options for the constant
!          ! concentration cases (2 without and 3 with profile). However, instead
!          ! of adding a fifth option, it seems more advisable to clean up the 
!          ! complete handling of radiation switches (including ighg), later.
!          zx_m = (gas_scenario+xp(1)*gas_scenario)*0.5_wp
!          zx_d = (gas_scenario-xp(1)*gas_scenario)*0.5_wp
!          gas_profile(1:kproma,:)=(1-(zx_d/zx_m)*TANH(LOG(pressure(1:kproma,:)   &
!             &                  /xp(2)) /xp(3))) * zx_m
!        ELSE
!          gas_profile(1:kproma,:) = gas_scenario
!        ENDIF
!        gas_initialized = .TRUE.
!      ELSE IF (PRESENT(gas_scenario_v)) THEN
!        gas_profile(1:kproma,:) = gas_scenario_v(1:kproma,:)
!        gas_initialized = .TRUE.
!      END IF
!    END SELECT
!    IF (.NOT. gas_initialized) &
!         CALL finish('radiation','gas_profile options not supported')
!
!  END FUNCTION gas_profile
  !-----------------------------------------------------------------------------
  !>
  !! @brief arranges input and calls rrtm sw and lw routines
  !! 
  !! @par Revision History
  !! Original Source Rewritten and renamed by B. Stevens (2009-08)
  !!
  !! @remarks
  !!   Because the RRTM indexes vertical levels differently than ECHAM a chief
  !!   function of thise routine is to reorder the input in the vertical.  In 
  !!   addition some cloud physical properties are prescribed, which are 
  !!   required to derive cloud optical properties
  !!
  !! @par The gases are passed into RRTM via two multi-constituent arrays: 
  !!   zwkl and wx_r. zwkl has JPINPX species and  wx_r has JPXSEC species  
  !!   The species are identifed as follows.  
  !!     ZWKL [#/cm2]          WX_R [#/cm2]
  !!    index = 1 => H20     index = 1 => n/a
  !!    index = 2 => CO2     index = 2 => CFC11
  !!    index = 3 =>  O3     index = 3 => CFC12
  !!    index = 4 => N2O     index = 4 => n/a
  !!    index = 5 => n/a
  !!    index = 6 => CH4
  !!    index = 7 => O2
  !


  SUBROUTINE rrtm_interface( &
       & iaero           ,kproma          ,kbdim           ,klev            ,&
       & krow            ,ktrac           ,ktype           ,nb_sw           ,&
       & laland          ,laglac          ,cemiss          ,diff            ,&
       & pmu0            ,pgeom1          ,alb_vis         ,alb_nir         ,&
       & alb_vis_dir     ,alb_nir_dir     ,alb_vis_dif     ,alb_nir_dif     ,&
       & pp_fl           ,pp_hl           ,pp_sfc          ,tk_fl           ,&
       & tk_hl           ,tk_sfc          ,xm_vap          ,xm_sat          ,&
       & xm_liq          ,xm_ice          ,cdnc            ,icnc            ,&
       & cld_frc         ,&
       & cld_cvr         ,xm_o3           ,xm_co2          ,xm_ch4          ,&
       & xm_n2o          ,xm_cfc          ,xm_o2           ,pxtm1           ,&
       & flx_lw_net      ,flx_sw_net      ,flx_lw_net_clr  ,flx_sw_net_clr  ,&
       & flx_uplw_sfc    ,flx_upsw_sfc    ,flx_uplw_sfc_clr,flx_upsw_sfc_clr,&
       & emiss_sfc       ,sw_irr_toa      ,nir_sfc         ,nir_dff_sfc     ,&
       & vis_sfc         ,vis_dff_sfc     ,dvis_sfc        ,dpar_sfc        ,&
       & par_dff_sfc                                                         )

    INTEGER,INTENT(IN)  ::             &
         iaero,                        & !< aerosol control
         kproma,                       & !< number of longitudes
         kbdim,                        & !< first dimension of 2-d arrays
         krow,                         & !< first dimension of 2-d arrays
         klev,                         & !< number of levels
         ktrac,                        & !< number of tracers
         ktype(kbdim),                 & !< type of convection
         nb_sw                           !< number of shortwave bands

    LOGICAL,INTENT(IN) ::                  &
         laland(kbdim),                    & !< land sea mask, land=.true.
         laglac(kbdim)                       !< glacier mask, glacier=.true.

    REAL(WP),INTENT(IN)  ::            &
         cemiss,                       & !< surface emissivity
         diff,                         & !< solar diffusivity
         pmu0(kbdim),                  & !< mu0 for solar zenith angle
         pgeom1(kbdim,klev),           & !< geopotential above ground
         alb_vis(kbdim),               & !< surface albedo for visible range
                                !<  (direct and diffuse combined)
         alb_nir(kbdim),               & !< surface albedo for NIR range
                                !<  (direct and diffuse combined)
         alb_vis_dir(kbdim),           & !< surface albedo for vis range and dir light
         alb_nir_dir(kbdim),           & !< surface albedo for NIR range and dir light
         alb_vis_dif(kbdim),           & !< surface albedo for vis range and dif light
         alb_nir_dif(kbdim),           & !< surface albedo for NIR range and dif light
         pp_fl(kbdim,klev),            & !< full level pressure in Pa
         pp_hl(kbdim,klev+1),          & !< half level pressure in Pa
         pp_sfc(kbdim),                & !< surface pressure in Pa
         tk_fl(kbdim,klev),            & !< full level temperature in K
         tk_hl(kbdim,klev+1),          & !< half level temperature in K
         tk_sfc(kbdim),                & !< surface temperature in K
         xm_vap(kbdim,klev),           & !< specific humidity in g/g
         xm_sat(kbdim,klev),           & !< satur. specific humidity
         xm_liq(kbdim,klev),           & !< specific liquid water content
         xm_ice(kbdim,klev),           & !< specific ice content in g/g
         cdnc(kbdim,klev),             & !< cloud nuclei concentration
         icnc(kbdim,klev),             & !< ice crystal number concentration
         cld_frc(kbdim,klev),          & !< fractional cloud cover
         cld_cvr(kbdim),               & !< total cloud cover in m2/m2
         xm_o3(kbdim,klev),            & !< o3 mass mixing ratio
         xm_co2(kbdim,klev),           & !< co2 mass mixing ratio
         xm_ch4(kbdim,klev),           & !< ch4 mass mixing ratio
         xm_n2o(kbdim,klev),           & !< n2o mass mixing ratio
         xm_cfc(kbdim,klev,2),         & !< cfc volume mixing ratio
         xm_o2(kbdim,klev),            & !< o2  mass mixing ratio
         pxtm1(kbdim,klev,ktrac)         !< tracer mass mixing ratios

    REAL (wp), INTENT (OUT) ::            &
         flx_lw_net(kbdim,klev+1),        & !< net downward LW flux profile,
         flx_sw_net(kbdim,klev+1),        & !< net downward SW flux profile,
         flx_lw_net_clr(kbdim,klev+1),    & !< clrsky downward LW flux profile,
         flx_sw_net_clr(kbdim,klev+1),    & !< clrsky downward SW flux profile,
         flx_uplw_sfc(kbdim),             & !< sfc LW upward flux,
         flx_upsw_sfc(kbdim),             & !< sfc SW upward flux,
         flx_uplw_sfc_clr(kbdim),         & !< clrsky sfc LW upward flux,
         flx_upsw_sfc_clr(kbdim),         & !< clrsky sfc SW upward flux,
         emiss_sfc(kbdim),                & !< effective surface emissivity
         sw_irr_toa(kbdim),               & !< top of atmosphere solar irradiation
         nir_sfc(kbdim),                  & !< net surface NIR flux
         nir_dff_sfc(kbdim),              & !< fraction of diffuse NIR
         vis_sfc(kbdim),                  & !< net surface visible flux
         vis_dff_sfc(kbdim),              & !< fraction of diffuse visible
         dvis_sfc(kbdim),                 & !< surf. visible downward
         dpar_sfc(kbdim),                 & !< surf. PAR downw.
         par_dff_sfc(kbdim)                     !< fraction of diffuse PAR

    INTEGER  :: jk, jl, jp, jkb,          & !< loop indicies
         icldlyr(kbdim,klev)                !< index for clear or cloudy

    REAL(wp) ::                           &
         zsemiss(kbdim,jpband),           & !< LW surface emissivity by band
         clr_frc(kbdim),                  & !< clear sky fraction of total column
         ppd_hl(kbdim,klev),              & !< pressure thickness in Pa
         pm_sfc(kbdim),                   & !< surface pressure in mb
         dnir_sfc(kbdim),                 & !< surf. NIR downw.
         unir_toa(kbdim),                 & !< top NIR upw.
         uvis_toa(kbdim),                 & !< top vis. upw.
         dnir_sfc_cld(kbdim),             & !< surf. NIR downw. cloudy sky
         dvis_sfc_cld(kbdim),             & !< surf. vis. downw. cloudy sky
         unir_toa_cld(kbdim),             & !< top NIR upw. cloudy sky
         uvis_toa_cld(kbdim),             & !< top vis. upw. cloudy sky
         zsudu(kbdim),                    & !< direct beam SW
         zfuvf(kbdim),                    & !< dummy output array from sw
         amm,                             & !< molecular weight of moist air
         delta,                           & !< pressure thickness 
         zscratch                           !< scratch array
    !
    ! --- vertically reversed _vr variables
    !
    REAL(wp) ::                           &
         xcld_frc_vr(kbdim,klev),         & !< lyr cloud frac wrt clmn cld frac
         col_dry_vr(kbdim,klev),          & !< number of molecules/cm2 of
         pm_fl_vr(kbdim,klev),            & !< full level pressure [mb] 
         pm_hl_vr(kbdim,klev+1),          & !< half level pressure [mb] 
         tk_fl_vr(kbdim,klev),            & !< full level temperature [K]
         tk_hl_vr(kbdim,klev+1),          & !< half level temperature [K]
         cdnc_vr(kbdim,klev),             & !< cloud nuclei concentration
         icnc_vr(kbdim,klev),             & !< ice crystal number concentration
         cld_frc_vr(kbdim,klev),          & !< secure cloud fraction
         ziwgkg_vr(kbdim,klev),           & !< specific ice water content
         ziwc_vr(kbdim,klev),             & !< ice water content per volume
         ziwp_vr(kbdim,klev),             & !< ice water path in g/m2  
         zlwgkg_vr(kbdim,klev),           & !< specific liquid water content
         zlwp_vr(kbdim,klev),             & !< liquid water path in g/m2  
         zlwc_vr(kbdim,klev),             & !< liquid water content per
         zoz_vr(kbdim,klev),              & !< ozone mass mixing ratio *
         wkl_vr(kbdim,jpinpx,klev),       & !< number of molecules/cm2 of
         wx_vr(kbdim,jpxsec,klev),        & !< number of molecules/cm2 of
         cld_tau_lw_vr(kbdim,klev,jpband),& !< LW optical thickness of clouds
         cld_tau_sw_vr(kbdim,nb_sw,klev), & !< extincion
         cld_cg_sw_vr(kbdim,nb_sw,klev),  & !< asymmetry factor
         cld_piz_sw_vr(kbdim,nb_sw,klev), & !< single scattering albedo
         aer_tau_lw_vr(kbdim,klev,jpband),& !< LW optical thickness of aerosols
         aer_tau_sw_vr(kbdim,klev,nb_sw), & !< aerosol optical thickness
         aer_cg_sw_vr(kbdim,klev,nb_sw),  & !< aerosol asymmetry factor
         aer_piz_sw_vr(kbdim,klev,nb_sw), & !< aerosol single scattering albedo
         flx_uplw_vr(kbdim,klev+1),       & !< upward flux, total sky
         flx_uplw_clr_vr(kbdim,klev+1),   & !< upward flux, clear sky
         flx_dnlw_vr(kbdim,klev+1),       & !< downward flux, total sky
         flx_dnlw_clr_vr(kbdim,klev+1),   & !< downward flux, clear sky
         flx_upsw_vr(kbdim,klev+1),       & !< upward flux total sky
         flx_upsw_clr_vr(kbdim,klev+1),   & !< upward flux clear sky
         flx_dnsw_vr(kbdim,klev+1),       & !< downward flux total sky
         flx_dnsw_clr_vr(kbdim,klev+1),   & !< downward flux clear sky
         flx_upsw(kbdim,klev+1),          & !< upward flux total sky
         flx_upsw_clr(kbdim,klev+1),      & !< upward flux clear sky
         flx_dnsw(kbdim,klev+1),          & !< downward flux total sky
         flx_dnsw_clr(kbdim,klev+1),      & !< downward flux clear sky
         htngrt_vr(kbdim,klev),           & !< heating rate total sky
         htngrt_clr_vr(kbdim,klev),       & !< heating rate clear sky
         htngrt(kbdim,klev),              & !< heating rate total sky
         htngrt_clr(kbdim,klev),          & !< heating rate clear sky
         re_drop(kbdim,klev),             & !< effective radius of liquid
         re_cryst(kbdim,klev)               !< effective radius of ice

    REAL (wp), PARAMETER :: fcf = 2.0_wp*api*10000.0_wp !< flux conversion fact
    !
    ! 1.0 Constituent properties 
    !--------------------------------
    DO jk = 1, klev
      jkb = klev+1-jk
      DO jl = 1, kproma
        cld_frc_vr(jl,jk) = MAX(EPSILON(1.0_wp),cld_frc(jl,jkb))
        ziwgkg_vr(jl,jk)  = xm_ice(jl,jkb)*1000.0_wp/cld_frc_vr(jl,jk)
        zlwgkg_vr(jl,jk)  = xm_liq(jl,jkb)*1000.0_wp/cld_frc_vr(jl,jk)
      END DO
    END DO
    !
    ! --- control for zero, infintesimal or negative cloud fractions
    !
    WHERE (cld_frc_vr(1:kproma,:) > 2.0_wp*EPSILON(1.0_wp))
      icldlyr(1:kproma,:) = 1
    ELSEWHERE
      icldlyr(1:kproma,:)  = 0
      ziwgkg_vr(1:kproma,:) = 0.0_wp
      zlwgkg_vr(1:kproma,:) = 0.0_wp
    END WHERE
    !
    ! --- main constituent reordering
    !
    DO jl = 1, kproma
      pm_hl_vr(jl,klev+1) = 0.01_wp*pp_hl(jl,1)
      tk_hl_vr(jl,klev+1) = tk_hl(jl,1)
      pm_sfc(jl)          = 0.01_wp*pp_sfc(jl)
    END DO

    DO jk = 1, klev
      jkb = klev+1-jk
      DO jl = 1, kproma
        delta = pp_hl(jl,jkb+1)-pp_hl(jl,jkb)
        !
        ! --- thermodynamic arrays
        !
        pm_hl_vr(jl,jk) = 0.01_wp*pp_hl(jl,jkb+1)
        pm_fl_vr(jl,jk) = 0.01_wp*pp_fl(jl,jkb)
        tk_hl_vr(jl,jk) = tk_hl(jl,jkb+1)
        tk_fl_vr(jl,jk) = tk_fl(jl,jkb)
        !
        ! --- cloud properties
        !
        zscratch      = pp_fl(jl,jkb)/tk_fl(jl,jkb)
        ziwc_vr(jl,jk) = ziwgkg_vr(jl,jk)*zscratch/rd
        ziwp_vr(jl,jk) = ziwgkg_vr(jl,jk)*delta/g
        zlwc_vr(jl,jk) = zlwgkg_vr(jl,jk)*zscratch/rd
        zlwp_vr(jl,jk) = zlwgkg_vr(jl,jk)*delta/g
        cdnc_vr(jl,jk) = cdnc(jl,jkb)*1.e-6_wp
        icnc_vr(jl,jk) = icnc(jl,jkb)*1.e-6_wp
        !
        ! --- radiatively active gases
        !
        zoz_vr(jl,jk)     = xm_o3(jl,jkb)*delta*46.6968_wp/g
        wkl_vr(jl,:,jk)   = 0.0_wp
        wkl_vr(jl,1,jk)   = xm_vap(jl,jkb)*amd/amw
        wkl_vr(jl,2,jk)   = xm_co2(jl,jkb)*amd/amco2
        wkl_vr(jl,3,jk)   = xm_o3(jl,jkb) *amd/amo3
        wkl_vr(jl,4,jk)   = xm_n2o(jl,jkb)*amd/amn2o
        wkl_vr(jl,6,jk)   = xm_ch4(jl,jkb)*amd/amch4
        wkl_vr(jl,7,jk)   = xm_o2(jl,jkb)*amd/amo2
        amm               = (1.0_wp-wkl_vr(jl,1,jk))*amd + wkl_vr(jl,1,jk)*amw
        col_dry_vr(jl,jk) = (0.01_wp*delta)*10.0_wp*avo/g/amm                &
             &              / (1.0_wp+wkl_vr(jl,1,jk))
        !
        ! --- alternate treatment for cfcs
        !
        wx_vr(jl,:,jk) = 0.0_wp
        wx_vr(jl,2,jk) = col_dry_vr(jl,jk)*xm_cfc(jl,jkb,1)*1.e-20_wp
        wx_vr(jl,3,jk) = col_dry_vr(jl,jk)*xm_cfc(jl,jkb,2)*1.e-20_wp
      END DO
    END DO
    DO jp = 1, 7
      wkl_vr(1:kproma,jp,:)=col_dry_vr(1:kproma,:)*wkl_vr(1:kproma,jp,:)
    END DO
    !
    ! 2.0 Surface Properties
    ! --------------------------------
    zsemiss(1:kproma,:) = cemiss 
    !
    ! 3.0 Particulate Optical Properties
    ! --------------------------------
    ppd_hl(1:kproma,:) = pp_hl(1:kproma,2:klev+1)-pp_hl(1:kproma,1:klev)

    SELECT CASE (iaero)
    CASE (0)
      aer_tau_lw_vr(:,:,:) = 0.0_wp
      aer_tau_sw_vr(:,:,:) = 0.0_wp
      aer_piz_sw_vr(:,:,:) = 0.0_wp
      aer_cg_sw_vr(:,:,:)  = 0.0_wp
!LT    CASE (2)
!      CALL aerosol_optics_tanre( &
!           & krow            ,kproma          ,kbdim           ,klev           ,&
!           & jpband          ,nb_sw           ,cetah           ,ppd_hl         ,&
!           & pp_fl           ,tk_hl           ,aer_tau_lw_vr   ,aer_tau_sw_vr  ,&
!           & aer_cg_sw_vr    ,aer_piz_sw_vr   )
!LT
    CASE (3)
      CALL set_aop_kinne( &
           & kproma           ,kbdim                 ,klev             ,&
           & krow             ,jpband                ,nb_sw            ,&
           & aer_tau_lw_vr    ,aer_tau_sw_vr         ,aer_piz_sw_vr    ,&
           & aer_cg_sw_vr     ,ppd_hl                ,pp_fl            ,&
           & tk_fl            ,pgeom1 )
!    CASE (5)
!      CALL set_aop_kinne( &
!           & kproma           ,kbdim                 ,klev             ,&
!           & krow             ,jpband                ,nb_sw            ,&
!           & aer_tau_lw_vr    ,aer_tau_sw_vr         ,aer_piz_sw_vr    ,&
!           & aer_cg_sw_vr     ,ppd_hl                ,pp_fl            ,&
!           & tk_fl            ,pgeom1 )
!      CALL add_aop_volc( &
!           & kproma           ,kbdim                 ,klev             ,&
!           & krow             ,jpband                ,nb_sw            ,&
!           & aer_tau_lw_vr    ,aer_tau_sw_vr         ,aer_piz_sw_vr    ,&
!           & aer_cg_sw_vr     ,ppd_hl                ,pp_fl            ,&
!           & tk_fl )

    CASE DEFAULT
    END SELECT

    IF (l_newoptics .AND. l_srtm) THEN
      CALL newcld_optics(                                                       &
           & laglac          ,laland          ,kproma          ,kbdim          ,&
           & klev            ,ktype           ,jpband          ,nb_sw          ,&
           & icldlyr         ,zlwp_vr         ,ziwp_vr         ,zlwc_vr        ,&
           & ziwc_vr         ,cdnc_vr         ,icnc_vr         ,cld_tau_lw_vr  ,&
           & cld_tau_sw_vr   ,cld_piz_sw_vr   ,cld_cg_sw_vr    ,re_drop        ,&
           & re_cryst       )
    ELSE       
      CALL cloud_optics(                                                        &
           & laglac          ,laland          ,kproma          ,kbdim          ,&
           & klev            ,ktype           ,jpband          ,nb_sw          ,&
           & diff            ,icldlyr         ,zlwp_vr         ,ziwp_vr        ,&
           & zlwc_vr         ,ziwc_vr         ,cdnc_vr         ,cld_tau_lw_vr  ,&
           & cld_tau_sw_vr   ,cld_piz_sw_vr   ,cld_cg_sw_vr    ,re_cryst       ,&
           & re_drop                                                           )
    END IF

    clr_frc(1:kproma) = 1.0_wp-cld_cvr(1:kproma)
    DO jl = 1, kproma
      IF (cld_cvr(jl) < EPSILON(1.0_wp)) THEN
        xcld_frc_vr(jl,:) = 0.0_wp
      ELSE
        xcld_frc_vr(jl,:) = cld_frc_vr(jl,:)/cld_cvr(jl)
      END IF
    END DO

!cms++
!  IF (locosp) THEN
!    DO jk=1,klev
!      jkb = klev+1-jk
!     DO jl=1,kproma
!        cosp_reffl(jl,jk,krow) = re_drop(jl,jkb)
!        cosp_reffi(jl,jk,krow) = re_cryst(jl,jkb)
!        cosp_f3d(jl,jk,krow) = cld_frc(jl,jk)
!     END DO
!    END DO
!     IF ( Lisccp_sim ) THEN   
        !don't really need vert. re-arrange ...
!      IF (l_newoptics .and. l_srtm) THEN
!        DO jk=1,klev
!          jkb = klev+1-jk
!          DO jl=1,kproma
!            cisccp_cldtau3d(jl,jk,krow) = &
!               cld_tau_sw_vr(jl, 9  ,jkb)  ! band 9 is 625 - 778 nm, needed is 670 nm
!          END DO
!        END DO
!      ELSE
!        DO jk=1,klev
!          jkb = klev+1-jk
!          DO jl=1,kproma
!            cisccp_cldtau3d(jl,jk,krow) = &
!               cld_tau_sw_vr(jl,  3 ,jkb)  ! band 3 is 440 - 690 nm, needed is 670 nm
!          END DO
!        END DO
!      END IF
!        DO jk=1,klev
!          jkb = klev+1-jk
!          DO jl=1,kproma
           !this is ONLY o.k. as long as wp equals dp, else conversion needed
!             cisccp_cldemi3d(jl,jk,krow) = &
!               1._wp - exp(-1._wp*cld_tau_lw_vr(jl,jkb,6)) ! band 6 is 820 - 980 cm-1, 
!                                                                   !  we need 10.5 µm channel
!          END DO
!        END DO
!     END IF
!  END IF
!cms--

    !
    ! 3.5 Interface for submodels that provide aerosol and/or cloud radiative properties:
    ! -----------------------------------------------------------------------------------
!    IF (lanysubmodel) THEN
!      CALL radiation_subm_1(kproma           ,kbdim            ,klev         ,krow  ,&
!                            ktrac            ,iaero            ,jpband       ,nb_sw ,&
!                            aer_tau_sw_vr    ,aer_piz_sw_vr    ,aer_cg_sw_vr        ,&
!                            aer_tau_lw_vr                                           ,&
!                            ppd_hl           ,pxtm1                                  )
!    END IF
    !
    ! 4.0 Radiative Transfer Routines
    ! --------------------------------
    IF (l_lrtm) THEN
      CALL lrtm(                                                                 &
           & kproma, kbdim, klev              ,pm_fl_vr        ,pm_sfc          ,&
           & tk_fl_vr        ,tk_hl_vr        ,tk_sfc          ,wkl_vr          ,&
           & wx_vr           ,col_dry_vr      ,zsemiss         ,cld_frc_vr      ,&
           & cld_tau_lw_vr   ,aer_tau_lw_vr   ,flx_uplw_vr     ,flx_dnlw_vr     ,&
           & flx_uplw_clr_vr, flx_dnlw_clr_vr )
    ELSE  
      CALL rrtm_rrtm_140gp(                                                      &
           & kproma          ,kbdim           ,klev            ,pm_fl_vr        ,& 
           & tk_fl_vr        ,tk_hl_vr        ,tk_sfc          ,zsemiss         ,&
           & cld_frc_vr      ,icldlyr         ,col_dry_vr      ,wkl_vr          ,&
           & wx_vr           ,cld_tau_lw_vr   ,aer_tau_lw_vr   ,flx_uplw_vr     ,&
           & flx_uplw_clr_vr ,flx_dnlw_vr     ,flx_dnlw_clr_vr ,emiss_sfc       )
      flx_uplw_clr_vr(1:kproma,:) = fcf * flx_uplw_clr_vr(1:kproma,:)
      flx_dnlw_clr_vr(1:kproma,:) = fcf * flx_dnlw_clr_vr(1:kproma,:)
      flx_uplw_vr(1:kproma,:) = fcf * flx_uplw_vr(1:kproma,:)
      flx_dnlw_vr(1:kproma,:) = fcf * flx_dnlw_vr(1:kproma,:)
    END IF
    IF (nb_sw == 14) THEN
      CALL srtm_srtm_224gp(                                                     &
           & kproma          ,kbdim           ,klev            ,nb_sw          ,&
           & alb_vis_dir     ,alb_nir_dir     ,alb_vis_dif     ,alb_nir_dif    ,&
           & pm_fl_vr        ,tk_fl_vr                                         ,&
           & ppd_hl          ,pmu0            ,col_dry_vr      ,wkl_vr         ,&
           & cld_frc_vr      ,cld_tau_sw_vr   ,cld_cg_sw_vr    ,cld_piz_sw_vr  ,&
           & aer_tau_sw_vr   ,aer_cg_sw_vr    ,aer_piz_sw_vr   ,psctm          ,&
           & ssi_factor      ,zfuvf           ,htngrt          ,flx_dnsw       ,&
           & flx_upsw        ,htngrt_clr      ,flx_dnsw_clr    ,flx_upsw_clr   ,&
           & dnir_sfc        ,dvis_sfc        ,unir_toa        ,uvis_toa       ,&
           & dnir_sfc_cld    ,dvis_sfc_cld    ,unir_toa_cld    ,uvis_toa_cld   ,&
           & zsudu           ,nir_sfc         ,nir_dff_sfc     ,vis_sfc        ,&
           & vis_dff_sfc     ,dpar_sfc        ,par_dff_sfc)
    ELSE
      CALL sw(                                                                  &
           & 1               ,kproma          ,kbdim           ,klev           ,&
           & psctm           ,pmu0            ,pp_sfc          ,pm_hl_vr       ,&
           & tk_fl_vr        ,alb_vis         ,alb_nir         ,xm_co2         ,&
           & xm_vap          ,xm_sat          ,clr_frc         ,xcld_frc_vr    ,&
           & ppd_hl          ,zoz_vr          ,cld_cg_sw_vr    ,cld_piz_sw_vr  ,&
           & cld_tau_sw_vr   ,aer_tau_sw_vr   ,aer_piz_sw_vr   ,aer_cg_sw_vr   ,&
           & htngrt_vr       ,flx_dnsw_vr     ,flx_upsw_vr     ,htngrt_clr_vr  ,&
           & flx_dnsw_clr_vr ,flx_upsw_clr_vr ,dnir_sfc        ,dvis_sfc       ,&
           & unir_toa        ,uvis_toa        ,dnir_sfc_cld    ,dvis_sfc_cld   ,&
           & unir_toa_cld    ,uvis_toa_cld    ,zsudu           ,nir_sfc        ,&
           & nir_dff_sfc     ,vis_sfc         ,vis_dff_sfc                     )
    END IF
    !
    ! 5.0 Post Processing
    ! --------------------------------
    IF (nb_sw == 14) THEN
      DO jk = 1, klev+1
        jkb = klev+2-jk
        DO jl = 1, kproma
          flx_lw_net(jl,jk)     = flx_dnlw_vr(jl,jkb)-flx_uplw_vr(jl,jkb) 
          flx_lw_net_clr(jl,jk) = flx_dnlw_clr_vr(jl,jkb)-flx_uplw_clr_vr(jl,jkb) 
          flx_sw_net(jl,jk)     = flx_dnsw(jl,jk) - flx_upsw(jl,jk)
          flx_sw_net_clr(jl,jk) = flx_dnsw_clr(jl,jk)-flx_upsw_clr(jl,jk)
        END DO
      END DO
      flx_uplw_sfc(1:kproma)     = flx_uplw_vr(1:kproma,1)   
      flx_uplw_sfc_clr(1:kproma) = flx_uplw_clr_vr(1:kproma,1)   
      flx_upsw_sfc(1:kproma)     = flx_upsw(1:kproma,klev+1)
      flx_upsw_sfc_clr(1:kproma) = flx_upsw_clr(1:kproma,klev+1)
      sw_irr_toa(1:kproma)       = flx_dnsw(1:kproma,1) 
    ELSE
      DO jk = 1, klev+1
        jkb = klev+2-jk
        DO jl=1,kproma
          flx_lw_net(jl,jk)     = flx_dnlw_vr(jl,jkb)-flx_uplw_vr(jl,jkb) 
          flx_lw_net_clr(jl,jk) = flx_dnlw_clr_vr(jl,jkb)-flx_uplw_clr_vr(jl,jkb) 
          flx_sw_net(jl,jk)     = flx_dnsw_vr(jl,jkb) - flx_upsw_vr(jl,jkb)
          flx_sw_net_clr(jl,jk) = flx_dnsw_clr_vr(jl,jkb)-flx_upsw_clr_vr(jl,jkb)
        END DO
      END DO
      flx_uplw_sfc(1:kproma)     = flx_uplw_vr(1:kproma,1)   
      flx_uplw_sfc_clr(1:kproma) = flx_uplw_clr_vr(1:kproma,1)   
      flx_upsw_sfc(1:kproma)     = flx_upsw_vr(1:kproma,1)
      flx_upsw_sfc_clr(1:kproma) = flx_upsw_clr_vr(1:kproma,1)
      sw_irr_toa(1:kproma)       = flx_dnsw_vr(1:kproma,klev+1) 
    END IF
    !
    ! 6.0 Interface for submodel diagnosics after radiation calculation:
    ! ------------------------------------------------------------------
!LT
!    IF (lanysubmodel) THEN
!      CALL radiation_subm_2(kproma, kbdim, krow, klev, &
!                            ktrac,  iaero,             &
!                            pxtm1                      )
!    END IF

  END SUBROUTINE rrtm_interface

  SUBROUTINE setup_radiation_new(ih2o_in,ico2_in,ich4_in,io3_in, &
       & io2_in,in2o_in,iaero_in)

    USE mo_aero_kinne,

    INTEGER, INTENT(IN) :: ih2o_in
    INTEGER, INTENT(IN) :: ico2_in
    INTEGER, INTENT(IN) :: ich4_in
    INTEGER, INTENT(IN) :: io3_in
    INTEGER, INTENT(IN) :: io2_in
    INTEGER, INTENT(IN) :: in2o_in
    INTEGER, INTENT(IN) :: iaero_in


    ! This replaces reading the namelist
    ih2o=ih2o_in
    ico2=ico2_in
    ich4=ich4_in
    io3=io3_in
    io2=io2_in
    in2o=in2o_in
    iaero=iaero_in


    IF (lrad) THEN
       IF (l_srtm) THEN
          nb_sw = jpsw
          CALL setup_srtm
       ELSE
          nb_sw = nsw
          CALL su_sw4
       END IF
       IF (l_lrtm) THEN
          CALL lrtm_setup
       ELSE
          CALL su_rrtm
       END IF

       SELECT CASE (ico2)
       CASE(0)
          co2mmr = co2vmr*amco2/amd  
       CASE(2)
          co2mmr = co2vmr*amco2/amd
       CASE default
       END SELECT


       SELECT CASE (ich4)
       CASE(2)
          ch4mmr = ch4vmr*amch4/amd
       CASE default
       END SELECT

       SELECT CASE (ich4)
       CASE(2)
          ch4mmr = ch4vmr*amch4/amd
       CASE default
       END SELECT

       SELECT CASE (in2o)
       CASE(2)
          n2ommr = n2ovmr*amn2o/amd
       CASE default
       END SELECT

       SELECT CASE (io2)
       CASE(2)
          o2mmr = o2vmr*amo2/amd
       CASE default
       END SELECT

       SELECT CASE (iaero)
       CASE(0)
          ! CALL message('','iaero= 0 --> no aerosol in radiation')
       CASE(3)
          ! CALL message('','iaero= 3 --> Kinne climatology')
          CALL su_aero_kinne(nb_sw)
          CALL read_aero_kinne(nb_sw)
       CASE default
       END SELECT

       IF (l_newoptics .AND. l_srtm) THEN
          CALL setup_newcld_optics
       ELSE 
          CALL setup_cloud_optics(lcouple, lipcc, nn, cetah)
       END IF

    END IF


  END SUBROUTINE setup_radiation_new
  
  

  !---------------------------------------------------------------------------
  !>
  !! @brief  Sets up (initializes) radation routines
  !! 
  !! @remarks
  !!   Modify preset variables of module MO_RADIATION which control the 
  !!   configuration of the radiation scheme.
  !
  SUBROUTINE setup_radiation

    USE mo_aero_kinne,       ONLY: su_aero_kinne, read_aero_kinne
!    USE mo_aero_volc,        ONLY: su_aero_volc
!    USE mo_solar_irradiance, ONLY: su_solar_irradiance

    INTEGER :: ierr
    !
    ! 1.0 Read radctl namelist to modify mo_radiation
    ! --------------------------------
!LT
!  READ (nin, radctl)  	
!LT
!    IF (p_parallel_io) THEN
!      !
!      ! --- In case of coupled runs initially set basic values to 1850 values
!      ! 
!      IF(lcouple .OR. lipcc) THEN
!        co2vmr    = 284.7e-06_wp !< 1850 concentration
!        ch4vmr    = 791.0e-09_wp !< 1850 concentration
!        n2ovmr    = 275.4e-09_wp !< 1850 concentration
!        cfcvmr(1) = 0.0_wp
!        cfcvmr(2) = 0.0_wp
!      ENDIF
!      !
!      ! --- Change default behavior for non Middle Atmosphere runs
!      ! 
!      IF(.NOT. lmidatm) THEN
!        ich4 = 2
!        in2o = 2
!      ENDIF
!      !
!      ! --- Read NAMELIST
!      ! 
!      CALL position_nml ('RADCTL', status=ierr)
!      SELECT CASE (ierr)
!      CASE (POSITIONED)
!        READ (nnml, radctl)
!      END SELECT
!      IF(ich4 == 4 .OR. in2o == 4 .OR. ico2 == 4 .OR. icfc == 4 ) ighg = 1
!    ENDIF
!
!    !
!    ! 2.0 Broadcast NAMELIST variables
!    ! --------------------------------
!    IF (p_parallel) THEN
!      CALL p_bcast (nmonth, p_io)
!      CALL p_bcast (isolrad, p_io)
!      CALL p_bcast (ldiur, p_io)
!      CALL p_bcast (lradforcing, p_io)
!      CALL p_bcast_event (trigrad, p_io)
!      CALL p_bcast (ih2o, p_io)
!      CALL p_bcast (ico2, p_io)
!      CALL p_bcast (ich4, p_io)
!      CALL p_bcast (io3, p_io)
!      CALL p_bcast (io2, p_io)
!      CALL p_bcast (in2o, p_io)
!      CALL p_bcast (icfc, p_io)
!      CALL p_bcast (ighg, p_io)
!      CALL p_bcast (iaero, p_io)
!      CALL p_bcast (co2vmr, p_io)
!      CALL p_bcast (ch4vmr, p_io)
!      CALL p_bcast (n2ovmr, p_io)
!      CALL p_bcast (cfcvmr, p_io)
!      CALL p_bcast (o2vmr, p_io)
!      CALL p_bcast (cecc, p_io)
!      CALL p_bcast (cobld, p_io)
!      CALL p_bcast (clonp, p_io)
!      CALL p_bcast (yr_perp, p_io)
!      CALL p_bcast (l_srtm, p_io)
!      CALL p_bcast (l_lrtm, p_io)
!      CALL p_bcast (l_newoptics, p_io)
!    ENDIF
!    !
!    IF (nmonth > 0 .AND. get_calendar_type() /= JULIAN) THEN
!      CALL finish('setup_radiation', &
!           &      ' ly360=.TRUE. cannot run perpetual month setup (nmonth > 0).')
!    ENDIF
    !
    ! 3.0 If radiation is active check NAMELIST variable conformance
    ! --------------------------------
    IF (lrad) THEN

      IF (l_srtm) THEN
        nb_sw = jpsw
        CALL setup_srtm
      ELSE
        nb_sw = nsw
        CALL su_sw4
      END IF
      IF (l_lrtm) THEN
        CALL lrtm_setup
      ELSE
        CALL su_rrtm
      END IF
      !
!      CALL message('','lrad = .TRUE.  --> ECMWF Cy21R4 radiation')
      !
      ! --- Check  H2O
      !
!      SELECT CASE (ih2o)
!      CASE(0)
!        CALL message('','ih2o = 0 --> no H2O(gas,liquid,ice) in radiation')
!      CASE(1)
!        CALL message('','ih2o = 1 --> prognostic H2O(gas,liquid,ice)')
!      CASE default
!        WRITE (message_text, '(a,i2,a)') &
!             'ih2o =', ih2o, ' in radctl namelist is not supported'
!        CALL message('', message_text)
!        CALL finish('setup_radiation','Run terminated ih2o')
!      END SELECT
      !
      ! --- Check  CO2
      ! 
       SELECT CASE (ico2)
       CASE(0)
!        CALL message('','ico2 = 0 --> no CO2 in radiation')
         co2mmr=co2vmr*amco2/amd   ! Necessary for use with lco2=.TRUE.
!      CASE(1)  
!        IF (lco2) THEN
!          WRITE (message_text, '(a,f12.8)') &
!               'ico2 = 1 --> Initial CO2 volume mixing ratio=', co2vmr
!          CALL message('',message_text)
!          co2mmr = co2vmr*amco2/amd
!        ELSE
!          CALL finish('setup_radiation','ico2=1 (interactive CO2) not '// &
!               &      'a valid choice for lco2=.false.')
!        END IF
       CASE(2)
!        WRITE (message_text, '(a,f12.8)') &
!             'ico2 = 2 --> CO2 volume mixing ratio=', co2vmr
!        CALL message('',message_text)
         co2mmr = co2vmr*amco2/amd
!      CASE(4)
!        CALL message('','ico2 = 4 --> CO2 volume mixing ratio from scenario')
!        co2mmr = co2vmr*amco2/amd    ! This is only a dummy value for the first
                                     ! initialization of co2m1 in the co2-module. 
                                     ! co2m1 will be overwritten with the correct
                                     ! values as soon as the ghg-data are
                                     ! interpolated to the right date.
       CASE default
!        WRITE (message_text, '(a,i2,a)') &
!             'ico2 = ', ico2, ' in radctl namelist is not supported'
!        CALL message('',message_text)
!        CALL finish('setup_radiation','Run terminated ico2')
       END SELECT
      !
      ! --- Check CH4
      ! 
       SELECT CASE (ich4)
!      CASE(0)
!        CALL message('','ich4 = 0 --> no CH4 in radiation')
!      CASE(1)
!        CALL message('','ich4 = 1 --> transported CH4 is not yet implemented')
!        CALL finish('setup_radiation','Run terminated ich4')
       CASE(2)
!        WRITE (message_text, '(a,f12.8)') &
!             'ich4 = 2 --> CH4 volume mixing ratio=', ch4vmr
!        CALL message('',message_text)
         ch4mmr = ch4vmr*amch4/amd
!      CASE(3)
!        WRITE (message_text, '(a,f12.8)') &
!             'ich4 = 3 --> CH4 (trop) volume mixing ratio =', ch4vmr
!        CALL message('',message_text)
!        ch4mmr = ch4vmr*amch4/amd
!      CASE(4)
!        CALL message('','ich4 = 4 --> CH4 volume mixing ratio from scenario')
       CASE default
!        WRITE (message_text, '(a,i2,a)') &
!             'ich4 =', ich4, ' in radctl namelist is not supported'
!        CALL message('',message_text)
!        CALL finish('setup_radiation','Run terminated ich4')
       END SELECT
      !
      ! --- Check O3
      ! 
!      SELECT CASE (io3)
!      CASE(0)
!        CALL message('','io3  = 0 --> no O3 in radiation')
!      CASE(1)
!        CALL message('','io3  = 1 --> transported O3 is not yet implemented')
!        CALL finish('setup_radiation','Run terminated io3')
!      CASE(2)
!        CALL message('','io3  = 2 --> spectral O3 climatology (ECHAM4)')
!      CASE(3)
!        CALL message('','io3  = 3 --> gridpoint O3 climatology from NetCDF file')
!      CASE(4)
!        CALL message('','io3  = 4 --> gridpoint O3 climatology from IPCC-NetCDF file')
!      CASE default
!        WRITE (message_text, '(a,i2,a)') &
!             'io3  =', io3, ' in radctl namelist is not supported'
!        CALL message('',message_text)
!        CALL finish('setup_radiation','Run terminated io3')
!      END SELECT
      !
      ! --- Check N2O
      ! 
       SELECT CASE (in2o)
!      CASE(0)
!        CALL message('','in2o = 0 --> no N2O in radiation')
!      CASE(1)
!        CALL message('','in2o = 1 --> transported N2O is not yet implemented')
!        CALL finish('setup_radiation','Run terminated in2o')
       CASE(2)
!        WRITE (message_text, '(a,f12.8)') &
!             'in2o = 2 --> N2O volume mixing ratio=', n2ovmr
!        CALL message('',message_text)
         n2ommr = n2ovmr*amn2o/amd
!      CASE(3)
!        WRITE (message_text, '(a,f12.8)') &
!             'in2o = 3 --> N2O (trop) volume mixing ratio=', n2ovmr
!        CALL message('',message_text)
!        n2ommr = n2ovmr*amn2o/amd
!      CASE(4)
!        CALL message('','in2o = 4 --> N2O volume mixing ratio from scenario')
       CASE default
!        WRITE (message_text, '(a,i2,a)') &
!             'in2o =',in2o,' in radctl namelist is not supported'
!        CALL message('',message_text)
!        CALL finish('setup_radiation','Run terminated in2o')
       END SELECT
      !
      ! --- Check CFCs
      ! 
!      SELECT CASE (icfc)
!      CASE(0)
!        CALL message('','icfc = 0 --> no CFCs in radiation')
!      CASE(1)
!        CALL message('','icfc = 1 --> transported CFCs not yet implemented')
!        CALL finish('setup_radiation','Run terminated icfc')
!      CASE(2)
!        WRITE (message_text, '(a,f12.8)') &
!             'icfc = 2 --> CFC11    volume mixing ratio=', cfcvmr(1)
!        CALL message('',message_text)
!        WRITE (message_text, '(a,f12.8)') &
!             '             CFC12    volume mixing ratio=', cfcvmr(2)
!        CALL message('',message_text)
!      CASE(4)
!        CALL message('','icfc = 4 --> CFC volume mixing ratio from scenario')
!      CASE default
!        WRITE (message_text, '(a,i2,a)') &
!             'icfc=', icfc, ' in radctl namelist is not supported'
!        CALL message('',message_text)
!        CALL finish('setup_radiation','Run terminated icfc')
!      END SELECT
      !
      ! --- Check Scenario
      ! 
!      SELECT CASE (ighg)
!      CASE(0)
!        CALL message('','ighg = 0 --> no scenario, fixed greenhouse gases and/or cfc')
!      CASE(1)
!        CALL message('','ighg = 1 --> greenhouse gases from scenario, check setting of switches')
!      END SELECT
      !
      ! --- Check O2
      ! 
       SELECT CASE (io2)
!      CASE(0)
!        CALL message('','io2  = 0 --> no O2  in radiation')
       CASE(2)
!        WRITE (message_text, '(a,f12.8)') &
!             'io2  = 2 --> O2    volume mixing ratio=', o2vmr
!        CALL message('',message_text)
         o2mmr = o2vmr*amo2/amd
       CASE default
!        WRITE (message_text, '(a,i2,a)') &
!             'io2 =', io2, ' in radctl namelist is not supported'
!        CALL message('',message_text)
!        CALL finish('setup_radiation','Run terminated io2')
       END SELECT
      !
      ! --- Check aerosol
      ! 
      SELECT CASE (iaero)
      CASE(0)
        CALL message('','iaero= 0 --> no aerosol in radiation')
!      CASE(1)
!        CALL message('','iaero= 1 --> prognostic aerosol (sub model)')
!      CASE(2)
!        CALL message('','iaero= 2 --> Tanre aerosol climatology')
      CASE(3)
        CALL message('','iaero= 3 --> Kinne climatology')
        CALL su_aero_kinne(nb_sw)
        CALL read_aero_kinne(nb_sw)
!      CASE(5)
!        CALL message('','iaero= 5 --> Kinne climatology + Stenchikov volcanic aerosol')
!        CALL su_aero_kinne(nb_sw)
!        CALL su_aero_volc(nb_sw)
      CASE default
        WRITE (message_text, '(a,i2,a)') &
             'iaero=', iaero, ' in radctl namelist is not supported'
        CALL message('',message_text)
        CALL finish('setup_radiation','Run terminated iaero')
      END SELECT
      !
      ! --- Check annual cycle
      ! 
!      SELECT CASE (nmonth)
!      CASE(0)
!        CALL message('','nmonth=0 --> annual cycle on')
!      CASE(1:12)
!        WRITE (message_text, '(a,i2.2,a)') &
!             'nmonth = ', nmonth, ' --> perpetual month'
!        CALL message('',message_text)
!      CASE default
!        WRITE (message_text, '(a,i2,a)') &
!             'nmonth=', nmonth, ' in radctl namelist is not supported'
!        CALL message('',message_text)
!        CALL finish('setup_radiation','Run terminated nmonth')
!      END SELECT
      !
      ! --- Check Shortwave Model
      ! 
!      IF (l_srtm) THEN
!        CALL message('','l_srtm =.TRUE.  --> USE AER RRTM Shortwave Model')
!      ELSE
!        CALL message('','l_srtm =.FALSE. --> USE ECHAM5 Shortwave')
!      ENDIF
      !
      ! --- Check Longwave Model
      ! 
!      IF (l_lrtm) THEN
!        CALL message('','l_lrtm =.TRUE.  --> USE New (V4) LRTM Model')
!      ELSE
!        CALL message('','l_lrtm =.FALSE. --> USE ECHAM5 RRTM')
!      ENDIF
      !
      ! --- Check New Cloud Optics
      ! 
!      IF (l_newoptics) THEN
!        CALL message('','l_newoptics =.TRUE.  --> USE New Cloud Optics')
!      ELSE
!        CALL message('','l_newoptics =.FALSE. --> USE ECHAM5 Cloud Optical Properties')
!      ENDIF
      !
      ! --- Check solar constant
      !
!      SELECT CASE (isolrad)
!      CASE (0) 
!        CALL message('','isolrad = 0 --> standard rrtm solar constant')
!      CASE (1) 
!        CALL message('','isolrad = 1 --> time dependent spectrally resolved solar constant read from file')
!        CALL su_solar_irradiance(nb_sw)
!      CASE (2) 
!        CALL message('','isolrad = 2 --> preindustrial solar constant')
!      CASE (3) 
!        CALL message('','isolrad = 3 --> solar constant for amip runs')
!      CASE default 
!        WRITE (message_text, '(a,i3,a)') &
!             'Run terminated isolrad = ', isolrad, ' not supported'
!        CALL message('',message_text)
!        CALL finish('setup_radiation', message_text)
!      END SELECT
      !
      ! --- Check diurnal cycle
      ! 
!      IF (ldiur) THEN
!        CALL message('','ldiur =.TRUE.  --> diurnal cycle on')
!      ELSE
!        CALL message('','ldiur =.FALSE. --> diurnal cycle off')
!      ENDIF
      !
      ! --- Check for diagnosis of instantaneous aerosol radiative forcing
      ! 
!      CALL message('','instantaneous forcing diagnostic:')
!      WRITE (message_text,'(a16,L3,a18,L3)')       &
!           ' solar radiation: ',   lradforcing(1), &
!           ' thermal radiation: ', lradforcing(2)
!      CALL message('',message_text)
      !
      ! --- Check perpetual orbit
      ! 
!      IF (yr_perp.NE.-99999)  lyr_perp = .TRUE.
!      CALL p_bcast (lyr_perp, p_io)

!      IF (lyr_perp) THEN
!        IF (l_orbvsop87) THEN
!          WRITE (message_text, '(a,i0,a)') &
!               'yr_perp=', yr_perp, ' --> perpetual year for orbit'
!          CALL message('',message_text)
!        ELSE
!          WRITE (message_text, '(a,i0,a,l1,a)') &
!               'yr_perp = ', yr_perp, ' l_orbvsop87 = ',l_orbvsop87,' not allowed!'
!          CALL message('',message_text)
!          CALL finish('setup_radiation', &
!               ' yr_perp.ne.-99999 cannot run  PCMDI-orbit (l_orbvsop87=.F.).')
!        END IF
!      END IF
      !
      ! 4.0 Initialization for radiation
      ! -------------------------------
      !
      ! --- resolution/run dependent cloud optical parameters (tuning)
      !
       IF (l_newoptics .AND. l_srtm) THEN
         CALL setup_newcld_optics
       ELSE 
         CALL setup_cloud_optics(lcouple, lipcc, nn, cetah)
       END IF
      !
      ! --- Ozone climatology
      ! 
!      IF (io3==3) CALL su_o3clim_3
      !
!    ELSE
!      CALL message('','lrad = .FALSE. --> no radiation')
     ENDIF

  END SUBROUTINE setup_radiation

END MODULE mo_radiation
