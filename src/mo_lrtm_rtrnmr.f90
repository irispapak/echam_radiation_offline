!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_lw/src/rrtmg_lw_rtrnmr.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.7 $
!     created:   $Date: 2009/11/12 20:52:26 $
!
#ifdef __xlC__
@PROCESS HOT
@PROCESS NOSTRICT
#endif
module mo_lrtm_rtrnmr

  !  --------------------------------------------------------------------------
  ! |                                                                          |
  ! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).  |
  ! |  This software may be used, copied, or redistributed as long as it is    |
  ! |  not sold and this copyright notice is reproduced on each copy made.     |
  ! |  This model is provided as is without any express or implied warranties. |
  ! |                       (http://www.rtweb.aer.com/)                        |
  ! |                                                                          |
  !  --------------------------------------------------------------------------

  ! ------- Modules -------

  use mo_kind,      only: wp, widx
  use mo_constants, only: api
  use parrrtm,      only: mg, nbndlw, ngptlw, delwave, ngs

  implicit none

  real(wp), parameter :: fluxfac = 2.0e+04_wp * api

  integer(widx), parameter    :: ntbl = 10000
  real(wp), parameter         :: tblint = 10000.0_wp
  real(wp), dimension(0:ntbl) :: tau_tbl
  real(wp), dimension(0:ntbl) :: exp_tbl
  real(wp), dimension(0:ntbl) :: tfn_tbl
  real(wp)                    :: bpade

  public lrtm_rtrnmr

contains

  !-----------------------------------------------------------------------------
  subroutine lrtm_rtrnmr(kproma, kbdim, nlayers, &
       istart, iend, iout, semiss, ncbands, &
       cldfrac, taucloud, planklay, planklev, plankbnd, &
       pwvcm, fracs, taut, & 
       totuflux, totdflux, fnet, &
       totuclfl, totdclfl, fnetc, &
       idrv, dplankbnd_dt, dtotuflux_dt, dtotuclfl_dt )
    !-----------------------------------------------------------------------------
    !
    !  Original version:   E. J. Mlawer, et al. RRTM_V3.0
    !  Revision for GCMs:  Michael J. Iacono; October, 2002
    !  Revision for F90:  Michael J. Iacono; June, 2006
    !  Revision for dFdT option: M. J. Iacono and E. J. Mlawer, November 2009
    !
    !  This program calculates the upward fluxes, downward fluxes, and
    !  heating rates for an arbitrary clear or cloudy atmosphere.  The input
    !  to this program is the atmospheric profile, all Planck function
    !  information, and the cloud fraction by layer.  A variable diffusivity 
    !  angle (SECDIFF) is used for the angle integration.  Bands 2-3 and 5-9 
    !  use a value for SECDIFF that varies from 1.50 to 1.80 as a function of 
    !  the column water vapor, and other bands use a value of 1.66.  The Gaussian 
    !  weight appropriate to this angle (WTDIFF=0.5) is applied here.  Note that 
    !  use of the emissivity angle for the flux integration can cause errors of 
    !  1 to 4 W/m2 within cloudy layers.  
    !  Clouds are treated with a maximum-random cloud overlap method.
    !  This subroutine also provides the optional capability to calculate
    !  the derivative of upward flux respect to surface temperature using
    !  the pre-tabulated derivative of the Planck function with respect to 
    !  temperature integrated over each spectral band.
    !***************************************************************************

    ! ------- Declarations -------

    ! ----- Input -----
    integer(widx), intent(in) :: kproma          ! number of columns
    integer(widx), intent(in) :: kbdim           ! memory extent
    integer(widx), intent(in) :: nlayers         ! total number of layers
    integer(widx), intent(in) :: istart          ! beginning band of calculation
    integer(widx), intent(in) :: iend            ! ending band of calculation
    integer, intent(in) :: iout            ! output option flag

    ! Atmosphere
    !    Dimensions: (0:nlayers)
    real(wp), intent(in) :: pwvcm(:)                         ! precipitable water vapor (cm)
    real(wp), intent(in) :: semiss(kbdim,nbndlw)             ! lw surface emissivity
    !    Dimensions: (nbndlw)
    real(wp), intent(in) :: planklay(kproma,nlayers,nbndlw)  ! 
    !    Dimensions: (nlayers,nbndlw)
    real(wp), intent(in) :: planklev(kproma,0:nlayers,nbndlw)! 
    !    Dimensions: (0:nlayers,nbndlw)
    real(wp), intent(in) :: plankbnd(kproma,nbndlw)          ! 
    !    Dimensions: (nbndlw)
    real(wp), intent(in) :: fracs(kproma,nlayers,ngptlw)     ! 
    !    Dimensions: (kproma,nlayers,ngptlw)
    real(wp), intent(in) :: taut(kproma,nlayers,ngptlw)      ! gaseous + aerosol optical depths
    !    Dimensions: (kproma,nlayers,ngptlw)

    ! Clouds
    integer, intent(in) :: ncbands(:)                        ! number of cloud spectral bands
    real(wp), intent(in) :: cldfrac(kbdim,nlayers)           ! layer cloud fraction
    !    Dimensions: (nlayers)
    real(wp), intent(in) :: taucloud(kproma,nlayers,nbndlw)  ! layer cloud optical depth
    !    Dimensions: (nlayers,nbndlw)
    integer, intent(in) :: idrv                              ! flag for calculation of dF/dt from 
    ! Planck derivative [0=off, 1=on]
    real(wp), intent(in) :: dplankbnd_dt(kproma,nbndlw)      ! derivative of Planck function wrt temp
    !    Dimensions: (nbndlw)

    ! ----- Output -----
    real(wp), intent(out) :: totuflux(kproma,0:nlayers)      ! upward longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    real(wp), intent(out) :: totdflux(kproma,0:nlayers)      ! downward longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    real(wp), intent(out) :: fnet(:,0:)                      ! net longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    real(wp), intent(out) :: totuclfl(kproma,0:nlayers)      ! clear sky upward longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    real(wp), intent(out) :: totdclfl(kproma,0:nlayers)      ! clear sky downward longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    real(wp), intent(out) :: fnetc(kproma,0:nlayers)         ! clear sky net longwave flux (w/m2)
    !    Dimensions: (0:nlayers)
    real(wp), intent(out) :: dtotuflux_dt(kproma,0:nlayers)  ! change in upward longwave flux (w/m2/k)
    ! with respect to surface temperature
    !    Dimensions: (0:nlayers)
    real(wp), intent(out) :: dtotuclfl_dt(kproma,0:nlayers)  ! change in upward longwave flux (w/m2/k)
    ! with respect to surface temperature
    !    Dimensions: (0:nlayers)

    ! ----- Local -----
    ! Declarations for radiative transfer
    real(wp) :: atot(kproma,nlayers)
    real(wp) :: atrans(kproma,nlayers)
    real(wp) :: bbugas(kproma,nlayers)
    real(wp) :: bbutot(kproma,nlayers)
    real(wp) :: clrurad(kproma,0:nlayers)
    real(wp) :: clrdrad(kproma,0:nlayers)
    real(wp) :: uflux
    real(wp) :: dflux
    real(wp) :: urad(kproma,0:nlayers)
    real(wp) :: drad(kproma,0:nlayers)
    real(wp) :: uclfl
    real(wp) :: dclfl
    real(wp) :: odcld(kproma,nlayers,nbndlw)

    real(wp) :: secdiff(kproma,nbndlw)          ! secant of diffusivity angle
    real(wp) :: a0(nbndlw),a1(nbndlw),a2(nbndlw)! diffusivity angle adjustment coefficients
    real(wp) :: wtdiff, rec_6
    real(wp) :: radld(kproma), radclrd(kproma), plfrac, blay, dplankup(kproma), dplankdn(kproma)
    real(wp) :: odepth(kproma), odtot(kproma), odtot_rec, gassrc(kproma), odepth_rec, ttot
    real(wp) :: tblind, tfactot, bbd(kproma), bbdtot(kproma), tfacgas, transc, tausfac
    real(wp) :: rad0, reflect, radlu(kproma), radclru(kproma)

    real(wp) :: duflux_dt
    real(wp) :: duclfl_dt
    real(wp) :: d_urad_dt(kproma,0:nlayers)
    real(wp) :: d_clrurad_dt(kproma,0:nlayers)
    real(wp) :: d_rad0_dt, d_radlu_dt(kproma), d_radclru_dt(kproma)

    integer(widx) :: icldidx(kproma,nlayers),iclridx(kproma,nlayers)             ! flag for cloud in layer
    integer(widx) :: icldctr(nlayers),iclrctr(nlayers)

    integer(widx) :: ibnd, ib(kproma), iband, lay, lev    ! loop indices
    integer(widx) :: igc                          ! g-point interval counter
    integer :: iclddn(kproma)               ! flag for cloud in down path
    integer(widx) :: ittot, itgas, itr            ! lookup table indices
    integer(widx) :: ipat(16,0:2)
    integer(widx) :: jl,n,m
    integer(widx) :: icnt1,icnt2,icnt3
    integer(widx) :: idx1(kproma),idx2(kproma),idx3(kproma)

    ! Declarations for cloud overlap adjustment
    real(wp) :: faccld1(kproma,nlayers+1),faccld2(kproma,nlayers+1)
    real(wp) :: facclr1(kproma,nlayers+1),facclr2(kproma,nlayers+1)
    real(wp) :: faccmb1(kproma,nlayers+1),faccmb2(kproma,nlayers+1)
    real(wp) :: faccld1d(kproma,0:nlayers),faccld2d(kproma,0:nlayers)
    real(wp) :: facclr1d(kproma,0:nlayers),facclr2d(kproma,0:nlayers)
    real(wp) :: faccmb1d(kproma,0:nlayers),faccmb2d(kproma,0:nlayers)
    !--------------------------------------------------------------------------
    !
    ! Maximum/Random cloud overlap variables
    ! for upward radiative transfer
    !  facclr2  fraction of clear radiance from previous layer that needs to
    !           be switched to cloudy stream
    !  facclr1  fraction of the radiance that had been switched in the previous
    !           layer from cloudy to clear that needs to be switched back to
    !           cloudy in the current layer
    !  faccld2  fraction of cloudy radiance from previous layer that needs to
    !           be switched to clear stream
    !  faccld1  fraction of the radiance that had been switched in the previous
    !           layer from clear to cloudy that needs to be switched back to
    !           clear in the current layer
    ! for downward radiative transfer
    !  facclr2d fraction of clear radiance from previous layer that needs to
    !           be switched to cloudy stream
    !  facclr1d fraction of the radiance that had been switched in the previous
    !           layer from cloudy to clear that needs to be switched back to
    !           cloudy in the current layer
    !  faccld2d fraction of cloudy radiance from previous layer that needs to
    !           be switched to clear stream
    !  faccld1d fraction of the radiance that had been switched in the previous
    !           layer from clear to cloudy that needs to be switched back to
    !           clear in the current layer
    !
    !--------------------------------------------------------------------------

    real(wp) :: fmax, fmin, rat1(kproma), rat2(kproma), zodcld(kproma)
    real(wp) :: clrradd(kproma), cldradd(kproma), clrradu(kproma), cldradu(kproma), oldclr, oldcld
    real(wp) :: rad(kproma), cldsrc, radmod

    integer(widx) :: istcld(kproma,nlayers+1),istcldd(kproma,0:nlayers)
    integer(widx) :: ilookup1(kproma),ilookup2(kproma)
    ! ------- Definitions -------
    ! input
    !    nlayers                      ! number of model layers
    !    ngptlw                       ! total number of g-point subintervals
    !    nbndlw                       ! number of longwave spectral bands
    !    ncbands                      ! number of spectral bands for clouds
    !    secdiff                      ! diffusivity angle
    !    wtdiff                       ! weight for radiance to flux conversion
    !    pavel                        ! layer pressures (mb)
    !    tavel                        ! layer temperatures (k)
    !    tz                           ! level (interface) temperatures(mb)
    !    tbound                       ! surface temperature (k)
    !    cldfrac                      ! layer cloud fraction
    !    taucloud                     ! layer cloud optical depth
    !    itr                          ! integer look-up table index
    !    icldidx                      ! index for cloudy columns/layers
    !    cclridx                      ! index for clear columns/layers
    !    iclddn                       ! flag for cloud in column at any layer
    !    semiss                       ! surface emissivities for each band
    !    reflect                      ! surface reflectance
    !    bpade                        ! 1/(pade constant)
    !    tau_tbl                      ! clear sky optical depth look-up table
    !    exp_tbl                      ! exponential look-up table for transmittance
    !    tfn_tbl                      ! tau transition function look-up table

    ! local
    !    atrans                       ! gaseous absorptivity
    !    atot                         ! combined gaseous and cloud absorptivity
    !    odclr                        ! clear sky (gaseous) optical depth
    !    odcld                        ! cloud optical depth
    !    odtot                        ! optical depth of gas and cloud
    !    tfacgas                      ! gas-only pade factor, used for planck fn
    !    tfactot                      ! gas and cloud pade factor, used for planck fn
    !    bbdgas                       ! gas-only planck function for downward rt
    !    bbugas                       ! gas-only planck function for upward rt
    !    bbdtot                       ! gas and cloud planck function for downward rt
    !    bbutot                       ! gas and cloud planck function for upward calc.
    !    gassrc                       ! source radiance due to gas only
    !    radlu                        ! spectrally summed upward radiance 
    !    radclru                      ! spectrally summed clear sky upward radiance 
    !    urad                         ! upward radiance by layer
    !    clrurad                      ! clear sky upward radiance by layer
    !    radld                        ! spectrally summed downward radiance 
    !    radclrd                      ! spectrally summed clear sky downward radiance 
    !    drad                         ! downward radiance by layer
    !    clrdrad                      ! clear sky downward radiance by layer
    !    d_radlu_dt                   ! spectrally summed upward radiance 
    !    d_radclru_dt                 ! spectrally summed clear sky upward radiance 
    !    d_urad_dt                    ! upward radiance by layer
    !    d_clrurad_dt                 ! clear sky upward radiance by layer

    ! output
    !    totuflux                     ! upward longwave flux (w/m2)
    !    totdflux                     ! downward longwave flux (w/m2)
    !    fnet                         ! net longwave flux (w/m2)
    !    totuclfl                     ! clear sky upward longwave flux (w/m2)
    !    totdclfl                     ! clear sky downward longwave flux (w/m2)
    !    fnetc                        ! clear sky net longwave flux (w/m2)
    !    dtotuflux_dt                 ! change in upward longwave flux (w/m2/k)
    !                                 ! with respect to surface temperature
    !    dtotuclfl_dt                 ! change in clear sky upward longwave flux (w/m2/k)
    !                                 ! with respect to surface temperature


    ! These arrays indicate the spectral 'region' (used in the 
    ! calculation of ice cloud optical depths) corresponding
    ! to each spectral band.  See cldprop.f for more details.
    data ipat /1,1,1,1,1,1,1,1,1, 1, 1, 1, 1, 1, 1, 1, &
         1,2,3,3,3,4,4,4,5, 5, 5, 5, 5, 5, 5, 5, &
         1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/

    ! This secant and weight corresponds to the standard diffusivity 
    ! angle.  This initial value is redefined below for some bands.
    data wtdiff /0.5_wp/
    data rec_6 /0.166667_wp/

    ! Reset diffusivity angle for Bands 2-3 and 5-9 to vary (between 1.50
    ! and 1.80) as a function of total column water vapor.  The function
    ! has been defined to minimize flux and cooling rate errors in these bands
    ! over a wide range of precipitable water values.
    data a0 / 1.66_wp,  1.55_wp,  1.58_wp,  1.66_wp, &
         1.54_wp, 1.454_wp,  1.89_wp,  1.33_wp, &
         1.668_wp,  1.66_wp,  1.66_wp,  1.66_wp, &
         1.66_wp,  1.66_wp,  1.66_wp,  1.66_wp /
    data a1 / 0.00_wp,  0.25_wp,  0.22_wp,  0.00_wp, &
         0.13_wp, 0.446_wp, -0.10_wp,  0.40_wp, &
         -0.006_wp,  0.00_wp,  0.00_wp,  0.00_wp, &
         0.00_wp,  0.00_wp,  0.00_wp,  0.00_wp /
    data a2 / 0.00_wp, -12.0_wp, -11.7_wp,  0.00_wp, &
         -0.72_wp,-0.243_wp,  0.19_wp,-0.062_wp, &
         0.414_wp,  0.00_wp,  0.00_wp,  0.00_wp, &
         0.00_wp,  0.00_wp,  0.00_wp,  0.00_wp /


    do ibnd = 1,nbndlw
      if (ibnd.eq.1 .or. ibnd.eq.4 .or. ibnd.ge.10) then
         secdiff(1:kproma,ibnd) = 1.66_wp
      else
         secdiff(1:kproma,ibnd) = a0(ibnd) + a1(ibnd)*exp(a2(ibnd)*pwvcm(1:kproma))
         secdiff(1:kproma,ibnd) = MAX(MIN(secdiff(1:kproma,ibnd),1.80_wp),1.50_wp)
      endif
    enddo

    urad(:,0) = 0.0_wp
    drad(:,0) = 0.0_wp
    totuflux(:,0) = 0.0_wp
    totdflux(:,0) = 0.0_wp
    clrurad(:,0) = 0.0_wp
    clrdrad(:,0) = 0.0_wp
    totuclfl(:,0) = 0.0_wp
    totdclfl(:,0) = 0.0_wp


    if (idrv .eq. 1) then
       d_urad_dt(:,0) = 0.0_wp
       d_clrurad_dt(:,0) = 0.0_wp
       dtotuflux_dt(:,0) = 0.0_wp
       dtotuclfl_dt(:,0) = 0.0_wp
    endif

    faccld1(:,:) = 0.0_wp
    faccld2(:,:) = 0.0_wp
    facclr1(:,:) = 0.0_wp
    facclr2(:,:) = 0.0_wp
    faccmb1(:,:) = 0.0_wp
    faccmb2(:,:) = 0.0_wp
    faccld1d(:,:) = 0.0_wp
    faccld2d(:,:) = 0.0_wp
    facclr1d(:,:) = 0.0_wp
    facclr2d(:,:) = 0.0_wp
    faccmb1d(:,:) = 0.0_wp
    faccmb2d(:,:) = 0.0_wp

    do lay = 1, nlayers
       icldctr(lay) = 0
       iclrctr(lay) = 0

       urad(:,lay) = 0.0_wp
       drad(:,lay) = 0.0_wp
       totuflux(:,lay) = 0.0_wp
       totdflux(:,lay) = 0.0_wp
       clrurad(:,lay) = 0.0_wp
       clrdrad(:,lay) = 0.0_wp
       totuclfl(:,lay) = 0.0_wp
       totdclfl(:,lay) = 0.0_wp

       if (idrv .eq. 1) then
          d_urad_dt(:,lay) = 0.0_wp
          d_clrurad_dt(:,lay) = 0.0_wp
          dtotuflux_dt(:,lay) = 0.0_wp
          dtotuclfl_dt(:,lay) = 0.0_wp
       endif

       do jl = 1, kproma  ! loop over columns
          if (cldfrac(jl,lay) .ge. 1.e-6_wp) then
             icldctr(lay) = icldctr(lay) + 1
             icldidx(icldctr(lay),lay) = jl
             do iband = 1, ncbands(jl)
                odcld(jl,lay,iband) = secdiff(jl,iband) * taucloud(jl,lay,iband)
             enddo
          else
             iclrctr(lay) = iclrctr(lay) + 1
             iclridx(iclrctr(lay),lay) = jl
             do iband = 1, ncbands(jl)
                odcld(jl,lay,iband) = 0.0_wp
             enddo
          endif
       enddo
       
    enddo

    ! Maximum/Random cloud overlap parameter

    istcld(:,1) = 1
    istcldd(:,nlayers) = 1

    do lev = 1, nlayers
       m = icldctr(lev)
       do n = 1,m ! loop over cloudy kproma indices
          jl = icldidx(n,lev)
          ! Maximum/random cloud overlap
          istcld(jl,lev+1) = 0
          if (lev .eq. nlayers) then
             faccld1(jl,lev+1) = 0._wp
             faccld2(jl,lev+1) = 0._wp
             facclr1(jl,lev+1) = 0._wp
             facclr2(jl,lev+1) = 0._wp
             faccmb1(jl,lev+1) = 0._wp
             faccmb2(jl,lev+1) = 0._wp
          elseif (cldfrac(jl,lev+1) .ge. cldfrac(jl,lev)) then
             faccld1(jl,lev+1) = 0._wp
             faccld2(jl,lev+1) = 0._wp
             if (istcld(jl,lev) .eq. 1) then
                facclr1(jl,lev+1) = 0._wp
                facclr2(jl,lev+1) = 0._wp
                if (cldfrac(jl,lev) .lt. 1._wp) facclr2(jl,lev+1) = &
                     (cldfrac(jl,lev+1)-cldfrac(jl,lev))/(1._wp-cldfrac(jl,lev))
                facclr2(jl,lev) = 0._wp
                faccld2(jl,lev) = 0._wp
             else
                fmax = max(cldfrac(jl,lev),cldfrac(jl,lev-1))
                if (cldfrac(jl,lev+1) .gt. fmax) then
                   facclr1(jl,lev+1) = rat2(jl)
                   facclr2(jl,lev+1) = (cldfrac(jl,lev+1)-fmax)/(1._wp-fmax)
                elseif (cldfrac(jl,lev+1) .lt. fmax) then
                   facclr1(jl,lev+1) = (cldfrac(jl,lev+1)-cldfrac(jl,lev))/ &
                        (cldfrac(jl,lev-1)-cldfrac(jl,lev))
                   facclr2(jl,lev+1) = 0._wp
                else
                   facclr1(jl,lev+1) = rat2(jl)
                   facclr2(jl,lev+1) = 0._wp
                endif
             endif
             if (facclr1(jl,lev+1).gt.0._wp .or. facclr2(jl,lev+1).gt.0._wp) then
                rat1(jl) = 1._wp
                rat2(jl) = 0._wp
             else
                rat1(jl) = 0._wp
                rat2(jl) = 0._wp
             endif
          else
            facclr1(jl,lev+1) = 0._wp
            facclr2(jl,lev+1) = 0._wp
            if (istcld(jl,lev) .eq. 1) then
              faccld1(jl,lev+1) = 0._wp
              faccld2(jl,lev+1) = (cldfrac(jl,lev)-cldfrac(jl,lev+1))/cldfrac(jl,lev)

              facclr2(jl,lev) = 0._wp
              faccld2(jl,lev) = 0._wp
            else
              fmin = min(cldfrac(jl,lev),cldfrac(jl,lev-1))
              if (cldfrac(jl,lev+1) .le. fmin) then
                faccld1(jl,lev+1) = rat1(jl)
                faccld2(jl,lev+1) = (fmin-cldfrac(jl,lev+1))/fmin
              else
                faccld1(jl,lev+1) = (cldfrac(jl,lev)-cldfrac(jl,lev+1))/(cldfrac(jl,lev)-fmin)
                faccld2(jl,lev+1) = 0._wp
              endif
            end if
            if (faccld1(jl,lev+1).gt.0._wp .or. faccld2(jl,lev+1).gt.0._wp) then
              rat1(jl) = 0._wp
              rat2(jl) = 1._wp
            else
              rat1(jl) = 0._wp
              rat2(jl) = 0._wp
           endif
         endif
          if (lev == 1) then 
            faccmb1(jl,lev+1) = 0._wp
            faccmb2(jl,lev+1) = faccld1(jl,lev+1) * facclr2(jl,lev)
          else
            faccmb1(jl,lev+1) = facclr1(jl,lev+1) * faccld2(jl,lev) * cldfrac(jl,lev-1) 
            faccmb2(jl,lev+1) = faccld1(jl,lev+1) * facclr2(jl,lev) * (1._wp - cldfrac(jl,lev-1)) 
         endif
       end do

       m = iclrctr(lev)
!IBM* ASSERT(NODEPS)
       do n = 1,m ! loop over clear kproma inidices
          jl = iclridx(n,lev)
          istcld(jl,lev+1) = 1
       enddo
    enddo

    do lev = nlayers, 1, -1

       m = icldctr(lev)
       do n = 1,m ! loop over cloudy columns indices
          jl = icldidx(n,lev)

          istcldd(jl,lev-1) = 0
          if (lev .eq. 1) then
             faccld1d(jl,lev-1) = 0._wp
             faccld2d(jl,lev-1) = 0._wp
             facclr1d(jl,lev-1) = 0._wp
             facclr2d(jl,lev-1) = 0._wp
             faccmb1d(jl,lev-1) = 0._wp
             faccmb2d(jl,lev-1) = 0._wp
          elseif (cldfrac(jl,lev-1) .ge. cldfrac(jl,lev)) then
             faccld1d(jl,lev-1) = 0._wp
             faccld2d(jl,lev-1) = 0._wp
             if (istcldd(jl,lev) .eq. 1) then
                facclr1d(jl,lev-1) = 0._wp
                facclr2d(jl,lev-1) = 0._wp
                if (cldfrac(jl,lev) .lt. 1._wp) facclr2d(jl,lev-1) = &
                     (cldfrac(jl,lev-1)-cldfrac(jl,lev))/(1._wp-cldfrac(jl,lev))
                facclr2d(jl,lev) = 0._wp
                faccld2d(jl,lev) = 0._wp
             else
                fmax = max(cldfrac(jl,lev),cldfrac(jl,lev+1))
                if (cldfrac(jl,lev-1) .gt. fmax) then
                   facclr1d(jl,lev-1) = rat2(jl)
                   facclr2d(jl,lev-1) = (cldfrac(jl,lev-1)-fmax)/(1._wp-fmax)
                elseif (cldfrac(jl,lev-1) .lt. fmax) then
                   facclr1d(jl,lev-1) = (cldfrac(jl,lev-1)-cldfrac(jl,lev))/ &
                        (cldfrac(jl,lev+1)-cldfrac(jl,lev))
                   facclr2d(jl,lev-1) = 0.
                else
                   facclr1d(jl,lev-1) = rat2(jl)
                   facclr2d(jl,lev-1) = 0._wp
                endif
             endif
             if (facclr1d(jl,lev-1).gt.0._wp .or. facclr2d(jl,lev-1).gt.0._wp)then
                rat1(jl) = 1._wp
                rat2(jl) = 0._wp
             else
                rat1(jl) = 0._wp
                rat2(jl) = 0._wp
             endif
          else
             facclr1d(jl,lev-1) = 0._wp
             facclr2d(jl,lev-1) = 0._wp
             if (istcldd(jl,lev) .eq. 1) then
                faccld1d(jl,lev-1) = 0._wp
                faccld2d(jl,lev-1) = (cldfrac(jl,lev)-cldfrac(jl,lev-1))/cldfrac(jl,lev)
                facclr2d(jl,lev) = 0._wp
                faccld2d(jl,lev) = 0._wp
             else
                fmin = min(cldfrac(jl,lev),cldfrac(jl,lev+1))
                if (cldfrac(jl,lev-1) .le. fmin) then
                   faccld1d(jl,lev-1) = rat1(jl)
                   faccld2d(jl,lev-1) = (fmin-cldfrac(jl,lev-1))/fmin
                else
                   faccld1d(jl,lev-1) = (cldfrac(jl,lev)-cldfrac(jl,lev-1))/(cldfrac(jl,lev)-fmin)
                   faccld2d(jl,lev-1) = 0._wp
                endif
             endif
             if (faccld1d(jl,lev-1).gt.0._wp .or. faccld2d(jl,lev-1).gt.0._wp)then
                rat1(jl) = 0._wp
                rat2(jl) = 1._wp
             else
                rat1(jl) = 0._wp
                rat2(jl) = 0._wp
             endif
          endif
          if (lev == nlayers) then
             faccmb1d(jl,lev-1) = 0._wp
             faccmb2d(jl,lev-1) = faccld1d(jl,lev-1) * facclr2d(jl,lev)
          else  
             faccmb1d(jl,lev-1) = facclr1d(jl,lev-1) * faccld2d(jl,lev) * cldfrac(jl,lev+1) 
             faccmb2d(jl,lev-1) = faccld1d(jl,lev-1) * facclr2d(jl,lev) * (1._wp - cldfrac(jl,lev+1))
          endif
       end do

       m = iclrctr(lev)
!IBM* ASSERT(NODEPS)
       do n = 1,m ! loop over clear columns
          jl = iclridx(n,lev)
          istcldd(jl,lev-1) = 1
       enddo
    enddo

      igc = 1
      ! Loop over frequency bands.
      do iband = istart, iend

        ! Reinitialize g-point counter for each band if output for each band is requested.
        if (iout.gt.0.and.iband.ge.2) igc = ngs(iband-1)+1
        do jl = 1,kproma
           if (ncbands(jl) .eq. 1) then
              ib(jl) = ipat(iband,0)
           elseif (ncbands(jl) .eq.  5) then
              ib(jl) = ipat(iband,1)
           elseif (ncbands(jl) .eq. 16) then
              ib(jl) = ipat(iband,2)
           endif
        end do

        ! Loop over g-channels.
1000    continue

        ! Radiative transfer starts here.
        radld(:) = 0._wp
        radclrd(:) = 0._wp
        iclddn(:) = 0

        ! Downward radiative transfer loop.  

        do lev = nlayers, 1, -1
           do jl = 1,kproma  ! loop over columns

              blay = planklay(jl,lev,iband)
              dplankup(jl) = planklev(jl,lev,iband) - blay
              dplankdn(jl) = planklev(jl,lev-1,iband) - blay
              odepth(jl) = max(0.0_wp,secdiff(jl,iband) * taut(jl,lev,igc))
              zodcld(jl)  = odcld(jl,lev,ib(jl))
           end do

           ! Cloudy layer
           icnt1 = 0
           icnt2 = 0
           icnt3 = 0
           m = icldctr(lev)
           do n = 1,m ! loop over cloudy kproma indices
              jl = icldidx(n,lev)
              
              iclddn(jl) = 1
              odtot(jl) = odepth(jl) + zodcld(jl)
              if (odtot(jl) .lt. 0.06_wp) then
                 icnt1 = icnt1 + 1
                 idx1(icnt1) = jl
              elseif (odepth(jl) .le. 0.06_wp) then
                 icnt2 = icnt2 + 1
                 idx2(icnt2) = jl
              else
                 icnt3 = icnt3 + 1
                 idx3(icnt3) = jl
              end if
           end do
              
!IBM* ASSERT(NODEPS)
           do n = 1,icnt1
              jl = idx1(n)
              atrans(jl,lev) = odepth(jl) - 0.5_wp*odepth(jl)*odepth(jl)
              atot(jl,lev)   =  odtot(jl) - 0.5_wp*odtot(jl)*odtot(jl)

              odepth_rec = rec_6*odepth(jl)
              odtot_rec  = rec_6*odtot(jl)

              plfrac = fracs(jl,lev,igc)
              blay   = planklay(jl,lev,iband)

              bbd(jl)    = plfrac * (blay + dplankdn(jl) * odepth_rec)
              bbdtot(jl) = plfrac * (blay + dplankdn(jl) * odtot_rec)
              gassrc(jl) = atrans(jl,lev) * bbd(jl)
                 
              bbugas(jl,lev) = plfrac * (blay + dplankup(jl) * odepth_rec)
              bbutot(jl,lev) = plfrac * (blay + dplankup(jl) * odtot_rec)
           end do

           if (icnt2 > 0) then ! trick compiler to not fuse the following loops
!IBM NOVECTOR
!IBM* ASSERT(NODEPS)
              do n = 1,icnt2
                 jl = idx2(n)

                 odtot(jl) = odepth(jl) + zodcld(jl)

                 ! table lookups here
                 tblind = odtot(jl)/(bpade+odtot(jl))
                 ilookup1(n) = INT(tblint*tblind + 0.5_wp,widx)
              end do
           end if
           
!IBM* ASSERT(NODEPS)
           do n = 1,icnt2
              jl = idx2(n)
              ittot = ilookup1(n)
              tfactot = tfn_tbl(ittot)
              atot(jl,lev) = 1._wp - exp_tbl(ittot)
              
              odepth_rec = rec_6*odepth(jl)
              atrans(jl,lev) = odepth(jl) - 0.5_wp*odepth(jl)*odepth(jl)

              plfrac = fracs(jl,lev,igc)
              blay   = planklay(jl,lev,iband)

              bbd(jl)    = plfrac * (blay + odepth_rec * dplankdn(jl))
              bbdtot(jl) = plfrac * (blay +    tfactot * dplankdn(jl))
              gassrc(jl) = atrans(jl,lev) * bbd(jl)

              bbugas(jl,lev) = plfrac * (blay + odepth_rec * dplankup(jl))
              bbutot(jl,lev) = plfrac * (blay +    tfactot * dplankup(jl))
           end do

!IBM NOVECTOR
!IBM* ASSERT(NODEPS)
           do n = 1,icnt3
              jl = idx3(n)

              ! table lookups here, separate index computation from use (P6 specific)
              tblind = odepth(jl)/(bpade+odepth(jl))
              ilookup1(n) = INT(tblint*tblind+0.5_wp,widx)
           end do

           if (icnt3 > 0) then ! trick the compiler not to fuse the loops
!IBM NOVECTOR
!IBM* ASSERT(NODEPS)
              do n = 1,icnt3
                 jl = idx3(n)
                 itgas = ilookup1(n)
                 odepth(jl) = tau_tbl(itgas)
                 odtot(jl) = odepth(jl) + zodcld(jl)

                 ! more table lookups here
                 tblind = odtot(jl)/(bpade+odtot(jl))
                 ilookup2(n) = INT(tblint*tblind + 0.5_wp,widx)
              end do
           endif

!IBM* ASSERT(NODEPS)
           do n = 1,icnt3
              jl = idx3(n)
              itgas = ilookup1(n)
              ittot = ilookup2(n)
              tfacgas = tfn_tbl(itgas)
              tfactot = tfn_tbl(ittot)
              atrans(jl,lev) = 1._wp - exp_tbl(itgas)
              atot(jl,lev)   = 1._wp - exp_tbl(ittot)

              plfrac = fracs(jl,lev,igc)
              blay   = planklay(jl,lev,iband)

              bbd(jl)    = plfrac * (blay + tfacgas * dplankdn(jl))
              bbdtot(jl) = plfrac * (blay + tfactot * dplankdn(jl))
              gassrc(jl) = atrans(jl,lev) * bbd(jl)

              bbugas(jl,lev) = plfrac * (blay + tfacgas * dplankup(jl))
              bbutot(jl,lev) = plfrac * (blay + tfactot * dplankup(jl))
           end do

           m = icldctr(lev)
!IBM* ASSERT(NODEPS)
           do n = 1,m ! loop over cloudy kproma indices
              jl = icldidx(n,lev)

              if (istcldd(jl,lev) .eq. 1) then
                 cldradd(jl) = cldfrac(jl,lev) * radld(jl)
                 clrradd(jl) = radld(jl) - cldradd(jl)
                 rad(jl) = 0._wp
              endif
           end do

!IBM* ASSERT(NODEPS)
           do n = 1,m ! loop over cloudy kproma indices
              jl = icldidx(n,lev)

              ttot = 1._wp - atot(jl,lev)
              cldsrc = bbdtot(jl) * atot(jl,lev)
              cldradd(jl) = cldradd(jl) * ttot + cldfrac(jl,lev) * cldsrc
              clrradd(jl) = clrradd(jl) * (1._wp-atrans(jl,lev)) + (1._wp - cldfrac(jl,lev)) * gassrc(jl)
              radld(jl) = cldradd(jl) + clrradd(jl)
              drad(jl,lev-1) = drad(jl,lev-1) + radld(jl)

              radmod = rad(jl) * &
                   (facclr1d(jl,lev-1) * (1.-atrans(jl,lev)) + &
                   faccld1d(jl,lev-1) *  ttot) - &
                   faccmb1d(jl,lev-1) * gassrc(jl) + &
                   faccmb2d(jl,lev-1) * cldsrc
              
              oldcld = cldradd(jl) - radmod
              oldclr = clrradd(jl) + radmod
              rad(jl) = -radmod + facclr2d(jl,lev-1)*oldclr - faccld2d(jl,lev-1)*oldcld
              cldradd(jl) = cldradd(jl) + rad(jl)
              clrradd(jl) = clrradd(jl) - rad(jl)
           end do

           ! Clear layer
           icnt2 = 0
           icnt3 = 0
           m = iclrctr(lev)
           do n = 1,m ! loop over cloudy kproma indices
              jl = iclridx(n,lev)
              
              if (odepth(jl) .le. 0.06_wp) then
                 icnt2 = icnt2 + 1
                 idx2(icnt2) = jl
              else
                 icnt3 = icnt3 + 1
                 idx3(icnt3) = jl
              end if
           end do
           
!IBM* ASSERT(NODEPS)
           do n = 1,icnt2
              jl = idx2(n)

              plfrac = fracs(jl,lev,igc)
              blay   = planklay(jl,lev,iband)

              odepth_rec = rec_6*odepth(jl)
              atrans(jl,lev) = odepth(jl)-0.5_wp*odepth(jl)*odepth(jl)

              bbd(jl)        = plfrac * (blay + dplankdn(jl) * odepth_rec)
              bbugas(jl,lev) = plfrac * (blay + dplankup(jl) * odepth_rec)

              radld(jl) = radld(jl) + (bbd(jl)-radld(jl))*atrans(jl,lev)
              drad(jl,lev-1) = drad(jl,lev-1) + radld(jl)
           end do

!IBM NOVECTOR
!IBM* ASSERT(NODEPS)
           do n = 1,icnt3
              jl = idx3(n)

              ! table lookups here
              tblind = odepth(jl)/(bpade+odepth(jl))
              ilookup1(n) = INT(tblint*tblind+0.5_wp,widx)
           end do

!IBM* ASSERT(NODEPS)
           do n = 1,icnt3
              jl = idx3(n)
              itr = ilookup1(n)
              transc = exp_tbl(itr)
              tausfac = tfn_tbl(itr)
              atrans(jl,lev) = 1._wp-transc

              plfrac = fracs(jl,lev,igc)
              blay   = planklay(jl,lev,iband)

              bbd(jl)        = plfrac * (blay + tausfac * dplankdn(jl))
              bbugas(jl,lev) = plfrac * (blay + tausfac * dplankup(jl))

              radld(jl)      = radld(jl) + (bbd(jl)-radld(jl))*atrans(jl,lev)
              drad(jl,lev-1) = drad(jl,lev-1) + radld(jl)
           end do

           !  Set clear sky stream to total sky stream as long as layers
           !  remain clear.  Streams diverge when a cloud is reached (iclddn=1),
           !  and clear sky stream must be computed separately from that point.

           icnt1 = 0
           icnt2 = 0
           do jl = 1,kproma   ! loop over columns
              if (iclddn(jl).eq.1) then
                 icnt1 = icnt1 + 1
                 idx1(icnt1) = jl
              else
                 icnt2 = icnt2 + 1
                 idx2(icnt2) = jl
              end if
           end do

!IBM* ASSERT(NODEPS)
           do n = 1,icnt1
              jl = idx1(n)
              radclrd(jl) = radclrd(jl) + (bbd(jl)-radclrd(jl)) * atrans(jl,lev) 
              clrdrad(jl,lev-1) = clrdrad(jl,lev-1) + radclrd(jl)
           end do

!IBM* ASSERT(NODEPS)
           do n = 1,icnt2
              jl = idx2(n)
              radclrd(jl) = radld(jl)
              clrdrad(jl,lev-1) = drad(jl,lev-1)
           end do
        enddo

        ! Spectral emissivity & reflectance
        !  Include the contribution of spectrally varying longwave emissivity
        !  and reflection from the surface to the upward radiative transfer.
        !  Note: Spectral and Lambertian reflection are identical for the
        !  diffusivity angle flux integration used here.
        !  Note: The emissivity is applied to plankbnd and dplankbnd_dt when 
        !  they are defined in subroutine setcoef. 

        do jl=1,kproma
           rad0 = fracs(jl,1,igc) * plankbnd(jl,iband)
           !  Add in reflection of surface downward radiance.
           reflect = 1._wp - semiss(jl,iband)
           radlu(jl) = rad0 + reflect * radld(jl)
           radclru(jl) = rad0 + reflect * radclrd(jl)
           
           ! Upward radiative transfer loop.

           urad(jl,0) = urad(jl,0) + radlu(jl)
           clrurad(jl,0) = clrurad(jl,0) + radclru(jl)
        end do

        if (idrv .eq. 1) then
           do jl=1,kproma
              d_rad0_dt = fracs(jl,1,igc) * dplankbnd_dt(jl,iband)
              d_radlu_dt(jl) = d_rad0_dt
              d_urad_dt(jl,0) = d_urad_dt(jl,0) + d_radlu_dt(jl)
              d_radclru_dt(jl) = d_rad0_dt
              d_clrurad_dt(jl,0) = d_clrurad_dt(jl,0) + d_radclru_dt(jl)
           end do
        endif

        do lev = 1, nlayers
          ! Cloudy layer
           m = icldctr(lev)
!IBM* ASSERT(NODEPS)
           do n = 1,m ! loop over cloudy kproma indices
              jl = icldidx(n,lev)
              
              gassrc(jl) = bbugas(jl,lev) * atrans(jl,lev)
              if (istcld(jl,lev) .eq. 1) then
                 cldradu(jl) = cldfrac(jl,lev) * radlu(jl)
                 clrradu(jl) = radlu(jl) - cldradu(jl)
                 rad(jl) = 0._wp
              endif
           end do

!IBM* ASSERT(NODEPS)
           do n = 1,m ! loop over cloudy kproma indices
              jl = icldidx(n,lev)
              ttot = 1._wp - atot(jl,lev)
              cldsrc = bbutot(jl,lev) * atot(jl,lev)
              cldradu(jl) = cldradu(jl) * ttot + cldfrac(jl,lev) * cldsrc
              clrradu(jl) = clrradu(jl) * (1.0_wp - atrans(jl,lev)) + (1._wp - cldfrac(jl,lev)) * gassrc(jl)
              ! Total sky radiance
              radlu(jl) = cldradu(jl) + clrradu(jl)
              urad(jl,lev) = urad(jl,lev) + radlu(jl)
              radmod = rad(jl) * &
                   (facclr1(jl,lev+1)*(1.0_wp - atrans(jl,lev))+ &
                   faccld1(jl,lev+1) *  ttot) - &
                   faccmb1(jl,lev+1) * gassrc(jl) + &
                   faccmb2(jl,lev+1) * cldsrc
              oldcld = cldradu(jl) - radmod
              oldclr = clrradu(jl) + radmod
              rad(jl) = -radmod + facclr2(jl,lev+1)*oldclr - faccld2(jl,lev+1)*oldcld
              cldradu(jl) = cldradu(jl) + rad(jl)
              clrradu(jl) = clrradu(jl) - rad(jl)
           end do

           if (idrv .eq. 1) then
!IBM* ASSERT(NODEPS)
              do n = 1,m ! loop over cloudy kproma indices
                 jl = icldidx(n,lev)
                 d_radlu_dt(jl) = d_radlu_dt(jl) * cldfrac(jl,lev) * (1.0_wp - atot(jl,lev)) + &
                      d_radlu_dt(jl) * (1.0_wp - cldfrac(jl,lev)) * (1.0_wp - atrans(jl,lev))
                 d_urad_dt(jl,lev) = d_urad_dt(jl,lev) + d_radlu_dt(jl)
              end do
           endif

           ! Clear layer
           m = iclrctr(lev)
!IBM* ASSERT(NODEPS)
           do n = 1,m ! loop over clear columns
              jl = iclridx(n,lev)
              
              radlu(jl) = radlu(jl) + (bbugas(jl,lev)-radlu(jl))*atrans(jl,lev)
              urad(jl,lev) = urad(jl,lev) + radlu(jl)
           end do
           if (idrv .eq. 1) then
!IBM* ASSERT(NODEPS)
              do n = 1,m ! loop over clear columns
                 jl = iclridx(n,lev)
                 d_radlu_dt(jl) = d_radlu_dt(jl) * (1.0_wp - atrans(jl,lev))
                 d_urad_dt(jl,lev) = d_urad_dt(jl,lev) + d_radlu_dt(jl)
              end do
           endif

           !  Set clear sky stream to total sky stream as long as all layers
           !  are clear (iclddn=0).  Streams must be calculated separately at 
           !  all layers when a cloud is present (iclddn=1), because surface 
           !  reflectance is different for each stream.
           icnt1 = 0
           icnt2 = 0
           do jl = 1,kproma   ! loop over columns
              if (iclddn(jl).eq.1) then
                 icnt1 = icnt1 + 1
                 idx1(icnt1) = jl
              else
                 icnt2 = icnt2 + 1
                 idx2(icnt2) = jl
              end if
           end do

!IBM* ASSERT(NODEPS)
           do n = 1,icnt1
              jl = idx1(n)
              radclru(jl) = radclru(jl) + (bbugas(jl,lev)-radclru(jl))*atrans(jl,lev) 
              clrurad(jl,lev) = clrurad(jl,lev) + radclru(jl)
           end do

!IBM* ASSERT(NODEPS)
           do n = 1,icnt2
              jl = idx2(n)
              radclru(jl) = radlu(jl)
              clrurad(jl,lev) = urad(jl,lev)
           end do
           
           if (idrv .eq. 1) then
!IBM* ASSERT(NODEPS)
              do n = 1,icnt1
                 jl = idx1(n)
                 d_radclru_dt(jl) = d_radclru_dt(jl) * (1.0_wp - atrans(jl,lev))
                 d_clrurad_dt(jl,lev) = d_clrurad_dt(jl,lev) + d_radclru_dt(jl)
              end do
!IBM* ASSERT(NODEPS)
              do n = 1,icnt2
                 jl = idx2(n)
                 d_radclru_dt(jl) = d_radlu_dt(jl)
                 d_clrurad_dt(jl,lev) = d_urad_dt(jl,lev)
              enddo
           endif
        enddo

        ! Increment g-point counter
        igc = igc + 1
        ! Return to continue radiative transfer for all g-channels in present band
        if (igc .le. ngs(iband)) go to 1000

        ! Process longwave output from band.
        ! Calculate upward, downward, and net flux.
        do lev = nlayers, 0, -1
           do jl=kproma,1,-1
              uflux = urad(jl,lev)*wtdiff
              dflux = drad(jl,lev)*wtdiff
              urad(jl,lev) = 0.0_wp
              drad(jl,lev) = 0.0_wp
              totuflux(jl,lev) = totuflux(jl,lev) + uflux * delwave(iband)
              totdflux(jl,lev) = totdflux(jl,lev) + dflux * delwave(iband)
              uclfl = clrurad(jl,lev)*wtdiff
              dclfl = clrdrad(jl,lev)*wtdiff
              clrurad(jl,lev) = 0.0_wp
              clrdrad(jl,lev) = 0.0_wp
              totuclfl(jl,lev) = totuclfl(jl,lev) + uclfl * delwave(iband)
              totdclfl(jl,lev) = totdclfl(jl,lev) + dclfl * delwave(iband)
           end do
        enddo

        ! Calculate total change in upward flux wrt surface temperature
        if (idrv .eq. 1) then
          do lev = nlayers, 0, -1
             do jl=kproma,1,-1
                duflux_dt = d_urad_dt(jl,lev) * wtdiff
                d_urad_dt(jl,lev) = 0.0_wp
                dtotuflux_dt(jl,lev) = dtotuflux_dt(jl,lev) + duflux_dt * delwave(iband) * fluxfac
                
                duclfl_dt = d_clrurad_dt(jl,lev) * wtdiff
                d_clrurad_dt(jl,lev) = 0.0_wp
                dtotuclfl_dt(jl,lev) = dtotuclfl_dt(jl,lev) + duclfl_dt * delwave(iband) * fluxfac
             enddo
          enddo
        endif

        ! End spectral band loop
     enddo

    ! Calculate fluxes at surface
    do jl = 1, kproma  ! loop over columns
      totuflux(jl,0) = totuflux(jl,0) * fluxfac
      totdflux(jl,0) = totdflux(jl,0) * fluxfac
      fnet(jl,0) = totuflux(jl,0) - totdflux(jl,0)

      totuclfl(jl,0) = totuclfl(jl,0) * fluxfac
      totdclfl(jl,0) = totdclfl(jl,0) * fluxfac
      fnetc(jl,0) = totuclfl(jl,0) - totdclfl(jl,0)
    enddo

    ! Calculate fluxes at model levels
    do lev = 1, nlayers
      do jl = 1, kproma  ! loop over columns
        totuflux(jl,lev) = totuflux(jl,lev) * fluxfac
        totdflux(jl,lev) = totdflux(jl,lev) * fluxfac
        fnet(jl,lev) = totuflux(jl,lev) - totdflux(jl,lev)
        totuclfl(jl,lev) = totuclfl(jl,lev) * fluxfac
        totdclfl(jl,lev) = totdclfl(jl,lev) * fluxfac
        fnetc(jl,lev) = totuclfl(jl,lev) - totdclfl(jl,lev)
      enddo
    enddo

  end subroutine lrtm_rtrnmr

end module mo_lrtm_rtrnmr


