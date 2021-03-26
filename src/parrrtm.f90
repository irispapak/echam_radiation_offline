      module parrrtm

      use mo_kind, only : wp
      implicit none
      save

!------------------------------------------------------------------
! rrtmg_lw main parameters
!
! Initial version:  JJMorcrette, ECMWF, Jul 1998
! Revised: MJIacono, AER, Jun 2006
! Revised: MJIacono, AER, Aug 2007
! Revised: MJIacono, AER, Aug 2008
!------------------------------------------------------------------

!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! mxlay  :  integer: maximum number of layers
! mg     :  integer: number of original g-intervals per spectral band
! nbndlw :  integer: number of spectral bands
! maxxsec:  integer: maximum number of cross-section molecules
!                    (e.g. cfcs)
! maxinpx:  integer: 
! ngptlw :  integer: total number of reduced g-intervals for rrtmg_lw
! ngNN   :  integer: number of reduced g-intervals per spectral band
! ngsNN  :  integer: cumulative number of g-intervals per band
!------------------------------------------------------------------

      integer, parameter :: mxlay  = 203
      integer, parameter :: mg     = 16
      integer, parameter :: nbndlw = 16
      integer, parameter :: maxxsec= 4
      integer, parameter :: mxmol  = 38
      integer, parameter :: maxinpx= 38
      integer, parameter :: nmol   = 7

      integer, parameter :: ngptlw = 140 ! set to 256 for 256 gpt model
      integer, parameter :: ng1  = 10
      integer, parameter :: ng2  = 12
      integer, parameter :: ng3  = 16
      integer, parameter :: ng4  = 14
      integer, parameter :: ng5  = 16
      integer, parameter :: ng6  = 8
      integer, parameter :: ng7  = 12
      integer, parameter :: ng8  = 8
      integer, parameter :: ng9  = 12
      integer, parameter :: ng10 = 6
      integer, parameter :: ng11 = 8
      integer, parameter :: ng12 = 8
      integer, parameter :: ng13 = 4
      integer, parameter :: ng14 = 2
      integer, parameter :: ng15 = 2
      integer, parameter :: ng16 = 2

      integer, parameter :: ngs1  = 10
      integer, parameter :: ngs2  = 22
      integer, parameter :: ngs3  = 38
      integer, parameter :: ngs4  = 52
      integer, parameter :: ngs5  = 68
      integer, parameter :: ngs6  = 76
      integer, parameter :: ngs7  = 88
      integer, parameter :: ngs8  = 96
      integer, parameter :: ngs9  = 108
      integer, parameter :: ngs10 = 114
      integer, parameter :: ngs11 = 122
      integer, parameter :: ngs12 = 130
      integer, parameter :: ngs13 = 134
      integer, parameter :: ngs14 = 136
      integer, parameter :: ngs15 = 138

!------------------------------------------------------------------
! rrtmg_lw spectral information

! Initial version:  JJMorcrette, ECMWF, jul1998
! Revised: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!------------------------------------------------------------------

!  name     type     purpose
! -----  :  ----   : ----------------------------------------------
! ng     :  integer: Number of original g-intervals in each spectral band
! nspa   :  integer: For the lower atmosphere, the number of reference
!                    atmospheres that are stored for each spectral band
!                    per pressure level and temperature.  Each of these
!                    atmospheres has different relative amounts of the 
!                    key species for the band (i.e. different binary
!                    species parameters).
! nspb   :  integer: Same as nspa for the upper atmosphere
!wavenum1:  real   : Spectral band lower boundary in wavenumbers
!wavenum2:  real   : Spectral band upper boundary in wavenumbers
! delwave:  real   : Spectral band width in wavenumbers
! totplnk:  real   : Integrated Planck value for each band; (band 16
!                    includes total from 2600 cm-1 to infinity)
!                    Used for calculation across total spectrum
!totplk16:  real   : Integrated Planck value for band 16 (2600-3250 cm-1)
!                    Used for calculation in band 16 only if 
!                    individual band output requested
!totplnkderiv: real: Integrated Planck function derivative with respect
!                    to temperature for each band; (band 16
!                    includes total from 2600 cm-1 to infinity)
!                    Used for calculation across total spectrum
!totplk16deriv:real: Integrated Planck function derivative with respect
!                    to temperature for band 16 (2600-3250 cm-1)
!                    Used for calculation in band 16 only if 
!                    individual band output requested
!
! ngc    :  integer: The number of new g-intervals in each band
! ngs    :  integer: The cumulative sum of new g-intervals for each band
! ngm    :  integer: The index of each new g-interval relative to the
!                    original 16 g-intervals in each band
! ngn    :  integer: The number of original g-intervals that are 
!                    combined to make each new g-intervals in each band
! ngb    :  integer: The band index for each new g-interval
! wt     :  real   : RRTM weights for the original 16 g-intervals
! rwgt   :  real   : Weights for combining original 16 g-intervals 
!                    (256 total) into reduced set of g-intervals 
!                    (140 total)
! nxmol  :  integer: Number of cross-section molecules
! ixindx :  integer: Flag for active cross-sections in calculation
!------------------------------------------------------------------

      integer :: ng(nbndlw)
      integer :: nspa(nbndlw)
      integer :: nspb(nbndlw)

      real(wp) :: wavenum1(nbndlw)
      real(wp) :: wavenum2(nbndlw)
      real(wp) :: delwave(nbndlw)

      real(wp) :: totplnk(181,nbndlw)
      real(wp) :: totplk16(181)

      real(wp) :: totplnkderiv(181,nbndlw)
      real(wp) :: totplk16deriv(181)

      integer :: ngc(nbndlw)
      integer :: ngs(nbndlw)
      integer :: ngn(ngptlw)
      integer :: ngb(ngptlw)
      integer :: ngm(nbndlw*mg)

      real(wp) :: wt(mg)
      real(wp) :: rwgt(nbndlw*mg)

      integer :: nxmol
      integer :: ixindx(maxinpx)

      end module parrrtm

