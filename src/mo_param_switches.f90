MODULE mo_param_switches

  IMPLICIT NONE

  ! M.A. Giorgetta, March 2000, lrad added
  !
  ! ----------------------------------------------------------------
  !
  ! module *mo_param_switches* switches related to the parameterisations of
  !                            diabatic processes except for radiation.
  !
  ! ----------------------------------------------------------------

  LOGICAL :: lphys        !   *true for parameterisation of diabatic processes.
!LT
  LOGICAL :: lrad = .TRUE.        !   *true for radiation.
!LT
  LOGICAL :: lvdiff       !   *true for vertical diffusion.
  LOGICAL :: lcond        !   *true for large scale condensation scheme.
  LOGICAL :: lsurf        !   *true for surface exchanges.
  LOGICAL :: lcover       !   *true for prognostic cloud cover scheme
  LOGICAL :: lconv        !   *true to allow convection
  LOGICAL :: lgwdrag      !   *true for gravity wave drag scheme
  LOGICAL :: lice         !   *true for sea-ice temperature calculation
  LOGICAL :: lconvmassfix !   *false for switching off aerosol mass fixer in conv

  INTEGER :: iconv = 1 !   *1,2,3 for different convection schemes

!++mgs : new switches for interactive cloud scheme
  LOGICAL :: lcdnc_progn  !   true for prognostic cloud droplet activation
  INTEGER :: ncd_activ    !   type of cloud droplet activation scheme (see physctl.inc)
  LOGICAL :: lice_supersat!   true for saturation over ice for cirrus clouds (former icnc=2)
  INTEGER :: nauto        !   autoconversion scheme    (1,2)
  INTEGER :: ncvmicro     !   conv microphysics scheme (0,1)

END MODULE mo_param_switches
