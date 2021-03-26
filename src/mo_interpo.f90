MODULE mo_interpo

  ! Weighting factors and indices for time interpolation

  USE mo_kind,  ONLY: dp

  IMPLICIT NONE

  PUBLIC

  REAL(dp):: wgt1, wgt2
  INTEGER :: nmw1, nmw2, nmw1cl, nmw2cl

  REAL(dp):: wgtd1, wgtd2
  INTEGER :: ndw1, ndw2

  ! weightings and indicies for radiation calculation time step

  REAL(dp), PARAMETER:: wgt1_m=0.5_dp, wgt2_m=0.5_dp
  INTEGER, PARAMETER :: nmw1_m=1, nmw2_m=1
 
END MODULE mo_interpo
