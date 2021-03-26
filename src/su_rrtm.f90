#ifdef __xlC__
@PROCESS STRICT
#endif
SUBROUTINE su_rrtm

  !============================================================================
  !
  !- Description:
  !
  !   Initialisation of RRTM modules
  !
  !   M.A. Giorgetta, MPI, March 2000
  !
  !============================================================================

  USE mo_rrtab, ONLY: surrtab

  IMPLICIT NONE

  EXTERNAL :: surrtpk, surrtrf, surrtftr, surrtbg2, surrta

  ! init RRTM modules
  ! -----------------
  CALL surrtab        ! mo_rrtab
  CALL surrtpk        ! mo_rrtwn
  CALL surrtrf        ! mo_rrtrf
  CALL surrtftr       ! mo_rrtftr
  CALL surrtbg2       ! mo_rrtbg2
  CALL surrta         ! mo_rrtaN (N=1:16) from file 

END SUBROUTINE su_rrtm
