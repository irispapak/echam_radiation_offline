#ifdef __xlC__
@PROCESS STRICT
#endif
MODULE mo_rrtab

  ! J.-J. Morcrette, ECMWF, 1998-07-14, original version
  ! L. Kornblueh, MPI, 2008-11-17, optimizations

  USE mo_kind  , ONLY : dp

  IMPLICIT NONE

  SAVE

  REAL(dp), PARAMETER :: bpade = 1.0_dp/0.278_dp

  REAL(dp) :: trans(0:5000)
!cdir duplicate(trans,256)

CONTAINS

  SUBROUTINE surrtab

    ! J.-J. Morcrette, ECMWF, 1998-07-14, original version
    ! M. Giorgetta, MPI, 2000-02-25, comment added
    ! L. Kornblueh, MPI, removed definition of bpade 
    
    USE mo_kind,  ONLY : dp
    
    !     local integer scalars
    INTEGER :: itr
  
    !     local real scalars
    REAL(dp) :: ztau, ztfn
    
    trans(0)    = 1.0_dp
    trans(5000) = 0.0_dp
    DO itr = 1, 4999
      ztfn = REAL(itr,dp)/5000.0_dp
      ztau = bpade*ztfn/(1.0_dp-ztfn)
      trans(itr) = EXP(-ztau)
    ENDDO
!CDIR DU_UPDATE(TRANS)
    
  END SUBROUTINE surrtab
  
END MODULE mo_rrtab
