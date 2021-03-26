MODULE mo_memory_g3b

  USE mo_kind,          ONLY: dp
  USE mo_decomposition, ONLY: ldc=>local_decomposition

  REAL(dp), POINTER, PUBLIC :: geosp(:,:,:)

CONTAINS

  SUBROUTINE construct_g3b

    ALLOCATE(geosp(ldc%nlon,ldc%nlat,ldc%ntime))

  END SUBROUTINE construct_g3b
END MODULE mo_memory_g3b
