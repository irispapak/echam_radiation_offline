MODULE mo_hyb

  USE mo_kind, ONLY: dp

  IMPLICIT NONE

  ! ------------------------------------------------------------------
  !
  ! module *mo_hyb* - *loop indices and surface-pressure independent
  !        variables associated with the vertical finite-difference scheme.
  !
  ! a.j.simmons     e.c.m.w.f.     16/11/81.
  !
  ! ------------------------------------------------------------------


  INTEGER :: nlevm1         !  (number of levels)-1.
  INTEGER :: nplev          !  *number of pressure levels.
  INTEGER :: nplvp1         !  *nplev+1.*
  INTEGER :: nplvp2         !  *nplev+2.*
  INTEGER :: nplvpa         !  *nplvp1,* or 2 if *nplev=0.*
  INTEGER :: nlmsgl         !  *nlev* - (number of sigma levels).
  INTEGER :: nlmslp         !  *nlmsgl+1.*
  INTEGER :: nlmsla         !  *nlmslp,* or 2 if *nlmslp=1.*

  REAL(dp) :: apzero            !  *reference pressure for computation of the
                            !   hybrid vertical levels.
  REAL(dp) :: rpr               !  *reciprocal of reference surface pressure.
  REAL(dp) :: rdtr              !   rd*(reference temperature).
  REAL(dp) :: apsurf            !   fixed global mean of surface pressure
  REAL(dp) :: t0icao            !  *surface temperatur of reference atmosphere
  REAL(dp) :: tsticao           !  *stratospheric temperature of reference atmosphere
  REAL(dp) :: rdtstic           !  *rd*tsticao
  REAL(dp) :: rdlnp0i           !  *rd*ln(surface pressure) of reference atmosphere
  REAL(dp) :: alrrdic           !  *lapse-rate parameter of reference atmosphere
  REAL(dp) :: rdt0ral           !  *rd*t0icao/alphaic
  REAL(dp) :: ptricao           !  *tropopause pressure of reference atmosphere
  REAL(dp) :: rdlnpti           !  *rd*ln(ptricao)
  REAL(dp) :: gsticao           !  *constant used in geopotential calculation
  REAL(dp), ALLOCATABLE :: ralpha(:)    !   rd*alpha at pressure and sigma levels.
  REAL(dp), ALLOCATABLE :: rlnpr(:)     !   rd*ln(p(k+.5)/p(k-.5)) at pressure and sigma levels.
  REAL(dp), ALLOCATABLE :: dela(:)      !   a(k+.5)-a(k-.5).
  REAL(dp), ALLOCATABLE :: delb(:)      !   b(k+.5)-b(k-.5).
  REAL(dp), ALLOCATABLE :: rddelb(:)    !   rd*delb.
  REAL(dp), ALLOCATABLE :: cpg(:)       !   a(k+.5)*b(k-.5)-b(k+.5)*a(k-.5).
  REAL(dp), ALLOCATABLE :: delpr(:)     !   p(k+.5)-p(k-.5) for reference surface pressure.
  REAL(dp), ALLOCATABLE :: rdelpr(:)    !  *reciprocal of *delpr.*
  REAL(dp), ALLOCATABLE :: ralphr(:)    !  *constant array for use by pgrad.
  REAL(dp), ALLOCATABLE :: alpham(:)    !  *constant array for use by dyn.
  REAL(dp), ALLOCATABLE :: ardprc(:)    !  *constant array for use by dyn.
  REAL(dp), ALLOCATABLE :: rlnmar(:)    !  *constant array for use by pgrad.
  REAL(dp), ALLOCATABLE :: aktlrd(:)    !  *constant array for use by conteq.
  REAL(dp), ALLOCATABLE :: altrcp(:)    !  *constant array for use by conteq.
  REAL(dp), ALLOCATABLE :: ceta(:)      !  *full hybrid vertical levels.
!LT: feste vertikale Auflösung
  REAL(dp) :: cetah(48)    !  *half hybrid vertical levels.
!  REAL(dp), ALLOCATABLE :: cetah(:)    !  *half hybrid vertical levels.
  REAL(dp), ALLOCATABLE :: bb(:,:) !  *gravity wave matrix
  REAL(dp) :: vct_a(48)
  REAL(dp) :: vct_b(48)

CONTAINS

  !+ initializes constants for vertical coordinate calculations.
  !+ $Id$

  SUBROUTINE inihyb

    ! Description:
    !
    ! Initializes constants for vertical coordinate calculations.
    !
    ! Method:
    !
    ! Compute loop indices and surface-pressure independent
    ! variables associated with the vertical finite-difference scheme.
    !
    ! Output is in module *mo_hyb*
    !
    ! Authors:
    !
    ! A. J. Simmons, ECMWF, November 1981, original source
    ! L. Kornblueh, MPI, May 1998, f90 rewrite
    ! U. Schulzweida, MPI, May 1998, f90 rewrite
    ! A. Rhodin, MPI, Jan 1999, subroutine inihyb -> module mo_hyb
    ! 
    ! for more details see file AUTHORS
    !
    
    USE mo_constants, ONLY: g, rcpd, rd
    USE mo_control,   ONLY: nlevp1, nvclev, vct

    !  Local scalars: 
    REAL(dp) :: za, zb, zetam, zetap, zp, zp0icao, zpp, zrd, zs, zsm
    INTEGER :: ilev, ilevp1, iplev, iplvp1, is, ism, ist, jk, jlev
    INTEGER, PARAMETER 	:: nlev=47

    !  Intrinsic functions 
    INTRINSIC EXP, LOG


vct(1) = 0._dp
vct(2) = 1.98918533325195_dp
vct(3) = 6.57208919525146_dp
vct(4) = 15.6739025115967_dp
vct(5) = 30.6242828369141_dp
vct(6) = 54.5457153320312_dp
vct(7) = 92.558837890625_dp
vct(8) = 150.504699707031_dp
vct(9) = 235.327453613281_dp
vct(10) = 356.100341796875_dp
vct(11) = 523.91943359375_dp
vct(12) = 751.04296875_dp
vct(13) = 1051.13720703125_dp
vct(14) = 1438.98852539062_dp
vct(15) = 1930.17724609375_dp
vct(16) = 2540.69702148438_dp
vct(17) = 3286.55297851562_dp
vct(18) = 4199.57421875_dp
vct(19) = 5303.95703125_dp
vct(20) = 6624.703125_dp
vct(21) = 8187.18359375_dp
vct(22) = 9976.13671875_dp
vct(23) = 11820.5390625_dp
vct(24) = 13431.390625_dp
vct(25) = 14736.359375_dp
vct(26) = 15689.2109375_dp
vct(27) = 16266.609375_dp
vct(28) = 16465._dp
vct(29) = 16297.62109375_dp
vct(30) = 15791.6015625_dp
vct(31) = 14985.26953125_dp
vct(32) = 13925.51953125_dp
vct(33) = 12665.2890625_dp
vct(34) = 11261.23046875_dp
vct(35) = 9771.40625_dp
vct(36) = 8253.2109375_dp
vct(37) = 6761.33984375_dp
vct(38) = 5345.9140625_dp
vct(39) = 4050.71801757812_dp
vct(40) = 2911.56909179688_dp
vct(41) = 1954.80493164062_dp
vct(42) = 1195.88989257812_dp
vct(43) = 638.14892578125_dp
vct(44) = 271.62646484375_dp
vct(45) = 72.0635986328125_dp
vct(46) = 0._dp
vct(47) = 0._dp
vct(48) = 0._dp
vct(49) = 0._dp
vct(50) = 0._dp
vct(51) = 0._dp
vct(52) = 0._dp
vct(53) = 0._dp
vct(54) = 0._dp
vct(55) = 0._dp
vct(56) = 0._dp
vct(57) = 0._dp
vct(58) = 0._dp
vct(59) = 0._dp
vct(60) = 0._dp
vct(61) = 0._dp
vct(62) = 0._dp
vct(63) = 0._dp
vct(64) = 0._dp
vct(65) = 0._dp
vct(66) = 0._dp
vct(67) = 0._dp
vct(68) = 0._dp
vct(69) = 0._dp
vct(70) = 0.000400000018998981_dp
vct(71) = 0.00289999996311963_dp
vct(72) = 0.00919999927282333_dp
vct(73) = 0.0203000009059906_dp
vct(74) = 0.0370000004768372_dp
vct(75) = 0.0595000013709068_dp
vct(76) = 0.0878999829292297_dp
vct(77) = 0.121999979019165_dp
vct(78) = 0.161400020122528_dp
vct(79) = 0.205699980258942_dp
vct(80) = 0.254199981689453_dp
vct(81) = 0.30620002746582_dp
vct(82) = 0.361100018024445_dp
vct(83) = 0.418200016021729_dp
vct(84) = 0.476700007915497_dp
vct(85) = 0.535899996757507_dp
vct(86) = 0.595099985599518_dp
vct(87) = 0.653599977493286_dp
vct(88) = 0.710600018501282_dp
vct(89) = 0.765399992465973_dp
vct(90) = 0.817200005054474_dp
vct(91) = 0.865000009536743_dp
vct(92) = 0.907700002193451_dp
vct(93) = 0.944199979305267_dp
vct(94) = 0.97299998998642_dp
vct(95) = 0.992299973964691_dp
vct(96) = 1._dp

nvclev = 48
nlevp1 = 48

vct_a = vct(1:nvclev)
vct_b = vct(nvclev+1:2*nvclev)


    !  Executable statements 

!LT
       ALLOCATE (ralpha(nlev))
       ALLOCATE (rlnpr(nlev))
       ALLOCATE (dela(nlev))
       ALLOCATE (delb(nlev))
       ALLOCATE (rddelb(nlev))
       ALLOCATE (cpg(nlev))
       ALLOCATE (delpr(nlev))
       ALLOCATE (rdelpr(nlev))
       ALLOCATE (ralphr(nlev))
       ALLOCATE (alpham(nlev))
       ALLOCATE (ardprc(nlev))
       ALLOCATE (rlnmar(nlev))
       ALLOCATE (aktlrd(nlev))
       ALLOCATE (altrcp(nlev))
       ALLOCATE (ceta(nlev))
!       ALLOCATE (cetah(nlevp1))
       ALLOCATE (bb(nlev,nlev))

    !-- 1. Initialize variables
    
    apzero = 101325.0_dp
    zrd = rd
    ralpha(1) = zrd*LOG(2.0_dp)
    rlnpr(1) = 2.0_dp*ralpha(1)
    ilev = nlev
    ilevp1 = ilev + 1
    nlevp1 = ilevp1
    nlevm1 = ilev - 1
    iplev = 0
    iplvp1 = 1
    is = nvclev + ilevp1
    ism = is - 1
    zpp = vct(1)
    zsm = vct(is)
    apsurf = 98550.0_dp

    t0icao = 288.0_dp
    tsticao = 216.5_dp
    zp0icao = 101320.0_dp
    rdlnp0i = rd*LOG(zp0icao)
    rdtstic = rd*tsticao
    alrrdic = 0.0065_dp/g
    rdt0ral = t0icao/alrrdic
    rdlnpti = rdlnp0i + (LOG(tsticao/t0icao))/alrrdic
    ptricao = EXP(rdlnpti/rd)
    gsticao = tsticao*(rdlnpti-1.0_dp/alrrdic)

    !-- 2. Calculate pressure-level values

10  CONTINUE

    zb = vct(nvclev+iplvp1+1)
    IF (zb > 0.0_dp) THEN
      nplev = iplev
      nplvp1 = iplvp1
      nplvp2 = iplvp1 + 1
      IF (iplev==0) THEN
        nplvpa = 2
      ELSE
        nplvpa = iplvp1
      END IF
      GO TO 20
    ELSE
      iplev = iplvp1
      iplvp1 = iplev + 1
      IF (iplvp1==ilevp1) GO TO 40
      zp = zpp
      zpp = vct(iplvp1)
      delpr(iplev) = zpp - zp
      rdelpr(iplev) = 1.0_dp/delpr(iplev)
      IF (iplev>1) THEN
        rlnpr(iplev) = zrd*LOG(zpp/zp)
        ralpha(iplev) = zrd - zp*rlnpr(iplev)/delpr(iplev)
      END IF
      alpham(iplev) = ralpha(iplev)*rcpd
      ardprc(iplev) = rlnpr(iplev)*rdelpr(iplev)*rcpd
      GO TO 10
    END IF

    !-- 3. Calculate sigma-level values

20  CONTINUE

    za = vct(ism-nvclev)
    IF (za > 0.0_dp) THEN
      nlmsgl = ism - nvclev
      nlmslp = nlmsgl + 1
      nlmsla = nlmslp
      GO TO 30
    ELSE
      is = ism
      ism = is - 1
      ist = is - nvclev
      zs = zsm
      zsm = vct(is)
      IF (ist==1) THEN
        nlmsgl = 0
        nlmslp = 1
        nlmsla = 2
        GO TO 30
      ELSE
        rlnpr(ist) = zrd*LOG(zs/zsm)
        ralpha(ist) = zrd - zsm*rlnpr(ist)/(zs-zsm)
      END IF
      GO TO 20
    END IF

    !-- 4. Calculate dela, delb, rddelb, cpg, and complete alphdb

30  CONTINUE

    DO jk = 1, nlev
      dela(jk) = vct(jk+1) - vct(jk)
      delb(jk) = vct(nvclev+jk+1) - vct(nvclev+jk)
      rddelb(jk) = rd*delb(jk)
      cpg(jk) = vct(nvclev+jk)*vct(jk+1) - vct(nvclev+jk+1)*vct(jk)
    END DO

    DO jk = nlmslp, nlev
      alpham(jk) = ralpha(jk)*delb(jk)
    END DO

    !-- 5. Compute full level values of the hybrid coordinate

40  CONTINUE

!LT
    zetam = vct(1)/apzero + vct(nvclev+1)
    cetah(1) = zetam

    DO jlev = 1, nlev
      zetap = vct(jlev+1)/apzero + vct(nvclev+1+jlev)
      ceta(jlev) = (zetam+zetap)*0.5_dp
      cetah(jlev+1) = zetap
      zetam = zetap
    END DO
    
  END SUBROUTINE inihyb

END MODULE mo_hyb
