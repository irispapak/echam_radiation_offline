#ifdef __xlC__
@PROCESS HOT
#else
#define SWDIV_NOCHK(a,b) ((a)/(b))
#endif

MODULE MO_ECHAM5_SW

  USE MO_KIND  , ONLY : DP
  USE MO_CONSTANTS,    ONLY : RG=>G,RD, NDAYLEN=>DAYLEN ! LT
!LT
!  USE MO_TIME_CONTROL, ONLY : NDAYLEN

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: nsw, sw, su_sw4

  SAVE

  ! PARAMETERS:
  !
  !  nsw      : number of bands of SW scheme
  !
  !  novlp    : index for cloud overlap assumption in radiation computation
  !             1 : maximum-random overlap
  !             2 : maximum overlap
  !             3 : random overlap
  !
  !        - C. Cagnazzo, June 2005:
  !        - from 4 to 6 bands, following 26Rcycle Morcrette Code

  INTEGER, PARAMETER :: NSW=6
  INTEGER, PARAMETER :: NOVLP=1


  !     ------------------------------------------------------------------
  !*    ** *MO_SW* - COEFFICIENTS FOR SHORTWAVE RADIATION TRANSFER
  !     ------------------------------------------------------------------

  REAL(DP):: APAD(6,3,7)
  REAL(DP):: BPAD(6,3,7)
  REAL(DP):: RRAY(6,6)
  REAL(DP):: RSUN(6)
  REAL(DP):: RPDH1
  REAL(DP):: RPDU1
  REAL(DP):: RPNH
  REAL(DP):: RPNU
  REAL(DP):: RSWCE(6)
  REAL(DP):: RSWCP(6)
  REAL(DP):: RTDH2O
  REAL(DP):: RTDUMG
  REAL(DP):: RTH2O
  REAL(DP):: RTUMG
  REAL(DP):: D(6,3)
  REAL(DP):: REXPO3(6,2,7)
  INTEGER :: NEXPO3(6)

  REAL(DP):: RYFWCA(6)
  REAL(DP):: RYFWCB(6)
  REAL(DP):: RYFWCC(6)
  REAL(DP):: RYFWCD(6)
  REAL(DP):: RYFWCE(6)
  REAL(DP):: RYFWCF(6)

  REAL(DP):: REBCUA(6)
  REAL(DP):: REBCUB(6)
  REAL(DP):: REBCUC(6)
  REAL(DP):: REBCUD(6)
  REAL(DP):: REBCUE(6)
  REAL(DP):: REBCUF(6)
  REAL(DP):: REBCUG(16)
  REAL(DP):: REBCUH(16)
  REAL(DP):: REBCUI(6)
  REAL(DP):: REBCUJ(6)

  REAL(DP):: RASWCA(6)
  REAL(DP):: RASWCB(6)
  REAL(DP):: RASWCC(6)
  REAL(DP):: RASWCD(6)
  REAL(DP):: RASWCE(6)
  REAL(DP):: RASWCF(6)

  REAL(DP):: RFULIO(16,3)
  REAL(DP):: RFLAA0(6)
  REAL(DP):: RFLAA1(6)
  REAL(DP):: RFLBB0(6)
  REAL(DP):: RFLBB1(6)
  REAL(DP):: RFLBB2(6)
  REAL(DP):: RFLBB3(6)
  REAL(DP):: RFLCC0(6)
  REAL(DP):: RFLCC1(6)
  REAL(DP):: RFLCC2(6)
  REAL(DP):: RFLCC3(6)
  REAL(DP):: RFLDD0(6)
  REAL(DP):: RFLDD1(6)
  REAL(DP):: RFLDD2(6)
  REAL(DP):: RFLDD3(6)

  REAL(DP):: RHSAVI(16,3)

  REAL(DP):: RSUSHE(6)
  REAL(DP):: RSUSHF(6)
  REAL(DP):: RSUSHH(6)
  REAL(DP):: RSUSHK(6)
  REAL(DP):: RSUSHA(6)
  REAL(DP):: RSUSHG(6)
  REAL(DP):: RSUSHFA(4)
  REAL(DP):: RSUSHC
  REAL(DP):: RSUSHD

  REAL(DP):: REFFIA
  REAL(DP):: REFFIB
  REAL(DP):: RTIW
  REAL(DP):: RRIW
  REAL(DP):: RROMA(6)
  REAL(DP):: RROMB(6)
  REAL(DP):: RRASY(6)

  REAL(DP):: RHSRA(6)
  REAL(DP):: RHSRB(6)
  REAL(DP):: RHSRC(6)
  REAL(DP):: RHSRD(6)
  REAL(DP):: RHSRE(6)
  REAL(DP):: RHSRF(6)
  REAL(DP):: RHSRTA
  REAL(DP):: RHSRTB

  REAL(DP):: RTAUA(6,6)
  REAL(DP):: RPIZA(6,6)
  REAL(DP):: RCGA(6,6)
  REAL(DP):: RAER(6,6)
  REAL(DP):: RSNOALB(6)
  REAL(DP):: RSNOMEL(6)
  REAL(DP):: RWEIGS(6)
  REAL(DP):: RWEIGV(6)


  !        * E.C.M.W.F. PHYSICS PACKAGE *

  !     J.-J. MORCRETTE       E.C.M.W.F.      89/07/14

  !  NAME     TYPE     PURPOSE
  !  ----  :  ----   : ---------------------------------------------------
  !  APAD  :  REAL     PADE APPROXIMANTS NUMERATOR
  !  BPAD  :  REAL     PADE APPROXIMANTS DENOMINATOR
  !  D     :  REAL     TRANSMISSION LIMIT FOR INFINITE ABSORBER AMOUNT
  !  RRAY  :  REAL     RAYLEIGH SCATTERING COEFFICIENTS
  !  RSUN  :  REAL     SOLAR FRACTION IN SPECTRAL INTERVALS
  !  RPDH1 :  1 + EXPONENT PRESSURE DEPENDENCE H2O
  !  RPDU1 :  1 + EXPONENT PRESSURE DEPENDENCE UNIFORMLY MIXED GASES
  !  RPNH  :  REFERENCE PRESSURE FACTOR FOR H2O
  !  RPNU  :  REFERENCE PRESSURE FACTOR FOR UNIFORMLY MIXED GASES
  !  RSWCE :  E-TYPE, H2O CONTINUUM ABSORPTION COEFFICIENT 
  !  RSWCP :  P-TYPE, H2O CONTINUUM ABSORPTION COEFFICIENT 
  !  RTDH2O:  EXPONENT TEMPERATURE DEPENDENCE H2O
  !  RTDUMG:  EXPONENT TEMPERATURE DEPENDENCE UNIFORMLY MIXED GASES
  !  RTH2O :  REFERENCE TEMPERATURE H2O
  !  RTUMG :  REFERENCE TEMPERATURE UNIFORMLY MIXED GASES
  !     -----------------------------------------------------------------

  !        * E.C.M.W.F. PHYSICS PACKAGE *

  !     J.-J. MORCRETTE       E.C.M.W.F.      89/07/14

  !  NAME     TYPE     PURPOSE
  !  ----  :  ----   : ---------------------------------------------------
  !*    FOUQUART (1987) WATER CLOUD OPTICAL PROPERTIES

  ! RYFWCA :  REAL   : C1 IN OPTICAL THICKNESS FORMULA
  ! RYFWCB :  REAL   : C2 IN OPTICAL THICKNESS FORMULA
  ! RYFWCC :  REAL   : SINGLE SCATTERING ALBEDO PARAMETER
  ! RYFWCD :  REAL   : SINGLE SCATTERING ALBEDO PARAMETER
  ! RYFWCE :  REAL   : SINGLE SCATTERING ALBEDO PARAMETER
  ! RYFWCF :  REAL   : ASSYMETRY FACTOR

  !*    SLINGO (1989) WATER CLOUD OPTICAL PROPERTIES

  ! RASWCA :  REAL   : C1 IN OPTICAL THICKNESS FORMULA
  ! RASWCB :  REAL   : C2 IN OPTICAL THICKNESS FORMULA
  ! RASWCC :  REAL   : SINGLE SCATTERING ALBEDO PARAMETER
  ! RASWCD :  REAL   : SINGLE SCATTERING ALBEDO PARAMETER
  ! RASWCE :  REAL   : SINGLE SCATTERING ALBEDO PARAMETER
  ! RASWCF :  REAL   : ASSYMETRY FACTOR

  !*   SAVIJARVI (1998) WATER CLOUD OPTICAL PROPERTIES (RRTM)

  ! RHSAVI : REAL    : MASS ABSORPTION COEFFICIENTS (POLYNOMIAL DEVELOPM)

  !*    ICE CLOUD OPTICAL PROPERTIES DERIVED FROM EBERT-CURRY (1992)

  ! REBCUA :  REAL   : C1 IN OPTICAL THICKNESS FORMULA
  ! REBCUB :  REAL   : C2 IN OPTICAL THICKNESS FORMULA
  ! REBCUC :  REAL   : 1-C3  IN SINGLE SCATTERING ALBEDO FORMULA
  ! REBCUD :  REAL   : C4 IN SINGLE SCATTERING ALBEDO FORMULA
  ! REBCUE :  REAL   : C5 IN ASSYMETRY FACTOR FORMULA
  ! REBCUF :  REAL   : C6 IN ASSYMETRY FACTOR FORMULA
  ! REBCUG :  REAL   : C7 IN MASS ABSORPTION COEFFICIENT FORMULA
  ! REBCUH :  REAL   : C8 IN MASS ABSORPTION COEFFICIENT FORMULA
  ! REBCUI :  REAL   : C7 IN MASS ABSORPTION COEFFICIENT SPECTRAL FORMULA
  ! REBCUJ :  REAL   : C8 IN MASS ABSORPTION COEFFICIENT SPECTRAL FORMULA

  !*    ICE CLOUD OPTICAL PROPERTIES DERIVED FROM SUN-SHINE (1995)

  ! RSHSUE :  REAL   : E IN SINGLE SCATTERING ALBEDO FORMULA
  ! RSHSUF :  REAL   : F IN SINGLE SCATTERING ALBEDO FORMULA
  ! RSHSUH :  REAL   : H IN ASSYMETRY FACTOR FORMULA
  ! RSHSUK :  REAL   : K IN ASSYMETRY FACTOR FORMULA
  ! RSHSUA :  REAL   : ALPHA IN SSA CORRECTION FACTOR FORMULA
  ! RSHSUG :  REAL   : GAMMA IN ASSYMETRY CORRECTION FACTOR FORMULA
  ! RSHSUFA:  REAL   : COEFFICIENTS IN TEMPERATURE CORRECTION FACTOR

  ! REFFIA :  REAL   : C9  IN EFFECTIVE RADIUS FORMULA
  ! REFFIB :  REAL   : C10 IN EFFECTIVE RADIUS FORMULA

  !*    ICE CLOUD OPTICAL PROPERTIES DERIVED FROM FU-LIOU (1993)

  ! RFULIO :  REAL   : COEFFICIENTS IN EXPRESSION FOR LW EXTINCTION COEFF.
  ! RFLAA  :  REAL   : COEFFICIENTS IN EXPRESSION FOR SW EXTINCTION COEFF.
  ! RFLBB  :  REAL   : COEFFICIENTS IN EXPRESSION FOR SW SINGLE SCATT.ALB.
  ! RFLCC  :  REAL   : COEFFICIENTS IN EXPRESSION FOR SW ASSYMETRY FACTOR
  ! RFLDD  :  REAL   : COEFFICIENTS IN EXPRESSION FOR SW ASSYMETRY FACTOR

  !*    TRANSITION BETWEEN LIQUID AND SOLID WATER

  ! RTIW   :  REAL   : TEMPERATURE THRESHOLD
  ! RRIW   :  REAL   : TRANSITION RANGE

  !*    RAIN OPTICAL PROPERTIES FROM SAVIJARVI (1996)

  ! RROMA  :  REAL   : COEFFICIENTS FOR SINGLE SCATTERING ALBEDO
  ! RROMB  :  REAL   : COEFFICIENTS FOR SINGLE SCATTERING ALBEDO
  ! RRASY  :  REAL   : COEFFICIENTS FOR ASSYMETRY FACTOR
  ! RHSRA  :  REAL   : COEFFICIENTS FOR OPTICAL THICKNESS
  ! RHSRB  :  REAL   : COEFFICIENTS FOR OPTICAL THICKNESS
  ! RHSRC  :  REAL   : COEFFICIENTS FOR SINGLE SCATTERING ALBEDO
  ! RHSRD  :  REAL   : COEFFICIENTS FOR SINGLE SCATTERING ALBEDO
  ! RHSRE  :  REAL   : COEFFICIENTS FOR ASSYMETRY FACTOR 
  ! RHSRF  :  REAL   : COEFFICIENTS FOR ASSYMETRY FACTOR
  ! RHSRTA :  REAL   : COEFFICIENTS FOR OPTICAL THICKNESS
  ! RHSRTB :  REAL   : COEFFICIENTS FOR OPTICAL THICKNESS
  !     -----------------------------------------------------------------

  !        * E.C.M.W.F. PHYSICS PACKAGE *

  !     J.-J. MORCRETTE       E.C.M.W.F.      89/07/14

  !  NAME     TYPE     PURPOSE
  !  ----  :  ----   : -------
  !  RTAUA :  REAL     S.W. NORMALIZED OPTICAL THICKNESS AT 0.55 MICRON
  !  RPIZA :  REAL     S.W. SINGLE SCATTERING ALBEDO
  !  RCGA  :  REAL     S.W. ASSYMETRY FACTOR
  !  RAER  :  REAL     L.W. ABSORPTION COEFFICIENTS
  !     -----------------------------------------------------------------

  !        * E.C.M.W.F. PHYSICS PACKAGE *

  !     J.-J. MORCRETTE       E.C.M.W.F.      89/07/14

  !  NAME     TYPE     PURPOSE
  !  ----  :  ----   : -------
  ! RSNOALB:  REAL     S.W. SPECTRAL ALBEDO (Fresh Snow) after WARREN
  ! RSNOMEL:  REAL     S.W. SPECTRAL ALBEDO (Aging Snow) after WARREN

  ! RWEIGS :  REAL     S.W. SPECTR WEIGHT for soil (Briegleb, Ramanathan)
  ! RWEIGV :  REAL     S.W. SPECTR WEIGHT for vegetation (BR86)
  !     -----------------------------------------------------------------
CONTAINS 
  SUBROUTINE SW &
       &( KIDIA, KFDIA , KBDIM  , KLEV &
       &, PSCT , PRMU0, PPSOL, PPMB, PTAVE &
       &, PALB_VIS , PALB_NIR , PCARDI, PWV , PQS  &
       &, PCLEAR, PCLDFRAC, PDP, POZ, PCG, POMEGA, PTAU  &
       &, PSWEXT, PSWSSA, PSWASY &
       &, PHEAT, PFDOWN, PFUP, PCEAT, PCDOWN, PCUP &
       &, PFDNN, PFDNV , PFUPN, PFUPV, PCDNN, PCDNV , PCUPN, PCUPV &
       &, PSUDU, P_NIR_NET_SURFACE, P_NIR_FRACT_DIFFUSE &
       &, P_VIS_NET_SURFACE, P_VIS_FRACT_DIFFUSE &
       &)

    !**** *SW* - COMPUTES THE SHORTWAVE RADIATION FLUXES.

    !     PURPOSE.
    !     --------

    !          THIS ROUTINE COMPUTES THE SHORTWAVE RADIATION FLUXES IN TWO
    !     SPECTRAL INTERVALS FOLLOWING FOUQUART AND BONNEL (1980).

    !**   INTERFACE.
    !     ----------
    !          *SW* IS CALLED FROM *RADLSW*

    !        IMPLICIT ARGUMENTS :
    !        --------------------

    !     ==== INPUTS ===
    !     ==== OUTPUTS ===

    !     METHOD.
    !     -------

    !          1. COMPUTES ABSORBER AMOUNTS                 (SWU)
    !          2. COMPUTES FLUXES IN U.V./VISIBLE  SPECTRAL INTERVAL (SW1S)
    !          3. COMPUTES FLUXES IN NEAR-INFRARED SPECTRAL INTERVAL (SWNI)

    !     EXTERNALS.
    !     ----------

    !          *SWU*, *SW1S*, *SWNI*

    !     REFERENCE.
    !     ----------

    !        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
    !        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

    !     AUTHOR.
    !     -------
    !        JEAN-JACQUES MORCRETTE  *ECMWF*

    !     MODIFICATIONS.
    !     --------------
    !        ORIGINAL : 89-07-14
    !        95-01-01   J.-J. MORCRETTE  Direct/Diffuse Albedo
    !        95-12-07   J.-J. MORCRETTE  Near-Infrared in nsw-1 Intervals
    !        990128     JJMorcrette      sunshine duration
    !        99-05-25   JJMorcrette      Revised aerosols
    !        C. Cagnazzo, June 2005: 
    !                from 4 to 6 bands, following 26R and 28Rcycle Morcrette Code
    !     ------------------------------------------------------------------


    REAL(DP) :: RCDAY

    !     DUMMY INTEGER SCALARS
    INTEGER :: KFDIA
    INTEGER :: KIDIA
    INTEGER :: KLEV
    INTEGER :: KBDIM

    !     DUMMY REAL SCALARS
    REAL(DP) :: PSCT

    !     ------------------------------------------------------------------
    !*       0.1   ARGUMENTS
    !              ---------

    REAL(DP):: PPSOL(KBDIM), PRMU0(KBDIM)&
         &,  PWV(KBDIM,KLEV),PQS(KBDIM,KLEV)

    REAL(DP):: PALB_VIS(KBDIM) , PALB_NIR(KBDIM)                        &
         &,  PCG(KBDIM,NSW,KLEV)   , PCLEAR(KBDIM), PCLDFRAC(KBDIM,KLEV)&
         &,  PDP(KBDIM,KLEV), PCARDI(KBDIM,KLEV)                        &
         &,  POMEGA(KBDIM,NSW,KLEV), POZ(KBDIM,KLEV)                    &
         &,  PPMB(KBDIM,KLEV+1)                                         &
         &,  PTAU(KBDIM,NSW,KLEV)  , PTAVE(KBDIM,KLEV)

    REAL(DP) :: PHEAT(KBDIM,KLEV), PFDOWN(KBDIM,KLEV+1), PFUP(KBDIM,KLEV+1),&
         &PFUPV(KBDIM), PFUPN(KBDIM), PFDNV(KBDIM), PFDNN(KBDIM)&
         &,  PCEAT(KBDIM,KLEV), PCDOWN(KBDIM,KLEV+1), PCUP(KBDIM,KLEV+1)&
         &,  PCUPV(KBDIM), PCUPN(KBDIM), PCDNV(KBDIM), PCDNN(KBDIM)&
         &,  PSUDU(KBDIM)

    REAL(DP) :: PSWEXT(KBDIM,KLEV,NSW), PSWSSA(KBDIM,KLEV,NSW), &
         &         PSWASY(KBDIM,KLEV,NSW)

    REAL(DP) ::  P_NIR_NET_SURFACE(KBDIM), P_NIR_FRACT_DIFFUSE(KBDIM)
    REAL(DP) ::  P_VIS_NET_SURFACE(KBDIM), P_VIS_FRACT_DIFFUSE(KBDIM)

    !     ------------------------------------------------------------------

    !*       0.2   LOCAL ARRAYS
    !              ------------

    REAL(DP) :: ZAKI(KBDIM,2,NSW)                   & 
         &,  ZDSIG(KBDIM,KLEV)   , ZFACT(KBDIM)     &
         &,  ZFD(KBDIM,KLEV+1)   , ZCD(KBDIM,KLEV+1)&
         &,  ZCDOWN(KBDIM,KLEV+1), ZCDNIR(KBDIM,KLEV+1), ZCDUVS(KBDIM,KLEV+1) &
         &,  ZFDOWN(KBDIM,KLEV+1), ZFDNIR(KBDIM,KLEV+1), ZFDUVS(KBDIM,KLEV+1) &
         &,  ZFU(KBDIM,KLEV+1)   , ZCU(KBDIM,KLEV+1)&
         &,  ZCUP(KBDIM,KLEV+1)  , ZCUNIR(KBDIM,KLEV+1),  ZCUUVS(KBDIM,KLEV+1)&
         &,  ZFUP(KBDIM,KLEV+1)  , ZFUNIR(KBDIM,KLEV+1),  ZFUUVS(KBDIM,KLEV+1)&
         &,  ZRMU(KBDIM)         , ZSEC(KBDIM)      &
         &,  ZSUDU1(KBDIM)       , ZSUDU2(KBDIM)    &
         &,  ZSUDU2T(KBDIM),  ZSUDU1T(KBDIM)        &
         &,  ZUD(KBDIM,5,KLEV+1)

    REAL(DP) :: P_VIS(KBDIM), P_VIS_DIFFUSE(KBDIM), P_VIS_DIFFUSE_SUM(KBDIM)
    REAL(DP) :: P_NIR(KBDIM), P_NIR_DIFFUSE(KBDIM), P_NIR_DIFFUSE_SUM(KBDIM)

    !     LOCAL INTEGER SCALARS
    INTEGER :: JK, JKL, JL, JNU,INUVS, INUIR

    !     LOCAL REAL SCALARS
    REAL(DP):: ZDCNET, ZDFNET

    !     ------------------------------------------------------------------
    !*         1.     ABSORBER AMOUNTS AND OTHER USEFUL QUANTITIES
    !                 --------------------------------------------

    RCDAY = REAL(NDAYLEN,DP)*RG/3.5_dp/RD

    CALL SWU ( KIDIA,KFDIA ,KBDIM ,KLEV &
         &, PSCT ,PCARDI,PPMB ,PPSOL &
         &, PRMU0,PTAVE ,PWV &
         &, ZAKI ,ZDSIG,ZFACT,ZRMU,ZSEC,ZUD )

    !     ------------------------------------------------------------------
    !*         2.     INTERVAL (0.185/0.25-0.68 MICRON): U.V. AND VISIBLE
    !                 ---------------------------------------------------
    INUVS=1
    INUIR=4

    DO JK = 1 , KLEV+1
      DO JL = KIDIA,KFDIA
        ZFD(JL,JK) =0._DP
        ZFU(JL,JK) =0._DP
        ZCD(JL,JK) =0._DP
        ZCU(JL,JK) =0._DP
      ENDDO
    ENDDO
    DO JL = KIDIA,KFDIA
      ZSUDU1T(JL)=0._DP
      P_VIS_NET_SURFACE(JL) = 0._DP
      P_VIS_DIFFUSE_SUM(JL) = 0._DP
    ENDDO

    DO JNU = INUVS , INUIR-1
      CALL SW1S &
           &( KIDIA, KFDIA, KBDIM, KLEV, JNU, PSWEXT, PSWSSA, PSWASY &
           &,  PALB_VIS , PCG  , PCLDFRAC , PCLEAR                         &
           &,  ZDSIG, POMEGA, POZ  , ZRMU , ZSEC , PTAU  , ZUD             &
           &,  ZFDUVS,ZFUUVS, ZCDUVS,ZCUUVS, ZSUDU1                        &
           &,  P_VIS, P_VIS_DIFFUSE )

      DO JK = 1 , KLEV+1
        DO JL = KIDIA,KFDIA
          ZFD(JL,JK)=ZFD(JL,JK)+ZFDUVS(JL,JK)
          ZFU(JL,JK)=ZFU(JL,JK)+ZFUUVS(JL,JK)
          ZCD(JL,JK)=ZCD(JL,JK)+ZCDUVS(JL,JK)
          ZCU(JL,JK)=ZCU(JL,JK)+ZCUUVS(JL,JK)
        ENDDO
      ENDDO
      DO JL = KIDIA,KFDIA
        ZSUDU1T(JL)=ZSUDU1T(JL)+ZSUDU1(JL)
        P_VIS_NET_SURFACE(JL) = P_VIS_NET_SURFACE(JL) + P_VIS(JL)
        P_VIS_DIFFUSE_SUM(JL) = P_VIS_DIFFUSE_SUM(JL) + P_VIS(JL) &
             * P_VIS_DIFFUSE(JL)
      ENDDO
    ENDDO

    DO JL = KIDIA,KFDIA
      P_VIS_FRACT_DIFFUSE(JL) = P_VIS_DIFFUSE_SUM(JL) /           &
           (P_VIS_NET_SURFACE(JL) + EPSILON(1._dp))
    ENDDO

    !     ------------------------------------------------------------------
    !*         3.     INTERVAL (0.68-4.00 MICRON): NEAR-INFRARED
    !                 ------------------------------------------
    DO JK = 1 , KLEV+1
      DO JL = KIDIA,KFDIA
        ZFDOWN(JL,JK) = 0._DP
        ZFUP  (JL,JK) = 0._DP
        ZCDOWN(JL,JK) = 0._DP
        ZCUP  (JL,JK) = 0._DP
      ENDDO
    ENDDO
    DO JL = KIDIA,KFDIA
      ZSUDU2T(JL) = 0._DP
      P_NIR_NET_SURFACE(JL) = 0._DP
      P_NIR_DIFFUSE_SUM(JL) = 0._DP
    ENDDO

    DO JNU = INUIR , NSW
      CALL SWNI &
           &(  KIDIA, KFDIA, KBDIM, KLEV, JNU, PSWEXT, PSWSSA, PSWASY &
           &,  ZAKI  , PALB_NIR , PCG  , PCLDFRAC, PCLEAR                   &
           &,  ZDSIG ,POMEGA, POZ  , ZRMU , ZSEC , PTAU    , ZUD            &
           &,  PWV   ,PQS                                                   &
           &,  ZFDNIR,ZFUNIR,ZCDNIR,ZCUNIR,ZSUDU2                           &
           &,  P_NIR, P_NIR_DIFFUSE )

      DO JK = 1 , KLEV+1
        DO JL = KIDIA,KFDIA
          ZFDOWN(JL,JK)=ZFDOWN(JL,JK)+ZFDNIR(JL,JK)
          ZFUP  (JL,JK)=ZFUP  (JL,JK)+ZFUNIR(JL,JK)
          ZCDOWN(JL,JK)=ZCDOWN(JL,JK)+ZCDNIR(JL,JK)
          ZCUP  (JL,JK)=ZCUP  (JL,JK)+ZCUNIR(JL,JK)
        ENDDO
      ENDDO
      DO JL = KIDIA,KFDIA
        ZSUDU2T(JL) = ZSUDU2T(JL) + ZSUDU2(JL)
        P_NIR_NET_SURFACE(JL) = P_NIR_NET_SURFACE(JL) + P_NIR(JL)
        P_NIR_DIFFUSE_SUM(JL) = P_NIR_DIFFUSE_SUM(JL) + P_NIR(JL) &
             * P_NIR_DIFFUSE(JL)
      ENDDO
    ENDDO
    DO JL = KIDIA,KFDIA
      P_NIR_FRACT_DIFFUSE(JL) = P_NIR_DIFFUSE_SUM(JL)             &
           / (P_NIR_NET_SURFACE(JL) + EPSILON(1._dp))
    ENDDO

    !     ------------------------------------------------------------------
    !*         4.     FILL THE DIAGNOSTIC ARRAYS
    !                 --------------------------
    DO JL = KIDIA,KFDIA
      PFDNN(JL)=ZFDOWN(JL,1)*ZFACT(JL)
      PFDNV(JL)=ZFD(JL,1)*ZFACT(JL)
      PFUPN(JL)=ZFUP(JL,KLEV+1)*ZFACT(JL)
      PFUPV(JL)=ZFU(JL,KLEV+1)*ZFACT(JL)

      PCDNN(JL)=ZCDOWN(JL,1)*ZFACT(JL)
      PCDNV(JL)=ZCD(JL,1)*ZFACT(JL)
      PCUPN(JL)=ZCUP(JL,KLEV+1)*ZFACT(JL)
      PCUPV(JL)=ZCU(JL,KLEV+1)*ZFACT(JL)

      PSUDU(JL)=(ZSUDU1T(JL)+ZSUDU2T(JL))*ZFACT(JL)

      P_VIS_NET_SURFACE(JL) = P_VIS_NET_SURFACE(JL) * ZFACT(JL)
      P_NIR_NET_SURFACE(JL) = P_NIR_NET_SURFACE(JL) * ZFACT(JL)
    ENDDO

    DO JK = 1 , KLEV+1
      DO JL = KIDIA,KFDIA
        PFUP(JL,JK)   = (ZFUP(JL,JK)   + ZFU(JL,JK)) * ZFACT(JL)
        PFDOWN(JL,JK) = (ZFDOWN(JL,JK) + ZFD(JL,JK)) * ZFACT(JL)
        PCUP(JL,JK)   = (ZCUP(JL,JK)   + ZCU(JL,JK)) * ZFACT(JL)
        PCDOWN(JL,JK) = (ZCDOWN(JL,JK) + ZCD(JL,JK)) * ZFACT(JL)
      ENDDO
    ENDDO

    DO JKL = 1 , KLEV
      JK = KLEV+1 - JKL
      DO JL = KIDIA,KFDIA
        ZDFNET = PFUP(JL,JK+1) - PFDOWN(JL,JK+1)-PFUP(JL,JK  ) + PFDOWN(JL,JK  )
        PHEAT(JL,JK) = RCDAY * ZDFNET / PDP(JL,JKL)
        ZDCNET = PCUP(JL,JK+1) - PCDOWN(JL,JK+1)-PCUP(JL,JK  ) + PCDOWN(JL,JK  )
        PCEAT(JL,JK) = RCDAY * ZDCNET / PDP(JL,JKL)
      ENDDO
    ENDDO

  END SUBROUTINE SW

  SUBROUTINE SW1S &
       &( KIDIA , KFDIA , KBDIM , KLEV , KNU &
       &, PSWEXT, PSWSSA, PSWASY      & 
       &, PALB, PCG  , PCLD , PCLEAR &
       &, PDSIG , POMEGA, POZ  , PRMU , PSEC , PTAU  , PUD  &
       &, PFD   , PFU   , PCD  , PCU  , PSUDU1 &
       &, P_VIS , P_VIS_FRACT_DIFFUSE &
       &)

    !**** *SW1S* - SHORTWAVE RADIATION, FIRST SPECTRAL INTERVAL

    !     PURPOSE.
    !     --------

    !          THIS ROUTINE COMPUTES THE SHORTWAVE RADIATION FLUXES IN TWO
    !     SPECTRAL INTERVALS FOLLOWING FOUQUART AND BONNEL (1980).

    !**   INTERFACE.
    !     ----------

    !          *SW1S* IS CALLED FROM *SW*.


    !        IMPLICIT ARGUMENTS :
    !        --------------------

    !     ==== INPUTS ===
    !     ==== OUTPUTS ===

    !     METHOD.
    !     -------

    !          1. COMPUTES QUANTITIES FOR THE CLEAR-SKY FRACTION OF THE
    !     COLUMN
    !          2. COMPUTES UPWARD AND DOWNWARD FLUXES CORRESPONDING TO
    !     CONTINUUM SCATTERING
    !          3. MULTIPLY BY OZONE TRANSMISSION FUNCTION

    !     EXTERNALS.
    !     ----------

    !          *SWCLR*, *SWR*, *SWTT*, *SWUVO3*

    !     REFERENCE.
    !     ----------

    !        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
    !        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

    !     AUTHOR.
    !     -------
    !        JEAN-JACQUES MORCRETTE  *ECMWF*

    !     MODIFICATIONS.
    !     --------------
    !        ORIGINAL : 89-07-14
    !        94-11-15   J.-J. MORCRETTE    DIRECT/DIFFUSE ALBEDO
    !        96-01-15   J.-J. MORCRETTE    SW in nsw SPECTRAL INTERVALS 
    !        990128     JJMorcrette        sunshine duration
    !        99-05-25   JJMorcrette        Revised aerosols
    !        June 2005:  CCagnazzo 6 spectral bands from Morcrette code CY26R3

    !     ------------------------------------------------------------------



    !     DUMMY INTEGER SCALARS
    INTEGER :: KFDIA
    INTEGER :: KIDIA
    INTEGER :: KLEV
    INTEGER :: KBDIM
    INTEGER :: KNU



    !     ------------------------------------------------------------------

    !*       0.1   ARGUMENTS
    !              ---------

    REAL(DP):: &
         &   PALB(KBDIM) &
         &,  PCG(KBDIM,NSW,KLEV)   , PCLD(KBDIM,KLEV) &
         &,  PCLEAR(KBDIM)&
         &,  PDSIG(KBDIM,KLEV)&
         &,  POMEGA(KBDIM,NSW,KLEV), POZ(KBDIM,KLEV)&
         &,  PRMU(KBDIM)           , PSEC(KBDIM)&
         &,  PTAU(KBDIM,NSW,KLEV)  , PUD(KBDIM,5,KLEV+1)


    REAL(DP):: PFD(KBDIM,KLEV+1)     , PFU(KBDIM,KLEV+1)&
         &,  PCD(KBDIM,KLEV+1)     , PCU(KBDIM,KLEV+1)&
         &,  PSUDU1(KBDIM)

    REAL(DP) :: PSWEXT(KBDIM,KLEV,NSW), PSWSSA(KBDIM,KLEV,NSW) &
         &, PSWASY(KBDIM,KLEV,NSW)

    REAL(DP)::  P_VIS(KBDIM), P_VIS_FRACT_DIFFUSE(KBDIM)
    !     ------------------------------------------------------------------
    !*       0.2   LOCAL ARRAYS
    !              ------------

    INTEGER :: IIND(4)

    REAL(DP):: ZCGAZ(KBDIM,KLEV)&
         &,  ZDIFF(KBDIM)        , ZDIRF(KBDIM)        &
         &,  ZDIFT(KBDIM)        , ZDIRT(KBDIM)        &
         &,  ZPIZAZ(KBDIM,KLEV)&
         &,  ZRAYL(KBDIM), ZRAY1(KBDIM,KLEV+1), ZRAY2(KBDIM,KLEV+1)&
         &,  ZREFZ(KBDIM,2,KLEV+1)&
         &,  ZRJ(KBDIM,6,KLEV+1), ZRJ0(KBDIM,6,KLEV+1)&
         &,  ZRK(KBDIM,6,KLEV+1), ZRK0(KBDIM,6,KLEV+1)&
         &,  ZRMUE(KBDIM,KLEV+1), ZRMU0(KBDIM,KLEV+1)&
         &,  ZR(KBDIM,4)&			    ! new
         &,  ZTAUAZ(KBDIM,KLEV)&
         &,  ZTRA1(KBDIM,KLEV+1), ZTRA2(KBDIM,KLEV+1)&
         &,  ZTRCLD(KBDIM)      , ZTRCLR(KBDIM)&
         &,  ZW(KBDIM,4),  ZO(KBDIM,2) ,ZT(KBDIM,2) & ! new
         &,  Z_FRACT_DIFFUSE_CLEAR(KBDIM), Z_FRACT_DIFFUSE_CLOUDY(KBDIM)

    REAL(DP)::ZTAUAZ_N(KBDIM,KLEV),ZPIZAZ_N(KBDIM,KLEV),ZCGAZ_N(KBDIM,KLEV)

    INTEGER :: IKL, IKM1, JAJ, JK, JL


    !     ------------------------------------------------------------------
    !*         1.     FIRST SPECTRAL INTERVAL (0.25-0.68 MICRON)
    !                 ----------------------- ------------------


    !*         1.1    OPTICAL THICKNESS FOR RAYLEIGH SCATTERING
    !                 -----------------------------------------


    DO JL = KIDIA,KFDIA
      ZRAYL(JL) =  RRAY(KNU,1) + PRMU(JL) * (RRAY(KNU,2) + PRMU(JL)&
           &* (RRAY(KNU,3) + PRMU(JL) * (RRAY(KNU,4) + PRMU(JL)&
           &* (RRAY(KNU,5) + PRMU(JL) *  RRAY(KNU,6)       ))))
    ENDDO


    !     ------------------------------------------------------------------
    !*         2.    CONTINUUM SCATTERING CALCULATIONS
    !                ---------------------------------


    !*         2.1   CLEAR-SKY FRACTION OF THE COLUMN
    !                --------------------------------


    CALL SWCLR &
         &( KIDIA  , KFDIA , KBDIM  , KLEV , KNU &
         &, PSWEXT, PSWSSA, PSWASY &
         &, PALB , PDSIG , ZRAYL, PSEC &
         &, ZTAUAZ, ZCGAZ  , ZPIZAZ, ZRAY1 , ZRAY2, ZREFZ, ZRJ0 &
         &, ZRK0   , ZRMU0 , ZTRA1, ZTRA2, ZTRCLR &
         &, Z_FRACT_DIFFUSE_CLEAR &
         &, ZTAUAZ_N, ZPIZAZ_N, ZCGAZ_N)


    !*         2.2   CLOUDY FRACTION OF THE COLUMN
    !                -----------------------------

    CALL SWR &
         &( KIDIA ,KFDIA ,KBDIM  ,KLEV  , KNU &
         &, PALB ,PCG   ,PCLD  ,POMEGA, PSEC , PTAU &
         &, ZCGAZ ,ZPIZAZ,ZRAY1 ,ZRAY2 , ZREFZ, ZRJ  ,ZRK , ZRMUE &
         &, ZTAUAZ,ZTRA1 ,ZTRA2 ,ZTRCLD &
         &, Z_FRACT_DIFFUSE_CLOUDY &
         &, ZTAUAZ_N, ZPIZAZ_N, ZCGAZ_N)


    !     ------------------------------------------------------------------
    !*         3.2   SIX SPECTRAL INTERVALS
    !                ----------------------

    IIND(1)=1
    IIND(2)=2
    IIND(3)=1
    IIND(4)=2

    !*         3.1   DOWNWARD FLUXES
    !                ---------------
    JAJ = 2

    DO JL = KIDIA,KFDIA
      ZW(JL,1)=0._DP
      ZW(JL,2)=0._DP
      ZW(JL,3)=0._DP
      ZW(JL,4)=0._DP

      ZO(JL,1)=0._DP 
      ZO(JL,2)=0._DP 
      PFD(JL,KLEV+1)=((1._DP-PCLEAR(JL))*ZRJ(JL,JAJ,KLEV+1)&
           &+ PCLEAR(JL) *ZRJ0(JL,JAJ,KLEV+1)) * RSUN(KNU)
      PCD(JL,KLEV+1)= ZRJ0(JL,JAJ,KLEV+1) * RSUN(KNU)
    ENDDO
    DO JK = 1 , KLEV
      IKL = KLEV+1-JK
      DO JL = KIDIA,KFDIA
        ZW(JL,1)=ZW(JL,1)+PUD(JL,1,IKL)/ZRMUE(JL,IKL)
        ZW(JL,2)=ZW(JL,2)+PUD(JL,2,IKL)/ZRMUE(JL,IKL)
        ZW(JL,3)=ZW(JL,3)+PUD(JL,1,IKL)/ZRMU0(JL,IKL)
        ZW(JL,4)=ZW(JL,4)+PUD(JL,2,IKL)/ZRMU0(JL,IKL)

        ZO(JL,1)=ZO(JL,1)+POZ(JL,  IKL)/ZRMUE(JL,IKL)
        ZO(JL,2)=ZO(JL,2)+POZ(JL,  IKL)/ZRMU0(JL,IKL)
      ENDDO

      CALL SWTT1 ( KIDIA, KFDIA, KBDIM, KNU, 4, IIND, ZW, ZR )
      CALL SWUVO3( KIDIA, KFDIA, KBDIM, KNU, 2, ZO, ZT )

      DO JL = KIDIA,KFDIA
        ZDIFF(JL) = ZR(JL,1)*ZR(JL,2)*ZT(JL,1)*ZRJ(JL,JAJ,IKL) 
        ZDIRF(JL) = ZR(JL,3)*ZR(JL,4)*ZT(JL,2)*ZRJ0(JL,JAJ,IKL)
        PFD(JL,IKL) = ((1._DP-PCLEAR(JL)) * ZDIFF(JL)&
             &+PCLEAR(JL)  * ZDIRF(JL)) * RSUN(KNU)
        PCD(JL,IKL) = ZDIRF(JL) * RSUN(KNU)
      ENDDO
      IF (IKL == 1) THEN
        DO JL = KIDIA,KFDIA
          P_VIS_FRACT_DIFFUSE(JL) = ((1._DP-PCLEAR(JL)) * ZDIFF(JL)         &
               * Z_FRACT_DIFFUSE_CLOUDY(JL) + PCLEAR(JL) * ZDIRF(JL)        &
               * Z_FRACT_DIFFUSE_CLEAR(JL)) * RSUN(KNU)                     &
               / (PFD(JL,IKL) + EPSILON(1._DP))
        ENDDO
      END IF
    ENDDO

    DO JL=KIDIA,KFDIA
      ZDIFT(JL) = ZR(JL,1)*ZR(JL,2)*ZT(JL,1)*ZTRCLD(JL)
      ZDIRT(JL) = ZR(JL,3)*ZR(JL,4)*ZT(JL,2)*ZTRCLR(JL)
      PSUDU1(JL) = ((1._DP-PCLEAR(JL)) * ZDIFT(JL)            &
           &       + PCLEAR(JL) * ZDIRT(JL)) * RSUN(KNU)
    ENDDO


    !*         3.2   UPWARD FLUXES
    !                -------------

    DO JL = KIDIA,KFDIA
      PFU(JL,1) = ((1._DP-PCLEAR(JL))*ZDIFF(JL)*PALB(JL)&
           &+ PCLEAR(JL) * ZDIRF(JL) * PALB(JL))&
           &* RSUN(KNU)
      P_VIS(JL) = PFD(JL,1) - PFU(JL,1)
      PCU(JL,1) = ZDIRF(JL) * PALB(JL) * RSUN(KNU)
    ENDDO

    DO JK = 2 , KLEV+1
      IKM1=JK-1
      DO JL = KIDIA,KFDIA
        ZW(JL,1)=ZW(JL,1)+PUD(JL,1,IKM1)*1.66_DP
        ZW(JL,2)=ZW(JL,2)+PUD(JL,2,IKM1)*1.66_DP
        ZW(JL,3)=ZW(JL,3)+PUD(JL,1,IKM1)*1.66_DP
        ZW(JL,4)=ZW(JL,4)+PUD(JL,2,IKM1)*1.66_DP

        ZO(JL,1)=ZO(JL,1)+POZ(JL,  IKM1)*1.66_DP
        ZO(JL,2)=ZO(JL,2)+POZ(JL,  IKM1)*1.66_DP
      ENDDO

      CALL SWTT1 (KIDIA, KFDIA, KBDIM, KNU, 4, IIND, ZW, ZR )
      CALL SWUVO3(KIDIA, KFDIA, KBDIM, KNU, 2, ZO, ZT)

      DO JL = KIDIA,KFDIA
        ZDIFF(JL) = ZR(JL,1)*ZR(JL,2)*ZT(JL,1)*ZRK(JL,JAJ,JK)
        ZDIRF(JL) = ZR(JL,3)*ZR(JL,4)*ZT(JL,2)*ZRK0(JL,JAJ,JK)
        PFU(JL,JK) = ((1._DP-PCLEAR(JL)) * ZDIFF(JL)                &
             &       + PCLEAR(JL)  * ZDIRF(JL)) * RSUN(KNU)
        PCU(JL,JK) = ZDIRF(JL) * RSUN(KNU)
      ENDDO
    ENDDO

    RETURN
  END SUBROUTINE SW1S

  SUBROUTINE SWNI &
       &( KIDIA , KFDIA , KBDIM  , KLEV ,KNU &
       &, PSWEXT, PSWSSA, PSWASY &
       &, PAKI  , PALB, PCG &
       &, PCLD  , PCLEAR &
       &, PDSIG , POMEGA, POZ   , PRMU , PSEC, PTAU &
       &, PUD   , PWV   , PQS &
       &, PFDOWN, PFUP  , PCDOWN, PCUP , PSUDU2 &
       &, P_NIR, P_NIR_DIFFUSE &
       &)

    !**** *SWNI* - SHORTWAVE RADIATION, NEAR-INFRARED SPECTRAL INTERVALS

    !     PURPOSE.
    !     --------

    !          COMPUTES THE SHORTWAVE RADIATION FLUXES IN THE NEAR-INFRARED 
    !     SPECTRAL INTERVALS FOLLOWING FOUQUART AND BONNEL (1980).

    !**   INTERFACE.
    !     ----------

    !          *SWNI* IS CALLED FROM *SW*.


    !        IMPLICIT ARGUMENTS :
    !        --------------------

    !     ==== INPUTS ===
    !     ==== OUTPUTS ===

    !     METHOD.
    !     -------

    !          1. COMPUTES REFLECTIVITY/TRANSMISSIVITY CORRESPONDING TO
    !     CONTINUUM SCATTERING
    !          2. COMPUTES REFLECTIVITY/TRANSMISSIVITY CORRESPONDING FOR
    !     A GREY MOLECULAR ABSORPTION
    !          3. LAPLACE TRANSFORM ON THE PREVIOUS TO GET EFFECTIVE AMOUNTS
    !     OF ABSORBERS
    !          4. APPLY H2O AND U.M.G. TRANSMISSION FUNCTIONS
    !          5. MULTIPLY BY OZONE TRANSMISSION FUNCTION

    !     EXTERNALS.
    !     ----------

    !          *SWCLR*, *SWR*, *SWDE*, *SWTT*

    !     REFERENCE.
    !     ----------

    !        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
    !        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

    !     AUTHOR.
    !     -------
    !        JEAN-JACQUES MORCRETTE  *ECMWF*

    !     MODIFICATIONS.
    !     --------------
    !        ORIGINAL : 89-07-14
    !        94-11-15   J.-J. MORCRETTE    DIRECT/DIFFUSE ALBEDO
    !        95-12-07   J.-J. MORCRETTE    NEAR-INFRARED SW
    !        990128     JJMorcrette        Sunshine duration
    !        99-05-25   JJMorcrette        Revised aerosols

    !     ------------------------------------------------------------------


    !     DUMMY INTEGER SCALARS
    INTEGER :: KFDIA
    INTEGER :: KIDIA
    INTEGER :: KLEV
    INTEGER :: KBDIM
    INTEGER :: KNU

    !#include "yoeaer.h"
    !     ------------------------------------------------------------------

    !*       0.1   ARGUMENTS
    !              ---------

    REAL(DP):: PAKI(KBDIM,2,NSW)&
         &,  PALB(KBDIM) &
         &,  PCG(KBDIM,NSW,KLEV)   , PCLD(KBDIM,KLEV)&
         &,  PCLEAR(KBDIM)         , PDSIG(KBDIM,KLEV)&
         &,  POMEGA(KBDIM,NSW,KLEV), POZ(KBDIM,KLEV)&
         &,  PQS(KBDIM,KLEV)&
         &,  PRMU(KBDIM)           , PSEC(KBDIM)&
         &,  PTAU(KBDIM,NSW,KLEV)  , PUD(KBDIM,5,KLEV+1)&
         &,  PWV(KBDIM,KLEV)


    REAL(DP):: PFDOWN(KBDIM,KLEV+1)  , PFUP(KBDIM,KLEV+1)&
         &,  PCDOWN(KBDIM,KLEV+1)  , PCUP(KBDIM,KLEV+1)&
         &,  PSUDU2(KBDIM)

    REAL(DP) :: PSWEXT(KBDIM,KLEV,NSW), PSWSSA(KBDIM,KLEV,NSW) &
         &,        PSWASY(KBDIM,KLEV,NSW)

    REAL(DP) ::  P_NIR(KBDIM), P_NIR_DIFFUSE(KBDIM)

    !     ------------------------------------------------------------------

    !*       0.2   LOCAL ARRAYS
    !              ------------

    INTEGER :: IIND2(2), IIND3(3)
    REAL(DP):: ZCGAZ(KBDIM,KLEV)  , ZDIFF(KBDIM)         , ZDIRF(KBDIM)&
         &,  ZG(KBDIM)          , ZGG(KBDIM)&
         &,  ZPIZAZ(KBDIM,KLEV)&
         &,  ZRAYL(KBDIM)       , ZRAY1(KBDIM,KLEV+1)  , ZRAY2(KBDIM,KLEV+1)&
         &,  ZREF(KBDIM)        , ZREFZ(KBDIM,2,KLEV+1)&
         &,  ZRE1(KBDIM)        , ZRE2(KBDIM)&
         &,  ZRJ(KBDIM,6,KLEV+1), ZRJ0(KBDIM,6,KLEV+1)&
         &,  ZRK(KBDIM,6,KLEV+1), ZRK0(KBDIM,6,KLEV+1)&
         &,  ZRL(KBDIM,8)&
         &,  ZRMUE(KBDIM,KLEV+1), ZRMU0(KBDIM,KLEV+1)  , ZRMUZ(KBDIM)&
         &,  ZRNEB(KBDIM)       , ZRUEF(KBDIM,8)       , ZR1(KBDIM) &
         &,  ZR2(KBDIM,2)       , ZR3(KBDIM,3)         , ZR4(KBDIM)&
         &,  ZS(KBDIM)&
         &,  ZTAUAZ(KBDIM,KLEV) , ZTO1(KBDIM)          , ZTR(KBDIM,2,KLEV+1)&
         &,  ZTRA1(KBDIM,KLEV+1), ZTRA2(KBDIM,KLEV+1)&
         &,  ZTRCLD(KBDIM)      , ZTRCLR(KBDIM)&
         &,  ZTR1(KBDIM)        , ZTR2(KBDIM)&
         &,  ZW(KBDIM)          , ZW1(KBDIM)           , ZW2(KBDIM,2)&
         &,  ZW3(KBDIM,3)       , ZW4(KBDIM)           , ZW5(KBDIM) &
         &,  Z_FRACT_DIFFUSE_CLEAR(KBDIM), Z_FRACT_DIFFUSE_CLOUDY(KBDIM)

    REAL(DP) :: ZWC1(KBDIM), ZWC4(KBDIM), ZWC5(KBDIM),ZRC4(KBDIM), &
         &         ZRC1(KBDIM)
    REAL(DP) :: ZTAUAZ_N(KBDIM,KLEV), ZPIZAZ_N(KBDIM,KLEV), ZCGAZ_N(KBDIM,KLEV)
    REAL(DP) :: ZFDOWNC(KBDIM,KLEV+1)  , ZFUPC(KBDIM,KLEV+1)


    !     LOCAL INTEGER SCALARS
    INTEGER :: IABS, IKL, IKM1, JABS, JAJ, JAJP, JK, JKKI,&
         &JKKP4, JKL, JKLP1, JKM1, JL, JN, JN2J, JREF

    !     LOCAL REAL SCALARS
    REAL(DP):: ZAA, ZBB, ZCNEB, ZRE11, ZRKI, ZRMUM1, ZWH2O
    REAL(DP):: ZCHKS, ZCHKG, ZARGJ, ZARGK, ZARGJ0, ZARGK0

    REAL(DP) :: ZEPS
    !*         0.        Initializations
    !                       ---------------
    ZEPS=EPSILON(1.0_DP)

    !     ------------------------------------------------------------------

    !*         1.     NEAR-INFRARED SPECTRAL INTERVAL (0.68-4.00 MICRON)
    !                 --------------------------------------------------



    !*         1.1    OPTICAL THICKNESS FOR RAYLEIGH SCATTERING
    !                 -----------------------------------------


    DO JL = KIDIA,KFDIA
      ZRMUM1 = 1._DP - PRMU(JL)
      ZRAYL(JL) =  RRAY(KNU,1) + ZRMUM1   * (RRAY(KNU,2) + ZRMUM1 &
           &* (RRAY(KNU,3) + ZRMUM1   * (RRAY(KNU,4) + ZRMUM1 &
           &* (RRAY(KNU,5) + ZRMUM1   *  RRAY(KNU,6)     ))))
      ZRAYL(JL) =  MAX(ZRAYL(JL),0._DP)
    ENDDO


    !     ------------------------------------------------------------------

    !*         2.    CONTINUUM SCATTERING CALCULATIONS
    !                ---------------------------------


    !*         2.1   CLEAR-SKY FRACTION OF THE COLUMN
    !                --------------------------------


    CALL SWCLR &
         &( KIDIA , KFDIA , KBDIM , KLEV , KNU &
         &, PSWEXT, PSWSSA, PSWASY &
         &, PALB , PDSIG , ZRAYL,  PSEC &
         &, ZTAUAZ, ZCGAZ , ZPIZAZ, ZRAY1 , ZRAY2, ZREFZ, ZRJ0 &
         &, ZRK0  , ZRMU0 , ZTRA1, ZTRA2, ZTRCLR &
         &, Z_FRACT_DIFFUSE_CLEAR &
         &, ZTAUAZ_N, ZPIZAZ_N, ZCGAZ_N)

    !*         2.1.A   SCATTERING CALCULATIONS WITH GREY MOLECULAR ABSORPTION
    !                ------------------------------------------------------


    JN = 2

    DO JABS=1,2


      !*         A.1  SURFACE CONDITIONS
      !               ------------------


      DO JL = KIDIA,KFDIA
        ZREFZ(JL,2,1) = PALB(JL)
        ZREFZ(JL,1,1) = PALB(JL)
      ENDDO


      !*         A.2  INTRODUCING CLOUD EFFECTS
      !               -------------------------


      DO JK = 2 , KLEV+1
        JKM1 = JK - 1
        IKL=KLEV+1-JKM1
        DO JL = KIDIA,KFDIA
          ZAA=PUD(JL,JABS,JKM1)
          ZBB=ZAA
          ZRKI = PAKI(JL,JABS,KNU)
          ZCHKS= MIN(200._DP,ZRKI * ZAA * 1.66_DP)
          ZS(JL) = EXP(-ZCHKS)
          ZCHKG= MIN(200._DP,ZRKI * ZAA / ZRMU0(JL,JK))
          ZG(JL) = EXP(-ZCHKG)
          ZTR1(JL) = 0._DP
          ZRE1(JL) = 0._DP
          ZTR2(JL) = 0._DP
          ZRE2(JL) = 0._DP

          !-combining aerosol+gaseous absorption

          ZTO1(JL) =  ZTAUAZ_N(JL,JKM1) &
               &      +ZBB*ZRKI
          ZW(JL)   =  ZTAUAZ_N(JL,JKM1)*ZPIZAZ_N(JL,JKM1)
          ZGG(JL)  =  ZTAUAZ_N(JL,JKM1)*ZPIZAZ_N(JL,JKM1)*ZCGAZ_N(JL,JKM1)

          !-asymmetry parameter

          IF(ZW(JL)>ZEPS) THEN
            ZGG(JL)  = ZGG(JL)/ZW(JL)
          ELSE
            ZGG(JL)  = 0.0_DP
          ENDIF

          !-single scattering albedo

          IF(ZTO1(JL)>ZEPS) THEN
            ZW(JL)   = ZW(JL)/ZTO1(JL)
          ELSE
            ZW(JL)   = 1.0_DP
          ENDIF

          ZW(JL)=MIN(ZW(JL),1.0_DP-ZEPS)
          !---

          ZREF(JL) = ZREFZ(JL,1,JKM1)
          ZRMUZ(JL) = ZRMU0(JL,JK)

        ENDDO


        CALL SWDE ( KIDIA, KFDIA, KBDIM &
             &, ZGG  , ZREF , ZRMUZ, ZTO1, ZW &
             &, ZRE1 , ZRE2 , ZTR1 , ZTR2     )

        DO JL = KIDIA,KFDIA

          ZREFZ(JL,2,JK) = (ZRAY1(JL,JKM1)&
               &+ ZREFZ(JL,2,JKM1) * ZTRA1(JL,JKM1)&
               &* ZTRA2(JL,JKM1)) * ZG(JL) * ZS(JL)

          ZTR(JL,2,JKM1)=ZTRA1(JL,JKM1) * ZG(JL)

          ZREFZ(JL,1,JK)=(ZRAY1(JL,JKM1)&
               &+ZREFZ(JL,1,JKM1)*ZTRA1(JL,JKM1)*ZTRA2(JL,JKM1)&
               &/(1._DP-ZRAY2(JL,JKM1)*ZREFZ(JL,1,JKM1)))*ZG(JL)*ZS(JL)

          ZTR(JL,1,JKM1)=(ZTRA1(JL,JKM1)/(1._DP-ZRAY2(JL,JKM1)&
               &* ZREFZ(JL,1,JKM1)))&
               &* ZG(JL)
        ENDDO
      ENDDO


      !*         3.3  REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL
      !               -------------------------------------------------


      DO JREF=1,2

        JN = JN + 1


        DO JL = KIDIA,KFDIA
          ZRJ0(JL,JN,KLEV+1) = 1._DP
          ZRK0(JL,JN,KLEV+1) = ZREFZ(JL,JREF,KLEV+1)
        ENDDO

        DO JK = 1 , KLEV
          JKL = KLEV+1 - JK
          JKLP1 = JKL + 1
          DO JL = KIDIA,KFDIA
            ZRE11 = ZRJ0(JL,JN,JKLP1) * ZTR(JL,JREF,JKL)
            ZRJ0(JL,JN,JKL) = ZRE11
            ZRK0(JL,JN,JKL) = ZRE11 * ZREFZ(JL,JREF,JKL)
          ENDDO
        ENDDO
      ENDDO
    ENDDO


    !     ------------------------------------------------------------------

    !*         4.    INVERT GREY AND CONTINUUM FLUXES
    !                --------------------------------



    !*         4.1   UPWARD (ZRK) AND DOWNWARD (ZRJ) PSEUDO-FLUXES
    !                ---------------------------------------------



    DO JK = 1 , KLEV+1
      DO JAJ = 1 , 5 , 2
        JAJP = JAJ + 1
        DO JL = KIDIA,KFDIA
          ZRJ0(JL,JAJ,JK)=        ZRJ0(JL,JAJ,JK) - ZRJ0(JL,JAJP,JK)
          ZRK0(JL,JAJ,JK)=        ZRK0(JL,JAJ,JK) - ZRK0(JL,JAJP,JK)
          ZRJ0(JL,JAJ,JK)= MAX( ZRJ0(JL,JAJ,JK) , EPSILON(1._DP) )
          ZRK0(JL,JAJ,JK)= MAX( ZRK0(JL,JAJ,JK) , EPSILON(1._DP) )
        ENDDO
      ENDDO
    ENDDO

    DO JK = 1 , KLEV+1
      DO JAJ = 2 , 6 , 2
        DO JL = KIDIA,KFDIA
          ZRJ0(JL,JAJ,JK)= MAX( ZRJ0(JL,JAJ,JK) , EPSILON(1._DP) )
          ZRK0(JL,JAJ,JK)= MAX( ZRK0(JL,JAJ,JK) , EPSILON(1._DP) )
        ENDDO
      ENDDO
    ENDDO
    !
    !*         4.2    EFFECTIVE ABSORBER AMOUNTS BY INVERSE LAPLACE
    !                 ---------------------------------------------


    DO JK = 1 , KLEV+1
      JKKI = 1
      DO JAJ = 1 , 2
        IIND2(1)=JAJ
        IIND2(2)=JAJ
        DO JN = 1 , 2
          JN2J = JN + 2 * JAJ
          JKKP4 = JKKI + 4

          !*         4.2.1  EFFECTIVE ABSORBER AMOUNTS
          !                 --------------------------


          DO JL = KIDIA,KFDIA
            ZARGJ0     = MAX( ZRJ0(JL,JN,JK) / ZRJ0(JL,JN2J,JK),1._DP)
            ZARGK0     = MAX( ZRK0(JL,JN,JK) / ZRK0(JL,JN2J,JK),1._DP)
            ZW2(JL,1) = LOG( ZARGJ0 )/ PAKI(JL,JAJ,KNU)
            ZW2(JL,2) = LOG( ZARGK0 )/ PAKI(JL,JAJ,KNU)
          ENDDO
          !
          !*         4.2.2  TRANSMISSION FUNCTION
          !                 ---------------------


          CALL SWTT1 ( KIDIA,KFDIA,KBDIM, KNU, 2, IIND2 &
               &, ZW2 &
               &, ZR2                              )

          DO JL = KIDIA,KFDIA
            ZRL(JL,JKKI) = ZR2(JL,1)
            ZRUEF(JL,JKKI) = ZW2(JL,1)
            ZRL(JL,JKKP4) = ZR2(JL,2)
            ZRUEF(JL,JKKP4) = ZW2(JL,2)
          ENDDO

          JKKI=JKKI+1
        ENDDO
      ENDDO

      !*         4.3    UPWARD AND DOWNWARD FLUXES WITH H2O AND UMG ABSORPTION
      !                 ------------------------------------------------------


      DO JL = KIDIA,KFDIA
        ZFDOWNC(JL,JK) = ZRJ0(JL,1,JK) * ZRL(JL,1) * ZRL(JL,3)&
             &+ ZRJ0(JL,2,JK) * ZRL(JL,2) * ZRL(JL,4)
        ZFUPC(JL,JK)   = ZRK0(JL,1,JK) * ZRL(JL,5) * ZRL(JL,7)&
             &+ ZRK0(JL,2,JK) * ZRL(JL,6) * ZRL(JL,8)
      ENDDO
    ENDDO




    !*         2.2   CLOUDY FRACTION OF THE COLUMN
    !                -----------------------------


    CALL SWR &
         &( KIDIA , KFDIA , KBDIM , KLEV  , KNU &
         &, PALB , PCG   , PCLD , POMEGA, PSEC , PTAU &
         &, ZCGAZ , ZPIZAZ, ZRAY1, ZRAY2 , ZREFZ, ZRJ  , ZRK, ZRMUE &
         &, ZTAUAZ, ZTRA1 , ZTRA2, ZTRCLD &
         &, Z_FRACT_DIFFUSE_CLOUDY &
         &, ZTAUAZ_N, ZPIZAZ_N, ZCGAZ_N)
    !     ------------------------------------------------------------------

    !*         3.    SCATTERING CALCULATIONS WITH GREY MOLECULAR ABSORPTION
    !                ------------------------------------------------------


    JN = 2

    DO JABS=1,2


      !*         3.1  SURFACE CONDITIONS
      !               ------------------


      DO JL = KIDIA,KFDIA
        ZREFZ(JL,2,1) = PALB(JL)
        ZREFZ(JL,1,1) = PALB(JL)
      ENDDO


      !*         3.2  INTRODUCING CLOUD EFFECTS
      !               -------------------------

      DO JK = 2 , KLEV+1
        JKM1 = JK - 1
        IKL=KLEV+1-JKM1
        DO JL = KIDIA,KFDIA

          ZRNEB(JL) = PCLD(JL,JKM1)

          IF (JABS == 1.AND. ZRNEB(JL) > 2._DP*EPSILON(1._DP)) THEN
            ZWH2O=MAX(PWV(JL,IKL),EPSILON(1._DP))
            ZCNEB=MAX(EPSILON(1._DP),MIN(ZRNEB(JL),1._DP-EPSILON(1._DP)))
            ZBB=PUD(JL,JABS,JKM1)*PQS(JL,IKL)/ZWH2O
            ZAA=MAX((PUD(JL,JABS,JKM1)-ZCNEB*ZBB)/(1._DP-ZCNEB),EPSILON(1._DP))
          ELSE
            ZAA=PUD(JL,JABS,JKM1)
            ZBB=ZAA
          ENDIF
          ZRKI = PAKI(JL,JABS,KNU)
          ZCHKS= MIN(200._DP,ZRKI * ZAA * 1.66_DP)
          ZS(JL) = EXP(-ZCHKS)
          ZCHKG= MIN(200._DP,ZRKI * ZAA / ZRMUE(JL,JK))
          ZG(JL) = EXP(-ZCHKG)
          ZTR1(JL) = 0._DP
          ZRE1(JL) = 0._DP
          ZTR2(JL) = 0._DP
          ZRE2(JL) = 0._DP

!combining aerosols+clouds+gaseous absorption

          ZTO1(JL) =  ZTAUAZ_N(JL,JKM1) &
               &      +PTAU(JL,KNU,JKM1) &
               &      +ZBB*ZRKI
          ZW(JL)   =  ZTAUAZ_N(JL,JKM1)*ZPIZAZ_N(JL,JKM1) &
               &      +PTAU(JL,KNU,JKM1)*POMEGA(JL,KNU,JKM1)
          ZGG(JL)  =  ZTAUAZ_N(JL,JKM1)*ZPIZAZ_N(JL,JKM1)*ZCGAZ_N(JL,JKM1)&
               &      +PTAU(JL,KNU,JKM1)*POMEGA(JL,KNU,JKM1)*PCG(JL,KNU,JKM1)

          IF(ZW(JL)>ZEPS) THEN
            ZGG(JL)  = ZGG(JL)/ZW(JL)
          ELSE
            ZGG(JL)  = 0.0_DP
          ENDIF

          IF(ZTO1(JL)>ZEPS) THEN
            ZW(JL)   = ZW(JL)/ZTO1(JL)
          ELSE
            ZW(JL)=1.0_DP
          ENDIF

          ZW(JL)=MIN(ZW(JL),1.0_DP-ZEPS)

          !----
          ZREF(JL) = ZREFZ(JL,1,JKM1)
          ZRMUZ(JL) = ZRMUE(JL,JK)
        ENDDO

        CALL SWDE ( KIDIA, KFDIA, KBDIM &
             &, ZGG  , ZREF , ZRMUZ, ZTO1, ZW &
             &, ZRE1 , ZRE2 , ZTR1 , ZTR2     )

        DO JL = KIDIA,KFDIA


          ZREFZ(JL,2,JK) = (1._DP-ZRNEB(JL)) * (ZRAY1(JL,JKM1)&
               &+ ZREFZ(JL,2,JKM1) * ZTRA1(JL,JKM1)&
               &* ZTRA2(JL,JKM1) ) * ZG(JL) * ZS(JL)&
               &+ ZRNEB(JL) * ZRE1(JL)

          ZTR(JL,2,JKM1)=ZRNEB(JL)*ZTR1(JL)&
               &+ (ZTRA1(JL,JKM1)) * ZG(JL) * (1._DP-ZRNEB(JL))


          ZREFZ(JL,1,JK)=(1._DP-ZRNEB(JL))*(ZRAY1(JL,JKM1)&
               &+ZREFZ(JL,1,JKM1)*ZTRA1(JL,JKM1)*ZTRA2(JL,JKM1)&
               &/(1._DP-ZRAY2(JL,JKM1)*ZREFZ(JL,1,JKM1)))*ZG(JL)*ZS(JL)&
               &+ ZRNEB(JL) * ZRE2(JL)

          ZTR(JL,1,JKM1)= ZRNEB(JL)*ZTR2(JL)&
               &+ (ZTRA1(JL,JKM1)/(1._DP-ZRAY2(JL,JKM1)&
               &* ZREFZ(JL,1,JKM1)))&
               &* ZG(JL) * (1._DP -ZRNEB(JL))

        ENDDO
      ENDDO


      !*         3.3  REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL
      !               -------------------------------------------------


      DO JREF=1,2

        JN = JN + 1


        DO JL = KIDIA,KFDIA
          ZRJ(JL,JN,KLEV+1) = 1._DP
          ZRK(JL,JN,KLEV+1) = ZREFZ(JL,JREF,KLEV+1)
        ENDDO

        DO JK = 1 , KLEV
          JKL = KLEV+1 - JK
          JKLP1 = JKL + 1
          DO JL = KIDIA,KFDIA
            ZRE11 = ZRJ(JL,JN,JKLP1) * ZTR(JL,JREF,JKL)
            ZRJ(JL,JN,JKL) = ZRE11
            ZRK(JL,JN,JKL) = ZRE11 * ZREFZ(JL,JREF,JKL)
          ENDDO
        ENDDO
      ENDDO
    ENDDO


    !     ------------------------------------------------------------------

    !*         4.    INVERT GREY AND CONTINUUM FLUXES
    !                --------------------------------



    !*         4.1   UPWARD (ZRK) AND DOWNWARD (ZRJ) PSEUDO-FLUXES
    !                ---------------------------------------------



    DO JK = 1 , KLEV+1
      DO JAJ = 1 , 5 , 2
        JAJP = JAJ + 1
        DO JL = KIDIA,KFDIA
          ZRJ(JL,JAJ,JK)=        ZRJ(JL,JAJ,JK) - ZRJ(JL,JAJP,JK)
          ZRK(JL,JAJ,JK)=        ZRK(JL,JAJ,JK) - ZRK(JL,JAJP,JK)
          ZRJ(JL,JAJ,JK)= MAX( ZRJ(JL,JAJ,JK) , EPSILON(1._DP) )
          ZRK(JL,JAJ,JK)= MAX( ZRK(JL,JAJ,JK) , EPSILON(1._DP) )
        ENDDO
      ENDDO
    ENDDO

    DO JK = 1 , KLEV+1
      DO JAJ = 2 , 6 , 2
        DO JL = KIDIA,KFDIA
          ZRJ(JL,JAJ,JK)= MAX( ZRJ(JL,JAJ,JK) , EPSILON(1._DP) )
          ZRK(JL,JAJ,JK)= MAX( ZRK(JL,JAJ,JK) , EPSILON(1._DP) )
        ENDDO
      ENDDO
    ENDDO
    !
    !*         4.2    EFFECTIVE ABSORBER AMOUNTS BY INVERSE LAPLACE
    !                 ---------------------------------------------


    DO JK = 1 , KLEV+1
      JKKI = 1
      DO JAJ = 1 , 2
        IIND2(1)=JAJ
        IIND2(2)=JAJ
        DO JN = 1 , 2
          JN2J = JN + 2 * JAJ
          JKKP4 = JKKI + 4

          !*         4.2.1  EFFECTIVE ABSORBER AMOUNTS
          !                 --------------------------


          DO JL = KIDIA,KFDIA
            ZARGJ     = MAX( ZRJ(JL,JN,JK) / ZRJ(JL,JN2J,JK),1._DP)
            ZARGK     = MAX( ZRK(JL,JN,JK) / ZRK(JL,JN2J,JK),1._DP)
            ZW2(JL,1) = LOG( ZARGJ )/ PAKI(JL,JAJ,KNU)
            ZW2(JL,2) = LOG( ZARGK )/ PAKI(JL,JAJ,KNU)
          ENDDO
          !
          !*         4.2.2  TRANSMISSION FUNCTION
          !                 ---------------------


          CALL SWTT1 ( KIDIA,KFDIA,KBDIM, KNU, 2, IIND2 &
               &, ZW2 &
               &, ZR2                              )

          DO JL = KIDIA,KFDIA
            ZRL(JL,JKKI) = ZR2(JL,1)
            ZRUEF(JL,JKKI) = ZW2(JL,1)
            ZRL(JL,JKKP4) = ZR2(JL,2)
            ZRUEF(JL,JKKP4) = ZW2(JL,2)
          ENDDO

          JKKI=JKKI+1
        ENDDO
      ENDDO

      !*         4.3    UPWARD AND DOWNWARD FLUXES WITH H2O AND UMG ABSORPTION
      !                 ------------------------------------------------------


      DO JL = KIDIA,KFDIA
        PFDOWN(JL,JK) = ZRJ(JL,1,JK) * ZRL(JL,1) * ZRL(JL,3)&
             &+ ZRJ(JL,2,JK) * ZRL(JL,2) * ZRL(JL,4)
        PFUP(JL,JK)   = ZRK(JL,1,JK) * ZRL(JL,5) * ZRL(JL,7)&
             &+ ZRK(JL,2,JK) * ZRL(JL,6) * ZRL(JL,8)
      ENDDO
    ENDDO



    !     ------------------------------------------------------------------

    !*         5.    MOLECULAR ABSORPTION ON CLEAR-SKY FLUXES
    !                ----------------------------------------



    !*         5.1   DOWNWARD FLUXES
    !                ---------------


    !JAJ = 2
    IIND3(1)=1
    IIND3(2)=2
    IIND3(3)=3

    DO JL = KIDIA,KFDIA
      ZW3(JL,1)=0._DP
      ZW3(JL,2)=0._DP
      ZW3(JL,3)=0._DP
      ZW4(JL)  =0._DP
      ZW5(JL)  =0._DP
      ZR4(JL)  =1._DP
    ENDDO
    DO JK = 1 , KLEV
      IKL = KLEV+1-JK
      DO JL = KIDIA,KFDIA
        ZW3(JL,1)=ZW3(JL,1)+PUD(JL,1,IKL)/ZRMU0(JL,IKL)
        ZW3(JL,2)=ZW3(JL,2)+PUD(JL,2,IKL)/ZRMU0(JL,IKL)
        ZW3(JL,3)=ZW3(JL,3)+POZ(JL,  IKL)/ZRMU0(JL,IKL)
        ZW4(JL)  =ZW4(JL)  +PUD(JL,4,IKL)/ZRMU0(JL,IKL)
        ZW5(JL)  =ZW5(JL)  +PUD(JL,5,IKL)/ZRMU0(JL,IKL)
      ENDDO

      CALL SWTT1 ( KIDIA,KFDIA,KBDIM, KNU, 3, IIND3 &
           &, ZW3 &
           &, ZR3                              )

      DO JL = KIDIA,KFDIA
        ZR4(JL) = EXP(-RSWCE(KNU)*ZW4(JL)-RSWCP(KNU)*ZW5(JL))
      ENDDO
    ENDDO

    DO JL=KIDIA,KFDIA
      ZDIFF(JL) = ZR3(JL,1)*ZR3(JL,2)*ZR3(JL,3)*ZR4(JL)*ZTRCLD(JL)
      ZDIRF(JL) = ZR3(JL,1)*ZR3(JL,2)*ZR3(JL,3)*ZR4(JL)*ZTRCLR(JL)
      PSUDU2(JL) = ((1._DP-PCLEAR(JL)) * ZDIFF(JL)&
           &+PCLEAR(JL) * ZDIRF(JL)) * RSUN(KNU)
    ENDDO



    !*         6.     INTRODUCTION OF OZONE AND H2O CONTINUUM ABSORPTION
    !                 --------------------------------------------------

    IABS=3

    !*         6.1    DOWNWARD FLUXES
    !                 ---------------

    DO JL = KIDIA,KFDIA
      ZW1(JL)=0._DP
      ZW4(JL)=0._DP
      ZW5(JL)=0._DP
      ZR1(JL)=0._DP
      ZWC1(JL)=0._DP
      ZWC4(JL)=0._DP
      ZWC5(JL)=0._DP
      ZRC1(JL)=0._DP
      PFDOWN(JL,KLEV+1) = ((1._DP-PCLEAR(JL))*PFDOWN(JL,KLEV+1)&
           &+ PCLEAR(JL) * ZFDOWNC(JL,KLEV+1)) * RSUN(KNU)
      PCDOWN(JL,KLEV+1) = ZFDOWNC(JL,KLEV+1) * RSUN(KNU)
    ENDDO

    DO JK = 1 , KLEV
      IKL=KLEV+1-JK
      DO JL = KIDIA,KFDIA
        ZW1(JL) = ZW1(JL)+POZ(JL,  IKL)/ZRMUE(JL,IKL)
        ZW4(JL) = ZW4(JL)+PUD(JL,4,IKL)/ZRMUE(JL,IKL)
        ZW5(JL) = ZW5(JL)+PUD(JL,5,IKL)/ZRMUE(JL,IKL)
        ZR4(JL) = EXP(-RSWCE(KNU)*ZW4(JL)-RSWCP(KNU)*ZW5(JL))

        !gases absorption for clear sky fluxes

        ZWC1(JL) = ZWC1(JL)+POZ(JL,  IKL)/ZRMU0(JL,IKL)
        ZWC4(JL) = ZWC4(JL)+PUD(JL,4,IKL)/ZRMU0(JL,IKL)
        ZWC5(JL) = ZWC5(JL)+PUD(JL,5,IKL)/ZRMU0(JL,IKL)
        ZRC4(JL) = EXP(-RSWCE(KNU)*ZWC4(JL)-RSWCP(KNU)*ZWC5(JL))
      ENDDO

      CALL SWTT ( KIDIA,KFDIA,KBDIM, KNU, IABS, ZW1, ZR1 )

      !gases absorption for clear sky fluxes

      CALL SWTT ( KIDIA,KFDIA,KBDIM, KNU, IABS, ZWC1, ZRC1 )

      IF (IKL == 1) THEN
        DO JL = KIDIA,KFDIA
          P_NIR_DIFFUSE(JL) = ((1._DP-PCLEAR(JL))*ZR1(JL)*ZR4(JL)                       &
               & * PFDOWN(JL,IKL)*Z_FRACT_DIFFUSE_CLOUDY(JL)                                &
               & + PCLEAR(JL)*ZRC1(JL)*ZRC4(JL)*ZFDOWNC(JL,IKL)*Z_FRACT_DIFFUSE_CLEAR(JL))  &
               & * RSUN(KNU)
        ENDDO
      END IF
      DO JL = KIDIA,KFDIA
        PFDOWN(JL,IKL) = ((1._DP-PCLEAR(JL))*ZR1(JL)*ZR4(JL)*PFDOWN(JL,IKL)&
             &+PCLEAR(JL)*ZRC1(JL)*ZRC4(JL)*ZFDOWNC(JL,IKL)) * RSUN(KNU)
        PCDOWN(JL,IKL) = ZRC1(JL)*ZRC4(JL)*ZFDOWNC(JL,IKL) * RSUN(KNU)
      ENDDO
      IF (IKL == 1) THEN
        DO JL = KIDIA,KFDIA
          P_NIR_DIFFUSE(JL) = P_NIR_DIFFUSE(JL)/(PFDOWN(JL,IKL)+EPSILON(1._DP))
        ENDDO
      END IF
    ENDDO


    !*         6.2    UPWARD FLUXES
    !                 -------------

    DO JL = KIDIA,KFDIA
      PFUP(JL,1) = ((1._DP-PCLEAR(JL))*ZR1(JL)*ZR4(JL) * PFUP(JL,1)&
           &+PCLEAR(JL)*ZRC1(JL)*ZRC4(JL)*ZFDOWNC(JL,1)*PALB(JL)) * RSUN(KNU)
      P_NIR(JL) = PFDOWN(JL,1) - PFUP(JL,1)
      PCUP(JL,1) = ZRC1(JL)*ZRC4(JL)*ZFDOWNC(JL,1)*PALB(JL)* RSUN(KNU)
    ENDDO

    DO JK = 2 , KLEV+1
      IKM1=JK-1
      DO JL = KIDIA,KFDIA
        ZW1(JL) = ZW1(JL)+POZ(JL  ,IKM1)*1.66_DP
        ZW4(JL) = ZW4(JL)+PUD(JL,4,IKM1)*1.66_DP
        ZW5(JL) = ZW5(JL)+PUD(JL,5,IKM1)*1.66_DP
        ZR4(JL) = EXP(-RSWCE(KNU)*ZW4(JL)-RSWCP(KNU)*ZW5(JL))
      ENDDO

      CALL SWTT ( KIDIA,KFDIA,KBDIM, KNU, IABS, ZW1, ZR1 )

      DO JL = KIDIA,KFDIA
        PFUP(JL,JK) = ((1._DP-PCLEAR(JL))*ZR1(JL)*ZR4(JL) * PFUP(JL,JK)&
             &+PCLEAR(JL)*ZR1(JL)*ZR4(JL) * ZFUPC(JL,JK)) * RSUN(KNU)
        PCUP(JL,JK) = ZR1(JL)*ZR4(JL) * ZFUPC(JL,JK) * RSUN(KNU)
      ENDDO
    ENDDO

    !     ------------------------------------------------------------------

    RETURN
  END SUBROUTINE SWNI

  SUBROUTINE SWU &
       &( KIDIA, KFDIA , KBDIM  , KLEV &
       &, PSCT , PCARDI, PPMB , PPSOL, PRMU0, PTAVE, PWV &
       &, PAKI , PDSIG , PFACT, PRMU , PSEC , PUD &
       &)

    !**** *SWU* - SHORTWAVE RADIATION, ABSORBER AMOUNTS

    !     PURPOSE.
    !     --------
    !           COMPUTES THE ABSORBER AMOUNTS USED IN SHORTWAVE RADIATION
    !     CALCULATIONS

    !**   INTERFACE.
    !     ----------
    !          *SWU* IS CALLED BY *SW*


    !        IMPLICIT ARGUMENTS :
    !        --------------------

    !     ==== INPUTS ===
    !     ==== OUTPUTS ===

    !     METHOD.
    !     -------

    !          1. COMPUTES ABSORBER AMOUNTS WITH TEMPERATURE AND PRESSURE
    !     SCALING.

    !     EXTERNALS.
    !     ----------

    !          *SWTT*

    !     REFERENCE.
    !     ----------

    !        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
    !        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

    !     AUTHOR.
    !     -------
    !        JEAN-JACQUES MORCRETTE  *ECMWF*

    !     MODIFICATIONS.
    !     --------------
    !        ORIGINAL : 89-07-14
    !        June 2005: CCagnazzo from 4 to 6 bands, following 26Rcycle Morcrette Code

    !     ------------------------------------------------------------------


    !     DUMMY INTEGER SCALARS
    INTEGER :: KFDIA
    INTEGER :: KIDIA
    INTEGER :: KLEV
    INTEGER :: KBDIM

    !     DUMMY REAL SCALARS
    REAL(DP):: PSCT



    !     ------------------------------------------------------------------

    !*       0.1   ARGUMENTS
    !              ---------

    REAL(DP):: PPMB(KBDIM,KLEV+1), PPSOL(KBDIM)&
         &,  PRMU0(KBDIM)      , PTAVE(KBDIM,KLEV) , PWV(KBDIM,KLEV) &
         &,  PCARDI(KBDIM,KLEV)

    REAL(DP):: PAKI(KBDIM,2,NSW)&
         &,  PDSIG(KBDIM,KLEV) , PFACT(KBDIM)      , PRMU(KBDIM)&
         &,  PSEC(KBDIM)       , PUD(KBDIM,5,KLEV+1)

    !     ------------------------------------------------------------------

    !*       0.2   LOCAL ARRAYS
    !              ------------

    INTEGER :: IIND(2)
    REAL(DP):: ZN175(KBDIM), ZN190(KBDIM), ZO175(KBDIM)&
         &,  ZO190(KBDIM), ZSIGN(KBDIM)&
         &,  ZR(KBDIM,2) , ZSIGO(KBDIM), ZUD(KBDIM,2)

    !     LOCAL INTEGER SCALARS
    INTEGER :: JA, JK, JKL, JKLP1, JKP1, JL, JNU

    !     LOCAL REAL SCALARS
    REAL(DP):: ZDSCO2, ZDSH2O, ZFPPW, ZRTH, ZRTU, ZWH2O


    !     ------------------------------------------------------------------

    !*         1.     COMPUTES AMOUNTS OF ABSORBERS
    !                 -----------------------------


    IIND(1)=1
    IIND(2)=2


    !*         1.1    INITIALIZES QUANTITIES
    !                 ----------------------


    DO JL = KIDIA,KFDIA
      PUD(JL,1,KLEV+1)=0._DP
      PUD(JL,2,KLEV+1)=0._DP
      PUD(JL,3,KLEV+1)=0._DP
      PUD(JL,4,KLEV+1)=0._DP
      PUD(JL,5,KLEV+1)=0._DP
      PFACT(JL)= PRMU0(JL) * PSCT
      PRMU(JL)=SQRT(1224._DP* PRMU0(JL) * PRMU0(JL) + 1._DP) / 35._DP
      PSEC(JL)=1._DP/PRMU(JL)
    ENDDO

    !*          1.3    AMOUNTS OF ABSORBERS
    !                  --------------------


    DO JL= KIDIA,KFDIA
      ZUD(JL,1) = 0._DP
      ZUD(JL,2) = 0._DP
      ZO175(JL) = PPSOL(JL)** RPDU1
      ZO190(JL) = PPSOL(JL)** RPDH1
      ZSIGO(JL) = PPSOL(JL)
    ENDDO

    DO JK = 1 , KLEV
      JKP1 = JK + 1
      JKL = KLEV+1 - JK
      JKLP1 = JKL+1
      DO JL = KIDIA,KFDIA
        ZRTH=(RTH2O/PTAVE(JL,JK))**RTDH2O
        ZRTU=(RTUMG/PTAVE(JL,JK))**RTDUMG
        ZWH2O = MAX (PWV(JL,JKL) , EPSILON(1._DP) )
        ZSIGN(JL) = 100._DP * PPMB(JL,JKP1)
        PDSIG(JL,JK) = (ZSIGO(JL) - ZSIGN(JL))/PPSOL(JL)
        ZN175(JL) = ZSIGN(JL) ** RPDU1
        ZN190(JL) = ZSIGN(JL) ** RPDH1
        ZDSCO2 = ZO175(JL) - ZN175(JL)
        ZDSH2O = ZO190(JL) - ZN190(JL)
        PUD(JL,1,JK) = RPNH * ZDSH2O * ZWH2O  * ZRTH
        PUD(JL,2,JK) = RPNU * ZDSCO2 * PCARDI(JL,JKL) * ZRTU
        ZFPPW=1.6078_DP*ZWH2O/(1._DP+0.608_DP*ZWH2O)
        PUD(JL,4,JK)=PUD(JL,1,JK)*ZFPPW
        PUD(JL,5,JK)=PUD(JL,1,JK)*(1._DP-ZFPPW)
        ZUD(JL,1) = ZUD(JL,1) + PUD(JL,1,JK)
        ZUD(JL,2) = ZUD(JL,2) + PUD(JL,2,JK)
        ZSIGO(JL) = ZSIGN(JL)
        ZO175(JL) = ZN175(JL)
        ZO190(JL) = ZN190(JL)
      ENDDO
    ENDDO


    !*         1.4    COMPUTES CLEAR-SKY GREY ABSORPTION COEFFICIENTS
    !                 -----------------------------------------------


    DO JA = 1,2
      DO JL = KIDIA,KFDIA
        ZUD(JL,JA) = ZUD(JL,JA) * PSEC(JL)
      ENDDO
    ENDDO

    !! DO JNU= 2,NSW
    DO JNU=4,NSW   !! new

      CALL SWTT1 ( KIDIA,KFDIA,KBDIM, JNU, 2, IIND &
           &, ZUD &
           &, ZR                            )

      DO JA = 1,2
        DO JL = KIDIA,KFDIA
          PAKI(JL,JA,JNU) = -LOG( ZR(JL,JA) ) / ZUD(JL,JA)
        ENDDO
      ENDDO
    ENDDO


    !     ------------------------------------------------------------------

    RETURN
  END SUBROUTINE SWU

  SUBROUTINE SWCLR &
       &( KIDIA , KFDIA , KBDIM  , KLEV , KNU &
       &, PSWEXT, PSWSSA, PSWASY &
       &, PALB , PDSIG , PRAYL , PSEC &
       &, PTAUAZ, PCGAZ , PPIZAZ, PRAY1 , PRAY2 , PREFZ , PRJ  &
       &, PRK   , PRMU0 , PTRA1 , PTRA2 , PTRCLR &
       &, PDIFFUSE, PTAUAZ_N, PPIZAZ_N, PCGAZ_N)

    !**** *SWCLR* - CLEAR-SKY COLUMN COMPUTATIONS

    !     PURPOSE.
    !     --------
    !           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY IN CASE OF
    !     CLEAR-SKY COLUMN

    !**   INTERFACE.
    !     ----------

    !          *SWCLR* IS CALLED EITHER FROM *SW1S*
    !                                OR FROM *SWNI*

    !        IMPLICIT ARGUMENTS :
    !        --------------------

    !     ==== INPUTS ===
    !     ==== OUTPUTS ===

    !     METHOD.
    !     -------


    !     EXTERNALS.
    !     ----------

    !          NONE

    !     REFERENCE.
    !     ----------

    !        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
    !        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

    !     AUTHOR.
    !     -------
    !        JEAN-JACQUES MORCRETTE  *ECMWF*

    !     MODIFICATIONS.
    !     --------------
    !        ORIGINAL : 94-11-15
    !        Modified : 96-03-19 JJM-PhD (loop 107 in absence of aerosols)
    !        JJMorcrette 990128 : sunshine duration
    !        99-05-25   JJMorcrette    Revised aerosols
    !        June 2005: CCagnazzo from 4 to 6 bands, following 26Rcycle Morcrette
    !        Code
    !        08-10-28 : Marco Giorgetta: set PTAUAZ_N, PPIZAZ_N and PCGAZ_N
    !                   also if PPIZAZ<=EPSILON(1._DP)
    !     ------------------------------------------------------------------

    !     DUMMY INTEGER SCALARS
    INTEGER :: KFDIA
    INTEGER :: KIDIA
    INTEGER :: KLEV
    INTEGER :: KBDIM
    INTEGER :: KNU



    !     ------------------------------------------------------------------
    !*       0.1   ARGUMENTS
    !              ---------

    REAL(DP):: &
         &   PALB(KBDIM) &
         &,  PDSIG(KBDIM,KLEV)&
         &,  PRAYL(KBDIM)     &
         &,  PSEC(KBDIM)

    REAL(DP)::&
         &   PCGAZ(KBDIM,KLEV)     &
         &,  PPIZAZ(KBDIM,KLEV)&
         &,  PRAY1(KBDIM,KLEV+1)  , PRAY2(KBDIM,KLEV+1)&
         &,  PREFZ(KBDIM,2,KLEV+1), PRJ(KBDIM,6,KLEV+1)&
         &,  PRK(KBDIM,6,KLEV+1)  , PRMU0(KBDIM,KLEV+1)&
         &,  PTAUAZ(KBDIM,KLEV)&
         &,  PTRA1(KBDIM,KLEV+1)  , PTRA2(KBDIM,KLEV+1)&
         &,  PTRCLR(KBDIM)        , PDIFFUSE(KBDIM)

    ! aerosol properties
    REAL(DP) :: PSWEXT(KBDIM,KLEV,NSW), PSWSSA(KBDIM,KLEV,NSW) &
         &,        PSWASY(KBDIM,KLEV,NSW)

    !new set of optical properties without the delta transformation
    REAL(DP):: PTAUAZ_N(KBDIM,KLEV), PPIZAZ_N(KBDIM,KLEV), PCGAZ_N(KBDIM,KLEV)

    !     ------------------------------------------------------------------
    !*       0.2   LOCAL ARRAYS
    !              ------------

    REAL(DP):: ZC0I(KBDIM,KLEV+1)&
         &,  ZC0I2(KBDIM), ZCLEAR2(KBDIM) &
         &,  ZCLE0(KBDIM,KLEV), ZCLEAR(KBDIM) &
         &,  ZR21(KBDIM),  ZR23(KBDIM) , ZSS0(KBDIM) , ZSCAT(KBDIM)&
         &,  ZTR(KBDIM,2,KLEV+1)

    REAL(DP):: ZGG_1(KBDIM), ZTO1_1(KBDIM), ZW_1(KBDIM), ZREF(KBDIM)&
         &,       PR1(KBDIM), PT1(KBDIM), PR2(KBDIM), PT2(KBDIM)       &
         &,       ZRMUZ1(KBDIM), ZRMUZ(KBDIM)

    !     LOCAL INTEGER SCALARS
    INTEGER :: JAJ, JK, JKL, JKLP1, JKM1, JL

    !     LOCAL REAL SCALARS
    REAL(DP):: ZCORAE, ZFACOA, ZFF, ZMUE, ZRATIO, ZRE11, ZTRAY, PRMU1  &
         &, ZMUETO1, ZNORM

    !     ------------------------------------------------------------------
    !*         1.    OPTICAL PARAMETERS FOR AEROSOLS AND RAYLEIGH
    !                --------------------------------------------

    PRJ(:,:,:) = 0._DP
    PRK(:,:,:) = 0._DP

    DO JK = 1, KLEV
      DO JL = KIDIA,KFDIA
        PTAUAZ(JL,JK) = PSWEXT(JL,JK,KNU)
        PPIZAZ(JL,JK) = PSWSSA(JL,JK,KNU)
        PCGAZ(JL,JK)  = PSWASY(JL,JK,KNU)

        IF ( PPIZAZ(JL,JK)>EPSILON(1._DP) ) THEN
          PCGAZ(JL,JK)=PCGAZ(JL,JK)/PPIZAZ(JL,JK)
          PPIZAZ(JL,JK)=PPIZAZ(JL,JK)/PTAUAZ(JL,JK)

          ZTRAY = PRAYL(JL) * PDSIG(JL,JK)
          ZRATIO = ZTRAY / (ZTRAY + PTAUAZ(JL,JK))

          !new set of optical properties without delta transformation

          PTAUAZ_N(JL,JK)= ZTRAY + PTAUAZ(JL,JK)
          PCGAZ_N(JL,JK) = (PTAUAZ(JL,JK)*PPIZAZ(JL,JK)*PCGAZ(JL,JK))&
               &/(ZTRAY+PTAUAZ(JL,JK)*PPIZAZ(JL,JK))
          PPIZAZ_N(JL,JK)= (ZTRAY+PTAUAZ(JL,JK)*PPIZAZ(JL,JK))&
               &/PTAUAZ_N(JL,JK)

          !mixture of Rayleigh scatt and aerosols combined first by a weighted average
          !and then a delta transformation done to the mixture
          ZFF=PCGAZ_N(JL,JK)*PCGAZ_N(JL,JK) 
          PTAUAZ(JL,JK)=PTAUAZ_N(JL,JK)*(1._DP-PPIZAZ_N(JL,JK)*ZFF)
          PCGAZ(JL,JK)=PCGAZ_N(JL,JK)/(1._DP+PCGAZ_N(JL,JK))
          PPIZAZ(JL,JK) =(PPIZAZ_N(JL,JK)*(1._DP-ZFF))&
               &/ (1._DP - PPIZAZ_N(JL,JK) * ZFF)

        ELSE
          ZTRAY = PRAYL(JL) * PDSIG(JL,JK)
          PTAUAZ(JL,JK) = ZTRAY
          PCGAZ(JL,JK) = 0._DP
          PPIZAZ(JL,JK) = 1._DP-EPSILON(1._DP)
          PTAUAZ_N(JL,JK) = ZTRAY
          PCGAZ_N(JL,JK) = 0._DP
          PPIZAZ_N(JL,JK) =  1._DP-EPSILON(1._DP)
        ENDIF
      ENDDO
    ENDDO

    !     ------------------------------------------------------------------
    !*         2.    TOTAL EFFECTIVE CLOUDINESS ABOVE A GIVEN LEVEL
    !                ----------------------------------------------

    ZR23(:) = 0._DP
    ZC0I(:,KLEV+1) = 0._DP
    ZCLEAR(:) = 1._DP
    ZCLEAR2(:) = 1._DP
    ZSCAT(:) = 0._DP

    JK = 1
    JKL = KLEV+1 - JK
    JKLP1 = JKL + 1
    DO JL = KIDIA,KFDIA
      !repetition of applying the delta transformation to ptauaz
      ZFACOA = 1._DP 
      ZCORAE = ZFACOA * PTAUAZ(JL,JKL) * PSEC(JL)
      ZR21(JL) = EXP(-ZCORAE   )
      ZSS0(JL) = 1._DP-ZR21(JL)
      ZCLE0(JL,JKL) = ZSS0(JL)
    ENDDO

    IF (NOVLP == 1) THEN
      !* maximum-random      
      DO JL = KIDIA,KFDIA
        ZNORM = SWDIV_NOCHK(1._DP , (1._DP-MIN(ZSCAT(JL),1._DP-EPSILON(1._DP))))
        ZCLEAR(JL) = ZCLEAR(JL)&
             &*(1._DP-MAX(ZSS0(JL),ZSCAT(JL)))&
             &*ZNORM
        ZCLEAR2(JL) = ZCLEAR2(JL)*(1._DP-ZSS0(JL))
        ZC0I(JL,JKL) = 1._DP - ZCLEAR(JL)
        ZSCAT(JL) = ZSS0(JL)
      END DO
    ELSEIF (NOVLP == 2) THEN
      !* maximum
      DO JL = KIDIA,KFDIA
        ZSCAT(JL) = MAX( ZSS0(JL) , ZSCAT(JL) )
        ZCLEAR2(JL) = ZCLEAR2(JL)*(1._DP-ZSS0(JL))
        ZC0I(JL,JKL) = ZSCAT(JL)
      END DO
    ELSEIF (NOVLP == 3) THEN
      !* random
      DO JL = KIDIA,KFDIA
        ZCLEAR(JL)=ZCLEAR(JL)*(1._DP-ZSS0(JL))
        ZSCAT(JL) = 1._DP - ZCLEAR(JL)
        ZC0I(JL,JKL) = ZSCAT(JL)
      END DO
    ENDIF

    DO JK = 2 , KLEV
      JKL = KLEV+1 - JK
      JKLP1 = JKL + 1
      DO JL = KIDIA,KFDIA
        !repetition of applying the delta transformation to ptauaz
        ZFACOA = 1._DP 
        ZCORAE = ZFACOA * PTAUAZ(JL,JKL) * PSEC(JL)
        ZR21(JL) = EXP(-ZCORAE   )
        ZSS0(JL) = 1._DP-ZR21(JL)
        ZCLE0(JL,JKL) = ZSS0(JL)
        ZCLEAR2(JL) = ZCLEAR2(JL)*(1._DP-ZSS0(JL))
      END DO

      IF (NOVLP == 1) THEN
        !* maximum-random      
        DO JL = KIDIA,KFDIA
          ZNORM = SWDIV_NOCHK(1._DP , (1._DP-MIN(ZSCAT(JL),1._DP-EPSILON(1._DP))))
          ZCLEAR(JL) = ZCLEAR(JL)&
               &*(1._DP-MAX(ZSS0(JL),ZSCAT(JL)))&
               &*ZNORM
          ZC0I(JL,JKL) = 1._DP - ZCLEAR(JL)
          ZSCAT(JL) = ZSS0(JL)
          IF (JKL == 1) ZC0I2(JL) = 1._DP - ZCLEAR2(JL)
        END DO
      ELSEIF (NOVLP == 2) THEN
        !* maximum
        DO JL = KIDIA,KFDIA
          ZSCAT(JL) = MAX( ZSS0(JL) , ZSCAT(JL) )
          ZC0I(JL,JKL) = ZSCAT(JL)
          IF (JKL == 1) ZC0I2(JL) = 1._DP-ZCLEAR2(JL)
        END DO
      ELSEIF (NOVLP == 3) THEN
        !* random
        DO JL = KIDIA,KFDIA
          ZCLEAR(JL)=ZCLEAR(JL)*(1._DP-ZSS0(JL))
          ZSCAT(JL) = 1._DP - ZCLEAR(JL)
          ZC0I(JL,JKL) = ZSCAT(JL)
          IF (JKL == 1) ZC0I2(JL) = 1._DP - ZCLEAR(JL)
        END DO
      ENDIF
    ENDDO
    !     ------------------------------------------------------------------
    !*         3.    REFLECTIVITY/TRANSMISSIVITY FOR PURE SCATTERING
    !                -----------------------------------------------

    PRAY1(kidia:kfdia,KLEV+1) = 0._DP
    PRAY2(kidia:kfdia,KLEV+1) = 0._DP
    PREFZ(kidia:kfdia,2,1) = PALB(kidia:kfdia)
    PREFZ(kidia:kfdia,1,1) = PALB(kidia:kfdia)
    PTRA1(kidia:kfdia,KLEV+1) = 1._DP
    PTRA2(kidia:kfdia,KLEV+1) = 1._DP

    DO JK = 2 , KLEV+1
      JKM1 = JK-1
      DO JL = KIDIA,KFDIA
        PR1(JL) = 0._DP
        PT1(JL) = 0._DP
        PR2(JL) = 0._DP
        PT2(JL) = 0._DP

        !     ------------------------------------------------------------------
        !*         3.1  EQUIVALENT ZENITH ANGLE
        !               -----------------------


        ZMUE = (1._DP-ZC0I(JL,JK)) * PSEC(JL)+ ZC0I(JL,JK) * 1.66_DP
        PRMU0(JL,JK) = SWDIV_NOCHK(1._DP,ZMUE)
        PRMU1 = 0.5_DP

        !     ------------------------------------------------------------------
        !*         3.2  REFLECT./TRANSMISSIVITY DUE TO RAYLEIGH AND AEROSOLS
        !               ----------------------------------------------------
        !MANU--modified

        ZREF(JL) = PREFZ(JL,1,JKM1)
        ZRMUZ(JL) = PRMU0(JL,JK)
        ZRMUZ1(JL) = PRMU1
        ZGG_1(JL)=PCGAZ_N(JL,JKM1)
        ZTO1_1(JL)=PTAUAZ_N(JL,JKM1)
        ZW_1(JL) = PPIZAZ_N(JL,JKM1)

      ENDDO

      !pr1 and pt1 are reflectivity and transmissivity without reflection using
      !the equivalent zenith angle, prmu0

      CALL SWDE_WR ( KIDIA, KFDIA , KBDIM &
           &, ZGG_1  , ZREF  , ZRMUZ , ZTO1_1 , ZW_1 &
           &, PR1 , PT1 )

      !pr2 and pt2 are reflectivity and transmissivity without reflection using 
      !the fixed zenith angle, prmu1

      CALL SWDE_WR ( KIDIA, KFDIA , KBDIM &
           &, ZGG_1  , ZREF  , ZRMUZ1 , ZTO1_1 , ZW_1 &
           &, PR2 , PT2 )

      DO JL = KIDIA,KFDIA

        PRAY1(JL,JKM1) = PR1(JL)
        PTRA1(JL,JKM1) = PT1(JL)

        PRAY2(JL,JKM1) = PR2(JL)
        PTRA2(JL,JKM1) = PT2(JL)

        ZNORM = SWDIV_NOCHK(1._DP , (1._DP-PRAY2(JL,JKM1)*PREFZ(JL,1,JKM1)))
        PREFZ(JL,1,JK) = (PRAY1(JL,JKM1)&
             &+ PREFZ(JL,1,JKM1) * PTRA1(JL,JKM1)&
             &* PTRA2(JL,JKM1) * ZNORM)

        ZTR(JL,1,JKM1) = (PTRA1(JL,JKM1) * ZNORM)

        PREFZ(JL,2,JK) = (PRAY1(JL,JKM1)&
             &+ PREFZ(JL,2,JKM1) * PTRA1(JL,JKM1)&
             &* PTRA2(JL,JKM1) )

        ZTR(JL,2,JKM1) = PTRA1(JL,JKM1)

      ENDDO
    ENDDO

    DO JL = KIDIA,KFDIA
      ZMUE = (1._DP-ZC0I(JL,1))*PSEC(JL)+ZC0I(JL,1)*1.66_DP
      PRMU0(JL,1)=SWDIV_NOCHK(1._DP, ZMUE)
      PTRCLR(JL)=1._DP-ZC0I(JL,1)
      ZMUETO1 = SWDIV_NOCHK(1._DP, (1._DP - ZC0I2(JL) * 0.5_DP))
      PDIFFUSE(JL) = ZC0I2(JL) * 0.5_DP * ZMUETO1
    ENDDO


    !     ------------------------------------------------------------------
    !*         3.5    REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL
    !                 -------------------------------------------------


    IF (KNU <= 3) THEN 
      JAJ = 2
      DO JL = KIDIA,KFDIA
        PRJ(JL,JAJ,KLEV+1) = 1._DP
        PRK(JL,JAJ,KLEV+1) = PREFZ(JL, 1,KLEV+1)
      ENDDO

      DO JK = 1 , KLEV
        JKL = KLEV+1 - JK
        JKLP1 = JKL + 1
        DO JL = KIDIA,KFDIA
          ZRE11= PRJ(JL,JAJ,JKLP1) * ZTR(JL,  1,JKL)
          PRJ(JL,JAJ,JKL) = ZRE11
          PRK(JL,JAJ,JKL) = ZRE11 * PREFZ(JL,  1,JKL)
        ENDDO
      ENDDO

    ELSE

      DO JAJ = 1 , 2
        DO JL = KIDIA,KFDIA
          PRJ(JL,JAJ,KLEV+1) = 1._DP
          PRK(JL,JAJ,KLEV+1) = PREFZ(JL,JAJ,KLEV+1)
        ENDDO

        DO JK = 1 , KLEV
          JKL = KLEV+1 - JK
          JKLP1 = JKL + 1
          DO JL = KIDIA,KFDIA
            ZRE11= PRJ(JL,JAJ,JKLP1) * ZTR(JL,JAJ,JKL)
            PRJ(JL,JAJ,JKL) = ZRE11
            PRK(JL,JAJ,JKL) = ZRE11 * PREFZ(JL,JAJ,JKL)
          ENDDO
        ENDDO
      ENDDO

    ENDIF
    !     ------------------------------------------------------------------

    RETURN
  END SUBROUTINE SWCLR

  SUBROUTINE SWR &
       &( KIDIA , KFDIA , KBDIM , KLEV  , KNU &
       &, PALB , PCG   , PCLD , POMEGA, PSEC , PTAU &
       &, PCGAZ , PPIZAZ, PRAY1, PRAY2 , PREFZ, PRJ  , PRK , PRMUE &
       &, PTAUAZ, PTRA1 , PTRA2, PTRCLD &
       &, PDIFFUSE &
       &, PTAUAZ_N, PPIZAZ_N, PCGAZ_N)

    ! additional optical properties "*_N" -- without the delta transformations

    !**** *SWR* - CONTINUUM SCATTERING COMPUTATIONS

    !     PURPOSE.
    !     --------
    !           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY IN CASE OF
    !     CONTINUUM SCATTERING

    !**   INTERFACE.
    !     ----------

    !          *SWR* IS CALLED EITHER FROM *SW1S*
    !                              OR FROM *SWNI*

    !        IMPLICIT ARGUMENTS :
    !        --------------------

    !     ==== INPUTS ===
    !     ==== OUTPUTS ===

    !     METHOD.
    !     -------

    !          1. COMPUTES CONTINUUM FLUXES CORRESPONDING TO AEROSOL
    !     OR/AND RAYLEIGH SCATTERING (NO MOLECULAR GAS ABSORPTION)

    !     EXTERNALS.
    !     ----------

    !          *SWDE*

    !     REFERENCE.
    !     ----------

    !        SEE RADIATION'S PART OF THE ECMWF RESEARCH DEPARTMENT
    !        DOCUMENTATION, AND FOUQUART AND BONNEL (1980)

    !     AUTHOR.
    !     -------
    !        JEAN-JACQUES MORCRETTE  *ECMWF*

    !     MODIFICATIONS.
    !     --------------
    !        ORIGINAL : 89-07-14
    !        Ph. DANDIN Meteo-France 05-96 : Effect of cloud layer
    !        JJMorcrette 990128 : sunshine duration
    !        June 2005: CCagnazzo from 4 to 6 bands, following 26Rcycle Morcrette
    !        Code

    !     ------------------------------------------------------------------

    !     DUMMY INTEGER SCALARS
    INTEGER :: KFDIA
    INTEGER :: KIDIA
    INTEGER :: KLEV
    INTEGER :: KBDIM
    INTEGER :: KNU



    !     ------------------------------------------------------------------

    !*       0.1   ARGUMENTS
    !              ---------

    REAL(DP):: PALB(KBDIM)      , PCG(KBDIM,NSW,KLEV)&
         &,  PCLD(KBDIM,KLEV)&
         &,  POMEGA(KBDIM,NSW,KLEV)&
         &,  PSEC(KBDIM)           , PTAU(KBDIM,NSW,KLEV)

    REAL(DP):: PRAY1(KBDIM,KLEV+1)   , PRAY2(KBDIM,KLEV+1)&
         &,  PREFZ(KBDIM,2,KLEV+1) , PRJ(KBDIM,6,KLEV+1)&
         &,  PRK(KBDIM,6,KLEV+1)   , PRMUE(KBDIM,KLEV+1)&
         &,  PCGAZ(KBDIM,KLEV)     , PPIZAZ(KBDIM,KLEV)&
         &,  PTAUAZ(KBDIM,KLEV)&
         &,  PTRA1(KBDIM,KLEV+1)   , PTRA2(KBDIM,KLEV+1)&
         &,  PTRCLD(KBDIM)         , PDIFFUSE(KBDIM)

    !new optical properties without the delta transformations
    REAL(DP):: PTAUAZ_N(KBDIM,KLEV), PPIZAZ_N(KBDIM,KLEV), PCGAZ_N(KBDIM,KLEV)

    !     ------------------------------------------------------------------

    !*       0.2   LOCAL ARRAYS
    !              ------------

    REAL(DP):: ZC1I(KBDIM,KLEV+1)    , ZCLEQ(KBDIM,KLEV)&
         &,  ZC1I2(KBDIM)          , ZCLEAR2(KBDIM) &
         &,  ZCLOUD2(KBDIM) &
         &,  ZCLEAR(KBDIM)         , ZCLOUD(KBDIM) &
         &,  ZGG(KBDIM)            , ZREF(KBDIM)&
         &,  ZRE1(KBDIM)           , ZRE2(KBDIM)&
         &,  ZRMUZ(KBDIM)          , ZRNEB(KBDIM)&
         &,  ZR21(KBDIM)           , ZR22(KBDIM)&
         &,  ZR23(KBDIM)           , ZSS1(KBDIM)&
         &,  ZSS2(KBDIM)           , ZSS3(KBDIM)&
         &,  ZTO1(KBDIM)           , ZTR(KBDIM,2,KLEV+1)&
         &,  ZTR1(KBDIM)           , ZTR2(KBDIM)&
         &,  ZW(KBDIM)
    REAL(DP):: ZGG_1(KBDIM), ZTO1_1(KBDIM), ZW_1(KBDIM) &
         &,       PR1(KBDIM), PT1(KBDIM), PR2(KBDIM), PT2(KBDIM)&
         &,       ZRMUZ1(KBDIM)

    !     LOCAL INTEGER SCALARS
    INTEGER :: IKL, IKLP1, JA, JAJ, JK, JKM1, JL

    !     LOCAL REAL SCALARS
    REAL(DP):: ZBMU0, ZBMU1, ZCORAE, ZCORCD, ZDEN, ZDEN1,&
         &ZFACOA, ZFACOC, ZGAP, ZMU1, ZMUE, ZRE11, &
         &ZTO, ZWW, PRMU1, ZMUETO1, ZNORM, ZNORM2

    REAL(DP):: ZEPS

    !     ------------------------------------------------------------------

    !*         1.    INITIALIZATION
    !                --------------

    ZEPS=EPSILON(1.0_DP)

    PRJ(:,:,:) = 0._DP
    PRK(:,:,:) = 0._DP


    !     ------------------------------------------------------------------

    !*         2.    TOTAL EFFECTIVE CLOUDINESS ABOVE A GIVEN LEVEL
    !                ----------------------------------------------


    ZR23(:) = 0._DP
    ZC1I(:,KLEV+1) = 0._DP
    ZCLEAR(:) = 1._DP
    ZCLEAR2(:) = 1._DP
    ZCLOUD(:) = 0._DP
    ZCLOUD2(:) = 0._DP

    JK = 1
    IKL = KLEV+1 - JK
    IKLP1 = IKL + 1
    DO JL = KIDIA,KFDIA

      !avoiding the repetition of the delta transformation

      ZFACOA = 1._DP
      ZFACOC = 1._DP - POMEGA(JL,KNU,IKL) * PCG(JL,KNU,IKL)* PCG(JL,KNU,IKL)

      ZCORAE = ZFACOA * PTAUAZ(JL,IKL) * PSEC(JL)
      ZCORCD = ZFACOC * PTAU(JL,KNU,IKL) * PSEC(JL)
      ZCORAE = MIN(200._DP,ZCORAE)
      ZCORCD = MIN(200._DP,ZCORCD)
      ZR21(JL) = EXP(-ZCORAE   )
      ZR22(JL) = EXP(-ZCORCD   )
      ZSS1(JL) = PCLD(JL,IKL)*(1._DP-ZR21(JL)*ZR22(JL))&
           &+ (1._DP-PCLD(JL,IKL))*(1._DP-ZR21(JL))
      ZCLEQ(JL,IKL) = ZSS1(JL)
      ZSS2(JL) = 1._DP-ZR21(JL)
      ZSS3(JL) = PCLD(JL,IKL) * (1._DP-ZR22(JL))
    END DO

    IF (NOVLP == 1) THEN
      !* maximum-random      
      DO JL = KIDIA,KFDIA
        ZNORM = SWDIV_NOCHK(1._DP,(1._DP-MIN(ZCLOUD(JL),1._DP-EPSILON(1._DP))))
        ZNORM2= SWDIV_NOCHK(1._DP,(1._DP-MIN(ZCLOUD2(JL),1._DP-EPSILON(1._DP))))
        ZCLEAR(JL) = ZCLEAR(JL)&
             &*(1._DP-MAX(ZSS1(JL),ZCLOUD(JL)))&
             &*ZNORM
        ZCLEAR2(JL) = ZCLEAR2(JL)&
             &*(1._DP-MAX(ZSS3(JL),ZCLOUD2(JL)))&
             &*ZNORM2
        ZCLEAR2(JL) = ZCLEAR2(JL)*(1._DP - ZSS2(JL))
        ZC1I(JL,IKL) = 1._DP - ZCLEAR(JL)
        ZCLOUD(JL) = ZSS1(JL)
        ZCLOUD2(JL) = ZSS3(JL)
      END DO
    ELSEIF (NOVLP == 2) THEN
      !* maximum
      DO JL = KIDIA,KFDIA
        ZCLOUD(JL) = MAX( ZSS1(JL) , ZCLOUD(JL) )
        ZCLOUD2(JL) = MAX( ZSS3(JL) , ZCLOUD2(JL) )
        ZCLEAR2(JL) = ZCLEAR2(JL)*(1._DP - ZSS2(JL))
        ZC1I(JL,IKL) = ZCLOUD(JL)
      END DO
    ELSEIF (NOVLP == 3) THEN
      !* random
      DO JL = KIDIA,KFDIA
        ZCLEAR(JL) = ZCLEAR(JL)*(1._DP - ZSS1(JL))
        ZCLOUD(JL) = 1._DP - ZCLEAR(JL)
        ZC1I(JL,IKL) = ZCLOUD(JL)
      END DO
    ENDIF

    DO JK = 2 , KLEV
      IKL = KLEV+1 - JK
      IKLP1 = IKL + 1
      DO JL = KIDIA,KFDIA

        !avoiding the repetition of the delta transformation

        !    ZFACOA = 1._DP - PPIZAZ(JL,IKL)*PCGAZ(JL,IKL)*PCGAZ(JL,IKL)

        ZFACOA = 1._DP

        ZFACOC = 1._DP - POMEGA(JL,KNU,IKL) * PCG(JL,KNU,IKL)* PCG(JL,KNU,IKL)
        ZCORAE = ZFACOA * PTAUAZ(JL,IKL) * PSEC(JL)
        ZCORCD = ZFACOC * PTAU(JL,KNU,IKL) * PSEC(JL)
        ZCORAE = MIN(200._DP,ZCORAE)
        ZCORCD = MIN(200._DP,ZCORCD)
        ZR21(JL) = EXP(-ZCORAE   )
        ZR22(JL) = EXP(-ZCORCD   )
        ZSS1(JL) = PCLD(JL,IKL)*(1._DP-ZR21(JL)*ZR22(JL))&
             &+ (1._DP-PCLD(JL,IKL))*(1._DP-ZR21(JL))
        ZCLEQ(JL,IKL) = ZSS1(JL)
        ZSS2(JL) = 1._DP-ZR21(JL)
        ZSS3(JL) = PCLD(JL,IKL) * (1._DP-ZR22(JL))
      END DO

      IF (NOVLP == 1) THEN
        !* maximum-random      
        DO JL = KIDIA,KFDIA
          ZNORM = SWDIV_NOCHK(1._DP,(1._DP-MIN(ZCLOUD(JL),1._DP-EPSILON(1._DP))))
          ZNORM2= SWDIV_NOCHK(1._DP,(1._DP-MIN(ZCLOUD2(JL),1._DP-EPSILON(1._DP))))
          ZCLEAR(JL) = ZCLEAR(JL)&
               &*(1._DP-MAX(ZSS1(JL),ZCLOUD(JL)))&
               &*ZNORM
          ZCLEAR2(JL) = ZCLEAR2(JL)&
               &*(1._DP-MAX(ZSS3(JL),ZCLOUD2(JL)))&
               &*ZNORM2
          ZCLEAR2(JL) = ZCLEAR2(JL)*(1._DP - ZSS2(JL))
          ZC1I(JL,IKL) = 1._DP - ZCLEAR(JL)
!!$         ZC1I(JL,IKL) = 1._DP - ZCLEAR2(JL)
          IF (IKL == 1) ZC1I2(JL) = 1._DP - ZCLEAR2(JL)
          ZCLOUD(JL) = ZSS1(JL)
          ZCLOUD2(JL) = ZSS3(JL)
        END DO
      ELSEIF (NOVLP == 2) THEN
        !* maximum
        DO JL = KIDIA,KFDIA
          ZCLOUD(JL) = MAX( ZSS1(JL) , ZCLOUD(JL) )
          ZCLOUD2(JL) = MAX( ZSS3(JL) , ZCLOUD2(JL) )
          ZCLEAR2(JL) = ZCLEAR2(JL)*(1._DP - ZSS2(JL))
          ZC1I(JL,IKL) = ZCLOUD(JL)
          IF (IKL == 1) ZC1I2(JL) = 1._DP - ((1._DP - ZCLOUD2(JL)) * ZCLEAR2(JL))
        END DO
      ELSEIF (NOVLP == 3) THEN
        !* random
        DO JL = KIDIA,KFDIA
          ZCLEAR(JL) = ZCLEAR(JL)*(1._DP - ZSS1(JL))
          ZCLOUD(JL) = 1._DP - ZCLEAR(JL)
          ZC1I(JL,IKL) = ZCLOUD(JL)
          IF (IKL == 1) ZC1I2(JL) = 1._DP - ZCLEAR(JL)
        END DO
      ENDIF
    ENDDO

    !     ------------------------------------------------------------------

    !*         3.    REFLECTIVITY/TRANSMISSIVITY FOR PURE SCATTERING
    !                -----------------------------------------------


    PRAY1(kidia:kfdia,KLEV+1) = 0._DP
    PRAY2(kidia:kfdia,KLEV+1) = 0._DP
    PREFZ(kidia:kfdia,2,1) = PALB(kidia:kfdia)
    PREFZ(kidia:kfdia,1,1) = PALB(kidia:kfdia)
    PTRA1(kidia:kfdia,KLEV+1) = 1._DP
    PTRA2(kidia:kfdia,KLEV+1) = 1._DP

    DO JK = 2 , KLEV+1
      JKM1 = JK-1
!IBM* NOVECTOR
      DO JL = KIDIA,KFDIA
        PR1(JL) = 0._DP
        PT1(JL) = 0._DP
        PR2(JL) = 0._DP
        PT2(JL) = 0._DP


        !     ------------------------------------------------------------------

        !*         3.1  EQUIVALENT ZENITH ANGLE
        !               -----------------------


        ZMUE = (1._DP-ZC1I(JL,JK)) * PSEC(JL)+ ZC1I(JL,JK) * 1.66_DP
        PRMUE(JL,JK) = SWDIV_NOCHK(1._DP , ZMUE)
        PRMU1 = 0.5_DP


        !     ------------------------------------------------------------------

        !*         3.2  REFLECT./TRANSMISSIVITY DUE TO RAYLEIGH AND AEROSOLS
        !               ----------------------------------------------------

        ZREF(JL) = PREFZ(JL,1,JKM1)
        ZRMUZ(JL) = PRMUE(JL,JK)
        ZRMUZ1(JL) = PRMU1
        ZGG_1(JL)=PCGAZ_N(JL,JKM1)
        ZTO1_1(JL)=PTAUAZ_N(JL,JKM1)
        ZW_1(JL) = PPIZAZ_N(JL,JKM1)

      ENDDO

      !pr1 and pt1 are reflectivity and transmissivity without reflection using the
      !equivalent zenith angle, prmu0

      CALL SWDE_WR ( KIDIA, KFDIA , KBDIM &
           &, ZGG_1  , ZREF  , ZRMUZ , ZTO1_1 , ZW_1 &
           &, PR1 , PT1 )

      !pr2 and pt2 are reflectivity and transmissivity without reflection using the
      !fixed zenith angle, prmu1

      CALL SWDE_WR ( KIDIA, KFDIA , KBDIM &
           &, ZGG_1  , ZREF  , ZRMUZ1 , ZTO1_1 , ZW_1 &
           &, PR2 , PT2 )

      !     ------------------------------------------------------------------

      !*         3.3  EFFECT OF CLOUD LAYER
      !               ---------------------

      DO JL = KIDIA,KFDIA
        ZRNEB(JL)= PCLD(JL,JKM1)
        ZRE1(JL)=0._DP
        ZTR1(JL)=0._DP
        ZRE2(JL)=0._DP
        ZTR2(JL)=0._DP

        ZTO1(JL) =  PTAU(JL,KNU,JKM1) &
             &       +PTAUAZ_N(JL,JKM1)
        ZW(JL)   =  PTAUAZ_N(JL,JKM1) *PPIZAZ_N(JL,JKM1) &
             &       +PTAU(JL,KNU,JKM1) *POMEGA(JL,KNU,JKM1)
        ZGG(JL)  =  PTAUAZ_N(JL,JKM1)*PPIZAZ_N(JL,JKM1)*PCGAZ_N(JL,JKM1) &
             &       +PTAU(JL,KNU,JKM1)*POMEGA(JL,KNU,JKM1)*PCG(JL,KNU,JKM1)
        IF(ZW(JL)> ZEPS) THEN
          ZGG(JL)  =  ZGG(JL)/ZW(JL)
        ELSE
          ZGG(JL)  = 0.0_DP
        ENDIF

        IF(ZTO1(JL)>ZEPS) THEN
          ZW(JL)   =  ZW(JL)/ZTO1(JL)
        ELSE
          ZW(JL)   = 1.0_DP
        ENDIF

        ZW(JL)=MIN(ZW(JL), 1.0_DP-ZEPS)


        ZREF(JL) = PREFZ(JL,1,JKM1)
        ZRMUZ(JL) = PRMUE(JL,JK)

      ENDDO

      CALL SWDE ( KIDIA, KFDIA , KBDIM &
           &, ZGG  , ZREF  , ZRMUZ , ZTO1 , ZW &
           &, ZRE1 , ZRE2  , ZTR1  , ZTR2      )

!IBM* NOVECTOR
      DO JL = KIDIA,KFDIA

        PRAY1(JL,JKM1) = PR1(JL)
        PTRA1(JL,JKM1) = PT1(JL)
        PRAY2(JL,JKM1) = PR2(JL)
        PTRA2(JL,JKM1) = PT2(JL)

        ZNORM = SWDIV_NOCHK(1._DP , (1._DP-PRAY2(JL,JKM1)*PREFZ(JL,1,JKM1)))

        PREFZ(JL,1,JK) = (1._DP-ZRNEB(JL)) * (PRAY1(JL,JKM1)&
             &+ PREFZ(JL,1,JKM1) * PTRA1(JL,JKM1)&
             &* PTRA2(JL,JKM1) * ZNORM)&
             &+ ZRNEB(JL) * ZRE2(JL)

        ZTR(JL,1,JKM1) = ZRNEB(JL) * ZTR2(JL) + PTRA1(JL,JKM1)&
             &* ZNORM&
             &* (1._DP-ZRNEB(JL))

        PREFZ(JL,2,JK) = (1._DP-ZRNEB(JL)) * (PRAY1(JL,JKM1)&
             &+ PREFZ(JL,2,JKM1) * PTRA1(JL,JKM1)&
             &* PTRA2(JL,JKM1) )&
             &+ ZRNEB(JL) * ZRE1(JL)

        ZTR(JL,2,JKM1) = ZRNEB(JL) * ZTR1(JL)+ PTRA1(JL,JKM1) * (1._DP-ZRNEB(JL))

      ENDDO
    ENDDO
    DO JL = KIDIA,KFDIA
      ZMUE = (1._DP-ZC1I(JL,1))*PSEC(JL)+ZC1I(JL,1)*1.66_DP
      PRMUE(JL,1)=SWDIV_NOCHK(1._DP,ZMUE)
      PTRCLD(JL)=1._DP-ZC1I(JL,1)
      ZMUETO1 = SWDIV_NOCHK(1._DP, (1._DP - ZC1I2(JL) * 0.5_DP))
      PDIFFUSE(JL) = ZC1I2(JL) * 0.5_DP * ZMUETO1
    ENDDO


    !     ------------------------------------------------------------------

    !*         3.5    REFLECT./TRANSMISSIVITY BETWEEN SURFACE AND LEVEL
    !                 -------------------------------------------------

    IF (KNU <= 3) THEN   
      JAJ = 2
      DO JL = KIDIA,KFDIA
        PRJ(JL,JAJ,KLEV+1) = 1._DP
        PRK(JL,JAJ,KLEV+1) = PREFZ(JL, 1,KLEV+1)
      ENDDO

      DO JK = 1 , KLEV
        IKL = KLEV+1 - JK
        IKLP1 = IKL + 1
        DO JL = KIDIA,KFDIA
          ZRE11= PRJ(JL,JAJ,IKLP1) * ZTR(JL,  1,IKL)
          PRJ(JL,JAJ,IKL) = ZRE11
          PRK(JL,JAJ,IKL) = ZRE11 * PREFZ(JL,  1,IKL)
        ENDDO
      ENDDO

    ELSE

      DO JAJ = 1 , 2
        DO JL = KIDIA,KFDIA
          PRJ(JL,JAJ,KLEV+1) = 1._DP
          PRK(JL,JAJ,KLEV+1) = PREFZ(JL,JAJ,KLEV+1)
        ENDDO

        DO JK = 1 , KLEV
          IKL = KLEV+1 - JK
          IKLP1 = IKL + 1
          DO JL = KIDIA,KFDIA
            ZRE11= PRJ(JL,JAJ,IKLP1) * ZTR(JL,JAJ,IKL)
            PRJ(JL,JAJ,IKL) = ZRE11
            PRK(JL,JAJ,IKL) = ZRE11 * PREFZ(JL,JAJ,IKL)
          ENDDO
        ENDDO
      ENDDO

    ENDIF

    !     ------------------------------------------------------------------

    RETURN
  END SUBROUTINE SWR

  SUBROUTINE SWDE &
       &( KIDIA, KFDIA, KBDIM &
       &, PGG  , PREF , PRMUZ, PTO1, PW &
       &, PRE1 , PRE2 , PTR1 , PTR2 &
       &)

    !**** *SWDE* - DELTA-EDDINGTON IN A CLOUDY LAYER

    !     PURPOSE.
    !     --------
    !           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY OF A CLOUDY
    !     LAYER USING THE DELTA-EDDINGTON'S APPROXIMATION.

    !**   INTERFACE.
    !     ----------
    !          *SWDE* IS CALLED BY *SWR*, *SWNI*


    !        EXPLICIT ARGUMENTS :
    !        --------------------
    ! PGG    : (KBDIM)             ; ASSYMETRY FACTOR
    ! PREF   : (KBDIM)             ; REFLECTIVITY OF THE UNDERLYING LAYER
    ! PRMUZ  : (KBDIM)             ; COSINE OF SOLAR ZENITH ANGLE
    ! PTO1   : (KBDIM)             ; OPTICAL THICKNESS
    ! PW     : (KBDIM)             ; SINGLE SCATTERING ALBEDO
    !     ==== OUTPUTS ===
    ! PRE1   : (KBDIM)             ; LAYER REFLECTIVITY ASSUMING NO
    !                             ; REFLECTION FROM UNDERLYING LAYER
    ! PTR1   : (KBDIM)             ; LAYER TRANSMISSIVITY ASSUMING NO
    !                             ; REFLECTION FROM UNDERLYING LAYER
    ! PRE2   : (KBDIM)             ; LAYER REFLECTIVITY ASSUMING
    !                             ; REFLECTION FROM UNDERLYING LAYER
    ! PTR2   : (KBDIM)             ; LAYER TRANSMISSIVITY ASSUMING
    !                             ; REFLECTION FROM UNDERLYING LAYER

    !        IMPLICIT ARGUMENTS :   NONE
    !        --------------------

    !     METHOD.
    !     -------

    !          STANDARD DELTA-EDDINGTON LAYER CALCULATIONS.

    !     EXTERNALS.
    !     ----------

    !          NONE

    !     REFERENCE.
    !     ----------

    !        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
    !        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

    !     AUTHOR.
    !     -------
    !        JEAN-JACQUES MORCRETTE  *ECMWF*

    !     MODIFICATIONS.
    !     --------------
    !        ORIGINAL : 88-12-15
    !                   96-05-30 Michel Deque (security in EXP())
    !                   08-10-17 Marco Giorgetta (PW<=1-EPSILON)

    !     ------------------------------------------------------------------


    !     ------------------------------------------------------------------

    !*       0.1   ARGUMENTS
    !              ---------

    !     DUMMY INTEGER SCALARS
    INTEGER :: KFDIA
    INTEGER :: KIDIA
    INTEGER :: KBDIM

    REAL(DP):: PGG(KBDIM),PREF(KBDIM),PRMUZ(KBDIM),PTO1(KBDIM),PW(KBDIM)
    REAL(DP):: PRE1(KBDIM),PRE2(KBDIM),PTR1(KBDIM),PTR2(KBDIM)

    !     LOCAL INTEGER SCALARS
    INTEGER :: JL

    !     LOCAL REAL SCALARS
    REAL(DP):: ZA11, ZA12, ZA13, ZA21, ZA22, ZA23, ZALPHA,&
         &ZAM2B, ZAP2B, ZARG, ZARG2, ZB21, ZB22, ZB23, &
         &ZBETA, ZC1A, ZC1B, ZC2A, ZC2B, ZDENA, ZDENB, &
         &ZDT, ZEXKM, ZEXKP, ZEXMU0, ZFF, ZGP, ZRI0A, &
         &ZRI0B, ZRI0C, ZRI0D, ZRI1A, ZRI1B, ZRI1C, &
         &ZRI1D, ZRK, ZRM2, ZRP, ZTOP, ZWCP, ZWM, ZX1, &
         &ZX2, ZXM2P, ZXP2P



    !     ------------------------------------------------------------------

    !*         1.      DELTA-EDDINGTON CALCULATIONS


    DO JL   =   KIDIA,KFDIA

      !*         1.0     PW must be smaller than 1

      PW(JL) = MIN(PW(JL),1._DP-EPSILON(1._DP))

      !*         1.1     SET UP THE DELTA-MODIFIED PARAMETERS


      ZFF = PGG(JL)*PGG(JL)
      ZGP = PGG(JL)/(1._DP+PGG(JL))
      ZTOP = (1._DP- PW(JL) * ZFF) * PTO1(JL)
      ZWCP = (1-ZFF)* PW(JL) /(1._DP- PW(JL) * ZFF)
      ZDT = 2._DP/3._DP
      ZX1 = 1._DP-ZWCP*ZGP
      ZWM = 1._DP-ZWCP
      ZRM2 =  PRMUZ(JL) * PRMUZ(JL)
      ZRK = SQRT(3._DP*ZWM*ZX1)
      ZX2 = 4._DP*(1._DP-ZRK*ZRK*ZRM2)
      ZRP=ZRK/ZX1
      ZALPHA = 3._DP*ZWCP*ZRM2*(1._DP+ZGP*ZWM)/ZX2
      ZBETA = 3._DP*ZWCP* PRMUZ(JL) *(1._DP+3._DP*ZGP*ZRM2*ZWM)/ZX2
      !     ZARG=MIN(ZTOP/PRMUZ(JL),200.)
      ZARG=MAX(-200._DP,MIN(ZTOP/PRMUZ(JL),200._DP))
      ZEXMU0=EXP(-ZARG)
      ZARG2=MIN(ZRK*ZTOP,200._DP)
      ZEXKP=EXP(ZARG2)
      ZEXKM = 1._DP/ZEXKP
      ZXP2P = 1._DP+ZDT*ZRP
      ZXM2P = 1._DP-ZDT*ZRP
      ZAP2B = ZALPHA+ZDT*ZBETA
      ZAM2B = ZALPHA-ZDT*ZBETA

      !*         1.2     WITHOUT REFLECTION FROM THE UNDERLYING LAYER


      ZA11 = ZXP2P
      ZA12 = ZXM2P
      ZA13 = ZAP2B
      ZA22 = ZXP2P*ZEXKP
      ZA21 = ZXM2P*ZEXKM
      ZA23 = ZAM2B*ZEXMU0
      ZDENA = ZA11 * ZA22 - ZA21 * ZA12
      ZC1A = (ZA22*ZA13-ZA12*ZA23)/ZDENA
      ZC2A = (ZA11*ZA23-ZA21*ZA13)/ZDENA
      ZRI0A = ZC1A+ZC2A-ZALPHA
      ZRI1A = ZRP*(ZC1A-ZC2A)-ZBETA
      PRE1(JL) = (ZRI0A-ZDT*ZRI1A)/ PRMUZ(JL)
      ZRI0B = ZC1A*ZEXKM+ZC2A*ZEXKP-ZALPHA*ZEXMU0
      ZRI1B = ZRP*(ZC1A*ZEXKM-ZC2A*ZEXKP)-ZBETA*ZEXMU0
      PTR1(JL) = ZEXMU0+(ZRI0B+ZDT*ZRI1B)/ PRMUZ(JL)

      !*         1.3     WITH REFLECTION FROM THE UNDERLYING LAYER


      ZB21 = ZA21- PREF(JL) *ZXP2P*ZEXKM
      ZB22 = ZA22- PREF(JL) *ZXM2P*ZEXKP
      ZB23 = ZA23- PREF(JL) *ZEXMU0*(ZAP2B - PRMUZ(JL) )
      ZDENB = ZA11 * ZB22 - ZB21 * ZA12
      ZC1B = (ZB22*ZA13-ZA12*ZB23)/ZDENB
      ZC2B = (ZA11*ZB23-ZB21*ZA13)/ZDENB
      ZRI0C = ZC1B+ZC2B-ZALPHA
      ZRI1C = ZRP*(ZC1B-ZC2B)-ZBETA
      PRE2(JL) = (ZRI0C-ZDT*ZRI1C) / PRMUZ(JL)
      ZRI0D = ZC1B*ZEXKM + ZC2B*ZEXKP - ZALPHA*ZEXMU0
      ZRI1D = ZRP * (ZC1B*ZEXKM - ZC2B*ZEXKP) - ZBETA*ZEXMU0
      PTR2(JL) = ZEXMU0 + (ZRI0D + ZDT*ZRI1D) / PRMUZ(JL)

    ENDDO
    RETURN
  END SUBROUTINE SWDE

  SUBROUTINE SWDE_WR &
       &( KIDIA, KFDIA, KBDIM &
       &, PGG  , PREF , PRMUZ, PTO1, PW &
       &, PRE1 , PTR1 &
       &)

    !**** *SWDE_WR* - DELTA-EDDINGTON IN A CLOUDY LAYER

    !     PURPOSE.
    !     --------
    !           COMPUTES THE REFLECTIVITY AND TRANSMISSIVITY OF A CLOUDY
    !     LAYER USING THE DELTA-EDDINGTON'S APPROXIMATION.

    !**   INTERFACE.
    !     ----------
    !          *SWDE* IS CALLED BY *SWR*, *SWNI*


    !        EXPLICIT ARGUMENTS :
    !        --------------------
    ! PGG    : (KBDIM)             ; ASSYMETRY FACTOR
    ! PREF   : (KBDIM)             ; REFLECTIVITY OF THE UNDERLYING LAYER
    ! PRMUZ  : (KBDIM)             ; COSINE OF SOLAR ZENITH ANGLE
    ! PTO1   : (KBDIM)             ; OPTICAL THICKNESS
    ! PW     : (KBDIM)             ; SINGLE SCATTERING ALBEDO
    !     ==== OUTPUTS ===
    ! PRE1   : (KBDIM)             ; LAYER REFLECTIVITY ASSUMING NO
    !                             ; REFLECTION FROM UNDERLYING LAYER
    ! PTR1   : (KBDIM)             ; LAYER TRANSMISSIVITY ASSUMING NO
    !                             ; REFLECTION FROM UNDERLYING LAYER

    !        IMPLICIT ARGUMENTS :   NONE
    !        --------------------

    !     METHOD.
    !     -------

    !          STANDARD DELTA-EDDINGTON LAYER CALCULATIONS.

    !     EXTERNALS.
    !     ----------

    !          NONE

    !     REFERENCE.
    !     ----------

    !        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
    !        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

    !     AUTHOR.
    !     -------
    !        JEAN-JACQUES MORCRETTE  *ECMWF*

    !     MODIFICATIONS.
    !     --------------
    !        ORIGINAL : 88-12-15
    !                   96-05-30 Michel Deque (security in EXP())
    !                   08-10-17 Marco Giorgetta (PW<=1-EPSILON)

    !     ------------------------------------------------------------------


    !     ------------------------------------------------------------------

    !*       0.1   ARGUMENTS
    !              ---------

    !     DUMMY INTEGER SCALARS
    INTEGER :: KFDIA
    INTEGER :: KIDIA
    INTEGER :: KBDIM

    REAL(DP):: PGG(KBDIM),PREF(KBDIM),PRMUZ(KBDIM),PTO1(KBDIM),PW(KBDIM)
    REAL(DP):: PRE1(KBDIM),PTR1(KBDIM)

    !     LOCAL INTEGER SCALARS
    INTEGER :: JL

    !     LOCAL REAL SCALARS
    REAL(DP):: ZA11, ZA12, ZA13, ZA21, ZA22, ZA23, ZALPHA,&
         &ZAM2B, ZAP2B, ZARG, ZARG2, &
         &ZBETA, ZC1A, ZC2A, ZDENA,  &
         &ZDT, ZEXKM, ZEXKP, ZEXMU0, ZFF, ZGP, ZRI0A, &
         &ZRI0B, ZRI1A, ZRI1B, &
         &ZRK, ZRM2, ZRP, ZTOP, ZWCP, ZWM, ZX1, &
         &ZX2, ZXM2P, ZXP2P



    !     ------------------------------------------------------------------

    !*         1.      DELTA-EDDINGTON CALCULATIONS


    DO JL   =   KIDIA,KFDIA

      !*         1.0     PW must be smaller than 1

      PW(JL) = MIN(PW(JL),1._DP-EPSILON(1._DP))

      !*         1.1     SET UP THE DELTA-MODIFIED PARAMETERS


      ZFF = PGG(JL)*PGG(JL)
      ZGP = PGG(JL)/(1._DP+PGG(JL))
      ZTOP = (1._DP- PW(JL) * ZFF) * PTO1(JL)
      ZWCP = (1-ZFF)* PW(JL) /(1._DP- PW(JL) * ZFF)
      ZDT = 2._DP/3._DP
      ZX1 = 1._DP-ZWCP*ZGP
      ZWM = 1._DP-ZWCP
      ZRM2 =  PRMUZ(JL) * PRMUZ(JL)
      ZRK = SQRT(3._DP*ZWM*ZX1)
      ZX2 = 4._DP*(1._DP-ZRK*ZRK*ZRM2)
      ZRP=ZRK/ZX1
      ZALPHA = 3._DP*ZWCP*ZRM2*(1._DP+ZGP*ZWM)/ZX2
      ZBETA = 3._DP*ZWCP* PRMUZ(JL) *(1._DP+3._DP*ZGP*ZRM2*ZWM)/ZX2
      !     ZARG=MIN(ZTOP/PRMUZ(JL),200.)
      ZARG=MAX(-200._DP,MIN(ZTOP/PRMUZ(JL),200._DP))
      ZEXMU0=EXP(-ZARG)
      ZARG2=MIN(ZRK*ZTOP,200._DP)
      ZEXKP=EXP(ZARG2)
      ZEXKM = 1._DP/ZEXKP
      ZXP2P = 1._DP+ZDT*ZRP
      ZXM2P = 1._DP-ZDT*ZRP
      ZAP2B = ZALPHA+ZDT*ZBETA
      ZAM2B = ZALPHA-ZDT*ZBETA

      !*         1.2     WITHOUT REFLECTION FROM THE UNDERLYING LAYER


      ZA11 = ZXP2P
      ZA12 = ZXM2P
      ZA13 = ZAP2B
      ZA22 = ZXP2P*ZEXKP
      ZA21 = ZXM2P*ZEXKM
      ZA23 = ZAM2B*ZEXMU0
      ZDENA = ZA11 * ZA22 - ZA21 * ZA12
      ZC1A = (ZA22*ZA13-ZA12*ZA23)/ZDENA
      ZC2A = (ZA11*ZA23-ZA21*ZA13)/ZDENA
      ZRI0A = ZC1A+ZC2A-ZALPHA
      ZRI1A = ZRP*(ZC1A-ZC2A)-ZBETA
      PRE1(JL) = (ZRI0A-ZDT*ZRI1A)/ PRMUZ(JL)
      ZRI0B = ZC1A*ZEXKM+ZC2A*ZEXKP-ZALPHA*ZEXMU0
      ZRI1B = ZRP*(ZC1A*ZEXKM-ZC2A*ZEXKP)-ZBETA*ZEXMU0
      PTR1(JL) = ZEXMU0+(ZRI0B+ZDT*ZRI1B)/ PRMUZ(JL)

    ENDDO
    RETURN
  END SUBROUTINE SWDE_WR

  SUBROUTINE SWTT ( KIDIA, KFDIA, KBDIM, KNU, KA , PU, PTR)

    !**** *SWTT* - COMPUTES THE SHORTWAVE TRANSMISSION FUNCTIONS

    !     PURPOSE.
    !     --------
    !           THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR ALL THE
    !     ABSORBERS (H2O, UNIFORMLY MIXED GASES, AND O3) IN THE TWO SPECTRAL
    !     INTERVALS.

    !**   INTERFACE.
    !     ----------
    !          *SWTT* IS CALLED FROM *SW1S*, *SWNI*.


    !        EXPLICIT ARGUMENTS :
    !        --------------------
    ! KNU    :                     ; INDEX OF THE SPECTRAL INTERVAL
    ! KA     :                     ; INDEX OF THE ABSORBER
    ! PU     : (KBDIM)             ; ABSORBER AMOUNT
    !     ==== OUTPUTS ===
    ! PTR    : (KBDIM)             ; TRANSMISSION FUNCTION

    !        IMPLICIT ARGUMENTS :   NONE
    !        --------------------

    !     METHOD.
    !     -------

    !          TRANSMISSION FUNCTION ARE COMPUTED USING PADE APPROXIMANTS
    !     AND HORNER'S ALGORITHM.

    !     EXTERNALS.
    !     ----------

    !          NONE

    !     REFERENCE.
    !     ----------

    !        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
    !        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

    !     AUTHOR.
    !     -------
    !        JEAN-JACQUES MORCRETTE  *ECMWF*

    !     MODIFICATIONS.
    !     --------------
    !        ORIGINAL : 88-12-15

    !-----------------------------------------------------------------------

    !     DUMMY INTEGER SCALARS
    INTEGER :: KA
    INTEGER :: KFDIA
    INTEGER :: KIDIA
    INTEGER :: KBDIM
    INTEGER :: KNU



    !-----------------------------------------------------------------------

    !*       0.1   ARGUMENTS
    !              ---------

    REAL(DP):: PU(KBDIM), PTR(KBDIM)

    !-----------------------------------------------------------------------

    !*       0.2   LOCAL ARRAYS
    !              ------------

    REAL(DP):: ZR1(KBDIM), ZR2(KBDIM)

    !     LOCAL INTEGER SCALARS
    INTEGER :: JL


    !-----------------------------------------------------------------------

    !*         1.      HORNER'S ALGORITHM TO COMPUTE TRANSMISSION FUNCTION


!IBM* novector
    DO JL = KIDIA,KFDIA
      ZR1(JL) = APAD(KNU,KA,1) + PU(JL) * (APAD(KNU,KA,2) + PU(JL)&
           &* ( APAD(KNU,KA,3) + PU(JL) * (APAD(KNU,KA,4) + PU(JL)&
           &* ( APAD(KNU,KA,5) + PU(JL) * (APAD(KNU,KA,6) + PU(JL)&
           &* ( APAD(KNU,KA,7) ))))))

      ZR2(JL) = BPAD(KNU,KA,1) + PU(JL) * (BPAD(KNU,KA,2) + PU(JL)&
           &* ( BPAD(KNU,KA,3) + PU(JL) * (BPAD(KNU,KA,4) + PU(JL)&
           &* ( BPAD(KNU,KA,5) + PU(JL) * (BPAD(KNU,KA,6) + PU(JL)&
           &* ( BPAD(KNU,KA,7) ))))))


      !*         2.      ADD THE BACKGROUND TRANSMISSION



      PTR(JL) = (ZR1(JL) / ZR2(JL)) * (1._DP - D(KNU,KA)) + D(KNU,KA)
    ENDDO

    RETURN
  END SUBROUTINE SWTT

  SUBROUTINE SWTT1 ( KIDIA,KFDIA,KBDIM,KNU,KABS,KKIND, PU, PTR )

    !**** *SWTT1* - COMPUTES THE SHORTWAVE TRANSMISSION FUNCTIONS

    !     PURPOSE.
    !     --------
    !           THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR ALL THE
    !     ABSORBERS (H2O, UNIFORMLY MIXED GASES, AND O3) IN THE TWO SPECTRAL
    !     INTERVALS.

    !**   INTERFACE.
    !     ----------
    !          *SWTT1* IS CALLED FROM *SW1S*.


    !        EXPLICIT ARGUMENTS :
    !        --------------------
    ! KNU    :                     ; INDEX OF THE SPECTRAL INTERVAL
    ! KABS   :                     ; NUMBER OF ABSORBERS
    ! KKIND  : (KABS)              ; INDICES OF THE ABSORBERS
    ! PU     : (KBDIM,KABS)         ; ABSORBER AMOUNT
    !     ==== OUTPUTS ===
    ! PTR    : (KBDIM,KABS)         ; TRANSMISSION FUNCTION

    !        IMPLICIT ARGUMENTS :   NONE
    !        --------------------

    !     METHOD.
    !     -------

    !          TRANSMISSION FUNCTION ARE COMPUTED USING PADE APPROXIMANTS
    !     AND HORNER'S ALGORITHM.

    !     EXTERNALS.
    !     ----------

    !          NONE

    !     REFERENCE.
    !     ----------

    !        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
    !        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

    !     AUTHOR.
    !     -------
    !        JEAN-JACQUES MORCRETTE  *ECMWF*

    !     MODIFICATIONS.
    !     --------------
    !        ORIGINAL : 95-01-20

    !-----------------------------------------------------------------------

    !     DUMMY INTEGER SCALARS
    INTEGER :: KABS
    INTEGER :: KFDIA
    INTEGER :: KIDIA
    INTEGER :: KBDIM
    INTEGER :: KNU



    !-----------------------------------------------------------------------

    !*       0.1   ARGUMENTS
    !              ---------

    INTEGER :: KKIND(KABS)
    REAL(DP):: PU(KBDIM,KABS)
    REAL(DP):: PTR(KBDIM,KABS)

    !-----------------------------------------------------------------------

    !*       0.2   LOCAL ARRAYS
    !              ------------

    REAL(DP):: ZR1(KBDIM), ZR2(KBDIM), ZU(KBDIM)

    !     LOCAL INTEGER SCALARS
    INTEGER :: IA, JA, JL


    !-----------------------------------------------------------------------

    !*         1.      HORNER'S ALGORITHM TO COMPUTE TRANSMISSION FUNCTION


    DO JA = 1,KABS
      IA=KKIND(JA)
!IBM* novector
      DO JL = KIDIA,KFDIA
        ZU(JL) = PU(JL,JA)
        ZR1(JL) = APAD(KNU,IA,1) + ZU(JL) * (APAD(KNU,IA,2) + ZU(JL)&
             &* ( APAD(KNU,IA,3) + ZU(JL) * (APAD(KNU,IA,4) + ZU(JL)&
             &* ( APAD(KNU,IA,5) + ZU(JL) * (APAD(KNU,IA,6) + ZU(JL)&
             &* ( APAD(KNU,IA,7) ))))))

        ZR2(JL) = BPAD(KNU,IA,1) + ZU(JL) * (BPAD(KNU,IA,2) + ZU(JL)&
             &* ( BPAD(KNU,IA,3) + ZU(JL) * (BPAD(KNU,IA,4) + ZU(JL)&
             &* ( BPAD(KNU,IA,5) + ZU(JL) * (BPAD(KNU,IA,6) + ZU(JL)&
             &* ( BPAD(KNU,IA,7) ))))))


        !*         2.      ADD THE BACKGROUND TRANSMISSION


        PTR(JL,JA) = (ZR1(JL)/ZR2(JL)) * (1._DP-D(KNU,IA)) + D(KNU,IA)
      ENDDO
    ENDDO

    RETURN
  END SUBROUTINE SWTT1

  SUBROUTINE SWUVO3 &
       &( KIDIA,KFDIA,KLON,KNU,KABS &
       &, PU, PTR &
       & )

    !**** *SWUVO3* - COMPUTES THE SHORTWAVE TRANSMISSION FUNCTIONS

    !     PURPOSE.
    !     --------
    !           THIS ROUTINE COMPUTES THE TRANSMISSION FUNCTIONS FOR OZONE
    !     IN THE UV and VISIBLE SPECTRAL INTERVALS.

    !**   INTERFACE.
    !     ----------
    !          *SWUVO3* IS CALLED FROM *SW1S*.


    !        EXPLICIT ARGUMENTS :
    !        --------------------
    ! KNU    :                     ; INDEX OF THE SPECTRAL INTERVAL
    ! KABS   :                     ; NUMBER OF ABSORBERS
    ! PU     : (KLON,KABS)         ; ABSORBER AMOUNT
    !     ==== OUTPUTS ===
    ! PTR    : (KLON,KABS)         ; TRANSMISSION FUNCTION

    !        IMPLICIT ARGUMENTS :   NONE
    !        --------------------

    !     METHOD.
    !     -------

    !          TRANSMISSION FUNCTION ARE COMPUTED USING SUMS OF EXPONENTIALS

    !     EXTERNALS.
    !     ----------

    !          NONE

    !     REFERENCE.
    !     ----------
    !        SEE RADIATION'S PART OF THE MODEL'S DOCUMENTATION AND
    !        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

    !     AUTHOR.
    !     -------
    !        JEAN-JACQUES MORCRETTE  *ECMWF*
    !        Ph.DUBUISSON/B.BONNEL L.O.A.

    !     MODIFICATIONS.
    !     --------------
    !        ORIGINAL : 00-12-18

    !-----------------------------------------------------------------------


    !     DUMMY INTEGER SCALARS
    INTEGER :: KABS
    INTEGER :: KFDIA
    INTEGER :: KIDIA
    INTEGER :: KLON
    INTEGER :: KNU

    !-----------------------------------------------------------------------

    !*       0.1   ARGUMENTS
    !              ---------

    REAL(DP) :: PU(KLON,KABS)
    REAL(DP) :: PTR(KLON,KABS)

    !-----------------------------------------------------------------------

    !*       0.2   LOCAL ARRAYS
    !              ------------

    !     LOCAL INTEGER SCALARS
    INTEGER :: JA, JL, IEXP, JX


    IEXP=NEXPO3(KNU)
    !print *,'IEXP(',KNU,')=',IEXP
    !print *,(REXPO3(KNU,1,JX),JX=1,IEXP)
    !print *,(REXPO3(KNU,2,JX),JX=1,IEXP)

    DO JA = 1,KABS
      DO JL=KIDIA,KFDIA
        PTR(JL,JA)=0._DP
      END DO

      DO JX=1,IEXP
        DO JL = KIDIA,KFDIA
          PTR(JL,JA) = PTR(JL,JA) &
               &+REXPO3(KNU,1,JX)*EXP(-REXPO3(KNU,2,JX)*PU(JL,JA))
        END DO
      END DO
    ENDDO

    RETURN
  END SUBROUTINE SWUVO3

SUBROUTINE su_sw4

  !**** *SUSW*   - INITIALIZE COMMON MO_SW

  !     PURPOSE.
  !     --------
  !           INITIALIZE MO_SW, THE COMMON THAT CONTAINS COEFFICIENTS
  !           NEEDED TO RUN THE SHORTWAVE RADIATION SUBROUTINES

  !**   INTERFACE.
  !     ----------
  !        *CALL* *SUSW

  !        EXPLICIT ARGUMENTS :
  !        --------------------
  !        NONE

  !        IMPLICIT ARGUMENTS :
  !        --------------------
  !        COMMON MO_SW

  !     METHOD.
  !     -------
  !        SEE DOCUMENTATION

  !     EXTERNALS.
  !     ----------

  !     REFERENCE.
  !     ----------
  !        ECMWF RESEARCH DEPARTMENT DOCUMENTATION OF THE IFS

  !     AUTHOR.
  !     -------
  !        JEAN-JACQUES MORCRETTE *ECMWF*

  !     MODIFICATIONS.
  !     --------------
  !        ORIGINAL : 88-12-15
  !        97-04-16 JJ Morcrette  2 and 4 interval spectral resolution
  !        M.A. Giorgetta, MPI, June 2000:
  !        - simplified to 4 bands only,
  !        - setup code not relevant for ECHAM is uncommented by !!$
  !        C. Cagnazzo, June 2005:
  !        - from 4 to 6 bands, following 26Rcycle and 28Rcycle Morcrette Code

  !     ------------------------------------------------------------------


  !     ----------------------------------------------------------------
  REAL(DP):: ZAPAD6(6,3,7)  , ZBPAD6(6,3,7)  , ZD6(6,3)&
       & ,  ZRAY6(6,6)     , ZSUN6(6)       , ZSWCE6(6)  ,   ZSWCP6(6)


  !     LOCAL INTEGER SCALARS
  INTEGER :: JC3, JC6, JI, JJ, JW

  !     LOCAL REAL SCALARS
  REAL(DP):: ZPDU6IS, ZPRU6IS, ZPDH6IS, ZPRH6IS &
       & ,  ZTDH6IS, ZTDU6IS, ZTH6IS,  ZTU6IS  &
       & ,  ZPDUMG,  ZPRUMG,  ZPDH2O,  ZPRH2O  &
       & ,  ZH2O,    ZUMG

!!$  REAL(DP):: ZRTO1,   ZRTO2

  !     ----------------------------------------------------------------

  !*        1.  CLEAR-SKY ABSORPTION COEFFICIENTS FOR N SPECTRAL INTERVALS
  !             --------------------------------------------------------


  !*        1.2  COEFFICIENTS FOR FOUR SPECTRAL INTERVALS
  !              ----------------------------------------


  !* DERIVED FROM HITRAN APRIL 1992 with LOWTRAN P AND T SCALING

 ZTDH6IS = 0.450_DP
 ZTDU6IS = 0.375_DP
 ZTH6IS  = 273._DP
 ZTU6IS  = 273._DP
 ZPDH6IS = 0.90_DP
 ZPDU6IS = 0.75_DP
 ZPRH6IS = 101300._DP
 ZPRU6IS = 101300._DP

!* 1st spectral interval: U.V.  (0.185 - 0.25 Micron)

ZSUN6(1) = 0.001917_DP

ZD6(1,:)= (/ 1.000000000_DP, 1.000000000_DP, 0.000000000_DP /)

ZAPAD6(1, 1, :) = (/&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP/)
ZAPAD6(1, 2, :) = (/&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP/)
ZAPAD6(1, 3, :) = (/&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP/)


ZBPAD6(1, 1, :) = (/&
 &0.100000000E+01_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP/)
ZBPAD6(1, 2, :) = (/&
 &0.100000000E+01_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP/)
ZBPAD6(1, 3, :) = (/&
 &0.100000000E+01_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP/)

ZRAY6(1,:)= (/&
 &.43959002E+01_DP, .000000E+00_DP, .000000E+00_DP,&
 &.000000E+00_DP, .000000E+00_DP, .000000E+00_DP/)

ZSWCE6(1) = 0._DP
ZSWCP6(1) = 0._DP

 NEXPO3(1) = 7
 REXPO3(1, 1, :) =(/&
  & 0.051395E+00_DP, 0.048250E+00_DP, 0.112339E+00_DP,&
  & 0.101426E+00_DP, 0.007700E+00_DP, 0.441320E+00_DP,&
  & 0.237571E+00_DP /)
 REXPO3(1, 2, :) =(/&
  & 0.100022E+02_DP, 0.851159E+02_DP, 0.346737E+03_DP,&
  & 0.158501E+02_DP, 0.724223E+01_DP, 0.177828E+03_DP,&
  & 0.467708E+02_DP /)  

!* 2nd spectral interval: U.V.  (0.25 - 0.44 Micron)

ZSUN6(2) = 0.135708_DP

ZD6(2,:)= (/ 1.000000000_DP, 1.000000000_DP, 0.000000000_DP /)

ZAPAD6(2, 1, :) = (/&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP/)
ZAPAD6(2, 2, :) = (/&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP/)
ZAPAD6(2, 3, :) = (/&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP/)
 
ZBPAD6(2, 1, :) = (/&
 &0.100000000E+01_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP/)
ZBPAD6(2, 2, :) = (/&
 &0.100000000E+01_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP/)
ZBPAD6(2, 3, :) = (/&
 &0.100000000E+01_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP,&
 &0.000000000E-00_DP/)
 
ZRAY6(2,:)= (/&
 &.55503070E+00_DP, .000000E+00_DP, .000000E+00_DP,&
 &.000000E+00_DP, .000000E+00_DP, .000000E+00_DP/)

ZSWCE6(2) = 0._DP
ZSWCP6(2) = 0._DP

NEXPO3(2) = 7
REXPO3(2, 1, :) =(/&
  &0.043801E+00_DP, 0.078893E+00_DP, 0.036839E+00_DP,&
  &0.022503E+00_DP, 0.042333E+00_DP, 0.037870E+00_DP,&
  &0.737762E+00_DP /)
REXPO3(2, 2, :) =(/&
  &0.234249E+01_DP, 0.125170E+00_DP, 0.549527E+02_DP,&
  &0.257041E+03_DP, 0.476838E+00_DP, 0.911993E+01_DP,&
  &0.000000E+00_DP /)

!* 3rd spectral interval: Visible  (0.44 - 0.69 Micron)

ZSUN6(3) = 0.322135_DP

ZD6(3,:)= (/ 0.800000000_DP, 0.900000000_DP, 0.000000000_DP /)

ZAPAD6(3, 1, :) = (/&
 & 0.1762097E+03_DP,&
 & 0.1641762E+03_DP,&
 & 0.8687919E+02_DP,&
 & 0.0000000E-00_DP,&
 & 0.0000000E-00_DP,&
 & 0.0000000E-00_DP,&
 & 0.0000000E-00_DP/)
ZAPAD6(3, 2, :) = (/&
 & 0.5581224E+00_DP,&
 & 0.1748430E+03_DP,&
 & 0.1134123E+04_DP,&
 & 0.3490429E+03_DP,&
 & 0.0000000E-00_DP,&
 & 0.0000000E-00_DP,&
 & 0.0000000E-00_DP/)
ZAPAD6(3, 3, :) = (/&
 & 0.000000000E-00_DP,&
 & 0.000000000E-00_DP,&
 & 0.000000000E-00_DP,&
 & 0.000000000E-00_DP,&
 & 0.000000000E-00_DP,&
 & 0.000000000E-00_DP,&
 & 0.000000000E-00_DP/)

ZBPAD6(3, 1, :) = (/&
 & 0.1762097E+03_DP,&
 & 0.1663950E+03_DP,&
 & 0.8939724E+02_DP,&
 & 0.1000000E+01_DP,&
 & 0.0000000E-00_DP,&
 & 0.0000000E-00_DP,&
 & 0.0000000E-00_DP/)
ZBPAD6(3, 2, :) = (/&
 & 0.5581224E+00_DP,&
 & 0.1749251E+03_DP,&
 & 0.1159910E+04_DP,&
 & 0.3893268E+03_DP,&
 & 0.1000000E+01_DP,&
 & 0.0000000E-00_DP,&
 & 0.0000000E-00_DP/)
ZBPAD6(3, 3, :) = (/&
 & 0.100000000E+01_DP,&
 & 0.000000000E-00_DP,&
 & 0.000000000E-00_DP,&
 & 0.000000000E-00_DP,&
 & 0.000000000E-00_DP,&
 & 0.000000000E-00_DP,&
 & 0.000000000E-00_DP/)
 
ZRAY6(3,:)= (/&
 &.10528199E+00_DP, .000000E+00_DP, .000000E+00_DP,&
 &.000000E+00_DP, .000000E+00_DP, .000000E+00_DP/)

ZSWCE6(3) = 0._DP
ZSWCP6(3) = 0._DP

NEXPO3(3) = 6
REXPO3(3, 1, :) =(/&
  &0.063442E+00_DP, 0.058550E+00_DP, 0.237534E+00_DP,&
  &0.412292E+00_DP, 0.126141E+00_DP, 0.102041E+00_DP,&
  &0.000000E+00_DP /)
REXPO3(3, 2, :) =(/&
  &0.125170E+00_DP, 0.119209E-01_DP, 0.119209E+00_DP,&
  &0.417233E-01_DP, 0.894070E-01_DP, 0.000000E+00_DP,&
  &0.000000E+00_DP /)


!     ----------------------------------------------------------------

!* Near-Infrared (0.69 - 4.0 Microns) is sub-divided into:

!     ----------------------------------------------------------------

!* 0.69 - 1.19 Micron

!* UMG is O2 only

ZSUN6(4) = 0.326158_DP

ZD6(4,:)= (/ 0.000000000_DP, 0.900000000_DP, 1.000000000_DP /)

ZAPAD6(4, 1, :) = (/&
 & 0.1335726E+02_DP,&
 & 0.2939136E+04_DP,&
 & 0.4010585E+05_DP,&
 & 0.7195030E+05_DP,&
 & 0.1648338E+05_DP,&
 & 0.3373738E+03_DP,&
 & 0.0000000E+00_DP/)
ZAPAD6(4, 2, :) = (/&
 & 0.2001271E-01_DP,&
 & 0.2480831E+01_DP,&
 & 0.3444162E+02_DP,&
 & 0.4788946E+02_DP,&
 & 0.0000000E+00_DP,&
 & 0.0000000E+00_DP,&
 & 0.0000000E+00_DP/)
ZAPAD6(4, 3, :) = (/&
 & 0.000000000E+00_DP,&
 & 0.000000000E+00_DP,&
 & 0.000000000E+00_DP,&
 & 0.000000000E+00_DP,&
 & 0.000000000E+00_DP,&
 & 0.000000000E+00_DP,&
 & 0.000000000E+00_DP/)

ZBPAD6(4, 1, :) = (/&
 & 0.1335726E+02_DP,&
 & 0.2942327E+04_DP,&
 & 0.4077237E+05_DP,&
 & 0.7749017E+05_DP,&
 & 0.2123132E+05_DP,&
 & 0.6659687E+03_DP,&
 & 0.1000000E+01_DP/)
ZBPAD6(4, 2, :) = (/&
 & 0.2001271E-01_DP,&
 & 0.2549067E+01_DP,&
 & 0.3752433E+02_DP,&
 & 0.6276637E+02_DP,&
 & 0.1000000E+01_DP,&
 & 0.0000000E+00_DP,&
 & 0.0000000E+00_DP/)
ZBPAD6(4, 3, :) = (/&
 & 1.000000000E+00_DP,&
 & 0.000000000E+00_DP,&
 & 0.000000000E+00_DP,&
 & 0.000000000E+00_DP,&
 & 0.000000000E+00_DP,&
 & 0.000000000E+00_DP,&
 & 0.000000000E+00_DP/)
 
ZRAY6(4,:) = (/&
 &.16436996E-01_DP, .000000E+00_DP, .000000E+00_DP,&
 &.000000E+00_DP, .000000E+00_DP, .000000E+00_DP /)

ZSWCE6(4) = 0._DP
ZSWCP6(4) = 0._DP

NEXPO3(4) = 4
REXPO3(4, 1, :) =(/&
  &0.000074E+00_DP, 0.320194E+00_DP, 0.082915E+00_DP,&
  &0.596816E+00_DP, 0.000000E+00_DP, 0.000000E+00_DP,&
  &0.000000E+00_DP /)
REXPO3(4, 2, :) =(/&
  &0.232458E+00_DP, 0.119209E-01_DP, 0.178814E-01_DP,&
  &0.000000E+00_DP, 0.000000E+00_DP, 0.000000E+00_DP,&
  &0.000000E+00_DP /)


!     ----------------------------------------------------------------

!* 1.19 - 2.38 Microns

!* UMG is CO2 only

ZSUN6(5) = 0.180608_DP

ZD6(5,:)= (/ 0.000000000_DP, 0.800000000_DP, 1.000000000_DP /)

ZAPAD6(5, 1, :) = (/&
 & 0.3325841E-03_DP,&
 & 0.6194496E+00_DP,&
 & 0.1497138E+03_DP,&
 & 0.2314864E+04_DP,&
 & 0.2380109E+04_DP,&
 & 0.9553823E+02_DP,&
 & 0.0000000E+00_DP/)
ZAPAD6(5, 2, :) = (/&
 & 0.4552471E-03_DP,&
 & 0.4084154E+00_DP,&
 & 0.6114905E+01_DP,&
 & 0.7102540E+01_DP,&
 & 0.0000000E+00_DP,&
 & 0.0000000E+00_DP,&
 & 0.0000000E+00_DP/)
ZAPAD6(5, 3, :) = (/&
 & 0.000000000E+00_DP,&
 & 0.000000000E+00_DP,&
 & 0.000000000E+00_DP,&
 & 0.000000000E+00_DP,&
 & 0.000000000E+00_DP,&
 & 0.000000000E+00_DP,&
 & 0.000000000E+00_DP/)

ZBPAD6(5, 1, :) = (/&
 & 0.3325841E-03_DP,&
 & 0.6231947E+00_DP,&
 & 0.1553098E+03_DP,&
 & 0.2822458E+04_DP,&
 & 0.3885194E+04_DP,&
 & 0.2700235E+03_DP,&
 & 0.1000000E+01_DP/)
ZBPAD6(5, 2, :) = (/&
 & 0.4552471E-03_DP,&
 & 0.4088242E+00_DP,&
 & 0.6411905E+01_DP,&
 & 0.9444439E+01_DP,&
 & 0.1000000E+01_DP,&
 & 0.0000000E+00_DP,&
 & 0.0000000E+00_DP/)
ZBPAD6(5, 3, :) = (/&
 & 1.000000000E+00_DP,&
 & 0.000000000E+00_DP,&
 & 0.000000000E+00_DP,&
 & 0.000000000E+00_DP,&
 & 0.000000000E+00_DP,&
 & 0.000000000E+00_DP,&
 & 0.000000000E+00_DP/)

ZRAY6(5,:)= (/&
 &.18073079E-02_DP, .000000E+00_DP, .000000E+00_DP,&
 &.000000E+00_DP, .000000E+00_DP, .000000E+00_DP/)

ZSWCE6(5) = 0._DP
ZSWCP6(5) = 0._DP

NEXPO3(5) = 0
REXPO3(5, 1, :) =(/&
  &0.000000E+00_DP, 0.000000E+00_DP, 0.000000E+00_DP,&
  &0.000000E+00_DP, 0.000000E+00_DP, 0.000000E+00_DP,&
  &0.000000E+00_DP /)
REXPO3(5, 2, :) =(/&
  &0.000000E+00_DP, 0.000000E+00_DP, 0.000000E+00_DP,&
  &0.000000E+00_DP, 0.000000E+00_DP, 0.000000E+00_DP,&
  &0.000000E+00_DP /)

!     ----------------------------------------------------------------

!* 2.38 - 4.00 Microns

ZSUN6(6) = 0.033474_DP

ZD6(6,:)= (/ 0.000000000_DP, 0.000000000_DP, 0.000000000_DP /)

ZAPAD6(6, 1, :) = (/&
 &0.2122889E-06_DP,&
 &0.9030576E-03_DP,&
 &0.2431282E+00_DP,&
 &0.4901345E+01_DP,&
 &0.3996347E+01_DP,&
 &0.3910227E+01_DP,&
 &0.0000000E+00_DP/)
ZAPAD6(6, 2, :) = (/&
 & 0.1215163E-03_DP,&
 & 0.1222574E+00_DP,&
 & 0.9382420E+01_DP,&
 & 0.6875727E+02_DP,&
 & 0.2746421E+02_DP,&
 & 0.0000000E+00_DP,&
 & 0.0000000E+00_DP/)
ZAPAD6(6, 3, :) = (/&
 &0.263068898E+02_DP,&
 &0.146425875E+03_DP,&
 &0.860137809E+02_DP,&
 &0.000000000E+00_DP,&
 &0.000000000E+00_DP,&
 &0.000000000E+00_DP,&
 &0.000000000E+00_DP/)

ZBPAD6(6, 1, :) = (/&
 &0.2122889E-06_DP,&
 &0.9379083E-03_DP,&
 &0.2957335E+00_DP,&
 &0.8747190E+01_DP,&
 &0.1015794E+02_DP,&
 &0.1361277E+02_DP,&
 &0.1000000E+01_DP/)
ZBPAD6(6, 2, :) = (/&
 & 0.1215163E-03_DP,&
 & 0.1255648E+00_DP,&
 & 0.1060119E+02_DP,&
 & 0.8414439E+02_DP,&
 & 0.4299438E+02_DP,&
 & 0.1000000E+01_DP,&
 & 0.0000000E+00_DP/)
ZBPAD6(6, 3, :) = (/&
 &0.263068898E+02_DP,&
 &0.152569217E+03_DP,&
 &0.976791971E+02_DP,&
 &0.100000000E+01_DP,&
 &0.000000000E+00_DP,&
 &0.000000000E+00_DP,&
 &0.000000000E+00_DP/)

ZRAY6(6,:)= (/&
 &.13618247E-03_DP, .000000E+00_DP, .000000E+00_DP,&
 &.000000E+00_DP, .000000E+00_DP, .000000E+00_DP/)

ZSWCE6(6) = 0._DP
ZSWCP6(6) = 0._DP

NEXPO3(6) = 0
REXPO3(6, 1, :) =(/&
  &0.000000E+00_DP, 0.000000E+00_DP, 0.000000E+00_DP,&
  &0.000000E+00_DP, 0.000000E+00_DP, 0.000000E+00_DP,&
  &0.000000E+00_DP /)
REXPO3(6, 2, :) =(/&
  &0.000000E+00_DP, 0.000000E+00_DP, 0.000000E+00_DP,&
  &0.000000E+00_DP, 0.000000E+00_DP, 0.000000E+00_DP,&
  &0.000000E+00_DP /)


!!$  !=====================================================================
!!$
!!$  !*    2.2   SPECTRAL ALBEDO OF SNOW
!!$  !           Fresh and Melting Snow (Warren, 1982)
!!$
!!$  ZSNAL4(1:4)= (/ 0.920_DP, 0.798_DP , 0.159_DP, 0.010_DP /)
!!$  ZSNML4(1:4)= (/ 0.860_DP, 0.664_DP , 0.092_DP, 0.010_DP /)
!!$
!!$  !     ----------------------------------------------------------------
!!$
!!$  !*    2.3   WEIGHTS FOR SPECTRAL ALBEDO OF SOIL AND VEGETATION
!!$  !           a la Briegleb and Ramanathan (1982)
!!$
!!$  ZWEIS4(1:4)= (/ 1._DP, 2._DP, 2._DP, 1._DP /)
!!$  ZWEIV4(1:4)= (/ 1._DP, 4._DP, 2._DP, 1._DP /)
!!$
!!$  !     ----------------------------------------------------------------
!!$
!!$  !*    2.4   OPTICAL PARAMETERS FOR RAIN DROPS
!!$  !           Savijarvi et al. (1996)
!!$
!!$  ZRTO1 =  0.003_DP
!!$  ZRTO2 = -0.22_DP
!!$  ! CAUTION JUST TEMPORARY PARAMETERS      
!!$
!!$  ZROMA4(1:4)= (/ 0.00008_DP , 0.0105_DP , 0.264_DP  , 0.465_DP   /)
!!$  ZROMB4(1:4)= (/ 0.23_DP    , 0.22_DP   , 0.09_DP   , 0.001_DP   /)
!!$  ZRASY4(1:4)= (/ 0.88_DP    , 0.89_DP   , 0.94_DP   , 0.97_DP    /)
!!$
!!$  ZRA4(1:4)= (/ 1.5_DP     , 1.5_DP    , 1.5_DP    , 1.5_DP     /)
!!$  ZRB4(1:4)= (/ 0.50_DP    , 0.78_DP   , 1.13_DP   , 2.00_DP    /)
!!$  ZRC4(1:4)= (/ 5.58E-7_DP , 2.18E-5_DP, 8.55E-4_DP, 1.94E-1_DP /)
!!$  ZRD4(1:4)= (/ 1.25E-7_DP , 2.25E-5_DP, 1.28E-3_DP, 8.04E-3_DP /)
!!$  ZRE4(1:4)= (/ 0.841_DP   , 0.821_DP  , 0.786_DP  , 0.820_DP   /)
!!$  ZRF4(1:4)= (/ 2.08E-3_DP , 3.06E-3_DP, 5.32E-3_DP, 5.59E-3_DP /)
!!$
!!$
!!$  !=====================================================================

  !*       2.    SET VALUES.
  !              -----------


  ZPDH2O = ZPDH6IS
  ZPDUMG = ZPDU6IS
  ZPRH2O = ZPRH6IS
  ZPRUMG = ZPRU6IS
  RTDH2O = ZTDH6IS
  RTDUMG = ZTDU6IS
  RTH2O  = ZTH6IS
  RTUMG  = ZTU6IS

  RPDH1=ZPDH2O+1._DP
  RPDU1=ZPDUMG+1._DP
  ZH2O=1._DP/( 10._DP* RG * RPDH1 )
  ZUMG=1._DP/( 10._DP* RG * RPDU1 )
  RPNU = ZUMG/(ZPRUMG**ZPDUMG)
  RPNH = ZH2O/(ZPRH2O**ZPDH2O)

!!$  RHSRTA=ZRTO1
!!$  RHSRTB=ZRTO2
  DO JW=1,6
     RSUN (JW)=ZSUN6(JW)

     RSWCE(JW)=ZSWCE6(JW)
     RSWCP(JW)=ZSWCP6(JW)

!!$     RSNOALB(JW)=ZSNAL4(JW)
!!$     RSNOMEL(JW)=ZSNML4(JW)
!!$
!!$     RWEIGS(JW)=ZWEIS4(JW)
!!$     RWEIGV(JW)=ZWEIV4(JW)
!!$
!!$     RROMA(JW)=ZROMA4(JW)
!!$     RROMB(JW)=ZROMB4(JW)
!!$     RRASY(JW)=ZRASY4(JW)
!!$     RHSRA(JW)=ZRA4(JW)
!!$     RHSRB(JW)=ZRB4(JW)
!!$     RHSRC(JW)=ZRC4(JW)
!!$     RHSRD(JW)=ZRD4(JW)
!!$     RHSRE(JW)=ZRE4(JW)
!!$     RHSRF(JW)=ZRF4(JW)

     DO JC3=1,3
        D(JW,JC3)=ZD6(JW,JC3)
     ENDDO
     DO JC6=1,6
        RRAY(JW,JC6)=ZRAY6(JW,JC6)
     ENDDO
     DO JI=1,3
        DO JJ=1,7
           APAD(JW,JI,JJ)=ZAPAD6(JW,JI,JJ)
           BPAD(JW,JI,JJ)=ZBPAD6(JW,JI,JJ)
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE su_sw4

END MODULE MO_ECHAM5_SW
