!>
!! @brief Module to provide aerosol properties from climatology
!!
!! @par Description
!!   This module contains a global aerosol climatology following Tanre and
!!   the necessary routines to make use of it.  The climatology distinguishes
!!   spacial distributions of sea, land, urban, and desert aerosols. Further,
!!   the climatology contains constant background aerosols of tropospheric, 
!!   stratospheric.  There is no annual cycle in the climatology, and volcanic
!!   aerosols (formerly aerosol type 4) have been excluded.
!!
!! @par
!!   These aerosol consist of four types
!!   <ol>
!!    <li> land + desert + tropospheric background
!!    <li> sea
!!    <li> urban
!!    <li> stratospheric background
!!   </ol>
!!   whose optical properties are combined
!!
!! @author Bjorn Stevens, MPI-M, Hamburg (2009-09-19)
!!
!! $ID: n/a$
!!
!! @par Origin
!!   The original source was written by J. F. Geleyn, ECMWF, (1982-11) and
!!   modified over the years, including rewrites by L.Kornblueh (1998-05); 
!!   M.A. Giorgetta, MPI-M (1999-04);  U. Schulzweida MPI, (1998-05,2002-05).
!!   This version was rewritten for the ICON standard and to isolate aerosol
!!   functionality by the author.
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
MODULE mo_aero_tanre

  USE mo_kind,          ONLY: wp
!LT
!  USE mo_exception,     ONLY: finish
!  USE mo_geoloc,        ONLY: zaes_x, zael_x, zaeu_x, zaed_x
!  USE mo_gaussgrid,     ONLY: coslon, sinlon, gl_twomu
!  USE mo_decomposition, ONLY: ldc => local_decomposition
!  USE mo_transpose,     ONLY: reorder

  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: aerosol_optics_tanre, cleanup_aero_tanre

  LOGICAL :: initialized = .FALSE.
  !
  ! 1.0 Tables defining aerosol properties given as a horizontal distribution of
  !  sea, land, urban and desert aerosols given by a T10 spectral representatin
  !--------------------------------
  !
  INTEGER, PARAMETER :: naer=4 !< number of aerosol mixture types

  REAL(wp), PARAMETER, DIMENSION(66) :: caesc= & !< sea aero cosine coeffs
       (/ 0.6688e+00_wp, -0.1172e+00_wp, -0.1013e+00_wp,  0.1636e-01_wp, &
       & -0.3699e-01_wp,  0.1775e-01_wp, -0.9635e-02_wp,  0.1290e-02_wp, &
       &  0.4681e-04_wp, -0.9106e-04_wp,  0.9355e-04_wp, -0.7076e-01_wp, &
       & -0.1782e-01_wp,  0.1856e-01_wp,  0.1372e-01_wp,  0.8210e-04_wp, &
       &  0.2149e-02_wp,  0.4856e-03_wp,  0.2231e-03_wp,  0.1824e-03_wp, &
       &  0.1960e-05_wp,  0.2057e-01_wp,  0.2703e-01_wp,  0.2424e-01_wp, &
       &  0.9716e-02_wp,  0.1312e-02_wp, -0.8846e-03_wp, -0.3347e-03_wp, &
       &  0.6231e-04_wp,  0.6397e-04_wp, -0.3341e-02_wp, -0.1295e-01_wp, &
       & -0.4598e-02_wp,  0.3242e-03_wp,  0.8122e-03_wp, -0.2975e-03_wp, &
       & -0.7757e-04_wp,  0.7793e-04_wp,  0.4455e-02_wp, -0.1584e-01_wp, &
       & -0.2551e-02_wp,  0.1174e-02_wp,  0.1335e-04_wp,  0.5112e-04_wp, &
       &  0.5605e-04_wp,  0.7412e-04_wp,  0.1857e-02_wp, -0.1917e-03_wp, &
       &  0.4460e-03_wp,  0.1767e-04_wp, -0.5281e-04_wp, -0.5043e-03_wp, &
       &  0.2467e-03_wp, -0.2497e-03_wp, -0.2377e-04_wp, -0.3954e-04_wp, &
       &  0.2666e-03_wp, -0.8186e-03_wp, -0.1441e-03_wp, -0.1904e-04_wp, &
       &  0.3337e-03_wp, -0.1696e-03_wp, -0.2503e-04_wp,  0.1239e-03_wp, &
       & -0.9983e-04_wp, -0.5283e-04_wp/)

  REAL(wp), PARAMETER, DIMENSION(55) :: caess= & !< sea aero sine coeffs
       (/-0.3374e-01_wp, -0.3247e-01_wp, -0.1012e-01_wp,  0.6002e-02_wp, &
       &  0.5190e-02_wp,  0.7784e-03_wp, -0.1090e-02_wp,  0.3294e-03_wp, &
       &  0.1719e-03_wp, -0.5866e-05_wp, -0.4124e-03_wp, -0.3742e-01_wp, &
       & -0.5054e-02_wp,  0.3430e-02_wp,  0.5513e-03_wp, -0.6235e-03_wp, &
       &  0.2892e-03_wp, -0.9730e-04_wp,  0.7078e-04_wp, -0.3300e-01_wp, &
       &  0.5104e-03_wp, -0.2156e-02_wp, -0.3194e-02_wp, -0.5079e-03_wp, &
       & -0.5517e-03_wp,  0.4632e-04_wp,  0.5369e-04_wp, -0.2731e-01_wp, &
       &  0.5126e-02_wp,  0.2241e-02_wp, -0.5789e-03_wp, -0.3048e-03_wp, &
       & -0.1774e-03_wp,  0.1946e-05_wp, -0.8247e-02_wp,  0.2338e-02_wp, &
       &  0.1021e-02_wp,  0.1575e-04_wp,  0.2612e-05_wp,  0.1995e-04_wp, &
       & -0.1319e-02_wp,  0.1384e-02_wp, -0.4159e-03_wp, -0.2337e-03_wp, &
       &  0.5764e-04_wp,  0.1495e-02_wp, -0.3727e-03_wp,  0.6075e-04_wp, &
       & -0.4642e-04_wp,  0.5368e-03_wp, -0.7619e-04_wp,  0.3774e-04_wp, &
       &  0.1206e-03_wp, -0.4104e-06_wp,  0.2158e-04_wp/)

  REAL(wp), PARAMETER, DIMENSION(66) :: caelc= & !< land aero cosine coeffs
       (/ 0.1542e+00_wp,  0.8245e-01_wp, -0.1879e-03_wp,  0.4864e-02_wp, &
       & -0.5527e-02_wp, -0.7966e-02_wp, -0.2683e-02_wp, -0.2011e-02_wp, &
       & -0.8889e-03_wp, -0.1058e-03_wp, -0.1614e-04_wp,  0.4206e-01_wp, &
       &  0.1912e-01_wp, -0.9476e-02_wp, -0.6780e-02_wp,  0.1767e-03_wp, &
       & -0.5422e-03_wp, -0.7753e-03_wp, -0.2106e-03_wp, -0.9870e-04_wp, &
       & -0.1721e-04_wp, -0.9536e-02_wp, -0.9580e-02_wp, -0.1050e-01_wp, &
       & -0.5747e-02_wp, -0.1282e-02_wp,  0.2248e-03_wp,  0.1694e-03_wp, &
       & -0.4782e-04_wp, -0.2441e-04_wp,  0.5781e-03_wp,  0.6212e-02_wp, &
       &  0.1921e-02_wp, -0.1102e-02_wp, -0.8145e-03_wp,  0.2497e-03_wp, &
       &  0.1539e-03_wp, -0.2538e-04_wp, -0.3993e-02_wp,  0.9777e-02_wp, &
       &  0.4837e-03_wp, -0.1304e-02_wp,  0.2417e-04_wp, -0.1370e-04_wp, &
       & -0.3731e-05_wp,  0.1922e-02_wp, -0.5167e-03_wp,  0.4295e-03_wp, &
       & -0.1888e-03_wp,  0.2427e-04_wp,  0.4012e-04_wp,  0.1529e-02_wp, &
       & -0.2120e-03_wp,  0.8166e-04_wp,  0.2579e-04_wp,  0.3488e-04_wp, &
       &  0.2140e-03_wp,  0.2274e-03_wp, -0.3447e-05_wp, -0.1075e-04_wp, &
       & -0.1018e-03_wp,  0.2864e-04_wp,  0.3442e-04_wp, -0.1002e-03_wp, &
       &  0.7117e-04_wp,  0.2045e-04_wp/)

  REAL(wp), PARAMETER, DIMENSION(55) :: caels= & !< land aero sine coeffs
       (/ 0.1637e-01_wp,  0.1935e-01_wp,  0.1080e-01_wp,  0.2784e-02_wp, &
       &  0.1606e-03_wp,  0.1860e-02_wp,  0.1263e-02_wp, -0.2707e-03_wp, &
       & -0.2290e-03_wp, -0.9761e-05_wp, -0.7317e-02_wp,  0.2465e-01_wp, &
       &  0.6799e-02_wp, -0.1913e-02_wp,  0.1382e-02_wp,  0.6691e-03_wp, &
       &  0.1414e-03_wp,  0.3527e-04_wp, -0.5210e-04_wp,  0.1873e-01_wp, &
       &  0.2977e-02_wp,  0.4650e-02_wp,  0.2509e-02_wp,  0.3680e-03_wp, &
       &  0.1481e-03_wp, -0.6594e-04_wp, -0.5634e-04_wp,  0.1592e-01_wp, &
       & -0.1875e-02_wp, -0.1093e-02_wp,  0.3022e-03_wp,  0.2625e-03_wp, &
       &  0.3252e-04_wp, -0.3803e-04_wp,  0.4218e-02_wp, -0.1843e-02_wp, &
       & -0.1351e-02_wp, -0.2952e-03_wp, -0.8171e-05_wp, -0.1473e-04_wp, &
       &  0.9076e-03_wp, -0.1057e-02_wp,  0.2676e-03_wp,  0.1307e-03_wp, &
       & -0.3628e-04_wp, -0.9158e-03_wp,  0.4335e-03_wp,  0.2927e-04_wp, &
       &  0.6602e-04_wp, -0.3570e-03_wp,  0.5760e-04_wp, -0.3465e-04_wp, &
       & -0.8535e-04_wp, -0.2011e-04_wp,  0.6612e-06_wp/)

  REAL(wp), PARAMETER, DIMENSION(66) :: caeuc= & !< urban aero cosine coeffs
       (/ 0.8005e-01_wp,  0.7095e-01_wp,  0.2014e-01_wp, -0.1412e-01_wp, &
       & -0.2425e-01_wp, -0.1332e-01_wp, -0.2904e-02_wp,  0.5068e-03_wp, &
       &  0.9369e-03_wp,  0.4114e-03_wp,  0.7549e-04_wp,  0.1922e-01_wp, &
       &  0.2534e-01_wp,  0.2088e-01_wp,  0.1064e-01_wp,  0.1063e-02_wp, &
       & -0.2526e-02_wp, -0.2091e-02_wp, -0.9660e-03_wp, -0.2030e-03_wp, &
       &  0.3865e-04_wp, -0.9900e-02_wp, -0.5964e-02_wp,  0.2223e-02_wp, &
       &  0.4941e-02_wp,  0.3277e-02_wp,  0.1038e-02_wp, -0.1480e-03_wp, &
       & -0.2844e-03_wp, -0.1208e-03_wp,  0.3999e-02_wp,  0.6282e-02_wp, &
       &  0.2813e-02_wp,  0.1475e-02_wp,  0.4571e-03_wp, -0.1349e-03_wp, &
       & -0.9011e-04_wp, -0.1936e-04_wp,  0.1994e-02_wp,  0.3540e-02_wp, &
       &  0.8837e-03_wp,  0.1992e-03_wp,  0.3092e-04_wp, -0.7979e-04_wp, &
       & -0.2664e-04_wp, -0.5006e-04_wp,  0.6447e-03_wp,  0.5550e-03_wp, &
       &  0.1197e-03_wp,  0.6657e-04_wp,  0.1488e-04_wp, -0.9141e-04_wp, &
       & -0.2896e-03_wp, -0.1561e-03_wp, -0.6524e-04_wp, -0.1559e-04_wp, &
       & -0.1082e-03_wp, -0.4126e-03_wp, -0.1732e-03_wp, -0.8286e-04_wp, &
       & -0.1993e-04_wp,  0.3850e-04_wp,  0.2870e-04_wp,  0.4493e-04_wp, &
       &  0.4721e-04_wp,  0.1338e-04_wp/)

  REAL(wp), PARAMETER, DIMENSION(55) :: caeus= & !< urban aero sine coeffs
       (/ 0.6646e-02_wp,  0.8373e-02_wp,  0.5463e-02_wp,  0.4554e-02_wp, &
       &  0.3301e-02_wp,  0.5725e-03_wp, -0.7482e-03_wp, -0.6222e-03_wp, &
       & -0.2603e-03_wp, -0.5127e-04_wp, -0.3849e-04_wp,  0.9741e-02_wp, &
       &  0.8190e-02_wp,  0.5712e-02_wp,  0.3039e-02_wp,  0.5290e-03_wp, &
       & -0.2044e-03_wp, -0.2309e-03_wp, -0.1160e-03_wp,  0.9160e-02_wp, &
       &  0.1286e-01_wp,  0.1170e-01_wp,  0.5491e-02_wp,  0.1393e-02_wp, &
       & -0.6288e-04_wp, -0.2715e-03_wp, -0.1047e-03_wp,  0.4873e-02_wp, &
       &  0.3545e-02_wp,  0.3069e-02_wp,  0.1819e-02_wp,  0.6947e-03_wp, &
       &  0.1416e-03_wp, -0.1538e-04_wp, -0.4351e-03_wp, -0.1907e-02_wp, &
       & -0.5774e-03_wp, -0.2247e-03_wp,  0.5345e-04_wp,  0.9052e-04_wp, &
       & -0.3972e-04_wp, -0.9665e-04_wp,  0.7912e-04_wp, -0.1094e-04_wp, &
       & -0.6776e-05_wp,  0.2724e-03_wp,  0.1973e-03_wp,  0.6837e-04_wp, &
       &  0.4313e-04_wp, -0.7174e-05_wp,  0.8527e-05_wp, -0.2160e-05_wp, & 
       & -0.7852e-04_wp,  0.3453e-06_wp, -0.2402e-05_wp/)

  REAL(wp), PARAMETER, DIMENSION(66) :: caedc= & !< desert aero cosine coeffs
       (/ 0.2840e-01_wp,  0.1775e-01_wp, -0.1069e-01_wp, -0.1553e-01_wp, &
       & -0.3299e-02_wp,  0.3583e-02_wp,  0.2274e-02_wp,  0.5767e-04_wp, &
       & -0.3678e-03_wp, -0.1050e-03_wp,  0.2133e-04_wp,  0.2326e-01_wp, &
       &  0.1566e-01_wp, -0.3130e-02_wp, -0.8253e-02_wp, -0.2615e-02_wp, &
       &  0.1247e-02_wp,  0.1059e-02_wp,  0.1196e-03_wp, -0.1303e-03_wp, &
       & -0.5094e-04_wp,  0.1185e-01_wp,  0.7238e-02_wp, -0.1562e-02_wp, &
       & -0.3665e-02_wp, -0.1182e-02_wp,  0.4678e-03_wp,  0.4448e-03_wp, &
       &  0.8307e-04_wp, -0.3468e-04_wp,  0.5273e-02_wp,  0.3037e-02_wp, &
       & -0.4014e-03_wp, -0.1202e-02_wp, -0.4647e-03_wp,  0.5148e-04_wp, &
       &  0.1014e-03_wp,  0.2996e-04_wp,  0.2505e-02_wp,  0.1495e-02_wp, &
       &  0.2438e-03_wp, -0.1223e-03_wp, -0.7669e-04_wp, -0.1638e-04_wp, &
       &  0.1869e-05_wp,  0.1094e-02_wp,  0.6131e-03_wp,  0.1508e-03_wp, &
       &  0.1765e-04_wp,  0.1360e-05_wp, -0.7998e-06_wp,  0.4475e-03_wp, &
       &  0.2737e-03_wp,  0.6430e-04_wp, -0.6759e-05_wp, -0.6761e-05_wp, &
       &  0.1992e-03_wp,  0.1531e-03_wp,  0.4828e-04_wp,  0.5103e-06_wp, &
       &  0.7454e-04_wp,  0.5917e-04_wp,  0.2152e-04_wp,  0.9300e-05_wp, &
       &  0.9790e-05_wp, -0.8853e-05_wp/)

  REAL(wp), PARAMETER, DIMENSION(55) :: caeds= &  !< desert aero sine coeffs
       (/ 0.9815e-02_wp,  0.8436e-02_wp,  0.1087e-02_wp, -0.2717e-02_wp, &
       & -0.1755e-02_wp, -0.1559e-03_wp,  0.2367e-03_wp,  0.8808e-04_wp, &
       &  0.2001e-05_wp, -0.1244e-05_wp,  0.1041e-01_wp,  0.8039e-02_wp, &
       &  0.1005e-02_wp, -0.1981e-02_wp, -0.1090e-02_wp,  0.1595e-05_wp, &
       &  0.1787e-03_wp,  0.4644e-04_wp, -0.1052e-04_wp,  0.6593e-02_wp, &
       &  0.3983e-02_wp, -0.1527e-03_wp, -0.1235e-02_wp, -0.5078e-03_wp, &
       &  0.3649e-04_wp,  0.1005e-03_wp,  0.3182e-04_wp,  0.3225e-02_wp, &
       &  0.1672e-02_wp, -0.7752e-04_wp, -0.4312e-03_wp, -0.1872e-03_wp, &
       & -0.1666e-04_wp,  0.1872e-04_wp,  0.1133e-02_wp,  0.5643e-03_wp, &
       &  0.7747e-04_wp, -0.2980e-04_wp, -0.2092e-04_wp, -0.8590e-05_wp, &
       &  0.2988e-03_wp,  0.6714e-04_wp, -0.6249e-05_wp,  0.1052e-04_wp, &
       &  0.8790e-05_wp,  0.1569e-03_wp, -0.1175e-04_wp, -0.3033e-04_wp, &
       & -0.9777e-06_wp,  0.1101e-03_wp,  0.6827e-05_wp, -0.1023e-04_wp, &
       &  0.4231e-04_wp,  0.4905e-05_wp,  0.6229e-05_wp/)

  REAL(wp) :: taua(naer) = & !< otpical depth factors for four aerosol types
       (/ 0.730719_wp, 0.912819_wp, 0.725059_wp, 0.682188_wp /)
  REAL(wp) :: cga(naer)  = & !< asymmetry factor for four aerosol types
       (/ 0.647596_wp, 0.739002_wp, 0.580845_wp, 0.624246_wp /)
  REAL(wp) :: piza(naer) = & !< sngl sctr albedo for four aerosol types
       (/ 0.872212_wp, 0.982545_wp, 0.623143_wp, 0.997975_wp /)
  !
  !> longwave optical depth factor in five LW bands
  ! 
  REAL(wp), DIMENSION(5,naer) :: caer=RESHAPE((/&  
       & 0.038520_wp, 0.037196_wp, 0.040532_wp, 0.054934_wp, 0.038520_wp ,&
       & 0.126130_wp, 0.183130_wp, 0.103570_wp, 0.064106_wp, 0.126130_wp ,&
       & 0.012579_wp, 0.013649_wp, 0.018652_wp, 0.025181_wp, 0.012579_wp ,&
       & 0.013792_wp, 0.026810_wp, 0.052203_wp, 0.066338_wp, 0.013792_wp  &
       & /), SHAPE=(/5,naer/))
  !
  ! --- variables for vertical distribution of different aerosol types
  !
  REAL(wp), ALLOCATABLE :: cvdaes(:) !< normalized vert. dist of sea aerosols
  REAL(wp), ALLOCATABLE :: cvdael(:) !< normalized vert. dist of land aerosols
  REAL(wp), ALLOCATABLE :: cvdaeu(:) !< normalized vert. dist of urban aerosols
  REAL(wp), ALLOCATABLE :: cvdaed(:) !< normalized vert. dist of desert aerosols


  REAL(wp)  :: caeadk(3) ! for moisture adsorption of aerosols
  REAL(wp)  :: caeadm    ! for moisture adsorption of aerosols
  REAL(wp)  :: caeops    ! optical depths of sea aerosols
  REAL(wp)  :: caeopl    ! optical depths of land aerosols
  REAL(wp)  :: caeopu    ! optical depths of urban aerosols
  REAL(wp)  :: caeopd    ! optical depths of desert aerosols
  REAL(wp)  :: ctrpt     ! temperature exponent for strat. aerosol
  !
  ! --- Pressure normalized optical depths, here the surface pressure is set to
  !     101325 hPa and the tropopause is fixed at 19330 hPa
  !
  REAL(wp), PARAMETER ::   &
       ctrbga = 0.03_wp/(101325.0_wp-19330.0_wp), & !< troposphere background
       cstbga = 0.045_wp/19330.0_wp                 !< stratosphere background

CONTAINS
  !-----------------------------------------------------------------------------
  !>
  !! @brief provides aerosol optical properties for the TANRE climatology
  !
  SUBROUTINE aerosol_optics_tanre(krow, kproma, kbdim, klev, nb_lw, nb_sw, &
       etah, ppd_hl, pp_fl, tk_hl, tau_lw, tau_sw, cg_sw, piz_sw )

    INTEGER, PARAMETER :: jbands = 16        !< number of lw bands in rad trans
    INTEGER, PARAMETER :: map_band(jbands) & !< map between lw bands and model
         &                = (/1,5,2,2,2,3,4,3,1,1,1,1,1,1,1,1/)

    INTEGER, INTENT (IN)  :: kproma, kbdim, klev, krow, nb_sw, nb_lw
    REAL(wp), INTENT (IN)::          &
         etah(:),                    & !< normalized pressure height
         ppd_hl(:,:),                & !< pressure thickness at half level [hPa]
         pp_fl(:,:),                 & !< pressure at full level [hPa]
         tk_hl(:,:)                    !< temperature at half level [K]
    REAL(wp), INTENT (OUT)::         &
         tau_lw(kbdim,klev,nb_lw),   & !< aerosol optical depth
         tau_sw(kbdim,klev,nb_sw),   & !< aerosol optical depth
         cg_sw(kbdim,klev,nb_sw),    & !< aerosol asymmetry factor
         piz_sw(kbdim,klev,nb_sw)      !< aerosol single scattering abledo

    REAL(wp) :: &
         aer(kbdim,klev,naer),       & !< aerosol optical depth
         zaer(kbdim)                   !< scratch array for aerosol

    INTEGER :: jl, jk, jb, ja          !< do loop indicies
    !
    ! --- initialization
    !
!LT
!    IF (nb_lw /= jbands)  CALL finish('aerosol_tau','band mismatch')
!    IF (nb_sw /= 6 .and. nb_sw /= 14) CALL finish('aerosol_tau','band mismatch')
    IF (.not.initialized) CALL setup_aero_tanre(etah)

    aer(1:kproma,:,:) = aero_tanre(krow,kproma,ppd_hl,pp_fl,tk_hl)

    cg_sw(:,:,:)  = 0.0_wp
    piz_sw(:,:,:) = 0.0_wp
    tau_sw(:,:,:) = 0.0_wp
    tau_lw(:,:,:) = 0.0_wp
    !
    ! --- Optical properties of Tanre Aerosol.  Note that the aerosol has its
    ! optical properties defined in five lw bands, which are then mapped onto
    ! the 16 RRTM bands through 'map_band'.  Also note that the 'AER' aerosols
    ! are entered from top to bottom (hence top-down index)
    !
    DO ja = 1, naer
      DO jk = 1, klev
        zaer(1:kproma) = aer(1:kproma,klev+1-jk,ja)
        DO jb = 1, nb_lw
          tau_lw(1:kproma,jk,jb) = tau_lw(1:kproma,jk,jb) &
               + caer(map_band(jb),ja)*zaer(1:kproma)
        END DO
        DO jb = 1, nb_sw
          tau_sw(1:kproma,jk,jb) = tau_sw(1:kproma,jk,jb)+zaer(1:kproma) &
               & *taua(ja)
          piz_sw(1:kproma,jk,jb) = piz_sw(1:kproma,jk,jb)+zaer(1:kproma) &
               & *taua(ja)*piza(ja)
          cg_sw(1:kproma,jk,jb)  = cg_sw(1:kproma,jk,jb) +zaer(1:kproma) &
               & *taua(ja)*piza(ja)*cga(ja)
        ENDDO
      END DO
    END DO
    IF (nb_sw /= 6) THEN
      DO jk = 1, klev
        DO jb = 1, nb_sw
          DO jl = 1, kproma
            IF (piz_sw(jl,jk,jb) > EPSILON(1.0_wp)) THEN
              cg_sw(jl,jk,jb)  =  cg_sw(jl,jk,jb)/piz_sw(jl,jk,jb)
              piz_sw(jl,jk,jb) = piz_sw(jl,jk,jb)/tau_sw(jl,jk,jb)
            END IF
          END DO
        END DO
      END DO
    END IF

  END SUBROUTINE aerosol_optics_tanre

  PURE FUNCTION aero_tanre(krow,kproma,pdp,ppf,pth)

    INTEGER, INTENT (IN) :: krow, kproma
    REAL(wp), INTENT(IN) :: &
         pdp(:,:),          & ! pressure thickness hPa
         ppf(:,:),          & ! pressure on full levels
         pth(:,:)             ! temperature

    REAL(wp), DIMENSION(kproma,SIZE(pdp,2),naer) :: aero_tanre, zaero
    REAL(wp), DIMENSION(kproma)                  :: zaetrn, zaetro
    REAL(wp)  :: zaeqsn, zaeqln, zaequn, zaeqdn, zaeqso, zaeqlo, zaequo, &
         &       zaeqdo, zaetr

    INTEGER :: jl,jk,klev

    klev = SIZE(pdp,2)
    !
    ! --- output: aerosol longitude height sections in aero_tanre
    !
    DO jl = 1, kproma
      zaetro(jl) = 1._wp
    END DO

    DO jk = 1, klev
      DO jl = 1, kproma
        zaeqso = caeops*zaes_x(jl,krow)*cvdaes(jk)
        zaeqsn = caeops*zaes_x(jl,krow)*cvdaes(jk+1)
        zaeqlo = caeopl*zael_x(jl,krow)*cvdael(jk)
        zaeqln = caeopl*zael_x(jl,krow)*cvdael(jk+1)
        zaequo = caeopu*zaeu_x(jl,krow)*cvdaeu(jk)
        zaequn = caeopu*zaeu_x(jl,krow)*cvdaeu(jk+1)
        zaeqdo = caeopd*zaed_x(jl,krow)*cvdaed(jk)
        zaeqdn = caeopd*zaed_x(jl,krow)*cvdaed(jk+1)

        IF (ppf(jl,jk) < 999._wp) THEN
          zaetr= 1._wp ! above 10 hPa 
        ELSE
          zaetrn(jl) = zaetro(jl)*(MIN(1.0_wp,pth(jl,jk)/pth(jl,jk+1)))**ctrpt
          zaetr      = SQRT(zaetrn(jl)*zaetro(jl))
          zaetro(jl) = zaetrn(jl)
        END IF

        zaero(jl,jk,1) = (1._wp-zaetr)*(ctrbga*pdp(jl,jk)    &
             &           + zaeqln - zaeqlo + zaeqdn - zaeqdo)
        zaero(jl,jk,2) = (1._wp-zaetr)*(zaeqsn-zaeqso)
        zaero(jl,jk,3) = (1._wp-zaetr)*(zaequn-zaequo)
        zaero(jl,jk,4) = zaetr*cstbga*pdp(jl,jk)
      END DO
    END DO
    !
    ! ---optical thickness is not negative
    !
    aero_tanre(:,:,:) = MAX(zaero(:,:,1:naer),EPSILON(1.0_wp)) 

  END FUNCTION aero_tanre

  !-----------------------------------------------------------------------------
  !>
  !! @brief: Retrieves aerosol properties from TANRE climatology
  !
  PURE FUNCTION aero_tanre3(krow,kproma,ppd_hl,pp_fl,tk_hl)

    INTEGER, INTENT (IN) :: krow, kproma
    REAL(wp), INTENT(IN) ::   &
         ppd_hl(:,:),         & !< pressure thickness [hPa]
         pp_fl(:,:),          & !< pressure on full levels [hPa]
         tk_hl(:,:)             !< temperature at half levels [K]

    REAL(wp) :: aero_tanre3(kproma,SIZE(ppd_hl,2),naer)

    INTEGER  :: jl, jk
    REAL(wp) :: zaer(kproma,SIZE(ppd_hl,2),naer)
    REAL(wp) :: z1n(kproma), z1o(kproma), z1, zlnd, zurb, zsea, zdes

    z1o(:) = 1.0_wp

    DO jk = 1, SIZE(ppd_hl,2)
      DO jl = 1, kproma
        zlnd = caeopl*zael_x(jl,krow)*cvdael(jk+1) 
        zlnd = zlnd - caeopl*zael_x(jl,krow)*cvdael(jk)
        zurb = caeopu*zaeu_x(jl,krow)*cvdaeu(jk+1)
        zurb = zurb - caeopu*zaeu_x(jl,krow)*cvdaeu(jk)
        zsea = caeops*zaes_x(jl,krow)*cvdaes(jk+1)
        zsea = zsea - caeops*zaes_x(jl,krow)*cvdaes(jk)
        zdes = caeopd*zaed_x(jl,krow)*cvdaed(jk+1)
        zdes = zdes - caeopd*zaed_x(jl,krow)*cvdaed(jk)
        !
        ! --- set scaling factors based on height/temperature gradients
        !
        IF (pp_fl(jl,jk) < 999.0_wp) THEN
          z1= 1.0_wp 
        ELSE
          z1n(jl) = z1o(jl)*(MIN(1.0_wp,tk_hl(jl,jk)/tk_hl(jl,jk+1)))**ctrpt
          z1      = SQRT(z1n(jl)*z1o(jl))
          z1o(jl) = z1n(jl)
        END IF

        zaer(jl,jk,1) = (1.0_wp-z1)*(ctrbga*ppd_hl(jl,jk)  + zlnd + zdes)
        zaer(jl,jk,2) = (1.0_wp-z1)*zsea 
        zaer(jl,jk,3) = (1.0_wp-z1)*zurb
        zaer(jl,jk,4) = z1*cstbga*ppd_hl(jl,jk)
      END DO
    END DO
    !
    ! ---optical thickness is not negative
    !
    aero_tanre3(:,:,:) = MAX(zaer(:,:,1:naer),EPSILON(1.0_wp)) 

  END FUNCTION aero_tanre3
  !-----------------------------------------------------------------------------
  !>
  !! @brief sets up aerosol properties from TANRE climatology
  !! 
  !! @remarks
  !!   This routine computes the values of the coefficients cvdae[xx] where xx
  !!   xx can be s (for sea), l (for land), u (for urban) or d (for desert) of
  !!   a surface-normalised vertical distribution of aerosols' optical depths
  !!   from the argument petah which is the normalized pressure height or 
  !!   vertical coordinate at klevp1 levels. 
  !
  SUBROUTINE setup_aero_tanre(petah)

    REAL(wp), INTENT(in) :: petah(:)

    REAL(wp), PARAMETER  :: zx1 = 8434.0_wp/1000.0_wp
    REAL(wp), PARAMETER  :: zx2 = 8434.0_wp/3000.0_wp

    INTEGER :: klevp1
    INTEGER :: jmm, imm, imnc, imns, jnn, jl, nglon, nglat, nlof, jlat, jglat

    REAL(wp), DIMENSION(21) :: zfaes,zfael,zfaeu,zfaed
    REAL(wp), DIMENSION(66) :: zalp

    REAL(wp) :: zsin, &
         zcos1, zcos2, zcos3, zcos4, zcos5, zcos6, zcos7, zcos8, zcos9, zcos10,&
         zsin1, zsin2, zsin3, zsin4, zsin5, zsin6, zsin7, zsin8, zsin9, zsin10
!LT         
!    REAL(wp) :: zaes(ldc% nglon, ldc% nglat)
!    REAL(wp) :: zael(ldc% nglon, ldc% nglat)
!    REAL(wp) :: zaeu(ldc% nglon, ldc% nglat)
!    REAL(wp) :: zaed(ldc% nglon, ldc% nglat)

    EXTERNAL  legtri !< for the Legendre transform at a given latitude 
    !
    ! 1.0 Fix constants used in setting aerosol properties
    !--------------------------------
    !
    ! --- get vertical extension and allocate space
    !
    klevp1 = SIZE(petah)
    IF (.NOT. ALLOCATED(cvdaes)) ALLOCATE(cvdaes(klevp1))
    IF (.NOT. ALLOCATED(cvdael)) ALLOCATE(cvdael(klevp1))
    IF (.NOT. ALLOCATED(cvdaeu)) ALLOCATE(cvdaeu(klevp1))
    IF (.NOT. ALLOCATED(cvdaed)) ALLOCATE(cvdaed(klevp1))
    !
    ! --- set up constants
    !
    cvdaes(:) = petah(:)**MAX(1._wp,8434._wp/1000._wp)
    cvdael(:) = petah(:)**MAX(1._wp,8434._wp/1000._wp)
    cvdaeu(:) = petah(:)**MAX(1._wp,8434._wp/1000._wp)
    cvdaed(:) = petah(:)**MAX(1._wp,8434._wp/3000._wp)

    IF (petah(1) == 0.0_wp) THEN
      cvdaes(1) = 0.0_wp
      cvdael(1) = 0.0_wp
      cvdaeu(1) = 0.0_wp
      cvdaed(1) = 0.0_wp
    END IF

    caeops = 0.05_wp
    caeopl = 0.2_wp
    caeopu = 0.1_wp
    caeopd = 1.9_wp
    ctrpt  = 30._wp
    caeadm = 2.6E-10_wp
    caeadk(1) =  0.3876E-03_wp
    caeadk(2) =  0.6693E-02_wp
    caeadk(3) =  0.8563E-03_wp
!LT
    !
    ! 2.0 Set up transformations from T10 input to global climatology
    !--------------------------------
!    nglon = ldc% nglon   !< local number of longitudes
!    nglat = ldc% nglat   !< local number of latitudes
!
!    DO jlat = 1, nglat
!
!      jglat = ldc% glat(jlat)
!      nlof  = ldc% glon(jlat)  !< longitude offset to global field
!      zsin  = 0.5_wp*gl_twomu(jglat)
      !
      ! --- Prepare and execute Legendre transform at given latitute
      !
!      CALL legtri(zsin,11,zalp)

      zfaes(:) = 0.0_wp
      zfael(:) = 0.0_wp
      zfaeu(:) = 0.0_wp
      zfaed(:) = 0.0_wp
      imm  = 0
      imnc = 0
      imns = 0
!      DO jmm = 1, 11
!        imm = imm + 1
!        DO jnn = jmm, 11
!          imnc = imnc + 1
!          zfaes(imm) = zfaes(imm) + zalp(imnc)*caesc(imnc)
!          zfael(imm) = zfael(imm) + zalp(imnc)*caelc(imnc)
!          zfaeu(imm) = zfaeu(imm) + zalp(imnc)*caeuc(imnc)
!          zfaed(imm) = zfaed(imm) + zalp(imnc)*caedc(imnc)
!        END DO
!        IF (jmm/=1) THEN
!          imm = imm + 1
!          DO jnn = jmm, 11
!            imns = imns + 1
!            zfaes(imm) = zfaes(imm) + zalp(imns+11)*caess(imns)
!            zfael(imm) = zfael(imm) + zalp(imns+11)*caels(imns)
!            zfaeu(imm) = zfaeu(imm) + zalp(imns+11)*caeus(imns)
!            zfaed(imm) = zfaed(imm) + zalp(imns+11)*caeds(imns)
!          END DO
!        END IF
!      END DO
      !
      ! --- Fourier transform
      !
!      DO jl = 1, nglon
!        zcos1 = coslon(jl+nlof)
!        zsin1 = sinlon(jl+nlof)
!        zcos2 = zcos1*zcos1 - zsin1*zsin1
!        zsin2 = zsin1*zcos1 + zcos1*zsin1
!        zcos3 = zcos2*zcos1 - zsin2*zsin1
!        zsin3 = zsin2*zcos1 + zcos2*zsin1
!        zcos4 = zcos3*zcos1 - zsin3*zsin1
!        zsin4 = zsin3*zcos1 + zcos3*zsin1
!        zcos5 = zcos4*zcos1 - zsin4*zsin1
!        zsin5 = zsin4*zcos1 + zcos4*zsin1
!        zcos6 = zcos5*zcos1 - zsin5*zsin1
!        zsin6 = zsin5*zcos1 + zcos5*zsin1
!        zcos7 = zcos6*zcos1 - zsin6*zsin1
!        zsin7 = zsin6*zcos1 + zcos6*zsin1
!        zcos8 = zcos7*zcos1 - zsin7*zsin1
!        zsin8 = zsin7*zcos1 + zcos7*zsin1
!        zcos9 = zcos8*zcos1 - zsin8*zsin1
!        zsin9 = zsin8*zcos1 + zcos8*zsin1
!        zcos10= zcos9*zcos1 - zsin9*zsin1
!        zsin10= zsin9*zcos1 + zcos9*zsin1
!
!        zaes(jl,jlat) = zfaes(1) +                                             &
!             2.0_wp*(zfaes(2)*zcos1+zfaes(3)*zsin1+zfaes(4)*zcos2+             &
!             zfaes(5)*zsin2+zfaes(6)*zcos3+zfaes(7)*zsin3+zfaes(8)*zcos4+      &
!             zfaes(9)*zsin4+zfaes(10)*zcos5+zfaes(11)*zsin5+zfaes(12)*zcos6+   &
!             zfaes(13)*zsin6+zfaes(14)*zcos7+zfaes(15)*zsin7+zfaes(16)*zcos8+  &
!             zfaes(17)*zsin8+zfaes(18)*zcos9+zfaes(19)*zsin9+zfaes(20)*zcos10+ &
!             zfaes(21)*zsin10)
!        zael(jl,jlat) = zfael(1) +                                             &
!             2.0_wp*(zfael(2)*zcos1+zfael(3)*zsin1+zfael(4)*zcos2+             &
!             zfael(5)*zsin2+zfael(6)*zcos3+zfael(7)*zsin3+zfael(8)*zcos4+      &
!             zfael(9)*zsin4+zfael(10)*zcos5+zfael(11)*zsin5+zfael(12)*zcos6+   &
!             zfael(13)*zsin6+zfael(14)*zcos7+zfael(15)*zsin7+zfael(16)*zcos8+  &
!             zfael(17)*zsin8+zfael(18)*zcos9+zfael(19)*zsin9+zfael(20)*zcos10+ &
!             zfael(21)*zsin10)
!        zaeu(jl,jlat) = zfaeu(1) +                                             &
!             2.0_wp*(zfaeu(2)*zcos1+zfaeu(3)*zsin1+zfaeu(4)*zcos2+             &
!             zfaeu(5)*zsin2+zfaeu(6)*zcos3+zfaeu(7)*zsin3+zfaeu(8)*zcos4+      &
!             zfaeu(9)*zsin4+zfaeu(10)*zcos5+zfaeu(11)*zsin5+zfaeu(12)*zcos6+   &
!             zfaeu(13)*zsin6+zfaeu(14)*zcos7+zfaeu(15)*zsin7+zfaeu(16)*zcos8+  &
!             zfaeu(17)*zsin8+zfaeu(18)*zcos9+zfaeu(19)*zsin9+zfaeu(20)*zcos10+ &
!             zfaeu(21)*zsin10)
!        zaed(jl,jlat) = zfaed(1) +                                             &
!             2.0_wp*(zfaed(2)*zcos1+zfaed(3)*zsin1+zfaed(4)*zcos2+             &
!             zfaed(5)*zsin2+zfaed(6)*zcos3+zfaed(7)*zsin3+zfaed(8)*zcos4+      &
!             zfaed(9)*zsin4+zfaed(10)*zcos5+zfaed(11)*zsin5+zfaed(12)*zcos6+   &
!             zfaed(13)*zsin6+zfaed(14)*zcos7+zfaed(15)*zsin7+zfaed(16)*zcos8+  &
!             zfaed(17)*zsin8+zfaed(18)*zcos9+zfaed(19)*zsin9+zfaed(20)*zcos10+ &
!             zfaed(21)*zsin10)
!      END DO
!    END DO

!    CALL reorder (zaes_x, zaes)
!    CALL reorder (zael_x, zael)
!    CALL reorder (zaeu_x, zaeu)
!    CALL reorder (zaed_x, zaed)

    initialized = .TRUE.

  END SUBROUTINE setup_aero_tanre
  !-----------------------------------------------------------------------------
  !>
  !! @brief Cleans up Tanre Aerosol
  !
  SUBROUTINE cleanup_aero_tanre

    IF (ALLOCATED(cvdaes)) DEALLOCATE(cvdaes)
    IF (ALLOCATED(cvdael)) DEALLOCATE(cvdael)
    IF (ALLOCATED(cvdaeu)) DEALLOCATE(cvdaeu)
    IF (ALLOCATED(cvdaed)) DEALLOCATE(cvdaed)

  END SUBROUTINE cleanup_aero_tanre

END MODULE mo_aero_tanre
