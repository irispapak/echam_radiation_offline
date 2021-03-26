module yoesrta16

use mo_kind, only : wp

implicit none

save

!     -----------------------------------------------------------------
!*    ** *YOESRTA16* - SRTM COEFFICIENTS FOR INTERVAL 16
!     BAND 16:  2600-3250 cm-1 (low - H2O,CH4; high - CH4)
!     -----------------------------------------------------------------

integer, parameter :: jpg=16, ng16 = 16, ngs15 = 0

real(wp) :: ka(9,5,13,jpg) 
real(wp) :: kb(5,13:59,jpg)
real(wp) :: selfref(10,jpg),forref(3,jpg)
real(wp) :: sfluxref(jpg)
real(wp) :: rayl            ,strrat1
integer :: layreffr

real(wp) :: kac(9,5,13,ng16),absa(585,ng16)
real(wp) :: kbc(5,13:59,ng16),absb(235,ng16)
real(wp) :: selfrefc(10,ng16),forrefc(3,ng16)
real(wp) :: sfluxrefc(ng16)

!EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)), (KB(1,13,1),ABSB(1,1))
equivalence (kac(1,1,1,1),absa(1,1)), (kbc(1,13,1),absb(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
!     M. J. IACONO          AER             12/09/03

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! KA      : REAL 
! KB      : REAL    
! SELFREF : REAL 
! FORREF  : REAL   
! SFLUXREF: REAL
! RAYL    : REAL 
! STRRAT1 : REAL
! LAYREFFR: INTEGER
! KAC     : REAL     Reduced g-point array for KA
! KBC     : REAL     Reduced g-point array for KB
! SELFREFC: REAL     Reduced g-point array for SELFREF
! FORREFC : REAL     Reduced g-point array for FORREF
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
!     -----------------------------------------------------------------
end module yoesrta16

module yoesrta17

use mo_kind, only : wp

implicit none

save

!     -----------------------------------------------------------------
!*    ** *YOESRTA17* - SRTM COEFFICIENTS FOR INTERVAL 17
!     BAND 17:  3250-4000 cm-1 (low - H2O,CO2; high - H2O,CO2)
!     -----------------------------------------------------------------

integer, parameter ::  jpg = 16, ng17 = 16, ngs16=16

real(wp) :: ka(9,5,13,jpg)  
real(wp) :: kb(5,5,13:59,jpg)
real(wp) :: selfref(10,jpg)  ,forref(4,jpg)
real(wp) :: sfluxref(jpg,5)
real(wp) :: rayl              ,strrat
integer :: layreffr

real(wp) :: kac(9,5,13,ng17),absa(585,ng17)
real(wp) :: kbc(5,5,13:59,ng17),absb(1175,ng17)
real(wp) :: selfrefc(10,ng17),forrefc(4,ng17)
real(wp) :: sfluxrefc(ng17,5)

!EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)), (KB(1,1,13,1),ABSB(1,1))
equivalence (kac(1,1,1,1),absa(1,1)), (kbc(1,1,13,1),absb(1,1))
!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
!     M. J. IACONO          AER             12/09/03

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! KA      : REAL 
! KB      : REAL    
! SELFREF : REAL 
! FORREF  : REAL  
! SFLUXREF: REAL  
! RAYL    : REAL
! STRRAT  : REAL
! LAYREFFR: INTEGER 
! KAC     : REAL     Reduced g-point array for KA
! KBC     : REAL     Reduced g-point array for KB
! SELFREFC: REAL     Reduced g-point array for SELFREF
! FORREFC : REAL     Reduced g-point array for FORREF
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
!     -----------------------------------------------------------------
end module yoesrta17

module yoesrta18

use mo_kind, only : wp

implicit none

save

!     -----------------------------------------------------------------
!*    ** *YOESRTA18* - SRTM COEFFICIENTS FOR INTERVAL 16
!     BAND 18:  4000-4650 cm-1 (low - H2O,CH4; high - CH4)
!     -----------------------------------------------------------------

integer, parameter :: jpg = 16, ng18 = 16, ngs17=32

real(wp) :: ka(9,5,13,jpg) 
real(wp) :: kb(5,13:59,jpg)
real(wp) :: selfref(10,jpg),forref(3,jpg)
real(wp) :: sfluxref(jpg,9)
real(wp) :: rayl            ,strrat

integer :: layreffr

real(wp) :: kac(9,5,13,ng18) ,absa(585,ng18)
real(wp) :: kbc(5,13:59,ng18),absb(235,ng18)
real(wp) :: selfrefc(10,ng18),forrefc(3,ng18)
real(wp) :: sfluxrefc(ng18,9)

!EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)), (KB(1,13,1),ABSB(1,1))
equivalence (kac(1,1,1,1),absa(1,1)), (kbc(1,13,1),absb(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
!     M. J. IACONO          AER             12/09/03

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! KA      : REAL 
! KB      : REAL    
! SELFREF : REAL 
! FORREF  : REAL   
! SFLUXREF: REAL 
! RAYL    : REAL
! STRRAT  : REAL
! LAYREFFR: INTEGER
! KAC     : REAL     Reduced g-point array for KA
! KBC     : REAL     Reduced g-point array for KB
! SELFREFC: REAL     Reduced g-point array for SELFREF
! FORREFC : REAL     Reduced g-point array for FORREF
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
!     -----------------------------------------------------------------
end module yoesrta18

module yoesrta19

use mo_kind, only : wp

implicit none

save

!     -----------------------------------------------------------------
!*    ** *YOESRTA19* - SRTM COEFFICIENTS FOR INTERVAL 19
!     BAND 19:  4650-5150 cm-1 (low - H2O,CO2; high - CO2)
!     -----------------------------------------------------------------

integer, parameter :: jpg = 16, ng19 = 16

real(wp) :: ka(9,5,13,jpg) 
real(wp) :: kb(5,13:59,jpg)
real(wp) :: selfref(10,jpg),forref(3,jpg)
real(wp) :: sfluxref(jpg,9)
real(wp) :: rayl            ,strrat
integer :: layreffr

real(wp) :: kac(9,5,13,ng19) ,absa(585,ng19)
real(wp) :: kbc(5,13:59,ng19),absb(235,ng19)
real(wp) :: selfrefc(10,ng19),forrefc(3,ng19)
real(wp) :: sfluxrefc(ng19,9)

!EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)), (KB(1,13,1),ABSB(1,1))
equivalence (kac(1,1,1,1),absa(1,1)), (kbc(1,13,1),absb(1,1))
!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
!     M. J. IACONO          AER             12/09/03

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! KA      : REAL 
! KB      : REAL    
! SELFREF : REAL 
! FORREF  : REAL
! SFLUXREF: REAL  
! RAYL    : REAL  
! STRRAT  : REAL
! LAYREFFR: INTEGER
! KAC     : REAL     Reduced g-point array for KA
! KBC     : REAL     Reduced g-point array for KB
! SELFREFC: REAL     Reduced g-point array for SELFREF
! FORREFC : REAL     Reduced g-point array for FORREF
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
!     -----------------------------------------------------------------
end module yoesrta19

module yoesrta20

use mo_kind, only : wp

implicit none

save

!     -----------------------------------------------------------------
!*    ** *YOESRTA20* - SRTM COEFFICIENTS FOR INTERVAL 20
!     BAND 20:  5150-6150 cm-1 (low - H2O; high - H2O)
!     -----------------------------------------------------------------

integer, parameter :: jpg = 16, ng20 = 16

real(wp) :: ka(5,13,jpg)   
real(wp) :: kb(5,13:59,jpg)
real(wp) :: selfref(10,jpg),forref(4,jpg)
real(wp) :: sfluxref(jpg)  ,absch4(jpg)
real(wp) :: rayl
integer :: layreffr

real(wp) :: kac(5,13,ng20)   ,absa(65,ng20)
real(wp) :: kbc(5,13:59,ng20),absb(235,ng20)
real(wp) :: selfrefc(10,ng20),forrefc(4,ng20)
real(wp) :: sfluxrefc(ng20)  ,absch4c(ng20)

!EQUIVALENCE (KA(1,1,1),ABSA(1,1)), (KB(1,13,1),ABSB(1,1))
equivalence (kac(1,1,1),absa(1,1)), (kbc(1,13,1),absb(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
!     M. J. IACONO          AER             12/09/03

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! KA      : REAL 
! KB      : REAL    
! SELFREF : REAL 
! FORREF  : REAL  
! SFLUXREF: REAL
! ABSCH4  : REAL  
! RAYL    : REAL
! LAYREFFR: INTEGER
! KAC     : REAL     Reduced g-point array for KA
! KBC     : REAL     Reduced g-point array for KB
! SELFREFC: REAL     Reduced g-point array for SELFREF
! FORREFC : REAL     Reduced g-point array for FORREF
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
! ABSCH4C : REAL     Reduced g-point array for ABSCH4
!     -----------------------------------------------------------------
end module yoesrta20

module yoesrta21

use mo_kind, only : wp

implicit none

save

!     -----------------------------------------------------------------
!*    ** *YOESRTA21* - SRTM COEFFICIENTS FOR INTERVAL 21
!     BAND 21:  6150-7700 cm-1 (low - H2O,CO2; high - H2O,CO2)
!     -----------------------------------------------------------------

integer, parameter :: jpg = 16, ng21 = 16

real(wp) :: ka(9,5,13,jpg)   
real(wp) :: kb(5,5,13:59,jpg)
real(wp) :: selfref(10,jpg)  ,forref(4,jpg)
real(wp) :: sfluxref(jpg,9)
real(wp) :: rayl              ,strrat
integer :: layreffr

real(wp) :: kac(9,5,13,ng21)   ,absa(585,ng21)
real(wp) :: kbc(5,5,13:59,ng21),absb(1175,ng21)
real(wp) :: selfrefc(10,ng21)  ,forrefc(4,ng21)
real(wp) :: sfluxrefc(ng21,9)

equivalence (ka(1,1,1,1),absa(1,1)), (kb(1,1,13,1),absb(1,1))
!EQUIVALENCE (KAC(1,1,1,1),ABSA(1,1)), (KBC(1,1,13,1),ABSB(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
!     M. J. IACONO          AER             12/09/03

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! KA      : REAL 
! KB      : REAL    
! SELFREF : REAL 
! FORREF  : REAL    
! SFLUXREF: REAL
! RAYL    : REAL
! STRRAT  : REAL
! LAYREFFR: INTEGER
! KAC     : REAL     Reduced g-point array for KA
! KBC     : REAL     Reduced g-point array for KB
! SELFREFC: REAL     Reduced g-point array for SELFREF
! FORREFC : REAL     Reduced g-point array for FORREF
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
!     -----------------------------------------------------------------
end module yoesrta21

module yoesrta22

use mo_kind, only : wp

implicit none

save

!     -----------------------------------------------------------------
!*    ** *YOESRTA22* - SRTM COEFFICIENTS FOR INTERVAL 22
!     BAND 22:  7700-8050 cm-1 (low - H2O,O2; high - O2)
!     -----------------------------------------------------------------

integer, parameter :: jpg = 16, ng22 = 16

real(wp) :: ka(9,5,13,jpg) 
real(wp) :: kb(5,13:59,jpg)
real(wp) :: selfref(10,jpg),forref(3,jpg)
real(wp) :: sfluxref(jpg,9)
real(wp) :: rayl            ,strrat
integer :: layreffr

real(wp) :: kac(9,5,13,ng22) ,absa(585,ng22)
real(wp) :: kbc(5,13:59,ng22),absb(235,ng22)
real(wp) :: selfrefc(10,ng22),forrefc(3,ng22)
real(wp) :: sfluxrefc(ng22,9)

!EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)), (KB(1,13,1),ABSB(1,1))
equivalence (kac(1,1,1,1),absa(1,1)), (kbc(1,13,1),absb(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
!     M. J. IACONO          AER             12/09/03

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! KA      : REAL 
! KB      : REAL    
! SELFREF : REAL 
! FORREF  : REAL    
! SFLUXREF: REAL
! RAYL    : REAL
! STRRAT  : REAL
! LAYREFFR: INTEGER
! KAC     : REAL     Reduced g-point array for KA
! KBC     : REAL     Reduced g-point array for KB
! SELFREFC: REAL     Reduced g-point array for SELFREF
! FORREFC : REAL     Reduced g-point array for FORREF
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
!     -----------------------------------------------------------------
end module yoesrta22

module yoesrta23

use mo_kind, only : wp

implicit none

save

!     -----------------------------------------------------------------
!*    ** *YOESRTA23* - SRTM COEFFICIENTS FOR INTERVAL 23
!     BAND 23:  8050-12850 cm-1 (low - H2O; high - nothing)
!     -----------------------------------------------------------------

integer, parameter :: jpg = 16, ng23 = 16

real(wp) :: ka(5,13,jpg)   
real(wp) :: selfref(10,jpg),forref(3,jpg)
real(wp) :: sfluxref(jpg)  ,rayl(jpg)
real(wp) :: givfac
integer :: layreffr

real(wp) :: kac(5,13,ng23)   ,absa(65,ng23)
real(wp) :: selfrefc(10,ng23),forrefc(3,ng23)
real(wp) :: sfluxrefc(ng23)  ,raylc(ng23)

!EQUIVALENCE (KA(1,1,1),ABSA(1,1))
equivalence (kac(1,1,1),absa(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
!     M. J. IACONO          AER             12/09/03

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! KA      : REAL 
! SELFREF : REAL 
! FORREF  : REAL    
! SFLUXREF: REAL
! RAYL    : REAL
! GIVFAC  : REAL
! LAYREFFR: INTEGER
! KAC     : REAL     Reduced g-point array for KA
! SELFREFC: REAL     Reduced g-point array for SELFREF
! FORREFC : REAL     Reduced g-point array for FORREF
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
! RAYLC   : REAL     Reduced g-point array for RAYL
!     -----------------------------------------------------------------
end module yoesrta23

module yoesrta24

use mo_kind, only : wp

implicit none

save

!     -----------------------------------------------------------------
!*    ** *YOESRTA24* - SRTM COEFFICIENTS FOR INTERVAL 24
!     BAND 24: 12850-16000 cm-1 (low - H2O,O2; high - O2)
!     -----------------------------------------------------------------

integer, parameter :: jpg = 16, ng24 = 16

real(wp) :: ka(9,5,13,jpg) 
real(wp) :: kb(5,13:59,jpg)
real(wp) :: selfref(10,jpg),forref(3,jpg)
real(wp) :: sfluxref(jpg,9)
real(wp) :: abso3a(jpg), abso3b(jpg), rayla(jpg,9), raylb(jpg)
real(wp) :: strrat
integer :: layreffr

real(wp) :: kac(9,5,13,ng24) ,absa(585,ng24)
real(wp) :: kbc(5,13:59,ng24),absb(235,ng24)
real(wp) :: selfrefc(10,ng24),forrefc(3,ng24)
real(wp) :: sfluxrefc(ng24,9)
real(wp) :: abso3ac(ng24), abso3bc(ng24), raylac(ng24,9), raylbc(ng24)

!EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)), (KB(1,13,1),ABSB(1,1))
equivalence (kac(1,1,1,1),absa(1,1)), (kbc(1,13,1),absb(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
!     M. J. IACONO          AER             12/09/03

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! KA      : REAL 
! KB      : REAL    
! SELFREF : REAL 
! FORREF  : REAL 
! SFLUXREF: REAL
! ABSO3A  : REAL
! ABSO3B  : REAL
! RAYLA   : REAL
! RAYLB   : REAL   
! STRRAT  : REAL
! LAYREFFR: INTEGER
! KAC     : REAL     Reduced g-point array for KA
! KBC     : REAL     Reduced g-point array for KB
! SELFREFC: REAL     Reduced g-point array for SELFREF
! FORREFC : REAL     Reduced g-point array for FORREF
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
! ABSO3AC : REAL     Reduced g-point array for ABSO3A
! ABSO3BC : REAL     Reduced g-point array for ABSO3B
! RAYLAC  : REAL     Reduced g-point array for RAYLA
! RAYLBC  : REAL     Reduced g-point array for RAYLB
!     -----------------------------------------------------------------
end module yoesrta24

module yoesrta25

use mo_kind, only : wp

implicit none

save

!     -----------------------------------------------------------------
!*    ** *YOESRTA25* - SRTM COEFFICIENTS FOR INTERVAL 25
!     BAND 25: 16000-22650 cm-1 (low - H2O; high - nothing)
!     -----------------------------------------------------------------

integer, parameter :: jpg = 16, ng25 = 16

real(wp) :: ka(5,13,jpg) 
real(wp) :: sfluxref(jpg)
real(wp) :: rayl(jpg), abso3a(jpg), abso3b(jpg)
integer :: layreffr

real(wp) :: kac(5,13,ng25) ,absa(65,ng25)
real(wp) :: sfluxrefc(ng25)
real(wp) :: raylc(ng25), abso3ac(ng25), abso3bc(ng25)

!EQUIVALENCE (KA(1,1,1),ABSA(1,1))
equivalence (kac(1,1,1),absa(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
!     M. J. IACONO          AER             12/09/03

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! KA      : REAL 
! SFLUXREF: REAL
! RAYL    : REAL
! ABSO3A  : REAL
! ABSO3B  : REAL
! LAYREFFR: INTEGER
! KAC     : REAL     Reduced g-point array for KA
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
! RAYLC   : REAL     Reduced g-point array for RAYL
! ABSO3AC : REAL     Reduced g-point array for ABSO3A
! ABSO3BC : REAL     Reduced g-point array for ABSO3B
!     -----------------------------------------------------------------
end module yoesrta25

module yoesrta26

use mo_kind, only : wp

implicit none

save

!     -----------------------------------------------------------------
!*    ** *YOESRTA26* - SRTM COEFFICIENTS FOR INTERVAL 26
!     BAND 26: 22650-29000 cm-1 (low - nothing; high - nothing)
!     -----------------------------------------------------------------

integer, parameter :: jpg = 16, ng26 = 16

real(wp) :: sfluxref(jpg), rayl(jpg)

real(wp) :: sfluxrefc(ng26), raylc(ng26)

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
!     M. J. IACONO          AER             12/09/03

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! SFLUXREF: REAL    
! RAYL    : REAL 
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
! RAYLC   : REAL     Reduced g-point array for RAYL
!     -----------------------------------------------------------------
end module yoesrta26

module yoesrta27

use mo_kind, only : wp

implicit none

save

!     -----------------------------------------------------------------
!*    ** *YOESRTA27* - SRTM COEFFICIENTS FOR INTERVAL 27
!     BAND 27: 29000-38000 cm-1 (low - O3; high - O3)
!     -----------------------------------------------------------------

integer, parameter :: jpg = 16, ng27 = 16

real(wp) :: ka(5,13,jpg)   
real(wp) :: kb(5,13:59,jpg)
real(wp) :: sfluxref(jpg)  ,rayl(jpg)
real(wp) :: scalekur
integer :: layreffr

real(wp) :: kac(5,13,ng27)   ,absa(65,ng27)
real(wp) :: kbc(5,13:59,ng27),absb(235,ng27)
real(wp) :: sfluxrefc(ng27)  ,raylc(ng27)

!EQUIVALENCE (KA(1,1,1),ABSA(1,1)), (KB(1,13,1),ABSB(1,1))
equivalence (kac(1,1,1),absa(1,1)), (kbc(1,13,1),absb(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
!     M. J. IACONO          AER             12/09/03

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! FRACREFA: REAL    
! KA      : REAL 
! KB      : REAL    
! SFLUXREF: REAL 
! RAYL    : REAL    
! SCALEKUR: REAL
! LAYREFFR:INTEGER
! KAC     : REAL     Reduced g-point array for KA
! KBC     : REAL     Reduced g-point array for KB
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
! RAYLC   : REAL     Reduced g-point array for RAYL
!     -----------------------------------------------------------------
end module yoesrta27

module yoesrta28

use mo_kind, only : wp

implicit none

save

!     -----------------------------------------------------------------
!*    ** *YOESRTA28* - SRTM COEFFICIENTS FOR INTERVAL 28
!     BAND 28: 38000-50000 cm-1 (low - O3, O2; high - O3, O2)
!     -----------------------------------------------------------------

integer, parameter :: jpg = 16, ng28 = 16

real(wp) :: ka(9,5,13,jpg)   
real(wp) :: kb(5,5,13:59,jpg)
real(wp) :: sfluxref(jpg,5)
real(wp) :: rayl              ,strrat
integer :: layreffr

real(wp) :: kac(9,5,13,ng28)   ,absa(585,ng28)
real(wp) :: kbc(5,5,13:59,ng28),absb(1175,ng28)
real(wp) :: sfluxrefc(ng28,5)

!EQUIVALENCE (KA(1,1,1,1),ABSA(1,1)), (KB(1,1,13,1),ABSB(1,1))
equivalence (kac(1,1,1,1),absa(1,1)), (kbc(1,1,13,1),absb(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
!     M. J. IACONO          AER             12/09/03

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! KA      : REAL 
! KB      : REAL    
! SFLUXREF: REAL 
! RAYL    : REAL    
! STRRAT  : REAL
! LAYREFFR: INTEGER
! KAC     : REAL     Reduced g-point array for KA
! KBC     : REAL     Reduced g-point array for KB
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
!     -----------------------------------------------------------------
end module yoesrta28

module yoesrta29

use mo_kind, only : wp

implicit none

save

!     -----------------------------------------------------------------
!*    ** *YOESRTA29* - SRTM COEFFICIENTS FOR INTERVAL 29
!     BAND 29:  820-2600 cm-1 (low - H2O; high - CO2)
!     -----------------------------------------------------------------

integer, parameter :: jpg = 16, ng29 = 16

real(wp) :: ka(5,13,jpg)   
real(wp) :: kb(5,13:59,jpg)
real(wp) :: selfref(10,jpg),forref(4,jpg)
real(wp) :: sfluxref(jpg)  ,absh2o(jpg)  , absco2(jpg)
real(wp) :: rayl
integer :: layreffr

real(wp) :: kac(5,13,ng29)   ,absa(65,ng29)
real(wp) :: kbc(5,13:59,ng29),absb(235,ng29)
real(wp) :: selfrefc(10,ng29),forrefc(4,ng29)
real(wp) :: sfluxrefc(ng29)  ,absh2oc(ng29)  , absco2c(ng29)

!EQUIVALENCE (KA(1,1,1),ABSA(1,1)), (KB(1,13,1),ABSB(1,1))
equivalence (kac(1,1,1),absa(1,1)), (kbc(1,13,1),absb(1,1))

!     -----------------------------------------------------------------
!        * E.C.M.W.F. PHYSICS PACKAGE ** RRTM SW RADIATION **

!     J.-J. MORCRETTE       E.C.M.W.F.      02/10/29
!     M. J. IACONO          AER             12/09/03

!  NAME     TYPE     PURPOSE
!  ----   : ----   : ---------------------------------------------------
! KA      : REAL 
! KB      : REAL    
! SELFREF : REAL 
! FORREF  : REAL 
! SFLUXREF: REAL
! ABSH2O  : REAL
! ABSCO2  : REAL   
! RAYL    : REAL
! LAYREFFR: INTEGER
! KAC     : REAL     Reduced g-point array for KA
! KBC     : REAL     Reduced g-point array for KB
! SELFREFC: REAL     Reduced g-point array for SELFREF
! FORREFC : REAL     Reduced g-point array for FORREF
!SFLUXREFC: REAL     Reduced g-point array for SFLUXREF
! ABSH2OC : REAL     Reduced g-point array for ABSH2O
! ABSCO2C : REAL     Reduced g-point array for ABSCO2
!     -----------------------------------------------------------------
end module yoesrta29

