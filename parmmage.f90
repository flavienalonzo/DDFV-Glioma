module parmmage
  !---------------------------------------------------------------------
  ! On definit ici les parametres du maillage
  !---------------------------------------------------------------------
  !  Nbs                    : nombres des sommets du maillages
  !  Nbt                    : nombres des triangles du maillage
  !  NsInt                  : nombres des sommets interieurs du domaine
  !  Nbord                  : nombres des sommets au bord du domaine
  !  Nseg                   : nombres des segments
  !  NsegInt                : nombres des segments interieurs 
  !  Nsegbord               : nombres des segments au bord
  !----------------------------------------------------------------------
  !  CoordS(2,Nbs)          : les coordonnees des sommets du maillage
  !  CoordK(2,Nbt)          : les coordonnees des centroides du maillage
  !  CoordMseg(2,Nseg)      : les coordonnees des milieux des segments
  !  NuSoK(3,Nbt)           : les numeros des sommets des triangles
  !  NumTVoisSeg            : numero des Triangles voisins du segment
  !----------------------------------------------------------------------

  use longr
  implicit none

  INTEGER                 ::  Nbs, Nbt,  Nseg, NsInt, Nbord, NsegInt, Nsegbord,NbInc
  REAL (kind=long),   DIMENSION(:,:), ALLOCATABLE :: CoordS, CoordK, CoordMseg, NsigK, NsigeKe
  INTEGER,            DIMENSION(:,:), ALLOCATABLE :: NuSoK, TypS
  INTEGER,            DIMENSION(:,:), ALLOCATABLE :: NuVoisK
  INTEGER,            DIMENSION(:)  , ALLOCATABLE :: Ntyps, Ntypseg
  INTEGER,            DIMENSION(:,:), ALLOCATABLE :: NuSeg, NumTVoisSeg
  REAL (kind=long),   DIMENSION(:),   ALLOCATABLE :: SxxK, SyyK, SxyK,Msig,Msige
  REAL (kind=long),   DIMENSION(:),   ALLOCATABLE :: SxxKv, SyyKv, SxyKv
  REAL (kind=long),   DIMENSION(:),   ALLOCATABLE :: uTKL,uTKeLe,ueta,uetaSSe,cTKL,cTKeLe,cetaSSe,ceta,&
                                                    &eTKL,eTKeLe,eetaSSe,eeta,vTKL,vTKeLe,vetaSSe,veta
  REAL (kind=long),   DIMENSION(:),   ALLOCATABLE :: AireK, AireD, AireDSommet
  INTEGER,            DIMENSION(:),   ALLOCATABLE :: Num_points, NombvoisSeg
  !-----------------------------------
  ! stockage pour le maillage régulier
  !------------------------------------
  INTEGER ::  Nbitertemps, typegeomaille, Maillage_Format


  REAL (kind = long), DIMENSION(:), ALLOCATABLE ::  Gb 

  !------------------------------------------------------------------- 
  ! Stockage du systeme lineaire a resoudre
  !-------------------------------------------------------------------
  ! IndPL      : INDice du Premier non coef non nul de la Ligne i
  ! Indc       : INDice de la colonne du coef (j)
  ! TMat       : les coefficients non nuls de la matrice stokes alors 
  !            : par ligne et puis ordre croissant par colonne
  ! F          : le second membre physique 
  ! Bg         : le second membre du systeme a resoudre Tmat * x = bg
  !            : Bg = F + Gb   ;  Bg la donnee au bord physique  
  !-------------------------------------------------------------------
  TYPE MatCreux
     INTEGER, DIMENSION(:), ALLOCATABLE :: IndPL, Indc
     REAL (kind = long), DIMENSION(:), ALLOCATABLE    :: TMat
     REAL (kind = long), DIMENSION(:), ALLOCATABLE    :: BG, F
  END TYPE MatCreux
end module parmmage


