module longr
  implicit none
  integer, parameter     :: long = 8
  character (len=8)      :: prefix
  integer                :: uread, uprint, usave, iprint, uplotvtk, upoints
  integer                :: unode,uele, uedge
  integer                :: CdtinitialeU, CdtinitialeV,NormU,NormV
  REAL (kind = long)     :: Choixmesh,Choixtri,selec_mesh
  real (kind = long)     :: dt,pi,epsilon,alpha,beta,ubar
  real (kind = long)     :: delta, deltau, deltax, deltay, deltaxy, deltaxu, deltayu, deltaxyu
  real (kind = long)     :: CoefDiffV, CoefDiffuAdeg, CoefTranspMu
  CHARACTER (len = 80)   :: nom_mesh,nom_ele, nom_node, nom_edge
  CHARACTER (len = 80)   :: nom_usave, FPLOTVTK, nom_upoints, dossiermaillage
  integer                :: choixkscalaire, choixkscalaireu, choixanis, choixanisu
  Integer                :: choixgb, choixpb, choixAdeg, ChoixMu, ChoixCdtBord,Choixdt,ChoixBord,gamma
  integer, parameter     :: transparent=3, Neumann=2, Dirichlet=1
  integer, parameter     :: len_buffer=80, triangle=3, quadrangle=4
  Integer                :: affichage
  Integer :: ntarr, Nb_points, kitermax
  REAL (kind=long), dimension(:), allocatable :: Tempstock, Xpoints, Ypoints
  Integer :: CLgauche, CLdroite, CLbas, CLhaut
  Integer :: TypecondinitialeU, TypecondinitialeV, Nb_iter_visit
  REAL (kind = long)     :: Tolerencegradient, TolerenceNewton, CFDT
  
end module longr
