MODULE intgradc
  USE longr 
  USE imprime
  USE parmmage
  USE intmatvec
  IMPLICIT NONE

  INTERFACE OPERATOR (/)
     MODULE PROCEDURE gradconj
  END INTERFACE
CONTAINS
  FUNCTION gradconj(b , A)
    !     *--------------------------------------------------
    !     * Ce sous programme resoud le systeme lineaire 
    !     *    A x = b 
    !     * par la methode des gradients conjugues 
    !     * le resultat de l'operation x = b/A 
    !--------------------------------------------------------
    !--------------------------
    ! Declaration des arguments
    !--------------------------
    TYPE(MatCreux), INTENT(in)                                :: A
    REAL(kind=long), DIMENSION( SIZE(A%IndPL) -1 ),INTENT(in) :: b
    REAL(kind=long), DIMENSION( SIZE(A%IndPL) -1 )            :: gradconj
    !----------------------------------
    ! Declaration des variables locales
    !----------------------------------
    CHARACTER(len=6)      :: oldprf
    REAL(kind=long), DIMENSION( SIZE(A%IndPL) -1 ) ::  R, P, Q, X
    REAL (kind=long)                       :: alfa, bta, seuil, seuil9
    INTEGER                                :: ngrad, maxitergrad
    !-------------------
    ! Debut du programme
    !-------------------
    oldprf = prefix
    prefix = 'GRADCG '
    !
    R = 0.D0 ; P = 0.D0 ; Q = 0.D0
    ngrad = 0
    !---------------------------
    ! Initialisation de X  ?? attention il faut donner X0 quand il le faut
    !---------------------------
    X = 0.D0       
    R = b - A*X
    P = R

    seuil = DOT_PRODUCT(R, R)
    Maxitergrad = 2000
    DO ngrad = 1, Maxitergrad
                            ! nombre d'iteration
       Q =  A* P 

       IF (SQRT(seuil) < 1.D-10) EXIT

       alfa = seuil/DOT_PRODUCT(Q,P)
       X = X + alfa * P
       R = R - alfa * Q
       seuil9 = DOT_PRODUCT(R,R)
       bta = seuil9/seuil
       seuil = seuil9 
       P = R + bta * P

    END DO
    
     gradconj = X

    WRITE(6,*)'GRAD  NITER = ',ngrad,' SEUIL = ', SQRT(seuil) 
    !------------
    ! Impressions
    !------------
    !   if (iprint >= 3) then
    !      call prvarr(uprint, 'seuil',seuil)
    !   end if

    !-----------------
    ! Fin du programme
    !-----------------
    prefix = oldprf
    !
    RETURN
  END FUNCTION  gradconj


 function bigradient(A,b, X0,tolerence) result(X)
   !     *--------------------------------------------------
    !     * Ce sous programme resoud le systeme lineaire 
    !     *    A x = b 
    !     * par la methode des bigradients conjugues 
    !     * le resultat de l'operation x = b/A 
    !--------------------------------------------------------    
    ! Modules
    !--------
 
    implicit none

    !--------------------------
    ! Declaration des arguments
    !--------------------------
    TYPE(MatCreux), INTENT(in)                                :: A
    REAL(kind=long), DIMENSION( SIZE(A%IndPL) -1 ),INTENT(in) :: b
    REAL(kind=long), DIMENSION( SIZE(A%IndPL) -1 )            :: X , X0

     real(kind=long), intent(in) :: tolerence

    !----------------------------------
    ! Declaration des variables locales
    !----------------------------------
    character(len=6) :: oldprf
    real(kind=long), dimension (size(b))   :: Rst, Prj, Q, Rstb, Prjb,Qb
    real(kind=long)  :: alfa, bta,seuil,prscal, prscal9
    integer          :: nbigr, nbitermax
    ! ----------------------------
    ! Signification des variables
    ! ---------------------------
    ! Rst : reste = B - A * X
    ! Prj : projection
    ! Q   := A * Prj
    ! Qb  := transposee (A) * Prjb
    ! Rstb:= B - transposee (A) * X
    ! prscal : designe un produit scalaire
    !-------------------
    ! Debut du programme
    !-------------------
    oldprf = prefix
    prefix = 'BIGRAD '
    !
    Rst = 0. ;  Prj = 0. ; Q = 0.
    Rstb = 0. ; Prjb = 0. ; Qb = 0.
    nbigr = 0
    !
    X= X0
    Rst = B - A*X
    Rstb = B - transposee_mat_vec(A,X)
    Prj = Rst
    Prjb = Rstb
    prscal = DOT_PRODUCT(Rstb,Rst)

    nbitermax = size(X0) + 10000

    DO  nbigr = 1, nbitermax
       seuil = sqrt(DOT_PRODUCT(Rst, Rst))

       if (seuil < tolerence) exit

       Q =  A*Prj
       alfa = prscal/DOT_PRODUCT(Q,Prjb)
       X = X + alfa * Prj
       Rst = Rst - alfa * Q
       Qb = transposee_mat_vec(A,Prjb) 
       Rstb = Rstb - alfa * Qb
       prscal9 = DOT_PRODUCT(Rstb,Rst)
       bta = prscal9/prscal
       prscal = prscal9 
       Prj = Rst + bta * Prj
       Prjb = Rstb + bta * Prjb
       !================
       ! Impression
       !================
!!$       write(6,*)'nbigr =', nbigr
!!$       write(6,*)'prscal=',prscal
!!$       write(6,*)'alfa=',alfa
!!$       write(6,*)'bta=',bta
!!$       write(6,*)'seuil=',seuil , 'tolerence = ', tolerence   
    END DO
    write(6,*)'nbigr =', nbigr, 'seuil bigrad=',seuil
    !------------
    ! Impressions
    !------------
    !   if (iprint >= 3) then
    !      call prvarr(uprint, 'seuil',seuil)
    !   end if

    !-----------------
    ! Fin du programme
    !-----------------
    prefix = oldprf

  end function bigradient









END MODULE intgradc
