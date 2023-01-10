module intmatvec
  use longr 
  USE imprime
  USE parmmage
  implicit none

  interface operator (*)
     module procedure mat_vec
  end interface
contains
  FUNCTION mat_vec(A, X) result(Y)
    !     *--------------------------------------
    !     * Ce sous programme calcule sur un triangle donne
    !     * la matrice elementaire locale pour l'aproximation P_1
    !     * MATLOC fourni une matrice 3X3

    IMPLICIT NONE

    !--------------------------
    ! Declaration des arguments
    !--------------------------
    TYPE(MatCreux), intent(in)                        :: A
    REAL(kind=long), DIMENSION( SIZE(A%IndPL) -1 ),intent(in) :: X
    REAL(kind=long), DIMENSION( SIZE(A%IndPL) -1 )            :: Y
    !----------------------------------
    ! Declaration des variables locales
    !----------------------------------
    CHARACTER(len=6)      :: oldprf
    INTEGER               :: i, j, ndim

    !-------------------
    ! Debut du programme
    !-------------------
    oldprf = prefix
    prefix = 'MATVEC'

    !------
    ! Corps
    !------
    Y = 0.D0
    ndim = SIZE(A%IndPL) -1 
    DO i = 1, ndim
       DO j = A%IndPL(i) ,  A%IndPL(i+1) -1 
          Y(i) =  Y(i) + A%TMat( j )*X(A%Indc(j))
       END DO
    ENDDO

    !------------
    ! Impressions
    !------------
!!$  IF (iprint >= 2) THEN
!!$     WRITE(uprint,*)' SOLUTION DE MATVEC'
!!$     WRITE(uprint,100) (matvec(i), i=1,ndim)
!!$  END IF
100 FORMAT (10(E10.3, 2X))
    !-----------------
    ! Fin du programme
    !-----------------
    prefix = oldprf

    RETURN
  END FUNCTION mat_vec

  FUNCTION transposee_mat_vec(A, X) result(Y)
    !     *--------------------------------------
    !     * Ce sous programme calcule le produit 
    !     * de la transposee d'une matrice par un vecteur 
    !     * 
    !     * La matrice A est stockée sous formme creuse
    !     * 

    IMPLICIT NONE

    !--------------------------
    ! Declaration des arguments
    !--------------------------
    TYPE(MatCreux), intent(in)                                :: A
    REAL(kind=long), DIMENSION( SIZE(A%IndPL) -1 ),intent(in) :: X
    REAL(kind=long), DIMENSION( SIZE(A%IndPL) -1 )            :: Y
    !----------------------------------
    ! Declaration des variables locales
    !----------------------------------
    CHARACTER(len=6)      :: oldprf
    INTEGER               :: i, j 

    !-------------------
    ! Debut du programme
    !-------------------
    oldprf = prefix
    prefix = 'MATVEC'

    !------
    ! Corps
    !------
    Y = 0.D0
    DO i = 1, SIZE(A%IndPL) -1 
       DO j = A%IndPL(i) ,  A%IndPL(i+1) -1 
          Y( A%IndC(j) ) =  Y( A%IndC(j) ) + A%TMat( j ) * X( i )
       END DO
    ENDDO

    !------------
    ! Impressions
    !------------
!!$  IF (iprint >= 2) THEN
!!$     WRITE(uprint,*)' SOLUTION DE MATVEC'
!!$     WRITE(uprint,100) (matvec(i), i=1,ndim)
!!$  END IF
100 FORMAT (10(E10.3, 2X))
    !-----------------
    ! Fin du programme
    !-----------------
    prefix = oldprf

    RETURN
  END FUNCTION transposee_mat_vec

end module intmatvec
