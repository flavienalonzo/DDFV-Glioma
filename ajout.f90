!            **************************
!            **  SUBROUTINE AJOUT**
!            **************************
!******************************************************************************
!     *    Ecole Centrale de Nantes
!     *--------------------------------------
!     * Ce sous programme assemble la contribution de la
!     * la matrice elementaire locale Aloc du triangle jt
!     * Ajouter a la ligne I et la colone J la valeur Coef
!     *---------------------------
!     * Description des parametres
!     *---------------------------
!     * parm1  (E ,I*4)  ...
!     * ...
!     *----------------------------
!     * Fonctionnement du programme
!     *----------------------------
!     * Utilise les modules longr, imprime, parmmmage
!******************************************************************************
SUBROUTINE AJOUT(II, JJ, coefa,  A)
  !--------
  ! Modules
  !--------
  USE longr
  USE imprime
  USE parmmage

  IMPLICIT NONE

  !--------------------------
  ! Declaration des arguments
  !--------------------------
  INTEGER, INTENT(in)                  :: II , JJ
  REAL(kind=long), INTENT(in)          :: coefa
  TYPE(MatCreux)                        :: A
  !----------------------------------
  ! Declaration des variables locales
  !----------------------------------
  CHARACTER(len=6)      :: oldprf
  INTEGER               :: j
  LOGICAL               :: trouve

  !-------------------
  ! Debut du programme
  !-------------------
  oldprf = prefix
  prefix = 'AJOUT'

  !------
  ! Corps
  !------
  ! IndPL (II)  : le numero du premier elt non nul de la ligne II
  trouve = .FALSE.
  !print*,' A%IndPL( II )=', A%IndPL( II )
  !print*,' A%IndPL( II+1 )=', A%IndPL( II+1 )
  DO j = A%IndPL( II ) ,   A%IndPL( II+1 )  -1
     IF ( A%Indc(j) == JJ) THEN
        A%TMat(j) = A%TMat(j) + coefa 
        trouve = .TRUE.
        EXIT
     ENDIF
  ENDDO
  IF (.NOT.trouve ) THEN
     PRINT*, ' probleme d''assemblage'
     STOP
  ENDIF
  !------------
  ! Impressions
  !------------
!!$  IF (iprint >= 2) THEN
!!$     CALL prvari(uprint, ' Numero triangle = ', jt) 
!!$     WRITE (uprint,*) ( (matloc(iloc,jloc), jloc=1,3),iloc=1,3)
!!$  END IF

  !-----------------
  ! Fin du programme
  !-----------------
  prefix = oldprf

  RETURN
END SUBROUTINE ajout








