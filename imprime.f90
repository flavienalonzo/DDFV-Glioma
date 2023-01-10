module imprime
  use longr
contains
  !----------------------------------------------------------------------------
  ! Module des variables communes aux utilitaires d'impression.
  ! Ce module contient deux subroutines :
  !
  !  1>        subroutine prchar    : PRint CHaine
  !  2>        subroutine prvarr    : PRint VARiable Real
  !  3>        subroutine prvari    : PRint VARiable Integer
  !**************************************************************************
  !            **************************
  !            **  SUBROUTINE PRCHAR   **
  !            **************************
  !**************************************************************************
  !     * Sous programmes : PRCHAR, PRVARI, PRVARR
  !     *-----------------------------------------
  !     * Ces sous programmes impriment un texte suivi
  !     *     d'un entier pour 'PRVARI'
  !     *     d'un reel   pour 'PRVARR'
  !     *-----------------
  !     * Sequence d'appel
  !     *-----------------
  !     * call prchar(uprint,texte)
  !     * call prvari(uprint,texte,entier)
  !     * call prvarr(uprint,texte,reel)
  !     *---------------------------
  !     * Description des parametres
  !*************************************************************************
  FUNCTION lireligne(uread)

    !---------------------------------------------------------------------
    ! Role : lit une ligne du fichier d'unite uread en ignorant les lignes 
    !  vides les lignes commencant par * ou #
    !---------------------------------------------------------------------

    IMPLICIT NONE

    CHARACTER(len=80)    :: lireligne
    INTEGER, INTENT(in)  :: uread

    CHARACTER(len=80)    :: buffer
    INTEGER :: i

    boucle :  DO
       READ(uread,'(A)', err = 1001, END = 1000) buffer
  !!!! print*,'buffer = ',buffer
       buffer=ADJUSTL(buffer)

       SELECT CASE(buffer(:1))
       CASE ('#') 
          CONTINUE
       CASE ('*')
          CONTINUE
       CASE (' ')
          CONTINUE
       CASE ('!')
          CONTINUE
       CASE default
          EXIT
       END SELECT
    END DO boucle
    i=1
    DO 
       IF (buffer(i:i)=='=') THEN
          buffer(i:i)=' '
          EXIT
       END IF
       buffer(i:i)=' '
       i=i+1
    END DO

10  lireligne = TRIM(buffer)
    !print*,'sortie'
    RETURN
1000 PRINT *, 'LIRELIGNE : fin de fichier '
    STOP
1001 PRINT *, 'LIRELIGNE : Erreur, contenu du tampon ; ', buffer
    STOP
  END FUNCTION lireligne

  FUNCTION lireligne2(uread)

    !---------------------------------------------------------------------
    ! Role : lit une ligne du fichier d'unite uread en ignorant les lignes 
    !  vides les lignes commencant par * ou #
    !---------------------------------------------------------------------

    IMPLICIT NONE

    CHARACTER(len=80)    :: lireligne2
    INTEGER, INTENT(in)  :: uread

    CHARACTER(len=80)    :: buffer
    INTEGER :: i

    boucle :  DO
       READ(uread,'(A)', err = 1001, END = 1000) buffer
 !! print*,'buffer = ',buffer
       buffer=ADJUSTL(buffer)

       SELECT CASE(buffer(:1))
       CASE ('#') 
          CONTINUE
       CASE ('*')
          CONTINUE
       CASE (' ')
          CONTINUE
       CASE default
          EXIT
       END SELECT
    END DO boucle

10  lireligne2 = TRIM(buffer)
    !print*,'sortie'
    RETURN
1000 PRINT *, 'LIRELIGNE2 : fin de fichier '
    STOP
1001 PRINT *, 'LIRELIGNE2 : Erreur, contenu du tampon ; ', buffer
    STOP
  END FUNCTION lireligne2


  !!//////////////////////////

  subroutine prchar(uprint,texte)
    !========
    ! Modules
    !========
    !use imprime
    !===========================
    ! Declarations des arguments
    !===========================
    integer,          intent(in) :: uprint
    character(len=*), intent(in) :: texte
    !===================
    ! Debut du programme
    !===================
    write(uprint,fmt='(1X,A8,1X,A)') prefix, texte
    !=================
    ! Fin du programme
    !=================
    return
  end subroutine prchar
  !***************************************************************************
  !            **************************
  !            **  SUBROUTINE PRVARI   **
  !            **************************
  !***************************************************************************
  !     *
  !     *----------------------------------
  !     * Voir les commentaires dans PRCHAR
  !     *----------------------------------
  !***************************************************************************
  subroutine prvari(uprint,texte,entier)
    !========
    ! Modules
    !========
    !use imprime
    !===========================
    ! Declarations des arguments
    !===========================
    integer,          intent(in) :: uprint, entier
    character(len=*), intent(in) :: texte
    !===================
    ! Debut du programme
    !===================
    write(uprint,fmt='(1X,A8,1X,A,I10)') prefix, texte, entier
    !=================
    ! Fin du programme
    !=================
    return
  end subroutine prvari
  !*************************************************************************
  !            **************************
  !            **  SUBROUTINE PRVARR   **
  !            ****************************************************
  !     * 
  !     *----------------------------------
  !     * Voir les commentaires dans PRCHAR
  !     *----------------------------------
  !***************************************************************************
  subroutine prvarr(uprint,texte,reel)
    !========
    ! Modules
    !========
    !===========================
    ! Declarations des arguments
    !===========================
    integer,          intent(in) :: uprint
    real(kind=long),     intent(in) :: reel
    character(len=*), intent(in) :: texte
    !===================
    ! Debut du programme
    !===================
    write(uprint,fmt='(1X,A8,1X,A,E15.7)') prefix, texte, reel
    !=================
    ! Fin du programme
    !=================
    return
  end subroutine prvarr





end module imprime






