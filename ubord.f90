!            *****************************
!            **  SUBROUTINEUbord CLIMITE**
!            *****************************
!*****************************************************************************
!     *--------------------------------------
!     * Ce sous programme donne les valeurs des conditions 
!     * aux limites dirichlet aux sommets du bord
!******************************************************************************
SUBROUTINE UBORD(chxgb,temps)
  !--------
  ! Modules
  !--------
  USE longr
  USE imprime
  USE parmmage
  USE fsource

  IMPLICIT NONE

  !--------------------------
  ! Declaration des arguments
  !--------------------------
  Integer, intent(in) :: chxgb
  Real(kind=long)     :: temps
  !--------------------------------
  ! Declaration des variables locales
  !----------------------------------
  CHARACTER(len=6)      :: oldprf
  integer               :: is, iseg

  !-------------------
  ! Debut du programme
  !-------------------
  oldprf = prefix
  prefix = 'UBORD '

  !------
  ! Corps
  !------
  Gb = 0.D0
  ! On calcule en premier gbord sur les sommets du bord
  DO is = 1, Nbs
     Select case ( NtypS(is) )
     case (dirichlet)
        Gb(is) = gbord(temps, coordS(1,is), coordS(2,is),chxgb )
     case(0,neumann)
        !! On fait rien
     case default
        print*,'probleme TypSeg'
        stop
     End Select
  END DO
  ! Ensuite, on calcule gbord sur les milieux des segments du bord
  DO iseg = 1, Nseg 
     Select case ( NtypSeg(iseg) )
     case (dirichlet)
        Gb(iseg+Nbs) = gbord(temps, CoordMseg(1,iseg), CoordMseg(2,iseg),chxgb )
     case(0,neumann)
        !! On fait rien
     case default
        print*,'probleme TypSeg'
        stop
     End Select
  END DO
  print*,'max Gb dans ubord = ', maxval(Gb)
  print*,'min Gb dans ubord = ', minval(Gb)

  !-----------------
  ! Fin du programme
  !-----------------
  prefix = oldprf

  RETURN
END SUBROUTINE UBORD
