!            **************************
!            **  FUNCTION tenseurAMATLOC  **
!            **************************
!******************************************************************************
!     *-!!!!!!!!!@@@@@!!@@@@@!!@@@@@!!@@@@@!!!!!!!!!!!-*
!     *-!!!!!!!!!!!!!@!!@!!!@!!!!!!@!!@!!!@!!!!!!!!!!!-*
!     *-!!!!!!!!!@@@@@!!@!!!@!!@@@@@!!@!!!@!!!!!!!!!!!-*
!     *-!!!!!!!!!@!!!!!!@!!!@!!@@!!!!!@!!!@!!!!!!!!!!!-*
!     *-!!!!!!!!!@@@@@!!@@@@@!!@@@@@!!@@@@@!!!!!!!!!!!-*
!------------------------------------------------------------
! calcul sur chaque triangle le tenseur de diffusion 
! SxxK, SxyK, SyyK
!
!******************************************************************************
SUBROUTINE tenseur(p)
  !--------
  ! Modules
  !--------
  USE longr
  USE imprime
  USE parmmage

  IMPLICIT NONE

  !----------------------------------
  ! Declaration des variables locales
  !----------------------------------
  CHARACTER(len=6)                  :: oldprf
  Integer, intent(in)               :: p
  INTEGER                           :: iseg,js,ks, iloc, jloc
  INTEGER                           :: Kvois, Lvois
  REAL(kind=long), DIMENSION(2)     :: z
  REAL(kind=long), DIMENSION(4)     :: x, y
  REAL(kind=long)                   :: GEO,d1,d2,d3 
  REAL(kind=long) :: sally,a,kscalaire,merci1,merci2,merci4

  !-------------------
  ! Debut du programme
  !-------------------
  oldprf = prefix
  prefix = 'TENSEUR'
  !------
  ! Corps
  !------
  ALLOCATE(SxxK(Nseg),SyyK(Nseg),SxyK(Nseg))
  
  DO iseg = 1, Nseg
     Kvois = NumTVoisSeg(1,iseg)
     js = Nuseg(1,iseg) ; ks = Nuseg(2,iseg)
     !!
     x(1) = coordK(1,Kvois) ; y(1) = coordK(2,Kvois)
     x(2) = coordS(1,js) ; y(2) = coordS(2,js)
     x(3) = coordS(1,ks) ; y(3) = coordS(2,ks)
     !!
     z(1) = (x(1)+x(2)+x(3))/ 3.
     z(2) = (y(1)+y(2)+y(3))/ 3.
     IF (NombvoisSeg(iseg) == 2) THEN
        Lvois = NumTVoisSeg(2,iseg)
        js = NuSeg(1,iseg) ; ks = NuSeg(2,iseg)
        !!
        x(1) = coordK(1,Lvois) ; y(1) = coordK(2,Lvois)
        !! 
        z(1) = (3*z(1) + x(1))/4.
        z(2) = (3*z(2) + y(1))/4.
        !!
     END IF
     !------------------------
     ! anisotropie heterogene
     !------------------------
     sally= 1. /(z(1)**2 +z(2)**2)
     merci1 = sally*( deltau*(z(1)**2)+ z(2)**2)
     merci2 = -(1 - deltau)*sally* z(1)*z(2)
     merci4 = sally*(deltau*(z(2)**2)+ z(1)**2)
     !-------------------------------------------
     Select case (p)
     case(1) 
        !-----------
        ! S = Id
        !----------
        SxxK(iseg)= 1.D0
        SyyK(iseg) = 1.D0
        SxyK(iseg)= 0.D0
     case(2)
        !--------------------------
        !anisotropie exponentielle
        !--------------------------
        SxxK(iseg)= kscalaire(z(1),z(2),choixkscalaireu)
        SyyK(iseg) = SxxK(iseg)
        SxyK(iseg)= 0.D0

     case(3)
        !--------------------
        ! anisotropie moderee
        !--------------------
        SxxK(iseg) = deltaxu 
        SxyK(iseg) = deltaxyu
        SyyK(iseg) = deltayu

     case(4)
        !------------------------
        ! anisotropie heterogene
        !------------------------
        sally= 1. /(z(1)**2 +z(2)**2)
        SxxK(iseg) = sally*( deltau*(z(1)**2)+ z(2)**2)
        SxyK(iseg) = -(1 - deltau)*sally* z(1)*z(2)
        SyyK(iseg) = sally*(deltau*(z(2)**2)+ z(1)**2)
     end select
  end DO

  !------------
  ! Impressions
  !------------
  IF (iprint >= 5) THEN
     CALL prvari(uprint, ' Numero segment = ', iseg) 
     DO iloc=1,3
        WRITE (uprint,110) SxxK
     END DO
  END IF
110 FORMAT (3(E16.9, 2X))
  !-----------------
  ! Fin du programme
  !-----------------
  prefix = oldprf

END SUBROUTINE tenseur
