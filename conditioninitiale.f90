!            ************************************
!            **  SUBROUTINE conditioninitiale  **
!            ************************************
!****************************************************************
!     *-!!!!!!!!!@@@@@!!@@@@@!!@@@@@!!@@@@@!!!!!!!!!!!-*
!     *-!!!!!!!!!!!!!@!!@!!!@!!!!!!@!!@!!!@!!!!!!!!!!!-*
!     *-!!!!!!!!!@@@@@!!@!!!@!!@@@@@!!@!!!@!!!!!!!!!!!-*
!     *-!!!!!!!!!@!!!!!!@!!!@!!@@!!!!!@!!!@!!!!!!!!!!!-*
!     *-!!!!!!!!!@@@@@!!@@@@@!!@@@@@!!@@@@@!!!!!!!!!!!-*
!------------------------------------------------------------
!     * Ce sous programme lit le fichier uread
!     * le fichier d'entree
!****************************************************************
SUBROUTINE  conditioninitiale(CC,EE,UU,VV,Ndim)
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
  Integer , intent(in) :: Ndim
  REAL(kind=long), dimension(ndim), intent(inout) ::  CC,EE,UU,VV
  !----------------------------------
  ! Declaration des variables locales
  !----------------------------------
  CHARACTER(len=6)          :: oldprf
  REAL(kind = long)         :: xmin, xmax, ymin, ymax, uinit, vinit
  CHARACTER(len=len_buffer) :: buffer
  INTEGER                   :: Ndivu, Ndivv, jt, ii, is, ip, i, sizecdtU
  REAL(kind=long), dimension(2) :: S1, S2, S3, v0, v1, v2
  REAL(kind=long) :: dist, distmin
  LOGICAL :: trouver

  !-------------------
  ! Debut du programme
  !-------------------
  oldprf = prefix
  prefix = 'CdInit'

  !----------------------------------
  !  Lire les conditions initiales
  !----------------------------------
  UU = 0.D0
  CC = 0.D0
  EE = 0.D0
  VV = 0.D0
  
  !!----------------------------------------
  !! lecture de la condition initiale pour U
  !!----------------------------------------
  buffer=lireligne(uread)
  READ(buffer, *, err = 10) TypecondinitialeU
  CALL prvari (uprint, 'TypecondinitialeU       : ', TypecondinitialeU)
  !!
  Select case ( TypecondinitialeU )
     !
  case(0) 
     buffer=lireligne(uread)
     READ(buffer, *, err = 10) Ndivu
     CALL prvari (uprint, 'NdivU                : ', Ndivu)
     print*,'size(UU) cdt initiale =',size(UU)
     !!
     DO ii = 1, Ndivu 
        buffer = lireligne(uread)
        READ(buffer, *, err = 10) xmin, xmax, ymin, ymax, uinit
        print*, xmin, xmax, ymin, ymax, uinit
        DO is = 1, NsInt
           if (CoordS(1,is)> xmin .and. CoordS(1,is) < xmax &
                & .and. CoordS(2,is) > ymin .and. CoordS(2,is) < ymax) &
                & UU(is)= uinit
        END DO
        DO jt = 1, Nbt
           if (CoordK(1,jt) > xmin .and. CoordK(1,jt) < xmax &
                & .and. CoordK(2,jt) > ymin .and. CoordK(2,jt) < ymax) &
                & UU(jt + NsInt)= uinit
        END DO
     END DO

  case(1) 
     ! condition initiale faite par gbord
     !
     DO is = 1, NsInt
        UU(is) = gbord(0.D0, CoordS(1,is), CoordS(2,is), choixgb)
     END DO
     DO jt = 1, Nbt
        UU(jt + NsInt) = gbord(0.D0, CoordK(1,jt), CoordK(2,jt), choixgb)
     END DO

  case(2) ! condition initiale aleartoire 

     call random_number(UU)

  case(3) ! condition initiale faite par une condition existante
     DO jt = 1, NbInc
        READ(cdtinitialeU,"(f30.20)", err = 12) UU(jt)
     END DO
     print*,'Donne le numéro à partir duquel on numérote les résultat d''afficahge :'
     print*,'Ce numéro est accessible du fichier UPRINT derniere ligne :'
     read *, affichage
  end Select

  !!----------------------------------------
  !! lecture de la condition initiale pour V
  !!----------------------------------------
  buffer=lireligne(uread)
  READ(buffer, *, err = 10) TypecondinitialeV
  CALL prvari (uprint, 'TypecondinitialeU       : ', TypecondinitialeV)
  !!
  Select case ( TypecondinitialeV )
     !
  case(0) 
     buffer=lireligne(uread)
     READ(buffer, *, err = 10) Ndivv
     CALL prvari (uprint, 'NdivV                : ', Ndivv)
     print*,'size(VV) cdt initiale =',size(VV)
     !!
     DO ii = 1, Ndivv 
        buffer = lireligne(uread)
        READ(buffer, *, err = 10) xmin, xmax, ymin, ymax, vinit
        print*, xmin, xmax, ymin, ymax, vinit
        DO is = 1, NsInt
           if (CoordS(1,is)> xmin .and. CoordS(1,is) < xmax &
                & .and. CoordS(2,is) > ymin .and. CoordS(2,is) < ymax) &
                & VV(is)= vinit
        END DO
        DO jt = 1, Nbt
           if (CoordK(1,jt) > xmin .and. CoordK(1,jt) < xmax &
                & .and. CoordK(2,jt) > ymin .and. CoordK(2,jt) < ymax) &
                & VV(jt + NsInt)= vinit
        END DO
     END DO

  case(1) 
     ! condition initiale faite par gbord
     !
     DO is = 1, NsInt
        VV(is) = gbord(0.D0, CoordS(1,is), CoordS(2,is), choixgb)
     END DO
     DO jt = 1, Nbt
        VV(jt + NsInt) = gbord(0.D0, CoordK(1,jt), CoordK(2,jt), choixgb)
     END DO

  case(2) ! condition initiale aleartoire 

     call random_number(VV)

  case(3) ! condition initiale faite par une condition existante 
     DO jt = 1, NbInc
        READ(cdtinitialeV,"(f30.20)", err = 12) VV(jt)
     END DO
  end Select

  !=====================
  ! lecture des points 
  !====================

  buffer=lireligne(uread)
  READ(buffer, *, err = 10)  Nb_points
  CALL prvari (uprint, 'Nb_points       : ', Nb_points)
  Allocate ( Xpoints(Nb_points),Ypoints(Nb_points) ) 
  DO ii = 1, Nb_points
     buffer=lireligne(uread)
     READ(buffer, *, err = 10)  Xpoints(ii), Ypoints(ii)
     CALL prvarr (uprint, 'x points       : ',Xpoints(ii) )
     CALL prvarr (uprint, 'y points       : ',Ypoints(ii) )
  End DO
  !! -----------------------------------------
  !! calcul de la maille ou se trouve le point 
  !! -----------------------------------------
  Allocate (Num_points(Nb_points))
  Num_points = -99

  xmin = minval( coordS(1,:) )
  xmax = maxval( coordS(1,:) )
  ymin = minval( coordS(2,:) )
  ymax = maxval( coordS(2,:) )
  !! min du domaine et max du domaine
  ! ---------------------------------
  write(uprint,*)'limite domaine', xmin, xmax, ymin, ymax

  DO ii = 1, Nb_points
     i=1
     ip=1
     distmin = sqrt( (CoordS(1,i)-xpoints(ii))**2 + (CoordS(2,i)-ypoints(ii))**2 )
     do i=2, NsInt
        dist = sqrt( (CoordS(1,i)-xpoints(ii))**2 + (CoordS(2,i)-ypoints(ii))**2 )
        if (dist < distmin) then
           distmin = dist
           ip=i
        endif
     enddo
     Num_points(ii)=ip
  ENDDO

  write(uprint,*)'Num_points', Num_points
  do ii=1, Nb_points
     Write(uprint,*)'les points', CoordS(1,Num_points(ii)),  CoordS(2,Num_points(ii))
  end do


  !==============
  !==============
  buffer=lireligne(uread)
  READ(buffer, *, err = 10)  CLgauche
  CALL prvari (uprint, 'CLgauche                        : ', CLgauche)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10)  CLdroite
  CALL prvari (uprint, 'CLdroite                        : ', CLdroite)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10)  CLbas
  CALL prvari (uprint, 'CLbas                           : ', CLbas)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10)  CLhaut
  CALL prvari (uprint, 'CLhaut                          : ', CLhaut)

  close(uread)
  !-----------------
  ! Fin du programme
  !-----------------
  prefix = oldprf
  RETURN

10 PRINT*,"Erreur dans l'entree des parametres CI"
12 PRINT*,"Erreur dans la lecture de la condition initiale"
  close (uread)
  STOP


END SUBROUTINE CONDITIONINITIALE
