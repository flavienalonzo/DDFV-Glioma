!            *************************
!            **  PROGRAM GLIOMADDFV  **
!            *************************
!******************************************************************************
!     *-!!!!!!!!!@@@@@!!@@@@@!!@@@@@!!@@@@@!!!!!!!!!!!-*
!     *-!!!!!!!!!!!!!@!!@!!!@!!!!!!@!!!!!!@!!!!!!!!!!!-*
!     *-!!!!!!!!!@@@@@!!@!!!@!!@@@@@!!@@@@@!!!!!!!!!!!-*
!     *-!!!!!!!!!@!!!!!!@!!!@!!@@!!!!!!!!!@!!!!!!!!!!!-*
!     *-!!!!!!!!!@@@@@!!@@@@@!!@@@@@!!@@@@@!!!!!!!!!!!-*
!------------------------------------------------------------
!     * Programme resoud le syst�me de chimiotaxie anisotropic
!       sur un domaine qlq par la methode des DDFV:
!
!      Probleme  trait� est le suivant : 
!      Dt U - Div( S(x) nabla A(u) ) + div(S(x) X(u) nabla V) = 0    dans omega
!      Dt V - Div( D(x) nabla V) = a U-b V
!               DA(u).n = D V.n=0        sur gamma_1
!               u=udb, v=vdb             sur gamma_2
!
!----------------
! La m�thode utilise le maillage dual de Vernoi. 
!******************************************************************************
PROGRAM GLIOMADDFV
  
  !========
  ! Modules
  !========
  USE longr
  USE parmmage
  USE imprime
  USE fsource
  USE plotvtkmod

  IMPLICIT NONE
  !==========================
  ! Declaration des tableaux
  !==========================
  TYPE(MatCreux)       :: AC, AE, AU, AV
  !==================================
  ! Declaration des variables locales
  !==================================
  INTEGER         :: jt, i, j, is, ii, niter, y, m
  REAL(kind=long), DIMENSION(:),   Pointer  :: Uold, U, Vold, V, Cold, C, Eold, E, Um, Vm, Cm, Em
  REAL(kind = long)                         :: Tempsactuel, erreur_iteratif
  REAL(kind=long), DIMENSION(:),ALLOCATABLE :: NormL1U, NormL1C, NormL1E, NormL1V
  logical         :: cvge
  CHARACTER(20) :: chaine
  real :: start, finish
  !===================
  ! Debut du programme
  !===================
  call cpu_time(start)
  prefix = 'KelSeg '
  !
  !! coef devant la derivee en temps
  !! si CFDT = 0, pb stationnaire
  pi = 4.D0*atan(1.)
  !!
  CALL init                    !* Initialisation
  print*, 'init_fini'

  CALL maillage                !* lecture du maillage
  print*, 'maillage_fini'

  CALL meshtools               !* Mise � jour du maillage
  print*, 'meshtools_fini'
  !------------------------------------
  ! Nombre d'inconnues dans le probleme
  !------------------------------------
  IF (ChoixCdtBord == Neumann) NsInt = NbS
  NbInc = NsInt + Nbt
  print*,'NbInc = ',NbInc
  ALLOCATE(Uold(NbInc), U(NbInc), Vold(NbInc), V(NbInc), Um(NbInc), Vm(NbInc))
  ALLOCATE(Cold(NbInc), C(NbInc), Eold(NbInc), E(NbInc), Cm(NbInc), Em(NbInc))
  ALLOCATE( Gb(Nbs + Nseg) ) ! Gb est de taille  Nbs - NsInt
  !---------------------------------------------
  ! Condition initiale et conditions aux limites
  !---------------------------------------------
  affichage = 0
  CALL Conditioninitiale(C,E,U,V,NbInc)
 
  !---------------------------------------
  ! Allocation de la structure des donn�es
  !---------------------------------------
  IF (Choixdt == 0) dt = 0.1*(maxval(Msig))**2
  CALL matrixinitDDFV(AC)
  CALL matrixinitDDFV(AE)
  CALL matrixinitDDFV(AU)
  CALL matrixinitDDFV(AV)
  !Uold = U
  !Vold = V
  !----------------------------------------------
  ! Assemblagle du second membre de AU%F physique
  !----------------------------------------------
  !CALL scmem(AU, choixpb, dt)
  !-----------------------------------------
  ! Coefficients de transmissibilites pour u, c, e et v
  !-----------------------------------------
  CALL transmis(Choixanisu,Choixanisc,Choixanise,Choixanisv)

 
  !! -------------
  !! Temps initial 
  !! -------------
  Tempsactuel = 0.D0
  CALL plot_vtk (U,'U00000','U')
  CALL plot_vtk (V,'V00000','V')
  CALL plot_vtk (C,'C00000','C')
  CALL plot_vtk (E,'E00000','E')
  !
  Allocate(NormL1U(nbitertemps),NormL1V(nbitertemps),NormL1C(nbitertemps),NormL1E(nbitertemps))
  NormL1U = 0.D0;  NormL1V = 0.D0;  NormL1C = 0.D0;  NormL1E = 0.D0
  !
  !! =============================================
  !! %%%%%%%%%%%%%% Boucle en Temps %%%%%%%%%%%%%%
  !! =============================================
  niter = 1
  DO WHILE (niter  <= nbitertemps)  
     print*,' Boucle en TEMPS iter = ', niter     
     print*,' Max solution calculee U   : ', MAXVAL(U)
     print*,' Min solution calculee U   : ', MINVAL(U)
     print*,' Max solution calculee V   : ', MAXVAL(V)
     print*,' Min solution calculee V   : ', MINVAL(V)
     print*,' Max solution calculee C   : ', MAXVAL(C)
     print*,' Min solution calculee C   : ', MINVAL(C)
     print*,' Max solution calculee E   : ', MAXVAL(E)
     print*,' Min solution calculee E   : ', MINVAL(E)
     Tempsactuel = Tempsactuel + dt
     ! -------------------------------------
     ! Norme L1 de U et V � chaque iteration
     ! -------------------------------------
     DO is = 1, NsInt
        NormL1U(niter) = NormL1U(niter) + AireDsommet(is)*U(is)/2
        NormL1V(niter) = NormL1V(niter) + AireDsommet(is)*V(is)/2
        NormL1C(niter) = NormL1C(niter) + AireDsommet(is)*C(is)/2
        NormL1E(niter) = NormL1E(niter) + AireDsommet(is)*E(is)/2
      END DO
     !
     DO jt = 1, Nbt
        ii = jt + NsInt
        NormL1U(niter) = NormL1U(niter) + AireK(jt)*U(ii)/2
        NormL1V(niter) = NormL1V(niter) + AireK(jt)*V(ii)/2
        NormL1C(niter) = NormL1C(niter) + AireK(jt)*C(ii)/2
        NormL1E(niter) = NormL1E(niter) + AireK(jt)*E(ii)/2
     END DO
     ! Perform the surgery
     if (Abs(Tempsactuel-1.D0)<=dt/2.D0 .and. Use_surgery .and. Tempsactuel>=1.D0) then
         call surge(U,C,E,V,NbInc)
     end if
     !IF (niter == 2) epsilon = 1.D-12
     !IF (MINVAL(U)>=-1.D-12 .AND. MINVAL(U)< 0) U = 0.D0
     !IF (MINVAL(V)>=-1.D-12 .AND. MINVAL(V)< 0) V = 0.D0
     Uold = U; Vold = V; Cold = C; Eold = E
     !----------------------------------------------
     ! Assemblagle du second membre de AU%F physique
     !----------------------------------------------
     !CALL scmem(AU, choixpb, Tempsactuel) 
     !----------------------------------------------
     CALL ubord(Choixgb,Tempsactuel)

     ! Algorithme iteratif en m
     Um = Uold; Vm = Vold; Cm = Cold; Em = Eold
     cvge = .false.; m = 0; mitermax = 1000
     Do while ((cvge .neqv. .true.) .AND. (m < mitermax))
        !--------------------------------------------------
        ! Resolution du systeme nonlineaire pour calculer U
        !--------------------------------------------------
        ! Methode de Newton
        CALL NewtoncDDFV(AC,Cold,C,Em,Um,NbInc,choixpb,Tempsactuel)
        CALL NewtonuDDFV(AU,Uold,U,Cm,Em,NbInc,choixpb,Tempsactuel)
        CALL NewtonvDDFV(AV,Vold,V,Cm,Em,Um,NbInc,choixpb,Tempsactuel)
        CALL NewtoneDDFV(AE,Eold,E,Um,Vm,NbInc,choixpb,Tempsactuel)
        
        ! Calcul de l'erreur entre deux iterees,
        erreur_iteratif = sqrt(dot_product(Cm-C,Cm-C)) + sqrt(dot_product(Em-E,Em-E))&
             & + sqrt(dot_product(Um-U,Um-U))+sqrt(dot_product(Vm-V,Vm-V))
         print*, 'erreur relative = ',erreur_iteratif
        If (erreur_iteratif < Tolerenceiterative) cvge = .true.
        Um = U; Vm = V; Cm = C; Em = E
        m = m+1
     end Do

     print*,' Boucle en TEMPS iter = ', niter     
     print*,' Max solution calculee U   : ', MAXVAL(U)
     print*,' Min solution calculee U   : ', MINVAL(U)
     print*,' Max solution calculee V   : ', MAXVAL(V)
     print*,' Min solution calculee V   : ', MINVAL(V)
     print*,' Max solution calculee C   : ', MAXVAL(C)
     print*,' Min solution calculee C   : ', MINVAL(C)
     print*,' Max solution calculee E   : ', MAXVAL(E)
     print*,' Min solution calculee E   : ', MINVAL(E)

     !------------------------------------------
     ! Erreur L^infinie entre U approchee et Uex
     !------------------------------------------
     !DO is = 1, NsInt
     !   Uexacte(is) = gbord(Tempsactuel,coordS(1,is), CoordS(2,is), Choixgb) ! solution exacte
     !END DO
     !DO jt = 1, Nbt
     !   Uexacte(jt+ NsInt) = gbord(Tempsactuel,CoordK(1,jt), CoordK(2,jt), Choixgb) ! solution exacte
     !END DO
     !WRITE(* ,*)' Erreur infinie U-Uexacte =  : ', MAXVAL(ABS(U - Uexacte))
     !WRITE(* ,*)' Erreur infinie relative  =  : ', MAXVAL(ABS(U - Uexacte))/MAXVAL(ABS(Uexacte))

     !-------------------
     ! Affichage resultat
     !-------------------
     ! 
     IF (mod(niter+affichage,Nb_iter_visit) == 0) then
        Select case ( TypecondinitialeU )
           !
        case (3)
           write(chaine,'(i5.5)')niter+affichage
           CALL plot_vtk (U,'U'//trim(chaine),'U')
           CALL plot_vtk (C,'C'//trim(chaine),'C')
           CALL plot_vtk (E,'E'//trim(chaine),'E')
           CALL plot_vtk (V,'V'//trim(chaine),'V')
        case default
           write(chaine,'(i5.5)')niter
           CALL plot_vtk (U,'U'//trim(chaine),'U')
           CALL plot_vtk (C,'C'//trim(chaine),'C')
           CALL plot_vtk (E,'E'//trim(chaine),'E')
           CALL plot_vtk (V,'V'//trim(chaine),'V')
        END Select
     END IF
        
     niter = niter + 1
     !
  END DO
  !!
  CALL plot_vtk (U,'Ufinal','U')
  CALL plot_vtk (C,'Cfinal','C')
  CALL plot_vtk (E,'Efinal','E')
  CALL plot_vtk (V,'Vfinal','V')
  !!
  close (cdtinitialeU)
  close (cdtinitialeV)
  OPEN (cdtinitialeU,file = 'cdtinitialeU',status='unknown')
  OPEN (cdtinitialeV,file = 'cdtinitialeV',status='unknown')
  DO is = 1, NbInc 
     WRITE(cdtinitialeU,*) U(is)
     WRITE(cdtinitialeV,*) V(is)
  END DO
  !!
  CALL prvari (uprint, 'affichage       : ', niter+affichage-1)
  !!
  !DO i = 1, nbitertemps 
  !   WRITE(NormU,211) (i-1)*dt, NormL1U(i)
  !   WRITE(NormC,211) (i-1)*dt, NormL1C(i)
  !   WRITE(NormE,211) (i-1)*dt, NormL1E(i)
  !   WRITE(NormV,211) (i-1)*dt, NormL1V(i)
  !END DO

  !print*, AU%IndPL
  !print*, AU%Tmat
!!$  !------------------------------------------
!!$  ! Erreur L^infinie entre U approch�e et Uex
!!$  !-------------------------------------------
!!$  DO is = 1, NsInt
!!$     Uexacte(is) = gbord(dt,coordS(1,is),  CoordS(2,is), Choixgb) ! solution exacte
!!$  END DO
!!$  DO jt = 1, Nbt
!!$     Uexacte(jt+ NsInt) = gbord(dt,CoordK(1,jt),  CoordK(2,jt), Choixgb) ! solution exacte
!!$  END DO
!!$  WRITE(* ,*)' Erreur infinie V-Uexacte =  : ', MAXVAL(ABS(V - Uexacte))
!!$  WRITE(* ,*)' Erreur infinie relative  =  : ', MAXVAL(ABS(V - Uexacte))/MAXVAL(ABS(Uexacte))
100 FORMAT(10(E10.3,2x))
200 FORMAT(I6,2x, 6(E10.3,2x))
211 FORMAT (2(E12.6,2x))
  close (uprint)
  close (usave)
  close (upoints)
  Deallocate(uTKL,uTKeLe,uetaSSe,cTKL,cTKeLe,cetaSSe,eTKL,eTKeLe,eetaSSe,vTKL,vTKeLe,vetaSSe)
  deallocate(AireD,AireDsommet,U,Uold,C,Cold,E,Eold,V,Vold,Gb)
  Deallocate(NormL1U,NormL1C,NormL1E,NormL1V)
  Deallocate(AireK,NsigK,NsigeKe,Msig,Msige)
  call cpu_time(finish)
  print*,'finish',start,finish,finish-start
  PRINT '("Fin du travail Keller-Segel en ",f10.3," secondes.")',finish-start
END PROGRAM GLIOMADDFV