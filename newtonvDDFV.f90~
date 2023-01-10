!                ******************************
!                **  SUBROUTINE NewtonuDDFV  **
!                ******************************
! *****************************************************************************
!     *-!!!!!!!!!@@@@@!!@@@@@!!@@@@@!!@@@@@!!!!!!!!!!!-*
!     *-!!!!!!!!!!!!!@!!@!!!@!!!!!!@!!@!!!@!!!!!!!!!!!-*
!     *-!!!!!!!!!@@@@@!!@!!!@!!@@@@@!!@!!!@!!!!!!!!!!!-*
!     *-!!!!!!!!!@!!!!!!@!!!@!!@@!!!!!@!!!@!!!!!!!!!!!-*
!     *-!!!!!!!!!@@@@@!!@@@@@!!@@@@@!!@@@@@!!!!!!!!!!!-*
!------------------------------------------------------------
!  Dans ce sous-programme on discretise l'equation en U par un scheme  implicite 
!       pour le scheme DDFV auto 
!
! ******************************************************************************
! 
SUBROUTINE NewtoncDDFV(A,Cold,C,Em,Um,TKL,TKeLe,etaSSe,ndim,choixf,temps)
  !--------
  ! Modules
  !--------
  USE longr
  USE imprime
  USE parmmage
  USE intmatvec
  USE intgradc
  USE fsource

  IMPLICIT NONE

  !--------------------------
  ! Declaration des arguments
  !--------------------------
  TYPE(MatCreux)                                   :: A
  Integer, intent(in)                              :: ndim,choixf
  REAL(kind=long), dimension (ndim), intent(in)    :: Cold,Em,Um
  REAL(kind=long), dimension (Nseg), intent(in)    :: TKL,TKeLe,etaSSe
  REAL(kind=long), dimension (ndim), intent(out)   :: C
  REAL(kind=long), intent(in)                      :: temps

  !----------------------------------
  ! Declaration des variables locales
  !----------------------------------
  CHARACTER(len=6)      :: oldprf
  INTEGER               :: jt,is,js,ks,ii,iK,jL,iseg,kiter,jj,j
  REAL(kind=long)       :: CoefAjout, x1, y1 
  REAL(kind = long)     :: Adegenbord
  REAL(kind = long)     :: dVKLmoins, dVKLplus, FKL, XKjL
  REAL(kind = long)     :: dVKL,coef,coeta,coefe,Cbord,Cibord,Cjbord,RDVT
  REAL(kind = long), dimension(ndim)     :: Xk, Yk, X0, F
  LOGICAL               :: trouve
  !-------------------
  ! Debut du programme
  !-------------------
  oldprf = prefix
  prefix = 'NEWTON'

  !------
  ! Corps
  !------
  Xk = Cold
  print*,'max Xk initial',maxval(Xk)
  DO kiter = 1, kitermax
     !
     A%TMAT = 0.D0
     !!=============================
     !! calcul de F(Xk) et F'(Xk)=AU
     !!=============================
     !
     !--------------------
     !--1-- terme en temps
     !--------------------
     !
     DO is = 1, NsInt
        F(is) = AireDsommet(is)*(CFDT*( Xk(is)-Cold(is) )/dt - alpha2*Em(is) + (beta2+gamma2*Um(is))*Xk(is))
        !! derivee de F par rapport a Xk(is)
        call ajout(is,is, AireDSommet(is)*(CFDT/dt + beta2+gamma2*Um(is)), A )
     END DO
     !
     DO jt = 1, Nbt
        ii = jt + NsInt
        F(ii) = AireK(jt)*(CFDT*( Xk(ii)-Cold(ii) )/dt - alpha2*Em(ii) + (beta2+gamma2*Um(ii))*Xk(ii))
        !! derivee de F par rapport a Xk
        call ajout(ii,ii, AireK(jt)*(CFDT/dt + beta2+gamma2*Um(ii)), A )
     END DO
     !--------------------------------------------
     ! --2-- terme laplacien  -div(S(x) grad A(u))
     !--------------------------------------------      
     ! assemblage par segment
     Do iseg = 1,Nseg
        ii = (NtypSeg(iseg)) 
        Select case (ii) 
        case (0)   !! segment à l'intérieur  
           is = NuSeg(1,iseg); js = NuSeg(2,iseg)
           iK = NumTVoisSeg(1,iseg) + NsInt; jL = NumTVoisSeg(2,iseg) + NsInt
           !
           IF ( is <= NsInt) THEN ! sommet is à l'intérieur
              IF (js <= NsInt ) THEN ! sommet js à l'intérieur
                 coef = TKL(iseg) ; coeta = etaSSe(iseg)
                 RDVT = (Xk(is) + Xk(js) + Xk(iK) + Xk(jL))/4
                 !-------------------------------
                 ! 1. contribution du triangle iK
                 !-------------------------------
                 FKL = coef*(ln(Xk(iK)) - ln(Xk(jL))) + &
                      & coeta * ( ln(Xk(is)) - ln(Xk(js)) )
                 F(iK) = F(iK) + RDVT * FKL
                 !! derivee de F(iK) par rapport a Xk(iK)
                 CoefAjout = RDVT*coef*Derivln(Xk(iK)) + FKL/4
                 call ajout(iK,iK,CoefAjout, A )
                 !! derivee de F(iK) par rapport a Xk(jL)
                 CoefAjout = -RDVT*coef*Derivln(Xk(jL)) + FKL/4
                 call ajout(iK,jL,CoefAjout, A )
                 !! derivee de F(iK) par rapport a Xk(is)
                 CoefAjout = RDVT*coeta*Derivln(Xk(is)) + FKL/4
                 call ajout(iK,is, CoefAjout, A)
                 !! derivee de F(iK) par rapport a Xk(js)
                 CoefAjout = -RDVT*coeta*Derivln(Xk(js)) + FKL/4
                 call ajout(iK,js, CoefAjout, A)
                 !-------------------------------
                 ! 2. contribution du triangle jL
                 !-------------------------------
                 FKL = coef*(ln(Xk(jL)) - ln(Xk(iK))) + &
                      & coeta * ( ln(Xk(js)) - ln(Xk(is)) )
                 F(jL) = F(jL) + RDVT * FKL
                 !! derivee de F(jL) par rapport a Xk(jL)
                 CoefAjout = RDVT*coef*Derivln(Xk(jL)) + FKL/4
                 call ajout(jL,jL,CoefAjout, A )
                 !! derivee de F(jL) par rapport a Xk(iK)
                 CoefAjout = -RDVT*coef*Derivln(Xk(iK)) + FKL/4
                 call ajout(jL,iK,CoefAjout, A )
                 !! derivee de F(jL) par rapport a Xk(is)
                 CoefAjout = -RDVT*coeta*Derivln(Xk(is)) + FKL/4
                 call ajout(jL,is, CoefAjout, A)
                 !! derivee de F(jL) par rapport a Xk(js)
                 CoefAjout = RDVT*coeta*Derivln(Xk(js)) + FKL/4
                 call ajout(jL,js, CoefAjout, A)
                 !-----------------------------
                 ! 3. contribution du sommet is
                 !-----------------------------
                 coef = TKeLe(iseg)
                 FKL = coef*(ln(Xk(is)) - ln(Xk(js))) + &
                      & coeta * ( ln(Xk(iK)) - ln(Xk(jL)) )
                 F(is) =  F(is) + RDVT * FKL
                 !! derivee de F(is) par rapport a Xk(is)
                 CoefAjout = RDVT*coef*Derivln(Xk(is)) + FKL/4
                 call ajout(is,is,CoefAjout, A )
                 !! derivee de F(is) par rapport a Xk(js)
                 CoefAjout = -RDVT*coef*Derivln(Xk(js)) + FKL/4
                 call ajout(is,js,CoefAjout, A )
                 !! derivee de F(is) par rapport a Xk(iK)
                 CoefAjout = RDVT*coeta*Derivln(Xk(iK)) + FKL/4
                 call ajout(is,iK, CoefAjout, A)
                 !! derivee de F(is) par rapport a Xk(jL)
                 CoefAjout = -RDVT*coeta*Derivln(Xk(jL)) + FKL/4
                 call ajout(is,jL, CoefAjout, A)
                 !-----------------------------
                 ! 4. contribution du sommet js
                 !-----------------------------
                 FKL = coef*(ln(Xk(js)) - ln(Xk(is))) + &
                      & coeta * ( ln(Xk(jL)) - ln(Xk(iK)) )
                 F(js) = F(js) + RDVT * FKL
                 !! derivee de F(js) par rapport a Xk(js)
                 CoefAjout = RDVT*coef*Derivln(Xk(js)) + FKL/4
                 call ajout(js,js,CoefAjout, A )
                 !! derivee de F(js) par rapport a Xk(is)
                 CoefAjout = -RDVT*coef*Derivln(Xk(is)) + FKL/4
                 call ajout(js,is,CoefAjout, A )
                 !! derivee de F(js) par rapport a Xk(iK)
                 CoefAjout = -RDVT*coeta*Derivln(Xk(iK)) + FKL/4
                 call ajout(js,iK, CoefAjout, A)
                 !! derivee de F(js) par rapport a Xk(jL)
                 CoefAjout = RDVT*coeta*Derivln(Xk(jL)) + FKL/4
                 call ajout(js,jL, CoefAjout, A)
                 !
              ELSE
                 Select case ( Ntyps(js) ) 
                 case (dirichlet)
                    Cbord = Gb(js)
                    coef = TKL(iseg) ; coeta = etaSSe(iseg)
                    RDVT = (Xk(is) + Cbord + Xk(iK) + Xk(jL))/4
                    !-------------------------------
                    ! 1. contribution du triangle iK
                    !-------------------------------
                    FKL = coef*(ln(Xk(iK)) - ln(Xk(jL))) + &
                         & coeta * ( ln(Xk(is)) - ln(Cbord) )
                    F(iK) = F(iK) + RDVT * FKL
                    !! derivee de F(iK) par rapport a Xk(iK)
                    CoefAjout = RDVT*coef*Derivln(Xk(iK)) + FKL/4
                    call ajout(iK,iK,CoefAjout, A )
                    !! derivee de F(iK) par rapport a Xk(jL)
                    CoefAjout = -RDVT*coef*Derivln(Xk(jL)) + FKL/4
                    call ajout(iK,jL,CoefAjout, A )
                    !! derivee de F(iK) par rapport a Xk(is)
                    CoefAjout = RDVT*coeta*Derivln(Xk(is)) + FKL/4
                    call ajout(iK,is, CoefAjout, A)
                    !-------------------------------
                    ! 2. contribution du triangle jL
                    !-------------------------------
                    FKL = coef*(ln(Xk(jL)) - ln(Xk(iK))) + &
                         &coeta * ( ln(Cbord) - ln(Xk(is)) )
                    F(jL) = F(jL) + RDVT * FKL
                    !! derivee de F(jL) par rapport a Xk(jL)
                    CoefAjout = RDVT*coef*Derivln(Xk(jL)) + FKL/4
                    call ajout(jL,jL,CoefAjout, A )
                    !! derivee de F(jL) par rapport a Xk(iK)
                    CoefAjout = -RDVT*coef*Derivln(Xk(iK)) + FKL/4
                    call ajout(jL,iK,CoefAjout, A )
                    !! derivee de F(jL) par rapport a Xk(is)
                    CoefAjout = -RDVT*coeta*Derivln(Xk(is)) + FKL/4
                    call ajout(jL,is, CoefAjout, A)
                    !-----------------------------
                    ! 3. contribution du sommet is
                    !-----------------------------
                    coef = TKeLe(iseg)
                    FKL = coef*(ln(Xk(is)) - ln(Cbord)) + &
                         & coeta * ( ln(Xk(iK)) - ln(Xk(jL)) )
                    F(is) = F(is) + RDVT * FKL
                    !! derivee de F(is) par rapport a Xk(is)
                    CoefAjout = RDVT*coef*Derivln(Xk(is)) + FKL/4
                    call ajout(is,is,CoefAjout, A )
                    !! derivee de F(is) par rapport a Xk(iK)
                    CoefAjout = RDVT*coeta*Derivln(Xk(iK)) + FKL/4
                    call ajout(is,iK, CoefAjout, A)
                    !! derivee de F(is) par rapport a Xk(jL)
                    CoefAjout = -RDVT*coeta*Derivln(Xk(jL)) + FKL/4
                    call ajout(is,jL, CoefAjout, A)
                    !
                 case(neumann)
                    ! Déjà faite (is <= NsInt est la même chose que is <= Nbs ce qui est correcte)
                 End Select
              END IF
              !
           ELSE
              !
              Select case ( Ntyps(is) )
              case (dirichlet)
                 IF (js <= NsInt) THEN
                    Cbord= Gb(is)
                    coef = TKL(iseg) ; coeta = etaSSe(iseg)
                    RDVT = (Cbord + Xk(js) + Xk(iK) + Xk(jL))/4
                    !-------------------------------
                    ! 1. contribution du triangle iK
                    !-------------------------------
                    FKL = coef*(ln(Xk(iK)) - ln(Xk(jL))) + &
                         & coeta * ( ln(Cbord) - ln(Xk(js)) )
                    F(iK) = F(iK) + RDVT * FKL
                    !! derivee de F(iK) par rapport a Xk(iK)
                    CoefAjout = RDVT*coef*Derivln(Xk(iK)) + FKL/4
                    call ajout(iK,iK,CoefAjout, A )
                    !! derivee de F(iK) par rapport a Xk(jL)
                    CoefAjout = -RDVT*coef*Derivln(Xk(jL)) + FKL/4
                    call ajout(iK,jL,CoefAjout, A )
                    !! derivee de F(iK) par rapport a Xk(js)
                    CoefAjout = -RDVT*coeta*Derivln(Xk(js)) + FKL/4
                    call ajout(iK,js, CoefAjout, A)
                    !-------------------------------
                    ! 2. contribution du triangle jL
                    !-------------------------------
                    FKL = coef*(ln(Xk(jL)) - ln(Xk(iK))) + &
                         &coeta * ( ln(Xk(js)) - ln(Cbord) )
                    F(jL) = F(jL) + RDVT * FKL
                    !! derivee de F(jL) par rapport a Xk(jL)
                    CoefAjout = RDVT*coef*Derivln(Xk(jL)) + FKL/4
                    call ajout(jL,jL,CoefAjout, A )
                    !! derivee de F(jL) par rapport a Xk(iK)
                    CoefAjout = -RDVT*coef*Derivln(Xk(iK)) + FKL/4
                    call ajout(jL,iK,CoefAjout, A )
                    !! derivee de F(jL) par rapport a Xk(js)
                    CoefAjout = RDVT*coeta*Derivln(Xk(js)) + FKL/4
                    call ajout(jL,js, CoefAjout, A)
                    !-----------------------------
                    ! 3. contribution du sommet js
                    !-----------------------------
                    coef = TKeLe(iseg)
                    FKL = coef*(ln(Xk(js)) - ln(Cbord)) + &
                         & coeta * ( ln(Xk(jL)) - ln(Xk(iK)) )
                    F(js) = F(js) + RDVT * FKL
                    !! derivee de F(js) par rapport a Xk(js)
                    CoefAjout = RDVT*coef*Derivln(Xk(js)) + FKL/4
                    call ajout(js,js,CoefAjout, A )
                    !! derivee de F(js) par rapport a Xk(iK)
                    CoefAjout = -RDVT*coeta*Derivln(Xk(iK)) + FKL/4
                    call ajout(js,iK, CoefAjout, A)
                    !! derivee de F(js) par rapport a Xk(jL)
                    CoefAjout = RDVT*coeta*Derivln(Xk(jL)) + FKL/4
                    call ajout(js,jL, CoefAjout, A)
                 ELSE ! js se situe encore au bord (un triangle entier sur le bord)
                    Cibord= Gb(is) ; Cjbord= Gb(js)
                    coef = TKL(iseg) ; coeta = etaSSe(iseg)
                    RDVT = (Cibord + Cjbord + Xk(iK) + Xk(jL))/4
                    !-------------------------------
                    ! 1. contribution du triangle iK
                    !-------------------------------
                    FKL = coef*(ln(Xk(iK)) - ln(Xk(jL))) + &
                         & coeta * ( ln(Cibord) - ln(Cjbord) )
                    F(iK) = F(iK) + RDVT * FKL
                    !! derivee de F(iK) par rapport a Xk(iK)
                    CoefAjout = RDVT*coef*Derivln(Xk(iK)) + FKL/4
                    call ajout(iK,iK,CoefAjout, A )
                    !! derivee de F(iK) par rapport a Xk(jL)
                    CoefAjout = -RDVT*coef*Derivln(Xk(jL)) + FKL/4
                    call ajout(iK,jL,CoefAjout, A )
                    !-------------------------------
                    ! 2. contribution du triangle jL
                    !-------------------------------
                    FKL = coef*(ln(Xk(jL)) - ln(Xk(iK))) + &
                         & coeta * ( ln(Cjbord) - ln(Cibord) )
                    F(jL) = F(jL) + RDVT * FKL
                    !! derivee de F(jL) par rapport a Xk(jL)
                    CoefAjout = RDVT*coef*Derivln(Xk(jL)) + FKL/4
                    call ajout(jL,jL,CoefAjout, A )
                    !! derivee de F(jL) par rapport a Xk(iK)
                    CoefAjout = -RDVT*coef*Derivln(Xk(iK)) + FKL/4
                    call ajout(jL,iK,CoefAjout, A )
                 END IF
              case(neumann)
                 ! On fait rien
              End Select
           END IF
           !!
        case (dirichlet) !! segment au bord  (triangle dégénéré)
           !!
           is = NuSeg(1,iseg); js = NuSeg(2,iseg)
           iK = NumTVoisSeg(1,iseg) + NsInt; !jL = numero milieu iseg
           Cibord= Gb(is) ; Cjbord= Gb(js) ; Cbord= Gb(iseg + Nbs)
           coef = TKL(iseg) ; coeta = etaSSe(iseg)
           RDVT = (Cibord + Cjbord + Xk(iK) + Cbord)/4
           !-------------------------------
           ! 1. contribution du triangle iK
           !-------------------------------
           FKL = coef*(ln(Xk(iK)) - ln(Cbord)) + &
                &coeta * ( ln(Cibord) - ln(Cjbord) )
           F(iK) = F(iK) + RDVT * FKL
           !! derivee de F(iK) par rapport a Xk(iK)
           CoefAjout = RDVT*coef*Derivln(Xk(iK)) + FKL/4
           call ajout(iK,iK,CoefAjout, A )
           !!
        case (neumann)
           !!
           is = NuSeg(1,iseg); js = NuSeg(2,iseg)
           iK = NumTVoisSeg(1,iseg) + NsInt
           !-------------------------------
           ! 1. contribution du triangle iK
           !-------------------------------
           ! On fait rien
           !-----------------------------
           ! 2. contribution du sommet is
           !-----------------------------
           RDVT = (3*(Xk(is) + Xk(js))/2 + Xk(iK))/4
           coef = TKL(iseg) ; coeta = etaSSe(iseg) ; coefe = TKeLe(iseg)
           coefe = coefe - (coeta**2)/coef
           FKL = coefe*(ln(Xk(is)) - ln(Xk(js)))
           F(is) = F(is) + RDVT * FKL
           !! derivee de F(is) par rapport a Xk(iK)
           CoefAjout = FKL/4
           call ajout(is,iK,CoefAjout, A )
           !! derivee de F(is) par rapport a Xk(is)
           CoefAjout = RDVT*coefe*Derivln(Xk(is))+ 3*FKL/8
           call ajout(is,is,CoefAjout, A )
           !! derivee de F(is) par rapport a Xk(js)
           CoefAjout = -RDVT*coefe*Derivln(Xk(js))+ 3*FKL/8
           call ajout(is,js,CoefAjout, A )
           !-----------------------------
           ! 3. contribution du sommet js
           !-----------------------------
           FKL = coefe*(ln(Xk(js)) - ln(Xk(is)))
           F(js) = F(js) + RDVT * FKL
           !! derivee de F(js) par rapport a Xk(iK)
           CoefAjout = FKL/4
           call ajout(js,iK,CoefAjout, A )
           !! derivee de F(js) par rapport a Xk(js)
           CoefAjout = RDVT*coefe*Derivln(Xk(js))+ 3*FKL/8
           call ajout(js,js,CoefAjout, A )
           !! derivee de F(js) par rapport a Xk(is)
           CoefAjout = -RDVT*coefe*Derivln(Xk(is))+ 3*FKL/8
           call ajout(js,is,CoefAjout, A )
        END Select
     END Do
     !!
     !!============================
     !!============================
     !!print*,'debut sys lin'
     !Yk = -F/A
     !print*,'max F avant bigradient', maxval(F), 'min F', minval(F)
     X0 = 0.D0
     !print*,'A de Newton=',A%Tmat
     !Yk = gradconj(-F , A)
     Yk = bigradient(A,-F, X0,Tolerencegradient)
     CALL prvari(uprint,'Matrice A = ',Ndim )
     !write(*,*)'Yk=',Yk
     !write(*,*)'F(1)=',F(1)
     print*,'kiter = ', kiter,'erreur NEWTONV sqrt(sum(Yk*Yk))',sqrt(sum(Yk*Yk))
     If (sqrt(dot_product(Yk,Yk)) <TolerenceNewton) exit
     Xk=Xk+Yk
     !print*,'max Xk', maxval(Xk), 'min Xk', minval(Xk)
     !print*,'max F', maxval(F), 'min F', minval(F)
     !stop
     !! impression de A etde F
     !     DO i = 1, size(A%IndPL)-1
     !           print*,'iseg = ', i, 'sum line = ', SUM(A%Tmat(A%IndPL(i): A%IndPL(i+1)-1))
     !    END DO
  END Do
  C=Xk+Yk
  print*,'fin newtonv'
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


END SUBROUTINE NewtonvDDFV






