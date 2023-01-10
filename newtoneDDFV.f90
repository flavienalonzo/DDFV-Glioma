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
SUBROUTINE NewtoneDDFV(A,Eold,E,Um,Vm,TKL,TKeLe,etaSSe,ndim,choixf,temps)
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
  REAL(kind=long), dimension (ndim), intent(in)    :: Eold,Um,Vm
  REAL(kind=long), dimension (ndim), intent(out)   :: E
  REAL(kind=long), dimension (Nseg), intent(in)    :: TKL,TKeLe,etaSSe
  REAL(kind=long), intent(in)                      :: temps

  !----------------------------------
  ! Declaration des variables locales
  !----------------------------------
  CHARACTER(len=6)      :: oldprf
  INTEGER               :: jt,is,js,ks,ii,iK,jL,iseg,kiter,jj,j
  REAL(kind=long)       :: CoefAjout, x1, y1 
  REAL(kind = long)     :: Adegenbord,RDVT,XkjL,VjL,FKL
  REAL(kind = long)     :: dVKLmoins, dVKLplus
  REAL(kind = long)     :: dVKL,coef,coeta,coefe,Ubord,Uibord,Ujbord
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
  Xk = Eold
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
        F(is) = AireDsommet(is)*(CFDT*( Xk(is)-Eold(is) )/dt - rho3*fe(Xk(is),Um(is)) + beta3*Xk(is) ) 
        !! derivee de F par rapport a Xk(is)
        call ajout(is,is, AireDSommet(is)*(CFDT/dt - rho3*deriveefe(Xk(is),Um(is)) + beta3), A )
     END DO
     !
     DO jt = 1, Nbt
        ii = jt + NsInt
        F(ii) = CFDT*AireK(jt)*(( Xk(ii)-Eold(ii) )/dt - rho3*fe(Xk(ii),Um(ii)) + beta3*Xk(ii) )  
        !! derivee de F par rapport a Xk(ii)
        call ajout(ii,ii, AireK(jt)*(CFDT/dt - rho3*deriveefe(Xk(ii),Um(ii)) + beta3), A )
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
                 !-------------------------------
                 ! 1. contribution du triangle iK
                 !-------------------------------
                 dVKL = coeta * ( xi(Xk(is)) - xi(Xk(js)) )
                 dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                 F(iK) = F(iK) + coef*(Adegen(Xk(iK)) - Adegen(Xk(jL))) + &
                      & dVKLplus*( rmCroit(Xk(iK)) + rmDecroit(Xk(jL)) ) + &
                      & dVKLmoins*( rmCroit(Xk(jL)) + rmDecroit(Xk(iK)))
                 !! derivee de F(iK) par rapport a Xk(iK)
                 CoefAjout = coef*DerivAdeg(Xk(iK))+ dVKLplus*DerrmCroit(Xk(iK)) + dVKLmoins*DerrmDecroit(Xk(iK))
                 call ajout(iK,iK,CoefAjout, A )
                 !! derivee de F(iK) par rapport a Xk(jL)
                 CoefAjout = -coef*DerivAdeg(Xk(jL))+ dVKLplus*DerrmDeCroit(Xk(jL)) + dVKLmoins*Derrmcroit(Xk(jL))
                 call ajout(iK,jL,CoefAjout, A )
                 !! derivee de F(iK) par rapport a Xk(is)
                 CoefAjout = transfer((dVKL>0),1)*coeta* rm(Xk(is))*( rmCroit(Xk(iK)) + rmDecroit(Xk(jL)) ) &
                      & + transfer((dVKL<0),1)*coeta* rm(Xk(is))*( rmCroit(Xk(jL)) + rmDecroit(Xk(iK)) )
                 call ajout(iK,is, CoefAjout, A)
                 !! derivee de F(iK) par rapport a Xk(js)
                 CoefAjout = -transfer((dVKL>0),1)*coeta* rm(Xk(js))*( rmCroit(Xk(iK)) + rmDecroit(Xk(jL)) ) &
                      & - transfer((dVKL<0),1)*coeta* rm(Xk(js))*( rmCroit(Xk(jL)) + rmDecroit(Xk(iK)) )
                 call ajout(iK,js, CoefAjout, A)
                 !-------------------------------
                 ! 2. contribution du triangle jL
                 !-------------------------------
                 dVKL = coeta * ( xi(Xk(js)) - xi(Xk(is)) )
                 dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                 F(jL) = F(jL) + coef*(Adegen(Xk(jL)) - Adegen(Xk(iK))) + &
                      & dVKLplus*( rmCroit(Xk(jL)) + rmDecroit(Xk(iK)) ) + &
                      & dVKLmoins*( rmCroit(Xk(iK)) + rmDecroit(Xk(jL)))
                 !! derivee de F(jL) par rapport a Xk(iK)
                 CoefAjout = -coef*DerivAdeg(Xk(iK))+ dVKLplus*DerrmDeCroit(Xk(iK)) + dVKLmoins*Derrmcroit(Xk(iK))
                 call ajout(jL,iK,CoefAjout, A )
                 !! derivee de F(jL) par rapport a Xk(jL)
                 CoefAjout = coef*DerivAdeg(Xk(jL))+ dVKLplus*DerrmCroit(Xk(jL)) + dVKLmoins*DerrmDecroit(Xk(jL))
                 call ajout(jL,jL,CoefAjout, A )
                 !! derivee de F(jL) par rapport a Xk(is)
                 CoefAjout = -transfer((dVKL>0),1)*coeta* rm(Xk(is))*( rmCroit(Xk(jL)) + rmDecroit(Xk(iK)) ) &
                      & -transfer((dVKL<0),1)*coeta* rm(Xk(is))*( rmCroit(Xk(iK)) + rmDecroit(Xk(jL)) )
                 call ajout(jL,is, CoefAjout, A)
                 !! derivee de F(jL) par rapport a Xk(js)
                 CoefAjout = transfer((dVKL>0),1)*coeta* rm(Xk(js))*( rmCroit(Xk(jL)) + rmDecroit(Xk(iK)) ) &
                      & + transfer((dVKL<0),1)*coeta* rm(Xk(js))*( rmCroit(Xk(iK)) + rmDecroit(Xk(jL)) )
                 call ajout(jL,js, CoefAjout, A)
                 !-----------------------------
                 ! 3. contribution du sommet is
                 !-----------------------------
                 coef = TKeLe(iseg)
                 dVKL = coeta * ( xi(Xk(iK)) - xi(Xk(jL)) )
                 dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                 F(is) = F(is) + coef*(Adegen(Xk(is)) - Adegen(Xk(js))) + &
                      & dVKLplus*( rmCroit(Xk(is)) + rmDecroit(Xk(js)) ) + &
                      & dVKLmoins*( rmCroit(Xk(js)) + rmDecroit(Xk(is)))
                 !! derivee de F(is) par rapport a Xk(is)
                 CoefAjout = coef*DerivAdeg(Xk(is))+ dVKLplus*DerrmCroit(Xk(is)) + dVKLmoins*DerrmDecroit(Xk(is))
                 call ajout(is,is,CoefAjout, A )
                 !! derivee de F(is) par rapport a Xk(js)
                 CoefAjout = -coef*DerivAdeg(Xk(js))+ dVKLplus*DerrmDeCroit(Xk(js)) + dVKLmoins*Derrmcroit(Xk(js))
                 call ajout(is,js,CoefAjout, A )
                 !! derivee de F(is) par rapport a Xk(iK)
                 CoefAjout = transfer((dVKL>0),1)*coeta* rm(Xk(iK))*( rmCroit(Xk(is)) + rmDecroit(Xk(js)) ) &
                      & + transfer((dVKL<0),1)*coeta* rm(Xk(iK))*( rmCroit(Xk(js)) + rmDecroit(Xk(is)) )
                 call ajout(is,iK, CoefAjout, A)
                 !! derivee de F(is) par rapport a Xk(jL)
                 CoefAjout = -transfer((dVKL>0),1)*coeta* rm(Xk(jL))*( rmCroit(Xk(is)) + rmDecroit(Xk(js)) ) &
                      & - transfer((dVKL<0),1)*coeta* rm(Xk(jL))*( rmCroit(Xk(js)) + rmDecroit(Xk(is)) )
                 call ajout(is,jL, CoefAjout, A)
                 !-----------------------------
                 ! 4. contribution du sommet js
                 !-----------------------------
                 dVKL = coeta * ( xi(Xk(jL)) - xi(Xk(iK)) )
                 dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                 F(js) = F(js) + coef*(Adegen(Xk(js)) - Adegen(Xk(is))) + &
                      & dVKLplus*( rmCroit(Xk(js)) + rmDecroit(Xk(is)) ) + &
                      & dVKLmoins*( rmCroit(Xk(is)) + rmDecroit(Xk(js)))
                 !! derivee de F(js) par rapport a Xk(is)
                 CoefAjout = -coef*DerivAdeg(Xk(is))+ dVKLplus*DerrmDeCroit(Xk(is)) + dVKLmoins*Derrmcroit(Xk(is))
                 call ajout(js,is,CoefAjout, A )
                 !! derivee de F(js) par rapport a Xk(js)
                 CoefAjout = coef*DerivAdeg(Xk(js))+ dVKLplus*DerrmCroit(Xk(js)) + dVKLmoins*DerrmDecroit(Xk(js))
                 call ajout(js,js,CoefAjout, A )
                 !! derivee de F(js) par rapport a Xk(iK)
                 CoefAjout = -transfer((dVKL>0),1)*coeta* rm(Xk(iK))*( rmCroit(Xk(js)) + rmDecroit(Xk(is)) ) &
                      & - transfer((dVKL<0),1)*coeta* rm(Xk(iK))*( rmCroit(Xk(is)) + rmDecroit(Xk(js)) )
                 call ajout(js,iK, CoefAjout, A)
                 !! derivee de F(js) par rapport a Xk(jL)
                 CoefAjout = transfer((dVKL>0),1)*coeta* rm(Xk(jL))*( rmCroit(Xk(js)) + rmDecroit(Xk(is)) ) &
                      & + transfer((dVKL<0),1)*coeta* rm(Xk(jL))*( rmCroit(Xk(is)) + rmDecroit(Xk(js)) )
                 call ajout(js,jL, CoefAjout, A)
                 !
              ELSE
                 Select case ( Ntyps(js) ) 
                 case (dirichlet)
                    Ubord = Gb(js)
                    coef = TKL(iseg) ; coeta = etaSSe(iseg)
                    !-------------------------------
                    ! 1. contribution du triangle iK
                    !-------------------------------
                    dVKL = coeta * ( xi(Xk(is)) - xi(Ubord) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    F(iK) = F(iK) + coef*(Adegen(Xk(iK)) - Adegen(Xk(jL))) + &
                         & dVKLplus*( rmCroit(Xk(iK)) + rmDecroit(Xk(jL)) ) + &
                         & dVKLmoins*( rmCroit(Xk(jL)) + rmDecroit(Xk(iK)))
                    !! derivee de F(iK) par rapport a Xk(iK)
                    CoefAjout = coef*DerivAdeg(Xk(iK))+ dVKLplus*DerrmCroit(Xk(iK)) + dVKLmoins*DerrmDecroit(Xk(iK))
                    call ajout(iK,iK,CoefAjout, A )
                    !! derivee de F(iK) par rapport a Xk(jL)
                    CoefAjout = -coef*DerivAdeg(Xk(jL))+ dVKLplus*DerrmDeCroit(Xk(jL)) + dVKLmoins*Derrmcroit(Xk(jL))
                    call ajout(iK,jL,CoefAjout, A )
                    !! derivee de F(iK) par rapport a Xk(is)
                    CoefAjout = transfer((dVKL>0),1)*coeta* rm(Xk(is))*( rmCroit(Xk(iK)) + rmDecroit(Xk(jL)) ) &
                         & + transfer((dVKL<0),1)*coeta* rm(Xk(is))*( rmCroit(Xk(jL)) + rmDecroit(Xk(iK)) )
                    call ajout(iK,is, CoefAjout, A)
                    !-------------------------------
                    ! 2. contribution du triangle jL
                    !-------------------------------
                    dVKL = coeta * ( xi(Ubord) - xi(Xk(is)) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    F(jL) = F(jL) + coef*(Adegen(Xk(jL)) - Adegen(Xk(iK))) + &
                         & dVKLplus*( rmCroit(Xk(jL)) + rmDecroit(Xk(iK)) ) + &
                         & dVKLmoins*( rmCroit(Xk(iK)) + rmDecroit(Xk(jL)))
                    !! derivee de F(jL) par rapport a Xk(iK)
                    CoefAjout = -coef*DerivAdeg(Xk(iK))+ dVKLplus*DerrmDeCroit(Xk(iK)) + dVKLmoins*Derrmcroit(Xk(iK))
                    call ajout(jL,iK,CoefAjout, A )
                    !! derivee de F(jL) par rapport a Xk(jL)
                    CoefAjout = coef*DerivAdeg(Xk(jL))+ dVKLplus*DerrmCroit(Xk(jL)) + dVKLmoins*DerrmDecroit(Xk(jL))
                    call ajout(jL,jL,CoefAjout, A )
                    !! derivee de F(jL) par rapport a Xk(is)
                    CoefAjout = -transfer((dVKL>0),1)*coeta* rm(Xk(is))*( rmCroit(Xk(jL)) + rmDecroit(Xk(iK)) ) &
                         & - transfer((dVKL<0),1)*coeta* rm(Xk(is))*( rmCroit(Xk(iK)) + rmDecroit(Xk(jL)) )
                    call ajout(jL,is, CoefAjout, A)
                    !-----------------------------
                    ! 3. contribution du sommet is
                    !-----------------------------
                    coef = TKeLe(iseg)
                    dVKL = coeta * ( xi(Xk(iK)) - xi(Xk(jL)) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    F(is) = F(is) + coef*(Adegen(Xk(is)) - Adegen(Ubord)) + &
                         & dVKLplus*( rmCroit(Xk(is)) + rmDecroit(Ubord) ) + &
                         & dVKLmoins*( rmCroit(Ubord) + rmDecroit(Xk(is)))
                    !! derivee de F(is) par rapport a Xk(is)
                    CoefAjout = coef*DerivAdeg(Xk(is))+ dVKLplus*DerrmCroit(Xk(is)) + dVKLmoins*DerrmDecroit(Xk(is))
                    call ajout(is,is,CoefAjout, A )
                    !! derivee de F(is) par rapport a Xk(iK)
                    CoefAjout = transfer((dVKL>0),1)*coeta* rm(Xk(iK))*( rmCroit(Xk(is)) + rmDecroit(Ubord) ) &
                         & + transfer((dVKL<0),1)*coeta* rm(Xk(iK))*( rmCroit(Ubord) + rmDecroit(Xk(is)) )
                    call ajout(is,iK, CoefAjout, A)
                    !! derivee de F(is) par rapport a Xk(jL)
                    CoefAjout = -transfer((dVKL>0),1)*coeta* rm(Xk(jL))*( rmCroit(Xk(is)) + rmDecroit(Ubord) ) &
                         & - transfer((dVKL<0),1)*coeta* rm(Xk(jL))*( rmCroit(Ubord) + rmDecroit(Xk(is)) )
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
                    Ubord= Gb(is)
                    coef = TKL(iseg) ; coeta = etaSSe(iseg)
                    !-------------------------------
                    ! 1. contribution du triangle iK
                    !-------------------------------
                    dVKL = coeta * ( xi(Ubord) - xi(Xk(js)) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    F(iK) = F(iK) + coef*(Adegen(Xk(iK)) - Adegen(Xk(jL))) + &
                         & dVKLplus*( rmCroit(Xk(iK)) + rmDecroit(Xk(jL)) ) + &
                         & dVKLmoins*( rmCroit(Xk(jL)) + rmDecroit(Xk(iK)))
                    !! derivee de F(iK) par rapport a Xk(iK)
                    CoefAjout = coef*DerivAdeg(Xk(iK))+ dVKLplus*DerrmCroit(Xk(iK)) + dVKLmoins*DerrmDecroit(Xk(iK))
                    call ajout(iK,iK,CoefAjout, A )
                    !! derivee de F(iK) par rapport a Xk(jL)
                    CoefAjout = -coef*DerivAdeg(Xk(jL))+ dVKLplus*DerrmDeCroit(Xk(jL)) + dVKLmoins*Derrmcroit(Xk(jL))
                    call ajout(iK,jL,CoefAjout, A )
                    !! derivee de F(iK) par rapport a Xk(js)
                    CoefAjout = -transfer((dVKL>0),1)*coeta* rm(Xk(js))*( rmCroit(Xk(iK)) + rmDecroit(Xk(jL)) )&
                         & - transfer((dVKL<0),1)*coeta* rm(Xk(js))*( rmCroit(Xk(jL)) + rmDecroit(Xk(iK)) )
                    call ajout(iK,js, CoefAjout, A)
                    !-------------------------------
                    ! 2. contribution du triangle jL
                    !-------------------------------
                    dVKL = coeta * ( xi(Xk(js)) - xi(Ubord) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    F(jL) = F(jL) + coef*(Adegen(Xk(jL)) - Adegen(Xk(iK))) + &
                         & dVKLplus*( rmCroit(Xk(jL)) + rmDecroit(Xk(iK)) ) + &
                         & dVKLmoins*( rmCroit(Xk(iK)) + rmDecroit(Xk(jL)))
                    !! derivee de F(jL) par rapport a Xk(iK)
                    CoefAjout = -coef*DerivAdeg(Xk(iK))+ dVKLplus*DerrmDeCroit(Xk(iK)) + dVKLmoins*Derrmcroit(Xk(iK))
                    call ajout(jL,iK,CoefAjout, A )
                    !! derivee de F(jL) par rapport a Xk(jL)
                    CoefAjout = coef*DerivAdeg(Xk(jL))+ dVKLplus*DerrmCroit(Xk(jL)) + dVKLmoins*DerrmDecroit(Xk(jL))
                    call ajout(jL,jL,CoefAjout, A )
                    !! derivee de F(jL) par rapport a Xk(js)
                    CoefAjout = transfer((dVKL>0),1)*coeta* rm(Xk(js))*( rmCroit(Xk(jL)) + rmDecroit(Xk(iK)) ) &
                         & + transfer((dVKL<0),1)*coeta* rm(Xk(js))*( rmCroit(Xk(iK)) + rmDecroit(Xk(jL)) )
                    call ajout(jL,js, CoefAjout, A)
                    !-----------------------------
                    ! 3. contribution du sommet js
                    !-----------------------------
                    coef = TKeLe(iseg)
                    dVKL = coeta * ( xi(Xk(jL)) - xi(Xk(iK)) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    F(js) = F(js) + coef*(Adegen(Xk(js)) - Adegen(Ubord)) + &
                         & dVKLplus*( rmCroit(Xk(js)) + rmDecroit(Ubord) ) + &
                         & dVKLmoins*( rmCroit(Ubord) + rmDecroit(Xk(js)))
                    !! derivee de F(js) par rapport a Xk(js)
                    CoefAjout = coef*DerivAdeg(Xk(js))+ dVKLplus*DerrmCroit(Xk(js)) + dVKLmoins*DerrmDecroit(Xk(js))
                    call ajout(js,js,CoefAjout, A )
                    !! derivee de F(js) par rapport a Xk(iK)
                    CoefAjout = -transfer((dVKL>0),1)*coeta* rm(Xk(iK))*( rmCroit(Xk(js)) + rmDecroit(Ubord) ) &
                         & - transfer((dVKL<0),1)*coeta* rm(Xk(iK))*( rmCroit(Ubord) + rmDecroit(Xk(js)) )
                    call ajout(js,iK, CoefAjout, A)
                    !! derivee de F(js) par rapport a Xk(jL)
                    CoefAjout = transfer((dVKL>0),1)*coeta* rm(Xk(jL))*( rmCroit(Xk(js)) + rmDecroit(Ubord) ) &
                         & + transfer((dVKL<0),1)*coeta* rm(Xk(jL))*( rmCroit(Ubord) + rmDecroit(Xk(js)) )
                    call ajout(js,jL, CoefAjout, A)
                 ELSE ! js se situe encore au bord (un triangle entier sur le bord)
                    Uibord= Gb(is) ; Ujbord= Gb(js)
                    coef = TKL(iseg) ; coeta = etaSSe(iseg)
                    !-------------------------------
                    ! 1. contribution du triangle iK
                    !-------------------------------
                    dVKL = coeta * ( xi(Uibord) - xi(Ujbord) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    F(iK) = F(iK) + coef*(Adegen(Xk(iK)) - Adegen(Xk(jL))) + &
                         & dVKLplus*( rmCroit(Xk(iK)) + rmDecroit(Xk(jL)) ) + &
                         & dVKLmoins*( rmCroit(Xk(jL)) + rmDecroit(Xk(iK)))
                    !! derivee de F(iK) par rapport a Xk(iK)
                    CoefAjout = coef*DerivAdeg(Xk(iK))+ dVKLplus*DerrmCroit(Xk(iK)) + dVKLmoins*DerrmDecroit(Xk(iK))
                    call ajout(iK,iK,CoefAjout, A )
                    !! derivee de F(iK) par rapport a Xk(jL)
                    CoefAjout = -coef*DerivAdeg(Xk(jL))+ dVKLplus*DerrmDeCroit(Xk(jL)) + dVKLmoins*Derrmcroit(Xk(jL))
                    call ajout(iK,jL,CoefAjout, A )
                    !-------------------------------
                    ! 2. contribution du triangle jL
                    !-------------------------------
                    dVKL = coeta * ( xi(Ujbord) - xi(Uibord) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    F(jL) = F(jL) + coef*(Adegen(Xk(jL)) - Adegen(Xk(iK))) + &
                         & dVKLplus*( rmCroit(Xk(jL)) + rmDecroit(Xk(iK)) ) + &
                         & dVKLmoins*( rmCroit(Xk(iK)) + rmDecroit(Xk(jL)))
                    !! derivee de F(jL) par rapport a Xk(iK)
                    CoefAjout = -coef*DerivAdeg(Xk(iK))+ dVKLplus*DerrmDeCroit(Xk(iK)) + dVKLmoins*Derrmcroit(Xk(iK))
                    call ajout(jL,iK,CoefAjout, A )
                    !! derivee de F(jL) par rapport a Xk(jL)
                    CoefAjout = coef*DerivAdeg(Xk(jL))+ dVKLplus*DerrmCroit(Xk(jL)) + dVKLmoins*DerrmDecroit(Xk(jL))
                    call ajout(jL,jL,CoefAjout, A )
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
           Uibord= Gb(is) ; Ujbord= Gb(js) ; Ubord= Gb(iseg + Nbs)
           coef = TKL(iseg) ; coeta = etaSSe(iseg)
           !-------------------------------
           ! 1. contribution du triangle iK
           !-------------------------------
           dVKL = coeta * ( xi(Uibord) - xi(Ujbord) )
           dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
           F(iK) = F(iK) + coef*(Adegen(Xk(iK)) - Adegen(Ubord)) + &
                & dVKLplus*( rmCroit(Xk(iK)) + rmDecroit(Ubord) ) + &
                & dVKLmoins*( rmCroit(Ubord) + rmDecroit(Xk(iK)))
           !! derivee de F(iK) par rapport a Xk(iK)
           CoefAjout = coef*DerivAdeg(Xk(iK))+ dVKLplus*DerrmCroit(Xk(iK)) + dVKLmoins*DerrmDecroit(Xk(iK))
           call ajout(iK,iK,CoefAjout, A )
           !!
        case (neumann)
           !!
           is = NuSeg(1,iseg); js = NuSeg(2,iseg)
           iK = NumTVoisSeg(1,iseg) + NsInt
           coef = TKL(iseg) ; coeta = etaSSe(iseg) ; coefe = TKeLe(iseg)
           !-------------------------------
           ! 1. contribution du triangle iK
           !-------------------------------
           ! On fait rien
           !-----------------------------
           ! 2. contribution du sommet is
           !----------------------------- 
           dVKL = ((coeta**2)/coef) * ( xi(Xk(js)) - xi(Xk(is)) )
           dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
           F(is) = F(is) + coefe*(Adegen(Xk(is)) - Adegen(Xk(js))) + &
                & dVKLplus*( rmCroit(Xk(is)) + rmDecroit(Xk(js)) ) + &
                & dVKLmoins*( rmCroit(Xk(js)) + rmDecroit(Xk(is)))
           !! derivee de F(is) par rapport a Xk(is)
           CoefAjout = coefe*DerivAdeg(Xk(is))+ dVKLplus*DerrmCroit(Xk(is)) + dVKLmoins*DerrmDecroit(Xk(is)) &
                & - transfer((dVKL>0),1)*((coeta**2)/coef)* rm(Xk(is))*( rmCroit(Xk(is)) + rmDecroit(Xk(js)))&
                & - transfer((dVKL<0),1)*((coeta**2)/coef)* rm(Xk(is))*( rmCroit(Xk(js)) + rmDecroit(Xk(is)))
           call ajout(is,is,CoefAjout, A )
           !! derivee de F(is) par rapport a Xk(js)
           CoefAjout = -coefe*DerivAdeg(Xk(js))+ dVKLplus*DerrmDecroit(Xk(js)) + dVKLmoins*DerrmCroit(Xk(js))&
                & + transfer((dVKL>0),1)*((coeta**2)/coef)* rm(Xk(js))*( rmCroit(Xk(is)) + rmDecroit(Xk(js)))&
                & + transfer((dVKL<0),1)*((coeta**2)/coef)* rm(Xk(js))*( rmCroit(Xk(js)) + rmDecroit(Xk(is)))
           call ajout(is,js,CoefAjout, A )
           !-----------------------------
           ! 3. contribution du sommet js
           !-----------------------------
           dVKL = ((coeta**2)/coef) * ( xi(Xk(is)) - xi(Xk(js)) )
           dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
           F(js) = F(js) + coefe*(Adegen(Xk(js)) - Adegen(Xk(is))) + &
                & dVKLplus*( rmCroit(Xk(js)) + rmDecroit(Xk(is)) ) + &
                & dVKLmoins*( rmCroit(Xk(is)) + rmDecroit(Xk(js)))
           !! derivee de F(js) par rapport a Xk(is)
           CoefAjout = -coefe*DerivAdeg(Xk(is))+ dVKLplus*DerrmDecroit(Xk(is)) + dVKLmoins*DerrmCroit(Xk(is))&
                & + transfer((dVKL>0),1)*((coeta**2)/coef)* rm(Xk(is))*( rmCroit(Xk(js)) + rmDecroit(Xk(is)))&
                & + transfer((dVKL<0),1)*((coeta**2)/coef)* rm(Xk(is))*( rmCroit(Xk(is)) + rmDecroit(Xk(js)))
           call ajout(js,is,CoefAjout, A )
           !! derivee de F(js) par rapport a Xk(js)
           CoefAjout = coefe*DerivAdeg(Xk(js))+ dVKLplus*DerrmCroit(Xk(js)) + dVKLmoins*DerrmDecroit(Xk(js)) &
                & - transfer((dVKL>0),1)*((coeta**2)/coef)* rm(Xk(js))*( rmCroit(Xk(js)) + rmDecroit(Xk(is)))&
                & - transfer((dVKL<0),1)*((coeta**2)/coef)* rm(Xk(js))*( rmCroit(Xk(is)) + rmDecroit(Xk(js)))
           call ajout(js,js,CoefAjout, A )
        END Select
     END Do
     !!
     !------------------------------------------------
     ! --3-- terme de convection div(mu(u)S(x) grad v)
     !------------------------------------------------
     ! assemblage par segment
     Do iseg = 1, Nseg
        ii = (NtypSeg(iseg)) 
        Select case (ii) 
        case (0)   !! segment à l'intérieur  
           is = NuSeg(1,iseg); js = NuSeg(2,iseg)
           iK = NumTVoisSeg(1,iseg) + NsInt; jL = NumTVoisSeg(2,iseg) + NsInt
           !
           IF ( is <= NsInt) THEN ! sommet is à l'intérieur
              IF (js <= NsInt ) THEN ! sommet js à l'intérieur
                 coef = TKL(iseg) ; coeta = etaSSe(iseg)
                 RDVT = (V(is) + V(js) + V(iK) + V(jL))/4
                 !-------------------------------
                 ! 1. contribution du triangle iK
                 !-------------------------------
                 dVKL = coef*(ln(V(jL)) - ln(V(iK))) + &
                      & coeta*(ln(V(js)) - ln(V(is)))
                 dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                 FKL = dVKLplus*( MuCroit(Xk(iK)) + MuDecroit(Xk(jL)) ) + &
                      & dVKLmoins*( MuCroit(Xk(jL)) + MuDecroit(Xk(iK)) )
                 F(iK) = F(iK) + RDVT* FKL
                 !! derivee de F(iK) par rapport a Xk(iK)
                 CoefAjout = RDVT*(dVKLplus*DerivMuCroit(Xk(iK)) + dVKLmoins*DerivMuDecroit(Xk(iK)))
                 call ajout(iK,iK,CoefAjout, A )
                 !! derivee de F(iK) par rapport a Xk(jL)
                 CoefAjout = RDVT*(dVKLplus*DerivMuDecroit(Xk(jL)) + dVKLmoins*DerivMuCroit(Xk(jL)))
                 call ajout(iK,jL,CoefAjout, A )
                 !-------------------------------
                 ! 2. contribution du triangle jL
                 !-------------------------------
                 dVKL = coef*(ln(V(iK)) - ln(V(jL))) + &
                      & coeta * ( ln(V(is)) - ln(V(js)) )
                 dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                 FKL = dVKLplus*( MuCroit(Xk(jL)) + MuDecroit(Xk(iK)) ) + &
                      & dVKLmoins*( MuCroit(Xk(iK)) + MuDecroit(Xk(jL)) )
                 F(jL) = F(jL) + RDVT* FKL
                 !! derivee de F(jL) par rapport a Xk(jL)
                 CoefAjout = RDVT*(dVKLplus*DerivMuCroit(Xk(jL)) + dVKLmoins*DerivMuDecroit(Xk(jL)))
                 call ajout(jL,jL,CoefAjout, A )
                 !! derivee de F(jL) par rapport a Xk(iK)
                 CoefAjout = RDVT*(dVKLplus*DerivMuDecroit(Xk(iK)) + dVKLmoins*DerivMuCroit(Xk(iK)))
                 call ajout(jL,iK,CoefAjout, A )
                 !-----------------------------
                 ! 3. contribution du sommet is
                 !-----------------------------
                 coef = TKeLe(iseg)
                 dVKL = coef*(ln(V(js)) - ln(V(is))) + &
                      & coeta * ( ln(V(jL)) - ln(V(iK)) )
                 dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                 FKL = dVKLplus*( MuCroit(Xk(is)) + MuDecroit(Xk(js)) ) + &
                      & dVKLmoins*( MuCroit(Xk(js)) + MuDecroit(Xk(is)) )
                 F(is) = F(is) + RDVT* FKL
                 !! derivee de F(is) par rapport a Xk(is)
                 CoefAjout = RDVT*(dVKLplus*DerivMuCroit(Xk(is)) + dVKLmoins*DerivMuDecroit(Xk(is)))
                 call ajout(is,is,CoefAjout, A )
                 !! derivee de F(is) par rapport a Xk(js)
                 CoefAjout = RDVT*(dVKLplus*DerivMuDecroit(Xk(js)) + dVKLmoins*DerivMuCroit(Xk(js)))
                 call ajout(is,js,CoefAjout, A )
                 !----------------------- -----
                 ! 4. contribution du sommet js
                 !-----------------------------
                 dVKL = coef*(ln(V(is)) - ln(V(js))) + &
                      & coeta * ( ln(V(iK)) - ln(V(jL)) )
                 dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                 FKL = dVKLplus*( MuCroit(Xk(js)) + MuDecroit(Xk(is)) ) + &
                      & dVKLmoins*( MuCroit(Xk(is)) + MuDecroit(Xk(js)) )
                 F(js) = F(js) + RDVT* FKL
                 !! derivee de F(js) par rapport a Xk(js)
                 CoefAjout = RDVT*(dVKLplus*DerivMuCroit(Xk(js)) + dVKLmoins*DerivMuDecroit(Xk(js)))
                 call ajout(js,js,CoefAjout, A )
                 !! derivee de F(js) par rapport a Xk(is)
                 CoefAjout = RDVT*(dVKLplus*DerivMuDecroit(Xk(is)) + dVKLmoins*DerivMuCroit(Xk(is)))
                 call ajout(js,is,CoefAjout, A )
              ELSE
                 Select case ( Ntyps(js) ) 
                 case (dirichlet)
                    Ubord = Gb(js)
                    coef = TKL(iseg) ; coeta = etaSSe(iseg)
                    RDVT = (V(is) + Ubord + V(iK) + V(jL))/4
                    !-------------------------------
                    ! 1. contribution du triangle iK
                    !-------------------------------
                    dVKL = coef*(ln(V(jL)) - ln(V(iK))) + &
                         & coeta * ( ln(Ubord) - ln(V(is)) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    FKL = dVKLplus*( MuCroit(Xk(iK)) + MuDecroit(Xk(jL)) ) + &
                         & dVKLmoins*( MuCroit(Xk(jL)) + MuDecroit(Xk(iK)) )
                    F(iK) = F(iK) + RDVT* FKL
                    !! derivee de F(iK) par rapport a Xk(iK)
                    CoefAjout = RDVT*(dVKLplus*DerivMuCroit(Xk(iK)) + dVKLmoins*DerivMuDecroit(Xk(iK)))
                    call ajout(iK,iK,CoefAjout, A )
                    !! derivee de F(iK) par rapport a Xk(jL)
                    CoefAjout = RDVT*(dVKLplus*DerivMuDecroit(Xk(jL)) + dVKLmoins*DerivMuCroit(Xk(jL)))
                    call ajout(iK,jL,CoefAjout, A )
                    !-------------------------------
                    ! 2. contribution du triangle jL
                    !-------------------------------
                    dVKL = coef*(ln(V(iK)) - ln(V(jL))) + &
                         & coeta * ( ln(V(is)) - ln(Ubord) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    FKL = dVKLplus*( MuCroit(Xk(jL)) + MuDecroit(Xk(iK)) ) + &
                         & dVKLmoins*( MuCroit(Xk(iK)) + MuDecroit(Xk(jL)) )
                    F(jL) = F(jL) + RDVT* FKL
                    !! derivee de F(jL) par rapport a Xk(jL)
                    CoefAjout = RDVT*(dVKLplus*DerivMuCroit(Xk(jL)) + dVKLmoins*DerivMuDecroit(Xk(jL)))
                    call ajout(jL,jL,CoefAjout, A )
                    !! derivee de F(jL) par rapport a Xk(iK)
                    CoefAjout = RDVT*(dVKLplus*DerivMuDecroit(Xk(iK)) + dVKLmoins*DerivMuCroit(Xk(iK)))
                    call ajout(jL,iK,CoefAjout, A )
                    !-----------------------------
                    ! 3. contribution du sommet is
                    !-----------------------------
                    coef = TKeLe(iseg)
                    dVKL = coef*(ln(Ubord) - ln(V(is))) + &
                         & coeta * ( ln(V(jL)) - ln(V(iK)) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    FKL = dVKLplus*( MuCroit(Xk(is)) + MuDecroit(Ubord) ) + &
                         & dVKLmoins*( MuCroit(Ubord) + MuDecroit(Xk(is)) )
                    F(is) = F(is) + RDVT* FKL
                    !! derivee de F(is) par rapport a Xk(is)
                    CoefAjout = RDVT*(dVKLplus*DerivMuCroit(Xk(is)) + dVKLmoins*DerivMuDecroit(Xk(is)))
                    call ajout(is,is,CoefAjout, A )
                 case(neumann)
                    ! on fait rien
                 End Select
              END IF
              !
           ELSE
              !
              Select case ( Ntyps(is) )
              case (dirichlet)
                 IF (js <= NsInt) THEN
                    Ubord= Gb(is)
                    coef = TKL(iseg) ; coeta = etaSSe(iseg)
                    RDVT = (Ubord + V(js) + V(iK) + V(jL))/4
                    !-------------------------------
                    ! 1. contribution du triangle iK
                    !-------------------------------
                    dVKL = coef*(ln(V(jL)) - ln(V(iK))) + &
                         & coeta * ( ln(V(js)) - ln(Ubord) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    FKL = dVKLplus*( MuCroit(Xk(iK)) + MuDecroit(Xk(jL)) ) + &
                         & dVKLmoins*( MuCroit(Xk(jL)) + MuDecroit(Xk(iK)) )
                    F(iK) = F(iK) + RDVT* FKL
                    !! derivee de F(iK) par rapport a Xk(iK)
                    CoefAjout = RDVT*(dVKLplus*DerivMuCroit(Xk(iK)) + dVKLmoins*DerivMuDecroit(Xk(iK)))
                    call ajout(iK,iK,CoefAjout, A )
                    !! derivee de F(iK) par rapport a Xk(jL)
                    CoefAjout = RDVT*(dVKLplus*DerivMuDecroit(Xk(jL)) + dVKLmoins*DerivMuCroit(Xk(jL)))
                    call ajout(iK,jL,CoefAjout, A )
                    !-------------------------------
                    ! 2. contribution du triangle jL
                    !-------------------------------
                    dVKL = coef*(ln(V(iK)) - ln(V(jL))) + &
                         & coeta * ( ln(Ubord) - ln(V(js)) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    FKL = dVKLplus*( MuCroit(Xk(jL)) + MuDecroit(Xk(iK)) ) + &
                         & dVKLmoins*( MuCroit(Xk(iK)) + MuDecroit(Xk(jL)) )
                    F(jL) = F(jL) + RDVT* FKL
                    !! derivee de F(jL) par rapport a Xk(jL)
                    CoefAjout = RDVT*(dVKLplus*DerivMuCroit(Xk(jL)) + dVKLmoins*DerivMuDecroit(Xk(jL)))
                    call ajout(jL,jL,CoefAjout, A )
                    !! derivee de F(jL) par rapport a Xk(iK)
                    CoefAjout = RDVT*(dVKLplus*DerivMuDecroit(Xk(iK)) + dVKLmoins*DerivMuCroit(Xk(iK)))
                    call ajout(jL,iK,CoefAjout, A )
                    !-----------------------------
                    ! 3. contribution du sommet js
                    !-----------------------------
                    coef = TKeLe(iseg)
                    dVKL = coef*(ln(Ubord) - ln(V(js))) + &
                         & coeta * ( ln(V(iK)) - ln(V(jL)) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    FKL = dVKLplus*( MuCroit(Xk(js)) + MuDecroit(Ubord) ) + &
                         & dVKLmoins*( MuCroit(Ubord) + MuDecroit(Xk(js)) )
                    F(js) = F(js) + RDVT* FKL
                    !! derivee de F(js) par rapport a Xk(js)
                    CoefAjout = RDVT*(dVKLplus*DerivMuCroit(Xk(js)) + dVKLmoins*DerivMuDecroit(Xk(js)))
                    call ajout(js,js,CoefAjout, A )
                 ELSE ! js se situe encore au bord (un triangle entier sur le bord)
                    Uibord= Gb(is) ; Ujbord= Gb(js)
                    coef = TKL(iseg) ; coeta = etaSSe(iseg)
                    RDVT = (Uibord + Ujbord + V(iK) + V(jL))/4
                    !-------------------------------
                    ! 1. contribution du triangle iK
                    !-------------------------------
                    dVKL = coef*(ln(V(jL)) - ln(V(iK))) + &
                         & coeta * ( ln(Ujbord) - ln(Uibord) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    FKL = dVKLplus*( MuCroit(Xk(iK)) + MuDecroit(Xk(jL)) ) + &
                         & dVKLmoins*( MuCroit(Xk(jL)) + MuDecroit(Xk(iK)) )
                    F(iK) = F(iK) + RDVT* FKL
                    !! derivee de F(iK) par rapport a Xk(iK)
                    CoefAjout = RDVT*(dVKLplus*DerivMuCroit(Xk(iK)) + dVKLmoins*DerivMuDecroit(Xk(iK)))
                    call ajout(iK,iK,CoefAjout, A )
                    !! derivee de F(iK) par rapport a Xk(jL)
                    CoefAjout = RDVT*(dVKLplus*DerivMuDecroit(Xk(jL)) + dVKLmoins*DerivMuCroit(Xk(jL)))
                    call ajout(iK,jL,CoefAjout, A )
                    !-------------------------------
                    ! 2. contribution du triangle jL
                    !-------------------------------
                    dVKL = coef*(ln(V(iK)) - ln(V(jL))) + &
                         & coeta * ( ln(Uibord) - ln(Ujbord) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    FKL = dVKLplus*( MuCroit(Xk(jL)) + MuDecroit(Xk(iK)) ) + &
                         & dVKLmoins*( MuCroit(Xk(iK)) + MuDecroit(Xk(jL)) )
                    F(jL) = F(jL) + RDVT* FKL
                    !! derivee de F(jL) par rapport a Xk(jL)
                    CoefAjout = RDVT*(dVKLplus*DerivMuCroit(Xk(jL)) + dVKLmoins*DerivMuDecroit(Xk(jL)))
                    call ajout(jL,jL,CoefAjout, A )
                    !! derivee de F(jL) par rapport a Xk(iK)
                    CoefAjout = RDVT*(dVKLplus*DerivMuDecroit(Xk(iK)) + dVKLmoins*DerivMuCroit(Xk(iK)))
                    call ajout(jL,iK,CoefAjout, A )
                 END IF
              case(neumann)
                 ! On fait rien
              END Select
           END IF
           !!
        case (dirichlet) !! segment au bord  (triangle dégénéré)
           !!
           is = NuSeg(1,iseg); js = NuSeg(2,iseg)
           iK = NumTVoisSeg(1,iseg) + NsInt; !jL = numero milieu iseg
           Uibord= Gb(is) ; Ujbord= Gb(js) ; Ubord= Gb(iseg + Nbs)
           coef = TKL(iseg) ; coeta = etaSSe(iseg)
           RDVT = (Uibord + Ujbord + V(iK) + Ubord)/4
           !-------------------------------
           ! 1. contribution du triangle iK
           !-------------------------------
           dVKL = coef*(ln(Ubord) - ln(V(iK))) + &
                & coeta * ( ln(Ujbord) - ln(Uibord) )
           dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
           FKL = dVKLplus*( MuCroit(Xk(iK)) + MuDecroit(Ubord) ) + &
                & dVKLmoins*( MuCroit(Ubord) + MuDecroit(Xk(iK)) )
           F(iK) = F(iK) + RDVT* FKL
           !! derivee de F(iK) par rapport a Xk(iK)
           CoefAjout = RDVT*(dVKLplus*DerivMuCroit(Xk(iK)) + dVKLmoins*DerivMuDecroit(Xk(iK)))
           call ajout(iK,iK,CoefAjout, A )
           !!
        case (neumann)
           !!
           is = NuSeg(1,iseg); js = NuSeg(2,iseg)
           iK = NumTVoisSeg(1,iseg) + NsInt
           coef = TKL(iseg) ; coeta = etaSSe(iseg) ; coefe = TKeLe(iseg)
           RDVT = (3*(V(is) + V(js))/2 + V(iK))/4
           !-------------------------------
           ! 1. contribution du triangle iK
           !-------------------------------
           ! On fait rien
           !-----------------------------
           ! 2. contribution du sommet is
           !-----------------------------
           coefe = (coefe - (coeta**2)/coef)
           dVKL = coefe*(ln(V(js)) - ln(V(is)))
           dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
           FKL = dVKLplus*( MuCroit(Xk(is)) + MuDecroit(Xk(js)) ) + &
                & dVKLmoins*( MuCroit(Xk(js)) + MuDecroit(Xk(is)) )
           F(is) = F(is) + RDVT* FKL
           !! derivee de F(is) par rapport a Xk(is)
           CoefAjout = RDVT*(dVKLplus*DerivMuCroit(Xk(is)) + dVKLmoins*DerivMuDecroit(Xk(is)))
           call ajout(is,is,CoefAjout, A )
           !! derivee de F(is) par rapport a Xk(js)
           CoefAjout = RDVT*(dVKLplus*DerivMuDecroit(Xk(js)) + dVKLmoins*DerivMuCroit(Xk(js)))
           call ajout(is,js,CoefAjout, A )
           !-----------------------------
           ! 3. contribution du sommet js
           !-----------------------------
           dVKL = coefe*(ln(V(is)) - ln(V(js)))
           dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
           FKL = dVKLplus*( MuCroit(Xk(js)) + MuDecroit(Xk(is))) + &
                & dVKLmoins*( MuCroit(Xk(is)) + MuDecroit(Xk(js)))
           F(js) = F(js) + RDVT* FKL
           !! derivee de F(js) par rapport a Xk(js)
           CoefAjout = RDVT*(dVKLplus*DerivMuCroit(Xk(js)) + dVKLmoins*DerivMuDecroit(Xk(js)))
           call ajout(js,js,CoefAjout, A )
           !! derivee de F(js) par rapport a Xk(is)
           CoefAjout = RDVT*(dVKLplus*DerivMuDecroit(Xk(is)) + dVKLmoins*DerivMuCroit(Xk(is)))
           call ajout(js,is,CoefAjout, A )
        END Select
     END Do


     !!============================
     !!=============================
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
     print*,'kiter = ', kiter,'erreur NEWTON sqrt(sum(Yk*Yk))',sqrt(sum(Yk*Yk))
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
  E=Xk+Yk
  print*,'fin newton'
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


END SUBROUTINE NewtonuDDFV






