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
SUBROUTINE NewtoneDDFV(A,Eold,E,Um,Vm,ndim,choixf,temps)
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
                 coef = eTKL(iseg) ; coeta = eetaSSe(iseg)
                 !-------------------------------
                 ! 1. contribution du triangle iK
                 !-------------------------------
                 dVKL = coeta * ( xiE(Xk(is)) - xiE(Xk(js)) )
                 dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                 F(iK) = F(iK) + coef*(AdegenE(Xk(iK)) - AdegenE(Xk(jL))) + &
                      & dVKLplus*( rmCroitE(Xk(iK)) + rmDeCroitE(Xk(jL)) ) + &
                      & dVKLmoins*( rmCroitE(Xk(jL)) + rmDeCroitE(Xk(iK)))
                 !! derivee de F(iK) par rapport a Xk(iK)
                 CoefAjout = coef*DerivAdegE(Xk(iK))+ dVKLplus*DerrmCroitE(Xk(iK)) + dVKLmoins*DerrmDeCroitE(Xk(iK))
                 call ajout(iK,iK,CoefAjout, A )
                 !! derivee de F(iK) par rapport a Xk(jL)
                 CoefAjout = -coef*DerivAdegE(Xk(jL))+ dVKLplus*DerrmDeCroitE(Xk(jL)) + dVKLmoins*DerrmCroitE(Xk(jL))
                 call ajout(iK,jL,CoefAjout, A )
                 !! derivee de F(iK) par rapport a Xk(is)
                 CoefAjout = transfer((dVKL>0),1)*coeta* rmE(Xk(is))*( rmCroitE(Xk(iK)) + rmDeCroitE(Xk(jL)) ) &
                      & + transfer((dVKL<0),1)*coeta* rmE(Xk(is))*( rmCroitE(Xk(jL)) + rmDeCroitE(Xk(iK)) )
                 call ajout(iK,is, CoefAjout, A)
                 !! derivee de F(iK) par rapport a Xk(js)
                 CoefAjout = -transfer((dVKL>0),1)*coeta* rmE(Xk(js))*( rmCroitE(Xk(iK)) + rmDeCroitE(Xk(jL)) ) &
                      & - transfer((dVKL<0),1)*coeta* rmE(Xk(js))*( rmCroitE(Xk(jL)) + rmDeCroitE(Xk(iK)) )
                 call ajout(iK,js, CoefAjout, A)
                 !-------------------------------
                 ! 2. contribution du triangle jL
                 !-------------------------------
                 dVKL = coeta * ( xiE(Xk(js)) - xiE(Xk(is)) )
                 dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                 F(jL) = F(jL) + coef*(AdegenE(Xk(jL)) - AdegenE(Xk(iK))) + &
                      & dVKLplus*( rmCroitE(Xk(jL)) + rmDeCroitE(Xk(iK)) ) + &
                      & dVKLmoins*( rmCroitE(Xk(iK)) + rmDeCroitE(Xk(jL)))
                 !! derivee de F(jL) par rapport a Xk(iK)
                 CoefAjout = -coef*DerivAdegE(Xk(iK))+ dVKLplus*DerrmDeCroitE(Xk(iK)) + dVKLmoins*DerrmCroitE(Xk(iK))
                 call ajout(jL,iK,CoefAjout, A )
                 !! derivee de F(jL) par rapport a Xk(jL)
                 CoefAjout = coef*DerivAdegE(Xk(jL))+ dVKLplus*DerrmCroitE(Xk(jL)) + dVKLmoins*DerrmDeCroitE(Xk(jL))
                 call ajout(jL,jL,CoefAjout, A )
                 !! derivee de F(jL) par rapport a Xk(is)
                 CoefAjout = -transfer((dVKL>0),1)*coeta* rmE(Xk(is))*( rmCroitE(Xk(jL)) + rmDeCroitE(Xk(iK)) ) &
                      & -transfer((dVKL<0),1)*coeta* rmE(Xk(is))*( rmCroitE(Xk(iK)) + rmDeCroitE(Xk(jL)) )
                 call ajout(jL,is, CoefAjout, A)
                 !! derivee de F(jL) par rapport a Xk(js)
                 CoefAjout = transfer((dVKL>0),1)*coeta* rmE(Xk(js))*( rmCroitE(Xk(jL)) + rmDeCroitE(Xk(iK)) ) &
                      & + transfer((dVKL<0),1)*coeta* rmE(Xk(js))*( rmCroitE(Xk(iK)) + rmDeCroitE(Xk(jL)) )
                 call ajout(jL,js, CoefAjout, A)
                 !-----------------------------
                 ! 3. contribution du sommet is
                 !-----------------------------
                 coef = eTKeLe(iseg)
                 dVKL = coeta * ( xiE(Xk(iK)) - xiE(Xk(jL)) )
                 dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                 F(is) = F(is) + coef*(AdegenE(Xk(is)) - AdegenE(Xk(js))) + &
                      & dVKLplus*( rmCroitE(Xk(is)) + rmDeCroitE(Xk(js)) ) + &
                      & dVKLmoins*( rmCroitE(Xk(js)) + rmDeCroitE(Xk(is)))
                 !! derivee de F(is) par rapport a Xk(is)
                 CoefAjout = coef*DerivAdegE(Xk(is))+ dVKLplus*DerrmCroitE(Xk(is)) + dVKLmoins*DerrmDeCroitE(Xk(is))
                 call ajout(is,is,CoefAjout, A )
                 !! derivee de F(is) par rapport a Xk(js)
                 CoefAjout = -coef*DerivAdegE(Xk(js))+ dVKLplus*DerrmDeCroitE(Xk(js)) + dVKLmoins*DerrmCroitE(Xk(js))
                 call ajout(is,js,CoefAjout, A )
                 !! derivee de F(is) par rapport a Xk(iK)
                 CoefAjout = transfer((dVKL>0),1)*coeta* rmE(Xk(iK))*( rmCroitE(Xk(is)) + rmDeCroitE(Xk(js)) ) &
                      & + transfer((dVKL<0),1)*coeta* rmE(Xk(iK))*( rmCroitE(Xk(js)) + rmDeCroitE(Xk(is)) )
                 call ajout(is,iK, CoefAjout, A)
                 !! derivee de F(is) par rapport a Xk(jL)
                 CoefAjout = -transfer((dVKL>0),1)*coeta* rmE(Xk(jL))*( rmCroitE(Xk(is)) + rmDeCroitE(Xk(js)) ) &
                      & - transfer((dVKL<0),1)*coeta* rmE(Xk(jL))*( rmCroitE(Xk(js)) + rmDeCroitE(Xk(is)) )
                 call ajout(is,jL, CoefAjout, A)
                 !-----------------------------
                 ! 4. contribution du sommet js
                 !-----------------------------
                 dVKL = coeta * ( xiE(Xk(jL)) - xiE(Xk(iK)) )
                 dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                 F(js) = F(js) + coef*(AdegenE(Xk(js)) - AdegenE(Xk(is))) + &
                      & dVKLplus*( rmCroitE(Xk(js)) + rmDeCroitE(Xk(is)) ) + &
                      & dVKLmoins*( rmCroitE(Xk(is)) + rmDeCroitE(Xk(js)))
                 !! derivee de F(js) par rapport a Xk(is)
                 CoefAjout = -coef*DerivAdegE(Xk(is))+ dVKLplus*DerrmDeCroitE(Xk(is)) + dVKLmoins*DerrmCroitE(Xk(is))
                 call ajout(js,is,CoefAjout, A )
                 !! derivee de F(js) par rapport a Xk(js)
                 CoefAjout = coef*DerivAdegE(Xk(js))+ dVKLplus*DerrmCroitE(Xk(js)) + dVKLmoins*DerrmDeCroitE(Xk(js))
                 call ajout(js,js,CoefAjout, A )
                 !! derivee de F(js) par rapport a Xk(iK)
                 CoefAjout = -transfer((dVKL>0),1)*coeta* rmE(Xk(iK))*( rmCroitE(Xk(js)) + rmDeCroitE(Xk(is)) ) &
                      & - transfer((dVKL<0),1)*coeta* rmE(Xk(iK))*( rmCroitE(Xk(is)) + rmDeCroitE(Xk(js)) )
                 call ajout(js,iK, CoefAjout, A)
                 !! derivee de F(js) par rapport a Xk(jL)
                 CoefAjout = transfer((dVKL>0),1)*coeta* rmE(Xk(jL))*( rmCroitE(Xk(js)) + rmDeCroitE(Xk(is)) ) &
                      & + transfer((dVKL<0),1)*coeta* rmE(Xk(jL))*( rmCroitE(Xk(is)) + rmDeCroitE(Xk(js)) )
                 call ajout(js,jL, CoefAjout, A)
                 !
              ELSE
                 Select case ( Ntyps(js) ) 
                 case (dirichlet)
                    Ubord = Gb(js)
                    coef = eTKL(iseg) ; coeta = eetaSSe(iseg)
                    !-------------------------------
                    ! 1. contribution du triangle iK
                    !-------------------------------
                    dVKL = coeta * ( xiE(Xk(is)) - xiE(Ubord) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    F(iK) = F(iK) + coef*(AdegenE(Xk(iK)) - AdegenE(Xk(jL))) + &
                         & dVKLplus*( rmCroitE(Xk(iK)) + rmDeCroitE(Xk(jL)) ) + &
                         & dVKLmoins*( rmCroitE(Xk(jL)) + rmDeCroitE(Xk(iK)))
                    !! derivee de F(iK) par rapport a Xk(iK)
                    CoefAjout = coef*DerivAdegE(Xk(iK))+ dVKLplus*DerrmCroitE(Xk(iK)) + dVKLmoins*DerrmDeCroitE(Xk(iK))
                    call ajout(iK,iK,CoefAjout, A )
                    !! derivee de F(iK) par rapport a Xk(jL)
                    CoefAjout = -coef*DerivAdegE(Xk(jL))+ dVKLplus*DerrmDeCroitE(Xk(jL)) + dVKLmoins*DerrmCroitE(Xk(jL))
                    call ajout(iK,jL,CoefAjout, A )
                    !! derivee de F(iK) par rapport a Xk(is)
                    CoefAjout = transfer((dVKL>0),1)*coeta* rmE(Xk(is))*( rmCroitE(Xk(iK)) + rmDeCroitE(Xk(jL)) ) &
                         & + transfer((dVKL<0),1)*coeta* rmE(Xk(is))*( rmCroitE(Xk(jL)) + rmDeCroitE(Xk(iK)) )
                    call ajout(iK,is, CoefAjout, A)
                    !-------------------------------
                    ! 2. contribution du triangle jL
                    !-------------------------------
                    dVKL = coeta * ( xiE(Ubord) - xiE(Xk(is)) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    F(jL) = F(jL) + coef*(AdegenE(Xk(jL)) - AdegenE(Xk(iK))) + &
                         & dVKLplus*( rmCroitE(Xk(jL)) + rmDeCroitE(Xk(iK)) ) + &
                         & dVKLmoins*( rmCroitE(Xk(iK)) + rmDeCroitE(Xk(jL)))
                    !! derivee de F(jL) par rapport a Xk(iK)
                    CoefAjout = -coef*DerivAdegE(Xk(iK))+ dVKLplus*DerrmDeCroitE(Xk(iK)) + dVKLmoins*DerrmCroitE(Xk(iK))
                    call ajout(jL,iK,CoefAjout, A )
                    !! derivee de F(jL) par rapport a Xk(jL)
                    CoefAjout = coef*DerivAdegE(Xk(jL))+ dVKLplus*DerrmCroitE(Xk(jL)) + dVKLmoins*DerrmDeCroitE(Xk(jL))
                    call ajout(jL,jL,CoefAjout, A )
                    !! derivee de F(jL) par rapport a Xk(is)
                    CoefAjout = -transfer((dVKL>0),1)*coeta* rmE(Xk(is))*( rmCroitE(Xk(jL)) + rmDeCroitE(Xk(iK)) ) &
                         & - transfer((dVKL<0),1)*coeta* rmE(Xk(is))*( rmCroitE(Xk(iK)) + rmDeCroitE(Xk(jL)) )
                    call ajout(jL,is, CoefAjout, A)
                    !-----------------------------
                    ! 3. contribution du sommet is
                    !-----------------------------
                    coef = eTKeLe(iseg)
                    dVKL = coeta * ( xiE(Xk(iK)) - xiE(Xk(jL)) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    F(is) = F(is) + coef*(AdegenE(Xk(is)) - AdegenE(Ubord)) + &
                         & dVKLplus*( rmCroitE(Xk(is)) + rmDeCroitE(Ubord) ) + &
                         & dVKLmoins*( rmCroitE(Ubord) + rmDeCroitE(Xk(is)))
                    !! derivee de F(is) par rapport a Xk(is)
                    CoefAjout = coef*DerivAdegE(Xk(is))+ dVKLplus*DerrmCroitE(Xk(is)) + dVKLmoins*DerrmDeCroitE(Xk(is))
                    call ajout(is,is,CoefAjout, A )
                    !! derivee de F(is) par rapport a Xk(iK)
                    CoefAjout = transfer((dVKL>0),1)*coeta* rmE(Xk(iK))*( rmCroitE(Xk(is)) + rmDeCroitE(Ubord) ) &
                         & + transfer((dVKL<0),1)*coeta* rmE(Xk(iK))*( rmCroitE(Ubord) + rmDeCroitE(Xk(is)) )
                    call ajout(is,iK, CoefAjout, A)
                    !! derivee de F(is) par rapport a Xk(jL)
                    CoefAjout = -transfer((dVKL>0),1)*coeta* rmE(Xk(jL))*( rmCroitE(Xk(is)) + rmDeCroitE(Ubord) ) &
                         & - transfer((dVKL<0),1)*coeta* rmE(Xk(jL))*( rmCroitE(Ubord) + rmDeCroitE(Xk(is)) )
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
                    coef = eTKL(iseg) ; coeta = eetaSSe(iseg)
                    !-------------------------------
                    ! 1. contribution du triangle iK
                    !-------------------------------
                    dVKL = coeta * ( xiE(Ubord) - xiE(Xk(js)) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    F(iK) = F(iK) + coef*(AdegenE(Xk(iK)) - AdegenE(Xk(jL))) + &
                         & dVKLplus*( rmCroitE(Xk(iK)) + rmDeCroitE(Xk(jL)) ) + &
                         & dVKLmoins*( rmCroitE(Xk(jL)) + rmDeCroitE(Xk(iK)))
                    !! derivee de F(iK) par rapport a Xk(iK)
                    CoefAjout = coef*DerivAdegE(Xk(iK))+ dVKLplus*DerrmCroitE(Xk(iK)) + dVKLmoins*DerrmDeCroitE(Xk(iK))
                    call ajout(iK,iK,CoefAjout, A )
                    !! derivee de F(iK) par rapport a Xk(jL)
                    CoefAjout = -coef*DerivAdegE(Xk(jL))+ dVKLplus*DerrmDeCroitE(Xk(jL)) + dVKLmoins*DerrmCroitE(Xk(jL))
                    call ajout(iK,jL,CoefAjout, A )
                    !! derivee de F(iK) par rapport a Xk(js)
                    CoefAjout = -transfer((dVKL>0),1)*coeta* rmE(Xk(js))*( rmCroitE(Xk(iK)) + rmDeCroitE(Xk(jL)) )&
                         & - transfer((dVKL<0),1)*coeta* rmE(Xk(js))*( rmCroitE(Xk(jL)) + rmDeCroitE(Xk(iK)) )
                    call ajout(iK,js, CoefAjout, A)
                    !-------------------------------
                    ! 2. contribution du triangle jL
                    !-------------------------------
                    dVKL = coeta * ( xiE(Xk(js)) - xiE(Ubord) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    F(jL) = F(jL) + coef*(AdegenE(Xk(jL)) - AdegenE(Xk(iK))) + &
                         & dVKLplus*( rmCroitE(Xk(jL)) + rmDeCroitE(Xk(iK)) ) + &
                         & dVKLmoins*( rmCroitE(Xk(iK)) + rmDeCroitE(Xk(jL)))
                    !! derivee de F(jL) par rapport a Xk(iK)
                    CoefAjout = -coef*DerivAdegE(Xk(iK))+ dVKLplus*DerrmDeCroitE(Xk(iK)) + dVKLmoins*DerrmCroitE(Xk(iK))
                    call ajout(jL,iK,CoefAjout, A )
                    !! derivee de F(jL) par rapport a Xk(jL)
                    CoefAjout = coef*DerivAdegE(Xk(jL))+ dVKLplus*DerrmCroitE(Xk(jL)) + dVKLmoins*DerrmDeCroitE(Xk(jL))
                    call ajout(jL,jL,CoefAjout, A )
                    !! derivee de F(jL) par rapport a Xk(js)
                    CoefAjout = transfer((dVKL>0),1)*coeta* rmE(Xk(js))*( rmCroitE(Xk(jL)) + rmDeCroitE(Xk(iK)) ) &
                         & + transfer((dVKL<0),1)*coeta* rmE(Xk(js))*( rmCroitE(Xk(iK)) + rmDeCroitE(Xk(jL)) )
                    call ajout(jL,js, CoefAjout, A)
                    !-----------------------------
                    ! 3. contribution du sommet js
                    !-----------------------------
                    coef = eTKeLe(iseg)
                    dVKL = coeta * ( xiE(Xk(jL)) - xiE(Xk(iK)) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    F(js) = F(js) + coef*(AdegenE(Xk(js)) - AdegenE(Ubord)) + &
                         & dVKLplus*( rmCroitE(Xk(js)) + rmDeCroitE(Ubord) ) + &
                         & dVKLmoins*( rmCroitE(Ubord) + rmDeCroitE(Xk(js)))
                    !! derivee de F(js) par rapport a Xk(js)
                    CoefAjout = coef*DerivAdegE(Xk(js))+ dVKLplus*DerrmCroitE(Xk(js)) + dVKLmoins*DerrmDeCroitE(Xk(js))
                    call ajout(js,js,CoefAjout, A )
                    !! derivee de F(js) par rapport a Xk(iK)
                    CoefAjout = -transfer((dVKL>0),1)*coeta* rmE(Xk(iK))*( rmCroitE(Xk(js)) + rmDeCroitE(Ubord) ) &
                         & - transfer((dVKL<0),1)*coeta* rmE(Xk(iK))*( rmCroitE(Ubord) + rmDeCroitE(Xk(js)) )
                    call ajout(js,iK, CoefAjout, A)
                    !! derivee de F(js) par rapport a Xk(jL)
                    CoefAjout = transfer((dVKL>0),1)*coeta* rmE(Xk(jL))*( rmCroitE(Xk(js)) + rmDeCroitE(Ubord) ) &
                         & + transfer((dVKL<0),1)*coeta* rmE(Xk(jL))*( rmCroitE(Ubord) + rmDeCroitE(Xk(js)) )
                    call ajout(js,jL, CoefAjout, A)
                 ELSE ! js se situe encore au bord (un triangle entier sur le bord)
                    Uibord= Gb(is) ; Ujbord= Gb(js)
                    coef = eTKL(iseg) ; coeta = eetaSSe(iseg)
                    !-------------------------------
                    ! 1. contribution du triangle iK
                    !-------------------------------
                    dVKL = coeta * ( xiE(Uibord) - xiE(Ujbord) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    F(iK) = F(iK) + coef*(AdegenE(Xk(iK)) - AdegenE(Xk(jL))) + &
                         & dVKLplus*( rmCroitE(Xk(iK)) + rmDeCroitE(Xk(jL)) ) + &
                         & dVKLmoins*( rmCroitE(Xk(jL)) + rmDeCroitE(Xk(iK)))
                    !! derivee de F(iK) par rapport a Xk(iK)
                    CoefAjout = coef*DerivAdegE(Xk(iK))+ dVKLplus*DerrmCroitE(Xk(iK)) + dVKLmoins*DerrmDeCroitE(Xk(iK))
                    call ajout(iK,iK,CoefAjout, A )
                    !! derivee de F(iK) par rapport a Xk(jL)
                    CoefAjout = -coef*DerivAdegE(Xk(jL))+ dVKLplus*DerrmDeCroitE(Xk(jL)) + dVKLmoins*DerrmCroitE(Xk(jL))
                    call ajout(iK,jL,CoefAjout, A )
                    !-------------------------------
                    ! 2. contribution du triangle jL
                    !-------------------------------
                    dVKL = coeta * ( xiE(Ujbord) - xiE(Uibord) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    F(jL) = F(jL) + coef*(AdegenE(Xk(jL)) - AdegenE(Xk(iK))) + &
                         & dVKLplus*( rmCroitE(Xk(jL)) + rmDeCroitE(Xk(iK)) ) + &
                         & dVKLmoins*( rmCroitE(Xk(iK)) + rmDeCroitE(Xk(jL)))
                    !! derivee de F(jL) par rapport a Xk(iK)
                    CoefAjout = -coef*DerivAdegE(Xk(iK))+ dVKLplus*DerrmDeCroitE(Xk(iK)) + dVKLmoins*DerrmCroitE(Xk(iK))
                    call ajout(jL,iK,CoefAjout, A )
                    !! derivee de F(jL) par rapport a Xk(jL)
                    CoefAjout = coef*DerivAdegE(Xk(jL))+ dVKLplus*DerrmCroitE(Xk(jL)) + dVKLmoins*DerrmDeCroitE(Xk(jL))
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
           coef = eTKL(iseg) ; coeta = eetaSSe(iseg)
           !-------------------------------
           ! 1. contribution du triangle iK
           !-------------------------------
           dVKL = coeta * ( xiE(Uibord) - xiE(Ujbord) )
           dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
           F(iK) = F(iK) + coef*(AdegenE(Xk(iK)) - AdegenE(Ubord)) + &
                & dVKLplus*( rmCroitE(Xk(iK)) + rmDeCroitE(Ubord) ) + &
                & dVKLmoins*( rmCroitE(Ubord) + rmDeCroitE(Xk(iK)))
           !! derivee de F(iK) par rapport a Xk(iK)
           CoefAjout = coef*DerivAdegE(Xk(iK))+ dVKLplus*DerrmCroitE(Xk(iK)) + dVKLmoins*DerrmDeCroitE(Xk(iK))
           call ajout(iK,iK,CoefAjout, A )
           !!
        case (neumann)
           !!
           is = NuSeg(1,iseg); js = NuSeg(2,iseg)
           iK = NumTVoisSeg(1,iseg) + NsInt
           coef = eTKL(iseg) ; coeta = eetaSSe(iseg) ; coefe = eTKeLe(iseg)
           !-------------------------------
           ! 1. contribution du triangle iK
           !-------------------------------
           ! On fait rien
           !-----------------------------
           ! 2. contribution du sommet is
           !----------------------------- 
           dVKL = ((coeta**2)/coef) * ( xiE(Xk(js)) - xiE(Xk(is)) )
           dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
           F(is) = F(is) + coefe*(AdegenE(Xk(is)) - AdegenE(Xk(js))) + &
                & dVKLplus*( rmCroitE(Xk(is)) + rmDeCroitE(Xk(js)) ) + &
                & dVKLmoins*( rmCroitE(Xk(js)) + rmDeCroitE(Xk(is)))
           !! derivee de F(is) par rapport a Xk(is)
           CoefAjout = coefe*DerivAdegE(Xk(is))+ dVKLplus*DerrmCroitE(Xk(is)) + dVKLmoins*DerrmDeCroitE(Xk(is)) &
                & - transfer((dVKL>0),1)*((coeta**2)/coef)* rmE(Xk(is))*( rmCroitE(Xk(is)) + rmDeCroitE(Xk(js)))&
                & - transfer((dVKL<0),1)*((coeta**2)/coef)* rmE(Xk(is))*( rmCroitE(Xk(js)) + rmDeCroitE(Xk(is)))
           call ajout(is,is,CoefAjout, A )
           !! derivee de F(is) par rapport a Xk(js)
           CoefAjout = -coefe*DerivAdegE(Xk(js))+ dVKLplus*DerrmDeCroitE(Xk(js)) + dVKLmoins*DerrmCroitE(Xk(js))&
                & + transfer((dVKL>0),1)*((coeta**2)/coef)* rmE(Xk(js))*( rmCroitE(Xk(is)) + rmDeCroitE(Xk(js)))&
                & + transfer((dVKL<0),1)*((coeta**2)/coef)* rmE(Xk(js))*( rmCroitE(Xk(js)) + rmDeCroitE(Xk(is)))
           call ajout(is,js,CoefAjout, A )
           !-----------------------------
           ! 3. contribution du sommet js
           !-----------------------------
           dVKL = ((coeta**2)/coef) * ( xiE(Xk(is)) - xiE(Xk(js)) )
           dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
           F(js) = F(js) + coefe*(AdegenE(Xk(js)) - AdegenE(Xk(is))) + &
                & dVKLplus*( rmCroitE(Xk(js)) + rmDeCroitE(Xk(is)) ) + &
                & dVKLmoins*( rmCroitE(Xk(is)) + rmDeCroitE(Xk(js)))
           !! derivee de F(js) par rapport a Xk(is)
           CoefAjout = -coefe*DerivAdegE(Xk(is))+ dVKLplus*DerrmDeCroitE(Xk(is)) + dVKLmoins*DerrmCroitE(Xk(is))&
                & + transfer((dVKL>0),1)*((coeta**2)/coef)* rmE(Xk(is))*( rmCroitE(Xk(js)) + rmDeCroitE(Xk(is)))&
                & + transfer((dVKL<0),1)*((coeta**2)/coef)* rmE(Xk(is))*( rmCroitE(Xk(is)) + rmDeCroitE(Xk(js)))
           call ajout(js,is,CoefAjout, A )
           !! derivee de F(js) par rapport a Xk(js)
           CoefAjout = coefe*DerivAdegE(Xk(js))+ dVKLplus*DerrmCroitE(Xk(js)) + dVKLmoins*DerrmDeCroitE(Xk(js)) &
                & - transfer((dVKL>0),1)*((coeta**2)/coef)* rmE(Xk(js))*( rmCroitE(Xk(js)) + rmDeCroitE(Xk(is)))&
                & - transfer((dVKL<0),1)*((coeta**2)/coef)* rmE(Xk(js))*( rmCroitE(Xk(is)) + rmDeCroitE(Xk(js)))
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
                 coef = eTKL(iseg) ; coeta = eetaSSe(iseg)
                 RDVT = (Vm(is) + Vm(js) + Vm(iK) + Vm(jL))/4
                 !-------------------------------
                 ! 1. contribution du triangle iK
                 !-------------------------------
                 dVKL = coef*(ln(Vm(jL)) - ln(Vm(iK))) + &
                      & coeta*(ln(Vm(js)) - ln(Vm(is)))
                 dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                 FKL = dVKLplus*( MuCroitE(Xk(iK)) + MuDeCroitE(Xk(jL)) ) + &
                      & dVKLmoins*( MuCroitE(Xk(jL)) + MuDeCroitE(Xk(iK)) )
                 F(iK) = F(iK) + RDVT* FKL
                 !! derivee de F(iK) par rapport a Xk(iK)
                 CoefAjout = RDVT*(dVKLplus*DerivMuCroitE(Xk(iK)) + dVKLmoins*DerivMuDeCroitE(Xk(iK)))
                 call ajout(iK,iK,CoefAjout, A )
                 !! derivee de F(iK) par rapport a Xk(jL)
                 CoefAjout = RDVT*(dVKLplus*DerivMuDeCroitE(Xk(jL)) + dVKLmoins*DerivMuCroitE(Xk(jL)))
                 call ajout(iK,jL,CoefAjout, A )
                 !-------------------------------
                 ! 2. contribution du triangle jL
                 !-------------------------------
                 dVKL = coef*(ln(Vm(iK)) - ln(Vm(jL))) + &
                      & coeta * ( ln(Vm(is)) - ln(Vm(js)) )
                 dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                 FKL = dVKLplus*( MuCroitE(Xk(jL)) + MuDeCroitE(Xk(iK)) ) + &
                      & dVKLmoins*( MuCroitE(Xk(iK)) + MuDeCroitE(Xk(jL)) )
                 F(jL) = F(jL) + RDVT* FKL
                 !! derivee de F(jL) par rapport a Xk(jL)
                 CoefAjout = RDVT*(dVKLplus*DerivMuCroitE(Xk(jL)) + dVKLmoins*DerivMuDeCroitE(Xk(jL)))
                 call ajout(jL,jL,CoefAjout, A )
                 !! derivee de F(jL) par rapport a Xk(iK)
                 CoefAjout = RDVT*(dVKLplus*DerivMuDeCroitE(Xk(iK)) + dVKLmoins*DerivMuCroitE(Xk(iK)))
                 call ajout(jL,iK,CoefAjout, A )
                 !-----------------------------
                 ! 3. contribution du sommet is
                 !-----------------------------
                 coef = eTKeLe(iseg)
                 dVKL = coef*(ln(Vm(js)) - ln(Vm(is))) + &
                      & coeta * ( ln(Vm(jL)) - ln(Vm(iK)) )
                 dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                 FKL = dVKLplus*( MuCroitE(Xk(is)) + MuDeCroitE(Xk(js)) ) + &
                      & dVKLmoins*( MuCroitE(Xk(js)) + MuDeCroitE(Xk(is)) )
                 F(is) = F(is) + RDVT* FKL
                 !! derivee de F(is) par rapport a Xk(is)
                 CoefAjout = RDVT*(dVKLplus*DerivMuCroitE(Xk(is)) + dVKLmoins*DerivMuDeCroitE(Xk(is)))
                 call ajout(is,is,CoefAjout, A )
                 !! derivee de F(is) par rapport a Xk(js)
                 CoefAjout = RDVT*(dVKLplus*DerivMuDeCroitE(Xk(js)) + dVKLmoins*DerivMuCroitE(Xk(js)))
                 call ajout(is,js,CoefAjout, A )
                 !----------------------- -----
                 ! 4. contribution du sommet js
                 !-----------------------------
                 dVKL = coef*(ln(Vm(is)) - ln(Vm(js))) + &
                      & coeta * ( ln(Vm(iK)) - ln(Vm(jL)) )
                 dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                 FKL = dVKLplus*( MuCroitE(Xk(js)) + MuDeCroitE(Xk(is)) ) + &
                      & dVKLmoins*( MuCroitE(Xk(is)) + MuDeCroitE(Xk(js)) )
                 F(js) = F(js) + RDVT* FKL
                 !! derivee de F(js) par rapport a Xk(js)
                 CoefAjout = RDVT*(dVKLplus*DerivMuCroitE(Xk(js)) + dVKLmoins*DerivMuDeCroitE(Xk(js)))
                 call ajout(js,js,CoefAjout, A )
                 !! derivee de F(js) par rapport a Xk(is)
                 CoefAjout = RDVT*(dVKLplus*DerivMuDeCroitE(Xk(is)) + dVKLmoins*DerivMuCroitE(Xk(is)))
                 call ajout(js,is,CoefAjout, A )
              ELSE
                 Select case ( Ntyps(js) ) 
                 case (dirichlet)
                    Ubord = Gb(js)
                    coef = eTKL(iseg) ; coeta = eetaSSe(iseg)
                    RDVT = (Vm(is) + Ubord + Vm(iK) + Vm(jL))/4
                    !-------------------------------
                    ! 1. contribution du triangle iK
                    !-------------------------------
                    dVKL = coef*(ln(Vm(jL)) - ln(Vm(iK))) + &
                         & coeta * ( ln(Ubord) - ln(Vm(is)) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    FKL = dVKLplus*( MuCroitE(Xk(iK)) + MuDeCroitE(Xk(jL)) ) + &
                         & dVKLmoins*( MuCroitE(Xk(jL)) + MuDeCroitE(Xk(iK)) )
                    F(iK) = F(iK) + RDVT* FKL
                    !! derivee de F(iK) par rapport a Xk(iK)
                    CoefAjout = RDVT*(dVKLplus*DerivMuCroitE(Xk(iK)) + dVKLmoins*DerivMuDeCroitE(Xk(iK)))
                    call ajout(iK,iK,CoefAjout, A )
                    !! derivee de F(iK) par rapport a Xk(jL)
                    CoefAjout = RDVT*(dVKLplus*DerivMuDeCroitE(Xk(jL)) + dVKLmoins*DerivMuCroitE(Xk(jL)))
                    call ajout(iK,jL,CoefAjout, A )
                    !-------------------------------
                    ! 2. contribution du triangle jL
                    !-------------------------------
                    dVKL = coef*(ln(Vm(iK)) - ln(Vm(jL))) + &
                         & coeta * ( ln(Vm(is)) - ln(Ubord) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    FKL = dVKLplus*( MuCroitE(Xk(jL)) + MuDeCroitE(Xk(iK)) ) + &
                         & dVKLmoins*( MuCroitE(Xk(iK)) + MuDeCroitE(Xk(jL)) )
                    F(jL) = F(jL) + RDVT* FKL
                    !! derivee de F(jL) par rapport a Xk(jL)
                    CoefAjout = RDVT*(dVKLplus*DerivMuCroitE(Xk(jL)) + dVKLmoins*DerivMuDeCroitE(Xk(jL)))
                    call ajout(jL,jL,CoefAjout, A )
                    !! derivee de F(jL) par rapport a Xk(iK)
                    CoefAjout = RDVT*(dVKLplus*DerivMuDeCroitE(Xk(iK)) + dVKLmoins*DerivMuCroitE(Xk(iK)))
                    call ajout(jL,iK,CoefAjout, A )
                    !-----------------------------
                    ! 3. contribution du sommet is
                    !-----------------------------
                    coef = eTKeLe(iseg)
                    dVKL = coef*(ln(Ubord) - ln(Vm(is))) + &
                         & coeta * ( ln(Vm(jL)) - ln(Vm(iK)) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    FKL = dVKLplus*( MuCroitE(Xk(is)) + MuDeCroitE(Ubord) ) + &
                         & dVKLmoins*( MuCroitE(Ubord) + MuDeCroitE(Xk(is)) )
                    F(is) = F(is) + RDVT* FKL
                    !! derivee de F(is) par rapport a Xk(is)
                    CoefAjout = RDVT*(dVKLplus*DerivMuCroitE(Xk(is)) + dVKLmoins*DerivMuDeCroitE(Xk(is)))
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
                    coef = eTKL(iseg) ; coeta = eetaSSe(iseg)
                    RDVT = (Ubord + Vm(js) + Vm(iK) + Vm(jL))/4
                    !-------------------------------
                    ! 1. contribution du triangle iK
                    !-------------------------------
                    dVKL = coef*(ln(Vm(jL)) - ln(Vm(iK))) + &
                         & coeta * ( ln(Vm(js)) - ln(Ubord) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    FKL = dVKLplus*( MuCroitE(Xk(iK)) + MuDeCroitE(Xk(jL)) ) + &
                         & dVKLmoins*( MuCroitE(Xk(jL)) + MuDeCroitE(Xk(iK)) )
                    F(iK) = F(iK) + RDVT* FKL
                    !! derivee de F(iK) par rapport a Xk(iK)
                    CoefAjout = RDVT*(dVKLplus*DerivMuCroitE(Xk(iK)) + dVKLmoins*DerivMuDeCroitE(Xk(iK)))
                    call ajout(iK,iK,CoefAjout, A )
                    !! derivee de F(iK) par rapport a Xk(jL)
                    CoefAjout = RDVT*(dVKLplus*DerivMuDeCroitE(Xk(jL)) + dVKLmoins*DerivMuCroitE(Xk(jL)))
                    call ajout(iK,jL,CoefAjout, A )
                    !-------------------------------
                    ! 2. contribution du triangle jL
                    !-------------------------------
                    dVKL = coef*(ln(Vm(iK)) - ln(Vm(jL))) + &
                         & coeta * ( ln(Ubord) - ln(Vm(js)) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    FKL = dVKLplus*( MuCroitE(Xk(jL)) + MuDeCroitE(Xk(iK)) ) + &
                         & dVKLmoins*( MuCroitE(Xk(iK)) + MuDeCroitE(Xk(jL)) )
                    F(jL) = F(jL) + RDVT* FKL
                    !! derivee de F(jL) par rapport a Xk(jL)
                    CoefAjout = RDVT*(dVKLplus*DerivMuCroitE(Xk(jL)) + dVKLmoins*DerivMuDeCroitE(Xk(jL)))
                    call ajout(jL,jL,CoefAjout, A )
                    !! derivee de F(jL) par rapport a Xk(iK)
                    CoefAjout = RDVT*(dVKLplus*DerivMuDeCroitE(Xk(iK)) + dVKLmoins*DerivMuCroitE(Xk(iK)))
                    call ajout(jL,iK,CoefAjout, A )
                    !-----------------------------
                    ! 3. contribution du sommet js
                    !-----------------------------
                    coef = eTKeLe(iseg)
                    dVKL = coef*(ln(Ubord) - ln(Vm(js))) + &
                         & coeta * ( ln(Vm(iK)) - ln(Vm(jL)) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    FKL = dVKLplus*( MuCroitE(Xk(js)) + MuDeCroitE(Ubord) ) + &
                         & dVKLmoins*( MuCroitE(Ubord) + MuDeCroitE(Xk(js)) )
                    F(js) = F(js) + RDVT* FKL
                    !! derivee de F(js) par rapport a Xk(js)
                    CoefAjout = RDVT*(dVKLplus*DerivMuCroitE(Xk(js)) + dVKLmoins*DerivMuDeCroitE(Xk(js)))
                    call ajout(js,js,CoefAjout, A )
                 ELSE ! js se situe encore au bord (un triangle entier sur le bord)
                    Uibord= Gb(is) ; Ujbord= Gb(js)
                    coef = eTKL(iseg) ; coeta = eetaSSe(iseg)
                    RDVT = (Uibord + Ujbord + Vm(iK) + Vm(jL))/4
                    !-------------------------------
                    ! 1. contribution du triangle iK
                    !-------------------------------
                    dVKL = coef*(ln(Vm(jL)) - ln(Vm(iK))) + &
                         & coeta * ( ln(Ujbord) - ln(Uibord) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    FKL = dVKLplus*( MuCroitE(Xk(iK)) + MuDeCroitE(Xk(jL)) ) + &
                         & dVKLmoins*( MuCroitE(Xk(jL)) + MuDeCroitE(Xk(iK)) )
                    F(iK) = F(iK) + RDVT* FKL
                    !! derivee de F(iK) par rapport a Xk(iK)
                    CoefAjout = RDVT*(dVKLplus*DerivMuCroitE(Xk(iK)) + dVKLmoins*DerivMuDeCroitE(Xk(iK)))
                    call ajout(iK,iK,CoefAjout, A )
                    !! derivee de F(iK) par rapport a Xk(jL)
                    CoefAjout = RDVT*(dVKLplus*DerivMuDeCroitE(Xk(jL)) + dVKLmoins*DerivMuCroitE(Xk(jL)))
                    call ajout(iK,jL,CoefAjout, A )
                    !-------------------------------
                    ! 2. contribution du triangle jL
                    !-------------------------------
                    dVKL = coef*(ln(Vm(iK)) - ln(Vm(jL))) + &
                         & coeta * ( ln(Uibord) - ln(Ujbord) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    FKL = dVKLplus*( MuCroitE(Xk(jL)) + MuDeCroitE(Xk(iK)) ) + &
                         & dVKLmoins*( MuCroitE(Xk(iK)) + MuDeCroitE(Xk(jL)) )
                    F(jL) = F(jL) + RDVT* FKL
                    !! derivee de F(jL) par rapport a Xk(jL)
                    CoefAjout = RDVT*(dVKLplus*DerivMuCroitE(Xk(jL)) + dVKLmoins*DerivMuDeCroitE(Xk(jL)))
                    call ajout(jL,jL,CoefAjout, A )
                    !! derivee de F(jL) par rapport a Xk(iK)
                    CoefAjout = RDVT*(dVKLplus*DerivMuDeCroitE(Xk(iK)) + dVKLmoins*DerivMuCroitE(Xk(iK)))
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
           coef = eTKL(iseg) ; coeta = eetaSSe(iseg)
           RDVT = (Uibord + Ujbord + Vm(iK) + Ubord)/4
           !-------------------------------
           ! 1. contribution du triangle iK
           !-------------------------------
           dVKL = coef*(ln(Ubord) - ln(Vm(iK))) + &
                & coeta * ( ln(Ujbord) - ln(Uibord) )
           dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
           FKL = dVKLplus*( MuCroitE(Xk(iK)) + MuDeCroitE(Ubord) ) + &
                & dVKLmoins*( MuCroitE(Ubord) + MuDeCroitE(Xk(iK)) )
           F(iK) = F(iK) + RDVT* FKL
           !! derivee de F(iK) par rapport a Xk(iK)
           CoefAjout = RDVT*(dVKLplus*DerivMuCroitE(Xk(iK)) + dVKLmoins*DerivMuDeCroitE(Xk(iK)))
           call ajout(iK,iK,CoefAjout, A )
           !!
        case (neumann)
           !!
           is = NuSeg(1,iseg); js = NuSeg(2,iseg)
           iK = NumTVoisSeg(1,iseg) + NsInt
           coef = eTKL(iseg) ; coeta = eetaSSe(iseg) ; coefe = eTKeLe(iseg)
           RDVT = (3*(Vm(is) + Vm(js))/2 + Vm(iK))/4
           !-------------------------------
           ! 1. contribution du triangle iK
           !-------------------------------
           ! On fait rien
           !-----------------------------
           ! 2. contribution du sommet is
           !-----------------------------
           coefe = (coefe - (coeta**2)/coef)
           dVKL = coefe*(ln(Vm(js)) - ln(Vm(is)))
           dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
           FKL = dVKLplus*( MuCroitE(Xk(is)) + MuDeCroitE(Xk(js)) ) + &
                & dVKLmoins*( MuCroitE(Xk(js)) + MuDeCroitE(Xk(is)) )
           F(is) = F(is) + RDVT* FKL
           !! derivee de F(is) par rapport a Xk(is)
           CoefAjout = RDVT*(dVKLplus*DerivMuCroitE(Xk(is)) + dVKLmoins*DerivMuDeCroitE(Xk(is)))
           call ajout(is,is,CoefAjout, A )
           !! derivee de F(is) par rapport a Xk(js)
           CoefAjout = RDVT*(dVKLplus*DerivMuDeCroitE(Xk(js)) + dVKLmoins*DerivMuCroitE(Xk(js)))
           call ajout(is,js,CoefAjout, A )
           !-----------------------------
           ! 3. contribution du sommet js
           !-----------------------------
           dVKL = coefe*(ln(Vm(is)) - ln(Vm(js)))
           dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
           FKL = dVKLplus*( MuCroitE(Xk(js)) + MuDeCroitE(Xk(is))) + &
                & dVKLmoins*( MuCroitE(Xk(is)) + MuDeCroitE(Xk(js)))
           F(js) = F(js) + RDVT* FKL
           !! derivee de F(js) par rapport a Xk(js)
           CoefAjout = RDVT*(dVKLplus*DerivMuCroitE(Xk(js)) + dVKLmoins*DerivMuDeCroitE(Xk(js)))
           call ajout(js,js,CoefAjout, A )
           !! derivee de F(js) par rapport a Xk(is)
           CoefAjout = RDVT*(dVKLplus*DerivMuDeCroitE(Xk(is)) + dVKLmoins*DerivMuCroitE(Xk(is)))
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


END SUBROUTINE NewtoneDDFV






