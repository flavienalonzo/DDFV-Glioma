!                ******************************
!                **  SUBROUTINE NewtonuDDFV  **
!                ******************************
! *****************************************************************************
!     *-!!!!!!!!!@@@@@!!@@@@@!!@@@@@!!@@@@@!!!!!!!!!!!-*
!     *-!!!!!!!!!!!!!@!!@!!!@!!!!!!@!!!!!!@!!!!!!!!!!!-*
!     *-!!!!!!!!!@@@@@!!@!!!@!!@@@@@!!@@@@@!!!!!!!!!!!-*
!     *-!!!!!!!!!@!!!!!!@!!!@!!@@!!!!!!!!!@!!!!!!!!!!!-*
!     *-!!!!!!!!!@@@@@!!@@@@@!!@@@@@!!@@@@@!!!!!!!!!!!-*
!------------------------------------------------------------
!  Dans ce sous-programme on discretise l'equation en U par un scheme  implicite 
!       pour le scheme DDFV auto 
!
! ******************************************************************************
! 
SUBROUTINE NewtonuDDFV(A,Uold,U,Cm,Em,ndim,choixf,temps)
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
  REAL(kind=long), dimension (ndim), intent(in)    :: Uold,Cm,Em
  REAL(kind=long), dimension (ndim), intent(out)   :: U
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
  Xk = Uold
  
  DO kiter = 1, kitermax
     if ((minval(Xk)<0.D0) .or. (maxval(Xk)>1.D0)) then
          print*,count((Xk<0.D0))
          print*,'min max Xk initial pour U',minloc(Xk),minval(Xk),maxloc(Xk),maxval(Xk),NsInt
          !CALL plot_vtk (Xk,'U_pb','U')
          stop
     end if
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
        F(is) = AireDsommet(is)*(CFDT*( Xk(is)-Uold(is) )/dt - rho1*hu(Cm(is))*fu(Xk(is),Em(is)) + beta1*Xk(is) + Tu(temps,Xk(is)) ) 
        !! derivee de F par rapport a Xk(is)
        call ajout(is,is, AireDSommet(is)*(CFDT/dt - rho1*hu(Cm(is))*deriveefu(Xk(is),Em(is)) + beta1+deriveeTu(temps,Xk(is))), A)
     END DO
     !
     DO jt = 1, Nbt
        ii = jt + NsInt
        F(ii) = CFDT*AireK(jt)*(( Xk(ii)-Uold(ii) )/dt - rho1*hu(Cm(ii))*fu(Xk(ii),Em(ii)) + beta1*Xk(ii) + Tu(temps,Xk(ii))) 
        !! derivee de F par rapport a Xk(ii)
        call ajout(ii,ii, AireK(jt)*(CFDT/dt - rho1*hu(Cm(ii))*deriveefu(Xk(ii),Em(ii)) + beta1 + deriveeTu(temps,Xk(ii))), A )
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
                 coef = uTKL(iseg) ; coeta = uetaSSe(iseg)
                 !-------------------------------
                 ! 1. contribution du triangle iK
                 !-------------------------------
                 dVKL = coeta * ( xiU(Xk(is)) - xiU(Xk(js)) )
                 dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                 F(iK) = F(iK) + coef*(AdegenU(Xk(iK)) - AdegenU(Xk(jL))) + &
                      & dVKLplus*( rmCroitU(Xk(iK)) + rmDeCroitU(Xk(jL)) ) + &
                      & dVKLmoins*( rmCroitU(Xk(jL)) + rmDeCroitU(Xk(iK)))
                 !! derivee de F(iK) par rapport a Xk(iK)
                 CoefAjout = coef*DerivAdegU(Xk(iK))+ dVKLplus*DerrmCroitU(Xk(iK)) + dVKLmoins*DerrmDeCroitU(Xk(iK))
                 call ajout(iK,iK,CoefAjout, A )
                 !! derivee de F(iK) par rapport a Xk(jL)
                 CoefAjout = -coef*DerivAdegU(Xk(jL))+ dVKLplus*DerrmDeCroitU(Xk(jL)) + dVKLmoins*DerrmCroitU(Xk(jL))
                 call ajout(iK,jL,CoefAjout, A )
                 !! derivee de F(iK) par rapport a Xk(is)
                 CoefAjout = transfer((dVKL>0),1)*coeta* rmU(Xk(is))*( rmCroitU(Xk(iK)) + rmDeCroitU(Xk(jL)) ) &
                      & + transfer((dVKL<0),1)*coeta* rmU(Xk(is))*( rmCroitU(Xk(jL)) + rmDeCroitU(Xk(iK)) )
                 call ajout(iK,is, CoefAjout, A)
                 !! derivee de F(iK) par rapport a Xk(js)
                 CoefAjout = -transfer((dVKL>0),1)*coeta* rmU(Xk(js))*( rmCroitU(Xk(iK)) + rmDeCroitU(Xk(jL)) ) &
                      & - transfer((dVKL<0),1)*coeta* rmU(Xk(js))*( rmCroitU(Xk(jL)) + rmDeCroitU(Xk(iK)) )
                 call ajout(iK,js, CoefAjout, A)
                 !-------------------------------
                 ! 2. contribution du triangle jL
                 !-------------------------------
                 dVKL = coeta * ( xiU(Xk(js)) - xiU(Xk(is)) )
                 dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                 F(jL) = F(jL) + coef*(AdegenU(Xk(jL)) - AdegenU(Xk(iK))) + &
                      & dVKLplus*( rmCroitU(Xk(jL)) + rmDeCroitU(Xk(iK)) ) + &
                      & dVKLmoins*( rmCroitU(Xk(iK)) + rmDeCroitU(Xk(jL)))
                 !! derivee de F(jL) par rapport a Xk(iK)
                 CoefAjout = -coef*DerivAdegU(Xk(iK))+ dVKLplus*DerrmDeCroitU(Xk(iK)) + dVKLmoins*DerrmCroitU(Xk(iK))
                 call ajout(jL,iK,CoefAjout, A )
                 !! derivee de F(jL) par rapport a Xk(jL)
                 CoefAjout = coef*DerivAdegU(Xk(jL))+ dVKLplus*DerrmCroitU(Xk(jL)) + dVKLmoins*DerrmDeCroitU(Xk(jL))
                 call ajout(jL,jL,CoefAjout, A )
                 !! derivee de F(jL) par rapport a Xk(is)
                 CoefAjout = -transfer((dVKL>0),1)*coeta* rmU(Xk(is))*( rmCroitU(Xk(jL)) + rmDeCroitU(Xk(iK)) ) &
                      & -transfer((dVKL<0),1)*coeta* rmU(Xk(is))*( rmCroitU(Xk(iK)) + rmDeCroitU(Xk(jL)) )
                 call ajout(jL,is, CoefAjout, A)
                 !! derivee de F(jL) par rapport a Xk(js)
                 CoefAjout = transfer((dVKL>0),1)*coeta* rmU(Xk(js))*( rmCroitU(Xk(jL)) + rmDeCroitU(Xk(iK)) ) &
                      & + transfer((dVKL<0),1)*coeta* rmU(Xk(js))*( rmCroitU(Xk(iK)) + rmDeCroitU(Xk(jL)) )
                 call ajout(jL,js, CoefAjout, A)
                 !-----------------------------
                 ! 3. contribution du sommet is
                 !-----------------------------
                 coef = uTKeLe(iseg)
                 dVKL = coeta * ( xiU(Xk(iK)) - xiU(Xk(jL)) )
                 dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                 F(is) = F(is) + coef*(AdegenU(Xk(is)) - AdegenU(Xk(js))) + &
                      & dVKLplus*( rmCroitU(Xk(is)) + rmDeCroitU(Xk(js)) ) + &
                      & dVKLmoins*( rmCroitU(Xk(js)) + rmDeCroitU(Xk(is)))
                 !! derivee de F(is) par rapport a Xk(is)
                 CoefAjout = coef*DerivAdegU(Xk(is))+ dVKLplus*DerrmCroitU(Xk(is)) + dVKLmoins*DerrmDeCroitU(Xk(is))
                 call ajout(is,is,CoefAjout, A )
                 !! derivee de F(is) par rapport a Xk(js)
                 CoefAjout = -coef*DerivAdegU(Xk(js))+ dVKLplus*DerrmDeCroitU(Xk(js)) + dVKLmoins*DerrmCroitU(Xk(js))
                 call ajout(is,js,CoefAjout, A )
                 !! derivee de F(is) par rapport a Xk(iK)
                 CoefAjout = transfer((dVKL>0),1)*coeta* rmU(Xk(iK))*( rmCroitU(Xk(is)) + rmDeCroitU(Xk(js)) ) &
                      & + transfer((dVKL<0),1)*coeta* rmU(Xk(iK))*( rmCroitU(Xk(js)) + rmDeCroitU(Xk(is)) )
                 call ajout(is,iK, CoefAjout, A)
                 !! derivee de F(is) par rapport a Xk(jL)
                 CoefAjout = -transfer((dVKL>0),1)*coeta* rmU(Xk(jL))*( rmCroitU(Xk(is)) + rmDeCroitU(Xk(js)) ) &
                      & - transfer((dVKL<0),1)*coeta* rmU(Xk(jL))*( rmCroitU(Xk(js)) + rmDeCroitU(Xk(is)) )
                 call ajout(is,jL, CoefAjout, A)
                 !-----------------------------
                 ! 4. contribution du sommet js
                 !-----------------------------
                 dVKL = coeta * ( xiU(Xk(jL)) - xiU(Xk(iK)) )
                 dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                 F(js) = F(js) + coef*(AdegenU(Xk(js)) - AdegenU(Xk(is))) + &
                      & dVKLplus*( rmCroitU(Xk(js)) + rmDeCroitU(Xk(is)) ) + &
                      & dVKLmoins*( rmCroitU(Xk(is)) + rmDeCroitU(Xk(js)))
                 !! derivee de F(js) par rapport a Xk(is)
                 CoefAjout = -coef*DerivAdegU(Xk(is))+ dVKLplus*DerrmDeCroitU(Xk(is)) + dVKLmoins*DerrmCroitU(Xk(is))
                 call ajout(js,is,CoefAjout, A )
                 !! derivee de F(js) par rapport a Xk(js)
                 CoefAjout = coef*DerivAdegU(Xk(js))+ dVKLplus*DerrmCroitU(Xk(js)) + dVKLmoins*DerrmDeCroitU(Xk(js))
                 call ajout(js,js,CoefAjout, A )
                 !! derivee de F(js) par rapport a Xk(iK)
                 CoefAjout = -transfer((dVKL>0),1)*coeta* rmU(Xk(iK))*( rmCroitU(Xk(js)) + rmDeCroitU(Xk(is)) ) &
                      & - transfer((dVKL<0),1)*coeta* rmU(Xk(iK))*( rmCroitU(Xk(is)) + rmDeCroitU(Xk(js)) )
                 call ajout(js,iK, CoefAjout, A)
                 !! derivee de F(js) par rapport a Xk(jL)
                 CoefAjout = transfer((dVKL>0),1)*coeta* rmU(Xk(jL))*( rmCroitU(Xk(js)) + rmDeCroitU(Xk(is)) ) &
                      & + transfer((dVKL<0),1)*coeta* rmU(Xk(jL))*( rmCroitU(Xk(is)) + rmDeCroitU(Xk(js)) )
                 call ajout(js,jL, CoefAjout, A)
                 !
              ELSE
                 Select case ( Ntyps(js) ) 
                 case (dirichlet)
                    Ubord = Gb(js)
                    coef = uTKL(iseg) ; coeta = uetaSSe(iseg)
                    !-------------------------------
                    ! 1. contribution du triangle iK
                    !-------------------------------
                    dVKL = coeta * ( xiU(Xk(is)) - xiU(Ubord) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    F(iK) = F(iK) + coef*(AdegenU(Xk(iK)) - AdegenU(Xk(jL))) + &
                         & dVKLplus*( rmCroitU(Xk(iK)) + rmDeCroitU(Xk(jL)) ) + &
                         & dVKLmoins*( rmCroitU(Xk(jL)) + rmDeCroitU(Xk(iK)))
                    !! derivee de F(iK) par rapport a Xk(iK)
                    CoefAjout = coef*DerivAdegU(Xk(iK))+ dVKLplus*DerrmCroitU(Xk(iK)) + dVKLmoins*DerrmDeCroitU(Xk(iK))
                    call ajout(iK,iK,CoefAjout, A )
                    !! derivee de F(iK) par rapport a Xk(jL)
                    CoefAjout = -coef*DerivAdegU(Xk(jL))+ dVKLplus*DerrmDeCroitU(Xk(jL)) + dVKLmoins*DerrmCroitU(Xk(jL))
                    call ajout(iK,jL,CoefAjout, A )
                    !! derivee de F(iK) par rapport a Xk(is)
                    CoefAjout = transfer((dVKL>0),1)*coeta* rmU(Xk(is))*( rmCroitU(Xk(iK)) + rmDeCroitU(Xk(jL)) ) &
                         & + transfer((dVKL<0),1)*coeta* rmU(Xk(is))*( rmCroitU(Xk(jL)) + rmDeCroitU(Xk(iK)) )
                    call ajout(iK,is, CoefAjout, A)
                    !-------------------------------
                    ! 2. contribution du triangle jL
                    !-------------------------------
                    dVKL = coeta * ( xiU(Ubord) - xiU(Xk(is)) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    F(jL) = F(jL) + coef*(AdegenU(Xk(jL)) - AdegenU(Xk(iK))) + &
                         & dVKLplus*( rmCroitU(Xk(jL)) + rmDeCroitU(Xk(iK)) ) + &
                         & dVKLmoins*( rmCroitU(Xk(iK)) + rmDeCroitU(Xk(jL)))
                    !! derivee de F(jL) par rapport a Xk(iK)
                    CoefAjout = -coef*DerivAdegU(Xk(iK))+ dVKLplus*DerrmDeCroitU(Xk(iK)) + dVKLmoins*DerrmCroitU(Xk(iK))
                    call ajout(jL,iK,CoefAjout, A )
                    !! derivee de F(jL) par rapport a Xk(jL)
                    CoefAjout = coef*DerivAdegU(Xk(jL))+ dVKLplus*DerrmCroitU(Xk(jL)) + dVKLmoins*DerrmDeCroitU(Xk(jL))
                    call ajout(jL,jL,CoefAjout, A )
                    !! derivee de F(jL) par rapport a Xk(is)
                    CoefAjout = -transfer((dVKL>0),1)*coeta* rmU(Xk(is))*( rmCroitU(Xk(jL)) + rmDeCroitU(Xk(iK)) ) &
                         & - transfer((dVKL<0),1)*coeta* rmU(Xk(is))*( rmCroitU(Xk(iK)) + rmDeCroitU(Xk(jL)) )
                    call ajout(jL,is, CoefAjout, A)
                    !-----------------------------
                    ! 3. contribution du sommet is
                    !-----------------------------
                    coef = uTKeLe(iseg)
                    dVKL = coeta * ( xiU(Xk(iK)) - xiU(Xk(jL)) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    F(is) = F(is) + coef*(AdegenU(Xk(is)) - AdegenU(Ubord)) + &
                         & dVKLplus*( rmCroitU(Xk(is)) + rmDeCroitU(Ubord) ) + &
                         & dVKLmoins*( rmCroitU(Ubord) + rmDeCroitU(Xk(is)))
                    !! derivee de F(is) par rapport a Xk(is)
                    CoefAjout = coef*DerivAdegU(Xk(is))+ dVKLplus*DerrmCroitU(Xk(is)) + dVKLmoins*DerrmDeCroitU(Xk(is))
                    call ajout(is,is,CoefAjout, A )
                    !! derivee de F(is) par rapport a Xk(iK)
                    CoefAjout = transfer((dVKL>0),1)*coeta* rmU(Xk(iK))*( rmCroitU(Xk(is)) + rmDeCroitU(Ubord) ) &
                         & + transfer((dVKL<0),1)*coeta* rmU(Xk(iK))*( rmCroitU(Ubord) + rmDeCroitU(Xk(is)) )
                    call ajout(is,iK, CoefAjout, A)
                    !! derivee de F(is) par rapport a Xk(jL)
                    CoefAjout = -transfer((dVKL>0),1)*coeta* rmU(Xk(jL))*( rmCroitU(Xk(is)) + rmDeCroitU(Ubord) ) &
                         & - transfer((dVKL<0),1)*coeta* rmU(Xk(jL))*( rmCroitU(Ubord) + rmDeCroitU(Xk(is)) )
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
                    coef = uTKL(iseg) ; coeta = uetaSSe(iseg)
                    !-------------------------------
                    ! 1. contribution du triangle iK
                    !-------------------------------
                    dVKL = coeta * ( xiU(Ubord) - xiU(Xk(js)) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    F(iK) = F(iK) + coef*(AdegenU(Xk(iK)) - AdegenU(Xk(jL))) + &
                         & dVKLplus*( rmCroitU(Xk(iK)) + rmDeCroitU(Xk(jL)) ) + &
                         & dVKLmoins*( rmCroitU(Xk(jL)) + rmDeCroitU(Xk(iK)))
                    !! derivee de F(iK) par rapport a Xk(iK)
                    CoefAjout = coef*DerivAdegU(Xk(iK))+ dVKLplus*DerrmCroitU(Xk(iK)) + dVKLmoins*DerrmDeCroitU(Xk(iK))
                    call ajout(iK,iK,CoefAjout, A )
                    !! derivee de F(iK) par rapport a Xk(jL)
                    CoefAjout = -coef*DerivAdegU(Xk(jL))+ dVKLplus*DerrmDeCroitU(Xk(jL)) + dVKLmoins*DerrmCroitU(Xk(jL))
                    call ajout(iK,jL,CoefAjout, A )
                    !! derivee de F(iK) par rapport a Xk(js)
                    CoefAjout = -transfer((dVKL>0),1)*coeta* rmU(Xk(js))*( rmCroitU(Xk(iK)) + rmDeCroitU(Xk(jL)) )&
                         & - transfer((dVKL<0),1)*coeta* rmU(Xk(js))*( rmCroitU(Xk(jL)) + rmDeCroitU(Xk(iK)) )
                    call ajout(iK,js, CoefAjout, A)
                    !-------------------------------
                    ! 2. contribution du triangle jL
                    !-------------------------------
                    dVKL = coeta * ( xiU(Xk(js)) - xiU(Ubord) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    F(jL) = F(jL) + coef*(AdegenU(Xk(jL)) - AdegenU(Xk(iK))) + &
                         & dVKLplus*( rmCroitU(Xk(jL)) + rmDeCroitU(Xk(iK)) ) + &
                         & dVKLmoins*( rmCroitU(Xk(iK)) + rmDeCroitU(Xk(jL)))
                    !! derivee de F(jL) par rapport a Xk(iK)
                    CoefAjout = -coef*DerivAdegU(Xk(iK))+ dVKLplus*DerrmDeCroitU(Xk(iK)) + dVKLmoins*DerrmCroitU(Xk(iK))
                    call ajout(jL,iK,CoefAjout, A )
                    !! derivee de F(jL) par rapport a Xk(jL)
                    CoefAjout = coef*DerivAdegU(Xk(jL))+ dVKLplus*DerrmCroitU(Xk(jL)) + dVKLmoins*DerrmDeCroitU(Xk(jL))
                    call ajout(jL,jL,CoefAjout, A )
                    !! derivee de F(jL) par rapport a Xk(js)
                    CoefAjout = transfer((dVKL>0),1)*coeta* rmU(Xk(js))*( rmCroitU(Xk(jL)) + rmDeCroitU(Xk(iK)) ) &
                         & + transfer((dVKL<0),1)*coeta* rmU(Xk(js))*( rmCroitU(Xk(iK)) + rmDeCroitU(Xk(jL)) )
                    call ajout(jL,js, CoefAjout, A)
                    !-----------------------------
                    ! 3. contribution du sommet js
                    !-----------------------------
                    coef = uTKeLe(iseg)
                    dVKL = coeta * ( xiU(Xk(jL)) - xiU(Xk(iK)) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    F(js) = F(js) + coef*(AdegenU(Xk(js)) - AdegenU(Ubord)) + &
                         & dVKLplus*( rmCroitU(Xk(js)) + rmDeCroitU(Ubord) ) + &
                         & dVKLmoins*( rmCroitU(Ubord) + rmDeCroitU(Xk(js)))
                    !! derivee de F(js) par rapport a Xk(js)
                    CoefAjout = coef*DerivAdegU(Xk(js))+ dVKLplus*DerrmCroitU(Xk(js)) + dVKLmoins*DerrmDeCroitU(Xk(js))
                    call ajout(js,js,CoefAjout, A )
                    !! derivee de F(js) par rapport a Xk(iK)
                    CoefAjout = -transfer((dVKL>0),1)*coeta* rmU(Xk(iK))*( rmCroitU(Xk(js)) + rmDeCroitU(Ubord) ) &
                         & - transfer((dVKL<0),1)*coeta* rmU(Xk(iK))*( rmCroitU(Ubord) + rmDeCroitU(Xk(js)) )
                    call ajout(js,iK, CoefAjout, A)
                    !! derivee de F(js) par rapport a Xk(jL)
                    CoefAjout = transfer((dVKL>0),1)*coeta* rmU(Xk(jL))*( rmCroitU(Xk(js)) + rmDeCroitU(Ubord) ) &
                         & + transfer((dVKL<0),1)*coeta* rmU(Xk(jL))*( rmCroitU(Ubord) + rmDeCroitU(Xk(js)) )
                    call ajout(js,jL, CoefAjout, A)
                 ELSE ! js se situe encore au bord (un triangle entier sur le bord)
                    Uibord= Gb(is) ; Ujbord= Gb(js)
                    coef = uTKL(iseg) ; coeta = uetaSSe(iseg)
                    !-------------------------------
                    ! 1. contribution du triangle iK
                    !-------------------------------
                    dVKL = coeta * ( xiU(Uibord) - xiU(Ujbord) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    F(iK) = F(iK) + coef*(AdegenU(Xk(iK)) - AdegenU(Xk(jL))) + &
                         & dVKLplus*( rmCroitU(Xk(iK)) + rmDeCroitU(Xk(jL)) ) + &
                         & dVKLmoins*( rmCroitU(Xk(jL)) + rmDeCroitU(Xk(iK)))
                    !! derivee de F(iK) par rapport a Xk(iK)
                    CoefAjout = coef*DerivAdegU(Xk(iK))+ dVKLplus*DerrmCroitU(Xk(iK)) + dVKLmoins*DerrmDeCroitU(Xk(iK))
                    call ajout(iK,iK,CoefAjout, A )
                    !! derivee de F(iK) par rapport a Xk(jL)
                    CoefAjout = -coef*DerivAdegU(Xk(jL))+ dVKLplus*DerrmDeCroitU(Xk(jL)) + dVKLmoins*DerrmCroitU(Xk(jL))
                    call ajout(iK,jL,CoefAjout, A )
                    !-------------------------------
                    ! 2. contribution du triangle jL
                    !-------------------------------
                    dVKL = coeta * ( xiU(Ujbord) - xiU(Uibord) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    F(jL) = F(jL) + coef*(AdegenU(Xk(jL)) - AdegenU(Xk(iK))) + &
                         & dVKLplus*( rmCroitU(Xk(jL)) + rmDeCroitU(Xk(iK)) ) + &
                         & dVKLmoins*( rmCroitU(Xk(iK)) + rmDeCroitU(Xk(jL)))
                    !! derivee de F(jL) par rapport a Xk(iK)
                    CoefAjout = -coef*DerivAdegU(Xk(iK))+ dVKLplus*DerrmDeCroitU(Xk(iK)) + dVKLmoins*DerrmCroitU(Xk(iK))
                    call ajout(jL,iK,CoefAjout, A )
                    !! derivee de F(jL) par rapport a Xk(jL)
                    CoefAjout = coef*DerivAdegU(Xk(jL))+ dVKLplus*DerrmCroitU(Xk(jL)) + dVKLmoins*DerrmDeCroitU(Xk(jL))
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
           coef = uTKL(iseg) ; coeta = uetaSSe(iseg)
           !-------------------------------
           ! 1. contribution du triangle iK
           !-------------------------------
           dVKL = coeta * ( xiU(Uibord) - xiU(Ujbord) )
           dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
           F(iK) = F(iK) + coef*(AdegenU(Xk(iK)) - AdegenU(Ubord)) + &
                & dVKLplus*( rmCroitU(Xk(iK)) + rmDeCroitU(Ubord) ) + &
                & dVKLmoins*( rmCroitU(Ubord) + rmDeCroitU(Xk(iK)))
           !! derivee de F(iK) par rapport a Xk(iK)
           CoefAjout = coef*DerivAdegU(Xk(iK))+ dVKLplus*DerrmCroitU(Xk(iK)) + dVKLmoins*DerrmDeCroitU(Xk(iK))
           call ajout(iK,iK,CoefAjout, A )
           !!
        case (neumann)
           !!
           is = NuSeg(1,iseg); js = NuSeg(2,iseg)
           iK = NumTVoisSeg(1,iseg) + NsInt
           coef = uTKL(iseg) ; coeta = uetaSSe(iseg) ; coefe = uTKeLe(iseg)
           !-------------------------------
           ! 1. contribution du triangle iK
           !-------------------------------
           ! On fait rien
           !-----------------------------
           ! 2. contribution du sommet is
           !----------------------------- 
           dVKL = ((coeta**2)/coef) * ( xiU(Xk(js)) - xiU(Xk(is)) )
           dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
           F(is) = F(is) + coefe*(AdegenU(Xk(is)) - AdegenU(Xk(js))) + &
                & dVKLplus*( rmCroitU(Xk(is)) + rmDeCroitU(Xk(js)) ) + &
                & dVKLmoins*( rmCroitU(Xk(js)) + rmDeCroitU(Xk(is)))
           !! derivee de F(is) par rapport a Xk(is)
           CoefAjout = coefe*DerivAdegU(Xk(is))+ dVKLplus*DerrmCroitU(Xk(is)) + dVKLmoins*DerrmDeCroitU(Xk(is)) &
                & - transfer((dVKL>0),1)*((coeta**2)/coef)* rmU(Xk(is))*( rmCroitU(Xk(is)) + rmDeCroitU(Xk(js)))&
                & - transfer((dVKL<0),1)*((coeta**2)/coef)* rmU(Xk(is))*( rmCroitU(Xk(js)) + rmDeCroitU(Xk(is)))
           call ajout(is,is,CoefAjout, A )
           !! derivee de F(is) par rapport a Xk(js)
           CoefAjout = -coefe*DerivAdegU(Xk(js))+ dVKLplus*DerrmDeCroitU(Xk(js)) + dVKLmoins*DerrmCroitU(Xk(js))&
                & + transfer((dVKL>0),1)*((coeta**2)/coef)* rmU(Xk(js))*( rmCroitU(Xk(is)) + rmDeCroitU(Xk(js)))&
                & + transfer((dVKL<0),1)*((coeta**2)/coef)* rmU(Xk(js))*( rmCroitU(Xk(js)) + rmDeCroitU(Xk(is)))
           call ajout(is,js,CoefAjout, A )
           !-----------------------------
           ! 3. contribution du sommet js
           !-----------------------------
           dVKL = ((coeta**2)/coef) * ( xiU(Xk(is)) - xiU(Xk(js)) )
           dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
           F(js) = F(js) + coefe*(AdegenU(Xk(js)) - AdegenU(Xk(is))) + &
                & dVKLplus*( rmCroitU(Xk(js)) + rmDeCroitU(Xk(is)) ) + &
                & dVKLmoins*( rmCroitU(Xk(is)) + rmDeCroitU(Xk(js)))
           !! derivee de F(js) par rapport a Xk(is)
           CoefAjout = -coefe*DerivAdegU(Xk(is))+ dVKLplus*DerrmDeCroitU(Xk(is)) + dVKLmoins*DerrmCroitU(Xk(is))&
                & + transfer((dVKL>0),1)*((coeta**2)/coef)* rmU(Xk(is))*( rmCroitU(Xk(js)) + rmDeCroitU(Xk(is)))&
                & + transfer((dVKL<0),1)*((coeta**2)/coef)* rmU(Xk(is))*( rmCroitU(Xk(is)) + rmDeCroitU(Xk(js)))
           call ajout(js,is,CoefAjout, A )
           !! derivee de F(js) par rapport a Xk(js)
           CoefAjout = coefe*DerivAdegU(Xk(js))+ dVKLplus*DerrmCroitU(Xk(js)) + dVKLmoins*DerrmDeCroitU(Xk(js)) &
                & - transfer((dVKL>0),1)*((coeta**2)/coef)* rmU(Xk(js))*( rmCroitU(Xk(js)) + rmDeCroitU(Xk(is)))&
                & - transfer((dVKL<0),1)*((coeta**2)/coef)* rmU(Xk(js))*( rmCroitU(Xk(is)) + rmDeCroitU(Xk(js)))
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
                 coef = uTKL(iseg) ; coeta = uetaSSe(iseg)
                 RDVT = (Cm(is) + Cm(js) + Cm(iK) + Cm(jL))/4
                 !-------------------------------
                 ! 1. contribution du triangle iK
                 !-------------------------------
                 dVKL = coef*(ln(Cm(jL)) - ln(Cm(iK))) + &
                      & coeta*(ln(Cm(js)) - ln(Cm(is)))
                 dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                 FKL = dVKLplus*( MuCroitU(Xk(iK)) + MuDeCroitU(Xk(jL)) ) + &
                      & dVKLmoins*( MuCroitU(Xk(jL)) + MuDeCroitU(Xk(iK)) )
                 F(iK) = F(iK) + RDVT* FKL
                 !! derivee de F(iK) par rapport a Xk(iK)
                 CoefAjout = RDVT*(dVKLplus*DerivMuCroitU(Xk(iK)) + dVKLmoins*DerivMuDeCroitU(Xk(iK)))
                 call ajout(iK,iK,CoefAjout, A )
                 !! derivee de F(iK) par rapport a Xk(jL)
                 CoefAjout = RDVT*(dVKLplus*DerivMuDeCroitU(Xk(jL)) + dVKLmoins*DerivMuCroitU(Xk(jL)))
                 call ajout(iK,jL,CoefAjout, A )
                 !-------------------------------
                 ! 2. contribution du triangle jL
                 !-------------------------------
                 dVKL = coef*(ln(Cm(iK)) - ln(Cm(jL))) + &
                      & coeta * ( ln(Cm(is)) - ln(Cm(js)) )
                 dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                 FKL = dVKLplus*( MuCroitU(Xk(jL)) + MuDeCroitU(Xk(iK)) ) + &
                      & dVKLmoins*( MuCroitU(Xk(iK)) + MuDeCroitU(Xk(jL)) )
                 F(jL) = F(jL) + RDVT* FKL
                 !! derivee de F(jL) par rapport a Xk(jL)
                 CoefAjout = RDVT*(dVKLplus*DerivMuCroitU(Xk(jL)) + dVKLmoins*DerivMuDeCroitU(Xk(jL)))
                 call ajout(jL,jL,CoefAjout, A )
                 !! derivee de F(jL) par rapport a Xk(iK)
                 CoefAjout = RDVT*(dVKLplus*DerivMuDeCroitU(Xk(iK)) + dVKLmoins*DerivMuCroitU(Xk(iK)))
                 call ajout(jL,iK,CoefAjout, A )
                 !-----------------------------
                 ! 3. contribution du sommet is
                 !-----------------------------
                 coef = uTKeLe(iseg)
                 dVKL = coef*(ln(Cm(js)) - ln(Cm(is))) + &
                      & coeta * ( ln(Cm(jL)) - ln(Cm(iK)) )
                 dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                 FKL = dVKLplus*( MuCroitU(Xk(is)) + MuDeCroitU(Xk(js)) ) + &
                      & dVKLmoins*( MuCroitU(Xk(js)) + MuDeCroitU(Xk(is)) )
                 F(is) = F(is) + RDVT* FKL
                 !! derivee de F(is) par rapport a Xk(is)
                 CoefAjout = RDVT*(dVKLplus*DerivMuCroitU(Xk(is)) + dVKLmoins*DerivMuDeCroitU(Xk(is)))
                 call ajout(is,is,CoefAjout, A )
                 !! derivee de F(is) par rapport a Xk(js)
                 CoefAjout = RDVT*(dVKLplus*DerivMuDeCroitU(Xk(js)) + dVKLmoins*DerivMuCroitU(Xk(js)))
                 call ajout(is,js,CoefAjout, A )
                 !----------------------- -----
                 ! 4. contribution du sommet js
                 !-----------------------------
                 dVKL = coef*(ln(Cm(is)) - ln(Cm(js))) + &
                      & coeta * ( ln(Cm(iK)) - ln(Cm(jL)) )
                 dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                 FKL = dVKLplus*( MuCroitU(Xk(js)) + MuDeCroitU(Xk(is)) ) + &
                      & dVKLmoins*( MuCroitU(Xk(is)) + MuDeCroitU(Xk(js)) )
                 F(js) = F(js) + RDVT* FKL
                 !! derivee de F(js) par rapport a Xk(js)
                 CoefAjout = RDVT*(dVKLplus*DerivMuCroitU(Xk(js)) + dVKLmoins*DerivMuDeCroitU(Xk(js)))
                 call ajout(js,js,CoefAjout, A )
                 !! derivee de F(js) par rapport a Xk(is)
                 CoefAjout = RDVT*(dVKLplus*DerivMuDeCroitU(Xk(is)) + dVKLmoins*DerivMuCroitU(Xk(is)))
                 call ajout(js,is,CoefAjout, A )
              ELSE
                 Select case ( Ntyps(js) ) 
                 case (dirichlet)
                    Ubord = Gb(js)
                    coef = uTKL(iseg) ; coeta = uetaSSe(iseg)
                    RDVT = (Cm(is) + Ubord + Cm(iK) + Cm(jL))/4
                    !-------------------------------
                    ! 1. contribution du triangle iK
                    !-------------------------------
                    dVKL = coef*(ln(Cm(jL)) - ln(Cm(iK))) + &
                         & coeta * ( ln(Ubord) - ln(Cm(is)) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    FKL = dVKLplus*( MuCroitU(Xk(iK)) + MuDeCroitU(Xk(jL)) ) + &
                         & dVKLmoins*( MuCroitU(Xk(jL)) + MuDeCroitU(Xk(iK)) )
                    F(iK) = F(iK) + RDVT* FKL
                    !! derivee de F(iK) par rapport a Xk(iK)
                    CoefAjout = RDVT*(dVKLplus*DerivMuCroitU(Xk(iK)) + dVKLmoins*DerivMuDeCroitU(Xk(iK)))
                    call ajout(iK,iK,CoefAjout, A )
                    !! derivee de F(iK) par rapport a Xk(jL)
                    CoefAjout = RDVT*(dVKLplus*DerivMuDeCroitU(Xk(jL)) + dVKLmoins*DerivMuCroitU(Xk(jL)))
                    call ajout(iK,jL,CoefAjout, A )
                    !-------------------------------
                    ! 2. contribution du triangle jL
                    !-------------------------------
                    dVKL = coef*(ln(Cm(iK)) - ln(Cm(jL))) + &
                         & coeta * ( ln(Cm(is)) - ln(Ubord) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    FKL = dVKLplus*( MuCroitU(Xk(jL)) + MuDeCroitU(Xk(iK)) ) + &
                         & dVKLmoins*( MuCroitU(Xk(iK)) + MuDeCroitU(Xk(jL)) )
                    F(jL) = F(jL) + RDVT* FKL
                    !! derivee de F(jL) par rapport a Xk(jL)
                    CoefAjout = RDVT*(dVKLplus*DerivMuCroitU(Xk(jL)) + dVKLmoins*DerivMuDeCroitU(Xk(jL)))
                    call ajout(jL,jL,CoefAjout, A )
                    !! derivee de F(jL) par rapport a Xk(iK)
                    CoefAjout = RDVT*(dVKLplus*DerivMuDeCroitU(Xk(iK)) + dVKLmoins*DerivMuCroitU(Xk(iK)))
                    call ajout(jL,iK,CoefAjout, A )
                    !-----------------------------
                    ! 3. contribution du sommet is
                    !-----------------------------
                    coef = uTKeLe(iseg)
                    dVKL = coef*(ln(Ubord) - ln(Cm(is))) + &
                         & coeta * ( ln(Cm(jL)) - ln(Cm(iK)) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    FKL = dVKLplus*( MuCroitU(Xk(is)) + MuDeCroitU(Ubord) ) + &
                         & dVKLmoins*( MuCroitU(Ubord) + MuDeCroitU(Xk(is)) )
                    F(is) = F(is) + RDVT* FKL
                    !! derivee de F(is) par rapport a Xk(is)
                    CoefAjout = RDVT*(dVKLplus*DerivMuCroitU(Xk(is)) + dVKLmoins*DerivMuDeCroitU(Xk(is)))
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
                    coef = uTKL(iseg) ; coeta = uetaSSe(iseg)
                    RDVT = (Ubord + Cm(js) + Cm(iK) + Cm(jL))/4
                    !-------------------------------
                    ! 1. contribution du triangle iK
                    !-------------------------------
                    dVKL = coef*(ln(Cm(jL)) - ln(Cm(iK))) + &
                         & coeta * ( ln(Cm(js)) - ln(Ubord) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    FKL = dVKLplus*( MuCroitU(Xk(iK)) + MuDeCroitU(Xk(jL)) ) + &
                         & dVKLmoins*( MuCroitU(Xk(jL)) + MuDeCroitU(Xk(iK)) )
                    F(iK) = F(iK) + RDVT* FKL
                    !! derivee de F(iK) par rapport a Xk(iK)
                    CoefAjout = RDVT*(dVKLplus*DerivMuCroitU(Xk(iK)) + dVKLmoins*DerivMuDeCroitU(Xk(iK)))
                    call ajout(iK,iK,CoefAjout, A )
                    !! derivee de F(iK) par rapport a Xk(jL)
                    CoefAjout = RDVT*(dVKLplus*DerivMuDeCroitU(Xk(jL)) + dVKLmoins*DerivMuCroitU(Xk(jL)))
                    call ajout(iK,jL,CoefAjout, A )
                    !-------------------------------
                    ! 2. contribution du triangle jL
                    !-------------------------------
                    dVKL = coef*(ln(Cm(iK)) - ln(Cm(jL))) + &
                         & coeta * ( ln(Ubord) - ln(Cm(js)) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    FKL = dVKLplus*( MuCroitU(Xk(jL)) + MuDeCroitU(Xk(iK)) ) + &
                         & dVKLmoins*( MuCroitU(Xk(iK)) + MuDeCroitU(Xk(jL)) )
                    F(jL) = F(jL) + RDVT* FKL
                    !! derivee de F(jL) par rapport a Xk(jL)
                    CoefAjout = RDVT*(dVKLplus*DerivMuCroitU(Xk(jL)) + dVKLmoins*DerivMuDeCroitU(Xk(jL)))
                    call ajout(jL,jL,CoefAjout, A )
                    !! derivee de F(jL) par rapport a Xk(iK)
                    CoefAjout = RDVT*(dVKLplus*DerivMuDeCroitU(Xk(iK)) + dVKLmoins*DerivMuCroitU(Xk(iK)))
                    call ajout(jL,iK,CoefAjout, A )
                    !-----------------------------
                    ! 3. contribution du sommet js
                    !-----------------------------
                    coef = uTKeLe(iseg)
                    dVKL = coef*(ln(Ubord) - ln(Cm(js))) + &
                         & coeta * ( ln(Cm(iK)) - ln(Cm(jL)) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    FKL = dVKLplus*( MuCroitU(Xk(js)) + MuDeCroitU(Ubord) ) + &
                         & dVKLmoins*( MuCroitU(Ubord) + MuDeCroitU(Xk(js)) )
                    F(js) = F(js) + RDVT* FKL
                    !! derivee de F(js) par rapport a Xk(js)
                    CoefAjout = RDVT*(dVKLplus*DerivMuCroitU(Xk(js)) + dVKLmoins*DerivMuDeCroitU(Xk(js)))
                    call ajout(js,js,CoefAjout, A )
                 ELSE ! js se situe encore au bord (un triangle entier sur le bord)
                    Uibord= Gb(is) ; Ujbord= Gb(js)
                    coef = uTKL(iseg) ; coeta = uetaSSe(iseg)
                    RDVT = (Uibord + Ujbord + Cm(iK) + Cm(jL))/4
                    !-------------------------------
                    ! 1. contribution du triangle iK
                    !-------------------------------
                    dVKL = coef*(ln(Cm(jL)) - ln(Cm(iK))) + &
                         & coeta * ( ln(Ujbord) - ln(Uibord) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    FKL = dVKLplus*( MuCroitU(Xk(iK)) + MuDeCroitU(Xk(jL)) ) + &
                         & dVKLmoins*( MuCroitU(Xk(jL)) + MuDeCroitU(Xk(iK)) )
                    F(iK) = F(iK) + RDVT* FKL
                    !! derivee de F(iK) par rapport a Xk(iK)
                    CoefAjout = RDVT*(dVKLplus*DerivMuCroitU(Xk(iK)) + dVKLmoins*DerivMuDeCroitU(Xk(iK)))
                    call ajout(iK,iK,CoefAjout, A )
                    !! derivee de F(iK) par rapport a Xk(jL)
                    CoefAjout = RDVT*(dVKLplus*DerivMuDeCroitU(Xk(jL)) + dVKLmoins*DerivMuCroitU(Xk(jL)))
                    call ajout(iK,jL,CoefAjout, A )
                    !-------------------------------
                    ! 2. contribution du triangle jL
                    !-------------------------------
                    dVKL = coef*(ln(Cm(iK)) - ln(Cm(jL))) + &
                         & coeta * ( ln(Uibord) - ln(Ujbord) )
                    dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
                    FKL = dVKLplus*( MuCroitU(Xk(jL)) + MuDeCroitU(Xk(iK)) ) + &
                         & dVKLmoins*( MuCroitU(Xk(iK)) + MuDeCroitU(Xk(jL)) )
                    F(jL) = F(jL) + RDVT* FKL
                    !! derivee de F(jL) par rapport a Xk(jL)
                    CoefAjout = RDVT*(dVKLplus*DerivMuCroitU(Xk(jL)) + dVKLmoins*DerivMuDeCroitU(Xk(jL)))
                    call ajout(jL,jL,CoefAjout, A )
                    !! derivee de F(jL) par rapport a Xk(iK)
                    CoefAjout = RDVT*(dVKLplus*DerivMuDeCroitU(Xk(iK)) + dVKLmoins*DerivMuCroitU(Xk(iK)))
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
           coef = uTKL(iseg) ; coeta = uetaSSe(iseg)
           RDVT = (Uibord + Ujbord + Cm(iK) + Ubord)/4
           !-------------------------------
           ! 1. contribution du triangle iK
           !-------------------------------
           dVKL = coef*(ln(Ubord) - ln(Cm(iK))) + &
                & coeta * ( ln(Ujbord) - ln(Uibord) )
           dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
           FKL = dVKLplus*( MuCroitU(Xk(iK)) + MuDeCroitU(Ubord) ) + &
                & dVKLmoins*( MuCroitU(Ubord) + MuDeCroitU(Xk(iK)) )
           F(iK) = F(iK) + RDVT* FKL
           !! derivee de F(iK) par rapport a Xk(iK)
           CoefAjout = RDVT*(dVKLplus*DerivMuCroitU(Xk(iK)) + dVKLmoins*DerivMuDeCroitU(Xk(iK)))
           call ajout(iK,iK,CoefAjout, A )
           !!
        case (neumann)
           !!
           is = NuSeg(1,iseg); js = NuSeg(2,iseg)
           iK = NumTVoisSeg(1,iseg) + NsInt
           coef = uTKL(iseg) ; coeta = uetaSSe(iseg) ; coefe = uTKeLe(iseg)
           RDVT = (3*(Cm(is) + Cm(js))/2 + Cm(iK))/4
           !-------------------------------
           ! 1. contribution du triangle iK
           !-------------------------------
           ! On fait rien
           !-----------------------------
           ! 2. contribution du sommet is
           !-----------------------------
           coefe = (coefe - (coeta**2)/coef)
           dVKL = coefe*(ln(Cm(js)) - ln(Cm(is)))
           dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
           FKL = dVKLplus*( MuCroitU(Xk(is)) + MuDeCroitU(Xk(js)) ) + &
                & dVKLmoins*( MuCroitU(Xk(js)) + MuDeCroitU(Xk(is)) )
           F(is) = F(is) + RDVT* FKL
           !! derivee de F(is) par rapport a Xk(is)
           CoefAjout = RDVT*(dVKLplus*DerivMuCroitU(Xk(is)) + dVKLmoins*DerivMuDeCroitU(Xk(is)))
           call ajout(is,is,CoefAjout, A )
           !! derivee de F(is) par rapport a Xk(js)
           CoefAjout = RDVT*(dVKLplus*DerivMuDeCroitU(Xk(js)) + dVKLmoins*DerivMuCroitU(Xk(js)))
           call ajout(is,js,CoefAjout, A )
           !-----------------------------
           ! 3. contribution du sommet js
           !-----------------------------
           dVKL = coefe*(ln(Cm(is)) - ln(Cm(js)))
           dVKLplus = max(dVKL, 0.D0) ; dVKLmoins = min(dVKL, 0.D0)
           FKL = dVKLplus*( MuCroitU(Xk(js)) + MuDeCroitU(Xk(is))) + &
                & dVKLmoins*( MuCroitU(Xk(is)) + MuDeCroitU(Xk(js)))
           F(js) = F(js) + RDVT* FKL
           !! derivee de F(js) par rapport a Xk(js)
           CoefAjout = RDVT*(dVKLplus*DerivMuCroitU(Xk(js)) + dVKLmoins*DerivMuDeCroitU(Xk(js)))
           call ajout(js,js,CoefAjout, A )
           !! derivee de F(js) par rapport a Xk(is)
           CoefAjout = RDVT*(dVKLplus*DerivMuDeCroitU(Xk(is)) + dVKLmoins*DerivMuCroitU(Xk(is)))
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
     !write(*,*)'F =',F
     print*,'kiter = ', kiter,'erreur NEWTON sqrt(sum(Yk*Yk)) pour U',sqrt(sum(Yk*Yk))
     !print*,'max Xk', maxval(Xk), 'min Xk', minval(Xk)
     !print*,'max F', maxval(F), 'min F', minval(F)
     If (sqrt(dot_product(Yk,Yk)) <TolerenceNewton) exit
     Xk=Xk+Yk
     
     !stop
     !! impression de A etde F
     !     DO i = 1, size(A%IndPL)-1
     !           print*,'iseg = ', i, 'sum line = ', SUM(A%Tmat(A%IndPL(i): A%IndPL(i+1)-1))
     !    END DO
  END Do
  U=Xk+Yk
  print*,'fin newtonu'
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






