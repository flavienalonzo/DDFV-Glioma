!            **************************
!            **  FUNCTION transmis  **
!            **************************
!******************************************************************************
!     *-!!!!!!!!!@@@@@!!@@@@@!!@@@@@!!@@@@@!!!!!!!!!!!-*
!     *-!!!!!!!!!!!!!@!!@!!!@!!!!!!@!!!!!!@!!!!!!!!!!!-*
!     *-!!!!!!!!!@@@@@!!@!!!@!!@@@@@!!@@@@@!!!!!!!!!!!-*
!     *-!!!!!!!!!@!!!!!!@!!!@!!@@!!!!!!!!!@!!!!!!!!!!!-*
!     *-!!!!!!!!!@@@@@!!@@@@@!!@@@@@!!@@@@@!!!!!!!!!!!-*
!------------------------------------------------------------
! calcul sur chaque triangle le tenseur de diffusion 
! SxxK, SxyK, SyyK
!
!******************************************************************************
SUBROUTINE transmis(choixaniu,choixanic,choixanie,choixaniv)
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
    Integer, intent(in)               :: choixaniu,choixanic,choixanie,choixaniv
    INTEGER                           :: iseg,js,ks, iloc, jloc
    INTEGER                           :: Kvois, Lvois
    real(kind=long)                   :: AireDK, AireDL
    REAL(kind=long), DIMENSION(2)     :: z
    REAL(kind=long), DIMENSION(4)     :: x, y
    real(kind=long), dimension(Nseg)  :: uSxxK,uSyyK,uSxyK,cSxxK,cSyyK,cSxyK,eSxxK,eSyyK,eSxyK,vSxxK,vSyyK,vSxyK
  
    !-------------------
    ! Debut du programme
    !-------------------
    oldprf = prefix
    prefix = 'TENSEUR'
    !------
    ! Corps
    !------
    
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
       AireDK = ABS( (x(2)-x(1))* (y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1)) )/2.D0
       IF (NombvoisSeg(iseg) == 2) THEN
          Lvois = NumTVoisSeg(2,iseg)
          js = NuSeg(1,iseg) ; ks = NuSeg(2,iseg)
          !!
          x(4) = coordK(1,Lvois) ; y(4) = coordK(2,Lvois)
          !! 
          z(1) = (3*z(1) + x(4))/4.
          z(2) = (3*z(2) + y(4))/4.
          AireDL = ABS( (x(2)-x(4))* (y(3)-y(4)) - (x(3)-x(4))*(y(2)-y(4)) )/2.D0
          !!
       END IF
     !-------------------------------------------
     ! On approxime les matrices de diffusion par 
     ! S^D = 1/AireD*\int_D S(x)dx
     !----------------------------------------
     !
       Select case (choixaniu)
       case(1) 
          !-----------
          ! S = Id
          !----------
          uSxxK(iseg)= CoefDiffU*1.D0
          uSyyK(iseg) = CoefDiffU*1.D0
          uSxyK(iseg)= 0.D0

       case(2)
          !-----------
          ! S = d(x)*Id
          !----------
          if (TypS(1,NumTVoisSeg(1,iseg)) == 100) then 
               uSxxK(iseg)= CoefDiffU * CoefDiffGm *AireDK/AireD(iseg)
               uSyyK(iseg) = CoefDiffU * CoefDiffGm *AireDK/AireD(iseg)
               uSxyK(iseg)= 0.D0
          elseif (TypS(1,NumTVoisSeg(1,iseg)) == 200) then 
               uSxxK(iseg)= CoefDiffU * CoefDiffWm *AireDK/AireD(iseg)
               uSyyK(iseg) = CoefDiffU * CoefDiffWm *AireDK/AireD(iseg)
               uSxyK(iseg)= 0.D0
          elseif (Use_surgery) then 
               uSxxK(iseg)= CoefDiffU * CoefDiffPs *AireDK/AireD(iseg)
               uSyyK(iseg) = CoefDiffU * CoefDiffPs *AireDK/AireD(iseg)
               uSxyK(iseg)= 0.D0
          else 
               uSxxK(iseg)= CoefDiffU * CoefDiffWm *AireDK/AireD(iseg)
               uSyyK(iseg) = CoefDiffU * CoefDiffWm *AireDK/AireD(iseg)
               uSxyK(iseg)= 0.D0
          end if
          if (NombvoisSeg(iseg) == 2) then
               if (TypS(1,NumTVoisSeg(2,iseg)) == 100) then 
                    uSxxK(iseg)= uSxxK(iseg) + CoefDiffU * CoefDiffGm *AireDL/AireD(iseg)
                    uSyyK(iseg) = uSyyK(iseg) + CoefDiffU * CoefDiffGm *AireDL/AireD(iseg)
                    uSxyK(iseg) = uSxyK(iseg) + 0.D0
               elseif (TypS(1,NumTVoisSeg(2,iseg)) == 200) then 
                    uSxxK(iseg)= uSxxK(iseg) + CoefDiffU * CoefDiffWm *AireDL/AireD(iseg)
                    uSyyK(iseg) = uSyyK(iseg) + CoefDiffU * CoefDiffWm *AireDL/AireD(iseg)
                    uSxyK(iseg)= uSxyK(iseg) + 0.D0
               elseif (Use_surgery) then 
                    uSxxK(iseg)= uSxxK(iseg) + CoefDiffU * CoefDiffPs *AireDL/AireD(iseg)
                    uSyyK(iseg) = uSyyK(iseg) + CoefDiffU * CoefDiffPs *AireDL/AireD(iseg)
                    uSxyK(iseg)= uSxyK(iseg) + 0.D0
               else 
                    uSxxK(iseg)= uSxxK(iseg) + CoefDiffU * CoefDiffWm *AireDL/AireD(iseg)
                    uSyyK(iseg) = uSyyK(iseg) + CoefDiffU * CoefDiffWm *AireDL/AireD(iseg)
                    uSxyK(iseg)= uSxyK(iseg) + 0.D0
               end if
          end if
       end select

       Select case (choixanic)
       case(1) 
          !-----------
          ! S = Id
          !----------
          cSxxK(iseg)= 1.D0
          cSyyK(iseg) = 1.D0
          cSxyK(iseg)= 0.D0

       case(2) 
          !-----------
          ! S = D2* Id
          !----------
          cSxxK(iseg)= CoefDiffC
          cSyyK(iseg) = CoefDiffC
          cSxyK(iseg)= 0.D0
       end select

       Select case (choixanie)
       case(1) 
          !-----------
          ! S = Id
          !----------
          eSxxK(iseg)= CoefDiffE*1.D0
          eSyyK(iseg) = CoefDiffE*1.D0
          eSxyK(iseg)= 0.D0
     case(2)
          !-----------
          ! S = d(x)*Id
          !----------
          if (TypS(1,NumTVoisSeg(1,iseg)) == 100) then 
               eSxxK(iseg)= CoefDiffE * CoefDiffGm *AireDK/AireD(iseg)
               eSyyK(iseg) = CoefDiffE * CoefDiffGm *AireDK/AireD(iseg)
               eSxyK(iseg)= 0.D0
          elseif (TypS(1,NumTVoisSeg(1,iseg)) == 200) then 
               eSxxK(iseg)= CoefDiffE * CoefDiffWm *AireDK/AireD(iseg)
               eSyyK(iseg) = CoefDiffE * CoefDiffWm *AireDK/AireD(iseg)
               eSxyK(iseg)= 0.D0
          elseif (Use_surgery) then 
               eSxxK(iseg)= CoefDiffE * CoefDiffPs *AireDK/AireD(iseg)
               eSyyK(iseg) = CoefDiffE * CoefDiffPs *AireDK/AireD(iseg)
               eSxyK(iseg)= 0.D0
          else 
               eSxxK(iseg)= CoefDiffE * CoefDiffWm *AireDK/AireD(iseg)
               eSyyK(iseg) = CoefDiffE * CoefDiffWm *AireDK/AireD(iseg)
               eSxyK(iseg)= 0.D0
          end if
          if (NombvoisSeg(iseg) == 2) then
               if (TypS(1,NumTVoisSeg(2,iseg)) == 100) then 
                    eSxxK(iseg)= eSxxK(iseg) + CoefDiffE * CoefDiffGm *AireDL/AireD(iseg)
                    eSyyK(iseg) = eSyyK(iseg) + CoefDiffE * CoefDiffGm *AireDL/AireD(iseg)
                    eSxyK(iseg) = eSxyK(iseg) + 0.D0
               elseif (TypS(1,NumTVoisSeg(2,iseg)) == 200) then 
                    eSxxK(iseg)= eSxxK(iseg) + CoefDiffE * CoefDiffWm *AireDL/AireD(iseg)
                    eSyyK(iseg) = eSyyK(iseg) + CoefDiffE * CoefDiffWm *AireDL/AireD(iseg)
                    eSxyK(iseg)= eSxyK(iseg) + 0.D0
               elseif (Use_surgery) then 
                    eSxxK(iseg)= eSxxK(iseg) + CoefDiffE * CoefDiffPs *AireDL/AireD(iseg)
                    eSyyK(iseg) = eSyyK(iseg) + CoefDiffE * CoefDiffPs *AireDL/AireD(iseg)
                    eSxyK(iseg)= eSxyK(iseg) + 0.D0
               else 
                    eSxxK(iseg)= eSxxK(iseg) + CoefDiffE * CoefDiffWm *AireDL/AireD(iseg)
                    eSyyK(iseg) = eSyyK(iseg) + CoefDiffE * CoefDiffWm *AireDL/AireD(iseg)
                    eSxyK(iseg)= eSxyK(iseg) + 0.D0
               end if
          end if
       end select

       Select case (choixaniv)
       case(1) 
          !-----------
          ! S = Id
          !----------
          vSxxK(iseg)= 1.D0
          vSyyK(iseg) = 1.D0
          vSxyK(iseg)= 0.D0
       case(2) 
          !-----------
          ! S = D4* Id
          !----------
          vSxxK(iseg)= CoefDiffV
          vSyyK(iseg) = CoefDiffV
          vSxyK(iseg)= 0.D0
       end select
    end DO

    Allocate(uTKL(Nseg),uTKeLe(Nseg),uetaSSe(Nseg))
    Allocate(cTKL(Nseg),cTKeLe(Nseg),cetaSSe(Nseg))
    Allocate(eTKL(Nseg),eTKeLe(Nseg),eetaSSe(Nseg))
    Allocate(vTKL(Nseg),vTKeLe(Nseg),vetaSSe(Nseg))

    DO iseg = 1,Nseg
        !
        ! les transmissibiltés de diffusion TKL
        ! TKL = 0.5*Msig^2/AireD*S^D*NsigK.NsigK
        !--------------------------------------
        !
        uTKL(iseg) = ( (uSxxK(iseg)*NsigK(1,iseg)+uSxyK(iseg) *NsigK(2,iseg) )*NsigK(1,iseg)&
             & + ( uSxyK(iseg)*NsigK(1,iseg)+uSyyK(iseg)*NsigK(2,iseg) )*NsigK(2,iseg))* &
             &(Msig(iseg)**2)/(2*AireD(iseg))

        cTKL(iseg) = ( (cSxxK(iseg)*NsigK(1,iseg)+cSxyK(iseg) *NsigK(2,iseg) )*NsigK(1,iseg)&
             & + ( cSxyK(iseg)*NsigK(1,iseg)+cSyyK(iseg)*NsigK(2,iseg) )*NsigK(2,iseg))* &
             &(Msig(iseg)**2)/(2*AireD(iseg))

        eTKL(iseg) = ( (eSxxK(iseg)*NsigK(1,iseg)+eSxyK(iseg) *NsigK(2,iseg) )*NsigK(1,iseg)&
             & + ( eSxyK(iseg)*NsigK(1,iseg)+eSyyK(iseg)*NsigK(2,iseg) )*NsigK(2,iseg))* &
             &(Msig(iseg)**2)/(2*AireD(iseg))

        vTKL(iseg) = ( (vSxxK(iseg)*NsigK(1,iseg)+vSxyK(iseg) *NsigK(2,iseg) )*NsigK(1,iseg)&
             & + ( vSxyK(iseg)*NsigK(1,iseg)+vSyyK(iseg)*NsigK(2,iseg) )*NsigK(2,iseg))* &
             &(Msig(iseg)**2)/(2*AireD(iseg))
        !
        ! les transmissibiltés de diffusion TK*L*
        !  TKeLe = 0.5*Msige^2/AireD*S^D*NsigeKe.NsigeKe
        !----------------------------------------
        !
        uTKeLe(iseg) = ( (uSxxK(iseg)*NsigeKe(1,iseg)+uSxyK(iseg) *NsigeKe(2,iseg) )*NsigeKe(1,iseg)&
             & + ( uSxyK(iseg)*NsigeKe(1,iseg)+uSyyK(iseg)*NsigeKe(2,iseg) )*NsigeKe(2,iseg))* &
             &(Msige(iseg)**2)/(2*AireD(iseg))

        cTKeLe(iseg) = ( (cSxxK(iseg)*NsigeKe(1,iseg)+cSxyK(iseg) *NsigeKe(2,iseg) )*NsigeKe(1,iseg)&
             & + ( cSxyK(iseg)*NsigeKe(1,iseg)+cSyyK(iseg)*NsigeKe(2,iseg) )*NsigeKe(2,iseg))* &
             &(Msige(iseg)**2)/(2*AireD(iseg))

        eTKeLe(iseg) = ( (eSxxK(iseg)*NsigeKe(1,iseg)+eSxyK(iseg) *NsigeKe(2,iseg) )*NsigeKe(1,iseg)&
             & + ( eSxyK(iseg)*NsigeKe(1,iseg)+eSyyK(iseg)*NsigeKe(2,iseg) )*NsigeKe(2,iseg))* &
             &(Msige(iseg)**2)/(2*AireD(iseg))

        vTKeLe(iseg) = ( (vSxxK(iseg)*NsigeKe(1,iseg)+vSxyK(iseg) *NsigeKe(2,iseg) )*NsigeKe(1,iseg)&
             & + ( vSxyK(iseg)*NsigeKe(1,iseg)+vSyyK(iseg)*NsigeKe(2,iseg) )*NsigeKe(2,iseg))* &
             &(Msige(iseg)**2)/(2*AireD(iseg))
        !
        ! les transmissibiltés de convection eta_sigma,sigma*
        ! etaSSe = 0.5*Msig*Msige/AireD*S^D*NsigK.NsigeKe
        !----------------------------------------------------
        !
        uetaSSe(iseg) = ( (uSxxK(iseg)*NsigK(1,iseg)+uSxyK(iseg) *NsigK(2,iseg) )*NsigeKe(1,iseg)&
             & + ( uSxyK(iseg)*NsigK(1,iseg)+uSyyK(iseg)*NsigK(2,iseg) )*NsigeKe(2,iseg))* &
             &(Msig(iseg)*Msige(iseg))/(2*AireD(iseg))

        cetaSSe(iseg) = ( (cSxxK(iseg)*NsigK(1,iseg)+cSxyK(iseg) *NsigK(2,iseg) )*NsigeKe(1,iseg)&
             & + ( cSxyK(iseg)*NsigK(1,iseg)+cSyyK(iseg)*NsigK(2,iseg) )*NsigeKe(2,iseg))* &
             &(Msig(iseg)*Msige(iseg))/(2*AireD(iseg))

        eetaSSe(iseg) = ( (eSxxK(iseg)*NsigK(1,iseg)+eSxyK(iseg) *NsigK(2,iseg) )*NsigeKe(1,iseg)&
             & + ( eSxyK(iseg)*NsigK(1,iseg)+eSyyK(iseg)*NsigK(2,iseg) )*NsigeKe(2,iseg))* &
             &(Msig(iseg)*Msige(iseg))/(2*AireD(iseg))

        vetaSSe(iseg) = ( (vSxxK(iseg)*NsigK(1,iseg)+vSxyK(iseg) *NsigK(2,iseg) )*NsigeKe(1,iseg)&
             & + ( vSxyK(iseg)*NsigK(1,iseg)+vSyyK(iseg)*NsigK(2,iseg) )*NsigeKe(2,iseg))* &
             &(Msig(iseg)*Msige(iseg))/(2*AireD(iseg))
        !
     END DO
  
    !------------
    ! Impressions
    !------------

    !-----------------
    ! Fin du programme
    !-----------------
    prefix = oldprf
  
  END SUBROUTINE transmis
  