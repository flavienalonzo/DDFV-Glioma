Module fsource

Contains

   subroutine surge(UU,CC,EE,VV,Ndim)
      use longr
      USE imprime
      USE parmmage
      implicit none 
      Integer , intent(in) :: Ndim
      REAL(kind=long), DIMENSION(Ndim), intent(inout)  :: UU, CC, EE, VV
      integer, dimension(NsInt)   :: NbTpS
      integer :: jt, m, is
      NbTpS = 0
      DO jt = 1,Nbt
         if (TypS(1,jt) == 100) then
            UU(jt + NsInt) = UU(jt + NsInt)
            CC(jt + NsInt) = CC(jt + NsInt)
            EE(jt + NsInt) = EE(jt + NsInt)
            VV(jt + NsInt) = VV(jt + NsInt)
         elseif (TypS(1,jt) == 200) then
            UU(jt + NsInt) = UU(jt + NsInt)
            CC(jt + NsInt) = CC(jt + NsInt)
            EE(jt + NsInt) = EE(jt + NsInt)
            VV(jt + NsInt) = VV(jt + NsInt)
         elseif (TypS(1,jt) == 300) then
            UU(jt + NsInt) = 0.D0
            CC(jt + NsInt) = 0.D0
            EE(jt + NsInt) = 0.D0
            VV(jt + NsInt) = 0.D0
         elseif (TypS(1,jt) == 400) then
            UU(jt + NsInt) = 0.D0
            CC(jt + NsInt) = 0.D0
            EE(jt + NsInt) = 0.D0
            VV(jt + NsInt) = 0.D0
         elseif (TypS(1,jt) == 500) then
            UU(jt + NsInt) = 0.D0
            CC(jt + NsInt) = 0.D0
            EE(jt + NsInt) = 0.D0
            VV(jt + NsInt) = 0.D0
         end if
         Do m = 1, typegeomaille
            UU(NuSoK (m,jt)) = UU(NuSoK (m,jt)) + UU(jt + NsInt)
            CC(NuSoK (m,jt)) = CC(NuSoK (m,jt)) + CC(jt + NsInt)
            EE(NuSoK (m,jt)) = EE(NuSoK (m,jt)) + EE(jt + NsInt)
            VV(NuSoK (m,jt)) = VV(NuSoK (m,jt)) + VV(jt + NsInt)
            NbTpS(NuSoK (m,jt)) = NbTpS(NuSoK (m,jt)) + 1
         END Do
      END DO
      DO is = 1, NsInt
         UU(is) = UU(is)/NbTpS(is)
         CC(is) = CC(is)/NbTpS(is)
         EE(is) = EE(is)/NbTpS(is)
         VV(is) = VV(is)/NbTpS(is)
      END DO

   end subroutine surge

   function fu(u,e)
      use longr
      implicit none 
      REAL(kind=long), INTENT(in)     :: e,u
      REAL(kind=long)  :: fu

      if (u<0) then 
         fu = 0.D0
      elseif (u>1) then 
         fu = 0.D0
      elseif (u+e>1) then 
         fu = 0.D0
      else 
         fu = u*(1-u-e)
      end if

   end function

   function deriveefu(u,e)
      use longr
      implicit none 
      REAL(kind=long), INTENT(in)     :: e,u
      REAL(kind=long)  :: deriveefu

      if (u<0) then 
         deriveefu = 0.D0
      elseif (u>1) then 
         deriveefu = 0.D0
      elseif (u+e>1) then 
         deriveefu = 0.D0
      else 
         deriveefu = (1-e-2*u)
      end if

   end function

   function fe(e,u)
      use longr
      implicit none 
      REAL(kind=long), INTENT(in)     :: e,u
      REAL(kind=long)  :: fe

      if (e<0) then 
         fe = 0.D0
      elseif (e>1) then 
         fe = 0.D0
      elseif (u+e>1) then 
         fe = 0.D0
      else 
         fe = e*(1-u-e)
      end if

   end function

   function deriveefe(e,u)
      use longr
      implicit none 
      REAL(kind=long), INTENT(in)     :: e,u
      REAL(kind=long)  :: deriveefe

      if (e<0) then 
         deriveefe = 0.D0
      elseif (e>1) then 
         deriveefe = 0.D0
      elseif (u+e>1) then 
         deriveefe = 0.D0
      else 
         deriveefe = (1-u-2*e)
      end if

   end function

   function Tu(temps,u)
      use longr 
      implicit none 
      REAL(kind=long), INTENT(in)     :: temps,u
      REAL(kind=long)  :: Tu 
      if (Use_chemo .and. Use_radio) then 
         if (ABS(temps-16.5)<=2.5 .or. ABS(temps-23.5)<=2.5 .or. ABS(temps-30.5)<=2.5 .or. ABS(temps-37.5)<=2.5 .or. &
         & ABS(temps-44.5)<=2.5 .or. ABS(temps-51.5)<=2.5) then
            Tu = -(Dose_chemo+Gain_radio*Dose_radio)*u
         else 
            Tu = 0.D0
         end if
      elseif (.not.Use_chemo .and. Use_radio) then 
         if (ABS(temps-35)<=21) then
            Tu = -Dose_radio*u
         else 
            Tu = 0.D0
         end if
      elseif (Use_chemo .and. .not.Use_radio) then 
         if (ABS(temps-16.5)<=2.5 .or. ABS(temps-23.5)<=2.5 .or. ABS(temps-30.5)<=2.5 .or. ABS(temps-37.5)<=2.5 .or. &
         & ABS(temps-44.5)<=2.5 .or. ABS(temps-51.5)<=2.5) then
            Tu = -Dose_chemo*u
         else 
            Tu = 0.D0
         end if
      else
         Tu = 0.D0
      end if
   end function

   function deriveeTu(temps,u)
      use longr 
      implicit none 
      REAL(kind=long), INTENT(in)     :: temps,u
      REAL(kind=long)  :: deriveeTu 

      if (Use_chemo .and. Use_radio) then 
         if (ABS(temps-16.5)<=2.5 .or. ABS(temps-23.5)<=2.5 .or. ABS(temps-30.5)<=2.5 .or. ABS(temps-37.5)<=2.5 .or. &
         & ABS(temps-44.5)<=2.5 .or. ABS(temps-51.5)<=2.5) then
            deriveeTu = -(Dose_chemo+Gain_radio*Dose_radio)
         else 
            deriveeTu = 0.D0
         end if
      elseif (.not.Use_chemo .and. Use_radio) then 
         if (ABS(temps-35)<=21) then
            deriveeTu = -Dose_radio
         else 
            deriveeTu = 0.D0
         end if
      elseif (Use_chemo .and. .not.Use_radio) then 
         if (ABS(temps-16.5)<=2.5 .or. ABS(temps-23.5)<=2.5 .or. ABS(temps-30.5)<=2.5 .or. ABS(temps-37.5)<=2.5 .or. &
         & ABS(temps-44.5)<=2.5 .or. ABS(temps-51.5)<=2.5) then
            deriveeTu = -Dose_chemo
         else 
            deriveeTu = 0.D0
         end if
      else
         deriveeTu = 0.D0
      end if  
   end function

   function hu(c)
      use longr
      implicit none
      REAL(kind=long), INTENT(in)     :: c
      REAL(kind=long)  :: hu 
      if (c>=Chypo) then 
         hu = 1.D0
      else 
         hu = 0.D0
      end if
   end function

   function gv(c)
      use longr
      implicit none
      REAL(kind=long), INTENT(in)     :: c
      REAL(kind=long)  :: gv
      if (c>=Cnecro .and. c<=Chypo) then 
         gv = c 
      else 
         gv = 0.D0
      end if
   end function

  FUNCTION fsourceu(t,x,y,m) ! m represente le choix de probleme
    !--------
    ! Modules
    !--------
    USE longr

    IMPLICIT NONE

    !--------------------------
    ! Declaration des arguments
    !--------------------------
    REAL(kind=long), INTENT(in)     :: t, x , y
    Integer                         :: m
    REAL(kind=long)                 :: fsourceu, omega
    !----------------------------------
    ! Declaration des variables locales
    !----------------------------------
    CHARACTER(len=6)                :: oldprf

    !-------------------
    ! Debut du programme
    !-------------------
    oldprf = prefix
    prefix = 'FSOUR'

    !------
    ! Corps
    !------
    Select case (m)
    case (1)
       fsourceu  = 0.D0
    case (2)
       fsourceu = x*x - y*y
    case(3)
       omega    = 20.D0*atan(1.)
       fsourceu =  2* omega**2 *cos(omega*(x+y))
    case(66)! cas degenere avec a(u) =u(1-u)
       fsourceu = -2.D0*(1.D0-2.D0*(x+y))
    case(77)! cas degenere avec a(u) =u et A(u) = u**2/2, u =  x+y
       fsourceu = -2.D0
    case(99)
       fsourceu = x+y
    case default
       stop 'fsourceu'
    end Select
    !-----------------
    ! Fin du programme
    !-----------------
    prefix = oldprf

    RETURN
  END FUNCTION fsourceu
  !=====================================================
  !=====================================================
  FUNCTION gbord(t,x,y,m)
    !--------
    ! Modules
    !--------
    USE longr

    IMPLICIT NONE

    !--------------------------
    ! Declaration des arguments
    !--------------------------
    REAL(kind=long), INTENT(in)     ::  x, y,t
    integer :: m
    REAL(kind=long)                 :: gbord
    !----------------------------------
    ! Declaration des variables locales
    !----------------------------------
    CHARACTER(len=6)      :: oldprf

    !-------------------
    ! Debut du programme
    !-------------------
    oldprf = prefix
    prefix = 'Gbord'

    !------
    ! Corps
    !------
    select case(m)
    case(1)
       gbord  = x + y
    case(2) 
       gbord  = x*x - y*y
    case(3)
       gbord = sin(x)*exp(-y)
    case(4)
       gbord = sin(pi*x)*sin(pi*y)*exp(-2*t*pi**2)
    case(5)
       gbord = ( 1.D0 + (cos(pi*x))*exp(-t*pi**2) )/2.D0
    case(66)! cas degenere avec a(u) =u(1-u)
       gbord = epsilon
    case(77)! cas degenere avec a(u) =u
       gbord = -20.D0
    case(80)! cas spécifique du Laplacien avec Neumann
       gbord = 1.D0
    case(99)
       gbord = 0.D0
    case default
       stop 'gbord'
    end select
    !-----------------
    ! Fin du programme
    !-----------------
    prefix = oldprf

    RETURN
  END FUNCTION gbord

  !=======================
  function AdegenU(y) !Phi=FU
    !=======================
    USE longr
    Real(kind=long) :: y,AdegenU

    select case (ChoixAdegU)
    case(0)  !!test 0
       AdegenU= CoefdiffuAdeg * y
    case(1)  !!test 1
       AdegenU=CoefdiffuAdeg* (1.D0/2.D0 - y/3.D0 )*y**2
    case(2) !!test 2
       AdegenU=CoefdiffuAdeg* (1.D0/3.D0 -y/2.D0 +y*y/5.D0 )*y**3
    case(3) !!test 3
       AdegenU=CoefdiffuAdeg* (1.D0/5.D0 -2.D0*y/3.D0 + (6./7.D0)*y**2 -(y**3)/2.D0 +(y**4)/9.D0)*y**5
    case(4) !!test 4
       AdegenU=CoefdiffuAdeg*y**2
    case default   
       stop 'Adege'
    end select

  end function AdegenU

  !=============================
  function DerivAdegU(y) !dPhi=fu
    !===========================
    USE longr
    REAL( kind=long) :: y,DerivAdegU
    !
    select case (ChoixAdegU)
    case(0) !!test 0
       DerivAdegU = CoefdiffuAdeg
    case(1) !!test 1
       DerivAdegU = CoefdiffuAdeg * y*( 1.D0 - y )
    case(2) !!test 2
       DerivAdegU= CoefdiffuAdeg * ( y*(1.D0 - y) )**2
    case(3) !!test 3
       DerivAdegU= CoefdiffuAdeg * ( y*(1.D0 - y) )**4
    case(4) !!test 4
       DerivAdegU=CoefdiffuAdeg*2.*y
    case default  
       stop 'DerAdeg'

    end select
  end function DerivAdegU

  !==========================================================
  function xiU(y) ! Squared-Kirchoff transform: integral de rmU
    !========================================================
    USE longr
    REAL( kind=long) :: y, xiU, CoefrmU
    !
    CoefrmU = sqrt(CoefdiffuAdeg)
    select case (ChoixAdegU)
    case(0) !!test 0
       xiU = CoefrmU* y
    case(1)  !!test 1
       xiU = CoefrmU*( ((2*y-1)*sqrt(y*(1.D0 - y))/4) + (acos(-2*y+1.))/8)
    case(2) !!test 2
       xiU = CoefrmU * (1.D0/2.D0 -y/3.D0 )*y**2
    case(3) !!test 3
       xiU = CoefrmU * (1.D0/3.D0 -y/2.D0 +y*y/5.D0 )*y**3
    case(4) !!test 4
       xiU=sqrt(CoefdiffuAdeg*2.)*(2./3.)*y**(3./2.)
    case default  
       stop 'xiU'

    end select
  end function xiU

  !=======================================
  function rmU(y) ! racine carre de DerAdeg
    !=====================================
    USE longr
    REAL( kind=long) :: y, rmU, CoefrmU
    !
    CoefrmU = sqrt(CoefdiffuAdeg)
    select case (ChoixAdegU)
    case(0) !!test 0
       rmU = CoefrmU
    case(1)  !!test 1
       rmU = CoefrmU *sqrt(y*(1.D0 - y))
    case(2) !!test 2
       rmU = CoefrmU * y*(1.D0 - y)
    case(3) !!test 3
       rmU = CoefrmU * ( y*(1.D0 - y) )**2
    case(4) !!test 4
       rmU = sqrt(CoefdiffuAdeg*2.*y)
    case default  
       stop 'rmU'

    end select
  end function rmU

  !=============
  function ln(y)
    !===========
    USE longr
    REAL( kind=long) :: x,y, ln
    x = max(y,epsilon)
    ln = log(x)

  end function ln

  !=============
  function Derivln(y)
    !===========
    USE longr
    REAL( kind=long) :: x,y, Derivln
    x = max(y,epsilon)
    Derivln = 1.D0/(x)

  end function Derivln
  
  !=======================
  function rmCroitU(y)
    !=======================
    USE longr
    REAL( kind=long) :: y, z, rmCroitU, CoefrmU,umax
    !
    CoefrmU = sqrt(CoefdiffuAdeg)
    select case(ChoixAdegU)
    case(0) !!test 0
       rmCroitU = CoefrmU
    case(1)  !!test 1
       umax= 0.5D0
       z=min(y,umax)
       rmCroitU =  CoefrmU *sqrt( z*(1.D0-z))
    case(2)  !!test 2
       umax= 0.5D0
       z=min(y,umax)
       rmCroitU =  CoefrmU * z*(1.D0-z)
    case(3)  !!test 3
       umax= 0.5D0
       z=min(y,umax)
       rmCroitU=  CoefrmU * (z*(1.D0-z))**2
    case(4) !!test 4
       rmCroitU = sqrt(CoefdiffuAdeg*2.*y)
    case default
       stop 'pb rmCroitU'
    end select
  end function rmCroitU
  !=======================
  function DerrmCroitU(y)
    !=======================
    USE longr
    REAL( kind=long) :: y, DerrmCroitU, CoefrmU,umax
    !
    CoefrmU = sqrt(CoefdiffuAdeg)
    select case(ChoixAdegU) 
    case(0) !!test 0
       DerrmCroitU = 0.D0
    case(1)  !!test 1
       umax= 0.5D0
       y = max(y, epsilon)
       if (y <=umax) then
          DerrmCroitU =  CoefrmU * (1-2.D0*y)/(2*sqrt( y*(1.D0-y)))
       else
          DerrmCroitU = 0.D0
       end if
    case(2)  !!test 2
       umax= 0.5D0
       if (y <=umax) then
          DerrmCroitU =  CoefrmU * (1-2.D0*y)
       else
          DerrmCroitU = 0.D0
       end if
    case(3)  !!test 3
       umax= 0.5D0
       if (y <=umax) then
          DerrmCroitU =  CoefrmU * 2.D0*y*(1.D0 -3.D0*y +2.D0*y**2)
       else
          DerrmCroitU = 0.D0
       end if
    case(4) !!test 4
       DerrmCroitU = CoefrmU/(sqrt(2.*y))
    case default
       stop 'pb DerrmCroitU'
    end select
  end function DerrmCroitU

  !=======================
  function rmDeCroitU(y)
    !=======================
    USE longr
    REAL( kind=long) :: y, z, rmDeCroitU, CoefrmU,umax
    !
    CoefrmU = sqrt(CoefdiffuAdeg)
    select case(ChoixAdegU)
    case(0) !!test 0
       rmDeCroitU = CoefrmU
    case(1)  !!test 1
       umax= 0.5D0
       z = max(y,umax)
       rmDeCroitU =  CoefrmU *sqrt( z*(1.D0-z)) -sqrt( umax*(1.D0-umax) )
    case(2)  !!test 2
       umax= 0.5D0
       z=max(y,umax)
       rmDeCroitU =  CoefrmU * ( z*(1.D0-z) - umax*(1.D0-umax) )
    case(3)  !!test 3
       umax= 0.5D0
       z=max(y,umax)
       rmDeCroitU=  CoefrmU * ( (z*(1.D0-z))**2 - (umax*(1.D0-umax))**2 )
    case(4) !!test 4
       rmDeCroitU = 0.D0
    case default
       stop 'pb rmDeCroitU'
    end select
  end function rmDeCroitU

  !=======================
  function DerrmDeCroitU(y)
    !=====================
    USE longr
    REAL( kind=long) :: y, DerrmDeCroitU, CoefrmU
    !
    CoefrmU = sqrt(CoefdiffuAdeg)
    select case(ChoixAdegU)
    case(0) !!test 0
       DerrmDeCroitU = 0.D0
    case(1)  !!test 1
       umax= 0.5D0
       y = min(y,1.0-epsilon)
       if (y > umax) then
          DerrmDeCroitU =  CoefrmU * (1-2.D0*y)/(2*sqrt( y*(1.D0-y)))
       else
          DerrmDeCroitU = 0.D0
       end if
    case(2)  !!test 2
       umax= 0.5D0
       if (y > umax) then
          DerrmDeCroitU =  CoefrmU * (1-2.D0*y)
       else
          DerrmDeCroitU = 0.D0
       end if
    case(3)  !!test 3
       umax= 0.5D0
       if (y > umax) then
          DerrmDeCroitU =  CoefrmU * 2.D0*y*(1.D0 -3.D0*y +2.D0*y**2)
       else
          DerrmDeroit = 0.D0
       end if
    case(4) !!test 4
       DerrmDeroit = 0.D0
    case default
       stop 'pb DerrmDeCroitU'
    end select
  end function DerrmDeCroitU


  !=======================
  function DerrmU(y)
    !=======================
    USE longr
    REAL( kind=long) :: y, DerrmU, CoefrmU
    !
    CoefrmU = sqrt(CoefdiffuAdeg)
    select case(ChoixAdegU)
  
    case(1)  !!test 1
       DerrmU = 0.D0
    case(2)  !!test 2
       DerrmU = CoefrmU *(1.D0 - 2.D0*y)
    case(3)  !!test 3
       DerrmU= CoefrmU * 2.D0*y*(1.D0 -3.D0*y +2.D0*y**2)
    case default
       stop 'pb DerrmU'
    end select
  end function DerrmU

!=======================
  function MuU(y)
    !=======================
    USE longr
    REAL( kind=long) :: y, MuU


    select case(ChoixMuU)
    case(1)
       MuU = CoefTranspMuU 
    case(2)
       MuU = CoefTranspMuU *y
    case(3)
       MuU= CoefTranspMuU *y*(1.D0-y) 
    case(4)
       MuU= CoefTranspMuU *y*(1.D0-y)**2
    case(5)
       IF ( y <= ubar ) THEN
          MuU= CoefTranspMuU *y*( 1.D0-(y/ubar)**gamma)
       ELSE
          MuU= 0.D0
       END IF
       case(6)
          MuU= CoefTranspMuU *( y*(1.D0-y) )**2
    case default
       stop 'pb MuU'
    end select

  end function MuU

  !=======================
  function MuCroitU(y)
    !=======================
    USE longr
    implicit none
    REAL( kind=long) :: y, z, MuCroitU, umax
    !integer :: choix



    select case(ChoixMuU)
    case(1)
       MuCroitU = CoefTranspMuU 
    case(2)
       MuCroitU = CoefTranspMuU *y
    case(3)
       umax= 0.5D0
       z=min(y,umax)
       MuCroitU=  CoefTranspMuU * z*(1.D0-z)
    case(4)
       umax= 1.D0/3.D0
       z=min(y,umax)
       MuCroitU= CoefTranspMuU * z*(1.D0-z)**2
    case(5)
       umax= (ubar)/( (gamma+1.D0)**(1.D0/gamma) )
       z=min(y,umax)
       IF ( z <= ubar ) THEN
          MuCroitU= CoefTranspMuU * z*(1.D0-(z/ubar)**gamma)
       ELSE
          MuCroitU= 0.D0
       END IF
    case(6)
       umax= 0.5D0
       z=min(y,umax)
       MuCroitU=  CoefTranspMuU * (z*(1.D0-z))**2
    case default
       stop 'pb MuCroitU'
    end select

  end function MuCroitU

  !=======================
  function DerivMuCroitU(y)
    !=======================
    USE longr
    REAL( kind=long) :: y, DerivMuCroitU, umax
    !integer :: choix
    !

    select case(ChoixMuU)
    case(1)
       DerivMuCroitU = 0.D0 
    case(2)
       DerivMuCroitU = CoefTranspMuU
    case(3)
       umax= 0.5D0
       if (y <=umax) then
          DerivMuCroitU =  CoefTranspMuU * (1-2.D0*y)
       else
          DerivMuCroitU = 0.D0
       end if
    case(4)
       umax= 1.D0/3.D0
       if (y < umax) then
          DerivMuCroitU = CoefTranspMuU * ( (1.D0-y)**2 -2.D0*y*(1.D0-y) )
       else
          DerivMuCroitU = 0.D0
       end if
    case(5)
       umax= (ubar)/( (gamma+1.D0)**(1.D0/gamma) )
       IF (y < umax .AND. y< ubar) THEN !! Although we know that y<=ubar
          DerivMuCroitU = CoefTranspMuU * (1.D0 - ( gamma+1.D0 )*( (y/ubar)**gamma) )
       ELSE
          DerivMuCroitU = 0.D0
       END IF
    case(6)
       umax= 0.5D0
       if (y <=umax) then
          DerivMuCroitU =  CoefTranspMuU * 2.D0*y*(1.D0-y)*(1-2.D0*y)
       else
          DerivMuCroitU = 0.D0
       end if
    case default
       stop 'pb DerivMuCroitU'
    end select

  end function DerivMuCroitU

  !=======================
  function MuDeCroitU(y)
    !=======================
    USE longr
    REAL( kind=long) :: y, z, MuDeCroitU, umax
    !

    select case(ChoixMuU)
    case(1)
       MuDeCroitU = 0.D0 
    case(2)
       MuDeCroitU = 0.D0
    case(3)
       umax= 0.5D0
       z=max(y,umax)
       MuDeCroitU=  CoefTranspMuU * ( z*(1.D0-z) - umax*(1.D0-umax) )
    case(4)
       umax= 1.D0/3.D0
       z=max(y,umax)
       MuDeCroitU= CoefTranspMuU * ( z*(1.D0-z)**2 -umax*(1.D0-umax)**2 ) 
    case(5)
       umax= (ubar)/( (gamma+1.D0)**(1.D0/gamma) )
       z=max(y,umax)
       IF (z <= ubar ) THEN
          MuDeCroitU= CoefTranspMuU * ( z*(1.D0-(z/ubar)**gamma) )
       ELSE
          MuDeCroitU= 0.D0
       END IF
       IF ( umax <= ubar )  THEN
          MuDeCroitU= MuDeCroitU  - CoefTranspMuU*umax*( 1.D0-(umax/ubar)**gamma ) 
       END IF
    case(6)
       umax= 0.5D0
       z=max(y,umax)
       MuDeCroitU=  CoefTranspMuU * ( (z*(1.D0-z))**2 - (umax*(1.D0-umax))**2 )
    case default
       stop 'pb MuCroit'
    end select

  end function MuDeCroitU

  !=======================
  function DerivMuDeCroitU(y)
    !=======================
    USE longr
    REAL( kind=long) :: y, DerivMuDeCroitU, umax
    !

    select case(ChoixMuU)
    case(1)
       DerivMuDeCroitU = 0.D0 
    case(2)
       DerivMuDeCroitU = 0.D0
    case(3)
       umax= 0.5D0
       if (y > umax) then
          DerivMuDeCroitU =  CoefTranspMuU * (1.D0-2.D0*y)
       else
          DerivMuDeCroitU = 0.D0
       end if
    case(4)
       umax= 1.D0/3.D0
       if (y > umax) then
          DerivMuDeCroitU = CoefTranspMuU * ( (1.D0-y)**2 -2.D0*y*(1.D0-y) )
       else
          DerivMuDeCroitU = 0.D0
       end if
    case(5)
       umax= (ubar)/( (gamma+1.D0)**(1.D0/gamma) )
       IF (y > umax .AND. y < ubar) THEN
          DerivMuDeCroitU = CoefTranspMuU* ( 1.D0-(gamma+1.D0)*( (y/ubar)**gamma ) )
       ELSE
          DerivMuDeCroitU = 0.D0
       END IF
    case(6)
       umax= 0.5D0
       if (y > umax) then
          DerivMuDeCroitU =  CoefTranspMuU * 2.D0*y*(1.D0-y)*(1.D0-2.D0*y)
       else
          DerivMuDeCroitU = 0.D0
       end if
    case default
       stop 'pb DerivMuDeCroitU'
    end select

  end function DerivMuDeCroitU


  !=======================
  function DerivMuU(y)
    !=======================
    USE longr
    REAL( kind=long) :: y,DerivMuU


    select case(ChoixMuU)
    case(1)
       DerivMuU = 0.D0
    case(2)
       DerivMuU = CoefTranspMuU *1.D0
    case(3)
       DerivMuU = CoefTranspMuU * ( 1.D0 - 2.D0*y )
    case(4)
       DerivMuU = CoefTranspMuU *( 1.D0 - y ) * ( 1.D0 - 3.D0*y )
    case(5)
       IF ( y < ubar ) THEN
          DerivMuU= CoefTranspMuU *( 1.D0-( gamma+ 1.D0 )*( (y/ubar)**gamma ) )
       ELSE
          DerivMuU= 0.D0
       END IF
    case(6)
       DerivMuU = CoefTranspMuU * 2.D0 *y*(1.D0-y)*( 1.D0 - 2.D0*y )
    case default
       stop 'pb DerivMuU'
    end select
  end function DerivMuU

    !=======================
  function AdegenE(y) !Phi=FU
   !=======================
   USE longr
   Real(kind=long) :: y,AdegenE

   select case (ChoixAdegU)
   case(0)  !!test 0
      AdegenE= CoefdiffeAdeg * y
   case(1)  !!test 1
      AdegenE=CoefdiffeAdeg* (1.D0/2.D0 - y/3.D0 )*y**2
   case(2) !!test 2
      AdegenE=CoefdiffeAdeg* (1.D0/3.D0 -y/2.D0 +y*y/5.D0 )*y**3
   case(3) !!test 3
      AdegenE=CoefdiffeAdeg* (1.D0/5.D0 -2.D0*y/3.D0 + (6./7.D0)*y**2 -(y**3)/2.D0 +(y**4)/9.D0)*y**5
   case(4) !!test 4
      AdegenE=CoefdiffeAdeg*y**2
   case default   
      stop 'Adege'
   end select

 end function AdegenE

 !=============================
 function DerivAdegE(y) !dPhi=fu
   !===========================
   USE longr
   REAL( kind=long) :: y,DerivAdegE
   !
   select case (ChoixAdegE)
   case(0) !!test 0
      DerivAdegE = CoefdiffeAdeg
   case(1) !!test 1
      DerivAdegE = CoefdiffeAdeg * y*( 1.D0 - y )
   case(2) !!test 2
      DerivAdegE= CoefdiffeAdeg * ( y*(1.D0 - y) )**2
   case(3) !!test 3
      DerivAdegE= CoefdiffeAdeg * ( y*(1.D0 - y) )**4
   case(4) !!test 4
      DerivAdegE=CoefdiffeAdeg*2.*y
   case default  
      stop 'DerAdeg'

   end select
 end function DerivAdegE

 !==========================================================
 function xiE(y) ! Squared-Kirchoff transform: integral de rmU
   !========================================================
   USE longr
   REAL( kind=long) :: y, xiE, CoefrmE
   !
   CoefrmE = sqrt(CoefdiffeAdeg)
   select case (ChoixAdegE)
   case(0) !!test 0
      xiE = CoefrmE* y
   case(1)  !!test 1
      xiE = CoefrmE*( ((2*y-1)*sqrt(y*(1.D0 - y))/4) + (acos(-2*y+1.))/8)
   case(2) !!test 2
      xiE = CoefrmE * (1.D0/2.D0 -y/3.D0 )*y**2
   case(3) !!test 3
      xiE = CoefrmE * (1.D0/3.D0 -y/2.D0 +y*y/5.D0 )*y**3
   case(4) !!test 4
      xiE=sqrt(CoefdiffeAdeg*2.)*(2./3.)*y**(3./2.)
   case default  
      stop 'xiE'

   end select
 end function xiE

 !=======================================
 function rmE(y) ! racine carre de DerAdeg
   !=====================================
   USE longr
   REAL( kind=long) :: y, rmE, CoefrmE
   !
   CoefrmE = sqrt(CoefdiffeAdeg)
   select case (ChoixAdegE)
   case(0) !!test 0
      rmE = CoefrmE
   case(1)  !!test 1
      rmE = CoefrmE *sqrt(y*(1.D0 - y))
   case(2) !!test 2
      rmE = CoefrmE * y*(1.D0 - y)
   case(3) !!test 3
      rmE = CoefrmE * ( y*(1.D0 - y) )**2
   case(4) !!test 4
      rmE = sqrt(CoefdiffeAdeg*2.*y)
   case default  
      stop 'rmE'

   end select
 end function rmE
 
 !=======================
 function rmCroitE(y)
   !=======================
   USE longr
   REAL( kind=long) :: y, z, rmCroitE, CoefrmE,umax
   !
   CoefrmE = sqrt(CoefdiffeAdeg)
   select case(ChoixAdegE)
   case(0) !!test 0
      rmCroitE = CoefrmE
   case(1)  !!test 1
      umax= 0.5D0
      z=min(y,umax)
      rmCroitE =  CoefrmE *sqrt( z*(1.D0-z))
   case(2)  !!test 2
      umax= 0.5D0
      z=min(y,umax)
      rmCroitE =  CoefrmE * z*(1.D0-z)
   case(3)  !!test 3
      umax= 0.5D0
      z=min(y,umax)
      rmCroitE=  CoefrmE * (z*(1.D0-z))**2
   case(4) !!test 4
      rmCroitE = sqrt(CoefdiffeAdeg*2.*y)
   case default
      stop 'pb rmCroitE'
   end select
 end function rmCroitE
 !=======================
 function DerrmCroitE(y)
   !=======================
   USE longr
   REAL( kind=long) :: y, DerrmCroitE, CoefrmE,umax
   !
   CoefrmE = sqrt(CoefdiffeAdeg)
   select case(ChoixAdegE) 
   case(0) !!test 0
      DerrmCroitE = 0.D0
   case(1)  !!test 1
      umax= 0.5D0
      y = max(y, epsilon)
      if (y <=umax) then
         DerrmCroitE =  CoefrmE * (1-2.D0*y)/(2*sqrt( y*(1.D0-y)))
      else
         DerrmCroitE = 0.D0
      end if
   case(2)  !!test 2
      umax= 0.5D0
      if (y <=umax) then
         DerrmCroitE =  CoefrmE * (1-2.D0*y)
      else
         DerrmCroitE = 0.D0
      end if
   case(3)  !!test 3
      umax= 0.5D0
      if (y <=umax) then
         DerrmCroitE =  CoefrmE * 2.D0*y*(1.D0 -3.D0*y +2.D0*y**2)
      else
         DerrmCroitE = 0.D0
      end if
   case(4) !!test 4
      DerrmCroitE = CoefrmE/(sqrt(2.*y))
   case default
      stop 'pb DerrmCroitE'
   end select
 end function DerrmCroitE

 !=======================
 function rmDeCroitE(y)
   !=======================
   USE longr
   REAL( kind=long) :: y, z, rmDeCroitE, CoefrmE,umax
   !
   CoefrmE = sqrt(CoefdiffeAdeg)
   select case(ChoixAdegE)
   case(0) !!test 0
      rmDeCroitE = CoefrmE
   case(1)  !!test 1
      umax= 0.5D0
      z = max(y,umax)
      rmDeCroitE =  CoefrmE *sqrt( z*(1.D0-z)) -sqrt( umax*(1.D0-umax) )
   case(2)  !!test 2
      umax= 0.5D0
      z=max(y,umax)
      rmDeCroitE =  CoefrmE * ( z*(1.D0-z) - umax*(1.D0-umax) )
   case(3)  !!test 3
      umax= 0.5D0
      z=max(y,umax)
      rmDeCroitE=  CoefrmE * ( (z*(1.D0-z))**2 - (umax*(1.D0-umax))**2 )
   case(4) !!test 4
      rmDeCroitE = 0.D0
   case default
      stop 'pb rmDeCroitE'
   end select
 end function rmDeCroitE

 !=======================
 function DerrmDeCroitE(y)
   !=====================
   USE longr
   REAL( kind=long) :: y, DerrmDeCroitE, CoefrmE
   !
   CoefrmE = sqrt(CoefdiffeAdeg)
   select case(ChoixAdegE)
   case(0) !!test 0
      DerrmDeCroitE = 0.D0
   case(1)  !!test 1
      umax= 0.5D0
      y = min(y,1.0-epsilon)
      if (y > umax) then
         DerrmDeCroitE =  CoefrmE * (1-2.D0*y)/(2*sqrt( y*(1.D0-y)))
      else
         DerrmDeCroitE = 0.D0
      end if
   case(2)  !!test 2
      umax= 0.5D0
      if (y > umax) then
         DerrmDeCroitE =  CoefrmE * (1-2.D0*y)
      else
         DerrmDeCroitE = 0.D0
      end if
   case(3)  !!test 3
      umax= 0.5D0
      if (y > umax) then
         DerrmDeCroitE =  CoefrmE * 2.D0*y*(1.D0 -3.D0*y +2.D0*y**2)
      else
         DerrmDeroit = 0.D0
      end if
   case(4) !!test 4
      DerrmDeroit = 0.D0
   case default
      stop 'pb DerrmDeCroitE'
   end select
 end function DerrmDeCroitE


 !=======================
 function DerrmE(y)
   !=======================
   USE longr
   REAL( kind=long) :: y, DerrmE, CoefrmE
   !
   CoefrmE = sqrt(CoefdiffeAdeg)
   select case(ChoixAdegE)
 
   case(1)  !!test 1
      DerrmE = 0.D0
   case(2)  !!test 2
      DerrmE = CoefrmE *(1.D0 - 2.D0*y)
   case(3)  !!test 3
      DerrmE= CoefrmE * 2.D0*y*(1.D0 -3.D0*y +2.D0*y**2)
   case default
      stop 'pb DerrmE'
   end select
 end function DerrmE

!=======================
 function MuE(y)
   !=======================
   USE longr
   REAL( kind=long) :: y, MuE


   select case(ChoixMuE)
   case(1)
      MuE = CoefTranspMuE 
   case(2)
      MuE = CoefTranspMuE *y
   case(3)
      MuE= CoefTranspMuE *y*(1.D0-y) 
   case(4)
      MuE= CoefTranspMuE *y*(1.D0-y)**2
   case(5)
      IF ( y <= ubar ) THEN
         MuE= CoefTranspMuE *y*( 1.D0-(y/ubar)**gamma)
      ELSE
         MuE= 0.D0
      END IF
      case(6)
         MuE= CoefTranspMuE *( y*(1.D0-y) )**2
   case default
      stop 'pb MuE'
   end select

 end function MuE

 !=======================
 function MuCroitE(y)
   !=======================
   USE longr
   implicit none
   REAL( kind=long) :: y, z, MuCroitE, umax
   !integer :: choix



   select case(ChoixMuE)
   case(1)
      MuCroitE = CoefTranspMuE 
   case(2)
      MuCroitE = CoefTranspMuE *y
   case(3)
      umax= 0.5D0
      z=min(y,umax)
      MuCroitE=  CoefTranspMuE * z*(1.D0-z)
   case(4)
      umax= 1.D0/3.D0
      z=min(y,umax)
      MuCroitE= CoefTranspMuE * z*(1.D0-z)**2
   case(5)
      umax= (ubar)/( (gamma+1.D0)**(1.D0/gamma) )
      z=min(y,umax)
      IF ( z <= ubar ) THEN
         MuCroitE= CoefTranspMuE * z*(1.D0-(z/ubar)**gamma)
      ELSE
         MuCroitE= 0.D0
      END IF
   case(6)
      umax= 0.5D0
      z=min(y,umax)
      MuCroitE=  CoefTranspMuE * (z*(1.D0-z))**2
   case default
      stop 'pb MuCroitE'
   end select

 end function MuCroitE

 !=======================
 function DerivMuCroitE(y)
   !=======================
   USE longr
   REAL( kind=long) :: y, DerivMuCroitE, umax
   !integer :: choix
   !

   select case(ChoixMuE)
   case(1)
      DerivMuCroitE = 0.D0 
   case(2)
      DerivMuCroitE = CoefTranspMuE
   case(3)
      umax= 0.5D0
      if (y <=umax) then
         DerivMuCroitE =  CoefTranspMuE * (1-2.D0*y)
      else
         DerivMuCroitE = 0.D0
      end if
   case(4)
      umax= 1.D0/3.D0
      if (y < umax) then
         DerivMuCroitE = CoefTranspMuE * ( (1.D0-y)**2 -2.D0*y*(1.D0-y) )
      else
         DerivMuCroitE = 0.D0
      end if
   case(5)
      umax= (ubar)/( (gamma+1.D0)**(1.D0/gamma) )
      IF (y < umax .AND. y< ubar) THEN !! Although we know that y<=ubar
         DerivMuCroitE = CoefTranspMuE * (1.D0 - ( gamma+1.D0 )*( (y/ubar)**gamma) )
      ELSE
         DerivMuCroitE = 0.D0
      END IF
   case(6)
      umax= 0.5D0
      if (y <=umax) then
         DerivMuCroitE =  CoefTranspMuE * 2.D0*y*(1.D0-y)*(1-2.D0*y)
      else
         DerivMuCroitE = 0.D0
      end if
   case default
      stop 'pb DerivMuCroitE'
   end select

 end function DerivMuCroitE

 !=======================
 function MuDeCroitE(y)
   !=======================
   USE longr
   REAL( kind=long) :: y, z, MuDeCroitE, umax
   !

   select case(ChoixMuE)
   case(1)
      MuDeCroitE = 0.D0 
   case(2)
      MuDeCroitE = 0.D0
   case(3)
      umax= 0.5D0
      z=max(y,umax)
      MuDeCroitE=  CoefTranspMuE * ( z*(1.D0-z) - umax*(1.D0-umax) )
   case(4)
      umax= 1.D0/3.D0
      z=max(y,umax)
      MuDeCroitE= CoefTranspMuE * ( z*(1.D0-z)**2 -umax*(1.D0-umax)**2 ) 
   case(5)
      umax= (ubar)/( (gamma+1.D0)**(1.D0/gamma) )
      z=max(y,umax)
      IF (z <= ubar ) THEN
         MuDeCroitE= CoefTranspMuE * ( z*(1.D0-(z/ubar)**gamma) )
      ELSE
         MuDeCroitE= 0.D0
      END IF
      IF ( umax <= ubar )  THEN
         MuDeCroitE= MuDeCroitE  - CoefTranspMuE*umax*( 1.D0-(umax/ubar)**gamma ) 
      END IF
   case(6)
      umax= 0.5D0
      z=max(y,umax)
      MuDeCroitE=  CoefTranspMuE * ( (z*(1.D0-z))**2 - (umax*(1.D0-umax))**2 )
   case default
      stop 'pb MuCroit'
   end select

 end function MuDeCroitE

 !=======================
 function DerivMuDeCroitE(y)
   !=======================
   USE longr
   REAL( kind=long) :: y, DerivMuDeCroitE, umax
   !

   select case(ChoixMuE)
   case(1)
      DerivMuDeCroitE = 0.D0 
   case(2)
      DerivMuDeCroitE = 0.D0
   case(3)
      umax= 0.5D0
      if (y > umax) then
         DerivMuDeCroitE =  CoefTranspMuE * (1.D0-2.D0*y)
      else
         DerivMuDeCroitE = 0.D0
      end if
   case(4)
      umax= 1.D0/3.D0
      if (y > umax) then
         DerivMuDeCroitE = CoefTranspMuE * ( (1.D0-y)**2 -2.D0*y*(1.D0-y) )
      else
         DerivMuDeCroitE = 0.D0
      end if
   case(5)
      umax= (ubar)/( (gamma+1.D0)**(1.D0/gamma) )
      IF (y > umax .AND. y < ubar) THEN
         DerivMuDeCroitE = CoefTranspMuE* ( 1.D0-(gamma+1.D0)*( (y/ubar)**gamma ) )
      ELSE
         DerivMuDeCroitE = 0.D0
      END IF
   case(6)
      umax= 0.5D0
      if (y > umax) then
         DerivMuDeCroitE =  CoefTranspMuE * 2.D0*y*(1.D0-y)*(1.D0-2.D0*y)
      else
         DerivMuDeCroitE = 0.D0
      end if
   case default
      stop 'pb DerivMuDeCroitE'
   end select

 end function DerivMuDeCroitE


 !=======================
 function DerivMuE(y)
   !=======================
   USE longr
   REAL( kind=long) :: y,DerivMuE


   select case(ChoixMuE)
   case(1)
      DerivMuE = 0.D0
   case(2)
      DerivMuE = CoefTranspMuE *1.D0
   case(3)
      DerivMuE = CoefTranspMuE * ( 1.D0 - 2.D0*y )
   case(4)
      DerivMuE = CoefTranspMuE *( 1.D0 - y ) * ( 1.D0 - 3.D0*y )
   case(5)
      IF ( y < ubar ) THEN
         DerivMuE= CoefTranspMuE *( 1.D0-( gamma+ 1.D0 )*( (y/ubar)**gamma ) )
      ELSE
         DerivMuE= 0.D0
      END IF
   case(6)
      DerivMuE = CoefTranspMuE * 2.D0 *y*(1.D0-y)*( 1.D0 - 2.D0*y )
   case default
      stop 'pb DerivMuE'
   end select
 end function DerivMuE

END Module fsource
