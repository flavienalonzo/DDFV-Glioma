Module fsource

Contains

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

      Tu = -u 
   end function

   function deriveeTu(temps,u)
      use longr 
      implicit none 
      REAL(kind=long), INTENT(in)     :: temps,u
      REAL(kind=long)  :: deriveeTu 

      deriveeTu = -1
   end function

   function hu(c)
      use longr
      implicit none
      REAL(kind=long), INTENT(in)     :: c
      REAL(kind=long)  :: hu 

      hu = c 
   end function

   function gv(c)
      use longr
      implicit none
      REAL(kind=long), INTENT(in)     :: c
      REAL(kind=long)  :: gv

      gv = c 
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
  function Adegen(y) !Phi=FU
    !=======================
    USE longr
    Real(kind=long) :: y,Adegen

    select case (choixAdeg)
    case(0)  !!test 0
       Adegen= CoefdiffuAdeg * y
    case(1)  !!test 1
       Adegen=CoefdiffuAdeg* (1.D0/2.D0 - y/3.D0 )*y**2
    case(2) !!test 2
       Adegen=CoefdiffuAdeg* (1.D0/3.D0 -y/2.D0 +y*y/5.D0 )*y**3
    case(3) !!test 3
       Adegen=CoefdiffuAdeg* (1.D0/5.D0 -2.D0*y/3.D0 + (6./7.D0)*y**2 -(y**3)/2.D0 +(y**4)/9.D0)*y**5
    case(4) !!test 4
       Adegen=CoefdiffuAdeg*y**2
    case default   
       stop 'Adege'
    end select

  end function Adegen

  !=============================
  function DerivAdeg(y) !dPhi=fu
    !===========================
    USE longr
    REAL( kind=long) :: y,DerAdeg
    !
    select case (choixAdeg)
    case(0) !!test 0
       DerivAdeg = CoefdiffuAdeg
    case(1) !!test 1
       DerivAdeg = CoefdiffuAdeg * y*( 1.D0 - y )
    case(2) !!test 2
       DerivAdeg= CoefdiffuAdeg * ( y*(1.D0 - y) )**2
    case(3) !!test 3
       DerivAdeg= CoefdiffuAdeg * ( y*(1.D0 - y) )**4
    case(4) !!test 4
       DerivAdeg=CoefdiffuAdeg*2.*y
    case default  
       stop 'DerAdeg'

    end select
  end function DerivAdeg

  !==========================================================
  function xi(y) ! Squared-Kirchoff transform: integral de rm
    !========================================================
    USE longr
    REAL( kind=long) :: y, xi, Coefrm
    !
    Coefrm = sqrt(CoefdiffuAdeg)
    select case (choixAdeg)
    case(0) !!test 0
       xi = Coefrm* y
    case(1)  !!test 1
       xi = Coefrm*( ((2*y-1)*sqrt(y*(1.D0 - y))/4) + (acos(-2*y+1.))/8)
    case(2) !!test 2
       xi = Coefrm * (1.D0/2.D0 -y/3.D0 )*y**2
    case(3) !!test 3
       xi = Coefrm * (1.D0/3.D0 -y/2.D0 +y*y/5.D0 )*y**3
    case(4) !!test 4
       xi=sqrt(CoefdiffuAdeg*2.)*(2./3.)*y**(3./2.)
    case default  
       stop 'xi'

    end select
  end function xi

  !=======================================
  function rm(y) ! racine carre de DerAdeg
    !=====================================
    USE longr
    REAL( kind=long) :: y, rm, Coefrm
    !
    Coefrm = sqrt(CoefdiffuAdeg)
    select case (choixAdeg)
    case(0) !!test 0
       rm = Coefrm
    case(1)  !!test 1
       rm = Coefrm *sqrt(y*(1.D0 - y))
    case(2) !!test 2
       rm = Coefrm * y*(1.D0 - y)
    case(3) !!test 3
       rm = Coefrm * ( y*(1.D0 - y) )**2
    case(4) !!test 4
       rm = sqrt(CoefdiffuAdeg*2.*y)
    case default  
       stop 'rm'

    end select
  end function rm

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
  function rmCroit(y)
    !=======================
    USE longr
    REAL( kind=long) :: y, z, rmCroit, Coefrm,umax
    !
    Coefrm = sqrt(CoefdiffuAdeg)
    select case(choixAdeg)
    case(0) !!test 0
       rmCroit = Coefrm
    case(1)  !!test 1
       umax= 0.5D0
       z=min(y,umax)
       rmCroit =  Coefrm *sqrt( z*(1.D0-z))
    case(2)  !!test 2
       umax= 0.5D0
       z=min(y,umax)
       rmCroit =  Coefrm * z*(1.D0-z)
    case(3)  !!test 3
       umax= 0.5D0
       z=min(y,umax)
       rmCroit=  Coefrm * (z*(1.D0-z))**2
    case(4) !!test 4
       rmCroit = sqrt(CoefdiffuAdeg*2.*y)
    case default
       stop 'pb rmCroit'
    end select
  end function rmCroit
  !=======================
  function DerrmCroit(y)
    !=======================
    USE longr
    REAL( kind=long) :: y, DerrmCroit, Coefrm,umax
    !
    Coefrm = sqrt(CoefdiffuAdeg)
    select case(choixAdeg) 
    case(0) !!test 0
       DerrmCroit = 0.D0
    case(1)  !!test 1
       umax= 0.5D0
       y = max(y, epsilon)
       if (y <=umax) then
          DerrmCroit =  Coefrm * (1-2.D0*y)/(2*sqrt( y*(1.D0-y)))
       else
          DerrmCroit = 0.D0
       end if
    case(2)  !!test 2
       umax= 0.5D0
       if (y <=umax) then
          DerrmCroit =  Coefrm * (1-2.D0*y)
       else
          DerrmCroit = 0.D0
       end if
    case(3)  !!test 3
       umax= 0.5D0
       if (y <=umax) then
          DerrmCroit =  Coefrm * 2.D0*y*(1.D0 -3.D0*y +2.D0*y**2)
       else
          DerrmCroit = 0.D0
       end if
    case(4) !!test 4
       DerrmCroit = Coefrm/(sqrt(2.*y))
    case default
       stop 'pb DerrmCroit'
    end select
  end function DerrmCroit

  !=======================
  function rmDeCroit(y)
    !=======================
    USE longr
    REAL( kind=long) :: y, z, rmDeCroit, Coefrm,umax
    !
    Coefrm = sqrt(CoefdiffuAdeg)
    select case(choixAdeg)
    case(0) !!test 0
       rmDeCroit = Coefrm
    case(1)  !!test 1
       umax= 0.5D0
       z = max(y,umax)
       rmDeCroit =  Coefrm *sqrt( z*(1.D0-z)) -sqrt( umax*(1.D0-umax) )
    case(2)  !!test 2
       umax= 0.5D0
       z=max(y,umax)
       rmDeCroit =  Coefrm * ( z*(1.D0-z) - umax*(1.D0-umax) )
    case(3)  !!test 3
       umax= 0.5D0
       z=max(y,umax)
       rmDeCroit=  Coefrm * ( (z*(1.D0-z))**2 - (umax*(1.D0-umax))**2 )
    case(4) !!test 4
       rmDeCroit = 0.D0
    case default
       stop 'pb rmDeCroit'
    end select
  end function rmDeCroit

  !=======================
  function DerrmDeCroit(y)
    !=====================
    USE longr
    REAL( kind=long) :: y, DerrmDeCroit, Coefrm
    !
    Coefrm = sqrt(CoefdiffuAdeg)
    select case(choixAdeg)
    case(0) !!test 0
       DerrmDeCroit = 0.D0
    case(1)  !!test 1
       umax= 0.5D0
       y = min(y,1.0-epsilon)
       if (y > umax) then
          DerrmDeCroit =  Coefrm * (1-2.D0*y)/(2*sqrt( y*(1.D0-y)))
       else
          DerrmDeCroit = 0.D0
       end if
    case(2)  !!test 2
       umax= 0.5D0
       if (y > umax) then
          DerrmDeCroit =  Coefrm * (1-2.D0*y)
       else
          DerrmDeCroit = 0.D0
       end if
    case(3)  !!test 3
       umax= 0.5D0
       if (y > umax) then
          DerrmDecroit =  Coefrm * 2.D0*y*(1.D0 -3.D0*y +2.D0*y**2)
       else
          DerrmDeroit = 0.D0
       end if
    case(4) !!test 4
       DerrmDeroit = 0.D0
    case default
       stop 'pb DerrmDeCroit'
    end select
  end function DerrmDeCroit


  !=======================
  function Derrm(y)
    !=======================
    USE longr
    REAL( kind=long) :: y, Derrm, Coefrm
    !
    Coefrm = sqrt(CoefdiffuAdeg)
    select case(choixAdeg)
  
    case(1)  !!test 1
       Derrm = 0.D0
    case(2)  !!test 2
       Derrm = Coefrm *(1.D0 - 2.D0*y)
    case(3)  !!test 3
       Derrm= Coefrm * 2.D0*y*(1.D0 -3.D0*y +2.D0*y**2)
    case default
       stop 'pb Derrm'
    end select
  end function Derrm

!=======================
  function Mu(y)
    !=======================
    USE longr
    REAL( kind=long) :: y, Mu


    select case(ChoixMu)
    case(1)
       Mu = CoefTranspMu 
    case(2)
       Mu = CoefTranspMu *y
    case(3)
       Mu= CoefTranspMu *y*(1.D0-y) 
    case(4)
       Mu= CoefTranspMu *y*(1.D0-y)**2
    case(5)
       IF ( y <= ubar ) THEN
          Mu= CoefTranspMu *y*( 1.D0-(y/ubar)**gamma)
       ELSE
          Mu= 0.D0
       END IF
       case(6)
          Mu= CoefTranspMu *( y*(1.D0-y) )**2
    case default
       stop 'pb Mu'
    end select

  end function Mu

  !=======================
  function MuCroit(y)
    !=======================
    USE longr
    implicit none
    REAL( kind=long) :: y, z, MuCroit, umax
    !integer :: choix



    select case(ChoixMu)
    case(1)
       MuCroit = CoefTranspMu 
    case(2)
       MuCroit = CoefTranspMu *y
    case(3)
       umax= 0.5D0
       z=min(y,umax)
       MuCroit=  CoefTranspMu * z*(1.D0-z)
    case(4)
       umax= 1.D0/3.D0
       z=min(y,umax)
       MuCroit= CoefTranspMu * z*(1.D0-z)**2
    case(5)
       umax= (ubar)/( (gamma+1.D0)**(1.D0/gamma) )
       z=min(y,umax)
       IF ( z <= ubar ) THEN
          MuCroit= CoefTranspMu * z*(1.D0-(z/ubar)**gamma)
       ELSE
          MuCroit= 0.D0
       END IF
    case(6)
       umax= 0.5D0
       z=min(y,umax)
       MuCroit=  CoefTranspMu * (z*(1.D0-z))**2
    case default
       stop 'pb MuCroit'
    end select

  end function MuCroit

  !=======================
  function DerivMuCroit(y)
    !=======================
    USE longr
    REAL( kind=long) :: y, DerivMuCroit, umax
    !integer :: choix
    !

    select case(ChoixMu)
    case(1)
       DerivMuCroit = 0.D0 
    case(2)
       DerivMuCroit = CoefTranspMu
    case(3)
       umax= 0.5D0
       if (y <=umax) then
          DerivMuCroit =  CoefTranspMu * (1-2.D0*y)
       else
          DerivMuCroit = 0.D0
       end if
    case(4)
       umax= 1.D0/3.D0
       if (y < umax) then
          DerivMuCroit = CoefTranspMu * ( (1.D0-y)**2 -2.D0*y*(1.D0-y) )
       else
          DerivMuCroit = 0.D0
       end if
    case(5)
       umax= (ubar)/( (gamma+1.D0)**(1.D0/gamma) )
       IF (y < umax .AND. y< ubar) THEN !! Although we know that y<=ubar
          DerivMuCroit = CoefTranspMu * (1.D0 - ( gamma+1.D0 )*( (y/ubar)**gamma) )
       ELSE
          DerivMuCroit = 0.D0
       END IF
    case(6)
       umax= 0.5D0
       if (y <=umax) then
          DerivMuCroit =  CoefTranspMu * 2.D0*y*(1.D0-y)*(1-2.D0*y)
       else
          DerivMuCroit = 0.D0
       end if
    case default
       stop 'pb DerivMuCroit'
    end select

  end function DerivMuCroit

  !=======================
  function MuDeCroit(y)
    !=======================
    USE longr
    REAL( kind=long) :: y, z, MuDecroit, umax
    !

    select case(ChoixMu)
    case(1)
       MuDeCroit = 0.D0 
    case(2)
       MuDecroit = 0.D0
    case(3)
       umax= 0.5D0
       z=max(y,umax)
       MuDeCroit=  CoefTranspMu * ( z*(1.D0-z) - umax*(1.D0-umax) )
    case(4)
       umax= 1.D0/3.D0
       z=max(y,umax)
       MuDecroit= CoefTranspMu * ( z*(1.D0-z)**2 -umax*(1.D0-umax)**2 ) 
    case(5)
       umax= (ubar)/( (gamma+1.D0)**(1.D0/gamma) )
       z=max(y,umax)
       IF (z <= ubar ) THEN
          MuDeCroit= CoefTranspMu * ( z*(1.D0-(z/ubar)**gamma) )
       ELSE
          MuDeCroit= 0.D0
       END IF
       IF ( umax <= ubar )  THEN
          MuDeCroit= MuDeCroit  - CoefTranspMu*umax*( 1.D0-(umax/ubar)**gamma ) 
       END IF
    case(6)
       umax= 0.5D0
       z=max(y,umax)
       MuDeCroit=  CoefTranspMu * ( (z*(1.D0-z))**2 - (umax*(1.D0-umax))**2 )
    case default
       stop 'pb MuCroit'
    end select

  end function MuDeCroit

  !=======================
  function DerivMuDecroit(y)
    !=======================
    USE longr
    REAL( kind=long) :: y, DerivMuDecroit, umax
    !

    select case(ChoixMu)
    case(1)
       DerivMuDecroit = 0.D0 
    case(2)
       DerivMuDecroit = 0.D0
    case(3)
       umax= 0.5D0
       if (y > umax) then
          DerivMuDecroit =  CoefTranspMu * (1.D0-2.D0*y)
       else
          DerivMuDecroit = 0.D0
       end if
    case(4)
       umax= 1.D0/3.D0
       if (y > umax) then
          DerivMuDecroit = CoefTranspMu * ( (1.D0-y)**2 -2.D0*y*(1.D0-y) )
       else
          DerivMuDecroit = 0.D0
       end if
    case(5)
       umax= (ubar)/( (gamma+1.D0)**(1.D0/gamma) )
       IF (y > umax .AND. y < ubar) THEN
          DerivMuDeCroit = CoefTranspMu* ( 1.D0-(gamma+1.D0)*( (y/ubar)**gamma ) )
       ELSE
          DerivMuDeCroit = 0.D0
       END IF
    case(6)
       umax= 0.5D0
       if (y > umax) then
          DerivMuDecroit =  CoefTranspMu * 2.D0*y*(1.D0-y)*(1.D0-2.D0*y)
       else
          DerivMuDecroit = 0.D0
       end if
    case default
       stop 'pb DerivMuDecroit'
    end select

  end function DerivMuDecroit


  !=======================
  function DerivMu(y)
    !=======================
    USE longr
    REAL( kind=long) :: y,DerivMu


    select case(ChoixMu)
    case(1)
       DerivMu = 0.D0
    case(2)
       DerivMu = CoefTranspMu *1.D0
    case(3)
       DerivMu = CoefTranspMu * ( 1.D0 - 2.D0*y )
    case(4)
       DerivMu = CoefTranspMu *( 1.D0 - y ) * ( 1.D0 - 3.D0*y )
    case(5)
       IF ( y < ubar ) THEN
          DerivMu= CoefTranspMu *( 1.D0-( gamma+ 1.D0 )*( (y/ubar)**gamma ) )
       ELSE
          DerivMu= 0.D0
       END IF
    case(6)
       DerivMu = CoefTranspMu * 2.D0 *y*(1.D0-y)*( 1.D0 - 2.D0*y )
    case default
       stop 'pb derivMu'
    end select


  end function DerivMu


END Module fsource
