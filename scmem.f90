!            **********************
!            **  SUBROUTINE SCMEM**
!            **********************
!******************************************************************************
!     *-!!!!!!!!!@@@@@!!@@@@@!!@@@@@!!@@@@@!!!!!!!!!!!-*
!     *-!!!!!!!!!!!!!@!!@!!!@!!!!!!@!!@!!!@!!!!!!!!!!!-*
!     *-!!!!!!!!!@@@@@!!@!!!@!!@@@@@!!@!!!@!!!!!!!!!!!-*
!     *-!!!!!!!!!@!!!!!!@!!!@!!@@!!!!!@!!!@!!!!!!!!!!!-*
!     *-!!!!!!!!!@@@@@!!@@@@@!!@@@@@!!@@@@@!!!!!!!!!!!-* 
!------------------------------------------------------------
!     * Ce  sous programme calcul le second membre physique
!     * - D^2 u = F
!------------------------------------------------------------
!     * Utilise les modules longr, imprime, parmmmage
!******************************************************************************
SUBROUTINE SCMEM(A, choixf, temps)
  !--------
  ! Modules
  !--------
  USE longr
  USE parmmage
  USE fsource
  IMPLICIT NONE

  !--------------------------
  ! Declaration des arguments
  !--------------------------

  TYPE(MatCreux)                        :: A
  integer, intent(in)                   :: choixf
  REAL(kind=long), intent(in)           :: temps
  !----------------------------------
  ! Declaration des variables locales
  !----------------------------------
  CHARACTER(len=6)      :: oldprf
  integer               :: is,jt

  !-------------------
  ! Debut du programme
  !-------------------
  oldprf = prefix
  prefix = 'SCMEM '

  !------------------------------------------------------
  !  A%F   : est un tableau de taille NbInc = NsInt + Nbt  
  !        : NsInt le nombre des sommets à l'intérieur
  !        : Nbt le nombre des triangles
  !------------------------------------------------------
  
  A%F = 0.D0 
  do is = 1, NsInt
     A%F (is) = AireDsommet(is)*fsourceu( temps,CoordS(1,is), CoordS(2,is), choixf) 
  enddo
 do jt = 1, Nbt
     A%F (jt + NsInt) = AireK(jt)*fsourceu( temps,CoordK(1,jt), CoordK(2,jt), choixf ) 
  enddo
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

  RETURN
END SUBROUTINE SCMEM








