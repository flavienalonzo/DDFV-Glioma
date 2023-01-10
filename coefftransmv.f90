!******************************************************************************
!     *-!!!!!!!!!@@@@@!!@@@@@!!@@@@@!!@@@@@!!!!!!!!!!!-*
!     *-!!!!!!!!!!!!!@!!@!!!@!!!!!!@!!@!!!@!!!!!!!!!!!-*
!     *-!!!!!!!!!@@@@@!!@!!!@!!@@@@@!!@!!!@!!!!!!!!!!!-*
!     *-!!!!!!!!!@!!!!!!@!!!@!!@@!!!!!@!!!@!!!!!!!!!!!-*
!     *-!!!!!!!!!@@@@@!!@@@@@!!@@@@@!!@@@@@!!!!!!!!!!!-*
!------------------------------------------------------------
!     * Ce sous programme calcule sur un segment donné 
!     * la matrice élémentaire locale pour l'aproximation DDFV
!     * MATLOC fourni une matrice 2x2
!     *---------------------------
!     * Déscription des paramètres
!     *---------------------------
!******************************************************************************
SUBROUTINE coefftransmv
  !--------
  ! Modules
  !--------
  USE longr
  USE parmmage
  !----------------------------------
  ! Déclaration des variables locales
  !----------------------------------
  IMPLICIT NONE
  CHARACTER(len=6)                  :: oldprf
  INTEGER                           :: iseg
  !-------------------
  ! Début du programme
  !-------------------
  oldprf = prefix
  prefix = 'TRANSM'
  !------
  ! Corps
  !------
  Allocate(vTKL(Nseg),vTKeLe(Nseg),vetaSSe(Nseg))
  DO iseg = 1,Nseg
     !
     ! les transmissibiltés de diffusion TKL
     !--------------------------------------
     !
     vTKL(iseg) = CoefDiffV*( (SxxKv(iseg)*NsigK(1,iseg)+SxyKv(iseg) *NsigK(2,iseg) )*NsigK(1,iseg)&
          & + ( SxyKv(iseg)*NsigK(1,iseg)+SyyKv(iseg)*NsigK(2,iseg) )*NsigK(2,iseg))/(2*AireD(iseg))
     !
     ! les transmissibiltés de diffusion TK*L*
     !----------------------------------------
     !
     vTKeLe(iseg) = CoefDiffV*( (SxxKv(iseg)*NsigeKe(1,iseg)+SxyKv(iseg) *NsigeKe(2,iseg) )*NsigeKe(1,iseg)&
          & + ( SxyKv(iseg)*NsigeKe(1,iseg)+SyyKv(iseg)*NsigeKe(2,iseg) )*NsigeKe(2,iseg))/(2*AireD(iseg))
     !
     ! les transmissibiltés de convection eta_sigma,sigma*
     !----------------------------------------------------
     !
     vetaSSe(iseg) = CoefDiffV*( (SxxKv(iseg)*NsigK(1,iseg)+SxyKv(iseg) *NsigK(2,iseg) )*NsigeKe(1,iseg)&
          & + ( SxyKv(iseg)*NsigK(1,iseg)+SyyKv(iseg)*NsigK(2,iseg) )*NsigeKe(2,iseg))/(2*AireD(iseg))
     !
  END DO
  !------------
  ! Impressions
  !------------
  IF (iprint >= 5) THEN
     write(uprint,*)'iseg, vTKL(iseg), vTKeLe(iseg), vetaSSe(iseg)'
     DO iseg = 1,Nseg
        WRITE (uprint,110) iseg, vTKL(iseg), vTKeLe(iseg), vetaSSe(iseg)
     END DO
  END IF
110 FORMAT (i8,3(E16.9, 2X))
  !-----------------
  ! Fin du programme
  !-----------------
  prefix = oldprf

END SUBROUTINE coefftransmv
