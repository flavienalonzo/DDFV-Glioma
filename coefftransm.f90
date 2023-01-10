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
SUBROUTINE coefftransm
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
  Allocate(TKL(Nseg),TKeLe(Nseg),eta(Nseg))
  DO iseg = 1,Nseg
     !
     ! les transmissibiltés de diffusion TKL
     !--------------------------------------
     !
     TKL(iseg) = ( (SxxK(iseg)*NsigK(1,iseg)+SxyK(iseg) *NsigK(2,iseg) )*NsigK(1,iseg)&
          & + ( SxyK(iseg)*NsigK(1,iseg)+SyyK(iseg)*NsigK(2,iseg) )*NsigK(2,iseg))/(2*AireD(iseg))
     !
     ! les transmissibiltés de diffusion TK*L*
     !----------------------------------------
     !
     TKeLe(iseg) = ( (SxxK(iseg)*NsigeKe(1,iseg)+SxyK(iseg) *NsigeKe(2,iseg) )*NsigeKe(1,iseg)&
          & + ( SxyK(iseg)*NsigeKe(1,iseg)+SyyK(iseg)*NsigeKe(2,iseg) )*NsigeKe(2,iseg))/(2*AireD(iseg))
     !
     ! les transmissibiltés de convection eta_sigma,sigma*
     !----------------------------------------------------
     !
     eta(iseg) = ( (SxxK(iseg)*NsigK(1,iseg)+SxyK(iseg) *NsigK(2,iseg) )*NsigeKe(1,iseg)&
          & + ( SxyK(iseg)*NsigK(1,iseg)+SyyK(iseg)*NsigK(2,iseg) )*NsigeKe(2,iseg))/(2*AireD(iseg))
     !
  END DO
  !------------
  ! Impressions
  !------------
  IF (iprint >= 5) THEN
     write(uprint,*)'iseg, TKL(iseg), TKeLe(iseg), eta(iseg)'
     DO iseg = 1,Nseg
        WRITE (uprint,110) iseg, TKL(iseg), TKeLe(iseg), eta(iseg)
     END DO
  END IF
110 FORMAT (i8,3(E16.9, 2X))
  !-----------------
  ! Fin du programme
  !-----------------
  prefix = oldprf

END SUBROUTINE coefftransm
