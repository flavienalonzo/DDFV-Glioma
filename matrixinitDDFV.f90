!            *********************************
!            **  SUBROUTINE MATRIXINITDDFV  **
!            *********************************
!******************************************************************************
!     *-!!!!!!!!!@@@@@!!@@@@@!!@@@@@!!@@@@@!!!!!!!!!!!-*
!     *-!!!!!!!!!!!!!@!!@!!!@!!!!!!@!!@!!!@!!!!!!!!!!!-*
!     *-!!!!!!!!!@@@@@!!@!!!@!!@@@@@!!@!!!@!!!!!!!!!!!-*
!     *-!!!!!!!!!@!!!!!!@!!!@!!@@!!!!!@!!!@!!!!!!!!!!!-*
!     *-!!!!!!!!!@@@@@!!@@@@@!!@@@@@!!@@@@@!!!!!!!!!!!-*
!------------------------------------------------------------
!     * Ce sous programme donne la structure generale de 
!     * de la matrice A%TMAT ainsi le remplissage de 
!     * A%IndPL  
!     * A%Indc
!     * la matrice elementaire locale Aloc du triangle jt
!     * Ajouter a la ligne I et la colone J la valeur Coef
!     *----------------------------
!     * Fonctionnement du programme
!     *----------------------------
!     * Utilise les modules longr, imprime, parmmmage
!****************************************************************************
SUBROUTINE matrixinitDDFV(Mat)
  !--------
  ! Modules
  !--------
  USE longr
  USE parmmage
  USE imprime
  IMPLICIT NONE

  !--------------------------
  ! Declaration des arguments
  !--------------------------
  TYPE(MatCreux),         INTENT(OUT) :: Mat
  !----------------------------------
  ! Declaration des variables locales
  !----------------------------------
  CHARACTER(len=6)                    :: oldprf
  INTEGER                             :: iseg, is, js, i1, jv, kv,jmin,mloc
  INTEGER                             :: j,jt,iK,jL,NcoefMat
  INTEGER, DIMENSION(NbInc)           :: IndPL

  !-------------------
  ! Debut du programme
  !-------------------
  oldprf = prefix
  prefix = 'MATINIT'
  !------
  ! Corps
  !------
  ! Les termes diagonaux
  NcoefMat = NbInc
  ! contribution des triangles entre eux
  DO iseg = 1, Nseg
     IF(NtypSeg(iseg) == 0)   NcoefMat = NcoefMat + 2
  END DO
  ! contribution des centres triangles et sommets
  DO jt = 1, Nbt
     DO is = 1,Typegeomaille
        select case (ChoixCdtBord)
        case (Dirichlet)
           IF ( NuSoK(is,jt)  <= NsInt )  NcoefMat = NcoefMat + 2
        case (Neumann)
           NcoefMat = NcoefMat + 2
        END select
     END DO
  END DO
  ! contribution des sommets entre eux
  DO iseg = 1, Nseg
     is = NuSeg(1,iseg) ;   js = NuSeg(2,iseg)
     select case (ChoixCdtBord)
     case (Dirichlet)
        IF(is > NsInt .OR. js > NsInt) CYCLE
        NcoefMat = NcoefMat + 2
     case (Neumann)
        NcoefMat = NcoefMat + 2
     END select
  END DO
  !print*,'NcoefMat =', NcoefMat
  !!

  ALLOCATE( Mat%IndPL(1:NbInc+1), Mat%Indc(NcoefMat), Mat%TMat(NcoefMat) )
  ALLOCATE( Mat%F ( NbInc ), Mat%Bg ( NbInc ) )
  !
  Mat%TMat = 0.0D0 ;  Mat%F = 0.D0 ;  Mat%Bg = 0.D0; Mat%Indc = 0.D0
  !! La diagonale
  DO is = 1, NbInc+1
     Mat%IndPL(is) = is
  END DO
  !! KeLe
  DO iseg = 1, Nseg
     is = NuSeg(1,iseg) ; js = NuSeg(2, iseg)
     select case (ChoixCdtBord)
     case (Dirichlet)
        IF(is > NsInt .OR. js > NsInt) CYCLE
        Mat%IndPL(is+1:NbInc+1) = Mat%IndPL(is+1:NbInc+1) + 1
        Mat%IndPL(js+1:NbInc+1) = Mat%IndPL(js+1:NbInc+1) + 1
     case (Neumann)
           Mat%IndPL(is+1:NbInc+1) = Mat%IndPL(is+1:NbInc+1) + 1
           Mat%IndPL(js+1:NbInc+1) = Mat%IndPL(js+1:NbInc+1) + 1
     END select
  END DO
  !! KKe et KeK
  DO jt = 1, Nbt
     iK = jt + NsInt
     DO j = 1,Typegeomaille
        is = NuSoK(j,jt)
        IF (is <= NsInt )  THEN 
           Mat%IndPL(is+1:NbInc+1) = Mat%IndPL(is+1:NbInc+1) + 1
           Mat%IndPL(iK+1:NbInc+1) = Mat%IndPL(iK+1:NbInc+1) + 1
        END IF
     END DO
  END DO
  !! KL
  DO iseg = 1, Nseg
     IF (NtypSeg(iseg)==0) then 
        iK = NumTVoisSeg(1,iseg) + NsInt
        jL = NumTVoisSeg(2,iseg) + NsInt
        Mat%IndPL(iK+1:NbInc+1) = Mat%IndPL(iK+1:NbInc+1) + 1
        Mat%IndPL(jL+1:NbInc+1) = Mat%IndPL(jL+1:NbInc+1) + 1
     END IF
  END DO
  !print*,'IndPL =', Mat%IndPL

  !! Calcul des Indices des colonnes 
  IndPL(1:NbInc) = Mat%IndPL(1:NbInc)
  !! La diagonale
  DO is = 1, NbInc
     Mat%Indc( IndPL(is) ) = is
     IndPL(is)  = IndPL(is) + 1
  ENDDO
  !! KeLe
  DO iseg = 1, Nseg
     is = NuSeg(1,iseg) ;  js = NuSeg(2, iseg)
     !
     select case (ChoixCdtBord)
     case (Dirichlet)
        IF (is > NsInt .OR. js > NsInt) CYCLE
        Mat%Indc(IndPL(is)) = js
        Mat%Indc(IndPL(js)) = is
        IndPL(is)           = IndPL(is) + 1
        IndPL(js)           = IndPL(js) + 1
     case(Neumann)
           Mat%Indc(IndPL(is)) = js
           Mat%Indc(IndPL(js)) = is
           IndPL(is)           = IndPL(is) + 1
           IndPL(js)           = IndPL(js) + 1
     END select
  END DO
  !! KKe et KeK
  DO jt = 1, Nbt
     DO j = 1,Typegeomaille
        is = NuSoK(j,jt)
        select case (ChoixCdtBord)
        case (Dirichlet)
           IF (is  <= NsInt ) THEN 
              Mat%Indc(IndPL(jt+NsInt)) = is
              IndPL(jt+NsInt)           = IndPL(jt+NsInt) + 1
              Mat%Indc(IndPL(is))       = jt+NsInt
              IndPL(is)                 = IndPL(is) + 1
           END IF
        case (Neumann)
           Mat%Indc(IndPL(jt+Nbs))   = is
           IndPL(jt+Nbs)             = IndPL(jt+Nbs) + 1
           Mat%Indc(IndPL(is))       = jt+Nbs
           IndPL(is)                 = IndPL(is) + 1
        END select
     END DO
  END DO
  !!KL
  DO iseg = 1, Nseg
     IF (NtypSeg(iseg)==0) then 
        iK = NumTVoisSeg(1,iseg) + Nsint
        jL = NumTVoisSeg(2,iseg) + Nsint
        Mat%Indc(IndPL(iK)) = jL
        Mat%Indc(IndPL(jL)) = iK
        IndPL(iK)           = IndPL(iK) + 1
        IndPL(jL)           = IndPL(jL) + 1
     end IF
  END DO
  ! mettre les indices de colonnes par ordre croissant

  DO is = 1, NbInc 
     DO jv = Mat%IndPL(is), Mat%IndPL(is+1)-1 
        jmin = jv
        do   kv = jv +1,  Mat%IndPL(is+1)-1
           IF (Mat%Indc(kv) < Mat%Indc(jmin)) jmin = kv
        END DO
        mloc =  Mat%Indc(jmin )
        Mat%Indc(jmin ) = Mat%Indc(jv ) 
        Mat%Indc(jv ) = mloc 
     END DO
  END DO
  !------------
  ! Impressions
  !------------
  IF (iprint >=5) THEN
     CALL prvari(uprint,'NcoefMat = ', NcoefMat )
     WRITE(uprint,*) ( Mat%Tmat(is) , is=1, NcoefMat )
     WRITE(uprint,*) ( Mat%Indc(is) , is=1, NcoefMat )
     WRITE(uprint,*) ( Mat%IndPL(is), is=1, NbInc+1 )
  ENDIF
  !-----------------
  ! Fin du programme
  !-----------------
  prefix = oldprf

END SUBROUTINE matrixinitDDFV


  
