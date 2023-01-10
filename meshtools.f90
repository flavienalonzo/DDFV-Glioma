SUBROUTINE  meshtools
  !------------------------------------------------
  ! Nbs   : nombre des sommets total
  ! Nbt   : nombre des triangle
  ! NsInt : nombre des sommets interieur au domaine 
  !        nombre des sommets ou la solutions est calculee
  ! Nbord : nombre des sommets au bord
  ! NsInt = Nbs - Nbord
  ! NuMSeg : numero arete triangle, numero milieu segment
  ! CoordMseg : coordonnee du point milieu du segment
  !----------------------------------
  ! Declaration des variables locales
  !----------------------------------
  USE longr
  USE imprime
  USE parmmage
  IMPLICIT NONE 
  CHARACTER(len=6)      :: oldprf
  INTEGER               :: n, ntyp, k, ks,Kvois,Lvois
  INTEGER               :: isegmin, isegmax, jt, iseg, next,typeseg
  INTEGER               :: m, l,ii, ik, jL, is, js , ls, i , j , kt, ncentrevernoi
  INTEGER               :: Nbstriang, nusommet1, nusommet2, nextint,nextext

  REAL(kind=long) , DIMENSION(:,:), ALLOCATABLE :: CoordSbis
  Integer , DIMENSION(:,:), ALLOCATABLE         :: NuSegbis, NuSoKbis, NumTVoisSegbis 
  INTEGER, DIMENSION(:), ALLOCATABLE            :: PermutSegment,PermutSommet,ntypsbis,NTypSegbis

  real(kind=long), dimension(4) :: x, y
  real(kind=long) :: x1,y1,z
  CHARACTER(len=14)      :: chaine

  !-------------------
  ! Debut du programme
  !-------------------
  oldprf = prefix
  prefix = 'TOOLS'

  ! type geometrie de la maille (triangle)
  !! maillage triangulaire
  ! ----------------------------------------------------------
  ! liste des segments NtypSeg, NombVoisSeg
  ! calcul des aires des triangles
  ! ----------------------------------------------------------
  IF (Maillage_Format == 1) THEN
     DO iseg = 1,Nseg
        IF (NombVoisSeg(iseg) == 2) THEN
           NTypSeg(iseg) = 0
        ELSE
           NTypSeg(iseg) = ChoixCdtBord
        END IF
     END DO
  ELSE
     ALLOCATE( NombVoisSeg(NSeg))
     ALLOCATE (NTypSeg(Nseg))
     NTypSeg = ChoixCdtBord
     NombVoisSeg = 1
     DO iseg = 1,Nseg
        IF (NumTVoisSeg(2,iseg) /= 0) THEN
           NTypSeg(iseg) = 0
           NombVoisSeg(iseg) = 2
        END IF
     END DO
  END IF
  !! =============================================
  !! calcul de Ntyps dans le cas de maillage boyer 
  !! =============================================
  Allocate (Ntyps(Nbs) )  
  Ntyps = 0 
  DO iseg = 1,Nseg
     is = NuSeg(1,iseg) ; js = NuSeg(2,iseg)
     IF  (NTypSeg(iseg)/=0) then 
        Ntyps(is) = ChoixCdtBord ; Ntyps(js) = ChoixCdtBord
     Endif
  END DO
    
  IF (Maillage_Format == 1) THEN
     Nbord = 0
     DO is = 1,Nbs
        IF (Ntyps(is) == ChoixCdtBord) Nbord = Nbord + 1
     END DO
  END IF
  !-----------------------------------------------------------------
  ! calcul de NsInt le nombre de sommets ou la solution est calculee
  !-----------------------------------------------------------------
  NsInt = Nbs - Nbord
  write(*,*)'le nombre de sommets interieurs est:', NsInt
  write(*,*)'le nombre de sommets exterieurs est:', Nbord
  write(*,*)'le nombre de sommets total est:', Nbs

  CALL prvari(uprint,'nombre des sommets interieurs  (NsInt) ', NsInt )

  CALL prvari(uprint,'nombre des sommets sur le bord (Nbord) ', Nbord )
  !
  !! Calcule de Nsegint et Nsegbord
  !
  Nsegint = 0
  Do iseg =1, Nseg
     if(NtypSeg(iseg) == 0) Then
        Nsegint = Nsegint + 1
     end if
  end do

  Nsegbord= Nseg - Nsegint
  write(*,*)'le nombre de segments interieurs est:', Nsegint
  write(*,*)'le nombre de segments exterieurs est:', Nsegbord
  write(*,*)'le nombre de segments total est:', Nseg
  ! -------------------------------------- 
  ! -- reorientation positive des diamants
  ! --------------------------------------
  Do iseg = 1,Nseg
     is = NuSeg(1,iseg); js = NuSeg(2,iseg)
     iK = NumTVoisSeg(1,iseg)
     x(1) = coordS(1,is) ; y(1) = coordS(2,is)
     x(2) = coordS(1,js) ; y(2) = coordS(2,js)
     x(3) = coordK(1,iK) ; y(3) = coordK(2,iK)
     !
     IF ((x(2)-x(1))*(y(3)-y(1)) - (y(2)-y(1))*(x(3)-x(1)) < 0.D0) THEN
        !! ici on echange les sommets du segment si l'orientation du diamant n'est pas positive
        NuSeg(1,iseg) = js
        NuSeg(2,iseg) = is
     END IF
  END Do
  !!---------------------------------------------------
  !! renumerotation des sommets interieurs de 1 à NsInt
  !! --------------------------------------------------
  Allocate(PermutSommet(Nbs))
  Nextint = 1 ; Nextext = 1
  DO is = 1, Nbs
     IF( Ntyps(is) == 0) then 
        PermutSommet(is) = Nextint
        Nextint = Nextint+1
     else 
        PermutSommet(is) = Nsint + Nextext
        Nextext = Nextext + 1
     end IF
  Enddo
  !write(*,*)'PermutSommet=',PermutSommet

  !!----------------------------------------------------------------------
  !! permuter les coordonnees des sommets et type des sommets COORDS NTYPS
  !!----------------------------------------------------------------------

  Allocate(CoordSbis(2, Nbs),ntypSbis(Nbs))
  CoordSbis = CoordS
  NtypSbis  = NtypS
  DO is = 1, Nbs
     CoordS(1,PermutSommet(is))=  CoordSbis(1,is)
     CoordS(2,PermutSommet(is))=  CoordSbis(2,is)
     NtypS(PermutSommet(is)) = NtypSbis(is)
  EndDO
  !write(*,*)'CoordS avant permutation=',CoordSbis
  !write(*,*)'CoordS apres permutation=',CoordS

  !write(*,*)'ntyps avant permutation =', NtypSbis
  !write(*,*)'ntyps apres permutation =', NtypS
  DEALLOCATE(CoordSbis,Ntypsbis)

  !! --------------------------------------------------------
  !! Mise à jour des numeros des sommets des triangles: NuSoK
  !! --------------------------------------------------------
  Allocate(NuSoKbis(typegeomaille,Nbt))
  NuSoKbis = NuSoK
  DO kt = 1, Nbt
     Do m = 1, typegeomaille
        NuSoK (m,kt) = permutsommet(NuSoKbis(m,kt))
     END Do
  END DO
  !write(*,*)'NuSoK avant permutation=', NuSoKbis
  !write(*,*)'NuSoK apres permutation=', NuSoK
  Deallocate(NuSoKbis)

  !! --------------------------------------------------------
  !! Mise à jour des numeros des sommets des segments : NuSeg
  !! --------------------------------------------------------
  Allocate(NuSegbis(2,Nseg))
  NuSegbis = NuSeg
  Do iseg = 1,Nseg
     !! is = NuSegbis(1,iseg); js = NuSegbis(2,iseg)
     NuSeg (1,iseg) = permutsommet(NuSegbis(1,iseg))
     NuSeg (2,iseg) = permutsommet(NuSegbis(2,iseg))
  End do

  !write(*,*)'NuSeg avant permutation en Sommets =',NuSegbis
  !write(*,*)'NuSeg apres permutation en Sommets =',NuSeg
  Deallocate(NuSegbis)

  !! Les autres tableaux : 
  !! NuVoisK, Ntypseg, NumTVoisSeg, NombVoisSeg, CoordK : ces tableaux sont inchanges
  !!=================================================================================
  !!                 NUMEROTATION DES SEGMENTS A L'INTERIEUR
  !!=================================================================================
  !! A ce stade nous avons renumerote les inconnues à l'interieur 
  !! en gardant la nouvelle numerotation
  !! 
  !! Maintenant on va renumeroter les segments 
  !! a l'interieur et mettre à jour les tableau
  !!
  !! -----------------------------------------------------
  !! renumerotation des segmets interieurs de 1 a  Nsegint
  !! -----------------------------------------------------
  Allocate (PermutSegment(1:Nseg))
  Nextint = 1 ; Nextext = 1
  DO iseg = 1, Nseg
     IF( NtypSeg(iseg) == 0) then
        PermutSegment(iseg) = Nextint
        Nextint = Nextint + 1
     else
        PermutSegment(iseg) = Nextext + Nsegint
        Nextext = Nextext + 1
     end IF
  Enddo
  !write(*,*)'PermutSegment=',PermutSegment

  !!------------------------------------------------------
  !! Permuter les types des segment et numeros des sommets 
  !!------------------------------------------------------

  Allocate(NtypSegbis(Nseg), NuSegbis(2,Nseg))
  NtypSegbis  = NtypSeg 
  NuSegbis = NuSeg
  DO iseg = 1, Nseg
     Nuseg(1, PermutSegment(iseg)) =  Nusegbis(1,iseg)
     Nuseg(2, PermutSegment(iseg)) =  Nusegbis(2,iseg)
     NtypSeg(PermutSegment(iseg))  =  NtypSegbis(iseg)
  END DO
  !write(*,*)'NuSegbis avant permutation en Segments=', Nusegbis
  !write(*,*)'NuSeg apres permutation en segments=', NuSeg

  !write(*,*)'NtypSegbis avant  permutation =', NtypSegbis
  !write(*,*)'NtypSeg apres permutation =', NtypSeg

  DEALLOCATE(NtypSegbis)

  !--------------------------
  ! Mettre a jour nombVoisSeg
  !--------------------------
  DO iseg = 1, Nseg
     If (NTypSeg(iseg) == 0)  Then 
        NombVoisSeg(iseg) = 2
     Else 
        NombVoisSeg(iseg) = 1
     endif
  ENDDO
  ! ---------------------------------
  ! Mettre a jour NumTVoisSeg(2,Nseg)
  ! ---------------------------------

  Allocate(NumTVoisSegbis(2,Nseg))
  NumTVoisSegbis  = NumTVoisSeg
  DO iseg = 1, Nseg
     NumTVoisSeg(1:2,PermutSegment(iseg)) = NumTVoisSegbis(1:2,iseg)
  EndDO
  !write(*,*)'NumTVoisSegbis avant permutation =', NumTVoisSegbis
  !write(*,*)'NumTVoisSeg apres permutation =', NumTVoisSeg
  DEALLOCATE(NumTVoisSegbis)

  !----------------------------------------------------------------------------------
  !! Calcul des coordonnees aux milieux des segments
  !----------------------------------------------------------------------------------
  ALLOCATE(CoordMseg(2,Nseg))
  DO iseg = 1,Nseg
     is = NuSeg(1,iseg) ; js = NuSeg(2,iseg)
     CoordMseg(1,iseg)=(coordS(1,is)+coordS(1,js))/2.
     CoordMseg(2,iseg)=(coordS(2,is)+coordS(2,js))/2.
  End do

  ! ----------------------------------------------------------------------------------
  !                    Construction de Msig(Nseg) et Msige(Nseg)
  ! ----------------------------------------------------------------------------------
  ! Pour chaque segment iseg, Msig(Nseg) designe la distance entre les deux sommets du
  ! iseg, alors que Msige(Nseg) designe la distance entre les deux triangles voisins
  ! ----------------------------------------------------------------------------------
  ALLOCATE(Msig(Nseg),Msige(Nseg))
  DO iseg = 1, Nseg
     is = NuSeg(1,iseg); js = NuSeg(2,iseg)
     Msig(iseg) = sqrt((CoordS(1,is)-coordS(1,js))**2 +(CoordS(2,is)-coordS(2,js))**2)
     IF (NTypSeg(iseg) == 0)  THEN
        ik = NumTVoisSeg(1, iseg) ; jL = NumTVoisSeg(2, iseg)
        x(1) = coordK(1,ik) ; y(1) = coordK(2,ik)
        x(2) = coordK(1,jL) ; y(2) = coordK(2,jL)
        Msige(iseg) = sqrt((x(2)-x(1))**2 +(y(2)-y(1))**2)
     ELSE
        ik = NumTVoisSeg(1,iseg)
        x(1) = coordK(1,ik)    ; y(1) = coordK(2,ik)
        x(2) = coordMSeg(1,iseg) ; y(2) = coordMSeg(2,iseg)
        Msige(iseg) = sqrt((x(2)-x(1))**2 +(y(2)-y(1))**2)
     END IF
  END DO
  print*,'Le pas d''espace h = ',maxval(Msig)
  CALL prvarr(uprint, 'Max Msig  = ', maxval(Msig))
  CALL prvarr(uprint, 'Max Msige = ', maxval(Msige))
  ! -------------------------------------------------------------------------------------
  !                      Construction de NsigK(Nseg) et NsigeKe(Nseg)                   
  ! -------------------------------------------------------------------------------------
  ! Pour chaque segment iseg, NsigK(Nseg) designe le normal à l'interface primale sortant
  !  de K alors que NsigeKe(Nseg) designe le normal à l'interface duale sortant de Ke
  ! -------------------------------------------------------------------------------------
  ALLOCATE(NsigK(2,Nseg),NsigeKe(2,Nseg))
  !
  DO iseg = 1, Nseg
     is = NuSeg(1,iseg); js = NuSeg(2,iseg)  !! is,js numéro globale
     !! 1,2 numéro locale
     iK = NumTVoisSeg(1,iseg); jL = NumTVoisSeg(2,iseg) !! ik,jL numéro globale
     !! 1,2 numéro locale
     !!
     !
     IF (NTypSeg(iseg) == 0)  THEN
        !
        x(1) = coordS(1,is) ; y(1) = coordS(2,is)
        x(2) = coordS(1,js) ; y(2) = coordS(2,js)
        x(3) = coordK(1,iK) ; y(3) = coordK(2,iK)
        x(4) = coordK(1,jL) ; y(4) = coordK(2,jL)
        !---------------------------------------------------
        ! calcul du vecteur normal NsigK sortant de K vers L
        !---------------------------------------------------
        NsigK(1,iseg) = y(2)-y(1) ; NsigK(2,iseg) = -(x(2)-x(1))
        z = (x(4)-x(3))*NsigK(1,iseg) + (y(4)-y(3))*NsigK(2,iseg)
        IF ( z < 0.D0 ) THEN
           print*,'fausse orientation NsigK = ',iseg
           NsigK(1,iseg) = - NsigK(1,iseg)
           NsigK(2,iseg) = - NsigK(2,iseg)
        END IF
        !-------------------------------------------------------
        ! calcul du vecteur normal NsigeKe sortant de K* vers L*
        !-------------------------------------------------------
        NsigeKe(1,iseg) = -(y(4)-y(3)) ; NsigeKe(2,iseg) = x(4)-x(3)
        ! ici on effectue le prouduit scalaire pour vérifier la bonne
        ! direction de Nsigma*K*
        z = (x(2)-x(1))*NsigeKe(1,iseg) + (y(2)-y(1))*NsigeKe(2,iseg)
        !
        IF ( z <= 0.D0 ) THEN
           print*,' fausse orientation NsigeKe = ',iseg 
           NsigeKe(1,iseg) = - NsigeKe(1,iseg) 
           NsigeKe(2,iseg) = - NsigeKe(2,iseg)
        END IF
        !
     ELSE
        !
        x(1) = coordS(1,is)      ; y(1) = coordS(2,is)
        x(2) = coordS(1,js)      ; y(2) = coordS(2,js)
        x(3) = coordK(1,iK)      ; y(3) = coordK(2,iK)
        x(4) = CoordMseg(1,iseg) ; y(4) = CoordMseg(2,iseg)
        !-------------------------------------------------------
        ! calcul du vecteur normal NsigK au bord et sortant de K
        !-------------------------------------------------------
        NsigK(1,iseg) = y(2)-y(1) ; NsigK(2,iseg) = -(x(2)-x(1))
        z = (x(4)-x(3))*NsigK(1,iseg) + (y(4)-y(3))*NsigK(2,iseg)
        IF ( z <= 0.D0 ) THEN
           print*,'fausse orientation NsigK au bord = ',iseg
           NsigK(1,iseg) = - NsigK(1,iseg)
           NsigK(2,iseg) = - NsigK(2,iseg)
        END IF
        !-------------------------------------------------------
        ! calcul du vecteur normal NsigeKe sortant de K* vers L*
        !-------------------------------------------------------
        NsigeKe(1,iseg) = -(y(4)-y(3)) ; NsigeKe(2,iseg) = x(4)-x(3)
        ! ici on effectue le prouduit scalaire pour vérifier
        ! la direction de Nsigma*K*
        z = (x(2)-x(1))*NsigeKe(1,iseg) + (y(2)-y(1))*NsigeKe(2,iseg)
        !
        IF ( z <= 0.D0 ) THEN
           print*,' fausse orientation NsigeKe au bord = ',iseg 
           NsigeKe(1,iseg) = - NsigeKe(1,iseg) 
           NsigeKe(2,iseg) = - NsigeKe(2,iseg)
        END IF
     END IF
  END DO
  !! -------------------------------------------------------
  !!  calcul AireK des sommets pour le maillage triangulaire
  !!--------------------------------------------------------
  Allocate (AireK(Nbt))
  AireK = 0.D0
  Select case(Typegeomaille)
  case (triangle)
     Do jt = 1, Nbt
        i = NuSoK(1,jt) ; j = NuSoK(2,jt) ; k = NuSoK(3,jt) !! i,j,k numero globale

        x(1) = coordS(1,i) ; y(1) = coordS(2,i)
        x(2) = coordS(1,j) ; y(2) = coordS(2,j)
        x(3) = coordS(1,k) ; y(3) = coordS(2,k)
        !!
        AireK(jt) = ABS( (x(2)-x(1))* (y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1)) ) /2.
     End Do
  case (quadrangle) ! on calcule l'aire pour n'importe quelle ordre des sommets du quadrangle
     DO iseg = 1,Nseg
        Kvois = NumTVoisSeg(1,iseg)
        is = Nuseg(1,iseg) ; js = Nuseg(2,iseg)

        x(1) = coordK(1,Kvois) ; y(1) = coordK(2,Kvois)
        x(2) = coordS(1,is) ; y(2) = coordS(2,is)
        x(3) = coordS(1,js) ; y(3) = coordS(2,js)
        !!
        AireK(Kvois) = AireK(Kvois) + ABS( (x(2)-x(1))* (y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1)) ) /2.D0
        IF (NombvoisSeg(iseg)==2) then
           !!
           Lvois = NumTVoisSeg(2,iseg)
           is = NuSeg(1,iseg) ; js = NuSeg(2,iseg)
           !!
           x(1) = coordK(1,Lvois) ; y(1) = coordK(2,Lvois)
           x(2) = coordS(1,is) ; y(2) = coordS(2,is)
           x(3) = coordS(1,js) ; y(3) = coordS(2,js)
           !!
           AireK(Lvois) = AireK(Lvois) + &
                & ABS( (x(2)-x(1))* (y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1)) )/2.D0
        END IF
     END DO
  END Select
  CALL prvarr(uprint, 'Max aireK  = ', maxval(AireK))
  CALL prvarr(uprint, 'Min aireK  = ', minval(AireK))
  CALL prvarr(uprint, 'sum(AireK) = ', sum(AireK))
  !------------------------------
  ! Calcul des aires des diamants
  !------------------------------
  ALLOCATE(AireD(Nseg))
  Do iseg = 1, Nseg
     Kvois = NumTVoisSeg(1,iseg)
     is = Nuseg(1,iseg) ; js = Nuseg(2,iseg)

     x(1) = coordK(1,Kvois) ; y(1) = coordK(2,Kvois)
     x(2) = coordS(1,is) ; y(2) = coordS(2,is)
     x(3) = coordS(1,js) ; y(3) = coordS(2,js)
     !!
     AireD(iseg) = ABS( (x(2)-x(1))* (y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1)) ) /2.D0
     !!
     IF (NombvoisSeg(iseg)==2) then
        !!
        Lvois = NumTVoisSeg(2,iseg)
        is = NuSeg(1,iseg) ; js = NuSeg(2,iseg)
        !!
        x(1) = coordK(1,Lvois) ; y(1) = coordK(2,Lvois)
        x(2) = coordS(1,is) ; y(2) = coordS(2,is)
        x(3) = coordS(1,js) ; y(3) = coordS(2,js)
        !!
        AireD(iseg) = AireD(iseg)+ &
             & ABS( (x(2)-x(1))* (y(3)-y(1)) - (x(3)-x(1))*(y(2)-y(1)) )/2.D0
     END IF
  End Do
  CALL prvarr(uprint, 'Max aireD  = ', maxval(AireD))
  CALL prvarr(uprint, 'Min aireD  = ', minval(AireD))
  CALL prvarr(uprint, 'sum(AireD) = ', sum(AireD))
  ! ---------------------------------------------------------------
  ! Calcul AireDSommet associe a chaque sommet - maillage de Donald
  ! ---------------------------------------------------------------
  Allocate (AireDSommet(NbS))
  AireDSommet=0.D0
  Do jt = 1, Nbt
     i = NuSoK(1,jt) ; j = NuSoK(2,jt) ; k = NuSoK(3,jt) !! i,j,k numero globale

     AireDSommet(i) = AireDSommet(i) + AireK(jt)/Typegeomaille
     AireDSommet(j) = AireDSommet(j) + AireK(jt)/Typegeomaille
     AireDSommet(k) = AireDSommet(k) + AireK(jt)/Typegeomaille
     IF (Typegeomaille == quadrangle) THEN
        l = NuSoK(4,jt)
        AireDSommet(l) = AireDSommet(l) + AireK(jt)/Typegeomaille
     END IF
  END Do
  CALL prvarr(uprint, 'Max AireDSommet  = ', maxval(AireDSommet))
  CALL prvarr(uprint, 'Min AireDSommet  = ', minval(AireDSommet))
  CALL prvarr(uprint, 'sum(AireDSommet) = ', sum(AireDSommet))
  !====================================
  ! un test pour voir si ca bien marche
  !====================================
  Do iseg = 1, Nseg
     If(NTypSeg(iseg)==0 .and. NombVoisSeg(iseg)/= 2) then
        print*,'PROBLEME segment interieur iseg = ', iseg
        print*,'Nuseg', Nuseg(:,iseg)
        print*,'Ntypseg', Ntypseg(iseg)
        print*,'NumTVoisSeg',NumTVoisSeg(:, iseg)
        stop
     endif
     If(NTypSeg(iseg)/=0 .and. NombVoisSeg(iseg)/= 1) then
        print*,'PROBLEME segment bord = ',iseg
        stop
     endif
  END Do

  !------------------------
  ! Impressions eventuelles
  !------------------------
  IF (iprint >=5) THEN 
     CALL prvari(uprint, 'Nseg = ', Nseg)
     Write(uprint,*)'iseg, NombVoisSeg(iseg), NumTVoisSeg(1,iseg), NumTVoisSeg(2,iseg), Ntyp'
     DO iseg = 1 , Nseg 
        WRITE(uprint,200) iseg, NombVoisSeg(iseg), NumTVoisSeg(1,iseg), NumTVoisSeg(2,iseg),NTypSeg(iseg)
     ENDDO
     write(uprint,*)' is ,coordS(1,is),coordS(2,is),ntyps(is)'
     DO is =1,nbs
        WRITE(uprint,112) is ,coordS(1,is),coordS(2,is) ,ntyps(is)
     ENDDO
     write(uprint,*)' '
     write(uprint,*)'  kt , (NuSoK(j,kt),j=1,3 ou 4) '
     DO kt=1,nbt
        WRITE(uprint,114) kt , ( NuSoK(j,kt),j=1,4)
     END DO
     write(uprint,*)' '
     write(uprint,*)'  Nseg , (NuSeg(j,iseg),j=1,2), (NumTVoisSeg(j,iseg),j=1,2) '
     DO iseg=1,Nseg
        WRITE(uprint,114) iseg , (NuSeg(j,iseg),j=1,2), (NumTVoisSeg(j,iseg),j=1,2)
     END DO
     write(uprint,*)' '
     write(uprint,*)'  Nseg , (NsigK(j,iseg),j=1,2), (NsigeKe(j,iseg),j=1,2) '
     DO iseg=1,Nseg
        WRITE(uprint,115) iseg , (NsigK(j,iseg),j=1,2), (NsigeKe(j,iseg),j=1,2)
     END DO
  ENDIF

112 FORMAT (i8,2(E12.4,1x),i5)
114 FORMAT (5i8)
115 FORMAT (i8,4(E12.4,1x))
200 FORMAT (5(i8,2x))
300 FORMAT (2(E12.4,2x), I5)
  !-----------------
  ! Fin du programme
  !----------------
  prefix = oldprf

END SUBROUTINE meshtools



