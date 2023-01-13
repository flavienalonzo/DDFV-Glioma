SUBROUTINE   maillage
  !--------------------------------------------------------
  ! Nbs   : nombre des sommets total
  ! Nbt   : nombre des triangle
  ! NsInt : nombre des sommets interieur au domaine 
  !        nombre des sommets ou la solutions est calculee
  ! Nbord : nombre des sommets au bord
  ! NsInt = Nbs - Nbord
  !---------------------------------------------------------
  ! Declaration des variables locales
  !----------------------------------
  USE longr
  USE imprime
  USE parmmage
  IMPLICIT NONE 
  CHARACTER(len=6)                  :: oldprf
  CHARACTER(len=10)                 :: chaine
  INTEGER                           :: j,iseg,kt,is,js,iK,jL
  INTEGER                           :: bidon,bidon1,bidon2,bidon3,Nbt2,typet, nums,is1,isegmax,isegmin
  INTEGER                           :: ismax,ismin,jt,k,n,next,r1,r2,r3,is2,typeseg
  Integer                           :: Dimespace,Typedom,Typsommet 
  REAL(kind=long), DIMENSION(3)     :: x,y
  !-------------------
  ! Debut du programme
  !-------------------
  oldprf = prefix
  prefix = 'MAILLAGE'
  Select case (Maillage_Format)
  case (1)
     typegeomaille = triangle
     ! -------------------------------------------------------------- 
     !-- lecture des numeros des sommets de chaque triangle nom.1.ele
     ! --------------------------------------------------------------
     OPEN (uele,     file = nom_ele   ,  status = 'old')
     READ (uele,*) nbt, bidon, typet
     CALL prvari(uprint,'nombre des triangles   (Nbt) ', Nbt )
     ALLOCATE ( NuSoK(3,Nbt) )
     !
     if (typet==0) then 
     DO kt = 1,Nbt
        READ(uele,*) n, ( NuSoK(j,kt),j=1,3 )
     END DO
      else 
         allocate ( TypS(typet,Nbt) )
         DO kt = 1,Nbt
            READ(uele,*) n, ( NuSoK(j,kt),j=1,3 ), (TypS(k,kt),k=1,typet)
         END DO
      end if
     close(uele)
     ! -------------------------------------------------------- 
     ! -- lecture des coordonnees de chaque sommet : nom.1.node
     ! --------------------------------------------------------
     ! 
     OPEN (unode,     file = nom_node   ,  status = 'old')
     !!
     READ (unode,*) Nbs, Dimespace,Typedom,Typsommet  
     !!
     CALL prvari(uprint,'nombre des sommets     (Nbs) ', Nbs )
     !!
     ALLOCATE(coordS(2,Nbs))
     Select Case(Typedom)
     case(0)
        DO is = 1,nbs
           READ(unode,*) n,(coordS(j,is),j=1,2), r3
        END DO
     Case(1) 
        DO is = 1,nbs
           READ(unode,*) n,(coordS(j,is),j=1,2), r1, r3
        END DO
     Case(2) 
        DO is = 1,nbs
           READ(unode,*) n,(coordS(j,is),j=1,2), r1, r2, r3
        END DO
     end Select
     close(unode)
     ! ------------------------------------------------------------- 
     !-- lecture de la liste des segments et des types nom.1.edge
     ! -----------------------------------------------------------
     OPEN (uedge,     file = nom_edge   ,  status = 'old')
     !!
     READ (uedge,*) Nseg, typeseg
     !!
     CALL prvari(uprint,'nombre des segments   (Nseg) ', Nseg )
     !!
     ALLOCATE ( NuSeg(2,Nseg), NTypSeg(Nseg) )

     DO iseg = 1,Nseg
        READ(uedge,*) n, NuSeg(1,iseg), NuSeg(2,iseg), NTypSeg(iseg)
     END DO
     close(uedge)
     !===================================================================
     ! Construction  de la liste des deux triangles 
     !===================================================================
     ALLOCATE( NumTVoisSeg(2,NSeg))  ! Numero des triangles voisins
     ALLOCATE( NombVoisSeg(NSeg))    ! Nombre de voisins au plus 2; 
     !sinon 1 (c-a-d segment sur le bord)
     NombVoisSeg = 0 ; NumTVoisSeg = 0
     DO jt = 1, Nbt
        DO k = 1,typegeomaille
           next = MOD(k,typegeomaille)+1
           is1 = NuSoK(k,   jt)
           is2 = NuSoK(next,jt)
           ismin = MIN(is1,is2)
           ismax = MAX(is1,is2)
           !!
           DO iseg = 1, Nseg
              isegmin = MIN(NuSeg(1,iseg), NuSeg(2,iseg))
              isegmax = MAX(NuSeg(1,iseg), NuSeg(2,iseg))

              IF(ismin == isegmin.AND. ismax == isegmax) Then
                 NombVoisSeg(iseg) = NombVoisSeg(iseg)+1  !! Nombre de Voisins
                 NumTVoisSeg(NombVoisSeg(iseg), iseg)=jt            !! 
                 EXIT
              END IF
           ENDDO
        END DO
     END DO
     write(uprint,*)'  Nseg , (NuSeg(j,kt),j=1,2) '
     DO iseg = 1,Nseg
        WRITE(uprint,114) iseg , (NuSeg(j,iseg),j=1,2), (NumTVoisSeg(j,iseg),j=1,2)
     ENDDO
  case (2)
     ! fichier benchmark boyer
     ! -------------------------------------------
     ! -- lecture des coordonn�es de chaque sommet 
     ! -------------------------------------------
     !
     OPEN(unode,     file = nom_mesh  ,  status = 'old')
     READ(unode,*)
     READ(unode,*) Nbs
     !!
     CALL prvari(uprint,'nombre des sommets     (Nbs) ', Nbs )
     !!
     ALLOCATE(coordS(2,Nbs) )
     DO is = 1, Nbs
        READ(unode,*) coordS(1:2,is)
     END DO
     ! ---------------------------------------------------
     ! -- Lecture des sommets de chaque volume de controle 
     ! ---------------------------------------------------
     !
     READ(unode,*)
     READ(unode,*) Nbt
     !
     IF (Nbt == 0) THEN
        Typegeomaille = quadrangle
        READ(unode,*)
        READ(unode,*) Nbt
     ELSE
        Typegeomaille = triangle
     END IF
     !
     ALLOCATE(NuSoK(Typegeomaille,Nbt))
     !
     IF (Typegeomaille == triangle) THEN
        DO kt = 1,Nbt
           READ(unode,*) NuSoK(1:Typegeomaille,kt)
        END DO
        Do j = 1,3
           READ(unode,*)
        End Do
     ELSE
        DO kt = 1,Nbt
           READ(unode,*) NuSoK(1:Typegeomaille,kt)
        END DO
        READ(unode,*)
     END IF
     ! ----------------------------------- 
     ! -- Lecture de la liste des segments
     ! -----------------------------------
     READ(unode,*) Nbord
     Do is = 1,Nbord
        READ(unode,*)
     End Do
     READ(unode,*)
     READ(unode,*) Nseg
     ALLOCATE (NuSeg(2,Nseg),NumTVoisSeg(2,Nseg))
     Do iseg = 1, NSeg
        READ(unode,*) NuSeg(1:2,iseg), NumTVoisSeg(1:2,iseg)
     End Do
     !
     close(unode)
  END Select
  ! --------------------------------------------
  ! Ici on calcul les coordonnees barycentriques
  ! --------------------------------------------
  ALLOCATE(CoordK(2,Nbt))
  Select case (Typegeomaille)
  case(triangle)
     DO kt = 1,Nbt
        CoordK(1,kt) = (CoordS(1,NuSoK(1,kt)) + CoordS(1,NuSoK(2,kt)) + CoordS(1,NuSoK(3,kt)) )/Typegeomaille 
        CoordK(2,kt) = (CoordS(2,NuSoK(1,kt)) + CoordS(2,NuSoK(2,kt)) + CoordS(2,NuSoK(3,kt)) )/Typegeomaille
     END DO
  case(quadrangle)
     DO kt = 1,Nbt
        CoordK(1,kt) = (CoordS(1,NuSoK(1,kt)) + CoordS(1,NuSoK(2,kt)) + &
             & CoordS(1,NuSoK(3,kt)) + CoordS(1,NuSoK(4,kt)) )/Typegeomaille 
        CoordK(2,kt) = (CoordS(2,NuSoK(1,kt)) + CoordS(2,NuSoK(2,kt)) + &
             & CoordS(2,NuSoK(3,kt)) + CoordS(2,NuSoK(4,kt)) )/Typegeomaille
     END DO
  case default
     print*,'erreur dans la geometrie du maillage'
     stop
  END Select
 
  !------------------------
  ! Impressions eventuelles
  !------------------------
  IF (iprint >=2) THEN
     write(uprint,*)' '
     If (Typegeomaille == triangle) THEN
        write(uprint,*)' type du maillage = triangulation'
     ELSE
        write(uprint,*)' type du maillage = quadrangulation'
     END If
     write(uprint,*)' '
     write(uprint,*)' is ,coordS(1,is),coordS(2,is)'
     DO is =1,nbs
        WRITE(uprint,112) is ,coordS(1,is),coordS(2,is)
     ENDDO
     write(uprint,*)' '
     write(uprint,*)'  kt , (NuSoK(j,kt),j=1,3 ou 4) '
     DO kt=1,nbt
        WRITE(uprint,113) kt , ( NuSoK(j,kt),j=1,4)
     ENDDO
     write(uprint,*)' '
     write(uprint,*)'  Nseg , (NuSeg(j,kt),j=1,2) '
     DO kt=1,Nseg
        WRITE(uprint,114) kt , (NuSeg(j,kt),j=1,2), (NumTVoisSeg(j,kt),j=1,2)
     ENDDO
     write(uprint,*)' '
     write(uprint,*)'  kt , (CoordK(j,kt),j=1,2) '
     DO kt=1,nbt
        WRITE(uprint,112) kt , (CoordK(j,kt),j=1,2)
     ENDDO
  ENDIF
  
112 FORMAT (i6,2(E16.8,1x))
113 FORMAT (5i6)
114 FORMAT (i6,4i8)
  !-----------------
  ! Fin du programme
  !-----------------
  prefix = oldprf

END SUBROUTINE maillage



