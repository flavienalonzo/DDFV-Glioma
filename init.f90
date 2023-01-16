!            ***********************
!            **  SUBROUTINE INIT  **
!            ***********************
!***************************************************
!     * Ce sous programme lit le fichier uread
!     * le fichier d'entree
!***************************************************

SUBROUTINE  INIT
  !--------
  ! Modules
  !--------
  !
  USE longr
  USE imprime
  USE parmmage

  IMPLICIT NONE

  !--------------------------
  ! Declaration des arguments
  !--------------------------
  !----------------------------------
  ! Declaration des variables locales
  !----------------------------------
  CHARACTER(len=6)          :: oldprf
  CHARACTER(len=len_buffer) :: buffer,char_mesh
  INTEGER :: Ndivx, Ndivy, ii, jj, i, j , k, ieff, jeff, imx, jmy
  !-------------------
  ! Debut du programme
  !-------------------
  oldprf = prefix
  prefix = 'INIT'

  uread        =  9   ! fichier de donnee
  uprint       =  10  ! Unite de listage pour verifier les donnees
  usave        =  14  ! Unite de sauvegarde des solutions
  unode        =  15  ! unite de lecture du maillage
  uele         =  16  ! unite de lecture du maillage
  uedge        =  17  ! unite de lecture du maillage
  uplotvtk     =  18  ! ficher de sauvgarde
  upoints      =  19  ! ficher de sauvgarde
  NormU        =  20  ! ficher de sauvgarde
  NormV        =  21  ! ficher de sauvgarde
  cdtinitialeU =  22  ! fichier de sauvgarde a un instant donne
  cdtinitialeV =  23  ! fichier de sauvgarde a un instant donne

  print*,'UwU 48'
  OPEN (uread, file = 'UREAD'  ,  status = 'old') 
  OPEN (uprint, file = 'UPRINT',  status = 'unknown')
  OPEN (NormU, file = 'NormU'  ,  status = 'unknown')
  OPEN (NormV, file = 'NormV'  ,  status = 'unknown')
  OPEN (cdtinitialeU,file = 'cdtinitialeU',  status = 'unknown')
  OPEN (cdtinitialeV,file = 'cdtinitialeV',  status = 'unknown')

  buffer=lireligne(uread)
  READ(buffer, *, err = 10)  selec_mesh
  CALL prvarr (uprint, 'Selec_mesh                      : ', selec_mesh)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10)  Choixdt
  CALL prvari (uprint, 'Choix_dt                        : ', Choixdt)

  IF (Choixdt /= 0 .AND. Choixdt /= 1) THEN
     print*,'Erreur Choixdt'
     stop
  END IF

  IF (Choixdt == 1) THEN
     buffer=lireligne(uread)
     READ(buffer, *, err = 10)  dt
     CALL prvarr (uprint, 'dt                              : ', dt)
  END IF

  buffer=lireligne(uread)
  READ(buffer, *, err = 10)  ChoixBord
  CALL prvari (uprint, 'ChoixBord                       : ', ChoixBord)
  IF (ChoixBord /= 0 .AND. ChoixBord /= 1) THEN
     print*,'Erreur ChoixBord'
     print*, ChoixBord
     stop
  END IF
  IF (ChoixBord == 1) THEN
     buffer=lireligne(uread)
     READ(buffer, *, err = 10)  ChoixCdtBord
     CALL prvari (uprint, 'ChoixCdtBord                 : ', ChoixCdtBord)
  END IF
  buffer=lireligne(uread)
  READ(buffer, *, err = 10)  NbiterTemps
  CALL prvari (uprint, 'NbiterTemps                     : ', NbiterTemps)
  IF (ChoixBord == 0) THEN
     print*,' '
     print*,'========================================='
     print*,'******** CHOIX Condition au Bord ********'
     print*,'========================================='
     print*,'Dirichlet        = 1'
     print*,'Neumann homogene = 2'
     print*,'Choix Condition au Bord :'
     read *, ChoixCdtBord
     !
     do while (Int(ChoixCdtBord) /= ChoixCdtBord .OR. ChoixCdtBord < 1 .OR. ChoixCdtBord > 2)
        print*,'********* ERREUR choix Condition au Bord *********'
        print*,'La valeur doit etre soit 1 pour Dirichlet soit 2 pour Neumann'
        print*,' '
        print*,'Choix Condition au Bord :'
        read *,ChoixCdtBord
     end do
     !
  END IF
  buffer=lireligne(uread)
  READ(buffer, *, err = 10)  CFDT
  CALL prvarr (uprint, 'CFDT                            : ', cfdt)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10)  epsilon
  CALL prvarr (uprint, 'epsilon                         : ', epsilon)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10)  TolerenceGradient
  CALL prvarr (uprint, 'TolerenceGradient               : ', TolerenceGradient)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10)  TolerenceNewton
  CALL prvarr (uprint, 'TolerenceNewton                 : ', TolerenceNewton)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10)  Tolerenceiterative
  CALL prvarr (uprint, 'Tolerenceiterative              : ', Tolerenceiterative)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10)  kitermax
  CALL prvari (uprint, 'ItermaxNewton                   : ', kitermax)

  buffer=lireligne(uread)
  !print*,'dossier=',buffer
  ! READ(buffer, *, err = 10)  dossiermaillage
  dossiermaillage=trim(adjustl(buffer))
  !print*,'dossiermaillage=',dossiermaillage
  CALL prchar (uprint, 'Dossier Maillage                : '//dossiermaillage)
  print*, 'UwU 140'
  If (selec_mesh == 0) THEN
     buffer=lireligne(uread)
     READ(buffer, *, err = 10) Maillage_Format
     CALL prvari (uprint, 'Maillage_Format                  : ', Maillage_Format)

     Select case (Maillage_format)
     case (1)  ! maillage du logiciel triangle
        buffer=lireligne(uread)
        READ(buffer, *, err = 10 ) nom_ele
        nom_ele=trim(dossiermaillage)//trim(nom_ele)
        CALL prchar (uprint, 'Fichier des elements nom_ele   : '//nom_ele) 

        buffer=lireligne(uread)
        READ(buffer, *, err = 10 )  nom_node
        nom_node=trim(dossiermaillage)//trim(nom_node)
        CALL prchar (uprint, 'Fichier du centres             : '//nom_node) 

        buffer=lireligne(uread)
        READ(buffer, *, err = 10 )  nom_edge
        nom_edge=trim(dossiermaillage)//trim(nom_edge)
        CALL prchar (uprint, 'Fichier des segments           : '//nom_edge) 

        buffer=lireligne(uread)
     CASE (2)  ! maillage Boyer (Benchmark)
      buffer=lireligne(uread)
      buffer=lireligne(uread)
      buffer=lireligne(uread)
        buffer=lireligne(uread)
        READ(buffer, *, err = 10 )  nom_mesh
        nom_mesh=trim(dossiermaillage)//trim(nom_mesh)
        CALL prchar (uprint, 'Fichier du centres             : '//nom_mesh) 
     case default
        print*,'Erreur dans entree Maillage_Format'
        stop
     END Select
  ELSE IF (selec_mesh == 1) THEN
     print*,' '
     print*,'========================================='
     print*,'*********** CHOIX DU MAILLAGE ***********'
     print*,'========================================='
     print*,'Maillage triangulaire   = 1'
     print*,'Maillage quadrangulaire = 2'
     print*,'Choix du type du maillage :'
     read *, Choixmesh
     !
     do while (Int(Choixmesh) /= Choixmesh .OR. Choixmesh < 1 .OR. Choixmesh > 2)
        print*,'********* ERREUR dans le choix du maillage *********'
        print*,'La valeur doit etre soit 1 pour une triangulation soit 2 pour une quadrangulation'
        print*,' '
        print*,'Choix du type du maillage :'
        read *,Choixmesh
     end do
     !
     Maillage_Format = Choixmesh
     select case (INT(Choixmesh))
     case (1)
        print*,' '
        print*,'========================================='
        print*,'********* MAILLAGE TRIANGULAIRE *********'
        print*,'========================================='
        print*,'01. mesh_triangle_1 :  40 triangles'
        print*,'02. mesh_triangle_2 :  224 triangles'
        print*,'03. mesh_triangle_3 :  934 triangles'
        print*,'04. mesh_triangle_4 :  6422 triangles'
        print*,'05. mesh_triangle_5 :  25872 triangles'
        print*,'06. mesh_triangle_6 :  104420 triangles'
        print*,'******* Maillage orthogonale *******'
        print*,'07. mesh_VF_1.typ1  :  56 triangles'
        print*,'08. mesh_VF_2.typ1  :  224 triangles'
        print*,'09. mesh_VF_3.typ1  :  896 triangles'
        print*,'10. mesh_VF_4.typ1  :  3584 triangles'
        print*,'11. mesh_VF_5.typ1  :  14336 triangles'
        print*,'*** Mesh obtenu avec le logiciel triangle ***'
        print*,'12. mesh_Breast_1   :  898 triangles'
        print*,'13. mesh_Breast_2   :  1859 triangles'
        print*,'14. mesh_Breast_3   :  4522 triangles'
        print*,'15. mesh_Breast_4   :  9148 triangles'
        print*,'16. mesh_Carre_1    :  1547 triangles'
        print*,'17. mesh_Carre_2    :  3132 triangles'
        print*,'18. mesh_Carre_3    :  7845 triangles'
        print*,'Choix du type de triangulation :'
        read *, Choixtri
        !
        do while (Int(Choixtri) /= Choixtri .OR. Choixtri > 18 .OR. Choixtri < 1)
           print*,'********* ERREUR dans le choix de triangulation *********'
           print*,'La valeur doit etre un entier entre 1 et 18'
           print*,' '
           print*,'Choix du type de triangulation :'
           read *,Choixtri
        end do
        !
        If (Choixtri < 7) THEN
           Maillage_Format = 2
           write(char_mesh,'(I5)')INT(Choixtri)
           nom_mesh = 'mesh_tri_'//trim(adjustl(char_mesh))//'.typ1'
           nom_mesh = trim(dossiermaillage)//trim(nom_mesh)
           print*,'Fichier maillage est: ',nom_mesh
           CALL prchar (uprint, 'Fichier du centres              : '//nom_mesh)
        ELSE IF (Choixtri < 12) THEN
           Maillage_Format = 2
           write(char_mesh,'(I5)')INT(Choixtri)-6
           nom_mesh = 'mesh_VF_'//trim(adjustl(char_mesh))//'.typ1'
           nom_mesh = trim(dossiermaillage)//trim(nom_mesh)
           print*,'Fichier maillage est: ',nom_mesh
           CALL prchar (uprint, 'Fichier du centres              : '//nom_mesh)
        ELSE IF (Choixtri < 16) THEN
           write(char_mesh,'(I5)')INT(Choixtri)-11
           nom_ele = 'Breast.'//trim(adjustl(char_mesh))//'.ele'
           nom_ele = trim(dossiermaillage)//trim(nom_ele)
           nom_node = 'Breast.'//trim(adjustl(char_mesh))//'.node'
           nom_node = trim(dossiermaillage)//trim(nom_node)
           nom_edge = 'Breast.'//trim(adjustl(char_mesh))//'.edge'
           nom_edge = trim(dossiermaillage)//trim(nom_edge)
           nom_mesh = 'mesh_Breast_'//trim(adjustl(char_mesh))
           nom_mesh = trim(dossiermaillage)//trim(nom_mesh)
           print*,'Fichier maillage est: ',nom_mesh
           CALL prchar (uprint, 'Fichier du centres              : '//nom_mesh)
        ELSE
           write(char_mesh,'(I5)')INT(Choixtri)-15
           nom_ele = 'Carre.'//trim(adjustl(char_mesh))//'.ele'
           nom_ele = trim(dossiermaillage)//trim(nom_ele)
           nom_node = 'Carre.'//trim(adjustl(char_mesh))//'.node'
           nom_node = trim(dossiermaillage)//trim(nom_node)
           nom_edge = 'Carre.'//trim(adjustl(char_mesh))//'.edge'
           nom_edge = trim(dossiermaillage)//trim(nom_edge)
           nom_mesh = 'mesh_Carre_'//trim(adjustl(char_mesh))
           nom_mesh = trim(dossiermaillage)//trim(nom_mesh)
           print*,'Fichier maillage est: ',nom_mesh
           CALL prchar (uprint, 'Fichier du centres              : '//nom_mesh)
        END If
     case (2)
        print*,' '
        print*,'==========================================='
        print*,'********* MAILLAGE QUADRANGULAIRE *********'
        print*,'==========================================='
        print*,'01. mesh_quad_1.typ1 :  16 quadrangles'
        print*,'02. mesh_quad_2.typ1 :  64 quadrangles'
        print*,'03. mesh_quad_3.typ1 :  256 quadrangles'
        print*,'04. mesh_quad_4.typ1 :  1024 quadrangles'
        print*,'05. mesh_quad_5.typ1 :  4096 quadrangles'
        print*,'06. mesh_quad_6.typ1 :  16384 quadrangles'
        print*,'07. mesh_quad_7.typ1 :  65536 quadrangles'
        print*,'******* Maillage de type Kershaw *******'
        print*,'08. mesh4_1_1.typ1  :  289 quadrangles'
        print*,'09. mesh4_1_2.typ1  :  1156 quadrangles'
        print*,'10. mesh4_1_3.typ1  :  2601 quadrangles'
        print*,'11. mesh4_1_4.typ1  :  4624 quadrangles'
        print*,'12. mesh4_2_1.typ1  :  1089 quadrangles'
        print*,'13. mesh4_2_2.typ1  :  4356 quadrangles'
        print*,'14. mesh4_2_3.typ1  :  9801 quadrangles'
        print*,'Choix du type de quadrangulation :'
        read *, Choixtri
        do while (Int(Choixtri) /= Choixtri .OR. Choixtri > 14 .OR. Choixtri < 1)
           print*,'********* ERREUR dans le choix de triangulation *********'
           print*,'La valeur doit etre un entier entre 1 et 6'
           print*,' '
           print*,'Choix du type de triangulation :'
           read *,Choixtri
        end do
        If (Choixtri < 8) THEN
           write(char_mesh,'(I5)')INT(Choixtri)
           nom_mesh = 'mesh_quad_'//trim(adjustl(char_mesh))//'.typ1'
           nom_mesh = trim(dossiermaillage)//trim(nom_mesh)
           print*,'Fichier maillage est: ',nom_mesh
           CALL prchar (uprint, 'Fichier du centres             : '//nom_mesh)
        ELSEIF (Choixtri < 12) THEN
           write(char_mesh,'(I5)')INT(Choixtri)-7
           nom_mesh = 'mesh4_1_'//trim(adjustl(char_mesh))//'.typ1'
           nom_mesh = trim(dossiermaillage)//trim(nom_mesh)
           print*,'Fichier maillage est: ',nom_mesh
           CALL prchar (uprint, 'Fichier du centres             : '//nom_mesh)
        ELSE
           write(char_mesh,'(I5)')INT(Choixtri)-11
           nom_mesh = 'mesh4_2_'//trim(adjustl(char_mesh))//'.typ1'
           nom_mesh = trim(dossiermaillage)//trim(nom_mesh)
           print*,'Fichier maillage est: ',nom_mesh
           CALL prchar (uprint, 'Fichier du centres             : '//nom_mesh)
        END If
     END select
  ELSE
     PRINT*,'Erreur Selec_mesh'
     stop
  END If

  buffer=lireligne(uread)
  print*, buffer
  READ(buffer, *, err = 10)  CoefDiffU
  CALL prvarr (uprint, 'CoefDiffU                       : ', CoefDiffU)
  
  buffer=lireligne(uread)
  print*, buffer
  READ(buffer, *, err = 10)  CoefDiffC
  CALL prvarr (uprint, 'CoefDiffC                       : ', CoefDiffC)

  buffer=lireligne(uread)
  print*, buffer
  READ(buffer, *, err = 10)  CoefDiffE
  CALL prvarr (uprint, 'CoefDiffE                       : ', CoefDiffE)
  
  buffer=lireligne(uread)
  print*, buffer
  READ(buffer, *, err = 10)  CoefDiffV
  CALL prvarr (uprint, 'CoefDiffV                       : ', CoefDiffV)

  buffer=lireligne(uread)
  print*, buffer
  READ(buffer, *, err = 10)  CoefDiffWm
  CALL prvarr (uprint, 'CoefDiffWm                       : ', CoefDiffWm)

  buffer=lireligne(uread)
  print*, buffer
  READ(buffer, *, err = 10)  CoefDiffGm
  CALL prvarr (uprint, 'CoefDiffGm                       : ', CoefDiffGm)

  buffer=lireligne(uread)
  print*, buffer
  READ(buffer, *, err = 10)  CoefDiffPs
  CALL prvarr (uprint, 'CoefDiffPs                       : ', CoefDiffPs)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10)  rho1
  CALL prvarr (uprint, 'rho1                           : ', rho1)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10)  beta1
  CALL prvarr (uprint, 'beta1                            : ', beta1)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10)  alpha2
  CALL prvarr (uprint, 'alpha2                           : ', alpha2)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10)  beta2
  CALL prvarr (uprint, 'beta2                            : ', beta2)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10)  gamma2
  CALL prvarr (uprint, 'gamma2                            : ', gamma2)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10)  rho3
  CALL prvarr (uprint, 'rho3                           : ', rho3)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10)  beta3
  CALL prvarr (uprint, 'beta3                            : ', beta3)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10)  alpha4
  CALL prvarr (uprint, 'alpha4                           : ', alpha4)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10)  beta4
  CALL prvarr (uprint, 'beta4                            : ', beta4)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10)  gamma4
  CALL prvarr (uprint, 'gamma4                            : ', gamma4)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10)  Chypo
  CALL prvarr (uprint, 'Chypo                            : ', Chypo)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10)  Cnecro
  CALL prvarr (uprint, 'Cnecro                            : ', Cnecro)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10)  Use_surgery

  buffer=lireligne(uread)
  READ(buffer, *, err = 10)  Use_radio
  
  buffer=lireligne(uread)
  READ(buffer, *, err = 10)  Use_chemo

  print*, Use_surgery,Use_radio,Use_chemo
  buffer=lireligne(uread)
  READ(buffer, *, err = 10) Dose_radio

  buffer=lireligne(uread)
  READ(buffer, *, err = 10) Dose_chemo

  buffer=lireligne(uread)
  READ(buffer, *, err = 10) Gain_radio

  print*, Dose_radio,Dose_chemo,Gain_radio

  buffer=lireligne(uread)
  READ(buffer, *, err = 10) choixpb
  CALL prvari (uprint, 'choixpb                         : ', choixpb)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10) choixgb
  CALL prvari (uprint, 'choixgb                         : ', choixgb)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10) choixanisu
  CALL prvari (uprint, 'choixanisu                      : ', choixanisu)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10) choixanisc
  CALL prvari (uprint, 'choixanisc                      : ', choixanisc)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10) choixanise
  CALL prvari (uprint, 'choixanise                      : ', choixanise)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10) choixanisv
  CALL prvari (uprint, 'choixanisv                      : ', choixanisv)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10) ChoixAdegU
  CALL prvari (uprint, 'ChoixAdegU                       : ', ChoixAdegU)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10) CoefDiffuAdeg
  CALL prvarr (uprint, 'CoefDiffuAdeg                   : ', CoefDiffuAdeg)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10) ChoixMuU
  CALL prvari (uprint, 'ChoixMuU                         : ', ChoixMuU)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10) CoefTranspMuU
  CALL prvarr (uprint, 'CoefTranspMuU                    : ', CoefTranspMuU)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10) ChoixAdegE
  CALL prvari (uprint, 'ChoixAdegE                       : ', ChoixAdegE)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10) CoefDiffeAdeg
  CALL prvarr (uprint, 'CoefDiffeAdeg                   : ', CoefDiffeAdeg)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10) ChoixMuE
  CALL prvari (uprint, 'ChoixMuE                         : ', ChoixMuE)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10) CoefTranspMuE
  CALL prvarr (uprint, 'CoefTranspMuE                    : ', CoefTranspMuE)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10) ubar
  CALL prvarr (uprint, 'u_bar                           : ', ubar)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10)  gamma
  CALL prvari (uprint, 'gamma                           : ', gamma)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10)  iprint
  CALL prvari (uprint, 'IPRINT                          : ', iprint)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10) nom_usave
  CALL prchar (uprint, 'Fichier du maillage USAVE        : '//nom_usave) 

  buffer=lireligne(uread)
  READ(buffer, *, err = 10) nom_upoints
  CALL prchar (uprint, 'Fichier du maillage  UPOINTS     : '//nom_upoints) 
  OPEN (upoints, file = nom_upoints     ,  status = 'unknown')

  ! -----------------------------
  ! temps de sortie des solutions
  ! -----------------------------

  buffer=lireligne(uread)
  READ(buffer, *, err = 10) Nb_iter_visit
  call prvari (uprint, 'Nb_iter_visit                     : ',Nb_iter_visit)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10) ntarr
  call prvari (uprint, 'ntarr                             : ',ntarr)

  Allocate (Tempstock(ntarr))

  buffer=lireligne(uread)
  READ(buffer, *, err = 10) (Tempstock(i),i=1,ntarr)
  call prchar (uprint,'Temps de stockage des solutions')
  write (uprint,*)(Tempstock(i),i=1,ntarr)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10) TypecondinitialeU
  CALL prvari (uprint, 'TypecondinitialeU       : ', TypecondinitialeU)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10) TypecondinitialeC
  CALL prvari (uprint, 'TypecondinitialeC       : ', TypecondinitialeC)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10) TypecondinitialeE
  CALL prvari (uprint, 'TypecondinitialeE       : ', TypecondinitialeE)

  buffer=lireligne(uread)
  READ(buffer, *, err = 10) TypecondinitialeV
  CALL prvari (uprint, 'TypecondinitialeV       : ', TypecondinitialeV)

  !-----------------
  ! Fin du programme
  !-----------------
  prefix = oldprf
  RETURN

10 PRINT*,"Erreur dans l'entree des parametres INIT"

  STOP


END SUBROUTINE INIT








