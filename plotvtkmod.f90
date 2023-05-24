module plotvtkmod

  Use longr
  Use imprime
  Use parmmage

  Implicit none

contains

  Subroutine plot_vtk (vec,chaine,nomchamps)

    Real(kind=long), dimension(:),intent(in)   :: vec
    Character(len=*),intent(in)                :: chaine,nomchamps

    Integer                    :: kt, Lt, is, js, ks, iseg, countk, i, im, jm, km
    Real(kind=long), dimension(:), allocatable :: WW
    Integer, dimension(:,:), allocatable       :: NuSoDiam

    !print*,"Creation du fichier d'impression"
    uplotvtk = 16
    FPLOTVTK = chaine//'.vtk'

    open (unit=uplotvtk,file=FPLOTVTK,status='replace')
    write(uplotvtk,'(A)') '# vtk DataFile Version 3.0'
    write(uplotvtk,'(A)') 'LAPLACIEN  2D'
    write(uplotvtk,'(A)') 'ASCII'
    write(uplotvtk,'(A)') 'DATASET UNSTRUCTURED_GRID'

    !---------------------------------!
    ! maillage formé par les Diamants ! 
    !---------------------------------!

    write(uplotvtk,*) 'POINTS',Nbs+Nbt+Nsegbord,' float'

    ! 1. Construction de sous-mailles de Donald:
    !-------------------------------------------

    ! 1.1. Definir les coordonnes de sommets des sous-mailles: 
    DO is = 1, Nbs
       write (uplotvtk,*) CoordS(1,is), CoordS(2,is), 0.D0
    END DO

    DO kt = 1, Nbt
       write (uplotvtk,*) CoordK(1,kt), CoordK(2,kt), 0.D0
    END DO
    countk=0
    DO iseg = 1, Nseg
       Lt = NumTVoisSeg(2,iseg)
       !!write(*,*)'NumTVoisSeg(2,iseg)', NumTVoisSeg(2,iseg)
       if (Lt == 0) then
          !!write (uplotvtk,*) CoordMseg(1,iseg), CoordMseg(2,iseg), 0.D0
          is = NuSeg(1,iseg)
          write (uplotvtk,*) CoordS(1,is), CoordS(2,is), 0.D0
          countk = countk+1
       end if
    END DO
    !write(*,*)'nb seg bord via plot est ', countk
    ! 1.2. Definir le numero de sommets des diamants: NuSoDiam (1:4, Nseg+Nbt) 
    !--------------------------------------------------------------------------

    Allocate (NuSoDiam (1:4, Nseg) )
    countk = 1 
    do iseg = 1, Nseg
       is = NuSeg(1,iseg); js=NuSeg(2,iseg)
       kt = NumTVoisSeg(1,iseg); Lt= NumTVoisSeg(2,iseg)
       if (Lt > 0) then
          NuSoDiam (1,iseg)= is
          NuSoDiam (2,iseg)= kt + Nbs
          NuSoDiam (3,iseg)= js
          NuSoDiam (4,iseg)= Lt + Nbs            
       else
          NuSoDiam (1,iseg)= is
          NuSoDiam (2,iseg)= kt + Nbs
          NuSoDiam (3,iseg)= js
          NuSoDiam (4,iseg)= Nbt + Nbs + countk
          countk = countk+1 
       end if
    end do

    ! 1.3. Joindre les sommets de chaque sous-maille:
    !-------------------------------------------------
    write(uplotvtk,*) 'CELLS ',Nseg, 5*Nseg
    DO iseg = 1,Nseg          
       write (uplotvtk,*) 4, NuSoDiam(1,iseg)-1, NuSoDiam(2,iseg)-1, &
            & NuSoDiam(3,iseg)-1, NuSoDiam(4,iseg)-1
    END DO
    Deallocate (NuSoDiam)

    write(uplotvtk,*) 'CELL_TYPES ', Nseg

    DO iseg = 1, Nseg
       write(uplotvtk,*) 9
    END DO

    WRITE(uplotvtk,*) 'CELL_DATA',Nseg 
    WRITE(uplotvtk,*) 'SCALARS ',nomchamps,' float'
    WRITE(uplotvtk,*) 'LOOKUP_TABLE default'

    ! 2. Attribuer les valeurs aux diamands:
    !---------------------------------------
    Allocate (WW(Nseg))
    DO iseg = 1, Nseg
       is = NuSeg(1,iseg); js = NuSeg(2,iseg)
       Kt = NumTVoisSeg(1,iseg); Lt = NumTVoisSeg(2,iseg)
       Select case (NtypSeg(iseg))
       case (0) ! Pour les diamants ayant les sommets à l'intérieur
          IF (is <= NsInt) THEN
             IF (js <= NsInt) THEN
                WW(iseg) = (vec(is) + vec(js) + vec(Kt + NsInt) + vec(Lt + NsInt))/4.D0
             ELSE
                WW(iseg) = (vec(is) + vec(Kt + NsInt) + vec(Lt + NsInt))/3.D0
             END IF
          ELSE
             IF (js <= NsInt) THEN
                WW(iseg) = (vec(js) + vec(Kt + NsInt) + vec(Lt + NsInt))/3.D0
             ELSE
                WW(iseg) = (vec(Kt + NsInt) + vec(Lt + NsInt))/2.D0
             END IF
          END IF
       case (Dirichlet) ! Pour les diamants ayant des sommets sur le bord
          WW(iseg) = vec(Kt + NsInt)
       case (Neumann) ! Pour les diamants ayant des sommets sur le bord
          WW(iseg) = (3*(vec(is) + vec(js))/2 + vec(Kt + Nbs) )/4.D0
       END Select
    END DO
    ! 3. Creation des fichiers pour visit:
    !-------------------------------------

    DO iseg = 1, Nseg
       write (uplotvtk,500)  max(WW(iseg), 1.D-30)
    END DO
    deallocate(WW)
    close (uplotvtk)
    !print*,"OK"
    !print*," "

      500 format (E30.20)

  end subroutine plot_vtk


  Subroutine plot_vtk_primal (vec,chaine,nomchamps)

   Real(kind=long), dimension(:),intent(in)   :: vec
   Character(len=*),intent(in)                :: chaine,nomchamps

   Integer                    :: kt, Lt, is, js, ks, iseg, countk, i, im, jm, km
   Real(kind=long), dimension(:), allocatable :: WW
   Integer, dimension(:,:), allocatable       :: NuSoDiam

   !print*,"Creation du fichier d'impression"
   uplotvtk = 30
   FPLOTVTK = chaine//'.vtk'

   open (unit=uplotvtk,file=FPLOTVTK,status='replace')
   write(uplotvtk,'(A)') '# vtk DataFile Version 3.0'
   write(uplotvtk,'(A)') 'LAPLACIEN  2D'
   write(uplotvtk,'(A)') 'ASCII'
   write(uplotvtk,'(A)') 'DATASET UNSTRUCTURED_GRID'

   !---------------------------------!
   ! maillage primal (triangle) ! 
   !---------------------------------!

   write(uplotvtk,*) 'POINTS',Nbs,' float'

   ! 1. Construction des sous-mailles :
   !-------------------------------------------

   ! 1.1. Definir les coordonnes de sommets des sous-mailles: 

   DO is = 1, Nbs
      write (uplotvtk,*) CoordS(1,is), CoordS(2,is), 0.D0
   END DO

   ! 1.3. Joindre les sommets de chaque sous-maille:
   !-------------------------------------------------
   write(uplotvtk,*) 'CELLS ',Nbt, 4*Nbt
   DO kt = 1, Nbt        
      write (uplotvtk,*) 3, NuSoK (1,kt)-1, NuSoK (2,kt)-1 &
      & , NuSoK (3,kt)-1
   END DO

   write(uplotvtk,*) 'CELL_TYPES ', Nbt

   DO kt = 1, Nbt
      write(uplotvtk,*) 5  !7 = polygone
   END DO

   WRITE(uplotvtk,*) 'CELL_DATA',Nbt 
   WRITE(uplotvtk,*) 'SCALARS ',nomchamps,' float'
   WRITE(uplotvtk,*) 'LOOKUP_TABLE default'

   ! 3. Creation des fichiers pour visit:
   !-------------------------------------

   DO kt = 1, Nbt
      write (uplotvtk,500)  vec(kt)
   END DO

   close (uplotvtk)
   print*,"OK"
   print*," "

     500 format (E30.20)

 end subroutine plot_vtk_primal


 !Subroutine plot_vtk_dual (vec,chaine,nomchamps)
!
!    Real(kind=long), dimension(:),intent(in)   :: vec
!    Character(len=*),intent(in)                :: chaine,nomchamps
!
!    Integer                    :: kt, Lt, is, js, ks, iseg, countk, i, im, jm, km
!    Real(kind=long), dimension(:), allocatable :: WW
!    Integer, dimension(:,:), allocatable       :: NuSoDual
!    Integer, dimension(:), allocatable         :: NbSoSeg
!
!    print*,"Creation du fichier d'impression"
!    uplotvtk = 16
!    FPLOTVTK = chaine//'.vtk'
!
!    open (unit=uplotvtk,file=FPLOTVTK,status='replace')
!    write(uplotvtk,'(A)') '# vtk DataFile Version 3.0'
!    write(uplotvtk,'(A)') 'LAPLACIEN  2D'
!    write(uplotvtk,'(A)') 'ASCII'
!    write(uplotvtk,'(A)') 'DATASET UNSTRUCTURED_GRID'
!
!    !---------------------------------!
!    ! maillage dual ! 
!    !---------------------------------!
!
!    write(uplotvtk,*) 'POINTS',Nbt+Nsegbord,' float'
!
!    ! 1. Construction de sous-mailles de Donald:
!    !-------------------------------------------
!
!    ! 1.1. Definir les coordonnes de sommets des sous-mailles: 
!    
!    DO kt = 1, Nbt
!       write (uplotvtk,*) CoordK(1,kt), CoordK(2,kt), 0.D0
!    END DO
!
!    DO is = 1, Nbs
!      write (uplotvtk,*) CoordS(1,is), CoordS(2,is), 0.D0
!    END DO
!
!    countk=0
!    DO iseg = 1, Nseg
!       Lt = NumTVoisSeg(2,iseg)
!       !!write(*,*)'NumTVoisSeg(2,iseg)', NumTVoisSeg(2,iseg)
!       if (Lt == 0) then
!          write (uplotvtk,*) CoordMseg(1,iseg), CoordMseg(2,iseg), 0.D0
!          !is = NuSeg(1,iseg)
!          !write (uplotvtk,*) CoordS(1,is), CoordS(2,is), 0.D0
!          countk = countk+1
!       end if
!    END DO
!    write(*,*)'nb seg bord via plot est ', countk
!    ! 1.2. Definir le numero de sommets des diamants: NuSoDiam (1:4, Nseg+Nbt) 
!    !--------------------------------------------------------------------------
!
!    allocate(NbSoSeg(Nbs))
!    countk = 1 
!    NbSoSeg = 0
!    do iseg = 1, Nseg
!       is = NuSeg(1,iseg); js=NuSeg(2,iseg)
!       kt = NumTVoisSeg(1,iseg); Lt= NumTVoisSeg(2,iseg)
!       NbSoSeg(is) = NbSoSeg(is) + 1
!       NbSoSeg(js) = NbSoSeg(js) + 1
!    end do
!
!    Allocate (NuSoDual (maxval(NbSoSeg), Nbs) )
!    NuSoDual = 0
!    NbSoSeg = 0
!    countk = 1 
!    do iseg = 1, Nseg
!       is = NuSeg(1,iseg); js=NuSeg(2,iseg)
!       kt = NumTVoisSeg(1,iseg); Lt= NumTVoisSeg(2,iseg)
!       NbSoSeg(is) = NbSoSeg(is) + 1
!       NbSoSeg(js) = NbSoSeg(js) + 1
!       if (Lt > 0) then
!         NuSoDual (NbSoSeg(is),is)= kt + Nbs
!         NuSoDual (NbSoSeg(js),js)= Lt + Nbs
!         NuSoDual (3,iseg)= js
!         NuSoDual (4,iseg)= Lt + Nbs            
!       else
!         NuSoDual (1,iseg)= is
!         NuSoDual (2,iseg)= kt + Nbs
!         NuSoDual (3,iseg)= js
!         NuSoDual (4,iseg)= Nbt + Nbs + countk
!          countk = countk+1 
!       end if
!    end do
!
!    ! 1.3. Joindre les sommets de chaque sous-maille:
!    !-------------------------------------------------
!    write(uplotvtk,*) 'CELLS ',Nbs, 5*Nseg
!    DO iseg = 1,Nseg          
!       write (uplotvtk,*) 4, NuSoDiam(1,iseg)-1, NuSoDiam(2,iseg)-1, &
!            & NuSoDiam(3,iseg)-1, NuSoDiam(4,iseg)-1
!    END DO
!    Deallocate (NuSoDiam)
!    deallocate(NbSoSeg)
!
!    write(uplotvtk,*) 'CELL_TYPES ', Nbs
!
!    DO iseg = 1, Nbs
!       write(uplotvtk,*) 7
!    END DO
!
!    WRITE(uplotvtk,*) 'CELL_DATA',Nbs 
!    WRITE(uplotvtk,*) 'SCALARS ',nomchamps,' float'
!    WRITE(uplotvtk,*) 'LOOKUP_TABLE default'
!
!    
!    ! 2. Creation des fichiers pour visit:
!    !-------------------------------------
!
!    DO is = 1, Nbs
!       write (uplotvtk,500)  Vec(Nbt+is)
!    END DO
!    close (uplotvtk)
!    print*,"OK"
!    print*," "
!
!      500 format (E30.20)
!
!  end subroutine plot_vtk_dual


end module plotvtkmod
