      MODULE ocean_mod
      
!!======================================================================
!! Jully, 2009                                                         ! 
!!======================================================Ligia Pinto=====
!!                                                                     !
!! This module is adapted from ROMS (ocean_mod).                       !
!!                                                                     ! 
!!======================================================================   
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !

        !USE elfe_glbl, only: nvrt,nea,npa,mnei,ntracers
        use sed_param, only:r8
!        use sed_mod

        IMPLICIT NONE
        SAVE
        
        real(r8), allocatable :: Hz(:,:),Zob(:),bed(:,:,:),bed_frac(:,:,:), &
     &bed_mass(:,:,:,:),bottom(:,:),bedldu(:,:),bedldv(:,:)

#ifdef SED_MORPH
        real(r8), allocatable  :: bed_thick(:,:)
!        real(r8) :: bedldv(:,:)
#endif

#ifdef BEDLOAD
        real(r8), allocatable :: mcoefd(:,:)
#endif

      END MODULE ocean_mod
