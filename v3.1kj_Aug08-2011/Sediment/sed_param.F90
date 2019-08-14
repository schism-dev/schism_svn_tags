      MODULE sed_param
      
!!======================================================================
!! August, 2007                                                        ! 
!!======================================================Ligia Pinto=====
!!                                                                     !
!! This module is adapted from ROMS (mod_param.F + mod_scalars.F):     !
!!                                                                     ! 
!!======================================================================   
!
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!                                                                      !
!  Tracer parameters:                                                  !
!                                                                      !
!  NST        Number of non-cohesive (sand) sediment tracers.          !
!  Nbed       Number of sediment bed layers.                           !
!
!                                                                      !
!=======================================================================
!

        use elfe_glbl, only: ntracers
        IMPLICIT NONE
        SAVE
         
        integer, parameter :: r4 = 4 
        integer, parameter :: r8 = 8  
        integer, parameter :: c8 = selected_real_kind(6,30)        
!
!-----------------------------------------------------------------------
!  Sediment tracers parameters.
!-----------------------------------------------------------------------
!
!  Number of sediment bed layers.
!
        integer, parameter :: Nbed = 1 
!
! Error: set NST to ntracers? Can even read in NST, NCS, NNS from param.in
!
!  Number of non-cohesive (sand) sediments.
!
!        integer :: NST = 1 
!
!-----------------------------------------------------------------------
!  Tracer identification indices.
!-----------------------------------------------------------------------
!

! Temperature and salinity from SELFE
!        integer :: itemp              ! Potential temperature
!        integer :: isalt              ! Salinity

! Sediment	
        integer, allocatable :: isand(:)  ! Non-cohesive sediment

!-----------------------------------------------------------------------
!  Physical constants.   
!-----------------------------------------------------------------------
!
!    Cp            Specific heat for seawater (Joules/Kg/degC).
!    Csolar        Solar irradiantion constant (W/m2).
!    Eradius       Earth equatorial radius (m).
!    StefBo        Stefan-Boltzmann constant (W/m2/K4).
!    emmiss        Infrared emmissivity.
!    g             Acceleration due to gravity (m/s2).
!    gorho0        gravity divided by mean density anomaly.
!    rhow          fresh water density (kg/m3).
!    vonKar        von Karman constant.
!
!        real(r8) :: rhow = 998.0_r8             ! kg/m3         
        real(r8) :: g = 9.81_r8                 ! m/s2
        real(r8) :: vonKar = 0.40_r8            ! non-dimensional
!

!  Minimum and maximum threshold for transfer coefficient of momentum.
!
        real(r8) :: Cdb_min = 0.000001_r8
        real(r8) :: Cdb_max = 0.5_r8


!  Mean density (Kg/m3) used when the Boussinesq approximation is
!  inferred. In SELFE rho0 is the fresh water density
!
           real(r8) :: rhom = 1000.0_r8

      END MODULE sed_param
