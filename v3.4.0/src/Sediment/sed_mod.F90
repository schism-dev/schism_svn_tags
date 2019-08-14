      MODULE sed_mod
!--------------------------------------------------------------------!
! Variables declaration for 3D sediment model                        !
!                                                                    !
! This subroutine is adapted from ROMS routine ana_sediment.h        !
! Copyright (c) 2002-2007 The ROMS/TOMS Group                        !
!   Licensed under a MIT/X style license                             !
!   See License_ROMS.txt                                             !
!                                                                    !
! Author: Ligia Pinto                                                !
! Date: xx/08/2007                                                   !
!                                                                    !
! History: 2012/12 - F.Ganthy : form homogenisation of sediments     !
!                               routines                             !
!          2012/12 - F.Ganthy : modifications for Bottom Composition !
!                               Generation (BCG) purpose (added      !
!                               output files - bed median grain size,!
!                               bed fraction)                        !
!          2013/01 - F.Ganthy : Implementation of roughness predictor!
!          2013/01 - F.Ganthy : Implementation of avalanching        !
!          2013/03 - F.Ganthy : Implementation of wave-induced       !
!                               bedload transport                    !
!          2013/04 - F.Ganthy : Implementation of wave-current bottom!
!                               stress                               !
!          2013/05 - F.Ganthy : Updates related to ripple predictor  !
!          2013/05 - F.Ganthy : Added volume control area            !
!          2013/05 - F.Ganthy : Updates to the ripple predictor:     !
!                               - Changes on the total bedform       !
!                                 roughness computation              !
!                               - Add wave-ripple computation from   !
!                                 Nielsen (1992)                     !
!          2013/05 - F.Ganthy : Add different sediment behavior:     !
!                                - MUD-like or SAND-like             !
!          2013/05 - F.Ganthy : Re-introduction of the number of bed !
!                               sediment layers within sediment.in   !
!                                                                    !
!--------------------------------------------------------------------!
!                                                                    !
!  Parameters for sediment model:                                    !
!                                                                    !
!   Csed     Sediment concentration (kg/m3), used during analytical  !
!              initialization.                                       !
!   Erate    Surface erosion rate (kg/m2/s).                         !
!   Sd50     Median sediment grain diameter (m).                     !
!   Srho     Sediment grain density (kg/m3).                         !
!   SedIter  Maximum number of iterations.                           !
!   Wsed     Particle settling velocity (m/s).                       !
!   poros    Porosity (non-dimensional: 0.0-1.0):                    !
!              Vwater/(Vwater+Vsed).                                 !
!   tau_ce   Kinematic critical shear for erosion (m2/s2).           !
!   tau_cd   Kinematic critical shear for deposition (m2/s2).        !
!   morph_fac  Morphological scale factor (nondimensional).          !
!   newlayer_thick    New layer deposit thickness criteria (m).      !
!   bedload_coeff     Bedload rate coefficient (nondimensional).     !
!                                                                    !
!  BED properties indices:                                           !
!                                                                    !
!   MBEDP    Number of bed properties (array dimension).             !
!   idBmas   Sedimen mass index.                                     !
!   idSbed   IO indices for bed properties variables.                !
!   idfrac   sediment class fraction (non-dimensional).              !
!   ithck    Index into sediment layer thickness (m). 
!   iaged    Index into sediment layer age (s).                      
!   iporo    Index into sediment layer porosity (non-dimensional).    
!   idiff    Sediment layer bio-diffusivity (m2/s).                  !
!                                                                    !
!  BOTTOM properties indices:                                        !
!                                                                    !
!   MBOTP    Number of bottom properties (array dimension).          !
!   idBott   IO indices for bottom properties variables.             !
!   isd50    (Index into) Median sediment grain diameter (m).                     !
!   idens    Median sediment grain density (kg/m3).                  !
!   iwsed    Mean settling velocity (m/s).                           !
!   itauc    Mean critical erosion stress (m2/s2).                   !
!   irlen    Sediment ripple length (m).                             !
!   irhgt    Sediment ripple height (m).                             !
!   ibwav    Bed wave excursion amplitude (m).                       !
!   izdef    Default bottom roughness (m).                           !
!   izapp    Apparent bottom roughness (m).                          !
!   izNik    Nikuradse bottom roughness (m).                         !
!   izbio    Biological bottom roughness (m).                        !
!   izbfm    Bed form bottom roughness (m).                          !
!   izbld    Bed load bottom roughness (m).                          !
!   izwbl    Bottom roughness used wave BBL (m).                     ! 
!   iactv    Active layer thickness for erosive potential (m).       !
!   ishgt    Sediment saltation height (m).                          !
!   idoff    Offset for calculation of dmix erodibility profile (m)  !
!   idslp    Slope for calculation of dmix or erodibility profile    !
!   idtim    Time scale for restoring erodibility profile (s)        !
!   idbmx    bed biodifusivity max                                   !
!   idbmm    bed biodifusivity m                                     !
!   idbzs    bed biodifusivity zs                                    !
!   idbzm    bed biodifusivity zm                                    !
!   idbzp    bed biodifusivity phi                                   !
!                                                                    !
!--------------------------------------------------------------------!

       USE elfe_glbl, ONLY: rkind

       IMPLICIT NONE
       SAVE

!- Sediment variables -----------------------------------------------!

       INTEGER, PARAMETER :: MBEDP = 3  ! Bed Properties:
       INTEGER, PARAMETER :: ithck = 1  ! 1st property is layer thickness
       INTEGER, PARAMETER :: iaged = 2  ! 2nd property is layer age
       INTEGER, PARAMETER :: iporo = 3  ! 3rd property is layer porosity
!       INTEGER, PARAMETER :: idiff = 4  ! layer bio-diffusivity

       INTEGER, PARAMETER :: MBOTP = 12  ! Bottom Properties:
        
       INTEGER, PARAMETER :: isd50 = 1  ! 1st property is mean grain diameter
       INTEGER, PARAMETER :: idens = 2  ! 2nd property is mean grain density
       INTEGER, PARAMETER :: iwsed = 3  ! mean settle velocity
       INTEGER, PARAMETER :: itauc = 4  ! critical erosion stress
       INTEGER, PARAMETER :: iactv = 5  ! active layer thickness
       INTEGER, PARAMETER :: izdef = 6  ! default bottom roughness length (rough.gr3)
       INTEGER, PARAMETER :: izNik = 7  ! Nikuradse bottom roughness length (D50/12)
       INTEGER, PARAMETER :: izcr  = 8  ! Current ripple roughness length (Soulsby, 1997)
       INTEGER, PARAMETER :: izsw  = 9  ! Sand waves roughness length (Van Rijn, 1984)
       INTEGER, PARAMETER :: izwr  = 10 ! Wave ripples roughness length (Grant and Madsen, 1982 OR Nielsen, 1992)
       INTEGER, PARAMETER :: izbld = 11 ! Bed load bottom roughness length (Grant and Madsen, 1982 OR Nielsen, 1992)
       INTEGER, PARAMETER :: izapp = 12 ! Apparent (total) bottom roughnes
!       INTEGER, PARAMETER :: ibwav = 7  ! wave excursion amplitude
!       INTEGER, PARAMETER :: izbio = 11 ! biological bottom roughness
!       INTEGER, PARAMETER :: izwbl = 14 ! wave bottom roughness
!       INTEGER, PARAMETER :: ishgt = 16 ! saltation height

!       INTEGER  :: idBott(MBOTP)        ! bottom properties IDs
!       INTEGER  :: idSbed(MBEDP)        ! bed properties IDs

!       INTEGER, ALLOCATABLE :: idBmas(:)   ! class mass indices
!       INTEGER, ALLOCATABLE :: idfrac(:)   ! class fraction indices
       INTEGER :: sed_debug  ! used for sediment model debugging, outputs lots of variables to mirror out
       INTEGER :: bedload                ! activation key for bedload transport
       INTEGER :: suspended_load         ! activation key for suspended load tranport
       INTEGER :: slope_formulation      ! activation key for slope effects on beldload
       INTEGER :: bc_for_weno            ! activation of boundary condition for weno
       INTEGER :: sed_morph              ! activation key for morphodynamics
       INTEGER :: drag_formulation       ! key for drag fourmulation method
       INTEGER :: ddensed                ! activation key for sediment density effects on water density

       INTEGER :: comp_ws                ! activation of Soulsby settling velocity
       INTEGER :: comp_tauce             ! activation of Soulsby critical shear stress
       INTEGER :: bedforms_rough         ! activation of roughness predictor
       INTEGER :: iwave_ripple           ! wave ripple computed from Grant and Madsen (1982) or Nielsen (1992)
       INTEGER :: irough_bdld            ! activation of the bedload transport induced roughness
       INTEGER :: slope_avalanching      ! activation of the avalanching computation
       INTEGER :: bedmass_filter         ! activation of the bedmass filter


       INTEGER :: nstp,nnew

!----------------------
! coming from sed_param
       INTEGER, PARAMETER :: r4 = 4 
       INTEGER, PARAMETER :: r8 = 8  
       INTEGER, PARAMETER :: c8 = selected_real_kind(6,30)

!  Number of sediment bed layers.
!       INTEGER, PARAMETER :: Nbed = 1
       INTEGER :: Nbed !fixed thru'out run (by combining bottom-most 2 layers if necessary)
       INTEGER, ALLOCATABLE :: isand(:)  ! Non-cohesive sediment indices into tr_el; basically= 1:ntracers (isand(MAX(1,ntracers))

! user specified constants
! will be read in from sediment.in
       REAL(rkind) :: newlayer_thick          ! deposit thickness criteria
       REAL(rkind) :: bedload_coeff           ! bedload rate coefficient
       REAL(rkind) :: porosity                ! bed sediment porosity
       REAL(rkind) :: bdldiffu                ! bedload diffusivity coef. 
       REAL(rkind) :: dry_slope_cr            ! critical slope for dry nods
       REAL(rkind) :: wet_slope_cr            ! critical slope for wet nods
       REAL(rkind) :: bedmass_threshold       ! threshold for bedmass_filter

       REAL(rkind), ALLOCATABLE  :: bedthick_overall(:)        ! bed thickness; bedthick_overall(npa)
       REAL(rkind), ALLOCATABLE :: Csed(:)      ! initial concentration; not used
       REAL(rkind), ALLOCATABLE :: Erate(:)     ! erosion rate (>0); Erate(ntracers)
       REAL(rkind), ALLOCATABLE :: Sd50(:)      ! mediam grain diameter; Sd50(ntracers)
       REAL(rkind), ALLOCATABLE :: Srho(:)      ! grain density; Srho(ntracers)
       REAL(rkind), ALLOCATABLE :: Wsed(:)      ! settling velocity (>0); Wsed(ntracers)
       REAL(rkind), ALLOCATABLE :: poros(:)     ! porosity \in [0,1]; not used at the moment
       REAL(rkind), ALLOCATABLE :: tau_ce(:)    ! shear for erosion (>0); tau_ce(ntracers)
       REAL(rkind), ALLOCATABLE :: Sedtype(:)   ! Sediment type; Sedtype(ntracers)
!        REAL(rkind), ALLOCATABLE :: tau_cd(:)   ! shear for deposition - not used
       REAL(rkind), ALLOCATABLE :: morph_fac(:) ! morphological factor; morph_fac(ntracers)

! Model parameters for a single layer bed model - these are no longer used 
! all below variables have dimension of (nea,ntracers)
!       REAL(rkind), ALLOCATABLE :: bedthick(:,:) ! 
!       REAL(rkind), ALLOCATABLE :: bedfrac(:,:)  ! 
!       REAL(rkind), ALLOCATABLE :: bedmass(:,:)  ! 
!       REAL(rkind), ALLOCATABLE :: bedporo(:,:)  ! 

! Acceleration due to gravity (m/s2)
       REAL(rkind) :: g
! von Karman constant
       REAL(rkind) :: vonKar
! Minimum and maximum threshold for transfer coefficient of momentum.
       REAL(rkind) :: Cdb_min
       REAL(rkind) :: Cdb_max
! Mean density (Kg/m3) used when the Boussinesq approximation is
! inferred. In SELFE rho0 is the fresh water density
       REAL(rkind) :: rhom

       REAL(rkind), ALLOCATABLE :: Hz(:,:) !Hz(nvrt,nea)
       REAL(rkind), ALLOCATABLE :: Zob(:)
       REAL(rkind), ALLOCATABLE :: bed(:,:,:) !bed(Nbed,nea,MBEDP): properties
       REAL(rkind), ALLOCATABLE :: bed_frac(:,:,:) !bed_frac(Nbed,nea,ntracers) - sum over ntracers should =1
       REAL(rkind), ALLOCATABLE :: bed_mass(:,:,:,:) !bed_mass(Nbed,nea,2,ntracers)>=0; 
                                                     !'2' is used to store old and new time step alternately
                                                     ![kg/m^2] (= rho*thick)
       REAL(rkind), ALLOCATABLE :: bottom(:,:) !bottom(nea,MBOTP): bottom (water-sed interface) properties
       REAL(rkind), ALLOCATABLE :: bedldu(:,:) !bedldu(npa,ntracers)
       REAL(rkind), ALLOCATABLE :: bedldv(:,:) !bedldv(npa,ntracers)

       REAL(rkind), ALLOCATABLE :: bed_thick(:) !bed_thick(nea) - sum over all layers; not really used
       REAL(rkind), ALLOCATABLE :: mcoefd(:,:) !JCG matrix

       ! Variables shared by subroutines sediment, bedload_vr
       LOGICAL, ALLOCATABLE     :: lbc_sed(:) !b.c. flag for erosion eq., used in JCG
       REAL(rkind), ALLOCATABLE :: bc_sed(:)  !b.c. for erosion eq., used in JCG

       REAL(rkind) :: smgd
       REAL(rkind), ALLOCATABLE  :: FX_r(:) !bedload flux (FX_r(nea))
       REAL(rkind), ALLOCATABLE  :: FY_r(:)
       REAL(rkind) :: angleu,anglev
       REAL(rkind) :: dzdx, dzdy
       REAL(rkind) :: sed_angle

       ! Variable shared by sed_friction subroutines
       REAL(rkind), ALLOCATABLE  :: bustr(:)       !Current-induced bottom stress in direction x (m2.s-2)
       REAL(rkind), ALLOCATABLE  :: bvstr(:)       !Current-induced bottom stress in direction y (m2.s-2)
       REAL(rkind), ALLOCATABLE  :: tau_c(:)       !Current-induced bottom stres (m2.s.2)
       REAL(rkind), ALLOCATABLE  :: tau_w(:)       !Wave-induced bottom stres (m2.s.2)
       REAL(rkind), ALLOCATABLE  :: tau_wc(:)      !Wave-current mean bottom stress (m2.s-2)

       ! Bed variables defined at nodes
       REAL(rkind), ALLOCATABLE :: vc_area(:)      !Volume control area (vc_area(npa))
       REAL(rkind), ALLOCATABLE :: bed_d50n(:)     !Median grain size (m)
       REAL(rkind), ALLOCATABLE :: bed_fracn(:,:)  !Sediment fraction (0-1) for each class (bed_fracn(npa,ntracers))
       REAL(rkind), ALLOCATABLE :: bed_taun(:)     !Bottom shear stress (N.m-2)
       REAL(rkind), ALLOCATABLE :: bed_rough(:)    !Apparent Roughness length (bedform prediction)

       ! WWM variables defined at element centres
       REAL(rkind), ALLOCATABLE :: hs(:)     !Significant wave height from WWM (m) ; (hs(nea)
       REAL(rkind), ALLOCATABLE :: tp(:)     !Peak wave period from WWM (s); (tp(nea))
       REAL(rkind), ALLOCATABLE :: wlpeak(:) !Wave lenght associated with peak period from WWM (m); wlpeak(nea)
       REAL(rkind), ALLOCATABLE :: uorb(:)   !RMS Orbital velocity from WWM (m.s-1); uorb(nea)
       REAL(rkind), ALLOCATABLE :: uorbp(:)  !Peak orbital velocity from WWM (m.s-1); uorbp(nea) - not really used

!--------------------------------------------------------------------!
       END MODULE sed_mod
