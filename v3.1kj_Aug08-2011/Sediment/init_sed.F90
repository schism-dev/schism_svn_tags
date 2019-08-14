!=======================================================================
             
!=======================================================================
     SUBROUTINE initialize_scalars
!                                                                      !
!  This routine initializes several variables in module "mod_scalars"  !
!  for all nested grids.                                               !
!                                                                      !
!=======================================================================
!Ligia Pinto                                  	                       !
!  Kept only the ones with respect to sediment                         !
!                                                                      !
!=======================================================================

!
!  Local variable declarations.
!
      use sed_param, only: isand,Nbed
      use elfe_glbl, only: ntracers
      use elfe_msgp, only: parallel_abort
      implicit none  
      integer :: i, ic, j


      allocate ( isand(MAX(1,ntracers)), stat=i )
      if(i/=0) call parallel_abort('Main: allocation failure')
!
!---------------------------------------------------------------------
!  Set tracer identification indices.
!---------------------------------------------------------------------
!
! Marta Rodrigues
! In SELFE T and S are in diferent arrays from the tracers array...
       
       ic=0

!
!  Set noncohesive suspended sediment tracers indices
!
!      NST=ntracers
      DO i=1,ntracers
        ic=ic+1
        isand(i)=ic
      END DO

      END SUBROUTINE initialize_scalars    

      SUBROUTINE initialize_ocean
!
!=======================================================================
!                                                                      !
!  This routine initialize all variables in the module using first     !
!  touch distribution policy. In shared-memory configuration, this     !
!  operation actually performs propagation of the  "shared arrays"     !
!  across the cluster, unless another policy is specified to           !
!  override the default.                                               !
!                                                                      !
!=======================================================================
!
      use sed_param, only: Nbed,r8
      USE sed_mod, only: MBEDP,MBOTP
      USE ocean_mod
      USE elfe_glbl, only: nvrt,nea,npa,ntracers,mnei
      use elfe_msgp, only: myrank,parallel_abort
      implicit none
!

!  Local variable declarations.
!
      integer :: i, j
      integer :: itrc, k

      real(r8), parameter :: IniVal = 0.0_r8
!
!     Allocate arrays in ocean_mod
      allocate(Hz(nvrt,nea),stat=i)
      if(i/=0) call parallel_abort('main: Hz allocation failure')
      allocate(bed(Nbed,nea,MBEDP),stat=i)
      if(i/=0) call parallel_abort('main: bed allocation failure')
      allocate(bed_frac(Nbed,nea,ntracers),stat=i)
      if(i/=0) call parallel_abort('main: Nbed allocation failure')
      allocate(bed_mass(Nbed,nea,2,ntracers),stat=i)
      if(i/=0) call parallel_abort('main: bed_mass allocation failure')
      allocate(bottom(nea,MBOTP),stat=i)
      if(i/=0) call parallel_abort('main: bottom allocation failure')
      allocate(bedldu(npa,ntracers),stat=i)
      if(i/=0) call parallel_abort('main: bedldu allocation failure')
      allocate(bedldv(npa,ntracers),stat=i)
      if(i/=0) call parallel_abort('main: bedldv allocation failure')
#ifdef SED_MORPH 
     allocate(bed_thick(nea,2),stat=i)
     if(i/=0) call parallel_abort('main: bed_thick allocation failure')
#endif
     
!  Set array initialization range.
!
!
!-----------------------------------------------------------------------
!  Initialize module variables.
!-----------------------------------------------------------------------
        DO i=1,nea
            DO itrc=1,ntracers
              DO k=1,Nbed
                bed_frac(k,i,itrc) = IniVal
                bed_mass(k,i,1,itrc) = IniVal
                bed_mass(k,i,2,itrc) = IniVal
              END DO
            END DO
        END DO
        DO i=1,nea
            DO itrc=1,MBEDP
              DO k=1,Nbed
                 bed(k,i,itrc) = IniVal
              END DO
            END DO
       END DO
        DO i=1,nea
            DO itrc=1,MBOTP
               bottom(i,itrc) = IniVal
            END DO
        END DO
#ifdef BEDLOAD
          DO i=1,npa
            DO itrc=1,ntracers
              bedldu(i,itrc) = IniVal
              bedldv(i,itrc) = IniVal
            END DO
          END DO
#endif

      END SUBROUTINE initialize_ocean
