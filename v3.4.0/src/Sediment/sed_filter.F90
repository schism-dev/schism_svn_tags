!=====================================================================
!=====================================================================
! MORSELFE FILTER SUBROUTINES
!
! subroutine sed_avalanching
! subroutine sed_bedmass_filter
!
!=====================================================================
!=====================================================================

      SUBROUTINE sed_avalanching(dhnd)
!--------------------------------------------------------------------!
! This routine updates the bathymetry to account for avalanching     !
! The bed slopes are computed at each node. When a critical slope for!
! element is exceeded, bathymetry of each element node is modified   !
! in order to obtain a slope lower than critical threshold. The      !
! method used here conserves the volume.                             !
! Adapted from filter.f (SAND2D, A. Fortunato)                       !
!                                                                    !
! Currently it does not take account for modification bed sediment   !
! characteristics related to modification of bathymetry              !
!                                                                    !
! Author: florian ganthy (fganthy@lnec.pt ; florian.ganthy@gmail.com)!
! Date:   2013/01/16                                                 !
!                                                                    !
! History:                                                           !
! 2013/01 - F.Ganthy : Modification of wet/dry element consideration !
!                      to be more physical.                          !
! 2013/05 - F.Ganthy : Synchronized processors and mpi_reduced some  !
!                      variables for consistency.                    !
!                      Removed computation of volume control, now    ! 
!                      done in sed_init                              !
!                                                                    !
!--------------------------------------------------------------------!

      USE elfe_glbl, ONLY: dldxy,idry,nea,elnode,npa,rkind,area,xnd,ynd,    &
     &                     errmsg,dp,np
      USE elfe_msgp, ONLY: comm,exchange_p2d,ierr,itype,myrank,      &
     &                     nproc,parallel_abort
      USE sed_mod,   ONLY: dry_slope_cr,wet_slope_cr,vc_area

      IMPLICIT NONE

      INCLUDE 'mpif.h'

!- Arguments --------------------------------------------------------!
      REAL(rkind), DIMENSION(npa), INTENT(inout) :: dhnd
!- Local variables --------------------------------------------------!
      INTEGER :: i,j,iter,iflag,iflag_gb,n1,n2,n3,ndry
      REAL(rkind) :: slope_cr,slope,h1,h2,h3,vec1,vec2,vec3,         &
                     m11,m12,m13,m21,m22,m23,m31,m32,m33,det,        &
                     h1p,h2p,h3p
      REAL(rkind), DIMENSION(:)   :: dp0(npa),dp1(npa),area2(npa),   &
                                     dph(npa)
      REAL(rkind), DIMENSION(:,:) :: dpdxy_el(nea,2)
!- User-defined parameters ------------------------------------------      
      INTEGER,PARAMETER :: maxiter = 10 !Maximum number of iterations
      REAL(rkind), PARAMETER :: epsi  = 0.01d0

!- Start Statement --------------------------------------------------!

!YJZ: Error: this routine will violate bare rock limit
!--------------------------------------------------------------------!
! * Compute bed changes due to sediment transport
!--------------------------------------------------------------------!
      DO i=1,np
        dp0(i) = dp(i)+dhnd(i)
      ENDDO
      dp1 = dp0
      CALL exchange_p2d(dp0)
      CALL exchange_p2d(dp1)

!--------------------------------------------------------------------!
! * Start iterative procedure
!--------------------------------------------------------------------!
      iter     = 0
      iflag    = 1
      iflag_gb = 1
      
      DO WHILE (iflag_gb.EQ.1)
        iflag = 0
        iter  = iter+1
        dpdxy_el = 0.0d0
        DO i=1,nea
!--------------------------------------------------------------------!
! * Compute bed slope at element center
!--------------------------------------------------------------------!
          DO j=1,3
            dpdxy_el(i,1) = dpdxy_el(i,1)+dp1(elnode(j,i))*dldxy(j,1,i)
            dpdxy_el(i,2) = dpdxy_el(i,2)+dp1(elnode(j,i))*dldxy(j,2,i)
          ENDDO
          slope = dsqrt(dpdxy_el(i,1)*dpdxy_el(i,1)+                 &
          &             dpdxy_el(i,2)*dpdxy_el(i,2))
!--------------------------------------------------------------------!
! * Testing for critical slope value (dry or wet)
!   Here we consider that element is dry only if its three nodes are 
!   dry. In the former version, the test was done over idry_e (=1), 
!   but this potentially overestimated the slumping at the wet-dry
!   limit
!--------------------------------------------------------------------!
          ndry = idry(elnode(1,i))+idry(elnode(2,i))+idry(elnode(3,i))
          IF(ndry.EQ.3) THEN
            slope_cr = dry_slope_cr
          ELSE
            slope_cr = wet_slope_cr
          ENDIF

          IF(slope-slope_cr.LE.epsi) CYCLE
!--------------------------------------------------------------------!
! * Preparation of system equation
!--------------------------------------------------------------------!
          n1 = elnode(1,i)
          n2 = elnode(2,i)
          n3 = elnode(3,i)
          h1 = dp1(n1)
          h2 = dp1(n2)
          h3 = dp1(n3)
          m11 = vc_area(n1)
          m12 = vc_area(n2)
          m13 = vc_area(n3)
          m21 = ynd(n2)-ynd(n3)
          m22 = ynd(n3)-ynd(n1)
          m23 = ynd(n1)-ynd(n2)
          m31 = xnd(n3)-xnd(n2)
          m32 = xnd(n1)-xnd(n3)
          m33 = xnd(n2)-xnd(n1)
          vec1 = h1*m11 + h2*m12 + h3*m13
          !slope checked
          vec2 = slope_cr/slope * (h1*m21 + h2*m22 + h3*m23)
          vec3 = slope_cr/slope * (h1*m31 + h2*m32 + h3*m33)
!--------------------------------------------------------------------!
! * Solving the system by Cramer's rule
!--------------------------------------------------------------------!
          det = m11*(m22*m33-m32*m23)-                               &
          &     m12*(m21*m33-m31*m23)+                               &
          &     m13*(m21*m32-m31*m22)

          IF(det.EQ.0.0d0) THEN
            WRITE(errmsg,*)'SED_AVALANCHING: det=0.0'
            CALL parallel_abort(errmsg)
          ENDIF
!--------------------------------------------------------------------!
! * Compute new depth at nodes
!--------------------------------------------------------------------!
          h1p = (vec1*(m22*m33-m32*m23)-                             &
          &     m12*(vec2*m33-vec3*m23)+                             &
          &     m13*(vec2*m32-vec3*m22))/det

          h2p = (m11*(vec2*m33-vec3*m23)-                            &
          &     vec1*(m21*m33-m31*m23)+                              &
          &     m13*(m21*vec3-m31*vec2))/det

          h3p = (m11*(m22*vec3-m32*vec2)-                            &
          &     m12*(m21*vec3-m31*vec2)+                             &
          &     vec1*(m21*m32-m31*m22))/det
!--------------------------------------------------------------------!
! * Apply new depth to temporary bathymetry
!--------------------------------------------------------------------!
          dp1(n1) = h1p
          dp1(n2) = h2p
          dp1(n3) = h3p
          
          iflag = 1
        ENDDO ! End loop on nea
        CALL exchange_p2d(dp1)
        CALL mpi_allreduce(iflag,iflag_gb,1,itype,MPI_MAX,comm,ierr)

        IF (iter.GE.maxiter) THEN
          iflag_gb = 0 !reset flag for exit
          IF(myrank.EQ.0) WRITE(16,*)'Warning: max iterations     &
     &    number reached in sed_avalanching.'
        ENDIF
        
      ENDDO ! End do while on iflag

!--------------------------------------------------------------------!
! * Write number of iteration within mirror.out
!--------------------------------------------------------------------!
      IF(iter > 1.and.myrank==0) THEN
        WRITE(16,*)'# of iter. in sed_avalanching:',iter
      ENDIF
!--------------------------------------------------------------------!
! * Apply depth changes
!--------------------------------------------------------------------!
      DO i=1,np
        dhnd(i) = dhnd(i)+(dp1(i)-dp0(i))
        !IF(isnan(dhnd(i)).EQV..TRUE.) THEN
        IF(dhnd(i)/=dhnd(i)) THEN
          WRITE(errmsg,*) 'Avalanching: dhnd is NaN',myrank,dhnd(i), &
                          dp1(i),dp0(i)
          CALL parallel_abort(errmsg)
        ENDIF
      ENDDO
     CALL exchange_p2d(dhnd)


!--------------------------------------------------------------------!
!      IF(myrank.EQ.0) WRITE(16,*)'Leaving sed_avalanching'
!--------------------------------------------------------------------!
      END SUBROUTINE sed_avalanching

!=====================================================================
!=====================================================================

      SUBROUTINE sed_bedmass_filter(hdep)
!-------------------------------------------------------------------!
! Not quite working at the moment - YJZ
! This subroutine filters iteratively bed_mass for elements where   !
! median grain size (d50, locally computed) presents local extremum !
!(with respect to surrounding elements).                            !
!                                                                   !
! Adapted from filternl.f (SAND2D, A. Fortunato)                    !
!                                                                   !
! - the filter will compute the average with element with           !
!   the largest amplitude                                           !
!                                                                   !
!                     i                    i is a local maximum     !
!                     +                                             !
!                    / \                   the strong filter will   !
!              +----+   \     +-----+                               !
!                  i-1   \   /             use i+1  and the weak    !
!                         \ /                                       !
!                          +               filter will use i-1 to   !
!                          i+1                                      !
!                                          compute an averaged      !
!                                                                   !
!                                          value                    !
!                                                                   !
! Application of filter :                                           !
!  - element is dry                                  : no filter    !
!  - element is wet:                                                !
!             - element have only one wet neighbour  : no filter    !
!             - element have only two wet neighbours : weak filter  !
!             - element have three wet neighbours    : weak/strong  !
!                                                      filter       
!             - element have three wet neighbours    : weak/strong  !
!                                                      filter       !
!                                                                   !
! Bed mass is filtered from a basis of difference (i vs i+1 or i-1) !
! exeed a given threshold (bedmass_threshold, sediment.in)          !
! After filter application, the bed fraction and bed mass are       !
! updated and bed level changes due to these changes in bed mass    !
! are taken into account for morphodynamics                         !
!                                                                   !
! Author: guillaume dodet (gdodet01@univ-lr.fr, gdodet@lnec.pt)     !
! Date:   20/11/2012                                                !
!                                                                   !
! History:                                                          !
! 01/2013 - G.Dodet: Include this subroutine in module sed2d_filter !
! 03/2013 - F.Ganthy: Adapt this subroutine to Sediment (3D)        !
!                     environment, to filter sediment bed_mass      !
!                     at elements instead of depth at nodes from    !
!                     the basis of a threshold exces between a local!
!                     extrema and surrounding value of D50          !
!                                                                   !
!-------------------------------------------------------------------!
      USE elfe_glbl, ONLY : rkind,nea,ic3,idry_e,ntracers,area
      USE elfe_msgp, ONLY : myrank,exchange_e2d,parallel_abort
      USE sed_mod,   ONLY : bed_mass,bed,ithck,Srho,iporo,bed_frac, &
                            nnew,bedmass_filter,bedmass_threshold,  &
                            Sd50

      IMPLICIT NONE
      SAVE

!- Arguments --------------------------------------------------------!
      INTEGER, INTENT(INOUT) :: hdep(nea)

!- Local variables --------------------------------------------------!
      INTEGER :: i,j,iter,iflag,nwse,imin,imax,imin2,imax2,nij,      &
                 ised,jstart
      REAL(rkind) :: d50ij,d50min,d50max,d50min2,d50max2
      REAL(rkind) :: bmij,bmmin,bmmax,bmmed
      REAL(rkind) :: mtot,cff0,cff1
      REAL(rkind), DIMENSION(ntracers)     :: frac
      REAL(rkind), DIMENSION(nea)          :: old_thick,new_thick,d50
      REAL(rkind), DIMENSION(nea,ntracers) :: old_mass

      ! * Maximum number of iterations
      INTEGER, PARAMETER     :: maxiter = 25
      ! * Value to prevent division by zero
      REAL(rkind), PARAMETER :: eps     = 1.0d-14


!- Start Statement --------------------------------------------------!

!---------------------------------------------------------------------
! * Saving original bed_mass [kg.m-2]
!---------------------------------------------------------------------
      DO ised = 1,ntracers
        DO i = 1,nea
          old_mass(i,ised) = bed_mass(1,i,nnew,ised)
        ENDDO
      ENDDO
!---------------------------------------------------------------------
! * Look for local extremum and filter application
!---------------------------------------------------------------------
      iter  = 0
      iflag = 1
      DO WHILE(iflag.EQ.1)
        iflag = 0
!---------------------------------------------------------------------
! * Compute sediment fractions and total D50
!---------------------------------------------------------------------
        DO i = 1,nea
          mtot = 0.0d0
          frac(:)  = 0.0d0
          DO ised = 1,ntracers
            mtot = mtot + bed_mass(1,i,nnew,ised)
          ENDDO ! End loop ntracers
          if(mtot<=0) call parallel_abort('SED3D; sed_filter: mtot<=0')
          DO ised = 1,ntracers
            frac(ised) = bed_mass(1,i,nnew,ised) / mtot
          ENDDO ! End loop ntracers
          cff0 = 1.0d0
          cff1 = 0.0d0
          ! * Weighted geometric mean
          DO ised = 1,ntracers
            cff0 = cff0*Sd50(ised)**frac(ised)
            cff1 = cff1 + frac(ised)
          ENDDO ! End secondary loop ntracers
          if(cff1<=0) call parallel_abort('SED3D; sed_filter: cff1<=0')
          d50(i) = cff0**(1.0d0/cff1)
        ENDDO ! End loop nea
!---------------------------------------------------------------------
! * Starting loop on nea for surrounding elements search and filter
!---------------------------------------------------------------------
        DO i = 1,nea
          ! * DRY ELEMENT : NO FILTER
          IF (idry_e(i).EQ.1) CYCLE
          ! * WET ELEMENTS
          imin = i
          imax = imin
          ! * D50 [m]
          d50min = d50(i)
          d50max = d50min
          ! * NUMBER OF WET SURROUNDING ELEMENTS
          nwse = 0
!---------------------------------------------------------------------
! * Look surrounding element for the strong filter
!---------------------------------------------------------------------
          DO j = 1,3
            IF (ic3(j,i).GT.0) THEN
              nij  = ic3(j,i)
              ! * D50 [m]
              d50ij = d50(nij)
              IF (d50ij.LT.d50min) THEN
                ! * D50IJ IS A LOCAL MINIMUM (over i and surrounding)
                imin   = nij
                d50min = d50ij
              ELSEIF (d50ij.GT.d50max) THEN
                ! * D50IJ IS A LOCAL MAXIMUM  (over i and surrounding)
                imax   = nij
                d50max = d50ij
              ENDIF
              nwse   = nwse + 1 !Number of wet surrounding elements
            ENDIF ! End test on surrounding elements
          ENDDO ! End loop on neighbours
!---------------------------------------------------------------------
! * Look surrounding element for the weak filter
!---------------------------------------------------------------------
          ! * FOR SECOND MIN or MAXIMUM 
          IF (ic3(1,i).GT.0) THEN
            jstart = 1
          ELSEIF (ic3(2,i).GT.0) THEN
            jstart = 2
          ELSEIF (ic3(3,i).GT.0) THEN
            jstart = 3
          ENDIF
          imin2  = i
          imax2  = imin2
          ! * D50 [m]
          d50min2 = d50(i)
          d50max2 = d50min2
          DO j = jstart,3
            IF (ic3(j,i).GT.0) THEN
              nij  = ic3(j,i)
              ! * D50 [m]
              d50ij = d50(nij)
              IF (d50ij.LT.d50min2) THEN
                ! * d50ij IS A LOCAL MINIMUM (over surrounding only)
                imin2   = nij
                d50min2 = d50ij
              ELSEIF (d50ij.GT.d50max2) THEN
                ! * d50ij IS A LOCAL MAXIMUM (over surrounding only)
                imax2   = nij
                d50max2 = d50ij
              ENDIF
            ENDIF ! End test on surrounding elements
          ENDDO ! End loop on neighbours
!---------------------------------------------------------------------
! * Apply filters
!---------------------------------------------------------------------
          d50ij = d50(i)
          ! * TEST FOR FILTER APPLICATION
          IF ((nwse.GT.2).AND.(bedmass_filter.EQ.2)) THEN
!---------------------------------------------------------------------
! * Strong filter
!---------------------------------------------------------------------
            IF ((d50ij.EQ.d50min).AND.                               &
            &   (DABS(d50ij-d50max).GT.bedmass_threshold)) THEN
              ! * d50ij IS LOCAL MINIMUM and THRESHOLD IS EXCEEDED
              iflag = 1
              DO ised = 1,ntracers
                bmij  = bed_mass(1,i,nnew,ised)
                bmmax = bed_mass(1,imax,nnew,ised)
                ! Mass [kg] is conservated
                bmmed = (bmij  * DABS(area(i))+                      &
                &        bmmax * DABS(area(imax))) /                 &
                &       (DABS(area(i)) + DABS(area(imax)))
                ! * APPLICATION OF CALCULATED MEAN
                bed_mass(1,i,nnew,ised)    = bmmed
                bed_mass(1,imax,nnew,ised) = bmmed
              ENDDO ! End loop ntracers
            ENDIF ! End of local minimum
            IF ((d50ij.EQ.d50max).AND.                               &
            &   (DABS(d50ij-d50min).GT.bedmass_threshold)) THEN
              ! * d50ij IS LOCAL MAXIMUM and THRESHOLD IS EXCEEDED
              iflag = 1
              DO ised = 1,ntracers
                bmij  = bed_mass(1,i,nnew,ised)
                bmmin = bed_mass(1,imin,nnew,ised)
                ! Mass [kg] is conservated
                bmmed = (bmij  * DABS(area(i))+                      &
                &        bmmin * DABS(area(imin))) /                 &
                &       (DABS(area(i))+DABS(area(imin)))
                ! * APPLICATION OF CALCULATED MEAN
                bed_mass(1,i,nnew,ised)    = bmmed
                bed_mass(1,imin,nnew,ised) = bmmed
              ENDDO ! End loop ntracers
            ENDIF ! End of local maximum

          ELSEIF (nwse.GT.1) THEN
!---------------------------------------------------------------------
! * Weak filter
!---------------------------------------------------------------------
            IF ((d50ij.EQ.d50min).AND.                               &
            &   (DABS(d50ij-d50min2).GT.bedmass_threshold)) THEN
              ! * d50ij IS LOCAL MINIMUM and THRESHOLD IS EXCEEDED
              iflag = 1
              DO ised = 1,ntracers
                bmij  = bed_mass(1,i,nnew,ised)
                bmmin = bed_mass(1,imin2,nnew,ised)
                ! Mass [kg] is conservated
                bmmed = (bmij  * DABS(area(i))+                      &
                &        bmmin * DABS(area(imin2))) /                &
                &       (DABS(area(i))+DABS(area(imin2)))
                ! * APPLICATION OF CALCULATED MEAN
                bed_mass(1,i,nnew,ised)     = bmmed
                bed_mass(1,imin2,nnew,ised) = bmmed
              ENDDO ! End loop ntracers
            ENDIF ! End of local minimum
            IF ((d50ij.EQ.d50max).AND.                               &
            &   (DABS(d50ij-d50max2).GT.bedmass_threshold)) THEN
              ! * d50ij IS LOCAL MINIMUM and THRESHOLD IS EXCEEDED
              iflag = 1
              DO ised = 1,ntracers
                bmij  = bed_mass(1,i,nnew,ised)
                bmmax = bed_mass(1,imax2,nnew,ised)
                ! Mass [kg] is conservated
                bmmed = (bmij  * DABS(area(i))+                      &
                &        bmmax * DABS(area(imax2))) /                &
                &       (DABS(area(i)) + DABS(area(imax2)))
                ! * APPLICATION OF CALCULATED MEAN
                bed_mass(1,i,nnew,ised)     = bmmed
                bed_mass(1,imax2,nnew,ised) = bmmed
              ENDDO ! End loop ntracers
            ENDIF ! End of local maximum
          ENDIF ! End test for filter application

        ENDDO ! End loop nea
        IF (iflag.EQ.1) iter = iter +1
        IF (iter.GE.maxiter) iflag = 0
      ENDDO ! End while iflag
!---------------------------------------------------------------------
! * End of iterative loop for filter application
!---------------------------------------------------------------------
!      IF(myrank.EQ.0) WRITE(16,*) 'Number of iterations in           &
!     &sed_mass_filter', iter
      WRITE(12,*) 'SED3D: Number of iterations in sed_mass_filter=',iter

      DO ised = 1,ntracers
        CALL exchange_e2d(bed_mass(1,:,nnew,ised))
      ENDDO ! End loop on ntracers

      IF (iter.GT.0) THEN

!---------------------------------------------------------------------
! * Recalculate thickness and fractions for first bed layer
!---------------------------------------------------------------------
        old_thick = 0.0d0
        new_thick = 0.0d0
        DO i = 1,nea
          IF(idry_e(i).EQ.1) CYCLE
          !Get old thickness
          old_thick(i) = bed(1,i,ithck)
          !Calculate the new total mass
          mtot = 0.0d0
          DO ised = 1,ntracers
            mtot = mtot + bed_mass(1,i,nnew,ised)
          ENDDO ! End loop ntracers
          !IF(mtot.EQ.0.0d0) mtot = eps
          mtot=max(mtot,eps)
  
          if(bed(1,i,iporo)==1) call parallel_abort('SED3D; sed_filter: bed(1,i,iporo)==1')
          DO ised = 1,ntracers
            ! Update fraction
            bed_frac(1,i,ised) = bed_mass(1,i,nnew,ised)/mtot
            ! Update new thickness [m]=[kg.m-2]/[kg.m-3]
            bed(1,i,ithck)      = MAX(bed(1,i,ithck)+                  &
            &                         bed_mass(1,i,nnew,ised)/         &
            &                         (Srho(ised)*                     &
            &                          (1.d0-bed(1,i,iporo))),         &
            &                         0.0d0)
          ENDDO ! End loop ntracers
          new_thick(i)        =  bed(1,i,ithck)
!---------------------------------------------------------------------
! * Update bed level changes due to filter application
!---------------------------------------------------------------------
          !Caution, hdep is erosion height in reality so that:
          !hdep>0 if old_thick-new thick>0
          hdep(i) = hdep(i) + (old_thick(i)-new_thick(i))
        ENDDO ! End loop nea

      ENDIF ! End test on iter
!--------------------------------------------------------------------!
      END SUBROUTINE sed_bedmass_filter
