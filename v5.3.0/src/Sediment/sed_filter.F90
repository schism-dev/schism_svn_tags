!   Copyright 2014 College of William and Mary
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

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
! in order to obtain a slope=critical threshold. The      
! method used here conserves the volume.                             !
! Adapted from filter.f (SAND2D, A. Fortunato)                       !
!                                                                    !
! Currently it does not account for modification to bed sediment
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
! Details of algorithm/eqs.
! (0) split quads into pair of tri's
! (1) conservation of volume in an elem.
! \sum(Ai*xi)=\sum(Ai*hi) : Ai area of dual graph, hi is depth at node of elem,
!    and xi is the new hi;
! (2) Force slope to become slope_cr, and preserve the original direction
! \sum(xi*dldxy(,1:2,))=slope_cr*slope_[xy]/slope : slope_[xy] are \nabla h,
! and slope=|\nabla h|

      USE schism_glbl, ONLY: dldxy,idry,nea,i34,elnode,npa,rkind,area,    &
     &                     errmsg,dp,np,xel,yel,nxq
      USE schism_msgp, ONLY: comm,exchange_p2d,ierr,itype,myrank,      &
     &                     nproc,parallel_abort
      USE sed_mod,   ONLY: dry_slope_cr,wet_slope_cr,vc_area

      IMPLICIT NONE

      INCLUDE 'mpif.h'

      REAL(rkind) :: signa

!- Arguments --------------------------------------------------------!
      REAL(rkind), DIMENSION(npa), INTENT(inout) :: dhnd
!- Local variables --------------------------------------------------!
      INTEGER :: i,j,iter,iflag,iflag_gb,n1,n2,n3,ndry,m,jj,nwild(3),nwild2(3)
      REAL(rkind) :: slope_cr,slope,h1,h2,h3,vec1,vec2,vec3,         &
                     m11,m12,m13,m21,m22,m23,m31,m32,m33,det,        &
                     h1p,h2p,h3p,ar2,dldxy2(3,2)
      REAL(rkind), DIMENSION(:)   :: dp0(npa),dp1(npa),area2(npa),   &
                                     dph(npa)
      REAL(rkind), DIMENSION(:,:) :: dpdxy_el(nea,2)
!- User-defined parameters ------------------------------------------      
      INTEGER,PARAMETER :: maxiter = 10 !Maximum number of iterations
      !REAL(rkind), PARAMETER :: epsi  = 0.01d0 !too large??
      REAL(rkind), PARAMETER :: epsi  = 1.e-4

!- Start Statement --------------------------------------------------!

!YJZ: Error: this routine will violate bare rock limit
!--------------------------------------------------------------------!
! * Compute bed changes due to sediment transport
!--------------------------------------------------------------------!
      DO i=1,npa
        dp0(i) = dp(i)+dhnd(i)
      ENDDO
      dp1=dp0

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
          vec1=dot_product(dp1(elnode(1:i34(i),i)),dldxy(1:i34(i),1,i)) !/dx
          vec2=dot_product(dp1(elnode(1:i34(i),i)),dldxy(1:i34(i),2,i)) !/dy
          vec3=sqrt(vec1*vec1+vec2*vec2) !slope

!--------------------------------------------------------------------!
! * Testing for critical slope value (dry or wet)
!   Here we consider that element is dry only if all nodes are 
!   dry. In the former version, the test was done over idry_e (=1), 
!   but this potentially overestimated the slumping at the wet-dry
!   limit
!--------------------------------------------------------------------!
          ndry = sum(idry(elnode(1:i34(i),i)))
          IF(ndry==i34(i)) THEN
            slope_cr = dry_slope_cr
          ELSE
            slope_cr = wet_slope_cr
          ENDIF

          !IF(slope<=slope_cr+epsi) CYCLE
          IF(vec3<=slope_cr+epsi) CYCLE

!         Adjust depths
          iflag = 1
          do m=1,i34(i)-2 !split quads into 2 tri's
            !Calc derivatives of shape function if quads
            if(i34(i)==3) then
              nwild(1:3)=(/1,2,3/) !local indices
              nwild2(1:3)=elnode(nwild(1:3),i) !3 nodes
              dldxy2(1:3,1:2)=dldxy(1:3,1:2,i)
            else !quad
              if(m==1) then !node 1,2,3
                nwild(1:3)=(/1,2,3/) !local indices
              else !1,3,4
                nwild(1:3)=(/1,3,4/)
              endif !m

              nwild2(1:3)=elnode(nwild(1:3),i) !3 nodes of the tri
!              ar2=signa(xnd(nwild2(1)),xnd(nwild2(2)),xnd(nwild2(3)),ynd(nwild2(1)),ynd(nwild2(2)),ynd(nwild2(3)))
              ar2=signa(xel(nwild(1),i),xel(nwild(2),i),xel(nwild(3),i),yel(nwild(1),i),yel(nwild(2),i),yel(nwild(3),i))
              if(ar2<=0) call parallel_abort('SED_FILTER: ar2<=0')
              do jj=1,3
                !Elem. type is 3 not 4!!
                dldxy2(jj,1)=(yel(nwild(nxq(1,jj,3)),i)-yel(nwild(nxq(2,jj,3)),i))/2/ar2 !dL/dx
                dldxy2(jj,2)=(xel(nwild(nxq(2,jj,3)),i)-xel(nwild(nxq(1,jj,3)),i))/2/ar2 !dL/dy
              enddo !jj
            endif !i34

            dpdxy_el(i,1)=dot_product(dp1(nwild2(1:3)),dldxy2(1:3,1)) !dh/dx
            dpdxy_el(i,2)=dot_product(dp1(nwild2(1:3)),dldxy2(1:3,2))
!            DO j=1,i34(i)
!              dpdxy_el(i,1) = dpdxy_el(i,1)+dp1(elnode(j,i))*dldxy(j,1,i)
!              dpdxy_el(i,2) = dpdxy_el(i,2)+dp1(elnode(j,i))*dldxy(j,2,i)
!            ENDDO
            slope=sqrt(dpdxy_el(i,1)*dpdxy_el(i,1)+dpdxy_el(i,2)*dpdxy_el(i,2))

!--------------------------------------------------------------------!
! * Preparation of system equation
!--------------------------------------------------------------------!
            n1 = nwild2(1)
            n2 = nwild2(2)
            n3 = nwild2(3)
            h1 = dp1(n1)
            h2 = dp1(n2)
            h3 = dp1(n3)
            !Matrix
            m11 = vc_area(n1)
            m12 = vc_area(n2)
            m13 = vc_area(n3)
            m21 = dldxy2(1,1) !dL_1/dx
            m22 = dldxy2(2,1) 
            m23 = dldxy2(3,1) 
            m31 = dldxy2(1,2) 
            m32 = dldxy2(2,2) 
            m33 = dldxy2(3,2) 
            vec1 = h1*m11 + h2*m12 + h3*m13 !RHS
            !slope checked
            !vec2 = slope_cr/slope * (h1*m21 + h2*m22 + h3*m23)
            !vec3 = slope_cr/slope * (h1*m31 + h2*m32 + h3*m33)
            vec2 = slope_cr*dpdxy_el(i,1)/slope
            vec3 = slope_cr*dpdxy_el(i,2)/slope
!--------------------------------------------------------------------!
! * Solving the system by Cramer's rule
!--------------------------------------------------------------------!
            det=m11*(m22*m33-m32*m23)-m12*(m21*m33-m31*m23)+m13*(m21*m32-m31*m22)

            IF(det==0.0d0) THEN
              WRITE(errmsg,*)'SED_AVALANCHING: det=0.0'
              CALL parallel_abort(errmsg)
            ENDIF
!--------------------------------------------------------------------!
! * Compute new depth at nodes
!--------------------------------------------------------------------!
            dp1(n1)=(vec1*(m22*m33-m32*m23)-m12*(vec2*m33-vec3*m23)+m13*(vec2*m32-vec3*m22))/det
            dp1(n2)=(m11*(vec2*m33-vec3*m23)-vec1*(m21*m33-m31*m23)+m13*(m21*vec3-m31*vec2))/det
            dp1(n3)=(m11*(m22*vec3-m32*vec2)-m12*(m21*vec3-m31*vec2)+vec1*(m21*m32-m31*m22))/det
          enddo !m=1,i34(i)-2
        ENDDO !i=1,nea
        CALL exchange_p2d(dp1)
        CALL mpi_allreduce(iflag,iflag_gb,1,itype,MPI_MAX,comm,ierr)

        !No CPU dependency as _allreduce is a barrier
        IF(iter>=maxiter) THEN
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
      DO i=1,npa
        dhnd(i) = dhnd(i)+dp1(i)-dp0(i)
        IF(dhnd(i)/=dhnd(i)) THEN
          WRITE(errmsg,*) 'Avalanching: dhnd is NaN',myrank,dhnd(i), &
                          dp1(i),dp0(i)
          CALL parallel_abort(errmsg)
        ENDIF
      ENDDO
!      CALL exchange_p2d(dhnd)

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
      USE schism_glbl, ONLY : rkind,nea,ic3,idry_e,ntrs,area
      USE schism_msgp, ONLY : myrank,exchange_e2d,parallel_abort
      USE sed_mod,   ONLY : bed_mass,bed,ithck,Srho,iporo,bed_frac, &
                           &nnew,bedmass_filter,bedmass_threshold,  &
                           &Sd50,ntr_l

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
      REAL(rkind), DIMENSION(nea)          :: old_thick,new_thick,d50
!      REAL(rkind), DIMENSION(ntracers)     :: frac
!      REAL(rkind), DIMENSION(nea,ntracers) :: old_mass
      REAL(rkind), allocatable :: frac(:),old_mass(:,:)

      ! * Maximum number of iterations
      INTEGER, PARAMETER     :: maxiter = 25
      ! * Value to prevent division by zero
      REAL(rkind), PARAMETER :: eps     = 1.0d-14


!- Start Statement --------------------------------------------------!
      allocate(frac(ntr_l),old_mass(nea,ntr_l),stat=i)
      if(i/=0) call parallel_abort('SED: alloc failed')

!---------------------------------------------------------------------
! * Saving original bed_mass [kg.m-2]
!---------------------------------------------------------------------
      DO ised = 1,ntr_l
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
          DO ised = 1,ntr_l
            mtot = mtot + bed_mass(1,i,nnew,ised)
          ENDDO ! End loop ntr_l
          if(mtot<=0) call parallel_abort('SED3D; sed_filter: mtot<=0')
          DO ised = 1,ntr_l
            frac(ised) = bed_mass(1,i,nnew,ised) / mtot
          ENDDO ! End loop ntr_l
          cff0 = 1.0d0
          cff1 = 0.0d0
          ! * Weighted geometric mean
          DO ised = 1,ntr_l
            cff0 = cff0*Sd50(ised)**frac(ised)
            cff1 = cff1 + frac(ised)
          ENDDO ! End secondary loop ntr_l
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
              DO ised = 1,ntr_l
                bmij  = bed_mass(1,i,nnew,ised)
                bmmax = bed_mass(1,imax,nnew,ised)
                ! Mass [kg] is conservated
                bmmed = (bmij  * DABS(area(i))+                      &
                &        bmmax * DABS(area(imax))) /                 &
                &       (DABS(area(i)) + DABS(area(imax)))
                ! * APPLICATION OF CALCULATED MEAN
                bed_mass(1,i,nnew,ised)    = bmmed
                bed_mass(1,imax,nnew,ised) = bmmed
              ENDDO ! End loop ntr_l
            ENDIF ! End of local minimum
            IF ((d50ij.EQ.d50max).AND.                               &
            &   (DABS(d50ij-d50min).GT.bedmass_threshold)) THEN
              ! * d50ij IS LOCAL MAXIMUM and THRESHOLD IS EXCEEDED
              iflag = 1
              DO ised = 1,ntr_l
                bmij  = bed_mass(1,i,nnew,ised)
                bmmin = bed_mass(1,imin,nnew,ised)
                ! Mass [kg] is conservated
                bmmed = (bmij  * DABS(area(i))+                      &
                &        bmmin * DABS(area(imin))) /                 &
                &       (DABS(area(i))+DABS(area(imin)))
                ! * APPLICATION OF CALCULATED MEAN
                bed_mass(1,i,nnew,ised)    = bmmed
                bed_mass(1,imin,nnew,ised) = bmmed
              ENDDO ! End loop ntr_l
            ENDIF ! End of local maximum

          ELSEIF (nwse.GT.1) THEN
!---------------------------------------------------------------------
! * Weak filter
!---------------------------------------------------------------------
            IF ((d50ij.EQ.d50min).AND.                               &
            &   (DABS(d50ij-d50min2).GT.bedmass_threshold)) THEN
              ! * d50ij IS LOCAL MINIMUM and THRESHOLD IS EXCEEDED
              iflag = 1
              DO ised = 1,ntr_l
                bmij  = bed_mass(1,i,nnew,ised)
                bmmin = bed_mass(1,imin2,nnew,ised)
                ! Mass [kg] is conservated
                bmmed = (bmij  * DABS(area(i))+                      &
                &        bmmin * DABS(area(imin2))) /                &
                &       (DABS(area(i))+DABS(area(imin2)))
                ! * APPLICATION OF CALCULATED MEAN
                bed_mass(1,i,nnew,ised)     = bmmed
                bed_mass(1,imin2,nnew,ised) = bmmed
              ENDDO ! End loop ntr_l
            ENDIF ! End of local minimum
            IF ((d50ij.EQ.d50max).AND.                               &
            &   (DABS(d50ij-d50max2).GT.bedmass_threshold)) THEN
              ! * d50ij IS LOCAL MINIMUM and THRESHOLD IS EXCEEDED
              iflag = 1
              DO ised = 1,ntr_l
                bmij  = bed_mass(1,i,nnew,ised)
                bmmax = bed_mass(1,imax2,nnew,ised)
                ! Mass [kg] is conservated
                bmmed = (bmij  * DABS(area(i))+                      &
                &        bmmax * DABS(area(imax2))) /                &
                &       (DABS(area(i)) + DABS(area(imax2)))
                ! * APPLICATION OF CALCULATED MEAN
                bed_mass(1,i,nnew,ised)     = bmmed
                bed_mass(1,imax2,nnew,ised) = bmmed
              ENDDO ! End loop ntr_l
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

      DO ised = 1,ntr_l
        CALL exchange_e2d(bed_mass(1,:,nnew,ised))
      ENDDO ! End loop on ntr_l

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
          DO ised = 1,ntr_l
            mtot = mtot + bed_mass(1,i,nnew,ised)
          ENDDO ! End loop ntr_l
          !IF(mtot.EQ.0.0d0) mtot = eps
          mtot=max(mtot,eps)
  
          if(bed(1,i,iporo)==1) call parallel_abort('SED3D; sed_filter: bed(1,i,iporo)==1')
!'
          DO ised = 1,ntr_l
            ! Update fraction
            bed_frac(1,i,ised) = bed_mass(1,i,nnew,ised)/mtot
            ! Update new thickness [m]=[kg.m-2]/[kg.m-3]
            bed(1,i,ithck)      = MAX(bed(1,i,ithck)+                  &
            &                         bed_mass(1,i,nnew,ised)/         &
            &                         (Srho(ised)*                     &
            &                          (1.d0-bed(1,i,iporo))),         &
            &                         0.0d0)
          ENDDO ! End loop ntr_l
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

      deallocate(frac,old_mass)

      END SUBROUTINE sed_bedmass_filter
