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

SUBROUTINE sed2d_morpho(it)
!--------------------------------------------------------------------
! This subroutine computes bottom evolution based on Exner equation
!
! Adapted from sediment_v8.F90 (MORSELFE, L. Pinto) 
!
! Author: guillaume dodet (gdodet01@univ-lr.fr, gdodet@lnec.pt)
! Date:   06/12/2012
!
! History:
! 02/2013 - G.Dodet: - Removed arguments shared in schism_glbl; 
!                    - Removed useless ghost exchange; 
! 03/2013 - G.Dodet: - Corrected argument dimension for solve_jcg;
! 04/2013 - G.Dodet: - Added element-centered method;
!                    - Added information for debugging;
! 07/2013 - G.Dodet: - Flux time integration is now done in 
!                      sed2d_main based on dtsed2d
! 03/2014 - T.Guerin: - Adapted the routine for mutli-class multi-
!                       layer mode: computation of new bathymetry is
!                       done in sed2d_main_mcml for this mode
!--------------------------------------------------------------------

  USE schism_glbl, ONLY : area,dt,dp,eta2,idry,idry_e,indel,iplg,isdel,     &
                        isidenode,elside,nsa,moitn0,mxitn0,nea,i34,elnode,nne,np,&
                        npa,rkind,rtol0,xctr,xnd,yctr,ynd
  USE schism_msgp, ONLY : exchange_e2d,exchange_p2d,parallel_abort
  USE sed2d_mod, ONLY : bc_flag,bc_val,bed_delta,bed_del_e,h0_sed,   &
!YJZ: qtot_e not used; also consider removing imeth=2 altogether, and qtot_s
                        imeth,mcoefd,poro,qdt_e,qtot_s,nb_class

  IMPLICIT NONE

!- Arguments --------------------------------------------------------
  INTEGER, INTENT(IN) :: it          ! time step
!- Local variables --------------------------------------------------
  INTEGER :: i,iside,j  
  REAL(rkind) :: area_tot,flux,htot,x21,xp,y21,yp,tmp2 
  REAL(rkind), DIMENSION(nea,2) :: qtotdt_e
  REAL(rkind), DIMENSION(nsa,2) :: qtotdt_s
  REAL(rkind), DIMENSION(npa) :: qsan2d
  REAL(rkind), DIMENSION(np) :: tmp
!--------------------------------------------------------------------

  IF (imeth == 1) THEN !Node-centered method

!- Compute bottom update --------------------------------------------
!- Integrate transport flux (m^2) along node-centered control volume 
!- and solve equation system with jcg solver to obtain bed delta (m).
!--------------------------------------------------------------------
  qsan2d = 0.d0
  DO i = 1,nea
     htot = 0.d0
     DO j = 1,i34(i)
        htot = htot+(eta2(elnode(j,i))+dp(elnode(j,i)))/i34(i)
     ENDDO !j

     IF(idry_e(i) == 1.OR.htot.LE.h0_sed) CYCLE

!Error: YJZ - need to account for quads
     xp  = (xnd(elnode(2,i))+xnd(elnode(3,i)))/2.d0
     yp  = (ynd(elnode(2,i))+ynd(elnode(3,i)))/2.d0
     flux= qdt_e(i,1)*(yp-yctr(i))-qdt_e(i,2)*(xp-xctr(i))
     qsan2d(elnode(2,i)) = qsan2d(elnode(2,i))-flux
     qsan2d(elnode(3,i)) = qsan2d(elnode(3,i))+flux

     xp  = (xnd(elnode(3,i))+xnd(elnode(1,i)))/2.d0
     yp  = (ynd(elnode(3,i))+ynd(elnode(1,i)))/2.d0
     flux= qdt_e(i,1)*(yp-yctr(i))-qdt_e(i,2)*(xp-xctr(i))
     qsan2d(elnode(3,i)) = qsan2d(elnode(3,i))-flux
     qsan2d(elnode(1,i)) = qsan2d(elnode(1,i))+flux

     xp  = (xnd(elnode(1,i))+xnd(elnode(2,i)))/2.d0
     yp  = (ynd(elnode(1,i))+ynd(elnode(2,i)))/2.d0
     flux= qdt_e(i,1)*(yp-yctr(i))-qdt_e(i,2)*(xp-xctr(i))
     qsan2d(elnode(1,i)) = qsan2d(elnode(1,i))-flux
     qsan2d(elnode(2,i)) = qsan2d(elnode(2,i))+flux
  ENDDO !nea
  CALL exchange_p2d(qsan2d)

  DO i = 1,np
     tmp(i) = qsan2d(i)/(1.d0-poro)
  ENDDO !np

!- Solve Exner equation ---------------------------------------------
  bed_delta = 0.d0 
  CALL solve_jcg(it,moitn0,mxitn0,rtol0,mcoefd,bed_delta,tmp,bc_val,bc_flag)
  CALL exchange_p2d(bed_delta)

!- Compute new bathymetry -------------------------------------------
!YJZ: the following should be removed once two main routines r merged
  IF (nb_class==1) THEN
    DO i = 1,npa
       IF(idry(i)==1) CYCLE
       dp(i) = dp(i)+bed_delta(i)
       IF(dp(i)/=dp(i) .OR. bed_delta(i)>1.d0) THEN
         WRITE(12,*)'Warning (3) for it:',it,' and node: ', iplg(i)
         WRITE(12,*)'dp = ',dp(i)
         WRITE(12,*)'bed delta = ',bed_delta(i)
         CALL parallel_abort('Sed2d: depth is NaN or strange bed_delta')    
       ENDIF  
    ENDDO

    CALL exchange_p2d(dp)
  ENDIF

!--------------------------------------------------------------------
  ELSE !IMETH = 2 Element-centered method

!- Compute bottom update --------------------------------------------
!- Integrate transport fluxes (m^2) along element sides and divide by 
!- element surface area x porosity to obtain bed delta (m) at element
!- center.
!- Bed delta at nodes is obtained through node weighted average. 
!--------------------------------------------------------------------
  qtotdt_s(:,:) = qtot_s(:,:)*dt
  bed_del_e = 0.d0
  DO i = 1,nea
     tmp2 = 0.d0
     htot = 0.d0
     DO j = 1,i34(i)
        htot = htot+(eta2(elnode(j,i))+dp(elnode(j,i)))/i34(i)
     ENDDO !j
     IF((idry_e(i) == 1).OR.(htot.LE.h0_sed)) CYCLE
     DO j = 1,i34(i)
        iside = elside(i,j)
        IF(i == isdel(1,iside)) THEN
          x21 = xnd(isidenode(2,iside))-xnd(isidenode(1,iside))
          y21 = ynd(isidenode(2,iside))-ynd(isidenode(1,iside))
        ELSE
          x21 = xnd(isidenode(1,iside))-xnd(isidenode(2,iside))
          y21 = ynd(isidenode(1,iside))-ynd(isidenode(2,iside))
        ENDIF
        tmp2 = tmp2+qtotdt_s(iside,1)*y21-qtotdt_s(iside,2)*x21
     ENDDO
     bed_del_e(i) = tmp2/((area(i)*(1.d0-poro)))
  ENDDO !nea
  CALL exchange_e2d(bed_del_e)

  bed_delta = 0.d0
  DO i = 1,np
     area_tot = 0.d0
     DO j = 1,nne(i)
        bed_delta(i) = bed_delta(i)+area(indel(j,i))*bed_del_e(indel(j,i))
        area_tot = area_tot+area(indel(j,i))
     ENDDO
     IF(idry(i) == 1 .OR. bc_flag(i)) CYCLE
     dp(i) = dp(i)+bed_delta(i)/area_tot
!     IF(isnan(dp(i)) .EQV. .TRUE. .OR. bed_delta(i)>1.d0) THEN
       IF(dp(i) /= dp(i) .OR. bed_delta(i)>1.d0) THEN
       WRITE(12,*)'Warning (3) for it:',it,' and node: ', iplg(i)
       WRITE(12,*)'dp = ',dp(i)
       WRITE(12,*)'bed delta = ',bed_delta(i)
     ENDIF
  ENDDO
  CALL exchange_p2d(dp)

  ENDIF !IMETH

END SUBROUTINE sed2d_morpho
