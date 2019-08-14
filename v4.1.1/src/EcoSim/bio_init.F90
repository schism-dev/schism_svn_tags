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


    SUBROUTINE bio_init(trel0)

        
!!======================================================================
!! April/May, 2007 - Original code                                     ! 
!!======================================================Marta Rodrigues=
!!                                                                     !
!! This subroutine is from ROMS (where is called ana_biology.h):       !
!!                                                                     !
!! February, 2009 - Extended for oxygen cycle                          !
!!======================================================================	 
!!
!!======================================================================
!! Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!!   Licensed under a MIT/X style license                              !
!!   See License_ROMS.txt                                              !
!!                                                                     !
!=======================================================================
!                                                                      !
!  This routine sets initial conditions for biological tracer fields   !
!  using analytical expressions.                                       !
!                                                                      !
!=======================================================================
!
      USE bio_param
      USE schism_glbl, only : ntracers,nvrt,nea,tsel,idry_e
      USE biology 

      IMPLICIT NONE

!  Imported variable declarations.
      real(r8), dimension(ntracers,nvrt,nea), intent(out) :: trel0
!
!  Local variable declarations.
!
      integer :: i, iseco, itrc, j, k

      real(r8) :: cff1, cff2, cff3, cff4, cff5, cff6, cff7, cff8, cff9
      real(r8) :: cff10, cff11, cff12, cff13, cff14, cff15
      real(r8) :: salt, sftm, temp

!
!---------------------------------------------------------------------
!  EcoSim initial fields.
!---------------------------------------------------------------------
!
! Assumed maximum temperature gradient.
!
      cff3=1.0_r8/14.0_r8
      cff4=1.0_r8/16.0_r8
      cff5=32.0_r8
      cff7=1.0_r8/0.0157_r8
      cff8=1.0_r8/6.625_r8
      cff9=1.0_r8/16.0_r8
      cff10=1.0_r8/15.0_r8
      cff11=1.0_r8/8.0_r8
      cff12=1.0_r8/128.0_r8
      cff13=1.0_r8/1000.0_r8
      cff14=1.0_r8/12.0_r8
      cff15=cff5*cff8*cff14                  ! mole N : gram Chl

!Marta Rodrigues - For vertical and node      
      DO k= nvrt, 1, -1
        DO i=1, nea
          !if(idry_e(i)==1) cycle

! Initialization of surface chlorophyll.
!
!Marta Rodrigues - temperature and salinity values
            sftm=tsel(1,nvrt,i)	!sea surface temperature
            temp=tsel(1,k,i)
            salt=tsel(2,k,i)

	    	    
            cff1=-0.0827_r8*sftm+2.6386_r8
            cff2=MAX(0.00001_r8,cff1*(1.0_r8-(sftm-temp)*cff3))

!
! Initialization of nutrients.
!
            trel0(iNH4_,k,i)=0.053_r8*temp+0.7990_r8
            trel0(iNO3_,k,i)=8.5_r8-cff2*cff15-trel0(iNH4_,k,i)
            trel0(iPO4_,k,i)=(trel0(iNH4_,k,i)+trel0(iNO3_,k,i))*cff4
            IF(IRON==1) trel0(iFeO_,k,i)=1.0_r8  !MFR
!
! Assuming diatoms are 75% of initialized chlorophyll.
!
            trel0(iSiO_,k,i)=5.5_r8-(cff2*0.75_r8)*cff15*1.20_r8
            trel0(iDIC_,k,i)=2000.0_r8

! Marta Rodrigues - Initialization of DO and COD (oxygen)

            trel0(iDO_,k,i)=150.0_r8
            trel0(iCOD_,k,i)=0.0_r8

!
! Bacteria Initialization.
!    
              DO iseco=1,Nbac
                trel0(iBacC(iseco),k,i)=0.85_r8
                trel0(iBacN(iseco),k,i)=trel0(iBacC(iseco),k,i)*N2cBAC
                trel0(iBacP(iseco),k,i)=trel0(iBacC(iseco),k,i)*P2cBAC
                IF(IRON==1) trel0(iBacF(iseco),k,i)=trel0(iBacC(iseco),k,i) &
                           *Fe2cBAC  !MFR
              END DO
!

! Marta Rodrigues
! Zooplankton initialization.

            DO iseco=1,Nzoo
	      trel0(iZooC(iseco),k,i)=2_r8
	      trel0(iZooN(iseco),k,i)=0.2_r8
	      trel0(iZooP(iseco),k,i)=0.02_r8
	    END DO  

! Initialize phytoplankton populations.
!
            trel0(iPhyC(1),k,i)=MAX(0.02_r8,                            &
     &                              0.75_r8*0.75_r8*cff5*cff2*cff14)
            DO iseco=1,Nphy
              trel0(iPhyN(iseco),k,i)=trel0(iPhyC(iseco),k,i)*cff8
              trel0(iPhyP(iseco),k,i)=trel0(iPhyN(iseco),k,i)*cff4
              IF(IRON==1) trel0(iPhyF(iseco),k,i)=trel0(iPhyC(iseco),k,i)&
                          *cff13  !MFR
              IF (iPhyS(iseco).gt.0) THEN
                trel0(iPhyS(iseco),k,i)=trel0(iPhyN(iseco),k,i)*1.20_r8
              END IF
!  Initialize Pigments in ugrams/liter (not umole/liter).
!
              cff6=12.0_r8/cff5
              trel0(iPigs(iseco,1),k,i)=cff6*trel0(iPhyC(iseco),k,i)
!
!  Chlorophyll-b.
!
              cff6=cff5-b_C2Cl(iseco)
              IF (iPigs(iseco,2).gt.0) THEN
                 trel0(iPigs(iseco,2),k,i)=trel0(iPigs(iseco,1),k,i)*         &
     &                                  (b_ChlB(iseco)+                 &
     &                                   mxChlB(iseco)*cff6)
              END IF
!
!  Chlorophyll-c.
!
              IF (iPigs(iseco,3).gt.0) THEN
                 trel0(iPigs(iseco,3),k,i)=trel0(iPigs(iseco,1),k,i)*         &
     &                                  (b_ChlC(iseco)+                 &
     &                                   mxChlC(iseco)*cff6)
              END IF
!
!  Photosynthetic Carotenoids.
!
              IF (iPigs(iseco,4).gt.0) THEN
                 trel0(iPigs(iseco,4),k,i)=trel0(iPigs(iseco,1),k,i)*         &
     &                                  (b_PSC(iseco)+                  &
     &                                   mxPSC(iseco)*cff6)
              END IF
!
!  Photoprotective Carotenoids.
!
              IF (iPigs(iseco,5).gt.0) THEN
                 trel0(iPigs(iseco,5),k,i)=trel0(iPigs(iseco,1),k,i)*         &
     &                                  (b_PPC(iseco)+                  &
     &                                   mxPPC(iseco)*cff6)
              END IF
!
!  Low Urobilin Phycoeurythin Carotenoids.
!
              IF (iPigs(iseco,6).gt.0) THEN
                 trel0(iPigs(iseco,6),k,i)=trel0(iPigs(iseco,1),k,i)*         &
     &                                  (b_LPUb(iseco)+                 &
     &                                   mxLPUb(iseco)*cff6)
              END IF
!
!  High Urobilin Phycoeurythin Carotenoids.
!
              IF (iPigs(iseco,7).gt.0) THEN
                 trel0(iPigs(iseco,7),k,i)=trel0(iPigs(iseco,1),k,i)*         &
     &                                  (b_HPUb(iseco)+                 &
     &                                   mxHPUb(iseco)*cff6)
              END IF
            END DO
!
! DOC initialization.
!
            cff6=MAX(0.001_r8,-0.9833_r8*salt+33.411_r8)
            trel0(iDOMC(1),k,i)=0.1_r8
            trel0(iDOMN(1),k,i)=trel0(iDOMC(1),k,i)*cff8
            trel0(iDOMP(1),k,i)=trel0(iDOMN(1),k,i)*cff9
            IF(CDOC==1) trel0(iCDMC(1),k,i)=trel0(iDOMC(1),k,i)*cDOCfrac_c(1) !MFR
            IF(Ndom==2)THEN  !MFR
              trel0(iDOMC(2),k,i)=15.254_r8*cff6+70.0_r8
              trel0(iDOMN(2),k,i)=trel0(iDOMC(2),k,i)*cff10
              trel0(iDOMP(2),k,i)=0.0_r8
              IF(CDOC==1) trel0(iCDMC(2),k,i)=(0.243_r8*cff6+0.055_r8)*cff7   !MFR
            ENDIF 
!
! Fecal Initialization.
!
            DO iseco=1,Nfec
              trel0(iFecC(iseco),k,i)=0.002_r8
              trel0(iFecN(iseco),k,i)=trel0(iFecC(iseco),k,i)*cff11
              trel0(iFecP(iseco),k,i)=trel0(iFecC(iseco),k,i)*cff12
              IF(IRON==1) trel0(iFecF(iseco),k,i)=trel0(iFecC(iseco),k,i)*cff13 !MFR
              trel0(iFecS(iseco),k,i)=trel0(iFecC(iseco),k,i)*cff11
            END DO
        END DO
      END DO          
 
      RETURN
      END SUBROUTINE bio_init


      
      
      
