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
! MORSELFE INITIALIZATION SUBROUTINES
!
! subroutine sed_alloc
! subroutine sed_init
!
!=====================================================================
!=====================================================================

      SUBROUTINE sed_alloc()
!--------------------------------------------------------------------!
! This subroutine allocates and pre-initialize sediment model arrays !
!                                                                    !
! Adapted from former subroutines initialize_scalars and             !
! initialize_ocean (former init_sed.F90), other allocation formerly  !
! done within schism_init.F90 were also moved within this subroutine. !
!                                                                    !
! Author: florian ganthy (fganthy@lnec.pt ; florian.ganthy@gmail.com)!
! Date: 2013/06/07                                                   !
!                                                                    !
! History:                                                           !
!            ***  Previous history (former subroutines) ***          !
!                                                                    !
!   2012/12 - F.Ganthy : form homogenisation of sediments  routines  !
!   2012/12 - F.Ganthy : modifications for Bottom Composition        !
!                        Generation (BCG) purpose (added output      !
!                        files - bed median grain size, bed fraction)!
!   2013/01 - F.Ganthy : Implementation of roughness predictor       !
!   2013/03 - F.Ganthy : Implementation of wave-induced bedload      !
!                        transport                                   !
!   2013/04 - F.Ganthy : Implementation of wave-current bottom stress!
!   2013/05 - F.Ganthy : Updates related to ripple predictor         !
!   2013/05 - F.Ganthy : Added node-centered  volume control area    !
!   2013/05 - F.Ganthy : Updates to the ripple predictor:            !
!                         - Changes on the total bedform             !
!                           roughness computation                    !
!                         - Add wave-ripple computation from         !
!                           Nielsen (1992)                           !
!                                                                    !
!            ***  Current history (former subroutines) ***           !
!                                                                    !
!                                                                    !
!--------------------------------------------------------------------!

      USE sed_mod

      USE schism_glbl, ONLY: nea,np,npa,mnei,ntracers,rkind,dav,dave,  &
     &                     nvrt
      USE schism_msgp, ONLY: myrank,parallel_abort

      IMPLICIT NONE

!- Local variables --------------------------------------------------!

      INTEGER :: i

      REAL(rkind), PARAMETER :: IniVal = 0.0d0

!- Start Statement --------------------------------------------------!

      IF(myrank==0) write(16,*)'Entering sed_alloc'

!--------------------------------------------------------------------!
!* ARRAYS ALLOCATION
!--------------------------------------------------------------------!

      !--------------------------------------------------------------!
      !* 1D arrays defined on elements
      !--------------------------------------------------------------!
      ALLOCATE(Zob(nea),stat=i)
        IF(i/=0) CALL parallel_abort('Sed: Zob allocation failure')
      ALLOCATE(dave(nea),stat=i)
        IF(i/=0) CALL parallel_abort('Sed: dave allocation failure')
      ALLOCATE(FX_r(nea), stat=i )
        IF(i/=0) CALL parallel_abort('Sed: FX_r allocation failure')
      ALLOCATE(FY_r(nea), stat=i )
        IF(i/=0) CALL parallel_abort('Sed: FY_r allocation failure')
      ALLOCATE(bustr(nea), stat=i )
        IF(i/=0) CALL parallel_abort('Sed: bustr allocation failure')
      ALLOCATE(bvstr(nea), stat=i )
        IF(i/=0) CALL parallel_abort('Sed: bvstr allocation failure')
      ALLOCATE(bed_thick(nea),stat=i)
        IF(i/=0) CALL parallel_abort('Sed: bed_thick allocation failure')
      ALLOCATE(hs(nea),stat=i)
        IF(i/=0) CALL parallel_abort('Sed: hs allocation failure')
      ALLOCATE(tp(nea),stat=i)
        IF(i/=0) CALL parallel_abort('Sed: tp allocation failure')
      ALLOCATE(wlpeak(nea),stat=i)
        IF(i/=0) CALL parallel_abort('Sed: wlpeak allocation failure')
      ALLOCATE(uorb(nea),stat=i)
        IF(i/=0) CALL parallel_abort('Sed: uorb allocation failure')
      ALLOCATE(uorbp(nea),stat=i)
        IF(i/=0) CALL parallel_abort('Sed: uorbp allocation failure')
      ALLOCATE(tau_c(nea),stat=i)
        IF(i/=0) CALL parallel_abort('Sed: tau_c allocation failure')
      ALLOCATE(tau_w(nea),stat=i)
        IF(i/=0) CALL parallel_abort('Sed: tau_w allocation failure')
      ALLOCATE(tau_wc(nea),stat=i)
        IF(i/=0) CALL parallel_abort('Sed: tau_wc allocation failure')

      !--------------------------------------------------------------!
      !* 2D arrays defined on elements
      !--------------------------------------------------------------!
      ALLOCATE(mcoefd(0:(mnei+1),np),stat=i)
        IF (i/=0) CALL parallel_abort('Sed: mcoefd allocation failure')
      ALLOCATE(Hz(nvrt,nea),stat=i)
        IF(i/=0) CALL parallel_abort('Sed: Hz allocation failure')
      ALLOCATE(bottom(nea,MBOTP),stat=i)
        IF(i/=0) CALL parallel_abort('Sed: bottom allocation failure')

      !--------------------------------------------------------------!
      !* 3D arrays defined on elements
      !--------------------------------------------------------------!
      ALLOCATE(bed(Nbed,nea,MBEDP),stat=i)
        IF(i/=0) CALL parallel_abort('Sed: bed allocation failure')
      ALLOCATE(bed_frac(Nbed,nea,ntracers),stat=i)
        IF(i/=0) CALL parallel_abort('Sed: bed_frac allocation failure')

      !--------------------------------------------------------------!
      !* 4D arrays defined on elements
      !--------------------------------------------------------------!
      ALLOCATE(bed_mass(Nbed,nea,2,ntracers),stat=i)
        IF(i/=0) CALL parallel_abort('Sed: bed_mass allocation failure')

      !--------------------------------------------------------------!
      !* 1D arrays defined on nodes
      !--------------------------------------------------------------!
      ALLOCATE(vc_area(npa),stat=i)
        IF(i/=0) CALL parallel_abort('sed_init: vc_area allocation failure')
      ALLOCATE(bedthick_overall(npa),bed_d50n(npa),stat=i)
        IF(i/=0) CALL parallel_abort('sed_init: bed_d50n allocation failure')
      ALLOCATE(bed_taun(npa),stat=i)
        IF(i/=0) CALL parallel_abort('sed_init: bed_taun allocation failure')
      ALLOCATE(bed_rough(npa),stat=i)
        IF(i/=0) CALL parallel_abort('sed_init: bed_rough allocation failure')
      ALLOCATE(lbc_sed(npa),stat=i)
        IF(i/=0) CALL parallel_abort('main: lbc_sed allocation failure')
      ALLOCATE(bc_sed(npa),stat=i)
        IF(i/=0) CALL parallel_abort('main: bc_sed allocation failure')

      !--------------------------------------------------------------!
      !* 2D arrays defined on nodes
      !--------------------------------------------------------------!
      ALLOCATE(bedldu(npa,ntracers),stat=i)
        IF(i/=0) CALL parallel_abort('Sed: bedldu allocation failure')
      ALLOCATE(bedldv(npa,ntracers),stat=i)
        IF(i/=0) CALL parallel_abort('Sed: bedldv allocation failure')
      ALLOCATE(bed_fracn(npa,ntracers),stat=i)
        IF(i/=0) CALL parallel_abort('Sed: bed_fracn allocation failure')

      !--------------------------------------------------------------!
      !* Other 1D arrays
      !--------------------------------------------------------------!
      ALLOCATE ( isand(MAX(1,ntracers)), stat=i )
        IF(i/=0) CALL parallel_abort('Main: allocation failure')

!--------------------------------------------------------------------!
!* Pre-initialization
!--------------------------------------------------------------------!

      !--------------------------------------------------------------!
      !* 1D arrays defined on elements
      !--------------------------------------------------------------!
      Zob(:)       = IniVal
      dave(:)      = IniVal
      FX_r(:)      = IniVal
      FY_r(:)      = IniVal
      bustr(:)     = IniVal
      bvstr(:)     = IniVal
      bed_thick(:) = IniVal
      hs(:)        = IniVal
      tp(:)        = IniVal
      wlpeak(:)    = IniVal
      uorb(:)      = IniVal
      uorbp(:)     = IniVal
      tau_c(:)     = IniVal
      tau_w(:)     = IniVal
      tau_wc(:)    = IniVal

      !--------------------------------------------------------------!
      !* 2D arrays defined on elements
      !--------------------------------------------------------------!
      mcoefd(:,:) = IniVal
      Hz(:,:)     = IniVal
      bottom(:,:) = IniVal

      !--------------------------------------------------------------!
      !* 3D arrays defined on elements
      !--------------------------------------------------------------!
      bed(:,:,:)      = IniVal
      bed_frac(:,:,:) = IniVal

      !--------------------------------------------------------------!
      !* 3D arrays defined on elements
      !--------------------------------------------------------------!
      bed_mass(:,:,:,:) = IniVal

      !--------------------------------------------------------------!
      !* 1D arrays defined on nodes
      !--------------------------------------------------------------!
      vc_area(:)   = IniVal
      bed_d50n(:)  = IniVal
      bed_taun(:)  = IniVal
      bed_rough(:) = IniVal

      !Special initializations
      bc_sed(:)    = -9999
      lbc_sed(:)   = .FALSE.

      !--------------------------------------------------------------!
      !* 2D arrays defined on nodes
      !--------------------------------------------------------------!
      bedldu(:,:)    = IniVal
      bedldv(:,:)    = IniVal
      bed_fracn(:,:) = IniVal

      !--------------------------------------------------------------!
      !* Other 1D arrays
      !--------------------------------------------------------------!
      isand(:) = 1 ! Integer



      IF(myrank==0) write(16,*)'Leaving sed_alloc'
!--------------------------------------------------------------------!
      END SUBROUTINE sed_alloc      
      
      
!=====================================================================
!=====================================================================

      SUBROUTINE sed_init()
!--------------------------------------------------------------------!
! This subroutine initialize variables and arrays for the sediment   ! 
! model:                                                             !
!        - computes coefficients for the jgc solver                  !
!        - computes control volume area at nodes                     !
!        - read bed fraction file (bed_frac_x.ic)                    !
!        - mapping of bed fraction                                   !
!        - checking sum of bed fraction                              !
!        - initialize the bed model (layer thickness, age, porosity) !
!        - initialize the bed mass                                   !
!        - initialize exposed sediment layer properties              !
!        - initialize sediment roughness length                      !
!        - initialize total bed thickness                            !
!        - initialization of arrays defined at nodes                 !
!        - if debuging: write sediment initialization                !
!                                                                    !
! Adapted from former subroutines initialize_scalars,                !
! initialize_ocean (former init_sed.F90), and sed_init (former       !
! sed_init.F90) and other initializations formerly done within       !
!schism_init.F90 were also moved within this subroutine.              !
!                                                                    !
! The former subroutine sed_init.F90 was adapted from ROMS routine   !
! ana_sediment.h                                                     !
! Copyright (c) 2002-2007 The ROMS/TOMS Group                        !
!   Licensed under a MIT/X style license                             !
!   See License_ROMS.txt                                             !
!                                                                    !
! Author: florian ganthy (fganthy@lnec.pt ; florian.ganthy@gmail.com)!
! Date: 2013/06/07                                                   !
!                                                                    !
! History:                                                           !
!            ***  Previous history (former subroutines) ***          !
!                                                                    !
!   2012/12 - F.Ganthy : form homogenisation of sediments  routines  !
!   2012/12 - F.Ganthy : modifications for Bottom Composition        !
!                        Generation (BCG) purpose (added output      !
!                        files - bed median grain size, bed fraction)!
!   2013/01 - F.Ganthy : Implementation of roughness predictor       !
!   2013/03 - F.Ganthy : Implementation of wave-induced bedload      !
!                        transport                                   !
!   2013/04 - F.Ganthy : Implementation of wave-current bottom stress!
!   2013/05 - F.Ganthy : Updates related to ripple predictor         !
!   2013/05 - F.Ganthy : Added node-centered  volume control area    !
!   2013/05 - F.Ganthy : Updates to the ripple predictor:            !
!                         - Changes on the total bedform             !
!                           roughness computation                    !
!                         - Add wave-ripple computation from         !
!                           Nielsen (1992)                           !
!                                                                    !
!            ***  Current history (former subroutines) ***           !
!   2013/06 - F.Ganthy : Also removed old initializations currently  !
!                        which can be found in former svn revision   !
!
!                                                                    !
!                                                                    !
!--------------------------------------------------------------------!      

      USE sed_mod

      USE schism_glbl, ONLY: nea,npa,mnei,ntracers,ipgl,ielg,elnode,np_global,  &
     &                     ifile_char,ifile_len,area,np,nne,indel,     &
     &                     isbnd,rough_p,errmsg,ihot
      USE schism_msgp, ONLY: myrank,parallel_abort,exchange_p2d

      IMPLICIT NONE

!- Local variables --------------------------------------------------!

      INTEGER :: i,j,k,ie
      INTEGER :: ised,ic,itmp,istat

      REAL(rkind) :: aux1,aux2,xtmp,ytmp,tmp1,cff1,cff2,cff3,cff4,cff5
      REAL(rkind) :: bed_frac_sum

      REAL(rkind),DIMENSION(npa) :: bdfc

      !swild98 used for exchange (deallocate immediately afterwards)
      REAL(rkind),ALLOCATABLE :: swild98(:,:,:)

      CHARACTER(len=48) :: inputfile

!- Start Statement --------------------------------------------------!

      IF(myrank==0) WRITE(16,*)'Entering sed_init'

!--------------------------------------------------------------------!
! - Set tracer identification indices.
! In SELFE T and S are in diferent arrays from the tracers array...
!--------------------------------------------------------------------!
      ic = 0
      DO i = 1,ntracers
        ic = ic+1
        isand(i) = ic
      END DO

!--------------------------------------------------------------------!
! - Computes matrix coefficients for the JCG solver
! Used for the computation of depth variation induced by bedload
!--------------------------------------------------------------------!
      mcoefd = 0
      aux1 = 22.0d0/108.0d0
      aux2 = 7.0d0/108.0d0
      DO i = 1,np !residents
        DO j = 1,nne(i)
          ie = indel(j,i)
          mcoefd(0,i) = mcoefd(0,i)+area(ie)*aux1 !diagonal
          mcoefd(j,i) = mcoefd(j,i)+area(ie)*aux2
          IF((isbnd(1,i).EQ.0).AND.(j.EQ.nne(i))) THEN !internal ball
             mcoefd(1,i) = mcoefd(1,i)+area(ie)*aux2
          ELSE
             mcoefd(j+1,i) = mcoefd(j+1,i)+area(ie)*aux2
          ENDIF
        ENDDO ! END loop nne
      ENDDO ! END loop np

!--------------------------------------------------------------------!
! - Control volume at each node
!--------------------------------------------------------------------!
      vc_area = 0.0d0
      DO i=1,nea
        DO j=1,3
          vc_area(elnode(j,i)) = vc_area(elnode(j,i)) + (1.0d0/3.0d0)*area(i)
        ENDDO ! End loop j
      ENDDO ! End loop nea
      CALL exchange_p2d(vc_area)

!     Read in total bed thickness at nodes
      OPEN(10,FILE='bedthick.ic',STATUS='OLD')
      READ(10,*); READ(10,*)
      DO i = 1,np_global
        READ(10,*) itmp,xtmp,ytmp,tmp1
        IF(tmp1<=0) CALL parallel_abort('Sed_init: bed_thick<=0!')
        IF(ipgl(i)%rank==myrank) bedthick_overall(ipgl(i)%id)=tmp1
      ENDDO !i=1,np_global
      CLOSE(10)

!     For cold start only
      if(ihot==0) then
!========================================================
!--------------------------------------------------------------------!
! - Reading bed_frac files and convert bed_fraction from nodes to 
! elements.
! For instance, bed_fraction is applied to all bed layers, i.e. no vertical variation
!--------------------------------------------------------------------!
        ALLOCATE(swild98(Nbed,npa,ntracers),stat=istat)
        ! * READING bed_frac_x.ic
!        DO k = 1,Nbed
        DO ised = 1,ntracers
          WRITE(ifile_char,'(i03)')ised
          ifile_char=ADJUSTL(ifile_char); ifile_len=LEN_TRIM(ifile_char)
          inputfile='bed_frac_'//ifile_char(1:ifile_len)//'.ic'
          OPEN(10,FILE=inputfile,STATUS='OLD')
          READ(10,*) !read in first line, no need to store it
          READ(10,*) !read in second line, no need to store it
          DO i = 1,np_global
            READ(10,*) itmp,xtmp,ytmp,tmp1
            IF(tmp1.LT.0) CALL parallel_abort('Sed: bed_frac <0!')
            IF(tmp1.GT.1) CALL parallel_abort('Sed: bed_frac >1!')
            IF(ipgl(i)%rank==myrank) swild98(:,ipgl(i)%id,ised) = tmp1
          ENDDO !i=1,np_global
          CLOSE(10)
        ENDDO !ised=1,ntracers
!        ENDDO !k=1,Nbed

!--------------------------------------------------------------------!
! - Mapping bed fraction (conversion from node to elements)
!--------------------------------------------------------------------!
        DO ised = 1,ntracers
          DO k = 1,Nbed
            DO i = 1,nea
              bed_frac(k,i,ised) = (swild98(k,elnode(1,i),ised)+           &
              &                     swild98(k,elnode(2,i),ised)+           &
              &                     swild98(k,elnode(3,i),ised))/3.0d0
            ENDDO ! END loop nea
          ENDDO ! END loop Nbed
        ENDDO ! END loop ntracers
        DEALLOCATE(swild98,stat=istat)

!--------------------------------------------------------------------!
! - Check the sum of bed fractions 
!   Sum should be equal to 1 but not exactly possible due to rounding
!       - Sum must not exceed 1.01
!       - Sum must not be lower than 0.99
!--------------------------------------------------------------------!
        DO k=1,Nbed
          DO i=1,nea
            bed_frac_sum = 0.0d0
            DO ised=1,ntracers
              bed_frac_sum = bed_frac_sum+bed_frac(k,i,ised)
            ENDDO !End loop ntracers
            IF (bed_frac_sum.GT.1.01d0) THEN
              WRITE(errmsg,*)'SED: sum of bed_frac >1 at elem ',i,     &
              &              bed_frac_sum
              CALL parallel_abort(errmsg)
            ENDIF
            IF (bed_frac_sum.LT.0.99d0) THEN
              WRITE(errmsg,*)'SED: sum of bed_frac <1 at elem ',i,     &
              &              bed_frac_sum
              CALL parallel_abort(errmsg)
            ENDIF
          ENDDO ! End loop nea
        ENDDO ! End loop Nbed

!--------------------------------------------------------------------!
! - Initialize the bed model (layer thickness, age and porosity)
!--------------------------------------------------------------------!
        DO i=1,Nbed
          DO j=1,nea
            bed(i,j,ithck) = sum(bedthick_overall(elnode(1:3,j)))/3.d0/real(Nbed) !>0
            bed(i,j,iaged) = 0.0d0
            bed(i,j,iporo) = porosity
          ENDDO ! End loop Nbed
        ENDDO !End loop nea
!========================================================
      endif !ihot==0

!--------------------------------------------------------------------!
! - Initialize the bed mass ([kg/m2]=[m]*[kg.m-3]*[-])
!--------------------------------------------------------------------!
      DO k=1,Nbed
        DO j=1,nea
          DO ised=1,ntracers
             bed_mass(k,j,:,ised) = bed(k,j,ithck)*                  &
             &                      Srho(ised)*                      &
             &                      (1.0d0-bed(k,j,iporo))*          &
             &                      bed_frac(k,j,ised)
          ENDDO ! End loop ntracers
        ENDDO ! End loop nea
      END DO ! End loop Nbed

!--------------------------------------------------------------------!
! - Initialize exposed sediment layer properties:
!      - total D50 (m)
!      - total density (kg.m-3)
!      - total settling velocity (m.s-1)
!      - total critical shear stress (m2.s-2)
!--------------------------------------------------------------------!
      DO j=1,nea
        cff1 = 1.0d0
        cff2 = 1.0d0
        cff3 = 1.0d0
        cff4 = 1.0d0
        cff5 = 0.0d0
        ! Weighted geometric mean
        DO ised=1,ntracers
          cff1 = cff1*Sd50(ised)**bed_frac(1,j,ised)
          cff2 = cff2*Srho(ised)**bed_frac(1,j,ised)
          cff3 = cff3*Wsed(ised)**bed_frac(1,j,ised)
          cff4 = cff4*tau_ce(ised)**bed_frac(1,j,ised)
          cff5 = cff5+bed_frac(1,j,ised)
        ENDDO ! End loop ntracers

        if(cff5<0) then
          call parallel_abort('SED_INIT: cff5<0 (1)')
        else if(cff5==0) then !care-takers for all-eroded case
          WRITE(12,*)'SED_INIT: all eroded at elem. ',ielg(j)
          bottom(j,isd50) = Sd50(1)
          bottom(j,idens) = Srho(1)
          bottom(j,iwsed) = Wsed(1)
          bottom(j,itauc) = tau_ce(1)
        else
          bottom(j,isd50) = cff1**(1.0d0/cff5)
          bottom(j,idens) = cff2**(1.0d0/cff5)
          bottom(j,iwsed) = cff3**(1.0d0/cff5)
          bottom(j,itauc) = cff4**(1.0d0/cff5)
        endif !cff5
      ENDDO ! j

!--------------------------------------------------------------------!
! - Initialize sediment main roughness length and bedform predictor
! related roughness length
! Here, the initialization is the same if using or not the bedfrom
! predictor.
!--------------------------------------------------------------------!
      ! Current ripple roughness length (Soulsby, 1997)
      bottom(:,izcr)  = 0.0d0
      ! Sand waves roughness length (Van Rijn, 1984)
      bottom(:,izsw)  = 0.0d0
      ! Wave ripples roughness length (Grant and Madsen, 1982 OR Nielsen, 1992)
      bottom(:,izwr)  = 0.0d0
      ! Bed load bottom roughness length (Grant and Madsen, 1982 OR Nielsen, 1992)
      bottom(:,izbld) = 0.0d0

      DO i=1,nea
        ! Nikurasde roughness length
        bottom(i,izNik) = bottom(i,isd50)/12.0d0
        ! Default roughness
        bottom(i,izdef) = (rough_p(elnode(1,i))+                         &
        &                  rough_p(elnode(2,i))+                         &
        &                  rough_p(elnode(3,i)))/3.0d0
        ! Apparent initial roughness
        bottom(i,izapp) = bottom(i,izdef)
        ! Roughness length effectively used (even if bedform predictor is not used)
        Zob(i) = bottom(i,izdef)
        IF(zob(i).LE.0.0d0)THEN
          WRITE(errmsg,*) 'SED_INIT: Zob <=0 at elem, rank:', Zob(i),i,myrank
          CALL parallel_abort(errmsg)
        ENDIF
      ENDDO ! i

!--------------------------------------------------------------------!
! - Initialize total bed thickness by summing over all layers
!   This array is not actively updated after this
!--------------------------------------------------------------------!
      bed_thick(:)=0.0d0
      DO i=1,Nbed
        DO j=1,nea
          bed_thick(j) = bed_thick(j)+bed(i,j,ithck)
        ENDDO !End loop nea
      ENDDO ! End loop Nbed

!--------------------------------------------------------------------!
! - Initialization of variables defined at nodes
!--------------------------------------------------------------------!
      bed_d50n(:)    = 0.0d0
      bed_taun(:)    = 0.0d0
      bed_fracn(:,:) = 0.0d0
      bdfc(:)        = 0.0d0
      DO i=1,nea
        DO j=1,3
          DO ised=1,ntracers
            bed_fracn(elnode(j,i),ised) = bed_fracn(elnode(j,i),ised)+       &
            &                         bed_frac(1,i,ised)
          ENDDO ! END loop ntracers
          bed_d50n(elnode(j,i))  = bed_d50n(elnode(j,i))+bottom(i,isd50)
          bdfc(elnode(j,i))      = bdfc(elnode(j,i))+1.0d0
        ENDDO !END loop j
      ENDDO !END loop nea
      DO i=1,npa
        if(bdfc(i)==0) call parallel_abort('SED_INIT: bdfc(i)==0')
        DO ised=1,ntracers
          bed_fracn(i,ised) = bed_fracn(i,ised)/bdfc(i)
        ENDDO !END loop ntracers
        bed_d50n(i)    = bed_d50n(i)/bdfc(i)
        bed_rough(i)   = rough_p(i) !roughness length
      ENDDO !END loop npa
      DO ised=1,ntracers
        CALL exchange_p2d(bed_fracn(:,ised))
      ENDDO !END loop ntracers
      CALL exchange_p2d(bed_d50n)
      CALL exchange_p2d(bed_rough)

!--------------------------------------------------------------------!
! - Sediment debug outputs
!--------------------------------------------------------------------!
!      IF(sed_debug==1) THEN
!      IF(myrank==0) THEN
      DO i=1,Nbed
        DO j=1,nea
          DO k=1,ntracers
            WRITE(12,*)'Nbed:',i,' nea:',j,' ntracers:',k,       &
              &                    ' bed_frac.ic:',real(bed_frac(i,j,k))
          ENDDO !End loop ntracers
        ENDDO !End loop nea
      ENDDO !End loop Nbed

      DO i=1,Nbed
        DO j=1,nea
          WRITE(12,*)'Nbed:',i,' nea:',j,' bed_thick:',real(bed_thick(j)),  &
            &     ' bed(ithck):',real(bed(i,j,ithck)),' bed(iaged):',    &
            &     real(bed(i,j,iaged)),' bed(iporo):',real(bed(i,j,iporo))
        ENDDO !End loop nea
      ENDDO !End loop Nbed

      DO i=1,nea
        WRITE(12,*)'nea:',i,' bottom(itauc):',real(bottom(i,itauc)),        &
          &     ' bottom(isd50):',real(bottom(i,isd50)),                 &
          &     ' bottom(iwsed):',real(bottom(i,iwsed)),                 &
          &     ' bottom(idens):',real(bottom(i,idens))
      ENDDO !End loop nea
!      ENDIF !End myrank
!      ENDIF !End sed_debug

!--------------------------------------------------------------------!
      IF(myrank==0) write(16,*)'Leaving sed_init'
!--------------------------------------------------------------------!
      END SUBROUTINE sed_init



