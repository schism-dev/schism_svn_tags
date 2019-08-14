      subroutine sediment(it,moitn,mxitn,rtol,dave)
!--------------------------------------------------------------------!
!  This routine computes the sediment sources and sinks and adds     !
!  then the global sediment tracer fields. Currently, it includes    !
!  the following:                                                    !
!                                                                    !
!  * Vertical settling of sediment in the water column.              !
!  * Erosive and depositional flux interactions of sediment          !
!    between water column and the bed.                               !
!  * Transport of multiple grain sizes.                              !
!  * Bed layer stratigraphy.                                         !
!  * Bed morphology.                                                 !
!  * Bedload based on Meyer Peter Mueller or Van Rijn, 2007, Journal !
!    of Hydraulic Engineering                                        !
!  * Bedload slope term options: Damgaard et al, 1997, Journal       !
!    of Hydraulic Engineerig v123, p 1130-1138; Antunes do Carmo,    !
!    1995 PhD Thesis; Lesser et al, 2004, Coastal Engineering, v 51, !
!    p 883-915.                                                      !
!                                                                    !
!  * Seawater/sediment vertical level distribution:                  !
!                                                                    !
!         W-level  RHO-level                                         !
!                                                                    !
!            N     _________                                         !
!                 |         |                                        !
!                 |    N    |                                        !
!          N-1    |_________|  S                                     !
!                 |         |  E                                     !
!                 |   N-1   |  A                                     !
!            3    |_________|  W                                     !
!                 |         |  A                                     !
!                 |    3    |  T                                     !
!            2    |_________|  E                                     !
!                 |         |  R                                     !
!                 |    2    |                                        !
!            1    |_________|_____ bathymetry                        !
!                 |/////////|                                        !
!                 |    1    |                                        !
!            1    |_________|  S                                     !
!                 |         |  E                                     !
!                 |    2    |  D                                     !
!            2    |_________|  I                                     !
!                 |         |  M                                     !
!                 |  Nbed-1 |  E                                     !
!        Nbed-1   |_________|  N                                     !
!                 |         |  T                                     !
!                 |  Nbed   |                                        !
!         Nbed    |_________|                                        !
!                                                                    !
!                                                                    !
! This subroutine is adapted from ROMS routines                      !
! Copyright (c) 2002-2007 The ROMS/TOMS Group                        !
!   Licensed under a MIT/X style license                             !
!   See License_ROMS.txt                                             !
!                                                                    !
! Author: Ligia Pinto                                                !
! Date: xx/08/2007                                                   !
!                                                                    !
! History: 2012/12 - F.Ganthy : form homogenisation of sediments     !
!          routines                                                  !
!          2012/12 - F.Ganthy : modifications for Bottom Composition !
!                               Generation (BCG) purpose (added      !
!                               output files - bed median grain size,!
!                               bed fraction - and modification of   !
!                               morphodynamics management            !
!          2013/01 - F.Ganthy : Implementation of roughness predictor!
!          2013/01 - F.Ganthy : Implementation of avalanching        !
!          2013/03 - F.Ganthy : Implementation of bedmass filter     !
!          2013/03 - F.Ganthy : Implementation wave-induced bedload  !
!                               transport                            !
!          2013/04 - F.Ganthy : Implementation of wave-current bottom!
!                               stress                               !
!          2013/05 - F.Ganthy : Add different sediment behavior:     !
!                                - MUD-like or SAND-like             !
!          2013/06 - F.Ganthy : Modification for BCG related to the  !
!                               multiple bed layer model             !
!                                                                    !
!--------------------------------------------------------------------!

      USE sed_mod
      USE elfe_glbl, ONLY : rkind,nvrt,nea,npa,np,ntracers,idry_e,   &
     &                      idry,area,xnd,ynd,znl,dt,elnode,xctr,yctr,   &
     &                      kbe,ze,pi,nne,indel,tr_el,bdy_frc,mntr,    &
     &                      errmsg,ielg,iplg,nond_global,iond_global,&
     &                      ipgl,nope_global,np_global,dp,h0,dpe,    &
     &                      iegl,out_wwm,pi,eta2,dp00
      USE elfe_msgp

      IMPLICIT NONE
      SAVE

!- Local variables --------------------------------------------------!

      INTEGER,INTENT(IN)     :: it          ! time step
      INTEGER,INTENT(IN)     :: moitn,mxitn ! JCG solver int and max 
                                            ! iterations
      REAL(rkind),INTENT(IN) :: rtol        ! Relative tolerance 
      REAL(rkind),INTENT(IN) :: dave(nea)   ! Depth-averaged vel. at centroids


!      INTEGER :: nwild(nea+12)
!      INTEGER :: ibnd,isd,isd00,ind1,ind2,ndo,ndf(npa)
      INTEGER :: ie,nd,ip
      INTEGER :: Ksed,i,indx,ised,j,k,ks,l
      INTEGER :: nm1,nm2,nm3 !bnew
      INTEGER, PARAMETER :: top = 1      ! Top layer of bed
      INTEGER, PARAMETER :: mirror = 16  ! FD for mirror.out
      INTEGER, PARAMETER :: debug = 12   ! outputs/nonfatal_xxxx 
!      INTEGER , PARAMETER :: sf_dam = 1  ! Slope form. Damgaard
!      INTEGER , PARAMETER :: sf_del = 2  ! Slope form. Delft 
!      INTEGER , PARAMETER :: sf_car = 3  ! Slope form. Carmo 

      REAL(rkind), PARAMETER :: eps = 1.0d-14

      REAL(rkind) :: hdep(nea),hbed(npa) !depth change (in dt) due to suspended load and bedload
      REAL(rkind) :: hbed_ised(npa),hdep_nd(npa),ta,tmp
      REAL(rkind) :: time
      REAL(rkind) :: cff, cff1, cff2, cff3, cffL, cffR, dltL, dltR
      REAL(rkind) :: cu, ero_flux, cff4, cff6, aref, cff7, cff8, cff9
      REAL(rkind) :: thck_avail, thck_to_add
 
      ! - For suspended sediment
      INTEGER, dimension(nvrt,nea)     :: ksource
      REAL(rkind)                      :: Hz_inv3
      REAL(rkind), DIMENSION(nvrt)     :: Hz_inv
      REAL(rkind), DIMENSION(nvrt)     :: Hz_inv2
      REAL(rkind), DIMENSION(nvrt,nea) :: FC
      REAL(rkind), DIMENSION(nvrt,nea) :: qc
      REAL(rkind), DIMENSION(nvrt,nea) :: qR
      REAL(rkind), DIMENSION(nvrt,nea) :: qL
      REAL(rkind), DIMENSION(nvrt,nea) :: WR
      REAL(rkind), DIMENSION(nvrt,nea) :: WL
      ! - For bed load
      REAL(rkind) :: cff5
      REAL(rkind) :: bedld_mass
      REAL(rkind) :: smgdr, osmgd, Umag
      REAL(rkind) :: tauc0
      REAL(rkind) :: derx1,derx2,derx3,dery1,dery2,dery3
      REAL(rkind) :: yp,xp,flux      
      REAL(rkind), DIMENSION(nea,mntr) :: dep_mass

      ! - For MPM bed load
      REAL(rkind) :: alphas,tauc

      ! - For VR bed load
      REAL(rkind) :: tsta,dpar

      ! - For morphology
      INTEGER     :: kbed
      REAL(rkind) :: qsan(npa) 
      REAL(rkind) :: dhnd(npa) !total depth change in dt (=sus+bedload)

      ! - For waves
      REAL(rkind) :: htot
      REAL(rkind) :: kpeak


!- Start Statement --------------------------------------------------!

!      IF(myrank==0) WRITE(mirror,*)'SED: Entering sediment'

      ! Used to distinguish between old and new bedloads
      ! it-1 if initial time step=1
      ! nstp=1
      nstp = 1+MOD(it-1,2)
      nnew = 3-nstp
      !bnew = nnew

      ! Initialize  variables
      ks       = 0
      ksource  = 0
      Hz       = 0.d0
      qL       = 0.d0
      qR       = 0.d0
      qc       = 0.d0
      cffR     = 0.d0
      cffL     = 0.d0
      FC       = 0.d0
      cff      = 0.d0
      bdy_frc  = 0.d0
      dep_mass = 0.d0
      hdep_nd  = 0.d0
      hdep     = 0.d0
      hbed     = 0.d0
      angleu   = 0.d0
      anglev   = 0.d0
      bedldu   = 0.d0
      bedldv   = 0.d0

      DO i=1,nea
        FX_r(i)  = 0.d0
        FY_r(i)  = 0.d0
        bustr(i) = 0.d0
        bvstr(i) = 0.d0
      ENDDO

      ! Wave parameters initialisation
      hs     = 0.0d0
      tp     = 0.0d0
      wlpeak = 0.0d0
      uorb   = 0.0d0
      uorbp  = 0.0d0

      ! RUNTIME
      time=dt*it


!---------------------------------------------------------------------
! - Get wave parameters (defined at nodes) and converte to element
!   centres
!---------------------------------------------------------------------
#ifdef USE_WWM
      DO i = 1,nea
        IF (idry_e(i).EQ.1) CYCLE
        htot  = 0.0d0
        kpeak = 0.0d0
        DO j = 1,3
          htot      = htot + (dp(elnode(j,i))+eta2(elnode(j,i)))/3.0d0
          hs(i)     = hs(i) + out_wwm(elnode(j,i),1)/3.0d0
          tp(i)     = tp(i) + out_wwm(elnode(j,i),12)/3.0d0
          wlpeak(i) = wlpeak(i) + out_wwm(elnode(j,i),17)/3.0d0 !Peak wave length
          uorb(i)   = uorb(i) + out_wwm(elnode(j,i),22)/3.0d0
        ENDDO ! End loop 3
        ! * FG - Uorbp unused now
        !kpeak    = 2.0d0*pi/wlpeak(i)
        !uorbp(i) = pi*hs(i)/(tp(i)*DSINH(kpeak*htot))
      ENDDO ! End loop nea
#endif

!---------------------------------------------------------------------
! - Compute bottom stress
!   Test of whether WWM is used or not is performed within the subroutine
!---------------------------------------------------------------------
      CALL sed_wavecurrent_stress()

!---------------------------------------------------------------------
! - Sediment debug output
!---------------------------------------------------------------------

      IF(sed_debug==1) THEN
        CALL sed_write_debug(it)  
      ENDIF



!---------------------------------------------------------------------
!           *** Compute bedload sediment transport. ***
!                Adapted from [ROMS sed_bedload.F]
!---------------------------------------------------------------------
      IF(myrank==0) WRITE(16,*)'SED: start bedload...'

      ! Compute some constant bed slope parameters.
      ! Friction angle = 33º
      ! sed_angle in rad van rijn 35�
      sed_angle = DTAN(30.0_r8*pi/180.0_r8)


      ! Compute bedload boundary conditions
      IF (sed_morph.GE.1) THEN

        ! Identify open boundaries and impose flux of 0 for JCG
        bc_sed = -9999 !flags
        DO i=1,nope_global
          DO j=1,nond_global(i)
            nd = iond_global(i,j) !global
            IF(ipgl(nd)%rank==myrank) THEN
              ip = ipgl(nd)%id
              bc_sed(ip) = 0
            ENDIF !ipgl(nd)%rank==myrank
          ENDDO !End loop nond_global
        ENDDO !End loop nope_global

        ! Compute b.c. flag for all nodes for the matrix
        DO i=1,npa
          IF(bc_sed(i)>-9998) THEN
            lbc_sed(i)=.TRUE.
          ELSE
            lbc_sed(i)=.FALSE.
          ENDIF
        ENDDO ! End loop npa

      ENDIF ! sed_morph.GE.1

!---------------------------------------------------------------------
! - Loop over all non-cohesive sediment classes
!---------------------------------------------------------------------
      DO ised=1,ntracers 

        smgd  = (Srho(ised)/rhom-1.0d0)*g*Sd50(ised)
        if(smgd<=0) call parallel_abort('SED3D: smgd<=0')
        osmgd = 1.0d0/smgd
        smgdr = SQRT(smgd)*Sd50(ised)

!---------------------------------------------------------------------
! - Loop over elements
!---------------------------------------------------------------------
        DO i=1,nea

          IF (idry_e(i)==1) CYCLE

!---------------------------------------------------------------------
! - Compute bottom slopes
!   Jan check x-y-assignment
!jl. Derivatives of shape functions
!    An equivalent implementation is found in elfe_main for dl 
!    Dphi/Dx, Dphi/Dy
!---------------------------------------------------------------------
          nm1= elnode(1,i)
          nm2= elnode(2,i)
          nm3= elnode(3,i)
          derx1 = ynd(nm2)-ynd(nm3)
          derx2 = ynd(nm3)-ynd(nm1)
          derx3 = ynd(nm1)-ynd(nm2)
          dery1 = xnd(nm3)-xnd(nm2)
          dery2 = xnd(nm1)-xnd(nm3)
          dery3 = xnd(nm2)-xnd(nm1)

          ! Element area from SELFE
          dzdx = (dp(nm1)*derx1+dp(nm2)*derx2+dp(nm3)*derx3)/        &
          &      2.0d0/area(i)
          dzdy = (dp(nm1)*dery1+dp(nm2)*dery2+dp(nm3)*dery3)/        &
          &      2.0d0/area(i)

          ! Compute vector components of bottom stress induced by curents
          ! Here it is assumed that effects of waves on current direction
          ! are already taken into account
          IF (tau_c(i).EQ.0.d0) THEN
            angleu = 0.d0
            anglev = 0.d0
          ELSE
            angleu = bustr(i)/tau_c(i) !cos angle
            anglev = bvstr(i)/tau_c(i) !sin angle
          ENDIF

!---------------------------------------------------------------------
! - Computation of bedload
!   returns FX_r amd FY_r
!---------------------------------------------------------------------
          
          IF (bedload == 1) THEN

            IF(Sedtype(ised).EQ.0) THEN
              !* MUD-like sediment type, no bedload
              FX_r(i) = 0.0d0
              FY_r(i) = 0.0d0
            ELSEIF(Sedtype(ised).EQ.1) THEN
              !* SAND-like sediment (0.05 < D50 < 2.0 mm) 
              ! van Rijn bedload can be applied
              CALL sed_bedload_vr(ised,i,dave)
            ELSE !IF(Sedtype(ised).EQ.2) THEN
              !* GRAVEL-like sediment (>= 2.0 mm)
              ! Need a specific bedload transport subroutine
              ! NOT AVAILABLE YET
              WRITE(errmsg,*)'SED: GRAVEL type not available yet'
              CALL parallel_abort(errmsg)
            ENDIF

          ELSEIF (bedload == 2) THEN
            ! bedload mpm UNFINISHED
            !CALL sed_bedload_mpm()
            WRITE(errmsg,*)'SED: MPM bedload not ready yet'
            CALL parallel_abort(errmsg)
          ENDIF

!---------------------------------------------------------------------
! - Apply bedload transport rate coefficient (nondimensional).
! jl. Added bedload transport as implemented in ROMS 
!---------------------------------------------------------------------

          FX_r(i) = FX_r(i)*bedload_coeff
          FY_r(i) = FY_r(i)*bedload_coeff

!---------------------------------------------------------------------
! - Diffusive fluxes in flow direction (Fortunato et al., 2009; eq 2)
!---------------------------------------------------------------------
          FX_r(i) = FX_r(i)+bdldiffu*(1.0-porosity)*DABS(FX_r(i))*dzdx
          FY_r(i) = FY_r(i)+bdldiffu*(1.0-porosity)*DABS(FY_r(i))*dzdy

!---------------------------------------------------------------------
! - Consistency check
!---------------------------------------------------------------------
          
          IF (FX_r(i)/=FX_r(i)) THEN
            WRITE(errmsg,*)'FX_r is NaN',myrank,i,FX_r(i)
            CALL parallel_abort(errmsg)
          ENDIF
          IF (FY_r(i)/=FY_r(i)) THEN
            WRITE(errmsg,*)'FY_r is NaN',myrank,i,FY_r(i)
            CALL parallel_abort(errmsg)
          ENDIF

!---------------------------------------------------------------------
! - END Loop over elements
!---------------------------------------------------------------------
        ENDDO !End loop i=1,nea


!---------------------------------------------------------------------
! - Compute morphology/bed change (characteristics and dh) due to bed load
!---------------------------------------------------------------------

        IF (sed_morph.GE.1) THEN
          IF(myrank.EQ.0) WRITE(16,*)'SED: Entering bedchange_bedload:',ised,it
          CALL bedchange_bedload(ised,it,moitn,mxitn,rtol,qsan,      &
          &                      hbed,hbed_ised)
          IF(myrank.EQ.0) WRITE(16,*)'SED: leaving bedchange_bedload:',ised,it
        ENDIF

!-----------------------------------------------------------------------
!  Output bedload fluxes.
!-----------------------------------------------------------------------
        DO i=1,np
          IF(idry(i)==1) CYCLE
          ks = 0 !total count

          DO j=1,nne(i)
            k = indel(j,i)
            IF(idry_e(k)==1) CYCLE
            ks = ks+1
            ! compute bedload transport rate in [kg/m/s]
            bedldu(i,ised) = bedldu(i,ised) + FX_r(k)*Srho(ised)/dt
            bedldv(i,ised) = bedldv(i,ised) + FY_r(k)*Srho(ised)/dt
          ENDDO ! End loop nne

          if(ks==0) call parallel_abort('SED3D: (4)')
          bedldu(i,ised) = bedldu(i,ised)/ks
          bedldv(i,ised) = bedldv(i,ised)/ks
        ENDDO !End loop np


        ! Exchange ghosts
        CALL exchange_p2d(bedldu(:,ised))
        CALL exchange_p2d(bedldv(:,ised))

!---------------------------------------------------------------------
! - END Loop over all non-cohesive sediment classes
!---------------------------------------------------------------------
      ENDDO ! ised=1,ntracers

!---------------------------------------------------------------------
! - Update mean surface properties.
! Sd50 must be positive definite, due to BBL routines.
! Srho must be >1000, due to (s-1) in BBL routines.
!---------------------------------------------------------------------

      DO i=1,nea
! * FG. Here, I removed testing for dry element, because bed
!       characteristics have to be updated over the whole domain
!       IF (idry_e(i)==1) CYCLE

        ! update sediment fractions
        cff3 = 0.0d0
        DO ised=1,ntracers
          cff3 = cff3+bed_mass(top,i,nnew,ised)
        ENDDO
        ! Test to prevent div by 0
        cff3=max(cff3,eps)

        DO ised=1,ntracers
          bed_frac(1,i,ised) = bed_mass(top,i,nnew,ised)/cff3 !\in [0,1]
        ENDDO

        ! Weighted geometric mean 
        cff1 = 1.0d0
        cff2 = 1.0d0
        cff3 = 1.0d0
        cff4 = 1.0d0
        cff5 = 0.0d0
        DO ised=1,ntracers
          cff1 = cff1*tau_ce(ised)**bed_frac(top,i,ised)
          cff2 = cff2*Sd50(ised)**bed_frac(top,i,ised)
          cff3 = cff3*(Wsed(ised)+eps)**bed_frac(top,i,ised)
          cff4 = cff4*Srho(ised)**bed_frac(top,i,ised)
          cff5 = cff5+bed_frac(top,i,ised)
        ENDDO !End loop ntracer

        ! Update mean bottom properties
!YJZ Error: cff5 can =0 if all eroded
        if(cff5<0) then
          WRITE(errmsg,*)'SED3D: cff5<0 (2); ',cff5
          call parallel_abort(errmsg)
        else if(cff5==0) then
          WRITE(12,*)'SED3D: all eroded at elem. ',ielg(i),it
          !Care takers
          bottom(i,itauc) = tau_ce(1)
          bottom(i,isd50) = Sd50(1)
          bottom(i,iwsed) = Wsed(1)+eps
          bottom(i,idens) = Srho(1)
        else !cff5>0
          bottom(i,itauc) = cff1**(1.0d0/cff5)
!          bottom(i,isd50) = MIN(cff2,Zob(i))
          bottom(i,isd50) = MAX(MIN(cff2**(1.0d0/cff5),MAXVAL(Sd50(:))), &
     &MINVAL(Sd50(:)))
          bottom(i,iwsed) = cff3**(1.0d0/cff5)
          bottom(i,idens) = MAX(cff4**(1.0d0/cff5),1050.0d0)
        endif !cff5

! * FG. I don't understand why the grain diameter is the minimum
!       between the grain diameter and the roughness length.
!       Roughness length depends on the grain size, but
!       the grain size does not depends on the roughness...
!       However, total D50 should not be lower than the lowest
!       D50(ised) and higher than the highest D50(ised)
! * FG. I added the sum of weights (bed_frac) for the computation
!       of geometric mean, because the rounding of bed fractions
!       leads to total bed fraction slightly different than 1.0,
!       inducing significant errors for bed properties

      ENDDO !i=1,nea


!---------------------------------------------------------------------
!  End computing bedload sediment transport.
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!           *** Compute suspended sediment transport. ***
!
!---------------------------------------------------------------------
      IF (suspended_load == 1) THEN
        IF(myrank.EQ.0) WRITE(16,*)'SED: Entering suspended load...'

!---------------------------------------------------------------------
! - Add sediment Source/Sink terms.
! [ROMS sed_settling.F]
!---------------------------------------------------------------------

!---------------------------------------------------------------------
! - Compute thickness and actual depths in the middle of the volume of 
! control, Hz-layer thickness, z_w==ze - layer depth at RHO
! points and w points.
!---------------------------------------------------------------------

        DO i=1,nea
          IF (idry_e(i)==1) CYCLE

          DO k=kbe(i)+1,nvrt
             Hz(k,i) = ze(k,i)-ze(k-1,i)
             IF(Hz(k,i)<=0) CALL parallel_abort('SEDIMENT: (1)')
          ENDDO !End loop nvrt

        ENDDO !End loop nea


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        SED_LOOP: DO ised=1,ntracers
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!---------------------------------------------------------------------
! - Copy concentration of suspended sediment into scratch array "qc"
! (q-central, restrict it to be positive) which is hereafter
! interpreted as a set of grid-box averaged values for sediment
! concentration.
!---------------------------------------------------------------------

          indx=isand(ised) !=ised

          DO i=1,nea
            IF (idry_e(i)==1) CYCLE
            DO k=kbe(i)+1,nvrt
              qc(k,i) = tr_el(indx,k,i)
            ENDDO !End loop nvrt
          ENDDO !End loop nea

!---------------------------------------------------------------------
! - Vertical sinking of suspended sediment.
! Reconstruct vertical profile of suspended sediment "qc" in terms
! of a set of parabolic segments within each grid box. Then, compute
! semi-Lagrangian flux due to sinking.
!---------------------------------------------------------------------

          IF(myrank==0) WRITE(16,*)'SED: Reconstructing vertical:',ised

          DO i=1,nea
            IF (idry_e(i)==1) CYCLE
            Hz_inv2 = 0.d0

            DO k=kbe(i)+1,nvrt-1
               Hz_inv2(k) = 1.0d0/(Hz(k,i)+Hz(k+1,i))
!               Hz_inv2(k) = 2.0d0/(Hz(k,i)+Hz(k+1,i))
            ENDDO !End loop nvrt-1

            DO k=nvrt-1,kbe(i)+1,-1

              !Error: distance should be 1/2 of Hz_inv2?
              !FC = dqc/dz ? 
              FC(k,i) = (qc(k+1,i)-qc(k,i))*(Hz_inv2(k))

            ENDDO !End loop nvrt backwards
          ENDDO !End loop nea

          DO i=1,nea
            IF (idry_e(i)==1) CYCLE

            cffR = 0.d0
            cffL = 0.d0
            dltR = 0.d0
            dltL = 0.d0
            DO k=kbe(i)+2,nvrt-1
              dltR=Hz(k,i)*FC(k,i)
              dltL=Hz(k,i)*FC(k-1,i)
              cff=Hz(k-1,i)+2.0d0*Hz(k,i)+Hz(k+1,i)
              cffR=cff*FC(k,i)
              cffL=cff*FC(k-1,i)

!---------------------------------------------------------------------
! - Apply PPM monotonicity constraint to prevent oscillations within
! the grid box.
!---------------------------------------------------------------------

              IF (dltR*dltL.LE.0.0d0) THEN
                dltR = 0.0d0
                dltL = 0.0d0

              ELSEIF (ABS(dltR).GT.ABS(cffL)) THEN
                dltR = cffL

              ELSEIF (ABS(dltL).GT.ABS(cffR)) THEN
                dltL=cffR

              ENDIF

!---------------------------------------------------------------------
! - Compute right and left side values (qR,qL) of parabolic segments
! within grid box Hz(k); (WR,WL) are measures of quadratic variations. 
! NOTE: Although each parabolic segment is monotonic within its grid
!       box, monotonicity of the whole profile is not guaranteed,
!       because qL(k+1)-qR(k) may still have different sign than
!       qc(k+1)-qc(k).  This possibility is excluded, after qL and qR
!       are reconciled using WENO procedure.
!---------------------------------------------------------------------

              Hz_inv3 = 0.d0
              Hz_inv3 = 1.0d0/(Hz(k-1,i)+Hz(k,i)+Hz(k+1,i))
              cff = (dltR-dltL)*Hz_inv3
              dltR = dltR-cff*Hz(k+1,i)
              dltL = dltL+cff*Hz(k-1,i)
              qR(k,i) = qc(k,i)+dltR
              qL(k,i) = qc(k,i)-dltL
              WR(k,i) = (2.0d0*dltR-dltL)**2.0d0
              WL(k,i) = (dltR-2.0d08*dltL)**2.0d0

            ENDDO !End loop k=kbe(i)+2,nvrt-1

          ENDDO !End loop nea

          cff=1.0d-14

          DO i=1,nea
            IF (idry_e(i)==1) CYCLE
            DO k=kbe(i)+2,nvrt-2
              dltL = MAX(cff,WL(k,i))
              dltR = MAX(cff,WR(k+1,i))
              qR(k,i) = (dltR*qR(k,i)+dltL*qL(k+1,i))/(dltR+dltL)
              qL(k+1,i) = qR(k,i)
            ENDDO !End kbe(i)+2,nvrt-2
          ENDDO !End loop nea

          DO i=1,nea
            IF(idry_e(i)==1) CYCLE

            ! no-flux boundary condition
            FC(nvrt,i) = 0.0d0

            IF (bc_for_weno == 1) THEN
              ! Linear continuation
              qL(nvrt,i) = qR(nvrt-1,i)
              qR(nvrt,i) = 2.0d0*qc(nvrt,i)-qL(nvrt,i)
           ELSEIF (bc_for_weno == 2) THEN
              ! Neumann
              qL(nvrt,i) = qR(nvrt-1,i)
              qR(nvrt,i) = 1.5d0*qc(nvrt,i)-0.5d0*qL(nvrt,i)
            ELSE
              ! default strictly monotonic
              qR(nvrt,i)   = qc(nvrt,i)
              qL(nvrt,i)   = qc(nvrt,i)
              qR(nvrt-1,i) = qc(nvrt,i)
            ENDIF ! End test bc_for_weno


            IF (bc_for_weno == 1) THEN
              qR(kbe(i)+1,i) = qL(kbe(i)+2,i)
              qL(kbe(i)+1,i) = 2.0d0*qc(kbe(i)+1,i)-qR(kbe(i)+1,i)
            ELSEIF (bc_for_weno == 2) THEN
              qR(kbe(i)+1,i) = qL(kbe(i)+2,i)
              qL(kbe(i)+1,i) = 1.5d0*qc(kbe(i)+1,i)-                 &
              &                0.5d0*qR(kbe(i)+1,i)
            ELSE
              ! bottom grid boxes are re-assumed to be 
              ! piecewise constant.
              qL(kbe(i)+2,i) = qc(kbe(i)+1,i) 
              qR(kbe(i)+1,i) = qc(kbe(i)+1,i) 
              qL(kbe(i)+1,i) = qc(kbe(i)+1,i) 
             ENDIF ! End test bc_for_weno

          ENDDO !End loop nea

!---------------------------------------------------------------------
! - Apply monotonicity constraint again, since the reconciled
! interfacial values may cause a non-monotonic behavior of the 
! parabolic segments inside the grid box.
!---------------------------------------------------------------------

          DO i=1,nea
            IF (idry_e(i)==1) CYCLE

            DO k=kbe(i)+1,nvrt
              dltR = qR(k,i)-qc(k,i)
              dltL = qc(k,i)-qL(k,i)
              cffR = 2.0d0*dltR
              cffL = 2.0d0*dltL

              IF (dltR*dltL.LT.0.0d0) THEN
                dltR = 0.0d0
                dltL = 0.0d0
              ELSEIF (ABS(dltR).GT.ABS(cffL)) THEN
                dltR = cffL
              ELSEIF (ABS(dltL).GT.ABS(cffR)) THEN
                dltL = cffR
              ENDIF

              qR(k,i) = qc(k,i)+dltR
              qL(k,i) = qc(k,i)-dltL
            ENDDO !End loop kbe(i)+1,nvrt

          ENDDO !End loop nea


!---------------------------------------------------------------------
! - After this moment reconstruction is considered complete. The next
! stage is to compute vertical advective fluxes, FC. It is expected
! that sinking may occur relatively fast, the algorithm is designed
! to be free of CFL criterion, which is achieved by allowing
! integration bounds for semi-Lagrangian advective flux to use as
! many grid boxes in upstream direction as necessary.
!
! In the two code segments below, WL is the z-coordinate of the
! departure point for grid box interface z_w (ze) with the same 
! indices FC is the finite volume flux; ksource(k) is index of 
! vertical grid box which contains the departure point 
! (restricted by N(ng)). 
! During the search: also add in content of whole grid boxes
! participating in FC.
!---------------------------------------------------------------------

          cff = dt*ABS(Wsed(ised))

          DO i=1,nea
            IF(idry_e(i)==1) CYCLE

            DO k=kbe(i)+1,nvrt
              FC(k-1,i) = 0.0d0
              !layer depth+sinking(dt*wsed)
              WL(k,i) = ze(k-1,i)+cff
              !thickness*sediment concentration
              WR(k,i) = Hz(k,i)*qc(k,i)
              !index of v-grid box of departure point
              ksource(k,i) = k
            ENDDO !End loop kbe(i)+1,nvrt

            DO k=kbe(i)+1,nvrt
              DO ks=k,nvrt-1
                IF (WL(k,i).GT.ze(ks,i)) THEN
                  ksource(k,i) = ks+1
                  FC(k-1,i) = FC(k-1,i)+WR(ks,i)
                ENDIF
              ENDDO !End loop nvrt-1
            ENDDO !End loop kbe(i)+1,nvrt

          ENDDO !End loop i=1,nea

!---------------------------------------------------------------------
!  Finalize computation of flux: add fractional part.
!---------------------------------------------------------------------

          DO i=1,nea
            IF (idry_e(i)==1) CYCLE

            ! Compute inverse thicknessed to avoid repeated divisions.
            Hz_inv = 0.d0 
            cu = 0.d0
            DO k=kbe(i)+1,nvrt
              Hz_inv(k) = 1.0d0/Hz(k,i)
            ENDDO
            
            DO k=kbe(i)+1,nvrt
              ks = ksource(k,i)

              IF(ks<kbe(i)+1.OR.ks>nvrt) THEN
                WRITE(errmsg,*)'SEDIMENT: out of bound, ',ks,ielg(i),k
                CALL parallel_abort(errmsg)
              ENDIF

              cu = MIN(1.0d0,(WL(k,i)-ze(ks-1,i))*Hz_inv(ks))
              FC(k-1,i) = FC(k-1,i)+                                 &
              &           Hz(ks,i)*cu*                               &
              &          (qL(ks,i)+                                  &
              &           cu*(0.5d0*(qR(ks,i)-qL(ks,i))-             &
              &          (1.5d0-cu)*                                 &
              &          (qR(ks,i)+qL(ks,i)-2.0d0*qc(ks,i))))

            ENDDO !End loop nvrt

            DO k=kbe(i)+1,nvrt
               qc(k,i) = (FC(k,i)-FC(k-1,i))*Hz_inv(k)               
            ENDDO
       
          ENDDO !End loop nea

          IF(myrank==0) WRITE(16,*)'SED: done reconstructing vertical:',ised


!---------------------------------------------------------------------
! - Sediment deposition and resuspension near the bottom.
! Adapted from [ROMS sed_fluxes.F]
!
! The deposition and resuspension of sediment on the bottom "bed"
! is due to precepitation flux FC(kbe(i),:), already computed, and the
! resuspension (erosion, hence called ero_flux). The resuspension is
! applied to the bottom-most grid box value qc(:,kbe(i)+1) so the 
! total mass is conserved. Restrict "ero_flux" so that "bed" cannot go
! negative after both fluxes are applied.
!
! Formulation to transfer the sediment erosive flux between the bottom 
! (reference height) and the center of the bottom computational cell
! transference is made considering that the sediment concentration  
! profile between the reference height and the cell center follows the 
! Rouse concentration profile
!---------------------------------------------------------------------

          cff = 1.0d0/tau_ce(ised)

          ! Ref height above bottom 
          aref = 3.0d0*Sd50(ised)         

          DO i=1,nea
            IF (idry_e(i)==1) CYCLE

            cff6 = SQRT(ABS(bustr(i))+ABS(bvstr(i)))
            IF (cff6.EQ.0.0d0) THEN
              cff7 = 0.0d0
              cff8 = 0.0d0
            ELSE

!jl. The version I replaced it with is Eq. 33 in the paper.
!              cff7=Wsed(ised)/                                       &
!              &    ((1+2*(Wsed(ised)/cff6)**2.d0)*0.4d0*cff6)

              cff7 = Wsed(ised)/(vonKar*cff6)
              cff8 = (aref/Hz(kbe(i)+1,i)/2.d0)**cff7
            ENDIF

!---------------------------------------------------------------------
! - Compute erosion, ero_flux (kg/m2) following Ariathurai and 
! Arulanandan (1978)
!---------------------------------------------------------------------

            cff1 = (1.0d0-bed(top,i,iporo))*bed_frac(top,i,ised)
            cff2 = dt*Erate(ised)*cff1
            cff3 = Srho(ised)*cff1
            cff4 = bed_mass(top,i,nnew,ised)

            ero_flux = MIN(MAX(0.0d0,cff2*(cff*tau_wc(i)-1.0d0))* &
            &cff8,MIN(cff3*bottom(i,iactv),cff4)+FC(kbe(i),i))   


!jl. Any idea of why bed_frac sum would be checked here???
!
!              if (bed_frac(1,i,ised).ne. 1.d0)then
!                 write(mirror,*)time,i,myrank,bed_frac(1,i,ised)
!                 stop
!              endif

!---------------------------------------------------------------------
! - Update sediment concentration
!---------------------------------------------------------------------

            qc(kbe(i)+1,i) = qc(kbe(i)+1,i)+ero_flux/Hz(kbe(i)+1,i)


            IF (sed_morph.GE.1) THEN
!---------------------------------------------------------------------
! Addapted from [ROMS sed_bed.F]
! - Apply morphology factor to flux and settling...
!---------------------------------------------------------------------

              ero_flux = ero_flux*morph_fac(ised)
              FC(kbe(i),i) = FC(kbe(i),i)*morph_fac(ised)

!---------------------------------------------------------------------
! - Depth change due to erosion/deposition of suspended sediment
!---------------------------------------------------------------------
!YJZ Error: when poro~1, hdep would be very large???
              if(bed(top,i,iporo)==1) then
                WRITE(errmsg,*)'SED3D: bed(top,i,iporo)==1; ',top,i,iporo
                CALL parallel_abort(errmsg)
              endif
              hdep(i) = hdep(i)+((ero_flux-FC(kbe(i),i)) / &
              &(Srho(ised)*(1.0d0-bed(top,i,iporo))))


            ENDIF !End test sed_morph


            IF (ero_flux-FC(kbe(i),i).LT.0.0d0) THEN

!---------------------------------------------------------------------
! - If first time step of deposit, then store deposit material in
! temporary array, dep_mass.
!---------------------------------------------------------------------

              IF ((time.GT.(bed(top,i,iaged)+1.1d0*dt)).AND.         &
              &   (bed(top,i,ithck).GT.newlayer_thick)) THEN
                dep_mass(i,ised) = -(ero_flux-FC(kbe(i),i)) !>0  
              ENDIF
              bed(top,i,iaged) = time
            ENDIF


!---------------------------------------------------------------------
! - Update bed mass arrays.
!---------------------------------------------------------------------
!jl.  The whole nnew/bnew,nstp is way more confusing than need be

            bed_mass(top,i,nnew,ised)=MAX(bed_mass(top,i,nnew,ised)-(ero_flux-FC(kbe(i),i)),0.0d0)   
            DO k=2,Nbed
              bed_mass(k,i,nnew,ised)=bed_mass(k,i,nstp,ised)
            ENDDO
            !bed_mass(top,i,nnew,ised)=MAX(bed_mass(top,i,nnew,ised),0.0d0) 
          ENDDO !End loop i=1,nea

          IF (sed_morph.GE.1) THEN
            CALL exchange_e2d(hdep(:))
          ENDIF

!-----------------------------------------------------------------------
! - Update global tracer variables (m Tunits); convert effect of settling to a
!  body force.
!-----------------------------------------------------------------------
! divide by dt 

          DO i=1,nea
            IF (idry_e(i)==1) CYCLE

            DO k=kbe(i)+1,nvrt
              cff = qc(k,i)+tr_el(indx,k,i)
              cff1 = -eps
              IF(cff.LT.cff1) qc(k,i) = -tr_el(indx,k,i)
              bdy_frc(indx,k,i) = qc(k,i)/dt 
            ENDDO
          ENDDO ! End loop nea


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ENDDO SED_LOOP !ised=1,ntracers
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        IF(myrank==0) WRITE(16,*)'SED: out of sus. load loop...'

!---------------------------------------------------------------------
! - If first time step of deposit, create new layer and combine bottom
! two bed layers.
!---------------------------------------------------------------------

        DO i=1,nea
          IF (idry_e(i)==1) CYCLE

          cff = 0.0d0
        
          IF (Nbed.GT.1) THEN
            DO ised=1,ntracers
              cff = cff+dep_mass(i,ised)
            ENDDO

            IF (cff.GT.0.0d0) THEN !deposition ocurred

              ! Combine bottom layers.
              bed(Nbed,i,iporo) = 0.5d0*(bed(Nbed-1,i,iporo)+        &
              &                          bed(Nbed,i,iporo))
              bed(Nbed,i,iaged) = 0.5d0*(bed(Nbed-1,i,iaged)+        &
              &                          bed(Nbed,i,iaged))

              DO ised=1,ntracers
                bed_mass(Nbed,i,nnew,ised) =                         &
                &                   bed_mass(Nbed-1,i,nnew,ised)+    &
                &                   bed_mass(Nbed,i,nnew,ised)
              ENDDO

              ! Push layers down.
              DO k=Nbed-1,2,-1
                bed(k,i,iporo) = bed(k-1,i,iporo)
                bed(k,i,iaged) = bed(k-1,i,iaged)
                DO ised =1,ntracers
                  bed_mass(k,i,nnew,ised) = bed_mass(k-1,i,nnew,ised)
                ENDDO
              ENDDO

              ! Set new top layer parameters.
              DO ised=1,ntracers
                if(dep_mass(i,ised)<0) then
                  WRITE(errmsg,*)'SED3D: dep_mass(i,ised)<0; ',dep_mass(i,ised),i,ised
                  CALL parallel_abort(errmsg)
                endif
                bed_mass(2,i,nnew,ised)=MAX(bed_mass(2,i,nnew,ised)-dep_mass(i,ised),0.0d0) !>=0
                bed_mass(top,i,nnew,ised)=dep_mass(i,ised) !>=0
              ENDDO !ised

            ENDIF ! Deposition occurred (cff>0)
          ENDIF !Nbed>1

!---------------------------------------------------------------------
! If Nbed = 1
! - Recalculate thickness and fractions for all layers.
!   bed(ithck) = bed_mass/Srho
!     [m]      = [m.m-2] / [kg.m-2] --> Do not have to integrate over
!     area
!---------------------------------------------------------------------
          DO k=1,Nbed
            cff3=0.0d0
            DO ised=1,ntracers
              ! cff3 is overall bedmass
              cff3 = cff3+bed_mass(k,i,nnew,ised)
            ENDDO

            !IF (cff3.EQ.0.0d0) cff3=eps
            cff3=max(cff3,eps)
            bed(k,i,ithck) = 0.0d0

            if(bed(k,i,iporo)==1) call parallel_abort('SED3D: div. by 0 (6)')
!' YJZ Error: large velue when bed(k,i,iporo)~1??

            DO ised=1,ntracers
              bed_frac(k,i,ised) = bed_mass(k,i,nnew,ised)/cff3
              bed(k,i,ithck)=MAX(bed(k,i,ithck)+bed_mass(k,i,nnew,ised)/(Srho(ised)*(1.0d0-bed(k,i,iporo))),0.0d0)
            ENDDO !ised
          ENDDO !k=1,Nbed
        ENDDO !i=1,nea

        IF(myrank==0) WRITE(mirror,*)'SED: End of suspended sediment'
      ENDIF !End suspended_load
!---------------------------------------------------------------------
!         *** End of Suspended Sediment only section ***
!---------------------------------------------------------------------

!---------------------------------------------------------------------
! * Application of the bed_mass filter - not quite working so turn off
!---------------------------------------------------------------------
     IF(bedmass_filter.GE.1) THEN
       IF(myrank==0) WRITE(16,*)'SED: doing sed_bedmass_filter '
       CALL sed_bedmass_filter(hdep)
       IF(myrank==0) WRITE(16,*)'SED: done sed_bedmass_filter'
     ENDIF

     IF(myrank==0) WRITE(mirror,*)'SED: adjusting bed layers...'

!---------------------------------------------------------------------
! - Ensure top bed layer thickness is greater or equal than active 
! layer thickness. If need to add sed to top layer, then entrain from 
! lower levels. Create new layers at bottom to maintain Nbed.
!---------------------------------------------------------------------

      DO i=1,nea
        IF (idry_e(i)==1) CYCLE

!---------------------------------------------------------------------
! - Calculate active layer thickness, bottom(i,j,iactv). Based on the 
! relation of Harris and Wiberg (1997)
!---------------------------------------------------------------------

        bottom(i,iactv)=MAX(0.0d0,7.0d-3*(tau_wc(i)-bottom(i,itauc))*rhom)+6.0d0*bottom(i,isd50) !>0


!---------------------------------------------------------------------
! - Apply morphology factor
!jl. The application of morph_fac is arbitrary here.  This is not in  
!    need of immediate attention, but if morph_facs differ between  
!    sediment classes there should be some attention, but if  
!    morph_facs differ between sediment classes there should be some 
!    kind of averaging, or other method, to determine the morph_fac to
!    use.
!---------------------------------------------------------------------

        IF(sed_morph.GE.1) THEN
          bottom(i,iactv)=MAX(bottom(i,iactv)*morph_fac(1),bottom(i,iactv))
        ENDIF

        !IF(bottom(i,iactv).GT.bed(top,i,ithck)) THEN
        IF(bed(top,i,ithck)<bottom(i,iactv)) THEN
          IF(Nbed.EQ.1) THEN
            bottom(i,iactv) = bed(top,i,ithck) !possibly 0
          ELSE ! Nbded>1
            !Increase top bed layer - entrain from layers below
            thck_to_add = bottom(i,iactv)-bed(top,i,ithck) !>0
            thck_avail = 0.0d0
            Ksed=1 !initialize
            DO k=2,Nbed
              IF (thck_avail.LT.thck_to_add) THEN
                thck_avail = thck_avail+bed(k,i,ithck)
                Ksed=k
              ENDIF
            ENDDO !k

            IF(thck_avail<1.e-10) THEN
              write(12,*)'SED3D: not enough sed; likely all eroded:', &
     &ielg(i),thck_avail,thck_to_add,bed(:,i,ithck),it
              bed(:,i,ithck)=0
              bed_frac(:,i,:)=0
              bed_mass(:,i,nnew,:)=0
            ELSE !not bare rock

!---------------------------------------------------------------------
! - Catch here if there was not enough bed material
!---------------------------------------------------------------------
              IF (thck_avail.LT.thck_to_add) THEN
                bottom(i,iactv) = bed(top,i,ithck)+thck_avail
                thck_to_add = thck_avail
              ENDIF

!---------------------------------------------------------------------
! - Update bed mass of top layer and fractional layer
!---------------------------------------------------------------------
              cff2 = MAX(thck_avail-thck_to_add,0.0d0)/MAX(bed(Ksed,i,ithck),eps)
              DO ised=1,ntracers
                cff1=0.0d0
                DO k=1,Ksed
                  cff1 = cff1+bed_mass(k,i,nnew,ised)
                ENDDO
                cff3 = cff2*bed_mass(Ksed,i,nnew,ised) 
                bed_mass(top,i,nnew,ised) = MAX(0.d0,cff1-cff3)
                bed_mass(Ksed,i,nnew,ised) = cff3 !>=0
              ENDDO !ised

!---------------------------------------------------------------------
! - Update thickness of fractional layer ksource_sed
!---------------------------------------------------------------------
              bed(Ksed,i,ithck) = MAX(thck_avail-thck_to_add,0.0d0)

!---------------------------------------------------------------------
! - Upate bed fraction of top layer
!---------------------------------------------------------------------
              cff3=0.0d0
              DO ised=1,ntracers
                cff3 = cff3+bed_mass(top,i,nnew,ised)
              ENDDO

              cff3=max(cff3,eps)

              DO ised=1,ntracers
                bed_frac(top,i,ised) = bed_mass(top,i,nnew,ised)/cff3 !>=0
              ENDDO

!---------------------------------------------------------------------
! Upate bed thickness of top layer
!---------------------------------------------------------------------
              bed(top,i,ithck) = bottom(i,iactv)
            
!---------------------------------------------------------------------
! Pull all layers closer to the surface 
!---------------------------------------------------------------------
              ks = Ksed-2
              DO k=Ksed,Nbed
                bed(k-ks,i,ithck) = bed(k,i,ithck)
                bed(k-ks,i,iporo) = bed(k,i,iporo)
                bed(k-ks,i,iaged) = bed(k,i,iaged)

                DO ised=1,ntracers
                  bed_frac(k-ks,i,ised) = bed_frac(k,i,ised)
                  bed_mass(k-ks,i,nnew,ised) = bed_mass(k,i,nnew,ised)
                ENDDO
              ENDDO !End loop Nbed

!---------------------------------------------------------------------
! By now there are only Nbed-ks layers left
! Split what layer Nbed-ks to make Nbed layers in total and
! conserve total mass. (ks+1 is the number of new layers).
!---------------------------------------------------------------------
              !ks = Ksed-2
              if(ks+1==0) call parallel_abort('SED3D: ks+1==0')
              cff = 1.0d0/REAL(ks+1,r8)
              DO k=Nbed,Nbed-ks,-1
                bed(k,i,ithck) = bed(Nbed-ks,i,ithck)*cff
                bed(k,i,iaged) = bed(Nbed-ks,i,iaged)
                DO ised=1,ntracers
                  bed_frac(k,i,ised) = bed_frac(Nbed-ks,i,ised)
                  bed_mass(k,i,nnew,ised)=bed_mass(Nbed-ks,i,nnew,ised)*cff
                ENDDO
              ENDDO ! k

            ENDIF !thck_avail

          ENDIF  ! End test Nbed > 1
        ENDIF  ! End test increase top bed layer
      ENDDO ! End loop i=1,nea

!---------------------------------------------------------------------
! Update mean surface properties.
! write(*,*)'update mean surface prop' 
! Sd50 must be positive definite, due to BBL routines.
! Srho must be >1000, due to (s-1) in BBL routines
!---------------------------------------------------------------------

! * FG. Here, I removed testing for dry element, because bed
!       characteristics have to be updated over the whole domain
      DO i=1,nea
        cff1 = 1.0d0
        cff2 = 1.0d0
        cff3 = 1.0d0
        cff4 = 1.0d0
        cff5 = 0.0d0
        ! weighted geometric mean
        DO ised=1,ntracers
          cff1 = cff1*tau_ce(ised)**bed_frac(1,i,ised)
          cff2 = cff2*Sd50(ised)**bed_frac(1,i,ised)
          cff3 = cff3*(Wsed(ised)+eps)**bed_frac(1,i,ised)
          cff4 = cff4*Srho(ised)**bed_frac(1,i,ised)
          cff5 = cff5+bed_frac(top,i,ised)
        ENDDO !ntracers

        if(cff5<0) then
          WRITE(errmsg,*)'SED3D: cff5<0 (1); ',cff5
          call parallel_abort(errmsg)
        else if(cff5==0) then
          WRITE(12,*)'SED3D: all eroded at elem. (2):',ielg(i),it
          !Caretakers
          bottom(i,itauc) = tau_ce(1)
          bottom(i,isd50) = Sd50(1)
          bottom(i,iwsed) = Wsed(1)+eps
          bottom(i,idens) = Srho(1)
        else !cff5>0
          bottom(i,itauc) = cff1**(1.0d0/cff5)
!          bottom(i,isd50) = MIN(cff2,Zob(i))
          bottom(i,isd50) = MAX(MIN(cff2**(1.0d0/cff5),MAXVAL(Sd50(:))),MINVAL(Sd50(:)))
          bottom(i,iwsed) = cff3**(1.0d0/cff5)
          bottom(i,idens) = MAX(cff4**(1.0d0/cff5),1050.0d0)
        endif !cff5

! * FG. I don't understand why the grain diameter is the minimum
!       between the grain diameter and the roughness length.
!       Roughness length depends on the grain size, but
!       the grain size does not depends on the roughness...
!       However, total D50 should not be lower than the lowest
!       D50(ised) and higher than the highest D50(ised)
! * FG. I added the sum of weights (bed_frac) for the computation
!       of geometric mean, because the rounding of bed fractions
!       leads to total bed fraction slightly different than 1.0,
!       inducing significant errors for bed properties

        
      ENDDO !End loop i=1,nea

!---------------------------------------------------------------------
! Convert bed sediment properties from elements to node for outputs
!---------------------------------------------------------------------
      IF(myrank==0) WRITE(16,*)'SED: converting sed. arrays to nodes for outputs...'
      
      !' First, properties which cannot be equal to zero even if dry
      bed_d50n   = 0.0d0
      bed_fracn  = 0.0d0
      DO i = 1,np
        ta = 0.0d0
        DO j = 1,nne(i)
          ie = indel(j,i)
          ta = ta + area(ie)
          bed_d50n(i)   = bed_d50n(i)   + bottom(ie,isd50)*area(ie)
          DO ised = 1,ntracers
            bed_fracn(i,ised) = bed_fracn(i,ised) +                  &
            &                   bed_frac(1,ie,ised)*area(ie)
          ENDDO ! END loop ntracers
        ENDDO ! END loop nne
        IF(ta.EQ.0) THEN
          CALL parallel_abort('SEDIMENT: elem2nod (1)')
        ELSE
          bed_d50n(i)   = bed_d50n(i) / ta
          DO ised = 1,ntracers
            bed_fracn(i,ised) = bed_fracn(i,ised) / ta
          ENDDO ! END loop ntracers
        ENDIF
      ENDDO ! END loop i=1,np
      CALL exchange_p2d(bed_d50n(:))
      DO  ised=1,ntracers
        CALL exchange_p2d(bed_fracn(:,ised))
      ENDDO  ! END loop ntracers

      ! Then, bottom shear stress
      bed_taun  = 0.0d0
      DO i = 1,np
        IF(idry(i).EQ.1) CYCLE
        ta = 0.0d0
        DO j = 1,nne(i)
          ie = indel(j,i)
          IF(idry_e(ie).EQ.0)THEN
            ta = ta + area(ie)
            bed_taun(i) = bed_taun(i) + tau_wc(ie)*area(ie)
          ENDIF
        ENDDO ! END loop nne
        IF(ta.EQ.0)THEN
          CALL parallel_abort('SEDIMENT: elem2nod (2)')
        ELSE
          bed_taun(i) = bed_taun(i) / ta
        ENDIF
      ENDDO ! END loop i=1,np
      CALL exchange_p2d(bed_taun(:))

!---------------------------------------------------------------------
! - BCG: If only one bed layer, re-initialize bed thickness to initial
! thickness and update bed_mass according with.
!---------------------------------------------------------------------
      IF((sed_morph.EQ.2).AND.(Nbed.EQ.1))THEN
        k=1
        DO i=1,nea
          tmp=sum(bedthick_overall(elnode(1:3,i)))/3
          bed(1,i,ithck) = tmp !bedthick_overall
          DO ised=1,ntracers
            bed_mass(1,i,nnew,ised)=tmp*Srho(ised)*(1.0d0-bed(1,i,iporo))*  &
                                   &bed_frac(1,i,ised)
         
            if(bed_mass(1,i,nnew,ised)<0) then
              WRITE(errmsg,*)'SED3D: mass<0; ',bed_mass(1,i,nnew,ised),i,nnew,ised
              CALL parallel_abort(errmsg)
            endif
          ENDDO
        ENDDO
      ENDIF

!---------------------------------------------------------------------
! - Store old bed thickness.
!---------------------------------------------------------------------
!jl.  Do we not need to save the total bed thickness???
!     I didn't see where it is used in ROMS other than here although
!     it gets passed around a lot.
!
      IF (sed_morph.GE.1) THEN
!          do i=1,nea
!            bed_thick(i,nnew)=0.0_r8
!            do kbed=1,Nbed
!              bed_thick(i,nnew)=bed_thick(i,nnew)+bed(kbed,i,ithck)
!            end do
!          end do

!        Changing depth variation due to erosion/deposition from 
!        elements to nodes
        DO i=1,np
          IF(idry(i)==1) CYCLE
  
          tmp = 0
          ta = 0
          DO j=1,nne(i)
            ie = indel(j,i)
            IF(idry_e(ie)==0) THEN
              ta = ta+area(ie)
              tmp = tmp+hdep(ie)*area(ie)
            ENDIF
          ENDDO ! End loop nne

          IF(ta==0) THEN
            CALL parallel_abort('SEDIMENT: (4)')
          ELSE
            hdep_nd(i) = tmp/ta
          ENDIF

        ENDDO !End loop np

        ! Exchange
        CALL exchange_p2d(hdep_nd)
  
!---------------------------------------------------------------------
! Loop over all nodes to collect depth change due to bedload
! and suspended load:  dhnd=hdep_nd+hbed
!---------------------------------------------------------------------

        dhnd=0.0d0
        DO i=1,np
          IF (suspended_load == 1) THEN
            dhnd(i) = dhnd(i)+hdep_nd(i)
          ENDIF
          dhnd(i) = dhnd(i)+hbed(i)
        ENDDO

!---------------------------------------------------------------------
! Compute slope avalanching
!---------------------------------------------------------------------
        IF (slope_avalanching.EQ.1)THEN
          IF(myrank.EQ.0) WRITE(16,*)'start sed_avalanching'
          CALL sed_avalanching(dhnd)
          IF(myrank.EQ.0) WRITE(16,*)'done sed_avalanching'
        ENDIF

!---------------------------------------------------------------------
! Update depth according to bedload and suspended fluxes and 
! avalanching ;  update dp()
! Only if sed_morph==1, else, no application of depth changes:
! no morphology or BCG
!---------------------------------------------------------------------
        IF(sed_morph.EQ.1)THEN
          DO i=1,np
            dp(i)=dp(i)+dhnd(i)
!           Impose bare rock limit
            dp(i)=min(dp(i),dp00(i)+bedthick_overall(i))
          ENDDO !i 
!---------------------------------------------------------------------
! Exchange depths of nodes before going back to circulation code 
!---------------------------------------------------------------------
          CALL exchange_p2d(dp)
        ENDIF ! End test sed_morph == 1 (Not for BCG)

      ENDIF ! End sed_morph >=1 (Morpho or BCG)

!     Check interface arrays to main routine
!     Additonal arrays: rough_p (sed_friction.F90)
      cff=sum(dp)+sum(bdy_frc)
      if(cff/=cff) then
        WRITE(errmsg,*)'SED: has nan; ',sum(dp),sum(bdy_frc)
        CALL parallel_abort(errmsg)
      endif

!      IF(myrank==0) WRITE(mirror,*)'SED: Leaving sediment model'
!--------------------------------------------------------------------!      
      end SUBROUTINE sediment
     
