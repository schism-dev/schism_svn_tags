SUBROUTINE sed2d_main_mcml(it)
!--------------------------------------------------------------------
! This subroutine controls the workflow in sed2d for multi-class
! multi-layer approach (based on sed2d_main done by Guillaume Dodet)
!
! NB: this routine is adapted to node-centered method
!
! Author: Thomas Guerin (thomas.guerin@univ-lr.fr)
! Date: 18/03/2014
!
! History:
! 04/2014 - T.Guerin: Removed morphological time step (dtsed2d)
!                     because led to bug for multi-class mode (will
!                     be adapted in the future)
!--------------------------------------------------------------------

  USE elfe_glbl, ONLY : Cdp,dldxy,dp,dt,eta2,idry,idry_e,idry_s,ielg,&
                        ipgl,indel,iplg,isdel,isidenode,isidenei2,   &
                        elside,nea,nne,elnode,np,npa,ns,nsa,pi,      &
                        out_wwm,su2,sv2,timer_ns,rhosed,rkind,       &
                        ihydraulics,isblock_sd,errmsg
  USE elfe_msgp, ONLY : comm,exchange_e2d,exchange_p2d,exchange_s2d, &
                        ierr,myrank,parallel_abort,rtype,parallel_finalize
  USE hydraulic_structures, ONLY : nhtblocks
  USE sed2d_mod, ONLY : Cdsed,cflsed,d50,d90,diffac,dpdxy,dpdxy_e,   &
                        dtsed2d,h0_sed,idrag,idsed2d,ifilt,          &
                        imorpho,irough,iskip,islope,itrans,nskip,    &
                        poro,qav,qb,qb_e,qdt_e,qfilter,qramp,qs,qs_e,&
                        qsum_e,qtot,qtot_e,ufilter,vc_area,          &
                        transfac,z0_e,z0cr_e,z0sw_e,z0wr_e
  USE sed2d_mod, ONLY : u2_tmp,v2_tmp
  USE sed2d_mod, ONLY : bed_delta,d50moy,d90moy,F_class,h_inf,       &
                        h_lim_max,h_lim_min,h_top,h2,nb_class,       &
                        nb_layer
  USE sed2d_friction
  USE sed2d_transport
  USE sed2d_filter

  IMPLICIT NONE

  INCLUDE 'mpif.h'

!- Arguments --------------------------------------------------------
  INTEGER, INTENT(IN) :: it
!- Local variables --------------------------------------------------
  INTEGER :: i,id,iel,inode,iside,j,k,n,neigh,rank_tmp
  REAL(rkind) :: beta,D_star,d50_e,d90_e,hs,htot,kpeak,rampfac,Ucr_c,&
                 Ucr_w,u_star,tau_cr,theta_cr,tp,suru,surv,udir,ue,  &
                 uorb,uorbp,utmp,vtmp,ve,wdir,wdirc,wdirs,wlpeak,    &
                 wtmp,z0,tmp
  REAL(rkind) :: d50moy_e,d90moy_e,F_class_e,sum_dh,sum_F_class
  REAL(rkind), DIMENSION(nea,2) :: sum_qb_e,sum_qs_e,sum_qtot_e
  REAL(rkind), DIMENSION(npa) :: dp_filt,dp_tmp,dp1,neigh2
  REAL(rkind), DIMENSION(npa,nb_class) :: dh_class
  REAL(rkind), DIMENSION(nea) :: Cd_e,tmp_e
  REAL(rkind), DIMENSION(nsa) :: su2_tmp,su2_tmp1,su2_tmp2,sv2_tmp,  &
                                 sv2_tmp1,sv2_tmp2,tmp_s
  REAL(rkind), DIMENSION(npa,2) :: qdt
!- Parameter --------------------------------------------------------
  REAL(rkind), PARAMETER :: ONETHIRD = 1.d0/3.d0
!--------------------------------------------------------------------

! Sanity check
!  if(ishapiro/=1) call parallel_abort('Sed2d needs Shapiro filter')

!--------------------------------------------------------------------
!-                  Compute ramp factor                             -
!--------------------------------------------------------------------
  IF(it<iskip) THEN
    rampfac = 0.d0
  ELSE
    IF(qramp>0) THEN
      rampfac = TANH(2.d0*(it-iskip)*dt/qramp)
    ELSE
      rampfac = 1.d0
    ENDIF
  ENDIF
!--------------------------------------------------------------------
!-          Compute total transport at element center               -
!--------------------------------------------------------------------

#ifdef INCLUDE_TIMING
      wtmp = mpi_wtime() !start of timer
#endif

!- Side velocities filter before computing transport ----------------
  su2_tmp = su2(1,:)
  sv2_tmp = sv2(1,:)
  su2_tmp1 = su2(1,:)
  sv2_tmp1 = sv2(1,:)
  su2_tmp2 = su2(1,:)
  sv2_tmp2 = sv2(1,:)
  IF(ufilter /= 0) THEN
!- Skip if the element has a side at the wet-dry interface 
    DO i = 1,nsa
      IF(isdel(1,i)>0 .AND. isdel(2,i)>0) THEN !then resident and internal
        IF(idry_e(isdel(1,i))+idry_e(isdel(2,i))==1) THEN !wet/dry interface
          su2_tmp1(i) = 0.d0 
          sv2_tmp1(i) = 0.d0
        ELSE
          su2_tmp1(i) = su2(1,i) 
          sv2_tmp1(i) = sv2(1,i)
        ENDIF
      ENDIF
    ENDDO !nsa
    CALL exchange_s2d(su2_tmp1)
    CALL exchange_s2d(sv2_tmp1)
!- Iterative shapiro filter applied on the current field    
!JZ: lon/lat; 
    DO k = 1,ufilter
      DO i = 1,ns
        if(isdel(2,i)==0.or.idry_s(i)==1) CYCLE
        if(ihydraulics/=0.and.nhtblocks>0) then
          if(isblock_sd(1,i)/=0) cycle
        endif

        suru = 0.d0
        surv = 0.d0
        DO j = 1,4
          id = isidenei2(j,i)
          suru = suru+su2_tmp1(id)
          surv = surv+sv2_tmp1(id)
        ENDDO
        su2_tmp2(i) = su2_tmp1(i)+0.25d0*(suru-4.d0*su2_tmp1(i))
        sv2_tmp2(i) = sv2_tmp1(i)+0.25d0*(surv-4.d0*sv2_tmp1(i))
      ENDDO !ns
      DO i = 1,ns
        if(isdel(2,i)==0.or.idry_s(i)==1) CYCLE
        if(ihydraulics/=0.and.nhtblocks>0) then
          if(isblock_sd(1,i)/=0) cycle
        endif

        su2_tmp(i) = su2_tmp2(i)
        sv2_tmp(i) = sv2_tmp2(i)
      ENDDO
      CALL exchange_s2d(su2_tmp)
      CALL exchange_s2d(sv2_tmp)
      su2_tmp1 = su2_tmp
      sv2_tmp1 = sv2_tmp
    ENDDO !k=1,ufilter
    IF(myrank==0) WRITE(16,*)'done filtering velocities (sed2d)'
  ENDIF !ufilter/=0

#ifdef INCLUDE_TIMING
      timer_ns(11) = timer_ns(11)+mpi_wtime()-wtmp !end of timer 
#endif 

!- Initialization
  sum_qb_e   = 0.d0    !Sum of bed load transport over sediment classes (m3/s/m)
  sum_qs_e   = 0.d0    !Sum of suspended load transport over sediment classes (m3/s/m)
  sum_qtot_e = 0.d0    !Sum of total transport over sediment classes (m3/s/m)

!- Start loop over sediment classes
  DO k=1,nb_class

!- Initialize local arrays that are used to compute nodal outputs
    Cd_e       = 0.d0    !Drag coefficient (-)
    dpdxy_e    = 0.d0    !Bed slope (m/m)
    qb_e       = 0.d0    !Bed load transport (m3/s/m)
    qs_e       = 0.d0    !Suspended load transport (m3/s/m)
    qtot_e     = 0.d0    !Total transport (m3/s/m)
    z0_e       = 0.d0    !Roughness length z0 (m)
    d50_e = d50(k)    !D50 (m) corresponding to class k
    d90_e = d90(k)    !D90 (m) corresponding to class k

    DO i = 1,nea
#ifdef INCLUDE_TIMING
      wtmp = mpi_wtime()                         !start of timer
#endif 

!- Compute variables at element center: -----------------------------
      d50moy_e = 0.d0 !Mean D50 (m) for each sediment layer
      d90moy_e = 0.d0 !Mean D90 (m) for each sediment layer
      F_class_e = 0.d0 !Sediment class fraction (-)
      htot   = 0.d0  !Total depth (m)
      ue     = 0.d0  !X-velocity (m/s)
      ve     = 0.d0  !Y-velocity (m/s)
      uorb   = 0.d0  !Orbital velocity computed by wwm (m/s)
      hs     = 0.d0  !Significant wave height computed by wwm(m)
      tp     = 0.d0  !Peak period computed by wwm (s)
      wlpeak = 0.d0  !Wave length associated to peak period by wwm(m)
      uorbp  = 0.d0  !Peak orbital velocity (m/s)
      wdir   = 0.d0  !Mean average energy transport direction (rad)
      wdirc  = 0.d0 !cosinus value needed for calculating mean wave direction in each element
      wdirs  = 0.d0 !sinus value needed for calculating mean wave direction in each element

      dp1    = dp !Store depth before bed update for output purpose

      DO j = 1,3
        inode = elnode(j,i)
        iside = elside(j,i)

!- Total depth 
        htot = htot+(eta2(inode)+dp(inode))*ONETHIRD

!- Velocity 
        ue = ue+su2_tmp(iside)*ONETHIRD
        ve = ve+sv2_tmp(iside)*ONETHIRD

!- Mean D50 and D90
        d50moy_e = d50moy_e+d50moy(inode,1)*ONETHIRD
        d90moy_e = d90moy_e+d90moy(inode,1)*ONETHIRD

!- Sediment classes fraction
        F_class_e = F_class_e+F_class(inode,k,1)*ONETHIRD

!- bed slope  
        dpdxy_e(i,1) = dpdxy_e(i,1)+dp(inode)*dldxy(j,1,i)
        dpdxy_e(i,2) = dpdxy_e(i,2)+dp(inode)*dldxy(j,2,i)

!- Wave parameters 
#ifdef USE_WWM
        uorb = uorb+out_wwm(inode,22)*ONETHIRD
        hs = hs+out_wwm(inode,1)*ONETHIRD
        tp = tp+out_wwm(inode,12)*ONETHIRD
        wlpeak = wlpeak+out_wwm(inode,17)*ONETHIRD
        wdirc = wdirc+(DCOS((270.d0-out_wwm(inode,9))          &
                   *pi/180.d0))*ONETHIRD
        wdirs = wdirs+(DSIN((270.d0-out_wwm(inode,9))          &
                   *pi/180.d0))*ONETHIRD
#endif

!- Drag coefficient
        IF(idrag==1)THEN
          Cd_e(i) = Cd_e(i)+Cdp(inode)*ONETHIRD
        ENDIF
      ENDDO! j = 1,3

#ifdef USE_WWM
      wdir = DATAN2(wdirs,wdirc)
#endif

#ifdef INCLUDE_TIMING
      timer_ns(4) = timer_ns(4)+mpi_wtime()-wtmp !end of timer 
      wtmp = mpi_wtime()                         !start of timer
#endif 

!     Do rest for wet elem. only
      IF((idry_e(i) == 1).OR.(htot.LE.h0_sed)) CYCLE

!- Compute threshold parameters, roughness and drag coefficient  ----
! NB: threshold parameters are computed from d50(class), whereas
!     roughness and drag coefficient are computed from mean d50
!--------------------------------------------------------------------
      CALL compute_thresh(d50_e,d90_e,htot,tp,D_star,Ucr_c,Ucr_w,  &
                         &theta_cr,tau_cr)

      SELECT CASE(irough)
       CASE(1) !Skin-friction only
        z0_e(i) = d50moy_e/12.d0
       CASE(2) !Bedform associated roughness (Soulsby, 1997)
        CALL bedform_predictor_s97(d50moy_e,htot,ue,ve,uorb,tp,      &
            &tau_cr,theta_cr,z0cr_e(i),z0wr_e(i),z0sw_e(i),z0_e(i))
       CASE(3) !Bedform associated roughness (Van-Rijn, 2007)
        CALL bedform_predictor_vr07(d50moy_e,htot,ue,ve,uorb,        &
            &z0cr_e(i),z0wr_e(i),z0sw_e(i),z0_e(i))
      END SELECT

      IF (ABS(idrag).GT.1) THEN
        CALL compute_drag(d50moy_e,htot,z0_e(i),ue,ve,idrag,Cd_e(i))
      ENDIF

#ifdef INCLUDE_TIMING
      timer_ns(5) = timer_ns(5)+mpi_wtime()-wtmp !end of timer 
      wtmp = mpi_wtime()                         !start of timer
#endif 

!- Test if parameters are physically correct ------------------------
#ifdef USE_DEBUG
      IF(htot<0.d0 .OR. DSQRT(ue*ue+ve*ve)>7.5d0 .OR. d50_e>0.05d0.OR.&
        &d90_e>0.10d0 .OR.DSQRT(dpdxy_e(i,1)**2.d0+dpdxy_e(i,2)**2.d0)&
        &>5.d0 .OR. uorb>10.d0 .OR. hs>10.d0 .OR. tp>30.d0 .OR.       &
        &wlpeak > 1500.d0 .OR. Cd_e(i) > 10.d0 .OR.  z0 > 0.5d0) THEN
        WRITE(12,*)'Warning (1) for it:',it,' and el: ', ielg(i)
        WRITE(12,*)'Check following values: '
        WRITE(12,*)'total depth = ',htot
        WRITE(12,*)'current velocity = ',DSQRT(ue*ue+ve*ve)
        WRITE(12,*)'d50 = ',d50_e
        WRITE(12,*)'d90 = ',d90_e
        WRITE(12,*)'slope = ',DSQRT(dpdxy_e(i,1)**2.d0+dpdxy_e(i,2)**2.d0)
        WRITE(12,*)'orbital velocity = ',uorb
        WRITE(12,*)'hs = ',hs
        WRITE(12,*)'tp = ',tp
        WRITE(12,*)'wlpeak = ',wlpeak
        WRITE(12,*)'cd = ', cd_e(i)
        WRITE(12,*)'z0 = ',z0
!        CALL parallel_abort('Sed2d: chech value in outputs/nonfatal')
      ENDIF
#endif

!- Apply user-defined formula ---------------------------------------

!- Compute bed slope in streamwise direction used for slope effect on 
!  total transport (Soulsby, 1997). Positive if flows runs uphill.
      IF(islope == 1) THEN
        udir = DATAN2(ve,ue)
        beta = -(dpdxy_e(i,1)*DCOS(udir)+dpdxy_e(i,2)*DSIN(udir))
        beta = MIN(beta,0.6d0) !May lead negative transport otherwise
      ELSE
        beta = 0.d0
      ENDIF

      SELECT CASE(itrans)
       CASE(1) !Engelund-Hansen (1967)
        CALL eh67(Cd_e(i),d50_e,ue,ve,beta,qtot_e(i,:))

       CASE(2) !Ackers and White (1973)
        CALL aw73(Cd_e(i),d50_e,ue,ve,htot,beta,D_star,qtot_e(i,:))

       CASE(3) !Soulsby - Van Rijn (Soulsby, 1997)
        CALL svr97_bedl(Cd_e(i),d50_e,ue,ve,htot,dpdxy_e(i,:),Ucr_c, &
                       &tau_cr,uorb,qb_e(i,:))
        CALL svr97_susp(Cd_e(i),d50_e,ue,ve,D_star,Ucr_c,uorb,       &
                       &qs_e(i,:))
        qtot_e(i,:) = (qs_e(i,:)+qb_e(i,:))*(1.d0-1.6d0*beta) !last factor>0

       CASE(4) !Van-Rijn (2007) needs peak orbital velocity
!#ifdef USE_WWM
!       IF((tp.GT.0.d0) .AND. (wlpeak.GT.0.d0)) THEN
!         kpeak = 2.d0*pi/wlpeak
!         uorbp = pi*hs/(tp*DSINH(kpeak*htot))
!       ENDIF
!#endif
        CALL vr07_bedl(Cd_e(i),d50_e,ue,ve,htot,dpdxy_e(i,:),uorb,  &
                      &Ucr_c,Ucr_w,tau_cr,qb_e(i,:))
        CALL vr07_susp(d50_e,ue,ve,uorb,D_star,Ucr_c,Ucr_w,         &
                      &qs_e(i,:))
        qtot_e(i,:) = (qs_e(i,:)+qb_e(i,:))*(1.d0-1.6d0*beta)

       CASE(5) !Camenen and Larson (2011)
        CALL cl11_tot(d50_e,D_star,dpdxy_e(i,:),hs,htot,theta_cr,tp, &
                     &ue,ve,wdir,wlpeak,qb_e(i,:),qs_e(i,:))
        qtot_e(i,:) = (qs_e(i,:)+qb_e(i,:))*(1.d0-1.6d0*beta)
      END SELECT

!- Multiply total transport by sediment class fraction
      qtot_e(i,:) = F_class_e*qtot_e(i,:)

!- Apply ramp, diffusion and correction factor if required-----------
      qtot_e(i,1) = rampfac*transfac*(qtot_e(i,1)+(1.d0-poro)*diffac* &
                   &DABS(qtot_e(i,1))*dpdxy_e(i,1))
      qtot_e(i,2) = rampfac*transfac*(qtot_e(i,2)+(1.d0-poro)*diffac* &
                   &DABS(qtot_e(i,2))*dpdxy_e(i,2))

    ENDDO !nea

#ifdef INCLUDE_TIMING                                                                     
    timer_ns(6) = timer_ns(6)+mpi_wtime()-wtmp !end of timer                              
    wtmp = mpi_wtime()                         !start of timer                            
#endif

    CALL exchange_e2d(qtot_e(:,1))  !check if necessary
    CALL exchange_e2d(qtot_e(:,2))  !check if necessary
    CALL exchange_e2d(qb_e(:,1))    !check if necessary                                   
    CALL exchange_e2d(qb_e(:,2))    !check if necessary                                   
    CALL exchange_e2d(qs_e(:,1))    !check if necessary                                   
    CALL exchange_e2d(qs_e(:,2))    !check if necessary
    CALL exchange_e2d(dpdxy_e(:,1)) !check if necessary                                   
    CALL exchange_e2d(dpdxy_e(:,2)) !check if necessary                                   
    CALL exchange_e2d(Cd_e(:))      !check if necessary
    CALL exchange_e2d(z0_e(:))      !check if necessary

!- Apply diffusive filter on total transport if required ------------
    IF(qfilter == 1) THEN
      CALL sed2d_filter_diffu(qtot_e(:,1),tmp_e,nea)
      qtot_e(:,1) = tmp_e
      CALL sed2d_filter_diffu(qtot_e(:,2),tmp_e,nea)
      qtot_e(:,2) = tmp_e
      CALL exchange_e2d(qtot_e(:,1))  !check if necessary                                   
      CALL exchange_e2d(qtot_e(:,2))  !check if necessary
      IF(myrank==0) WRITE(16,*)'done filtering transport rate (sed2d):',k
    ENDIF

!- Sum over sediment classes (intermediate variables)
    sum_qb_e(:,:) = sum_qb_e(:,:) + qb_e(:,:)
    sum_qs_e(:,:) = sum_qs_e(:,:) + qs_e(:,:)
    sum_qtot_e(:,1) = sum_qtot_e(:,1) + qtot_e(:,1)
    sum_qtot_e(:,2) = sum_qtot_e(:,2) + qtot_e(:,2)

!    CALL exchange_e2d(sum_qtot_e(:,1)) !check if necessary
!    CALL exchange_e2d(sum_qtot_e(:,2)) !check if necessary

#ifdef INCLUDE_TIMING
    timer_ns(7) = timer_ns(7)+mpi_wtime()-wtmp !end of timer 
    wtmp = mpi_wtime()                         !start of timer
#endif 

!--------------------------------------------------------------------
!-                    Compute bed change at nodes                   -
!--------------------------------------------------------------------
    qdt_e(:,1) = qtot_e(:,1)*dt
    qdt_e(:,2) = qtot_e(:,2)*dt
!    CALL exchange_e2d(qdt_e(:,1)) !check if necessary
!    CALL exchange_e2d(qdt_e(:,2)) !check if necessary

    IF(imorpho == 1 .AND. it>iskip) THEN
      CALL sed2d_morpho(it)
      !Save bottom evolution for the current sediment class
      dh_class(:,k) = bed_delta(:)
      CALL exchange_p2d(dh_class(:,k)) !check if necessary
    ENDIF
  ENDDO !k=1,nb_class

!- Final values for bed-load, suspended load, and total transport
!  qb_e(:,:)   = sum_qb_e(:,:)
!  qs_e(:,:)   = sum_qs_e(:,:)
  qtot_e(:,1) = sum_qtot_e(:,1)
  qtot_e(:,2) = sum_qtot_e(:,2)
!  CALL exchange_e2d(qtot_e(:,1)) !check if necessary
!  CALL exchange_e2d(qtot_e(:,2)) !check if necessary

  IF(imorpho == 1 .AND. it>iskip) THEN
    DO i=1,npa
      IF (idry(i)==1) CYCLE

!- Compute total bottom evolution
      sum_dh = 0.d0 !need to be initialized
      DO k=1,nb_class
        sum_dh = sum_dh + dh_class(i,k)
      ENDDO

!- Compute new bathymetry
      dp(i) = dp(i) + sum_dh

!- Compute new sediment class fractions
      DO k=1,nb_class
        IF (sum_dh.GT.0.d0) THEN !EROSION case

!YJZ: this check may be too restrictive?
          IF(sum_dh.GE.h_lim_min) THEN
             write(errmsg,*)'negative thickness for layer #2: &
                 &need to increase h_lim_min, or decrease time step, or decrease dtsed2d:',sum_dh,h_lim_min
             CALL parallel_abort(errmsg)
          ENDIF

          !New fractions for surface layer; draw from layer #2 so top layer
          !thickness stays same
          F_class(i,k,1) = (h_top*F_class(i,k,1) +                 &
                    &sum_dh*F_class(i,k,2) - dh_class(i,k)) /h_top

          IF(F_class(i,k,1).LT.0.d0) F_class(i,k,1)=0.d0


        ELSEIF (sum_dh.LT.0.d0) THEN !DEPOSITION case

!YJZ: this check may be too restrictive?
          IF(ABS(sum_dh).GE.h_lim_max) THEN
            write(errmsg,*)'layer #2 thickness too large: &
              &need to increase h_lim_max, or decrease time step, or decrease dtsed2d:',sum_dh,h_lim_max
            CALL parallel_abort(errmsg)
          ENDIF

          !New fractions for surface layer
          !Deposit into layer #1 (-dh), and swap out same amount to layer #2
          F_class(i,k,1) = ((h_top+sum_dh)*F_class(i,k,1)            &
                     & - dh_class(i,k)) /h_top

          IF(F_class(i,k,1).LT.0.d0) F_class(i,k,1)=0.d0

          !New sub-surface layer fractions
          F_class(i,k,2) = (h2(i)*F_class(i,k,2)                     &
                     & - sum_dh*F_class(i,k,1)) /(h2(i)-sum_dh) !>=0
        ENDIF !sum_dh sign
      ENDDO !k=1,nb_class

!-    Update layer below 1
      IF (sum_dh.GT.0.d0) THEN !EROSION case
        !New sub-surface layer thickness
        h2(i) = h2(i) - sum_dh !should >0 b/cos sum_dh<=h_lim_min

        IF (h2(i).LT.h_lim_min) THEN !Merging sub-surface layer with layer #3
          !New sub-surface layer fractions and thickness
          if(h2(i)+h_inf==0) call parallel_abort('SED2D: (3)')
          F_class(i,:,2)=(h_inf*F_class(i,:,3)+h2(i)*F_class(i,:,2))/(h2(i)+h_inf) !>0
          h2(i) = h2(i) + h_inf

          !Adding new layer => updating fractions - implies infinite supply
          !of sediment
          DO n=3,nb_layer-1
            F_class(i,:,n) = F_class(i,:,n+1)
          ENDDO
        ENDIF !h2(i)

      ELSE !DEPOSITION
        !New sub-surface layer thickness
        h2(i) = h2(i) - sum_dh !sum_dh<=0

        IF (h2(i).GT.h_lim_max) THEN !Splitting sub-surface layer into 2 layers
          !New sub-surface layer thickness
          h2(i) = h2(i) - h_inf 
          if(h2(i)<h_lim_min) then
            write(errmsg,*)'SED2D: check h_lim_max etc: ',h2(i),h_lim_min,h_lim_max
            call parallel_abort(errmsg)      
          endif

          !Removing lowest layer => updating fractions
          DO n=1,nb_layer-2
            F_class(i,:,nb_layer-n+1) = F_class(i,:,nb_layer-n)
          ENDDO
        ENDIF !h2(i)>
      ENDIF !sum_dh

!- Compute sum of fractions and correct fractions if sum /= 1
      DO n=1,nb_layer
        sum_F_class = 0.d0
        DO k=1,nb_class
          sum_F_class = sum_F_class + F_class(i,k,n)
        ENDDO

        IF (sum_F_class/=1.d0) THEN
          if(sum_F_class==0) then
            write(errmsg,*)'SED2D: all eroded; ',n,iplg(i)
            call parallel_abort(errmsg)
          endif
          F_class(i,:,n) = F_class(i,:,n)/sum_F_class
        ENDIF
      ENDDO !nb_layer

!- Compute new mean d50 (weighted geometric mean)
      DO n=1,nb_layer
        d50moy(i,n) = d50(1)**F_class(i,1,n)
        DO k=2,nb_class
          d50moy(i,n) = d50moy(i,n) * d50(k)**F_class(i,k,n)
        ENDDO
      ENDDO !nb_layer

!- Compute new mean d90 (=2.5*d50moy approximation)
      d90moy(i,:) = 2.5d0*d50moy(i,:)
    ENDDO !i=1,npa
    CALL exchange_p2d(dp) !check if necessary
    CALL exchange_p2d(h2) !check if necessary
    DO n=1,nb_layer
       CALL exchange_p2d(d50moy(:,n))
       CALL exchange_p2d(d90moy(:,n))
       DO k=1,nb_class
          CALL exchange_p2d(F_class(:,k,n))
       ENDDO !nb_class
    ENDDO !nb_layer
  ENDIF !imorpho=1 .and. it>iskip

#ifdef INCLUDE_TIMING
  timer_ns(9) = timer_ns(9)+mpi_wtime()-wtmp !end of timer 
  wtmp = mpi_wtime()                         !start of timer
#endif 
!--------------------------------------------------------------------
!-                      Apply user-defined filter                   -
!--------------------------------------------------------------------
  IF ((imorpho==1).AND.(it>iskip).AND.(ifilt>0)) THEN
    SELECT CASE(ifilt)
     CASE(1)
       CALL sed2d_filter_extrema(dp,dp_filt)
     CASE(2)
       CALL sed2d_filter_slope(dp,dp_filt)
     CASE(3)
       CALL sed2d_filter_slope(dp,dp_tmp)
       CALL sed2d_filter_extrema(dp_tmp,dp_filt)
    END SELECT
    dp = dp_filt
    CALL exchange_p2d(dp) !check if necessary
    IF(myrank==0)WRITE(16,*)'done filtering new bathy. (sed2d)'
  ENDIF !ifilt>0

!--------------------------------------------------------------------
!-                      Check if new bathymetry /=NaN               -
!--------------------------------------------------------------------
  DO i = 1,npa
    IF(dp(i)/=dp(i)) THEN
      WRITE(errmsg,*)'Sed2d: NaN in depth at step:',it,' and node: ', iplg(i),';dp = ',dp(i)
      CALL parallel_abort(errmsg)
    ENDIF
  ENDDO !i


#ifdef INCLUDE_TIMING
  timer_ns(10) = timer_ns(10)+mpi_wtime()-wtmp !end of timer 
  wtmp = mpi_wtime()                         !start of timer
#endif 

!--------------------------------------------------------------------
!-       Interpolate variables at node for output purpose and also for 
!        passing back Cdsed to selfe_step 
!--------------------------------------------------------------------
  Cdsed = 0.d0
  cflsed = 0.d0
  dpdxy = 0.d0
  qav = 0.d0
  qdt = 0.d0
  qb = 0.d0
  qs = 0.d0
  qtot = 0.d0

  DO i = 1,np
    IF(idry(i) == 1) CYCLE

    neigh = 0
    DO j = 1,nne(i)
      iel = indel(j,i)
      Cdsed(i)   = Cdsed(i)+Cd_e(iel)
      dpdxy(i,:) = dpdxy(i,:)+dpdxy_e(iel,:)
      qav(i,:)   = qav(i,:)+qdt_e(iel,:)
      qdt(i,:)   = qdt(i,:)+qdt_e(iel,:)
      qb(i,:)    = qb(i,:)+sum_qb_e(iel,:) !qb_e(iel,:)
      qs(i,:)    = qs(i,:)+sum_qs_e(iel,:) !qs_e(iel,:)
      qtot(i,:)  = qtot(i,:)+qtot_e(iel,:)
      neigh = neigh+1
    ENDDO !j

    if(neigh==0.or.eta2(i)+dp1(i)==0) call parallel_abort('SED2D: div. by 0 (12)')
    Cdsed(i)   = Cdsed(i)/neigh
    dpdxy(i,1) = dpdxy(i,1)/neigh
    dpdxy(i,2) = dpdxy(i,2)/neigh
    qav(i,:)   = rhosed*qav(i,:)/(neigh*dt*dtsed2d) !(kg/m/s)
    qdt(i,:)   = qdt(i,:)/neigh                     !(m3/m)
    qb(i,:)    = rhosed*qb(i,:)/neigh               !(kg/m/s)
    qs(i,:)    = rhosed*qs(i,:)/neigh               !(kg/m/s)
    qtot(i,:)  = rhosed*qtot(i,:)/neigh             !(kg/m/s)
    cflsed(i) = 3.d0*DSQRT(qdt(i,1)*qdt(i,1)+qdt(i,2)*qdt(i,2))/  &
     &((eta2(i)+dp1(i))*DSQRT(vc_area(i)))
  ENDDO !np

  CALL exchange_p2d(Cdsed(:))
  DO i=1,2
    CALL exchange_p2d(dpdxy(:,i))
    CALL exchange_p2d(qav(:,i))
    CALL exchange_p2d(qdt(:,i))
    CALL exchange_p2d(qb(:,i))
    CALL exchange_p2d(qs(:,i))
    CALL exchange_p2d(qtot(:,i))
  ENDDO
  CALL exchange_p2d(cflsed(:))

  u2_tmp = 0.d0
  v2_tmp = 0.d0
  neigh2 = 0.d0
  DO i = 1,nsa
    DO j = 1,2
      inode = isidenode(j,i)
      IF(idry(inode) == 1) CYCLE
        u2_tmp(inode) = u2_tmp(inode)+su2_tmp(i)
        v2_tmp(inode) = v2_tmp(inode)+sv2_tmp(i)
        neigh2(inode) = neigh2(inode)+1
    ENDDO !j=1,2
  ENDDO !nsa
  CALL exchange_p2d(u2_tmp(:)) !check if necessary
  CALL exchange_p2d(v2_tmp(:)) !check if necessary
  CALL exchange_p2d(neigh2(:)) !check if necessary

  DO i = 1,npa
    IF(idry(i) == 1) CYCLE
    if(neigh2(i)==0) call parallel_abort('sed2d_main_node: impossible (1)')
    u2_tmp(i) = u2_tmp(i)/neigh2(i)
    v2_tmp(i) = v2_tmp(i)/neigh2(i)   
  ENDDO !npa
 
  tmp=sum(dp)
  WRITE(12,*)'SED2D, it=',it,' ; sum of depth=',tmp
  IF(myrank==0) WRITE(16,*)'done computing total sediment load (sed2d) and &
    &exiting sed2d_main_mcml; check nonfatal* for sum of depths'
  if(tmp/=tmp) call parallel_abort('sed2d_main_node: NaN in sum')

#ifdef INCLUDE_TIMING
  timer_ns(8) = timer_ns(8)+mpi_wtime()-wtmp !end of timer 
#endif 

  END SUBROUTINE sed2d_main_mcml
