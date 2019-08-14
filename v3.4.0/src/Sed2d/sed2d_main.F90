SUBROUTINE sed2d_main(it)
!--------------------------------------------------------------------
! This subroutine controls the workflow in sed2d. 
!
! Author: guillaume dodet (gdodet01@univ-lr.fr, gdodet@lnec.pt)
! Date: 19/02/2013
!
! History:
! 03/2013 - G.Dodet: - Corrected bugs in node centered FE-method; 
!                    - Modified filtering section;
!                    - Added space-variable d50 input and Cdsed 
!                      output;
! 04/2013 - G.Dodet: - Added nskip parameter to skip first iterations
!                      before updating bed level;
!                    - Split in 2 routines for each numerical method;
!                    - Added diffusive filter on total transport;
!--------------------------------------------------------------------

  USE sed2d_mod, ONLY : imeth

  IMPLICIT NONE

  INCLUDE 'mpif.h'

!- Arguments --------------------------------------------------------
  INTEGER, INTENT(IN) :: it
!--------------------------------------------------------------------
 
  IF(imeth == 1) THEN        !Node-centered method
    CALL sed2d_main_node(it)
  ELSE                       !Element-centered method
    CALL sed2d_main_el(it) 
  ENDIF

  CONTAINS

!--------------------------------------------------------------------
  SUBROUTINE sed2d_main_node(it)
!--------------------------------------------------------------------
! Main sed2d routine adapted to node-centered method (same as 
! previous sed2d_main before splitting).
!
! Author: guillaume dodet (gdodet01@univ-lr.fr, gdodet@lnec.pt)
! Date: 29/04/2013
!
! History:
! 06/2013 - G.Dodet:  - Added filter on velocity before computing 
!                       sediment fluxes to prevent development of
!                       oscillations along the coastline induced by 
!                       wave-force/barotropic pressure gradient 
!                       inconsistency (temporary solution);
!                     - Added diffusive filter on sediment fluxes
!                     - Added bedform predictors and new roughness;
!                     - Added timers;
! 07/2013 - G.Dodet:  - Converted output transport rate in kg/m/s;
!         - G.Dodet:  Added iterative shapiro filter for currents
!         - T.Guerin: Added Camenen and Larson (2011)
! 04/2014 - T.Guerin: - Removed part for getting offshore wave
!                       parameters (needed for Camemen and Larson
!                       formula and wave asymmetry calculation) and
!                       replaced it by local parameters computation
!                     - Added initialization of Cd_e, z0_e
!                     - Removed morphological time step (dtsed2d) to
!                       be consistent with multi-class mode (will be
!                       adapted in the future)
!--------------------------------------------------------------------


  USE elfe_glbl, ONLY : Cdp,dldxy,dp,dt,eta2,idry,idry_e,idry_s,ielg,   &
                        ipgl,indel,iplg,isdel,isidenode,isidenei2,elside,nea, &
                        nne,elnode,np,npa,ns,nsa,pi,out_wwm,su2,sv2,     &
                        timer_ns,rhosed,rkind
  USE elfe_msgp, ONLY : comm,exchange_e2d,exchange_p2d,exchange_s2d, &
                        ierr,myrank,parallel_abort,rtype
  USE sed2d_mod, ONLY : Cdsed,cflsed,d50,d90,diffac,dpdxy,dpdxy_e,   &
                        dtsed2d,h0_sed,idrag,idsed2d,ifilt,          &
                        imorpho,irough,iskip,islope,itrans,nskip,    &
                        poro,qav,qb,qb_e,qdt_e,qfilter,qramp,qs,qs_e,&
                        qsum_e,qtot,qtot_e,ufilter,vc_area,          &
                        transfac,z0_e,z0cr_e,z0sw_e,z0wr_e
  USE sed2d_mod, ONLY : u2_tmp,v2_tmp
  USE sed2d_friction
  USE sed2d_transport
  USE sed2d_filter

  IMPLICIT NONE

!- Arguments --------------------------------------------------------
  INTEGER, INTENT(IN) :: it
!- Local variables --------------------------------------------------
  INTEGER :: i,id,iel,inode,iside,j,k,neigh,rank_tmp
  REAL(rkind) :: beta,D_star,d50_e,d90_e,hs,htot,kpeak,rampfac,Ucr_c,&
                 Ucr_w,u_star,tau_cr,theta_cr,tp,suru,surv,udir,ue,  &
                 uorb,uorbp,utmp,vtmp,ve,wdir,wdirc,wdirs,wlpeak,    &
                 wtmp,z0
  REAL(rkind), DIMENSION(npa) :: dp_filt,dp_tmp,dp1,neigh2
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
!- If the element has a side at the wet-dry interface cycle
    DO i = 1,nsa
      IF(isdel(1,i)>0 .AND. isdel(2,i)>0) THEN !then resident
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
!JZ: lon/lat; hydraulics
    DO k = 1,ufilter
      DO i = 1,ns
        if(isdel(2,i)==0.or.idry_s(i)==1) CYCLE
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

!- Initialize local arrays that are used to compute nodal outputs
  Cd_e    = 0.d0    !Drag coefficient (-)
  dpdxy_e = 0.d0    !Bed slope (m/m)
  qb_e    = 0.d0    !Bed load transport (m3/s/m)
  qs_e    = 0.d0    !Suspended load transport (m3/s/m)
  qtot_e  = 0.d0    !Total transport (m3/s/m)
  z0_e    = 0.d0    !Roughness length z0 (m) 

  DO i = 1,nea
#ifdef INCLUDE_TIMING
      wtmp = mpi_wtime()                         !start of timer
#endif 

!- Compute variables at element center: -----------------------------
     d50_e  = 0.d0  !Median grain size (m)
     d90_e  = 0.d0  !D90 grain size (m)
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

!- Grain size
        d50_e = d50_e+d50(inode)*ONETHIRD
        d90_e = d90_e+d90(inode)*ONETHIRD

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

     IF((idry_e(i) == 1).OR.(htot.LE.h0_sed)) CYCLE

!- Compute threshold parameters, roughness and drag coefficient  ----
     CALL compute_thresh(d50_e,d90_e,htot,tp,D_star,Ucr_c,Ucr_w,     &
                         theta_cr,tau_cr)

     SELECT CASE(irough)
     CASE(1) !Skin-friction only
       z0_e(i) = d50_e/12.d0
     CASE(2) !Bedform associated roughness (Soulsby, 1997)
       CALL bedform_predictor_s97(d50_e,htot,ue,ve,uorb,tp,tau_cr,   &
            theta_cr,z0cr_e(i),z0wr_e(i),z0sw_e(i),z0_e(i))
     CASE(3) !Bedform associated roughness (Van-Rijn, 2007)
       CALL bedform_predictor_vr07(d50_e,htot,ue,ve,uorb,z0cr_e(i),  &
            z0wr_e(i),z0sw_e(i),z0_e(i))
     END SELECT

     IF(ABS(idrag).GT.1)CALL compute_drag(d50_e,htot,z0_e(i),ue,ve,  &
                                          idrag,Cd_e(i))

#ifdef INCLUDE_TIMING
      timer_ns(5) = timer_ns(5)+mpi_wtime()-wtmp !end of timer 
      wtmp = mpi_wtime()                         !start of timer
#endif 

!- Test if parameters are physically correct ------------------------
#ifdef USE_DEBUG
     IF(htot<0.d0 .OR. DSQRT(ue*ue+ve*ve)>7.5d0 .OR. d50_e>0.05d0.OR.&
        d90_e>0.10d0 .OR.DSQRT(dpdxy_e(i,1)**2.d0+dpdxy_e(i,2)**2.d0)&
        >5.d0 .OR. uorb>10.d0 .OR. hs>10.d0 .OR. tp>30.d0 .OR.       &
        wlpeak > 1500.d0 .OR. Cd_e(i) > 10.d0 .OR.  z0 > 0.5d0) THEN
        WRITE(12,*)'Warning (1) for it:',it,' and el: ', ielg(i)
        WRITE(12,*)'Check following values: '
        WRITE(12,*)'total depth = ',htot
        WRITE(12,*)'current velocity = ',DSQRT(ue*ue+ve*ve)
        WRITE(12,*)'d50 = ',d50_e
        WRITE(12,*)'d90 = ',d90_e
        WRITE(12,*)'slope = ',DSQRT(dpdxy_e(i,1)**2.d0+dpdxy_e(i,2)**&
                              2.d0)
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
       CALL svr97_bedl(Cd_e(i),d50_e,ue,ve,htot,dpdxy_e(i,:),Ucr_c,  &
                       tau_cr,uorb,qb_e(i,:))
       CALL svr97_susp(Cd_e(i),d50_e,ue,ve,D_star,Ucr_c,uorb,        &
                       qs_e(i,:))
       qtot_e(i,:) = (qs_e(i,:)+qb_e(i,:))*(1.d0-1.6d0*beta)

      CASE(4) !Van-Rijn (2007) needs peak orbital velocity
!#ifdef USE_WWM
!       IF((tp.GT.0.d0) .AND. (wlpeak.GT.0.d0)) THEN
!         kpeak = 2.d0*pi/wlpeak
!         uorbp = pi*hs/(tp*DSINH(kpeak*htot))
!       ENDIF
!#endif
       CALL vr07_bedl(Cd_e(i),d50_e,ue,ve,htot,dpdxy_e(i,:),uorb,   &
                      Ucr_c,Ucr_w,tau_cr,qb_e(i,:))
       CALL vr07_susp(d50_e,ue,ve,uorb,D_star,Ucr_c,Ucr_w,qs_e(i,:))
       qtot_e(i,:) = (qs_e(i,:)+qb_e(i,:))*(1.d0-1.6d0*beta)

      CASE(5) !Camenen and Larson (2011)
       CALL cl11_tot(d50_e,D_star,dpdxy_e(i,:),hs,htot,theta_cr,tp,  &
                     ue,ve,wdir,wlpeak,qb_e(i,:),qs_e(i,:))
       qtot_e(i,:) = (qs_e(i,:)+qb_e(i,:))*(1.d0-1.6d0*beta)
     END SELECT

!- Apply ramp, diffusion and correction factor if required-----------
     qtot_e(i,1) = rampfac*transfac*(qtot_e(i,1)+(1.d0-poro)*diffac* &
                   DABS(qtot_e(i,1))*dpdxy_e(i,1))
     qtot_e(i,2) = rampfac*transfac*(qtot_e(i,2)+(1.d0-poro)*diffac* &
                   DABS(qtot_e(i,2))*dpdxy_e(i,2))
  ENDDO !nea

#ifdef INCLUDE_TIMING
      timer_ns(6) = timer_ns(6)+mpi_wtime()-wtmp !end of timer
      wtmp = mpi_wtime()
#endif

  CALL exchange_e2d(qb_e(:,1))    !check if necessary                                     
  CALL exchange_e2d(qb_e(:,2))    !check if necessary                                     
  CALL exchange_e2d(qs_e(:,1))    !check if necessary                                     
  CALL exchange_e2d(qs_e(:,2))    !check if necessary
  CALL exchange_e2d(qtot_e(:,1))  !check if necessary
  CALL exchange_e2d(qtot_e(:,2))  !check if necessary
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
    IF(myrank==0) WRITE(16,*)'done filtering transport rate (sed2d)'
  ENDIF

  CALL exchange_e2d(qtot_e(:,1))  !check if necessary
  CALL exchange_e2d(qtot_e(:,2))  !check if necessary

#ifdef INCLUDE_TIMING
      timer_ns(7) = timer_ns(7)+mpi_wtime()-wtmp !end of timer 
      wtmp = mpi_wtime()                         !start of timer
#endif 

!--------------------------------------------------------------------
!-                    Compute bed change at nodes                   -
!--------------------------------------------------------------------
  qdt_e(:,1) = qtot_e(:,1)*dt
  qdt_e(:,2) = qtot_e(:,2)*dt
  CALL exchange_e2d(qdt_e(:,1)) !check if necessary
  CALL exchange_e2d(qdt_e(:,2)) !check if necessary

  IF(imorpho == 1 .AND. it>iskip) THEN
    CALL sed2d_morpho(it)
    IF(myrank==0)WRITE(16,*)'done computing new bathymetry (sed2d)'
#ifdef INCLUDE_TIMING
      timer_ns(9) = timer_ns(9)+mpi_wtime()-wtmp !end of timer 
      wtmp = mpi_wtime()                         !start of timer
#endif 
!--------------------------------------------------------------------
!-                      Apply user-defined filter                   -
!--------------------------------------------------------------------
    IF(ifilt>0 .AND. MOD(it,nskip*dtsed2d)==0) THEN
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
        CALL exchange_p2d(dp(:)) !check if necessary
        IF(myrank==0)WRITE(16,*)'done filtering new bathy. (sed2d)'
    ENDIF !ifilt>0
  ENDIF !imoprho == 1

!--------------------------------------------------------------------
!-                      Check if new bathymetry /=NaN               -
!--------------------------------------------------------------------
  DO i = 1,npa
    IF(dp(i)/=dp(i)) THEN
      WRITE(12,*)'Warning (2) for it:',it,' and node: ', iplg(i)
      WRITE(12,*)'dp = ',dp(i)
      CALL parallel_abort('Sed2d: depth is NaN')
    ENDIF
  ENDDO
#ifdef INCLUDE_TIMING
  timer_ns(10) = timer_ns(10)+mpi_wtime()-wtmp !end of timer 
  wtmp = mpi_wtime()                         !start of timer
#endif 

!--------------------------------------------------------------------
!-       Interpolate variables at node for output purpose only      -
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
          qb(i,:)    = qb(i,:)+qb_e(iel,:)
          qs(i,:)    = qs(i,:)+qs_e(iel,:)
          qtot(i,:)  = qtot(i,:)+qtot_e(iel,:)
          neigh = neigh+1
       ENDDO !j

       if(neigh==0.or.eta2(i)+dp1(i)==0) call parallel_abort('sed2d_main: neigh')
       Cdsed(i)   = Cdsed(i)/neigh
       dpdxy(i,1) = dpdxy(i,1)/neigh
       dpdxy(i,2) = dpdxy(i,2)/neigh
       qav(i,:)   = rhosed*qav(i,:)/(neigh*dt*dtsed2d) !(kg/m/s)
       qdt(i,:)   = qdt(i,:)/neigh                     !(m3/m)
       qb(i,:)    = rhosed*qb(i,:)/neigh               !(kg/m/s)
       qs(i,:)    = rhosed*qs(i,:)/neigh               !(kg/m/s)
       qtot(i,:)  = rhosed*qtot(i,:)/neigh             !(kg/m/s)
       cflsed(i) = 3.d0*DSQRT(qdt(i,1)*qdt(i,1)+qdt(i,2)*qdt(i,2))/  &
                   ((eta2(i)+dp1(i))*DSQRT(vc_area(i)))
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
 
  IF(myrank==0) WRITE(16,*)'done computing total sediment load (sed2d)'

#ifdef INCLUDE_TIMING
  timer_ns(8) = timer_ns(8)+mpi_wtime()-wtmp !end of timer 
#endif 

  END SUBROUTINE sed2d_main_node
!--------------------------------------------------------------------

!--------------------------------------------------------------------
  SUBROUTINE sed2d_main_el(it)
!--------------------------------------------------------------------
! Main sed2d routine adapted to element-centered method.
!
! Author: guillaume dodet (gdodet01@univ-lr.fr, gdodet@lnec.pt)
! Date: 29/04/2013
!
! History:
! 06/2013 - G.Dodet: - Added diffusive filter on sediment fluxes
!--------------------------------------------------------------------


  USE elfe_glbl, ONLY : Cdp,dldxy,dp,dt,eta2,idry,idry_s,iplg,isdel,       &
                        isidenode,islg,elnode,np,npa,ns,nsa,pi,out_wwm,  &
                        su2,sv2,rkind
  USE elfe_msgp, ONLY : exchange_s2d,myrank,parallel_abort
  USE sed2d_mod, ONLY : Cdsed,d50,d90,diffac,dpdxy,dpdxy_s,h0_sed,   &
                        idrag,ifilt,imorpho,idsed2d,irough,islope,   &
                        itrans,nskip,poro,qb,qb_s,qfilter,qramp,qs,  &
                        qs_s,qtot,qtot_s,transfac
  USE sed2d_friction
  USE sed2d_transport
  USE sed2d_filter

  IMPLICIT NONE

!- Arguments --------------------------------------------------------
  INTEGER, INTENT(IN) :: it
!- Local variables --------------------------------------------------
  INTEGER :: i,iel,inode,iside,j,k
  REAL(rkind) :: beta,D_star,d50_s,d90_s,hs,htot,kpeak,rampfac,Ucr_c,&
                 Ucr_w,u_star,tau_cr,theta_cr,tp,udir,us,uorb,uorbp, &
                 vs,wlpeak,z0
  REAL(rkind), DIMENSION(npa) :: dp_filt,dp_tmp,neigh2
  REAL(rkind), DIMENSION(nsa) :: Cd_s,tmp_s
!--------------------------------------------------------------------

!--------------------------------------------------------------------
!-                  Compute ramp factor                             -
!--------------------------------------------------------------------
  IF(qramp>0) THEN
    rampfac = TANH(2.d0*it*dt/qramp)
  ELSE
    rampfac = 1.d0
  ENDIF

!- Initialize arrays that are later used to compute nodal outputs
  Cd_s    = 0.d0    !Drag coefficient 
  dpdxy_s = 0.d0    !Bed slope (m/m)
  qb_s    = 0.d0
  qs_s    = 0.d0
  qtot_s  = 0.d0    !Total transport (m3.s-1.m-1)

  DO i = 1,ns

!- Compute variables at side centers --------------------------------
     d50_s  = 0.d0  !Median grain size (m)
     d90_s  = 0.d0  !D90 grain size (m)
     htot   = 0.d0  !Total depth (m)
     uorb   = 0.d0  !Orbital velocity computed by wwm (m/s)
     hs     = 0.d0  !Significant wave height computed by wwm(m)
     tp     = 0.d0  !Peak period computed by wwm (s)
     wlpeak = 0.d0  !Wave length associated to peak period by wwm(m)
     uorbp  = 0.d0  !Peak orbital velocity (m/s)

     DO j = 1,2
        inode = isidenode(j,i)
        iel = isdel(j,i)
!- Total depth 
        htot = htot+(eta2(inode)+dp(inode))/2.d0

!- Velocity 
        us = su2(1,i)
        vs = sv2(1,i)

!- Grain size
        d50_s = d50_s+d50(inode)/2.d0
        d90_s = d90_s+d90(inode)/2.d0

!- bed slope  
        IF(j==2 .AND. isdel(j,i)==0) THEN !boundary
          dpdxy_s(i,1) = 2.d0*dpdxy_s(i,1)
          dpdxy_s(i,2) = 2.d0*dpdxy_s(i,2)
        ELSE
          DO k = 1,3
             dpdxy_s(i,1) = dpdxy_s(i,1)+dp(elnode(k,iel))*dldxy(k,1,iel)/2.d0
             dpdxy_s(i,2) = dpdxy_s(i,2)+dp(elnode(k,iel))*dldxy(k,2,iel)/2.d0
          ENDDO
        ENDIF

!- Wave parameters 
#ifdef USE_WWM
        uorb = uorb+out_wwm(inode,22)/2.d0
        hs = hs+out_wwm(inode,1)/2.d0
        tp = tp+out_wwm(inode,12)/2.d0
        wlpeak = wlpeak+out_wwm(inode,17)/2.d0
#endif

!- Drag coefficient
        IF(idrag==1)THEN
          Cd_s(i) = Cd_s(i)+Cdp(inode)/2.d0
        ENDIF
     ENDDO! j = 1,2

     IF((idry_s(i) == 1).OR.(htot.LE.h0_sed)) CYCLE

!- Compute threshold parameters, roughness and drag coefficient  ----
     CALL compute_thresh(d50_s,d90_s,htot,tp,D_star,Ucr_c,Ucr_w,     &
                         theta_cr,tau_cr)

     z0 = d50_s/12.d0
     IF(idrag.GT.1)CALL compute_drag(d50_s,htot,z0,us,vs,idrag,      &
                                     Cd_s(i))

!- Test if parameters are physically correct ------------------------
#ifdef USE_DEBUG
     IF(htot<0.d0 .OR. DSQRT(us*us+vs*vs)>7.5d0 .OR. d50_s>0.05d0.OR.&
        d90_s>0.10d0 .OR.DSQRT(dpdxy_s(i,1)**2.d0+dpdxy_s(i,2)**2.d0)&
        >5.d0 .OR. uorb>10.d0 .OR. hs>10.d0 .OR. tp>30.d0 .OR.       &
        wlpeak > 1500.d0 .OR. Cd_s(i) > 10.d0 .OR.  z0 > 0.5d0) THEN
        WRITE(12,*)'Warning (1) for it:',it,' and side: ', islg(i)
        WRITE(12,*)'Check following values: '
        WRITE(12,*)'total depth = ',htot
        WRITE(12,*)'current velocity = ',DSQRT(us*us+vs*vs)
        WRITE(12,*)'d50 = ',d50_s
        WRITE(12,*)'d90 = ',d90_s
        WRITE(12,*)'slope = ',DSQRT(dpdxy_s(i,1)**2.d0+dpdxy_s(i,2)**&
                                    2.d0)
        WRITE(12,*)'orbital velocity = ',uorb
        WRITE(12,*)'hs = ',hs
        WRITE(12,*)'tp = ',tp
        WRITE(12,*)'wlpeak = ',wlpeak
        WRITE(12,*)'cd = ', cd_s(i)
        WRITE(12,*)'z0 = ',z0
!        CALL parallel_abort('Sed2d: chech value in outputs/nonfatal')        
     ENDIF
#endif

!- Apply user-defined formula ---------------------------------------

!- Compute bed slope in streamwise direction used for slope effect on 
!  total transport (Soulsby, 1997). Positive if flows runs uphill.
     IF(islope == 1) THEN
       udir = DATAN2(vs,us)
       beta = -(dpdxy_s(i,1)*DCOS(udir)+dpdxy_s(i,2)*DSIN(udir))
       beta = MIN(beta,0.6d0) !May lead negative transport otherwise
     ELSE
       beta = 0.d0
     ENDIF

     SELECT CASE(itrans)
      CASE(1) !Engelund-Hansen (1967)
       CALL eh67(Cd_s(i),d50_s,us,vs,beta,qtot_s(i,:))

      CASE(2) !Ackers and White (1973)
       CALL aw73(Cd_s(i),d50_s,us,vs,htot,beta,D_star,qtot_s(i,:))

      CASE(3) !Soulsby - Van Rijn (Soulsby, 1997)
       CALL svr97_bedl(Cd_s(i),d50_s,us,vs,htot,dpdxy_s(i,:),Ucr_c,  &
                       tau_cr,uorb,qb_s(i,:))
       CALL svr97_susp(Cd_s(i),d50_s,us,vs,D_star,Ucr_c,uorb,        &
                       qs_s(i,:))
       qtot_s(i,:) = (qs_s(i,:)+qb_s(i,:))*(1.d0-1.6d0*beta)

      CASE(4) !Van-Rijn (2007) needs peak orbital velocity
#ifdef USE_WWM
       IF((tp.GT.0.d0) .AND. (wlpeak.GT.0.d0)) THEN
         kpeak = 2.d0*pi/wlpeak
         uorbp = pi*hs/(tp*DSINH(kpeak*htot))
       ENDIF
#endif
       CALL vr07_bedl(Cd_s(i),d50_s,us,vs,htot,dpdxy_s(i,:),uorb,    &
                      Ucr_c,Ucr_w,tau_cr,qb_s(i,:))
       CALL vr07_susp(d50_s,us,vs,uorb,D_star,Ucr_c,Ucr_w,qs_s(i,:))
       qtot_s(i,:) = (qs_s(i,:)+qb_s(i,:))*(1.d0-1.6d0*beta)

      CASE(5) !Camenen and Larson (2008)
!       CALL cl08_bed(...,qb_e(i,:))
!       CALL cl08_susp(...,qs_e(i,:))
!       qtot_s(i,:) = qb_e(i,:)+qs_e(i,:)
     END SELECT

!- Apply diffusion and correction factor if required ----------------
     qtot_s(i,1) = rampfac*transfac*(qtot_s(i,1)+(1.d0-poro)*diffac* &
                   DABS(qtot_s(i,1))*dpdxy_s(i,1))
     qtot_s(i,2) = rampfac*transfac*(qtot_s(i,2)+(1.d0-poro)*diffac* &
                   DABS(qtot_s(i,2))*dpdxy_s(i,2))
  ENDDO !ns
  CALL exchange_s2d(qtot_s(:,1))
  CALL exchange_s2d(qtot_s(:,2))

!- Apply diffusive filter on total transport if required ------------
  IF(qfilter == 1) THEN
    CALL sed2d_filter_diffu(qtot_s(:,1),tmp_s,nsa)
    qtot_s(:,1) = tmp_s
    CALL sed2d_filter_diffu(qtot_s(:,2),tmp_s,nsa)
    qtot_s(:,2) = tmp_s
    IF(myrank==0) WRITE(16,*)'done filtering transport rate (sed2d)'
  ENDIF

  CALL exchange_s2d(qtot_s(:,1))  
  CALL exchange_s2d(qtot_s(:,2))  
  CALL exchange_s2d(dpdxy_s(:,1)) !check if necessary
  CALL exchange_s2d(dpdxy_s(:,2)) !check if necessary
  CALL exchange_s2d(Cd_s(:))      !check if necessary
  CALL exchange_s2d(qb_s(:,1))    !check if necessary
  CALL exchange_s2d(qb_s(:,2))    !check if necessary
  CALL exchange_s2d(qs_s(:,1))    !check if necessary
  CALL exchange_s2d(qs_s(:,2))    !check if necessary

!--------------------------------------------------------------------
!-       Interpolate variables at node for output purpose only      -
!--------------------------------------------------------------------
  Cdsed = 0.d0
  dpdxy(:,:) = 0.d0
  qb = 0.d0
  qs = 0.d0
  qtot = 0.d0
  neigh2 = 0
  DO i = 1,nsa
     DO j = 1,2
        inode = isidenode(j,i)
        IF(idry(inode) == 1) CYCLE
        Cdsed(inode) = Cdsed(inode)+Cd_s(i)
        dpdxy(inode,:) = dpdxy(inode,:)+dpdxy_s(i,:)
        qb(inode,:) = qb(inode,:)+qb_s(i,:)
        qs(inode,:) = qs(inode,:)+qs_s(i,:)
        qtot(inode,:) = qtot(inode,:)+qtot_s(i,:)
        neigh2(inode) = neigh2(inode)+1
     ENDDO !j=1,2
  ENDDO !nsa

  DO i = 1,np
     IF(idry(i) == 1) CYCLE
     if(neigh2(i)==0) call parallel_abort('sed2d_main: neigh2')
     Cdsed(i) = Cdsed(i)/neigh2(i)
     dpdxy(i,:) = dpdxy(i,:)/neigh2(i)
     qb(i,:) = qb(i,:)/neigh2(i)
     qs(i,:) = qs(i,:)/neigh2(i)
     qtot(i,:) = qtot(i,:)/neigh2(i)
  ENDDO !np

!Error: exchange to get ghosts rite

  IF(myrank==0)WRITE(16,*)'done computing total sediment load (sed2d)'

!--------------------------------------------------------------------
!-                    Compute bed change at nodes                   -
!--------------------------------------------------------------------
  IF(imorpho == 1) THEN
    CALL sed2d_morpho(it)
    IF(myrank==0)WRITE(16,*)'done computing bottom evolution (sed2d)'
  ENDIF

!--------------------------------------------------------------------
!-                      Apply user-defined filter                   -
!--------------------------------------------------------------------
  IF(ifilt>0 .AND. MOD(it,nskip)==0) THEN
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
    IF(myrank==0)WRITE(16,*)'done filtering new bathymetry (sed2d)'
  ENDIF

  DO i = 1,npa
     IF(dp(i)/=dp(i)) THEN
       WRITE(12,*)'Warning (2) for it:',it,' and node: ', iplg(i)
       WRITE(12,*)'dp = ',dp(i)
!       CALL parallel_abort('Sed2d: depth is NaN')
     ENDIF
  ENDDO

  END SUBROUTINE sed2d_main_el
!--------------------------------------------------------------------
END SUBROUTINE sed2d_main
