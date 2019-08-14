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

MODULE sed2d_transport
!--------------------------------------------------------------------
! This module contains several subroutines/formulae to compute the 
! sediment transport rate at one location:
!
! - eh67       Total transport from Engelund and Hansen (1967)        
! - aw73       Total transport from Ackers and White (1973)
! - svr97_bedl Bedload transport from Soulsby and Van-Rijn (1997)
! - svr97_susp Suspended transport from Soulsby and Van-Rijn (1997)
! - vr07_bedl  Bedload transport from Soulsby (2007)
! - vr07_susp  Suspended transport from Soulsby (2007)
! - cl11_tot   Total transport (bedload + suspended) from Camenen
!              and Larson (2011)
!
! and
!
! - lesser04_slope Coefficients for slope effect on bed load from 
!                  Lesser et al. (2004)
!                                                                
! Author: guillaume dodet (gdodet01@univ-lr.fr, gdodet@lnec.pt)    
! Date:   20/02/2013
!
! History:                                                
! 03/2013 - G.Dodet: - Decomposed svr97 into a bedload and a 
!                      suspended transport routine;
!                    - Added Van-Rijn (2007);
!                    - Put computation of coefficients for slope
!                      effects on bedload (Lesser et al. 2004) in a 
!                      distinct subroutine;
! 07/2013 - T.Guerin: - Added Camenen and Larson (2011) and wave 
!                       asymmetry computation based on Elfrink (2006)
!--------------------------------------------------------------------

  IMPLICIT NONE

CONTAINS

!--------------------------------------------------------------------
  SUBROUTINE eh67(Cd,d50,u,v,beta,qtot)
!--------------------------------------------------------------------
! This subroutine computes the total sediment transport rate based on 
! Engelund-Hansen (1967) see S97, p.175 
!
! NB: Cd should be determined by alluvial friction method (Engelund,
!     1966, see S97,p175)
!
! Author: guillaume dodet (gdodet01@univ-lr.fr, gdodet@lnec.pt)
! Date: 20/02/2013
!
! History:
! 03/2013 - G.Dodet: Added slope effect on total transport
!--------------------------------------------------------------------

  USE schism_glbl, ONLY : grav,rkind
  USE schism_msgp, ONLY: parallel_abort
  USE sed2d_mod, ONLY : s

  IMPLICIT NONE

!- Arguments --------------------------------------------------------
  REAL(rkind), INTENT(IN) :: beta,Cd,d50,u,v
  REAL(rkind), DIMENSION(2), INTENT(OUT) :: qtot
!- Local variables --------------------------------------------------
  REAL(rkind) :: C,unorm
!--------------------------------------------------------------------

  unorm = DSQRT(u*u+v*v)
  if(Cd<=0) call parallel_abort('eh67: Cd<=0')
  C = DSQRT(grav/Cd)

  qtot(1) = (1.d0-1.6d0*beta)*(0.05d0/((s-1.d0)**2.d0*d50*C**3.d0*   &
            DSQRT(grav)))*(unorm**4.d0)*u
  qtot(2) = (1.d0-1.6d0*beta)*(0.05d0/((s-1.d0)**2.d0*d50*C**3.d0*   &
            DSQRT(grav)))*(unorm**4.d0)*v

  END SUBROUTINE eh67

!--------------------------------------------------------------------
  SUBROUTINE aw73(Cd,d50,u,v,h,beta,D_star,qtot)
!--------------------------------------------------------------------
! This subroutine computes the total sediment transport rate based on
! Ackers and White (1973), S97, p.175
!
! NB: - with this formula, d50 in sed2d.in should be replaced by
!       the equivalent d35
!     - Cd should be determined by White et al.(1980) friction
!       alluvial method. see S97, p.175
!
! Author: guillaume dodet (gdodet01@univ-lr.fr, gdodet@lnec.pt)
! Date: 20/02/2013
!
! History:
! 03/2013 - G.Dodet: Added slope effect on total transport
!--------------------------------------------------------------------

  USE schism_glbl, ONLY : grav,rkind
  USE schism_msgp, ONLY : parallel_abort
  USE sed2d_mod, ONLY : s

  IMPLICIT NONE

!- Arguments --------------------------------------------------------
  REAL(rkind), INTENT(IN) :: beta,Cd,D_star,d50,h,u,v
  REAL(rkind), DIMENSION(2), INTENT(OUT) :: qtot
!- Local variables --------------------------------------------------
  REAL(rkind) :: Aaw,Caw,Faw,lslope,m,n,udir,unorm,u_star
!--------------------------------------------------------------------

  unorm = DSQRT(u*u+v*v)
  udir = DATAN2(v,u)
  qtot = 0.d0

  IF((D_star.GE.1.d0).AND.(D_star.LE.60.d0)) THEN
    n = 1.d0-0.243d0*DLOG(D_star)
    Aaw = 0.23d0/DSQRT(D_star)+0.14d0
    m = 6.83d0/D_star+1.67d0
    Caw = DEXP(2.79d0*DLOG(D_star)-0.426d0*(DLOG(D_star))**2.d0-7.97d0)
  ELSEIF(D_star.GT.60) THEN
    n = 0.d0
    Aaw = 0.17d0
    m = 1.78d0
    Caw = 0.025d0
  ELSE
    CALL parallel_abort('sed2d: D_star out of range in AW73')
  ENDIF

  IF(unorm.NE.0.d0) THEN
    u_star = DSQRT(Cd)*unorm
    Faw = (u_star**n/DSQRT(grav*(s-1.d0)*d50))*(unorm/(2.46d0*DLOG   &
          (10.d0*h/d50)))**(1.d0-n)
    IF(Faw.GE.Aaw) THEN
      qtot(1) = (1.d0-1.6d0*beta)*((Caw*unorm*d50*(unorm/u_star)**n)*&
                ((Faw-Aaw)/Aaw)**m)*DCOS(udir)
      qtot(2) = (1.d0-1.6d0*beta)*((Caw*unorm*d50*(unorm/u_star)**n)*&
                ((Faw-Aaw)/Aaw)**m)*DSIN(udir)
    ENDIF
  ENDIF

  END SUBROUTINE aw73

!--------------------------------------------------------------------
  SUBROUTINE svr97_bedl(Cd,d50,u,v,h,dpdxy,U_cr,tau_cr,urms,qb)
!--------------------------------------------------------------------
! This subroutine computes the bed-load transport rate based on 
! Soulsby - Van Rijn (S97,p. 183)
!
! NB: This formula was computed for rippled bed. A value of z0 = 6 mm
!     should be used (S97, p.184)
! Cd = (0.4d0/(1.d0+DLOG(0.006d0/htot)))**2.d0
!
! Author: guillaume dodet (gdodet01@univ-lr.fr, gdodet@lnec.pt)
! Date: 11/03/2013
!
! History:
!--------------------------------------------------------------------

  USE schism_glbl, ONLY : grav,rkind
  USE schism_msgp, ONLY : parallel_abort
  USE sed2d_mod, ONLY : s

  IMPLICIT NONE

!- Arguments --------------------------------------------------------
  REAL(rkind), INTENT(IN) :: Cd,d50,h,tau_cr,u,U_cr,urms,v
  REAL(rkind), DIMENSION(2), INTENT(IN) :: dpdxy
  REAL(rkind), DIMENSION(2), INTENT(OUT) :: qb
!- Local variables --------------------------------------------------
  REAL(rkind) :: alpha_n,alpha_s,Asb,aux,cff,udir,unorm
!--------------------------------------------------------------------

!- Slope effect on bedload transport as in Lesser et al (2004)
  CALL lesser04_slope(Cd,d50,u,v,dpdxy,tau_cr,alpha_s,alpha_n)

!- Transport rate
  if(Cd<=0) call parallel_abort('svr97_bedl: Cd<=0')
  qb = 0.d0
  unorm = SQRT(u*u+v*v)
  udir = ATAN2(v,u)
  Asb = (0.005d0*h*(d50/h)**1.2d0)/(((s-1.d0)*grav*d50)**1.2d0)
  aux = SQRT(unorm**2.d0+(0.018d0/Cd)*urms**2.d0)
  IF(((aux-U_cr).GE.0.d0)) THEN
    cff = Asb*unorm*(aux-U_cr)**2.4d0
    qb(1) = alpha_s*cff*COS(udir)-alpha_n*alpha_s*cff*SIN(udir)
    qb(2) = alpha_s*cff*SIN(udir)-alpha_n*alpha_s*cff*COS(udir)
  ENDIF

  END SUBROUTINE svr97_bedl

!--------------------------------------------------------------------
  SUBROUTINE svr97_susp(Cd,d50,u,v,D_star,U_cr,urms,qs)
!--------------------------------------------------------------------
! This subroutine computes the suspended and total sediment transport
! rate based on Soulsby - Van Rijn (S97,p. 183)
!
! NB: This formula was computed for rippled bed. A value of z0 = 6 mm
!     should be used (S97, p.184)
! Cd = (0.4d0/(1.d0+DLOG(0.006d0/htot)))**2.d0
!
! Author: guillaume dodet (gdodet01@univ-lr.fr, gdodet@lnec.pt)
! Date: 11/03/2013
!
! History:
!--------------------------------------------------------------------

  USE schism_glbl, ONLY : grav,rkind
  USE schism_msgp, ONLY : parallel_abort
  USE sed2d_mod, ONLY : s

  IMPLICIT NONE

!- Arguments --------------------------------------------------------
  REAL(rkind), INTENT(IN) :: Cd,d50,u,v,D_star,U_cr,urms
  REAL(rkind), DIMENSION(2), INTENT(OUT) :: qs
!- Local variables --------------------------------------------------
  REAL(rkind) :: Ass,aux,udir,unorm
!--------------------------------------------------------------------

  if(Cd<=0.or.D_star==0) call parallel_abort('svr97_susp: ')
  qs = 0.d0
  unorm = SQRT(u*u+v*v)
  udir = ATAN2(v,u)
  Ass = (0.012d0*d50*D_star**(-0.6d0))/(((s-1.d0)*grav*d50)**1.2d0)
  aux = SQRT(unorm**2.d0+(0.018d0/Cd)*urms**2.d0)

!- Transport rate
  IF(((aux-U_cr).GE.0.d0)) THEN
    qs(1) = (Ass*unorm*(aux-U_cr)**2.4d0)*COS(udir)
    qs(2) = (Ass*unorm*(aux-U_cr)**2.4d0)*SIN(udir)
  ENDIF

  END SUBROUTINE svr97_susp

!--------------------------------------------------------------------
  SUBROUTINE vr07_bedl(Cd,d50,u,v,h,dpdxy,uorb,Ucr_c,Ucr_w,tau_cr,qb)
!--------------------------------------------------------------------
! This subroutine computes the bed-load transport rate based on
! Van Rijn (2007)
!
! Author: guillaume dodet (gdodet01@univ-lr.fr, gdodet@lnec.pt)
! Date: 11/03/2013
!
! History:
! 04/2013 - G.Dodet: Added a maximum velocity value based on Van Rijn
!                    (2007) to compute qb.
! 09/2014 - T.Guerin: Added condition for unorm+uorb=0 case
!--------------------------------------------------------------------

  USE schism_glbl, ONLY : grav,rkind
  USE schism_msgp, ONLY : parallel_abort
  USE sed2d_mod, ONLY : s  

  IMPLICIT NONE

!- Arguments --------------------------------------------------------
  REAL(rkind), INTENT(IN) :: Cd,d50,h,tau_cr,u,Ucr_c,Ucr_w,uorb,v
  REAL(rkind), DIMENSION(2), INTENT(IN) :: dpdxy
  REAL(rkind), DIMENSION(2), INTENT(OUT) :: qb
!- Local variables --------------------------------------------------
  REAL(rkind) :: alpha_n,alpha_s,beta,cff,k,Me,Ucr,ueff,udir,unorm  
!- User-defined parameters ------------------------------------------
  REAL(rkind), PARAMETER :: gama = 0.4d0       !0.8 for regular waves 
!--------------------------------------------------------------------

!- Slope effect on bedload transport as in Lesser et al (2004)
  CALL lesser04_slope(Cd,d50,u,v,dpdxy,tau_cr,alpha_s,alpha_n)

!- Transport rate
  qb = 0.d0
  unorm = MIN(SQRT(u*u+v*v),1.8d0) ! Range of validity for the formula
  udir = ATAN2(v,u)
  ueff = unorm+gama*uorb

  IF (unorm+uorb==0.d0) THEN
    qb(:) = 0.d0
  ELSE
    beta = unorm/(unorm+uorb)
    Ucr = beta*Ucr_c+(1.d0-beta)*Ucr_w
    Me = (ueff-Ucr)/SQRT((s-1.d0)*grav*d50)

    IF(ueff-Ucr.GE.0.d0) THEN
      cff = 0.015d0*unorm*h*((d50/h)**1.2d0)*(Me**1.5d0)
      qb(1) = alpha_s*cff*COS(udir)-alpha_n*alpha_s*cff*SIN(udir)
      qb(2) = alpha_s*cff*SIN(udir)-alpha_n*alpha_s*cff*COS(udir)
    ENDIF
  ENDIF

  END SUBROUTINE vr07_bedl

!--------------------------------------------------------------------
  SUBROUTINE vr07_susp(d50,u,v,uorb,D_star,Ucr_c,Ucr_w,qs)
!--------------------------------------------------------------------
! This subroutine computes the suspended transport rate based on
! Van Rijn (2007)
!
! Author: guillaume dodet (gdodet01@univ-lr.fr, gdodet@lnec.pt)
! Date: 11/03/2013
!
! History:
! 04/2013 - G.Dodet: Added a maximum velocity value based on Van Rijn
!                    (2007) to compute qs.
! 09/2014 - T.Guerin: Added condition for unorm+uorb=0 case
!--------------------------------------------------------------------

  USE schism_glbl, ONLY : grav,rkind
  USE schism_msgp, ONLY : parallel_abort
  USE sed2d_mod, ONLY : s

  IMPLICIT NONE

!- Arguments --------------------------------------------------------
  REAL(rkind), INTENT(IN) :: D_star,d50,u,Ucr_c,Ucr_w,uorb,v
  REAL(rkind), DIMENSION(2), INTENT(OUT) :: qs
!- Local variables --------------------------------------------------
  REAL(rkind) :: beta,cff,k,Me,Ucr,udir,ueff,unorm
!- User-defined parameters ------------------------------------------
  REAL(rkind), PARAMETER :: gama = 0.4d0       !0.8 for regular waves
!--------------------------------------------------------------------

!- Transport rate
  qs = 0.d0
  unorm = MIN(SQRT(u*u+v*v),1.8d0) !Range of validity for the formula
  udir = DATAN2(v,u)
  ueff = unorm+gama*uorb

  IF (unorm+uorb==0.d0) THEN
    qs(:) = 0.d0
  ELSE
    beta = unorm/(unorm+uorb)
    Ucr = beta*Ucr_c+(1.d0-beta)*Ucr_w
    Me = (ueff-Ucr)/SQRT((s-1.d0)*grav*d50)

    IF(ueff-Ucr.GE.0.d0) THEN
      if(D_star<=0) call parallel_abort('vr07_susp: D_star<=0')
      cff = 0.012d0*unorm*d50*Me**2.4d0*D_star**(-0.6d0)
      qs(1) = cff*COS(udir)
      qs(2) = cff*SIN(udir)
    ENDIF
  ENDIF

  END SUBROUTINE vr07_susp

!--------------------------------------------------------------------
  SUBROUTINE cl11_tot(d50,D_star,dpdxy,hs,htot,theta_cr,tp,u,v,wdir, &
                     &wlpeak,qb,qs)
!--------------------------------------------------------------------
! This subroutine computes total sediment transport rate (m^2/s) 
! (bedload + suspended load) from Camenen and Larson (2011)
!
! Author: thomas guerin (thomas.guerin@univ-lr.fr)    
! Date: 24/04/2013
!
! History:
!--------------------------------------------------------------------

  USE schism_glbl, ONLY: grav,pi,rho0,rkind
  USE schism_msgp, ONLY : parallel_abort
  USE sed2d_mod, ONLY: iasym,s,wvisco

#ifdef USE_WWM
  USE DATAPOOL, ONLY: gammab =>BRHD
#endif

  IMPLICIT NONE

!- Arguments --------------------------------------------------------  
  REAL(rkind), INTENT(IN) :: d50,D_star,hs,htot,theta_cr,tp,u,v,wdir,&
                             wlpeak
  REAL(rkind), DIMENSION(2), INTENT(IN) :: dpdxy
  REAL(rkind), DIMENSION(2), INTENT(OUT) :: qb,qs
!- User-defined parameters ------------------------------------------  
  REAL(rkind), PARAMETER :: ksd_coef = 2.5d0 ! = 2 or 2.5
!- Constants --------------------------------------------------------  
  REAL(rkind), PARAMETER :: kappa = 0.4d0
  REAL(rkind), PARAMETER :: an = 12.d0
  REAL(rkind), PARAMETER :: b = 4.5d0
  REAL(rkind), PARAMETER :: pwr_cr = 100.d0
  REAL(rkind), PARAMETER :: Ad = 0.7d0
  REAL(rkind), PARAMETER :: d_ks_c = 0.0001d0
  INTEGER, PARAMETER :: ech = 50 !sampling of wave orbital velocity uw(t)
!- Local variables --------------------------------------------------
  REAL(rkind) :: a_b,A_epsi,A_w,AcR,alpha_off,alpha_on,alpha_plb,    &
                 alpha_pls,aux,aw,cR,D_b,D_c,D_w,delta_w,dpdxy_tot,  &
                 dpdxy_wdir,Dpl_off,Dpl_on,epsi,fc,fcw,fw,Hr,Hr_c,   &
                 Hr_w,Hrms,irribaren,k_b,k_c,k_cw_c,k_cw_w,          &
                 kpeak,ksd,ks_c,ks_w,ksf_c,ksf_w,kss_w,Lr_c,Lr_w,    &
                 omega,omega_off,omega_on,phi,psi_w,pwr,qb_n,qb_w,   &
                 qs_n,qs_w,R,rw,sigma_c,sigma_cw,sigma_w,slope_dir,  &
                 slope_wdir_angle,T_crest,T_trough,tau_c,tau_w,      &
                 theta_c,theta_cn,theta_cw,theta_cw_on,theta_cw_off, &
                 theta_cwm,theta_cx,theta_cy,theta_net,theta_w,      &
                 theta_wm,U_crest,U_trough,Uc,Uc_n,Uc_net,Uc_w,      &
                 Uc_x,Uc_y,Ucw_off,Ucw_on,Ucw2_off,Ucw2_on,udir,Uorb,&
                 ustar_c,ustar_w,Uw_crsf,v_off,v_on,wl0,Ws,Xt,Xv,Y, &
                 tmp,tmp2
  REAL(rkind), DIMENSION(ech-1) :: Uorb_Elfrink
  INTEGER :: i
!--------------------------------------------------------------------

  !- constant
  aux = 2.d0*(s-1.d0)*grav*d50 !>0

  !- settling velocity (Soulsby 1997)
  Ws = wvisco*((10.36d0**2+1.049d0*D_star**3)**0.5d0 - 10.36d0) /d50 !>0

  !- current amplitude and direction
  Uc = DSQRT(u*u+v*v)
  udir = DATAN2(v,u)

  IF (hs>0.AND.tp>0.AND.wlpeak>0) THEN
    !- wave angular frequency
    omega = 2.d0*pi/tp !>0

    !- peak wave number
    kpeak = 2.d0*pi/wlpeak !>0

    !- angle between current and waves directions
    phi = udir - wdir

    !- current velocity components in the wave frame
    Uc_w = Uc*DCOS(phi)
    Uc_n = Uc*DSIN(phi)
  ELSE
    Uc_x = Uc*DCOS(udir)
    Uc_y = Uc*DSIN(udir)
  ENDIF

  IF (hs>0.AND.tp>0.AND.wlpeak>0) THEN
    IF (iasym==0) THEN
      !- Airy waves
      tmp=DSINH(kpeak*htot)
      if(tmp==0) call parallel_abort('cl11_tot: tmp=0')
      Uorb = pi*hs/tp/tmp !>0
      U_crest = Uorb !>0
      T_crest = tp/2.d0 !>0
      T_trough = T_crest !>0

    ELSEIF (iasym==1) THEN
      !- Wave asymmetry treatment : compute wave orbital velocity
      !  (horizontal component) from Elfrink (2006)
      CALL wave_asymmetry_Elfrink(hs,tp,wdir,htot,-dpdxy(1),-dpdxy(2),&
           ech,Uorb_Elfrink,U_crest,U_trough,T_crest,T_trough)

      IF (U_crest<=0.OR.U_trough<=0.OR.T_crest<=0.OR.T_trough<=0) THEN
        !- Airy waves
        tmp=DSINH(kpeak*htot)
        if(tmp==0) call parallel_abort('cl11_tot: tmp=0')
        Uorb = pi*hs/tp/tmp !>0
        U_crest = Uorb !>0
        T_crest = tp/2.d0 !>0
        T_trough = T_crest !>0
      ELSE
        Uorb = (U_crest+U_trough)/2.d0 !>0
      ENDIF

    ENDIF !iasym
  ENDIF

  !- total bed roughness and friction factors -----------------------

  !- grain-related roughness
  ksd = ksd_coef*d50

  !--- total bed roughness for current alone
  Lr_c = 1000.d0*d50              !ripple length
  Hr_c = Lr_c/7.d0                !ripple height
  ksf_c = 7.5d0*Hr_c*Hr_c/Lr_c   !form-drag roughness

  !- iterative method
  if(htot<=0) call parallel_abort('cl11_tot: htot<=0')
  ks_c = d_ks_c
  DO WHILE ((-ks_c+ksd+ksf_c + 5.d0*kappa*kappa*Uc*Uc/               &
            &(grav*(s-1.d0)*(1.d0+DLOG(ks_c/30.d0/htot))**2.d0))      &
            &.GT. 0.d0)
    ks_c = ks_c + d_ks_c !total bed roughness for current alone       
    IF (ks_c .GE. (0.1d0*htot)) THEN
      ks_c = 0.1d0*htot !maximum value of ks_c is set to htot/10
      EXIT
    ENDIF
  ENDDO

  !- friction factor for current alone
  fc = 2.d0*(kappa/(1.d0+DLOG(ks_c/30.d0/htot)))**2.d0 !>0

  IF (hs>0.AND.tp>0.AND.wlpeak>0) THEN
    !--- total bed roughness for waves alone
    Aw = Uorb*tp/2.d0/pi
    psi_w = 2.d0*Uorb*Uorb/aux
    IF (psi_w .LT. 10.d0) THEN
      Hr_w = 0.22d0*Aw
      Lr_w = 1.25d0*Aw
    ELSEIF ((psi_w .GE. 10.d0) .AND. (psi_w .LT. 250.d0)) THEN
      Hr_w = 2.8d0*1.0e-13*(250d0-psi_w)**5d0*Aw
      Lr_w = 1.4d0*1.0e-6*(250d0-psi_w)**2.5d0*Aw
    ELSEIF (psi_w .GE. 250.d0) THEN
      Hr_w = 0.d0
      Lr_w = 0.d0
    ENDIF
    ksf_w = 7.5d0*Hr_w*Hr_w/(Lr_w + 1.0e-6)   !form-drag roughness

    !- calculation of friction factor for waves (Swart 1974)
    R = Aw/ksd
    fw = DEXP(5.21d0*R**(-0.19d0) - 6.d0)
    kss_w = 2.5d0*fw*Uorb*Uorb/grav   !sediment-related roughness
    ks_w = ksd + ksf_w + kss_w
    IF ((Aw/ks_w) .LE. 1.57d0) THEN
      fw = 0.3d0
    ENDIF

    !--- friction factor for waves and current
    Xv = Uc/(Uc+Uorb+1.0e-6)
    fcw = Xv*fc + (1.d0-Xv)*fw

    !- bed shear-stress for waves only
    tau_w = 0.5d0*rho0*fw*Uorb*Uorb !>=0
  ENDIF

  !- bed shear-stress for current only
  tau_c = 0.5d0*rho0*fc*Uc*Uc

  !--- Shields number
  !- current Shield nb
  theta_c = tau_c / (rho0*0.5d0*aux) !>=0

  IF (hs>0.AND.tp>0.AND.wlpeak>0) THEN
    !- current Shield nb relative to normal current component in the wave frame
    theta_cn = fc*Uc_n*Uc_n / aux !>=0

    !- maximum and mean wave Shield nb
    theta_w = tau_w / (rho0*0.5d0*aux) !>=0
    theta_wm = theta_w/2.d0   !valid for a sinusoidal wave profile (!)

    !- maximum and mean wave Shield nb due to wave-current interaction
    tmp=theta_c*theta_c + theta_w*theta_w+ 2.d0*theta_w*theta_c*DCOS(phi)
    tmp2=theta_c*theta_c + theta_wm*theta_wm+ 2.d0*theta_wm*theta_c*DCOS(phi)
    if(tmp<=0.or.tmp2<0) call parallel_abort('cl11_tot: tmp<=0.or.tmp2<0')

    theta_cw=sqrt(tmp)
    theta_cwm=sqrt(tmp2)
!  theta_cw = DSQRT(theta_c*theta_c + theta_w*theta_w                 &
!                   + 2.d0*theta_w*theta_c*DCOS(phi))
!  theta_cwm = DSQRT(theta_c*theta_c + theta_wm*theta_wm              &
!                    + 2.d0*theta_wm*theta_c*DCOS(phi))
  ELSE
    theta_cx = fc*Uc_x*Uc_x / aux !>=0
    theta_cy = fc*Uc_y*Uc_y / aux !>=0
  ENDIF

  !--- Bed load transport -------------------------------------------

  IF (hs>0.AND.tp>0.AND.wlpeak>0) THEN
    !- mean values of the instantaneous Shield nb over the wave half periods
    if(omega==0.or.T_trough<=0.or.T_crest<=0) then
       write(12,*) 'omega=',omega
       write(12,*) 'T_trough=',T_trough
       write(12,*) 'T_crest=',T_crest
       call parallel_abort('cl11_tot: omega==0')
    endif
    Ucw2_on = (Uc_w*Uc_w*T_crest                                       &
               + 2.d0/omega*Uc_w*Uorb*(1.d0-DCOS(omega*T_crest))       &
               + Uorb*Uorb/2.d0*(T_crest-DSIN(2.d0*omega*T_crest)      &
               /2.d0/omega)) /T_crest
    Ucw2_off = (Uc_w*Uc_w*T_trough + 2.d0/omega*Uc_w*Uorb              &
                *(DCOS(omega*T_crest)-DCOS(omega*tp))                  &
                + Uorb*Uorb/2.d0*(T_trough+(DSIN(2.d0*omega*T_crest)   &
                -DSIN(2.d0*omega*tp))/2.d0/omega)) /T_trough
    theta_cw_on = fcw*Ucw2_on / aux
    theta_cw_off = -fcw*Ucw2_off / aux

    !- calculation of the phase-lag effects coefficient (Camenen and Larson 2006)
    delta_w = DSQRT(wvisco*tp/pi)
    rw = U_crest/Uorb - 1.d0   !wave asymmetry parameter
    Uw_crsf = 8.35d0*DSQRT((s-1.d0)*grav*DSQRT(d50*delta_w))*(1.d0+rw)
    if(Ucw2_on<=0.or.Ucw2_off<=0) call parallel_abort('cl11_tot: Ucw2_on<0')
    Ucw_on = DSQRT(Ucw2_on)
    Ucw_off = DSQRT(Ucw2_off)
    alpha_on = wvisco**0.25d0*Ucw_on**0.5d0                            &
               *DEXP(-(Uw_crsf/Ucw_on)**2.d0) /(Ws*T_crest**0.75d0)
    alpha_off = wvisco**0.25d0*Ucw_off**0.5d0                          &
                *DEXP(-(Uw_crsf/Ucw_off)**2.d0) /(Ws*T_trough**0.75d0)
    alpha_plb = alpha_on-alpha_off

    !- net sediment transporting Shields nb
    theta_net = (1.d0-alpha_plb)*theta_cw_on                           &
                + (1.d0+alpha_plb)*theta_cw_off

    !--- bed load transport components 
    !- parallel and normal to the wave direction

    if(theta_c+theta_w==0) call parallel_abort('cl11_tot: theta_c+theta_w==0')
    Xt = theta_c/(theta_c+theta_w)
    aw = 6.d0+6.d0*Xt

    IF (theta_net .GE. 0.d0) THEN
        qb_w = DSQRT((s-1.d0)*grav*d50*d50*d50)*aw*DSQRT(theta_net)    &
               *theta_cwm*DEXP(-b*theta_cr/theta_cw)
    ELSE
        qb_w = -DSQRT((s-1.d0)*grav*d50*d50*d50)*aw*DSQRT(-theta_net)  &
               *theta_cwm*DEXP(-b*theta_cr/theta_cw)
    ENDIF
    IF (Uc_n .GE. 0.d0) THEN
        qb_n = DSQRT((s-1.d0)*grav*d50*d50*d50)*an*DSQRT(theta_cn)     &
               *theta_cwm*DEXP(-b*theta_cr/theta_cw)
    ELSE
        qb_n = -DSQRT((s-1.d0)*grav*d50*d50*d50)*an*DSQRT(theta_cn)    &
               *theta_cwm*DEXP(-b*theta_cr/theta_cw)
    ENDIF

    !- x,y components
    qb(1) = qb_w*DCOS(wdir) - qb_n*DSIN(wdir)
    qb(2) = qb_w*DSIN(wdir) + qb_n*DCOS(wdir)

  ELSE
    aw = 12.d0

    IF (theta_c .GT. 0.d0) THEN
      IF (Uc_x .GE. 0.d0) THEN
        qb(1) = DSQRT((s-1.d0)*grav*d50*d50*d50)*aw*DSQRT(theta_cx)    &
                 *theta_c*DEXP(-b*theta_cr/theta_c)
      ELSE
        qb(1) = -DSQRT((s-1.d0)*grav*d50*d50*d50)*aw*DSQRT(theta_cx)   &
                 *theta_c*DEXP(-b*theta_cr/theta_c)
      ENDIF
      IF (Uc_y .GE. 0.d0) THEN
        qb(2) = DSQRT((s-1.d0)*grav*d50*d50*d50)*an*DSQRT(theta_cy)    &
                 *theta_c*DEXP(-b*theta_cr/theta_c)
      ELSE
        qb(2) = -DSQRT((s-1.d0)*grav*d50*d50*d50)*an*DSQRT(theta_cy)   &
                 *theta_c*DEXP(-b*theta_cr/theta_c)
      ENDIF
    ENDIF
  ENDIF

  !- put zero if NaN value for qb (temporary solution)
!  IF ((qb(1) .NE. qb(1)) .OR. (qb(2) .NE. qb(2))) THEN
!    qb(:) = 0.d0
!  ENDIF

  !--- Suspended load transport -------------------------------------

  IF (hs>0.AND.tp>0.AND.wlpeak>0) THEN
    !- net mean current after a wave period
    Uc_net = Ucw_on - Ucw_off
  ENDIF

  !- bed reference sediment concentration
  AcR = 0.0035d0*DEXP(-0.3d0*D_star)
  IF (hs>0.AND.tp>0.AND.wlpeak>0) THEN
    IF (theta_cw .GT. 0.d0) THEN
      cR = AcR*theta_cwm*DEXP(-4.5d0*theta_cr/theta_cw)
    ELSE
      cR = 0.d0
    ENDIF
  ELSE
    IF (theta_c .GT. 0.d0) THEN
      cR = AcR*theta_c*DEXP(-4.5d0*theta_cr/theta_c)
    ELSE
      cR = 0.d0
    ENDIF
  ENDIF

  !--- calcul of sediment diffusivity

  !-- sediment diffusivity due to current alone
  ustar_c = DSQRT(tau_c/rho0)
  D_c = tau_c*ustar_c  !energy dissipation from bottom friction due to current
  !- Schmidt nb calculation
  IF (ustar_c.GT.0) THEN
    IF ((Ws/ustar_c) .LE. 1.d0) THEN
      sigma_c = 0.4d0 + 3.5d0*(DSIN(pi*Ws/2.d0/ustar_c))**2.d0
    ELSE
      sigma_c = 1.d0 + 2.9d0*(DSIN(pi*ustar_c/2.d0/Ws))**2.d0
      ! sigma_c = 1.d0 + 2.9d0*(DSIN(pi*Ws/2.d0/ustar_c))**2.d0
      ! sigma_c = 1.d0 + 2.9d0*(DSIN(pi*ustar_c/2.d0/Ws))**2.5d0
    ENDIF
    k_c = kappa*sigma_c/6.d0
  ELSE
    k_c = 0.d0
  ENDIF

  IF (hs>0.AND.tp>0.AND.wlpeak>0) THEN
    !-- sediment diffusivity due to waves alone
    ustar_w = DSQRT(tau_w/rho0)
    D_w = tau_w*ustar_w   !energy dissipation from bottom friction due to waves                    
    !- Schmidt nb calculation
    IF ((Ws/ustar_w) .LE. 1.d0) THEN
      if(ustar_w==0) call parallel_abort('cl11_tot: ustar_w==0')
      sigma_w = 0.15d0 + 1.5d0*(DSIN(pi*Ws/2.d0/ustar_w))**2.d0
    ELSE
      sigma_w = 1.d0 + 0.65d0*(DSIN(pi*ustar_w/2.d0/Ws))**2.d0
      ! sigma_w = 1.d0 + 0.65d0*(DSIN(pi*Ws/2.d0/ustar_w))**2.d0
      ! sigma_w = 1.d0 + 0.65d0*(DSIN(pi*ustar_w/2.d0/Ws))**2.5d0
    ENDIF

    !- sediment diffusivity due to waves and current
    if(theta_c+theta_w==0) call parallel_abort('cl11_tot: theta_c+theta_w==0 (2)')
    Y = theta_c/(theta_c+theta_w)
    sigma_cw = Y*sigma_c + (1.d0-Y)*sigma_w
    k_cw_c = kappa*sigma_cw/6.d0
    k_cw_w = kappa*sigma_cw/3.d0/pi

    !-- sediment diffusivity due to wave breaking
    !- calculation of the slope along the local wave direction
    dpdxy_tot = DSQRT(dpdxy(1)*dpdxy(1) + dpdxy(2)*dpdxy(2))
    slope_dir = DATAN2(-dpdxy(2),-dpdxy(1))
    slope_wdir_angle = DABS(wdir-slope_dir)
    IF (slope_wdir_angle .GT. pi) THEN
      slope_wdir_angle = 2.d0*pi - slope_wdir_angle
    ENDIF
    dpdxy_wdir = DCOS(slope_wdir_angle)*dpdxy_tot

    !- irribaren nb
    if(hs<0.or.wlpeak<=0) call parallel_abort('cl11_tot: hs<0')
    IF (dpdxy_wdir .GT. 0.d0) THEN
      irribaren = dpdxy_wdir/DSQRT(hs/wlpeak)
    ELSE
      irribaren = 0.001d0/DSQRT(hs/wlpeak)
    ENDIF

    !- calculation of energy dissipation from wave breaking
    Hrms = hs/DSQRT(2.d0)
    if(Hrms==0) call parallel_abort('cl11_tot: Hrms=0')
#ifdef USE_WWM
    a_b = DEXP(-(gammab*htot/Hrms)**2.d0)
#else
    a_b = 1.d0
#endif
    A_epsi = 2.d0*DTANH(5.d0*irribaren)
    tmp=tp*(4.d0*htot*htot-hs*hs)
    if(tmp==0) call parallel_abort('cl11_tot: tmp=0 (2)')
    D_b = a_b*A_epsi*rho0*grav*htot*hs**3.d0/tmp
!        / (tp*(4.d0*htot*htot-hs*hs))

    !- parametrization of the efficiency parameter for wave breaking (see Camenen 2007, p. 132)
    k_b = 0.01d0
    ! k_b = 0.012d0*(1.d0+DTANH(irribaren))
    ! k_b = 0.062d0*(1.d0-0.9d0*DTANH(0.25d0*ustar_w/Ws))

    !-- total sediment diffusivity
    epsi = htot*((D_b*k_b**3.d0 + D_c*k_cw_c**3.d0 + D_w*k_cw_w**3.d0) &
           /rho0)**(1.d0/3.d0)

    !--- suspended load transport components
    !- parallel and normal to the wave direction
    qs_w = Uc_net*cR*epsi/Ws * (1.d0-DEXP(-Ws*htot/epsi))
    qs_n = Uc_n*cR*epsi/Ws * (1.d0-DEXP(-Ws*htot/epsi))

    !- x,y components
    qs(1) = qs_w*DCOS(wdir) - qs_n*DSIN(wdir)
    qs(2) = qs_w*DSIN(wdir) + qs_n*DCOS(wdir)

  ELSE
    !-- total sediment diffusivity
    epsi = htot * ((D_c*k_c**3.d0)/rho0)**(1.d0/3.d0)

    !--- suspended load transport components
    IF (epsi.GT.0.d0) THEN
      qs(1) = Uc_x*cR*epsi/Ws * (1.d0-DEXP(-Ws*htot/epsi))
      qs(2) = Uc_y*cR*epsi/Ws * (1.d0-DEXP(-Ws*htot/epsi))
    ELSE
      qs(1) = 0.d0
      qs(2) = 0.d0
    ENDIF
  ENDIF

  !- put zero if NaN value for qs (temporary solution)
!  IF ((qs(1) .NE. qs(1)) .OR. (qs(2) .NE. qs(2))) THEN
!    qs(:) = 0.d0
!  ENDIF

  END SUBROUTINE cl11_tot

!--------------------------------------------------------------------
  SUBROUTINE lesser04_slope(Cd,d50,u,v,dpdxy,tau_cr,alpha_s,alpha_n)
!--------------------------------------------------------------------
! This subroutine computes the coefficients of Lesser et al. (2004) 
! to take into account the slope effect on bed-load transport
!
! Author: guillaume dodet (gdodet01@univ-lr.fr, gdodet@lnec.pt)
! Date: 11/03/2013
!
! History:
!--------------------------------------------------------------------

  USE schism_glbl, ONLY : pi,rho0,rkind
  USE schism_msgp, ONLY : parallel_abort
  USE sed2d_mod, ONLY : islope

  IMPLICIT NONE

!- Arguments --------------------------------------------------------
  REAL(rkind), INTENT(IN) :: Cd,d50,tau_cr,u,v
  REAL(rkind), DIMENSION(2), INTENT(IN) :: dpdxy
  REAL(rkind), INTENT(OUT) :: alpha_n,alpha_s
!- Local variables --------------------------------------------------
  REAL(rkind) :: cff,lslope,sed_angle,tau,tslope,udir,unorm
!--------------------------------------------------------------------

  unorm = SQRT(u*u+v*v)
  udir = ATAN2(v,u)
  lslope = dpdxy(1)*COS(udir)+dpdxy(2)*SIN(udir) !longitudinal slope
  tslope = dpdxy(2)*COS(udir)-dpdxy(1)*SIN(udir) !transverse slope

  IF(islope == 2) THEN
    sed_angle = TAN(30.d0*pi/180.d0)
    cff = MIN(ABS(lslope),0.9d0*sed_angle)*SIGN(1.d0,lslope)
    alpha_s = 1.d0+1.d0*((sed_angle/(COS(ATAN(cff))*(sed_angle-cff)))-1.d0)
    tau = rho0*Cd*unorm**2.d0
    IF(tau == 0.d0) THEN
       alpha_n = 0.d0
    ELSE
       alpha_n = 1.5d0*SQRT(tau_cr/tau)*tslope
    ENDIF
  ELSE
    alpha_n = 0.d0
    alpha_s = 1.d0
  ENDIF

  END SUBROUTINE lesser04_slope

!--------------------------------------------------------------------
  SUBROUTINE wave_asymmetry_Elfrink(H,wave_per,w_dir,depth,dhxi,dhyi,&
                     ech,Uorbi,Ucrest,Utrough,T_crest,T_trough)
!--------------------------------------------------------------------
! This subroutine computes wave asymmetry based on Elfrink et al. 
! (2006, Coastal Engineering)
!
! NB: offshore wave height and wave period are taken into account to
! compute irribaren number and offshore wave length
!
! Author: thomas guerin (thomas.guerin@univ-lr.fr)    
! Date: 26/04/2013
!--------------------------------------------------------------------
  
  USE schism_msgp, ONLY : parallel_abort
  IMPLICIT NONE

!- Arguments --------------------------------------------------------  
  REAL(8), INTENT(IN) :: H,wave_per,w_dir,depth,dhxi,dhyi
  INTEGER, INTENT(IN) :: ech
  REAL(8), INTENT(OUT) :: Ucrest,Utrough,T_crest,T_trough
  REAL(8), DIMENSION(ech-1), INTENT(OUT) :: Uorbi
!- Constants --------------------------------------------------------  
  REAL(8), PARAMETER :: g = 9.80665d0
  REAL(8), PARAMETER :: pi = 3.141592653589793d0 !DACOS(-1.d0)
  REAL(8), PARAMETER :: a1 = 0.38989d0
  REAL(8), PARAMETER :: a2 = -0.0145d0
  REAL(8), PARAMETER :: a3 = -0.0005d0
  REAL(8), PARAMETER :: a4 = 0.5028d0
  REAL(8), PARAMETER :: a5 = 0.9209d0
  REAL(8), PARAMETER :: b1 = 0.5366d0
  REAL(8), PARAMETER :: b2 = 1.16d0
  REAL(8), PARAMETER :: b3 = -0.2615d0
  REAL(8), PARAMETER :: b4 = 0.0958d0
  REAL(8), PARAMETER :: b5 = -0.5623d0
!- Local variables --------------------------------------------------
  REAL(8) :: a1bis,C1,C2,C3,C4,C5,D1,D2,D3,D4,D5,E1,E2,E3,F1,F2,F3,  &
             F4,G1,G2,G3,G4,G5,G6,G7,G8,Hadim,kh,L0,Ladim,P1,P2,P3,  &
             P4,P5,psi,Slope,T0,T1,T2,tv,U0,U1,U2,Uairy,uorb,Ur,     &
             Ustar,Zeta,tmp
  INTEGER :: t
!--------------------------------------------------------------------

!- Compute offshore wave length, local adimensional wave length, and
!- orbital velocity according to linear wave theory
  if(depth<=0) call parallel_abort('SED_TRANS: (1)')
  psi = 4.d0*pi**2.d0*depth/(g*wave_per**2.d0) !>0
  !kh>0
  IF (psi .LE. 1.d0) THEN
      kh = DSQRT(psi)*(1.d0 + 0.2d0*psi)
  ELSE
      kh = psi*(1.d0 + 0.2d0*DEXP(2.d0-2.d0*psi))
  ENDIF
  Ladim = 2.d0*pi/kh !>0
  Hadim = H/depth
  uorb = pi*Hadim*depth/wave_per/DSINH(kh)

!- Compute bed slope in the direction of wave propagation, irribaren 
!- number, and Ursell number
  Slope = -DSQRT(dhxi**2.d0+dhyi**2.d0)*DCOS(w_dir+DATAN2(dhyi,dhxi))
  if(Hadim<0.or.Ladim<=0) call parallel_abort('SED_TRANS: (2)')
  Zeta = DTAN(Slope)/DSQRT(Hadim/Ladim)
  Ur = Hadim*Ladim**2.d0

!- Compute the normalized maximal orbital velocity
  C1 = Ladim-10.d0
  C2 = DABS(C1-(Hadim-DABS(Zeta)))
  C3 = Zeta*(1.d0-C1)
  C4 = DTANH(DABS(C3-C2)/Ur)
  tmp=DABS(Zeta)+DTANH(C4)
  if(Hadim<=0.or.tmp<0) call parallel_abort('SED_TRANS: (3)')
  !Ur>0
  C5 = DSQRT(tmp) !DABS(Zeta)+DTANH(C4))
  P1 = DSQRT(Hadim)-C5*Hadim
  U1 = b1*P1+a1

!- Compute the velocity asymmetry parameter
  D1 = 3.d0*Zeta+2.d0*Ladim/Ur
  D2 = DSQRT(Ladim)-DTANH(DABS(D1))
  D3 = (2.d0*Zeta+DSQRT(Ladim/Ur))**2.d0
  D4 = Ur+Ladim/D3/Ur
!  if(D2/D4<0) call parallel_abort('SED_TRANS: (4)')
  IF (D2/D4<=0) THEN
    D5 = 0.d0
  ELSE
    D5 = DSQRT(D2/D4)
  ENDIF
  P2 = 1.2001d0*D5+0.4758d0
  U2 = b2*P2+a2

!- Compute the normalized phase of wave crest
  E1 = Hadim*Ladim*Zeta
  E2 = E1*(-9.8496d0*Zeta*Hadim)**2.d0
  E3 = DTANH(E2)+DTANH(E1)+Ladim-1.d0
  P3 = DTANH(-9.3852d0/E3)
  T1 = b3*P3+a3

!- Compute the normalized phase of zero down crossing
  F1 = 0.0113d0*Zeta*Ladim**2.d0
  F2 = 0.00035667d0*Zeta*Ladim**4.d0
  F3 = 0.1206d0*Ladim*DTANH(DTANH(Zeta))
  F4 = Hadim*DTANH(F2)/DTANH(F3)
  P4 = Hadim*DTANH(0.02899d0*Ladim*F1)-DTANH(F4)
  T0 = b4*P4+a4

!- Compute the normalized phase of wave trough
  G1 = Zeta+0.9206d0
  G2 = Ladim-DSQRT(Ur)+DSQRT(2.5185d0/Ladim)-4.6505d0
  G3 = DSQRT(DABS(G2/Hadim))
  G4 = DABS(Zeta+Ladim)-4.4995d0+Zeta
  G5 = DABS(G4+DABS(Zeta)-5.3981d0)
  G6 = DABS(Ladim+DSQRT(3.0176d0/Hadim)-5.2868d0+Hadim)
  G7 = DABS(Zeta+0.1950d0*(G6+Zeta))
  G8 = DABS(Zeta)+Ladim
  P5 = 4.1958d0/(G1+G3+G5+G7+G8)
  T2 = b5*P5+a5

!- Compute orbital velocity parameters and correct U0
  Uairy = uorb
  Ustar = 2.d0*U2*Uairy
  Ucrest = U1*Ustar
  Utrough = Ustar-Ucrest
  U0 = (Ucrest*T0-Utrough*(1.d0-T0))/(T0-T1)
  IF (U0 .GT. (0.25d0*Ucrest)) THEN
      U0 = 0.25d0*Ucrest
  ENDIF

  tmp=T0*(U0-Ucrest-Utrough)
  if(tmp==0) call parallel_abort('SED_TRANS: (5)')
  a1bis = (-Utrough+T1*U0)/tmp !(T0*(U0-Ucrest-Utrough))
  IF (a1bis .LT. 0.99d0) THEN
    T0 = 0.99d0*T0
    if(T0==1) call parallel_abort('SED_TRANS: (6)')
    Utrough = (-(T0-T1)*U0+T0*Ucrest)/(1.d0-T0)
  ELSE
    T0 = a1bis*T0
  ENDIF

!- Compute wave half-periods
  T_crest = T0*wave_per
  T_trough = wave_per-T_crest

!- Rebuilt of orbital velocity over one wave period
  if(ech==1) call parallel_abort('SED_TRANS: (7)')
  DO t=1,ech-1
    tv = (t-1.d0)/(ech-1.d0)
    IF ((tv .GE. 0.d0) .AND. (tv .LE. T1)) THEN
      if(T1==0) call parallel_abort('SED_TRANS: (8)')
      Uorbi(t) = Ucrest*DSIN(pi/2.d0*tv/T1)
    ELSEIF ((tv .GT. T1) .AND. (tv .LE. T0)) THEN
      if(T1==T0) call parallel_abort('SED_TRANS: (9)')
      Uorbi(t) = Ucrest*DCOS(pi/2.d0*(tv-T1)/(T0-T1))                &
                 - U0*DSIN(pi*(tv-T1)/(T0-T1))
    ELSEIF ((tv .GT. T0) .AND. (tv .LE. T2)) THEN
      if(T2==T0) call parallel_abort('SED_TRANS: (10)')
      Uorbi(t) = -Utrough*DSIN(pi/2.d0*(tv-T0)/(T2-T0))
    ELSEIF ((tv .GT. T2) .AND. (tv .LE. 1.d0)) THEN
      if(T2==1) call parallel_abort('SED_TRANS: (11)')
      Uorbi(t) = -Utrough*DCOS(pi/2.d0*(tv-T2)/(1.d0-T2))
    ENDIF
  ENDDO

  END SUBROUTINE wave_asymmetry_Elfrink
!--------------------------------------------------------------------
END MODULE sed2d_transport
