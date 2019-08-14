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

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!    Time loop part of SCHISM
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      subroutine schism_step(it)

      use schism_glbl
      use schism_msgp
      use misc_modules

#ifdef USE_GOTM
      use turbulence, only: do_turbulence, cde, tke1d => tke, eps1d => eps, L1d => L, num1d => num, nuh1d => nuh
!      use mtridiagonal, only: init_tridiagonal
#endif

#ifdef USE_ECO
      USE bio_param
      USE biology
      USE eclight
#endif

#ifdef USE_ICM
      USE icm_mod, only: iSun,iWQPS,nps,DTD,WWPRPOC,WWPLPOC, &
                          &WWPDOCA,WWPRPON,WWPLPON,WWPDON,WWPNH4,WWPNO3, &
                          &WWPRPOP,WWPLPOP,WWPDOP,WWPPO4t,WWPSU,WWPSAt, &
                          &WWPCOD,WWPDO,xPSQ,xPSK,PRPOC,PLPOC,PDOCA,PRPON, &
                          &PLPON,PDON,PNH4,PNO3,PRPOP,PLPOP,PDOP,PPO4t,PSU, &
                          &PSAt,PCOD,PDO,WMS,irea    !added by YC
      USE icm_sed_mod, only: sed_BENDO,CTEMP,BBM,CPOS,PO4T2TM1S,NH4T2TM1S,NO3T2TM1S, &
                            &HST2TM1S,CH4T2TM1S,CH41TM1S,SO4T2TM1S,SIT2TM1S,BENSTR1S,CPOP,CPON,CPOC,&
                            &NH41TM1S,NO31TM1S,HS1TM1S,SI1TM1S,PO41TM1S,PON1TM1S,PON2TM1S,PON3TM1S,POC1TM1S,POC2TM1S,&
                            &POC3TM1S,POP1TM1S,POP2TM1S,POP3TM1S,PSITM1S,BFORMAXS,ISWBENS,DFEEDM1S, &  !added by wangzg
                            &BENDOC,SED_BENNH4,SED_BENNO3,BENPO4,SED_BENCOD,BENSA
#endif

#ifdef USE_NAPZD
      USE biology_napzd
#endif

#ifdef USE_SED
       USE sed_mod, only : Wsed,Srho,Nbed,MBEDP,bedldu,bedldv,bed,bottom,    &
                          &bed_frac,mcoefd,bed_fracn,bed_d50n,bed_taun,&
                          &bedforms_rough,bed_rough,izcr,izsw,izwr,izbld, &
                          &bed,ithck,iaged
#endif

#ifdef USE_SED2D
      use sed2d_mod, only : Cdsed,cflsed,d50moy,dpdxy,nb_class,qav,  &
                           &qb,qs,qtot,z0cr_e,z0_e,z0sw_e,z0wr_e, &
                           &idrag_sed2d=>idrag
#endif

#ifdef USE_OIL
#endif

#ifdef USE_HA
      USE harm
#endif

      USE hydraulic_structures

      implicit none
!      implicit real(rkind)(a-h,o-z),integer(i-n)
!#ifndef USE_MPIMODULE
      include 'mpif.h'
!#endif

      integer, intent(in) :: it

!     External functions
      integer :: kronecker
      real(rkind) :: eqstate

!     Local variables
      integer :: istat,i,j,k,l,m,kk,mm,jj,ll,lll,nd,nd0,ie,ie0,iegb,icount, &
                 &icount1,icount2,icount3,jsj,k0,k1,k2,ipgb,ndgb1,ndgb2, &
                 &irank,jblock,jface,n1,n2,n3,ifl,isd,isd0,isd1,isd2,isd3, &
                 &ibnd,jfr,ncyc,iter,nlev,klev,kin,nqdim,limit,jmin, &
                 &ipsgb,iadvf,ifl_bnd,nnel,jlev,ndelt_min,ndelt_max,ii,id, &
                 &id2,id3,ip,ndelt,ibot_fl,ibelow,indx,nj,ind,ind2,lim,in1, &
                 &in2,irank_s,itmp,itmp1,itmp2,node1,node2,ndim,mk,nd_lam, &
                 &iee,idel,irow,icol,ieq,ij,kbb,lwrite,lit,ihot_len,IHOTSTP, &
                 &itmpf,ibt,mmk,ndo,n,ibtm
      real(rkind) :: cwtmp,cwtmp2,wtmp1,wtmp2,time,ramp,rampbc,rampwind,dzdx,dzdy, &
                     &dudz,dvdz,dudx,dudx2,dvdx,dvdx2,dudy,dudy2,dvdy,dvdy2, &
                     &dzz1,ta,wx2,wy2,wtratio,sum1,sum2,sum3,sum4,dragcmin, &
                     &dragcmax,wmag,vmag,vmag1,vmag2,dragcoef,tmp,tmp0,tmp1, &
                     &tmp2,theta,x1,stratio,rat,htot,ar,vnth0,arg,bthick, &
                     &taubx,tauby,ybm,wfr,wdir,z0b,fw,delta_wc,vmax,vmin, &
                     &drhodz,bvf,shear2,rich,u_taus,u_taub,ztmp,toth,z0s, &
                     &vts0,xctr2,yctr2,zctr2,dists,distb,fwall,q2fs,q2bot, &
                     &xlfs,xlbot,prod,buoy,diss,psi_n,psi_n1,q2l,upper, &
                     &xl_max,vd,td,qd1,qd2,t0,s0,rot_per,rot_f,xt,yt,zt, &
                     &xt4,yt4,zt4,uuint,vvint,wwint,vis_coe,suma,dtbk,eps, &
                     &time_rm,time_rm2,u1,u2,v1,v2,eic,eta_min,zmax,xn1,yn1, &
                     &xn2,yn2,x10,x20,y10,y20,bb1,bb2,rl10,rl20,delta, &
                     &sintheta,tau_x,tau_x2,tau_y,tau_y2,detadx,detady,dprdx, &
                     &dprdy,detpdx,detpdy,chigamma,ubstar,vbstar,hhat_bar, &
                     &h_bar,bigf1,bigf2,botf1,botf2,ub2,vb2,bigu1,bigu2,bigu3, &
                     &bigv1,bigv2,bigv3,av_elem_x,av_elem_y,sdbtu,sdbtv, &
                     &hat_gam_x,hat_gam_y,del,gam_x,gam_y,horx,hory,rs1,rs2, &
                     &bigfc1,bigfc2,dot1,dot2,dot3,hhatb,avg2,etam,tmpj,tmpj1, &
                     &fac,dep,ubed,vbed,wbed,dpdx,dpdy,vnorm,bigvn,vn1,vn2, &
                     &utmp,vtmp,ri3,con0,Unbar,ss,etatot,etatotl,tmpx,tmpy, &
                     &tmpxs,tmpys,tmpx1,tmpy1,tmpx2,tmpy2,tmpx3,tmpy3, &
                     &tmpx1s,tmpy1s,tmpx2s,tmpy2s,tmpx3s,tmpy3s,taux2,tauy2, &
                     &taux2s,tauy2s,uths,vths,vtan,suru,surv,dhdx,dhdy,ubar1, &
                     &ubar2,vbar1,vbar2,ubar,vbar,eta1_bar,eta2_bar,qnon_e1, &
                     &qnon_e2,xcon,ycon,zcon,vnor1,vnor2,bflux,bflux0,top, &
                     &deta_dx,deta_dy,hmin,dzds_av,css,dsigma,dgam0,dgam1, &
                     &hat_i0,dzds,dsdx,dsdy,dsig2,hat_ir,vol,dz,tmp_max, &
                     &tmp_max_gb,dia_min,dia_min_gb,df_max,qhat_e1,qhat_e2,dqdz,uvnu, &
                     &av_bdef1,av_bdef2,depth,zz1,rr,d_1,d_2,smin,smax,tmin, &
                     &tmax,vnn,vnf,snu,tnu,evap,precip,sflux_e,dp1,dp2,srad1, &
                     &srad2,bigv,tsel01,tsel02,snd_nu_e,tnd_nu_e,zrat,tt1,ss1, &
                     &cff1,cff2,difnum_max_l,total_loading,trnu, &
                     &av_df,vol1,tot_heat,tot_salt,tot_heat_gb, &
                     &tot_salt_gb,dav_mag,tvol,tmass,tpe,tkne,enerf,ener_ob, &
                     &av_dep,vel_m1,vel_m2,xtmp,ytmp,ftmp,tvol12,fluxbnd, &
                     &fluxchan,fluxchan1,fluxchan2,tot_s,flux_s,ah,ubm,ramp_ss,Cdmax, &
                     &wmag_e,wmag_factor, & !wmag_e and wmag_facotr addedy by wangzg
                     &bthick_ori,zsurf

!     Output handles
      character(len=72) :: it_char
      character(len=72) :: fgb  ! Processor specific global output file name
      integer :: lfgb       ! Length of processor specific global output file name
      real(4) :: floatout,floatout2

!     Inter-subdomain backtracking
      logical :: lbt(1),lbtgb(1)
!      logical :: lbt_l(1), lbtgb_l(1)
      integer :: nbtrk
      type(bt_type) :: btlist(mxnbt) !to avoid conflict with inter_btrack()

!     Solver arrays for TRIDAG
      real(rkind) :: alow(max(3,nvrt)),bdia(max(3,nvrt)),cupp(max(3,nvrt)),rrhs(100,nvrt), &
                    &soln(100,nvrt),gam(nvrt)

!     Non-hydrostatic arrays
      real(rkind),allocatable :: qhat(:,:),dqnon_dxy(:,:,:),qmatr(:,:,:,:),qir(:,:)

!     Misc 
      integer :: kbs_e(nsa),nwild(nea+12),nwild2(ne_global),nsubd(npa),icoef(npa+1), &
                 &jcoef(npa*(mnei+1)),ibt_p(npa),ibt_s(nsa)
      real(rkind) :: dfz(2:nvrt),dzz(2:nvrt),deta1_dx(nsa),deta1_dy(nsa),deta2_dx(nsa), &
                     &deta2_dy(nsa),dpr_dx(nsa),dpr_dy(nsa),detp_dx(nsa),detp_dy(nsa), &
                     &sne(3,nvrt),area_e(nvrt),srad_e(nea),qel(np),elbc(npa),hhat(nsa), &
                     &bubt(2,nea),bigu(2,nsa),ghat1(2,nea),etp(npa),h1d(0:nvrt),SS1d(0:nvrt), &
                     &NN1d(0:nvrt),q2tmp(nvrt),xltmp(nvrt),rzbt(nvrt),shearbt(2:nvrt), &
                     &xlmax(nvrt),cpsi3(2:nvrt),cpsi2p(2:nvrt),q2ha(2:nvrt),xlha(2:nvrt), &
                     &vsource(nea)
      real(rkind) :: swild(nsa+nvrt+12+ntracers),swild2(nvrt,12),swild10(max(3,nvrt),12), &
     &swild3(20+ntracers),swild4(2,3+2*ntracers)
      real(4) :: swild8(nvrt,2) !used in ST nudging
      logical :: lelbc(npa)
      
      real(4),allocatable :: swild9(:,:) !used in tracer nudging
      real(rkind),allocatable :: rwild(:,:) 
      real(rkind),allocatable :: swild99(:,:),swild98(:,:,:) !used for exchange (deallocate immediately afterwards)
      real(rkind),allocatable :: hp_int(:,:,:),buf1(:,:),buf2(:,:),buf3(:),msource(:,:)
      real(rkind),allocatable :: fluxes_vol(:),fluxes_vol_gb(:) !volume fluxes output between regions
      logical :: ltmp,ltmp1(1),ltmp2(1)

      logical,save :: first_call=.true.
      logical :: up_tvd
#ifdef DEBUG
      real(rkind),allocatable :: bpgr(:,:),wafo(:,:,:)
#endif
!     Tracers
!      integer :: flag_model,flag_ic
      real(rkind),allocatable :: Bio_bdefp(:,:),tr_tc(:,:),tr_tl(:,:)

#ifdef USE_WWM
      CHARACTER(LEN=3) :: RADFLAG
      REAL(rkind) ::  dJ_dx(nsa) ! modif AD
      real(rkind),allocatable :: stokes_w(:,:),stokes_w_nd(:,:), &
     &stokes_vel_sd(:,:,:) !,jpress(:),sbr(:,:),sbf(:,:),stokes_vel(:,:,:)
#endif /*USE_WWM*/

!     End of declarations

      if(nonhydro==1) then
        allocate(qhat(nvrt,npa),dqnon_dxy(2,nvrt,nsa),qmatr(nvrt,-1:1,0:(mnei+1),np), &
     &qir(nvrt,np),stat=istat)
        if(istat/=0) call parallel_abort('STEP: Nonhydro allocation failure')
!'
      endif

      allocate(hp_int(nvrt,nea,2),stat=istat)
      if(istat/=0) call parallel_abort('STEP: other allocation failure')

      if(ntracers>0.and.inu_tr==2) then
        allocate(swild9(ntracers,nvrt),stat=istat)
        if(istat/=0) call parallel_abort('STEP: alloc failure (3)')
      endif

#ifdef DEBUG
      allocate(bpgr(nsa,2),wafo(nvrt,nsa,2))
#endif
!     Source
      if(if_source==1) then
        allocate(msource(2+ntracers,nea),stat=istat)
        if(istat/=0) call parallel_abort('STEP: allocation failure (2)')
      endif !if_source

#ifdef USE_NAPZD
      allocate(Bio_bdefp(nvrt,np), stat=istat)
      if(istat/=0) call parallel_abort('STEP: NAPZD allocation failure')
#endif

#ifdef USE_SED
       allocate(tr_tc(ntracers,nea),tr_tl(ntracers,nea),stat=istat)
       if(istat/=0) call parallel_abort('STEP: sed. allocation failure')
#endif

#ifdef USE_WWM
       allocate(stokes_w(nvrt,nea),stokes_w_nd(nvrt,npa), &
     &stokes_vel_sd(2,nvrt,nsa),stat=istat)
       if(istat/=0) call parallel_abort('STEP: WWM allocation failure')
#endif

!      End alloc.

!      do it=iths+1,ntime

#ifdef INCLUDE_TIMING
      wtmp1=mpi_wtime() !Forcing preparation section
#endif

      time=it*dt 
     
!     Broadcast to global module
      time_stamp=time; it_main=it

!...  define ramp function for boundary elevation forcing, wind and pressure
!...  forcing and tidal potential forcing
!...
      if(ibc==0) then
        if(nrampbc/=0) then
          rampbc=tanh(2*time/86400/drampbc)
        else
          rampbc=1
        endif
      endif

      if(nws>0.and.nrampwind/=0) then
        rampwind=tanh(2*time/86400/drampwind)
      else
        rampwind=1
      endif

      if(nramp==1) then
        ramp=tanh(2*time/86400/dramp)
        !For sink/source
        ramp_ss=tanh(2*time/86400/dramp/2)
      else
        ramp=1
        ramp_ss=1
      endif

!...  Compute new bed deformation
      do i=1,npa
        bdef2(i)=bdef(i)/ibdef*min0(it,ibdef)
      enddo !i

!...  Horizontal viscosity; compute d2uv (incorporated hvis inside)
!     Pre-compute some derivatives in each element to speed up
!     Use sdbt(1:4,nvrt,nsa) (nsa>=nea checked) as temporary array.
      d2uv=0

      if(ihorcon/=0) then
        sdbt=0 !for dry and below bottom spots
        do i=1,nea
          if(idry_e(i)==1) cycle

!         Wet
          do k=kbe(i),nvrt
            k2=min(k+1,nvrt)
            k1=max(k-1,kbe(i))
            dzdx=0; dzdy=0; dudz=0; dvdz=0; dudx2=0; dudy2=0; dvdx2=0; dvdy2=0
            do j=1,3 !nodes
              nd=elnode(j,i)
              dzdx=dzdx+znl(max(k,kbp(nd)),nd)*dldxy(j,1,i)
              dzdy=dzdy+znl(max(k,kbp(nd)),nd)*dldxy(j,2,i)
              dzz1=znl(max(k2,kbp(nd)),nd)-znl(max(k1,kbp(nd)),nd)
!              if(k1==k2.or.dzz1<=0) then
              if(dzz1<0) then
                write(errmsg,*)'Impossible 91:',k1,k2,dzz1
                call parallel_abort(errmsg)
              endif
              if(dzz1/=0) then
                dudz=dudz+(uu2(k2,nd)-uu2(k1,nd))/dzz1/3
                dvdz=dvdz+(vv2(k2,nd)-vv2(k1,nd))/dzz1/3
              endif
              dudx2=dudx2+uu2(k,nd)*dldxy(j,1,i) !on S-plane for k>=kz
              dudy2=dudy2+uu2(k,nd)*dldxy(j,2,i)
              dvdx2=dvdx2+vv2(k,nd)*dldxy(j,1,i)
              dvdy2=dvdy2+vv2(k,nd)*dldxy(j,2,i)
            enddo !j
            if(k<kz) then !z-levels
              dzdx=0; dzdy=0
            endif
            sdbt(1,k,i)=dudx2-dudz*dzdx !dudx; whole level
            sdbt(2,k,i)=dudy2-dudz*dzdy !dudy
            sdbt(3,k,i)=dvdx2-dvdz*dzdx !dvdx
            sdbt(4,k,i)=dvdy2-dvdz*dzdy !dvdy
          enddo !k=kbe(i),nvrt
        enddo !i=1,nea

        allocate(swild98(2,nvrt,nsa),stat=istat)
        if(istat/=0) call parallel_abort('MAIN: fail to allocate swild98 (3)')
!'
        swild98=0
        do i=1,ns !residents
          if(idry_s(i)==1) cycle

          do k=1,nvrt
            dudx=0; dudy=0; dvdx=0; dvdy=0
            icount=0
            do j=1,2
              ie=isdel(j,i)
              if(ie/=0) then
                icount=icount+1
                dudx=dudx+sdbt(1,k,ie)
                dudy=dudy+sdbt(2,k,ie)
                dvdx=dvdx+sdbt(3,k,ie)
                dvdy=dvdy+sdbt(4,k,ie)
              endif !ie/=0
            enddo !j=1,2
            if(icount==0) then
              write(errmsg,*)'MAIN: Impossible 78'
              call parallel_abort(errmsg)
            endif
            dudx=dudx/icount
            dudy=dudy/icount
            dvdx=dvdx/icount
            dvdy=dvdy/icount
            swild98(1,k,i)=dudx*sframe(1,1,i)+dudy*sframe(2,1,i) !dudn; whole level
            swild98(2,k,i)=dvdx*sframe(1,1,i)+dvdy*sframe(2,1,i) !dvdn
            !Deal with free slip land bnd; this is just a back-up case as land bnd will
            !be dealt with separately below
            !if(isbs(i)==-1) swild98(:,k,i)=0 
          enddo !k=1,nvrt
        enddo !i=1,ns

!       Update ghost 
#ifdef INCLUDE_TIMING
        cwtmp=mpi_wtime()
#endif
        call exchange_s3d_2(swild98)
#ifdef INCLUDE_TIMING
        wtimer(3,2)=wtimer(3,2)+mpi_wtime()-cwtmp
#endif

!       Compute horizontal viscosity term: \nabla\cdot(\miu\nabla \vector{u})
        do j=1,ns !residents only
          do k=kbs(j)+1,nvrt !viscosity = 0 at bottom
            ta=0
            do l=1,2 !element
              ie=isdel(l,j)
              if(ie/=0) then
                ta=ta+area(ie)
                do i=1,3 !sides
                  jsj=elside(i,ie)
                  if(isbs(jsj)==-1) then !deal with land bnd
                    tmp=sqrt(su2(k,jsj)**2+sv2(k,jsj)**2)
                    d2uv(1,k,j)=d2uv(1,k,j)-distj(jsj)*cdh*tmp*su2(k,jsj)
                    d2uv(2,k,j)=d2uv(2,k,j)-distj(jsj)*cdh*tmp*sv2(k,jsj)
                  else if(isdel(2,j)==0.or.jsj/=j) then
                    d2uv(1:2,k,j)=d2uv(1:2,k,j)+ssign(i,ie)*distj(jsj)*hvis(k,ie)*swild98(1:2,k,jsj)
                  endif
                enddo !i
              endif !ie/=0
            enddo !l
            if(ta==0) then
              write(errmsg,*)'MAIN: Impossible 77'
              call parallel_abort(errmsg)
            endif
            d2uv(1:2,k,j)=d2uv(1:2,k,j)/ta
          enddo !k
        enddo !j=1,ns

!       Update ghost 
#ifdef INCLUDE_TIMING
        cwtmp=mpi_wtime()
#endif
        call exchange_s3d_2(d2uv)
#ifdef INCLUDE_TIMING
        wtimer(3,2)=wtimer(3,2)+mpi_wtime()-cwtmp
#endif
        deallocate(swild98)
      endif !ihorcon/=0

      if(myrank==0) write(16,*)'done hvis and bottom fric: '

!...  Earth tidal potential at nodes: pre-compute to save time
!...
      do i=1,npa
        etp(i)=0
        do j=1,ntip
          ncyc=int(tfreq(j)*time/2/pi)
          arg=tfreq(j)*time-ncyc*2*pi+jspc(j)*xlon(i)+tear(j)
          etp(i)=etp(i)+ramp*tamp(j)*tnf(j)*fun_lat(jspc(j),i)*cos(arg)
        enddo !j
      enddo !i

!...  process new wind info 
!...  Wind vectors always in lat/lon frame 
      if(nws==1) then
        if(time>=wtime2) then
          wtime1=wtime2
          wtime2=wtime2+wtiminc
          read(22,*)tmp,wx2,wy2
          windx1=windx2
          windy1=windy2
          windx2=wx2
          windy2=wy2
        endif

        wtratio=(time-wtime1)/wtiminc
        do i=1,npa
          windx(i)=windx1(i)+wtratio*(windx2(i)-windx1(i))
          windy(i)=windy1(i)+wtratio*(windy2(i)-windy1(i))
        enddo !i
      endif !nws=1

      if(nws==4) then
        if(time>=wtime2) then
          wtime1=wtime2
          wtime2=wtime2+wtiminc
          windx1=windx2
          windy1=windy2
          pr1=pr2
!          The large array for nws=4 option (may consider changing to
!          unformatted binary read)
          allocate(rwild(np_global,3),stat=istat)
          if(istat/=0) call parallel_abort('STEP: failed to alloc. (71)')
          read(22,*)tmp,rwild(:,:) 
          do i=1,np_global
            if(ipgl(i)%rank==myrank) then
              nd=ipgl(i)%id
              windx2(nd)=rwild(i,1)
              windy2(nd)=rwild(i,2)
              pr2(nd)=rwild(i,3)
            endif
          enddo !i
          deallocate(rwild)
        endif

        wtratio=(time-wtime1)/wtiminc
        do i=1,npa
          windx(i)=windx1(i)+wtratio*(windx2(i)-windx1(i))
          windy(i)=windy1(i)+wtratio*(windy2(i)-windy1(i))
          pr(i)=pr1(i)+wtratio*(pr2(i)-pr1(i))
        enddo !i
      endif !nws=4

#ifdef USE_ICM
! calculating WMS used for reareation,added by wangzg
      if(irea==1) then
        do i=1,nea
           if(idry_e(i)==1) cycle
            n1=elnode(1,i)
            n2=elnode(2,i)
            n3=elnode(3,i)
            wmag_e=sqrt(windx(n1)**2+windy(n1)**2)+sqrt(windx(n2)**2+windy(n2)**2)+sqrt(windx(n3)**2+windy(n3)**2)
            wmag_factor=(windfactor(n1)+windfactor(n2)+windfactor(n3))/3.
            wmag_e=wmag_e/3.
            WMS(i)=wmag_e !*wmag_factor !no windfactor for DO reareation
        enddo !i
      endif !irea=1
#endif

!     CORIE mode
      if(nws>=2.and.nws<=3) then
        if(time>=wtime2) then
!...      Heat budget & wind stresses
          if(ihconsv/=0) then
!#ifdef USE_SFLUX
            call surf_fluxes(wtime2,windx2,windy2,pr2,airt2,shum2,srad,fluxsu,fluxlu,hradu,hradd,tauxz,tauyz, &
#ifdef PREC_EVAP
     &                       fluxprc,fluxevp, &
#endif
     &                       nws,fluxsu00,srad00)
!#endif
            do i=1,npa
              sflux(i)=-fluxsu(i)-fluxlu(i)-(hradu(i)-hradd(i))
            enddo
            if(myrank==0) write(16,*)'heat budge model completes...'
          endif !ihconsv.ne.0

          wtime1=wtime2
          wtime2=wtime2+wtiminc
          do i=1,npa
            windx1(i)=windx2(i)
            windy1(i)=windy2(i)
            pr1(i)=pr2(i)
            airt1(i)=airt2(i)
            shum1(i)=shum2(i)
          enddo
          call get_wind(wtime2,windx2,windy2,pr2,airt2,shum2)
        endif !time>=wtime2

        wtratio=(time-wtime1)/wtiminc
        do i=1,npa
          windx(i)=windx1(i)+wtratio*(windx2(i)-windx1(i))
          windy(i)=windy1(i)+wtratio*(windy2(i)-windy1(i))
          pr(i)=pr1(i)+wtratio*(pr2(i)-pr1(i))
        enddo !i

!       Overwrite wind with wind.th
!        read(22,*)tmp,wx2,wy2
!        windx1=wx2; windy1=wy2
!        windx2=wx2; windy2=wy2
!        windx=wx2; windy=wy2
!       End

!       Read in new flux values for next step
        if(nws==3) read(23,*) tmp,fluxsu00,srad00
      endif !nws>=2

!...  Re-scale wind
      if(nws>0) then
        do i=1,npa
          windx(i)=windx(i)*windfactor(i)
          windy(i)=windy(i)*windfactor(i)
        enddo !i
      endif

!-------------------------------------------------------------------------------
!   Wind wave model (WWM)
!-------------------------------------------------------------------------------
#ifdef USE_WWM
      if(mod(it,nstep_wwm)==0) then
        wtmp1=mpi_wtime()
        if(myrank==0) write(16,*)'starting WWM'
        call WWM_II(it,icou_elfe_wwm,dt,nstep_wwm,RADFLAG)

!       Outputs (via datapool):
!       sbr(2,npa): momentum flux vector due to wave breaking (nearshore depth-induced breaking; see Bennis 2011)
!       sbf(2,npa): momentum lost by waves due to the bottom friction (not used for the moment)
!       stokes_vel(2,nvrt,npa): Stokes velocity
!       jpress(npa): waved-induced pressure
!       wwave_force(2,nvrt,nsa): =0 if icou_elfe_wwm=0. In [e,p]frame (not sframe!).
!       wwave_force(1:2,:,1:nsa) = Rsx, Rsy in my notes (the terms in momen. eq.)
!       and has a dimension of m/s/s. This is overwritten under Vortex
!       formulation later.
!       out_wwm_windpar(npa,10): 
!         1) = WINDXY(IP,1) ! wind vector u10,x
!         2) = WINDXY(IP,2) ! wind vector u10,y
!         3) = SQRT(WINDXY(IP,1)**2.+WINDXY(IP,2)**2.) ! wind magnitutde u10
!         4) = TAUW(IP)     ! wave stress from the discrete part of the spectra
!         5) = TAUHF(IP)    ! high freq. part of the waves.
!         6) = TAUTOT(IP)   ! total stress of the wave
!         7) = Z0(IP)       ! apparent rougnes lengths (m)
!         8) = UFRIC(IP)    ! ustar - frictional vel. (m/s)
!         9) = ALPHA_CH(IP) ! Charnock Parameter gz0/ustar**2
!        10) = CD(IP)       ! Drag Coefficient

!       out_wwm(npa,30): output variables from WWM (all 2D); see names in NVARS() in the routine
!                      BASIC_PARAMETER() in wwm_initio.F90 for details; below is a snapshot from there:
         !OUTPAR(1)   = HS       ! Significant wave height
         !OUTPAR(2)   = TM01     ! Mean average period
         !OUTPAR(3)   = TM02     ! Zero down crossing period for comparison with buoy.
         !OUTPAR(4)   = TM10     ! Average period of wave runup/overtopping ...
         !OUTPAR(5)   = KLM      ! Mean wave number
         !OUTPAR(6)   = WLM      ! Mean wave length
         !OUTPAR(7)   = ETOTS    ! Etot energy in y-direction
         !OUTPAR(8)   = ETOTC    ! Etot energy in x-direction
         !OUTPAR(9)   = DM       ! Mean average energy transport direction
         !OUTPAR(10)  = DSPR     ! Mean directional spreading
         !OUTPAR(11)  = FPP      ! Discrete peak period (sec)
         !OUTPAR(12)  = TPP      ! Continuous peak period based on higher order moments (sec) 
         !OUTPAR(13)  = CPP      ! Peak phase vel. (m/s)
         !OUTPAR(14)  = WNPP     ! Peak n-factor
         !OUTPAR(15)  = CGPP     ! Peak group vel.
         !OUTPAR(16)  = KPP      ! Peak wave number
         !OUTPAR(17)  = LPP      ! Peak wave length.
         !OUTPAR(18)  = PEAKD    ! Peak (dominant) direction (degr)
         !OUTPAR(19)  = PEAKDSPR ! Peak directional spreading
         !OUTPAR(20)  = DPEAK    ! Discrete peak direction
         !OUTPAR(21)  = UBOT     ! Orbital vel. (m/s)
         !OUTPAR(22)  = ORBITAL  ! RMS Orbital vel. (m/s)
         !OUTPAR(23)  = BOTEXPER ! Bottom excursion period.
         !OUTPAR(24)  = TMBOT    ! Bottom wave period (sec)
         !OUTPAR(25)  = URSELL   ! Uresell number based on peak period ...
         !OUTPAR(26)  = UFRIC(IP)    ! Friction velocity 
         !OUTPAR(27)  = Z0(IP)       ! Rougness length
         !OUTPAR(28)  = ALPHA_CH(IP) ! Charnock coefficient
         !OUTPAR(29)  = CD(IP)       ! Drag coefficient
         !OUTPAR(30)  = WINDXY(IP,1) ! windx
         !OUTPAR(31)  = WINDXY(IP,2) ! windy

        if(myrank==0) write(16,*)'WWM-RS part took (sec) ',mpi_wtime()-wtmp1

        !Check outputs from WWM
        sum1=sum(out_wwm_windpar(1:npa,1:10))
        sum2=sum(wwave_force)
        sum3=sum(out_wwm(:,1:35))
        if(sum1/=sum1.or.sum2/=sum2.or.sum3/=sum3) then
          if(sum1/=sum1) then
            do i=1,9
              write(12,*)'sum1:',i,sum(out_wwm_windpar(:,i))
            end do
          endif
          if(sum3/=sum3) then
            do i=1,31 
              sum4=sum(out_wwm(:,i))
              write(12,*)'sum4:',i,sum4
              if(sum4/=sum4) then
                do j=1,npa
                  write(12,*)i,j,out_wwm(j,i)
                enddo !j
              endif
            enddo !i
          endif !sum3
          write(errmsg,*)'NaN from WWM:',sum1,sum2,sum3
          call parallel_abort(errmsg)
        endif !sum
      endif !mod()

!     Caculate vortex force of Bennis (2011)
!     Overwrites wwave_force(2,nvrt,nsa) (in eframe if ics=2)
      if(RADFLAG.eq.'VOR') then
        wwave_force=0

        !Stokes vel. at side (in pframe if ics=2)
        stokes_vel_sd=0
        do i=1,nsa
          if(idry_s(i)==0) then
            n1=isidenode(1,i); n2=isidenode(2,i)
            stokes_vel_sd(:,:,i)=(stokes_vel(:,:,n1)+stokes_vel(:,:,n2))/2
          endif
        enddo !i

        !Vortex terms
        call hgrad_nodes(2,0,nvrt,npa,nsa,uu2,dr_dxy)
        do i=1,ns
          if(idry_s(i)==0) then
            n1=isidenode(1,i); n2=isidenode(2,i)
            !f*v_s-du/dy*v_s
            wwave_force(1,:,i)=wwave_force(1,:,i)+(cori(i)-dr_dxy(2,:,i))*stokes_vel_sd(2,:,i)
            !-f*u_s+du/dy*u_s
            wwave_force(2,:,i)=wwave_force(2,:,i)+(-cori(i)+dr_dxy(2,:,i))*stokes_vel_sd(1,:,i)
          endif
        enddo !i

        call hgrad_nodes(2,0,nvrt,npa,nsa,vv2,dr_dxy)
        do i=1,ns
          if(idry_s(i)==0) then
            n1=isidenode(1,i); n2=isidenode(2,i)
            !dv/dx*v_s
            wwave_force(1,:,i)=wwave_force(1,:,i)+dr_dxy(1,:,i)*stokes_vel_sd(2,:,i)
            !-dv/dx*u_s
            wwave_force(2,:,i)=wwave_force(2,:,i)-dr_dxy(1,:,i)*stokes_vel_sd(1,:,i)
          endif
        enddo !i

        !pressure term 
        do j=1,ns !resident
          if(idry_s(j)==1) cycle

          !Wet side
          icount1=0 
          tmp1=0; tmp2=0
          do l=1,2 !elements
            ie=isdel(l,j)
            if(ie/=0) then; if(idry_e(ie)==0) then
              icount1=icount1+1
              tmp1=tmp1+dot_product(jpress(elnode(:,ie)),dldxy(:,1,ie)) !in eframe
              tmp2=tmp2+dot_product(jpress(elnode(:,ie)),dldxy(:,2,ie))
            endif; endif 
          enddo !l

          if(icount1/=0) then
            tmp1=tmp1/icount1
            tmp2=tmp2/icount1
          endif
          dJ_dx(j)=tmp1
          wwave_force(1,:,j)=wwave_force(1,:,j)-tmp1
          wwave_force(2,:,j)=wwave_force(2,:,j)-tmp2
        enddo !j

!        do j=1,np
!          tmp=sqrt(sbr(1,j)**2+sbr(2,j)**2)
!        enddo
!JZ: what's this? nd is undefined
!        if(myrank==0) write(16,*)'norm of SBR',it,tmp/nd

        do j=1,ns !resident
          n1=isidenode(1,j); n2=isidenode(2,j)
          htot=max(h0,dps(j)+(eta2(n1)+eta2(n2))/2)
         
!          if((idry(n1)==1).or.(idry(n2)==1)) cycle
 
!!         if(idry_s(j)==1.or.isbs(j)>0) cycle          !Wet side

          if(lm2d) then
            if(idry_s(j)==1.or.isbs(j)>0) cycle          !Wet side
            wwave_force(1,:,j)=wwave_force(1,:,j)-sum(sbr(1,isidenode(:,j)))/(2.d0*htot) !/rho0!/grav
            wwave_force(2,:,j)=wwave_force(2,:,j)-sum(sbr(2,isidenode(:,j)))/(2.d0*htot) !/rho0!/grav
          else
             tmp0=sum(out_wwm(isidenode(:,j),1))/2.d0 !Hs
            if(idry_s(j)==1.or.isbs(j)>0.or.tmp0<=0.1) cycle 
!            if(idry_s(j)==1.or.isbs(j)>0) cycle 
!           !Wet side
            swild=0
            do k=kbs(j),nvrt
              swild(k)=cosh(5*sqrt(2.d0)*(zs(k,j)+dps(j)))/tmp0
            enddo !k
    
            sum1=0 !integral
            do k=kbs(j),nvrt-1
              sum1=sum1+(swild(k+1)+swild(k))*(zs(k+1,j)-zs(k,j))/2.d0 !*rho0
            enddo !k
          
            if(sum1==0) call parallel_abort('STEP: integral=0')
             wwave_force(1,:,j)=wwave_force(1,:,j)-swild(1:nvrt)/sum1*sum(sbr(1,isidenode(:,j)))/(2.d0*htot) !/rho0!/grav
             wwave_force(2,:,j)=wwave_force(2,:,j)-swild(1:nvrt)/sum1*sum(sbr(2,isidenode(:,j)))/(2.d0*htot) !/rho0!/grav
          endif !lm2d
        enddo !j

        call exchange_s3d_2(wwave_force)

!       Calculate Stokes w-vel.: stokes_w(nvrt,nea)
        stokes_w=0 !for dry and below bottom levels
        do i=1,nea
          if(idry_e(i)==1) cycle

!	      Wet elements with 3 wet nodes
!         Compute upward normals and areas @ all levels
          n1=elnode(1,i)
          n2=elnode(2,i)
          n3=elnode(3,i)
          if(kbe(i)==0) then
            write(errmsg,*)'Impossible 95'
            call parallel_abort(errmsg)
          endif
          do l=kbe(i),nvrt
            if(ics==1) then
              xcon=(ynd(n2)-ynd(n1))*(znl(max(l,kbp(n3)),n3)-znl(max(l,kbp(n1)),n1))-(ynd(n3)-ynd(n1))* &
     &(znl(max(l,kbp(n2)),n2)-znl(max(l,kbp(n1)),n1))
              ycon=(xnd(n3)-xnd(n1))*(znl(max(l,kbp(n2)),n2)-znl(max(l,kbp(n1)),n1))-(xnd(n2)-xnd(n1))* &
     &(znl(max(l,kbp(n3)),n3)-znl(max(l,kbp(n1)),n1))
              zcon=area(i)*2
            else !lat/lon
              !eframe
              call cross_product(xel(2,i)-xel(1,i),yel(2,i)-yel(1,i),znl(max(l,kbp(n2)),n2)-znl(max(l,kbp(n1)),n1), &
     &                           xel(3,i)-xel(1,i),yel(3,i)-yel(1,i),znl(max(l,kbp(n3)),n3)-znl(max(l,kbp(n1)),n1), &
     &                           xcon,ycon,zcon)
            endif !ics

            area_e(l)=sqrt(xcon**2+ycon**2+zcon**2)/2
            if(area_e(l)==0) then
              write(errmsg,*)'Zero area:',i,l
              call parallel_abort(errmsg)
            endif
            sne(1,l)=xcon/area_e(l)/2 !in eframe
            sne(2,l)=ycon/area_e(l)/2
            sne(3,l)=zcon/area_e(l)/2 !>0
          enddo !l

!         Bottom b.c.
          ubar=sum(stokes_vel(1,k,elnode(:,i)))/3 !average bottom hvel
          vbar=sum(stokes_vel(2,k,elnode(:,i)))/3
          dhdx=dp(n1)*dldxy(1,1,i)+dp(n2)*dldxy(2,1,i)+dp(n3)*dldxy(3,1,i) !eframe
          dhdy=dp(n1)*dldxy(1,2,i)+dp(n2)*dldxy(2,2,i)+dp(n3)*dldxy(3,2,i)
          stokes_w(kbe(i),i)=-dhdx*ubar-dhdy*vbar

          do l=kbe(i),nvrt-1
            sum1=0
            do j=1,3
              jsj=elside(j,i)
              tmp1=stokes_vel_sd(1,l,jsj)
              tmpx1=stokes_vel_sd(1,l+1,jsj)
              tmp2=stokes_vel_sd(2,l,jsj)
              tmpy2=stokes_vel_sd(2,l+1,jsj)
              if(ics==1) then
                vnor1=tmp1*sframe(1,1,jsj)+tmp2*sframe(2,1,jsj)
                vnor2=tmpx1*sframe(1,1,jsj)+tmpy2*sframe(2,1,jsj)
              else !lat/lon
                !vnor1=su2(l,jsj) !normal
                !vnor2=su2(l+1,jsj)
                call project_hvec(tmp1,tmp2,eframe(:,:,i),sframe(:,:,jsj),vnor1,x1)
                call project_hvec(tmpx1,tmpy2,eframe(:,:,i),sframe(:,:,jsj),vnor2,x1)
              endif !ics
              sum1=sum1+ssign(j,i)*(zs(max(l+1,kbs(jsj)),jsj)-zs(max(l,kbs(jsj)),jsj))*distj(jsj)*(vnor1+vnor2)/2
            enddo !j=1,3

            ubar=sum(stokes_vel(1,l,elnode(:,i)))/3 !level l
            vbar=sum(stokes_vel(2,l,elnode(:,i)))/3
            ubar1=sum(stokes_vel(1,l+1,elnode(:,i)))/3 !level l+1
            vbar1=sum(stokes_vel(2,l+1,elnode(:,i)))/3
!           Impose bottom no-flux b.c.
            if(l==kbe(i)) then
              bflux=0
            else
              bflux=ubar*sne(1,l)+vbar*sne(2,l)+stokes_w(l,i)*sne(3,l)
            endif

            stokes_w(l+1,i)=(-sum1-(ubar1*sne(1,l+1)+vbar1*sne(2,l+1))*area_e(l+1) + &
     &bflux*area_e(l))/sne(3,l+1)/area_e(l+1)
          enddo !l=kbe(i),nvrt-1
        enddo !i=1,nea

!       Convert to node
        stokes_w_nd=0
        do i=1,np !resident only
          if(idry(i)==1) cycle

!         Wet node
          do k=kbp(i),nvrt
            tmp0=0
            do j=1,nne(i)
              ie=indel(j,i)
              stokes_w_nd(k,i)=stokes_w_nd(k,i)+stokes_w(max(k,kbe(ie)),ie)*area(ie)
              tmp0=tmp0+area(ie)
            enddo !j
            stokes_w_nd(k,i)=stokes_w_nd(k,i)/tmp0
          enddo !k
        enddo !i

        call exchange_p3dw(stokes_w_nd)

      endif !RADFLAG.eq.'VOR'

#endif /*USE_WWM*/

!...  compute wind stress components (in lat/lon frame if ics=2; in map projection E-N direction if ics=1)
      dragcmin=1.0d-3*(0.61+0.063*6)
      dragcmax=1.0d-3*(0.61+0.063*50)
      tau=0 !init.
      do i=1,npa
        if(nws==0) then
          tau(1,i)=0
          tau(2,i)=0
        else if(nws==1.or.nws==4.or.nws>=2.and.ihconsv==0.or.iwind_form==-1) then
          wmag=sqrt(windx(i)**2+windy(i)**2)
          dragcoef=1.0d-3*(0.61+0.063*wmag)
          dragcoef=min(max(dragcoef,dragcmin),dragcmax)
          tau(1,i)=dragcoef*0.001293*wmag*windx(i)*rampwind
          tau(2,i)=dragcoef*0.001293*wmag*windy(i)*rampwind
        else !nws>=2 and ihconsv !=0 and iwind_form=0; tauxz and tauyz defined
          if(idry(i)==1) then
            tau(1,i)=0
            tau(2,i)=0
          else !rescale as well
            tau(1,i)=-tauxz(i)/rho0*rampwind*windfactor(i)**2 !sign and scale difference between stresses tauxz and tau
            tau(2,i)=-tauyz(i)/rho0*rampwind*windfactor(i)**2
          endif
        endif !nws
      enddo !i=1,npa

!     Overwrite by WWM values
#ifdef USE_WWM
      if(icou_elfe_wwm>0.and.iwind_form==-2) then 
        do i=1,npa
          if(idry(i)==1) then
            tau(1:2,i)=0
          else
            !stress=rho_air*ufric^2 [Pa]; scaled by rho_water
            tmp=1.293e-3*out_wwm_windpar(i,8)**2*rampwind 
            !Wind direction
            theta=atan2(windy(i),windx(i))
            tau(1,i)=tmp*cos(theta)
            tau(2,i)=tmp*sin(theta)
          endif
        enddo !i
      endif !icou_elfe_wwm
#endif

      if(myrank==0) write(16,*)'done adjusting wind stress ...'

!    VIMS mode added by YC
#ifdef USE_ICM 
      if(myrank==0) write(16,*)'start ICM point source..'
      if(iWQPS==2) then
        if(time>=npstime) then
          npstime=npstime+npstiminc
!org yc         npstime1=npstime2
!org yc         npstime2=npstime2+npstiminc
          PSQ(:)=0.
          PSK(:)=0
          WWPRPOC(:) = 0.
          WWPLPOC(:) = 0.
          WWPDOCA(:) = 0.
          WWPRPON(:) = 0.
          WWPLPON(:) = 0.
          WWPDON(:)  = 0.
          WWPNH4(:)  = 0.
          WWPNO3(:)  = 0.
          WWPRPOP(:) = 0.
          WWPLPOP(:) = 0.
          WWPDOP(:)  = 0.
          WWPPO4t(:) = 0.
          WWPSU(:)   = 0.
          WWPSAt(:)  = 0.
          WWPCOD(:)  = 0.
          WWPDO(:)   = 0.
          x1 = 1.0E3 !* DTD
          read(61,*)      !time
          do i=1,nps
            read(61,*) iegb,xPSK,xPSQ,PRPOC,PLPOC,PDOCA,PRPON,PLPON,PDON, &
           &             PNH4,PNO3,PRPOP,PLPOP,PDOP,PPO4t,PSU,PSAt,PCOD,PDO
            if(iegl(iegb)%rank==myrank) then
              PSQ(iegl(iegb)%id)     = xPSQ
              PSK(iegl(iegb)%id)     = xPSK
              WWPRPOC(iegl(iegb)%id) = PRPOC * x1  ! kg/d * 10^3  = g per day
              WWPLPOC(iegl(iegb)%id) = PLPOC * x1
              WWPDOCA(iegl(iegb)%id) = PDOCA * x1
              WWPRPON(iegl(iegb)%id) = PRPON * x1
              WWPLPON(iegl(iegb)%id) = PLPON * x1
              WWPDON(iegl(iegb)%id)  = PDON  * x1
              WWPNH4(iegl(iegb)%id)  = PNH4  * x1
              WWPNO3(iegl(iegb)%id)  = PNO3  * x1
              WWPRPOP(iegl(iegb)%id) = PRPOP * x1
              WWPLPOP(iegl(iegb)%id) = PLPOP * x1
              WWPDOP(iegl(iegb)%id)  = PDOP  * x1
              WWPPO4t(iegl(iegb)%id) = PPO4t * x1
              WWPSU(iegl(iegb)%id)   = PSU  * x1
              WWPSAt(iegl(iegb)%id)  = PSAt * x1
              WWPCOD(iegl(iegb)%id)  = PCOD * x1
              WWPDO(iegl(iegb)%id)   = PDO  * x1
            endif
          enddo
        endif !time>=npstime+npstiminc
      endif ! iWQPS=2
      if(myrank==0) write(16,*)'end ICM point source..'

!    VIMS surface temperature mode added by YC
      if(myrank==0) write(16,*)'doing ICM surface T..'
      if(iSun==2) then
        if(time>=surf_time2) then
          surf_time1=surf_time2
          surf_time2=surf_time2+86400.
          surf_t1=surf_t2
          read(62,*)
          do i=1,np_global
            read(62,*) ipgb,tmp
            if(ipgl(ipgb)%rank==myrank) then
              surf_t2(ipgl(ipgb)%id)=tmp
            endif
          enddo
        endif !time>=wtime2+idwindrv*86400.

        stratio=(time-surf_time1)/86400.

        do i=1,npa
          surf_t(i)=surf_t1(i)+stratio*(surf_t2(i)-surf_t1(i))
        enddo !i
      endif ! iSun=2
      if(myrank==0) write(16,*)'done ICM surface T..'
#endif /*USE_ICM*/
!End of VIMS surface temperature mode added by YC

!...  Read in temp. and salt for nudging
      if(inu_st==2) then
        if(time>time_nu) then
          irec_nu=irec_nu+1
          time_nu=time_nu+step_nu
          tnd_nu1=tnd_nu2
          snd_nu1=snd_nu2
          read(37)floatout
          read(35)floatout
          if(floatout/=time_nu) then
            write(errmsg,*)'Wrong nudging time:',floatout,time_nu
            call parallel_abort(errmsg)
          endif
          do i=1,np_global
            read(37)(swild8(j,1),j=1,nvrt)
            read(35)(swild8(j,2),j=1,nvrt)
            if(ipgl(i)%rank==myrank) then
              tnd_nu2(:,ipgl(i)%id)=swild8(1:nvrt,1)
              snd_nu2(:,ipgl(i)%id)=swild8(1:nvrt,2)
            endif
          enddo !i
        endif !time>time_nu

!       Compute S,T
        rat=(time_nu-time)/step_nu
        if(rat<0.or.rat>1) then
          write(errmsg,*)'Impossible 81:',rat
          call parallel_abort(errmsg)
        endif
        tnd_nu=tnd_nu1+(1-rat)*(tnd_nu2-tnd_nu1)
        snd_nu=snd_nu1+(1-rat)*(snd_nu2-snd_nu1)
      endif !nudging

!...  Read in tracer nudging
      if(ntracers>0.and.inu_tr==2) then
        if(time>time_nu_tr) then
          irec_nu_tr=irec_nu_tr+1
          time_nu_tr=time_nu_tr+step_nu_tr
          trnd_nu1=trnd_nu2
          read(45)floatout
          if(abs(floatout-time_nu_tr)>0.01) then
            write(errmsg,*)'Wrong nudging time (2):',floatout,time_nu
            call parallel_abort(errmsg)
          endif
          do i=1,np_global
            read(45)swild9
            if(ipgl(i)%rank==myrank) then
              trnd_nu2(:,:,ipgl(i)%id)=swild9
            endif
          enddo !i
        endif !time>time_nu

!       Compute S,T
        rat=(time_nu_tr-time)/step_nu_tr
        if(rat<0.or.rat>1) then
          write(errmsg,*)'Impossible 82:',rat
          call parallel_abort(errmsg)
        endif
        trnd_nu=trnd_nu1+(1-rat)*(trnd_nu2-trnd_nu1)
      endif !nudging

!...  Compute hydraulic transfer blocks together with reading in flux values
!...  in case the blocks are taken out
!...  isblock_nd(1:2,1:npa): this array does not change over time iteration (static).
!                            (1,:) points to the block # of either active or _inactive_ block or 0 (NEVER a part of a block)
!                            (2,:) points to the face # of either active or _inactive_ block or 0
!     isblock_el(1:nea): points to ACTIVE block #; 0 means it's either inactive or not part of a block;
!     isblock_sd(1:2,1:nsa): (1,:) points to ACTIVE block #; 0 means it's either on an 
!                            INACTIVE block or NEVER part of a block;
!                            (2,:) when the block is active, it points to the face # or -1 (inside the block);!                             0 means it's either on an inactive block or never part of a block.
!     q_block(1:nhtblocks): flow from face "1" to "2" at a step. If structures(istruct)%install is false, the block is deactivated; 
!                           otherwise the block is active.

      if(ihydraulics/=0.and.nhtblocks>0) then
        ! reads time varying parameters
        call read_struct_ts(time)

        !Message passing to get elev., vel. info for ref. node #2 for each block
        !do i=1,npa; eta2(i)=iplg(i); enddo !test
        block_refnd2_eta=-1.e6 !init. as flags
        !Post send
        do i=0,nproc-1
          if(nhtsend1(i)/=0) then
            !if(i==myrank) call parallel_abort('MAIN: illegal comm.(2)')
            call mpi_isend(eta2,1,htsend_type(i),i,601,comm,srqst(i),ierr)
            if(ierr/=MPI_SUCCESS) call parallel_abort('MAIN: send error (2)')
!'
          else
            srqst(i)=MPI_REQUEST_NULL
          endif
        enddo !i

        !Post recv
        do i=0,nproc-1
          if(nhtrecv1(i)/=0) then
            !if(i==myrank) call parallel_abort('MAIN: illegal comm.(2)')
            call mpi_irecv(block_refnd2_eta,1,htrecv_type(i),i,601,comm,rrqst(i),ierr)
            if(ierr/=MPI_SUCCESS) call parallel_abort('MAIN: recv error (2)')
!'
          else
            rrqst(i)=MPI_REQUEST_NULL
          endif
        enddo !i

        call mpi_waitall(nproc,rrqst,rstat,ierr)
        if(ierr/=MPI_SUCCESS) call parallel_abort('MAIN: mpi_waitall rrqst tag=601',ierr)
        call mpi_waitall(nproc,srqst,sstat,ierr)
        if(ierr/=MPI_SUCCESS) call parallel_abort('MAIN: mpi_waitall srqst tag=601',ierr)
!'

        !Compute fluxes by proc's that own ref. node #1 (as non-ghost)
        iq_block_lcl=0 !local
        q_block_lcl=0.d0
        do i=1,nhtblocks
          ndgb1=structures(i)%upnode
          if(ipgl(ndgb1)%rank==myrank) then;  if(ipgl(ndgb1)%id<=np) then
            ndgb2=structures(i)%downnode
            irank=ipgl(ndgb2)%rank
            if(irank/=myrank) then
              if(block_refnd2_eta(i)<-1.e6+1) then
                write(errmsg,*)'MAIN: htexchange not rite:',i,ndgb1,ndgb2,irank
                call parallel_abort(errmsg)
              !else
              !  write(12,*)'htex:',i,ndgb1,ndgb2,irank,block_refnd2_eta(i)
              endif
            else !node #2 inside myrank
              block_refnd2_eta(i)=eta2(ipgl(ndgb2)%id)
            endif !irank

            !Compute flux
            call calc_struc_flow(i,eta2(ipgl(ndgb1)%id),block_refnd2_eta(i),q_block_lcl(i))
            iq_block_lcl(i)=iq_block_lcl(i)+1
          endif; endif !ipgl
        enddo !i=1,nhtblocks

        !Broadcast flux to all proc's
        call mpi_allreduce(q_block_lcl,q_block,nhtblocks,rtype,MPI_SUM,comm,ierr)
        call mpi_allreduce(iq_block_lcl,iq_block,nhtblocks,itype,MPI_SUM,comm,ierr)
        do i=1,nhtblocks
          if(iq_block(i)<=0) then
            write(errmsg,*)'MAIN: q_block left out:',i,iq_block(i)
            call parallel_abort(errmsg)
          else
            q_block(i)=q_block(i)/iq_block(i)
          endif
        enddo !i

        !Compute flags for elements, sides on _active_ blocks
        isblock_el=0 !0 or active block #
        do i=1,nea
          jblock=minval(isblock_nd(1,elnode(1:3,i)))
          if(jblock>0) then; if(structures(jblock)%install) then
            isblock_el(i)=jblock
            !write(12,*)'Block elem:',jblock,ielg(i)
          endif; endif
        enddo !i
        !isblock_sd(1,1:nsa): active block #
        !isblock_sd(2,1:nsa): face # or -1 (inside block) or 0
        isblock_sd=0 !init
        !Compute cross sectional areas for 2 faces of each block
        allocate(buf1(nhtblocks,2),buf2(nhtblocks,2))
        buf1=0
        do i=1,nsa
          jblock=minval(isblock_nd(1,isidenode(1:2,i)))
          if(jblock>0) then; if(structures(jblock)%install) then
            n1=isidenode(1,i); n2=isidenode(2,i)
            if(isblock_nd(1,n1)==isblock_nd(1,n2)) then
              isblock_sd(1,i)=jblock !block #
              if(isblock_nd(2,n1)==isblock_nd(2,n2)) then !not internal; on same face
                 jface=isblock_nd(2,n1)
                 isblock_sd(2,i)=jface !face #

                 !Check
                 !write(12,*)'Block face side:',jblock,jface,iplg(isidenode(1:2,i)),i,ns

                 !For resident sides, compute local cross sectional area
                 if(i<=ns) then 
                    !Deal with interface sides
                    ifl=0 !logical flag
                    if(.not.associated(isgl(islg(i))%next)) ifl=1
                    if(associated(isgl(islg(i))%next)) then
                       if(isgl(islg(i))%next%rank>=myrank) ifl=1
                    endif

                    if(ifl==1) then
                       htot=max(h0,dps(i)+(eta2(n1)+eta2(n2))/2)
                       buf1(jblock,jface)=buf1(jblock,jface)+htot*distj(i)
                    endif !ifl
                 endif !i<=ns
              else
                 isblock_sd(2,i)=-1 !internal
              endif
            endif ! isblock_nd(1,n1)==isblock_nd(1,n2)
          endif; endif 
        enddo !i=1,nsa

#ifdef INCLUDE_TIMING
        cwtmp=mpi_wtime()
#endif
        call mpi_allreduce(buf1,buf2,2*nhtblocks,rtype,MPI_SUM,comm,ierr)
#ifdef INCLUDE_TIMING
        wtimer(3,2)=wtimer(3,2)+mpi_wtime()-cwtmp
#endif

        !Check
        !write(12,*)'Block face area:',it,(real(buf2(i,1:2)),i=1,nhtblocks) !,(real(buf1(i,1:2)),i=1,nhtblocks)
        !write(12,*)it,(real(buf2(i,1:2)),i=1,nhtblocks) !,(real(buf1(i,1:2)),i=1,nhtblocks)

        !Compute (uniform) normal vel. at faces for each block
        !Positive is from face '1' to '2' (given in dir_block())
        vnth_block=-99 !flag
        do i=1,nhtblocks
          if(structures(i)%install) then
            !Active block
            do j=1,2 !face
              ar=buf2(i,j)
              if(ar<=0) then
                write(errmsg,*) 'MAIN: Block areas<=0:',i,j,ar,it
                call parallel_abort(errmsg)
              endif
              !Test
              !vnth_block(j,i)=q_block(i)
              
              !add ramp??
              vnth_block(j,i)=q_block(i)/ar !positive from face 1 to 2
            enddo !j; face
          endif !q_block
        enddo !i=1,nhtblocks

        deallocate(buf1,buf2)
      endif !ihydraulics/=0 and nhtblocks>0

!...  Get new time series values from *.th
!...
      if(nettype>0) then
        if(time>th_time(1,2,1)) then !not '>=' to avoid last step
          ath(:,1,1,1)=ath(:,1,2,1)
          read(50,*) tmp,ath(1:nettype,1,2,1)
          th_time(1,1,1)=th_time(1,2,1)
          th_time(1,2,1)=th_time(1,2,1)+th_dt(1,1)
        endif !time
!        if(it==iths_main+1.and.abs(tmp-time)>1.e-4) then
!          write(errmsg,*)'Starting time wrong for eta',it,tmp
!          call parallel_abort(errmsg)
!        endif
      
        rat=(time-th_time(1,1,1))/th_dt(1,1)
        if(rat<-small1.or.rat>1+small1) then
          write(errmsg,*) 'STEP: rat out in elev.th:',rat,time,th_time(1,1:2,1)
          call parallel_abort(errmsg)
        endif
        icount=0
        do k=1,nope_global
          if(iettype(k)==1) then
            icount=icount+1
            if(icount>nettype) call parallel_abort('Wrong counting 1')
            eth(1,k)=(1-rat)*ath(icount,1,1,1)+rat*ath(icount,1,2,1)
          endif
        enddo 
      endif !nettype

      if(nfltype>0) then
        if(time>th_time(1,2,2)) then
          ath(:,1,1,2)=ath(:,1,2,2)
          read(51,*) tmp,ath(1:nfltype,1,2,2)
          th_time(1,1,2)=th_time(1,2,2)
          th_time(1,2,2)=th_time(1,2,2)+th_dt(1,2)
        endif !time
!        if(it==iths_main+1.and.abs(tmp-time)>1.e-4) then
!          write(errmsg,*)'Starting time wrong for flux',it,tmp,time
!          call parallel_abort(errmsg)
!        endif

        rat=(time-th_time(1,1,2))/th_dt(1,2)
        if(rat<-small1.or.rat>1+small1) then
          write(errmsg,*) 'STEP: rat out in flux.th:',rat,time,th_time(1,1:2,2)
          call parallel_abort(errmsg)
        endif
        icount=0
        do k=1,nope_global
          if(ifltype(k)==1) then
            icount=icount+1
            if(icount>nfltype) call parallel_abort('Wrong counting 2')
            qthcon(k)=(1-rat)*ath(icount,1,1,2)+rat*ath(icount,1,2,2)
          endif
        enddo !k
      endif !nfltype

      if(ntetype>0) then
        if(time>th_time(1,2,3)) then
          ath(:,1,1,3)=ath(:,1,2,3)
          read(52,*) tmp,ath(1:ntetype,1,2,3)
          th_time(1,1,3)=th_time(1,2,3)
          th_time(1,2,3)=th_time(1,2,3)+th_dt(1,3)
        endif !time
!        if(it==iths_main+1.and.abs(tmp-time)>1.e-4) then
!          write(errmsg,*)'Starting time wrong for temp',it,tmp
!          call parallel_abort(errmsg)
!        endif

        rat=(time-th_time(1,1,3))/th_dt(1,3)
        if(rat<-small1.or.rat>1+small1) then
          write(errmsg,*) 'STEP: rat out in temp.th:',rat,time,th_time(1,1:2,3)
          call parallel_abort(errmsg)
        endif
        icount=0
        do k=1,nope_global
          if(itetype(k)==1) then
            icount=icount+1
            if(icount>ntetype) call parallel_abort('Wrong counting 3')
            tth(1,1,k)=(1-rat)*ath(icount,1,1,3)+rat*ath(icount,1,2,3)
          endif
        enddo !k
      endif !ntetype

      if(nsatype>0) then
        if(time>th_time(1,2,4)) then
          ath(:,1,1,4)=ath(:,1,2,4)
          read(53,*) tmp,ath(1:nsatype,1,2,4)
          th_time(1,1,4)=th_time(1,2,4)
          th_time(1,2,4)=th_time(1,2,4)+th_dt(1,4)
        endif !time
!        if(it==iths_main+1.and.abs(tmp-time)>1.e-4) then
!          write(errmsg,*)'Starting time wrong for salt',it,tmp
!          call parallel_abort(errmsg)
!        endif

        rat=(time-th_time(1,1,4))/th_dt(1,4)
        if(rat<-small1.or.rat>1+small1) then
          write(errmsg,*) 'STEP: rat out in salt.th:',rat,time,th_time(1,1:2,4)
          call parallel_abort(errmsg)
        endif
        icount=0
        do k=1,nope_global
          if(isatype(k)==1) then
            icount=icount+1
            if(icount>nsatype) call parallel_abort('Wrong counting 4')
            sth(1,1,k)=(1-rat)*ath(icount,1,1,4)+rat*ath(icount,1,2,4)
          endif
        enddo !k
      endif

      if(ntrtype>0) then !type I
        do m=1,ntracers
          if(time>th_time(m,2,5)) then
            ath(:,m,1,5)=ath(:,m,2,5)
            read(300+m,*) tmp,ath(1:ntrtype,m,2,5)
            th_time(m,1,5)=th_time(m,2,5)
            th_time(m,2,5)=th_time(m,2,5)+th_dt(m,5)
          endif !time
!          if(it==iths_main+1.and.abs(tmp-time)>1.e-4) then
!            write(errmsg,*)'Starting time wrong for tracer',it,tmp
!            call parallel_abort(errmsg)
!          endif

          rat=(time-th_time(m,1,5))/th_dt(m,5)
          if(rat<-small1.or.rat>1+small1) then
            write(errmsg,*) 'STEP: rat out in htr_.th:',rat,time,th_time(m,1:2,5)
            call parallel_abort(errmsg)
          endif
          icount=0
          do k=1,nope_global
            if(itrtype(k)==1) then
              icount=icount+1
              if(icount>ntrtype) call parallel_abort('Wrong counting 5')
              trth(m,1,1,k)=(1-rat)*ath(icount,m,1,5)+rat*ath(icount,m,2,5)
            endif
          enddo !k
        enddo !# of tracers
      endif

      if(nettype2>0) then
        if(time>th_time2(2,1)) then        
          ath2(:,:,:,1,1)=ath2(:,:,:,2,1)
          irec_th(1)=irec_th(1)+1
          read(54,rec=irec_th(1)) floatout,ath2(1,1,1:nnode_et,2,1)
          th_time2(1,1)=th_time2(2,1)
          th_time2(2,1)=th_time2(2,1)+th_dt2(1)
        endif !time
!        if(it==iths_main+1.and.abs(floatout-time)>1.e-4) then
!          write(errmsg,*)'Starting time wrong for eta 2',it,floatout
!          call parallel_abort(errmsg)
!        endif

        rat=(time-th_time2(1,1))/th_dt2(1)
        if(rat<-small1.or.rat>1+small1) then
          write(errmsg,*) 'STEP: rat out in elev2D.th:',rat,time,th_time2(1:2,1)
          call parallel_abort(errmsg)
        endif
        icount=0
        icount2=0
        do k=1,nope_global
          if(iettype(k)==4) then
            icount=icount+1
            if(icount>nettype2) call parallel_abort('Wrong counting 7')
            do j=1,nond_global(k)
!              nd=iond_global(k,j)
              icount2=icount2+1
              if(icount2>nnode_et) call parallel_abort('Wrong counting nodes')
!'
              eth(j,k)=(1-rat)*ath2(1,1,icount2,1,1)+rat*ath2(1,1,icount2,2,1)
            enddo !j
          endif
        enddo !k
      endif !nettype2

      if(nfltype2>0) then
        if(time>th_time2(2,2)) then
          ath2(:,:,:,1,2)=ath2(:,:,:,2,2)
          irec_th(2)=irec_th(2)+1
          read(58,rec=irec_th(2)) floatout,ath2(1:2,1:nvrt,1:nnode_fl,2,2)
          th_time2(1,2)=th_time2(2,2)
          th_time2(2,2)=th_time2(2,2)+th_dt2(2)
        endif !time
!        if(it==iths_main+1.and.abs(floatout-time)>1.e-4) then
!          write(errmsg,*)'Starting time wrong for flux 2',it,floatout
!          call parallel_abort(errmsg)
!        endif

        rat=(time-th_time2(1,2))/th_dt2(2)
        if(rat<-small1.or.rat>1+small1) then
          write(errmsg,*) 'STEP: rat out in uv3D.th:',rat,time,th_time2(1:2,2)
          call parallel_abort(errmsg)
        endif
        icount=0
        icount2=0
        do k=1,nope_global
          if(iabs(ifltype(k))==4) then
            icount=icount+1
            if(icount>nfltype2) call parallel_abort('Wrong counting 6')
            do j=1,nond_global(k)
              icount2=icount2+1
              if(icount2>nnode_fl) call parallel_abort('Wrong counting vel')
!'
              uthnd(1:nvrt,j,k)=(1-rat)*ath2(1,1:nvrt,icount2,1,2)+rat*ath2(1,1:nvrt,icount2,2,2) !ll frame if ics=2
              vthnd(1:nvrt,j,k)=(1-rat)*ath2(2,1:nvrt,icount2,1,2)+rat*ath2(2,1:nvrt,icount2,2,2)
            enddo !j
          endif
        enddo !k
      endif !nfltype2

      if(ntetype2>0) then
        if(time>th_time2(2,3)) then
          ath2(:,:,:,1,3)=ath2(:,:,:,2,3)
          irec_th(3)=irec_th(3)+1
          read(56,rec=irec_th(3)) floatout,ath2(1,1:nvrt,1:nnode_te,2,3)
          th_time2(1,3)=th_time2(2,3)
          th_time2(2,3)=th_time2(2,3)+th_dt2(3)
        endif !time
!        if(it==iths_main+1.and.abs(floatout-time)>1.e-4) then
!          write(errmsg,*)'Starting time wrong for temp. 2',it,floatout
!          call parallel_abort(errmsg)
!        endif

        rat=(time-th_time2(1,3))/th_dt2(3)
        if(rat<-small1.or.rat>1+small1) then
          write(errmsg,*) 'STEP: rat out in temp3D.th:',rat,time,th_time2(1:2,3)
          call parallel_abort(errmsg)
        endif
        icount=0
        icount2=0
        do k=1,nope_global
          if(iabs(itetype(k))==4) then
            icount=icount+1
            if(icount>ntetype2) call parallel_abort('Wrong counting 8')
            do j=1,nond_global(k)
!              nd=iond_global(k,j)
              icount2=icount2+1
              if(icount2>nnode_te) call parallel_abort('Wrong counting te')
!'
              tth(1:nvrt,j,k)=(1-rat)*ath2(1,1:nvrt,icount2,1,3)+rat*ath2(1,1:nvrt,icount2,2,3)
            enddo !j
          endif
        enddo !k
      endif !ntetype2

      if(nsatype2>0) then
        if(time>th_time2(2,4)) then
          ath2(:,:,:,1,4)=ath2(:,:,:,2,4)
          irec_th(4)=irec_th(4)+1
          read(57,rec=irec_th(4)) floatout,ath2(1,1:nvrt,1:nnode_sa,2,4)
          th_time2(1,4)=th_time2(2,4)
          th_time2(2,4)=th_time2(2,4)+th_dt2(4)
        endif !time
!        if(it==iths_main+1.and.abs(floatout-time)>1.e-4) then
!          write(errmsg,*)'Starting time wrong for salt 2',it,floatout
!          call parallel_abort(errmsg)
!        endif
 
        rat=(time-th_time2(1,4))/th_dt2(4)
        if(rat<-small1.or.rat>1+small1) then
          write(errmsg,*) 'STEP: rat out in salt3D.th:',rat,time,th_time2(1:2,4)
          call parallel_abort(errmsg)
        endif
        icount=0
        icount2=0
        do k=1,nope_global
          if(iabs(isatype(k))==4) then
            icount=icount+1
            if(icount>nsatype2) call parallel_abort('Wrong counting 9')
            do j=1,nond_global(k)
!              nd=iond_global(k,j)
              icount2=icount2+1
              if(icount2>nnode_sa) call parallel_abort('Wrong counting sa')
!'
              sth(1:nvrt,j,k)=(1-rat)*ath2(1,1:nvrt,icount2,1,4)+rat*ath2(1,1:nvrt,icount2,2,4)
            enddo !j
          endif
        enddo !k
      endif !nsatype2

!     Tracers
      if(ntrtype2>0) then
        if(time>th_time2(2,5)) then
          ath2(:,:,:,1,5)=ath2(:,:,:,2,5)
          irec_th(5)=irec_th(5)+1
          read(59,rec=irec_th(5)) floatout,ath2(1:ntracers,1:nvrt,1:nnode_tr,2,5)
          th_time2(1,5)=th_time2(2,5)
          th_time2(2,5)=th_time2(2,5)+th_dt2(5)
        endif !time
!        if(it==iths_main+1.and.abs(floatout-time)>1.e-4) then
!          write(errmsg,*)'Starting time wrong for tracers 2',it,floatout
!          call parallel_abort(errmsg)
!        endif

        rat=(time-th_time2(1,5))/th_dt2(5)
        if(rat<-small1.or.rat>1+small1) then
          write(errmsg,*) 'STEP: rat out in tr3D.th:',rat,time,th_time2(1:2,5)
          call parallel_abort(errmsg)
        endif
        icount=0
        icount2=0
        do k=1,nope_global
          if(itrtype(k)==4) then
            icount=icount+1
            if(icount>ntrtype2) call parallel_abort('Wrong counting 10')
            do j=1,nond_global(k)
              icount2=icount2+1
              if(icount2>nnode_tr) call parallel_abort('Wrong counting tr')
!'
              trth(1:ntracers,1:nvrt,j,k)=(1-rat)*ath2(1:ntracers,1:nvrt,icount2,1,5)+ &
     &rat*ath2(1:ntracers,1:nvrt,icount2,2,5)
            enddo !j
          endif !itetype
        enddo !k
      endif !ntrtype2

!     Read in volume/mass sources/sinks
      vsource=0 !init; dimension [m^3/s]; includes sinks as well
      if(if_source==1) then
        msource=0 !init; dimension same as concentration (psu etc)
        if(nsources>0) then
          if(time>th_time3(2,1)) then !not '>=' to avoid last step
            ath3(:,1,1,1)=ath3(:,1,2,1)
            read(63,*)tmp,ath3(1:nsources,1,2,1)
            th_time3(1,1)=th_time3(2,1)
            th_time3(2,1)=th_time3(2,1)+th_dt3(1)
          endif !time

          if(time>th_time3(2,3)) then !not '>=' to avoid last step
            ath3(:,:,1,3)=ath3(:,:,2,3)  
            !do j=1,2+ntracers
            read(65,*)tmp,ath3(1:nsources,1:2+ntracers,2,3)
            !enddo !j
            th_time3(1,3)=th_time3(2,3)
            th_time3(2,3)=th_time3(2,3)+th_dt3(3)
          endif !time
 
          rat=(time-th_time3(1,1))/th_dt3(1)
          if(rat<-small1.or.rat>1+small1) then
            write(errmsg,*) 'STEP: rat out in vsource.th:',rat,time,th_time3(1:2,1)
            call parallel_abort(errmsg)
          endif

          do i=1,nsources
            if(ath3(i,1,1,1)<0.or.ath3(i,1,2,1)<0) then
              write(errmsg,*)'STEP: wrong vsource',it,i,ath3(i,1,1:2,1)
              call parallel_abort(errmsg)
            endif

            if(iegl(ieg_source(i))%rank==myrank) then
              ie=iegl(ieg_source(i))%id
              vsource(ie)=((1-rat)*ath3(i,1,1,1)+rat*ath3(i,1,2,1))*ramp_ss
            endif !ielg
          enddo !i

          rat=(time-th_time3(1,3))/th_dt3(3)
          if(rat<-small1.or.rat>1+small1) then
            write(errmsg,*) 'STEP: rat out in msource.th:',rat,time,th_time3(1:2,3)
            call parallel_abort(errmsg)
          endif

          do j=1,2+ntracers
            do i=1,nsources
              if(iegl(ieg_source(i))%rank==myrank) then
                ie=iegl(ieg_source(i))%id
                msource(j,ie)=(1-rat)*ath3(i,j,1,3)+rat*ath3(i,j,2,3) !swild(i)
              endif !ielg
            enddo !i
          enddo !j
        endif !nsources

        if(nsinks>0) then
          if(time>th_time3(2,2)) then !not '>=' to avoid last step
            ath3(:,1,1,2)=ath3(:,1,2,2)
            read(64,*)tmp,ath3(1:nsinks,1,2,2)
            th_time3(1,2)=th_time3(2,2)
            th_time3(2,2)=th_time3(2,2)+th_dt3(2)
          endif !time

          rat=(time-th_time3(1,2))/th_dt3(2)
          if(rat<-small1.or.rat>1+small1) then
            write(errmsg,*) 'STEP: rat out in vsink.th:',rat,time,th_time3(1:2,2)
            call parallel_abort(errmsg)
          endif

          do i=1,nsinks
            if(ath3(i,1,1,2)>0.or.ath3(i,1,2,2)>0) then
              write(errmsg,*)'STEP: wrong vsink',it,i,ath3(i,1,1:2,2)
              call parallel_abort(errmsg)
            endif

            if(iegl(ieg_sink(i))%rank==myrank) then
              ie=iegl(ieg_sink(i))%id
              vsource(ie)=vsource(ie)+((1-rat)*ath3(i,1,1,2)+rat*ath3(i,1,2,2))*ramp_ss
            endif !ielg
          enddo !i
        endif !nsinks
      endif !if_source

!     Calcualtion of cross-section areas for flow b.c.
      if(lflbc) then
        allocate(buf1(nope_global,1),buf2(nope_global,1)); buf1=0d0;
        do k=1,nope
          kk=iopelg(k) !global segment #
          if(ifltype(kk)/=0) then
            do i=1,nond(k)-1
              n1=iond(k,i)
              n2=iond(k,i+1)
              !Find a local side
              isd0=0
              loop01: do j=1,nne(n1)
                ie=indel(j,n1)
                if(ie>0) then
                  do l=1,3
                    isd=elside(l,ie)
                    if((isidenode(1,isd)==n1.or.isidenode(2,isd)==n1).and. &
                       (isidenode(1,isd)==n2.or.isidenode(2,isd)==n2)) then
                       isd0=isd; exit loop01
                    endif
                  enddo !l
                endif !ie>0
              end do loop01 !j=1,nne(n1)

              if(isd0==0.or.isd0>ns) cycle !skip ghost to avoid duplication

              htot=dps(isd0)+(eta2(n1)+eta2(n2))/2
              if(htot<=h0) then
                write(errmsg,*)'Dry bnd side: h_tot',htot, &
     &'open boundary',kk,'node',i,'node index',iplg(n1)
                call parallel_abort(errmsg)
              endif
              buf1(kk,1)=buf1(kk,1)+htot*distj(isd0)
            enddo !i=1,nond(k)-1
          endif
        enddo !k=1,nope

#ifdef INCLUDE_TIMING
        cwtmp=mpi_wtime()
#endif
        call mpi_allreduce(buf1,buf2,nope_global,rtype,MPI_SUM,comm,ierr)
#ifdef INCLUDE_TIMING
        wtimer(3,2)=wtimer(3,2)+mpi_wtime()-cwtmp
#endif
        carea=0
        do k=1,nope_global
          if(ifltype(k)/=0) carea(k)=buf2(k,1)
        enddo
        deallocate(buf1,buf2)
      endif !lflbc

!      if(myrank==8) write(99,*)carea

!...  Compute new vel. for flow b.c. (uth,vth)
!     For ics=1, uth, vth are in global frame
!     For ics=2, they are in lat/lon frame (even at poles) 
      do i=1,nsa
        n1=isidenode(1,i)
        n2=isidenode(2,i)
        ibnd=isbs(i) !global bnd #
        if(ibnd<=0) cycle

!       Open bnds
!       ll frame at side
        swild10(1:3,1:3)=(pframe(:,:,n1)+pframe(:,:,n2))/2

!       Find bnd node indices for n1,n2
        nwild(1:2)=0
        do j=1,2
          do jj=1,2
            if(isbnd(jj,isidenode(j,i))==ibnd) then
              nwild(j)=isbnd(-jj,isidenode(j,i)) !global index
              exit
            endif
          enddo !jj
          if(nwild(j)==0) then
            write(errmsg,*)'Open bnd side has non-bnd node:',i,ibnd,iplg(n1),iplg(n2)
            call parallel_abort(errmsg)
          endif
        enddo !j

        if(ifltype(ibnd)==1.or.ifltype(ibnd)==2) then
          if(carea(ibnd)==0) then
            write(errmsg,*)'Dry bnd side 2',ibnd,carea(ibnd)
            call parallel_abort(errmsg)
          endif
          vnth0=qthcon(ibnd)*ramp/carea(ibnd)
          if(ics==1) then
            uth(:,i)=vnth0*sframe(1,1,i)
            vth(:,i)=vnth0*sframe(2,1,i)
          else !lat/lon
            call project_hvec(vnth0,0.d0,sframe(:,:,i),swild10(1:3,1:3),tmp1,tmp2)
            uth(:,i)=tmp1
            vth(:,i)=tmp2
          endif !ics

        else if(ifltype(ibnd)==-1) then !Flather 1
!         uthnd is the normal vel.; no ramp up
          do k=1,nvrt
            if(uthnd(k,nwild(1),ibnd)<-98.or.uthnd(k,nwild(2),ibnd)<-98) then
              write(errmsg,*)'MAIN: Problem with Flather:',iplg(n1),iplg(n2)
              call parallel_abort(errmsg)
            endif
            tmp=(uthnd(k,nwild(1),ibnd)+uthnd(k,nwild(2),ibnd))/2
            if(ics==1) then
              uth(k,i)=tmp*sframe(1,1,i)
              vth(k,i)=tmp*sframe(2,1,i) 
            else !lat/lon
              call project_hvec(tmp,0.d0,sframe(:,:,i),swild10(1:3,1:3),tmp1,tmp2)
              uth(k,i)=tmp1
              vth(k,i)=tmp2
            endif !ics
          enddo !k

        else if(ifltype(ibnd)==3) then
          vnth0=0 !normal vel.
          do jfr=1,nbfr
            ncyc=int(amig(jfr)*time/2/pi)
            arg=amig(jfr)*time-ncyc*2*pi+face(jfr)-vfa(ibnd,jfr)
            vnth0=vnth0+ramp*ff(jfr)*vmo(ibnd,jfr)*cos(arg)
          enddo !jfr=1,nbfr
          if(ics==1) then
            uth(:,i)=vnth0*sframe(1,1,i)
            vth(:,i)=vnth0*sframe(2,1,i)
          else !lat/lon
            call project_hvec(vnth0,0.d0,sframe(:,:,i),swild10(1:3,1:3),tmp1,tmp2)
            uth(:,i)=tmp1
            vth(:,i)=tmp2
          endif !ics
        else if(iabs(ifltype(ibnd))==4) then
          do k=1,nvrt
            if(uthnd(k,nwild(1),ibnd)<-98.or.uthnd(k,nwild(2),ibnd)<-98.or. &
               vthnd(k,nwild(1),ibnd)<-98.or.vthnd(k,nwild(2),ibnd)<-98) then
              write(errmsg,*)'Wrong time series of vel.'
              call parallel_abort(errmsg)
            endif
            uth(k,i)=ramp*(uthnd(k,nwild(1),ibnd)+uthnd(k,nwild(2),ibnd))/2
            vth(k,i)=ramp*(vthnd(k,nwild(1),ibnd)+vthnd(k,nwild(2),ibnd))/2
          enddo !k
        endif
      enddo !i=1,nsa

      if(myrank==0) write(16,*)'done flow b.c.'

#ifdef INCLUDE_TIMING
!     End forcing preparation section
      wtmp2=mpi_wtime()
      wtimer(3,1)=wtimer(3,1)+wtmp2-wtmp1
!     Start btrack
      wtmp1=wtmp2
#endif


!...  Bottom drag coefficients for nchi=-1 or 1; Cd and Cdp for nchi=0 already read in
      if(nchi==-1) then !2D
        Cdp=0; Cd=0 !for dry pts
!       Drag at nodes
        do i=1,npa
          if(idry(i)==1) cycle
!         Wet node
          htot=max(hmin_man,dp(i)+eta2(i)) !>0
          Cdp(i)=grav*rmanning(i)*rmanning(i)/htot**0.333
#ifdef USE_SED2D
          if(idrag_sed2d<-1) then
            Cdp(i)=Cdsed(i)
            if(Cdp(i)/=Cdp(i)) call parallel_abort('SED2D: NaN for Cd')
          endif
#endif
        enddo !i
      endif !nchi==-1

      if(nchi==1) then !idrag=2; 3D
#ifdef USE_SED
        !Roughness predictor
        if(bedforms_rough.GE.1) THEN
          IF(myrank==0) WRITE(16,*)'start sed_roughness'
          CALL sed_roughness
          IF(myrank==0) WRITE(16,*)'done sed_roughness'
          !Check
          tmp=sum(rough_p)
          if(tmp/=tmp) call parallel_abort('SED3D gave NaN from sed_roughness')
        endif
#endif
        Cdp=0; Cd=0 !for dry pts
        !Cdmax=-1 !max. Cd at node for this process (info only)
!'      Drag at nodes
        ltmp1(1)=.false. !for WBL iteration
        do i=1,npa
          if(idry(i)==1) cycle

!         Wet node
          htot=dp(i)+eta2(i)
          if(rough_p(i)<=0) then !time-independent Cd
            Cdp(i)=abs(rough_p(i))
          else !roughness >0
            bthick_ori=znl(kbp(i)+1,i)-znl(kbp(i),i)  !thickness of bottom bnd layer
            bthick=max(dzb_min,bthick_ori)
            if(bthick<=rough_p(i)) then
              !if(ifort12(5)==0) then
              !  ifort12(5)=1
              !  write(12,*)'BL too fine (2):',i,bthick,rough_p(i),htot
              !endif
              !Cdp(i)=Cdmax
              write(errmsg,*)'MAIN: dzb_min <= roughness at node ',iplg(i),dzb_min,rough_p(i)
              call parallel_abort(errmsg)
            else
              Cdp(i)=1/(2.5*log(bthick/rough_p(i)))**2 

              if(dzb_decay/=0.and.bthick_ori<bthick) then !dzb_decay=0 leads to no decay
                Cdp(i)=Cdp(i)*exp(dzb_decay*(1.-bthick_ori/bthick))
              endif
              !WBL
#ifdef USE_WWM
              if(iwbl==1) then
                vmag=sqrt(uu2(kbp(i)+1,i)**2+vv2(kbp(i)+1,i)**2)
                taubx=Cdp(i)*vmag*uu2(kbp(i)+1,i)
                tauby=Cdp(i)*vmag*vv2(kbp(i)+1,i)
                ubm=out_wwm(i,22) !orbital vel.
                wfr=2*pi/max(0.1,out_wwm(i,12)) !angular freq.; out_wwm is not real*8
                wdir=out_wwm(i,18) !wave direction
                call wbl_GM(taubx,tauby,rough_p(i),ubm,wfr,wdir,z0b,fw,delta_wc,iter,ifl)
                if(ifl==2) ltmp1(1)=.true.                
                if(iter>iwbl_itmax) iwbl_itmax=iter
                !Debug
!                if(it==1000) write(12,*)'WBL:',iplg(i),dp(i),ubm,out_wwm(i,5),wdir, &
!     &rough_p(i),z0b,iter,ifl,delta_wc

                if(bthick<=z0b) then
                  write(errmsg,*)'MAIN: dzb_min <= z0b at node ',iplg(i),dzb_min,z0b
                  call parallel_abort(errmsg)
                else
                  Cdp(i)=1/(2.5*log(bthick/z0b))**2
                  if(dzb_decay/=0.and.bthick_ori<bthick) then
                    Cdp(i)=Cdp(i)*exp(dzb_decay*(1.-bthick_ori/bthick))
                  endif
                endif
              endif !iwbl             
#endif /*USE_WWM*/

            endif !bthick
          endif !rough_p

          !if(Cdp(i)>Cdmax) Cdmax=Cdp(i)
        enddo !i=1,npa
        if(it==iths_main+1) write(12,*)'Cd min/max at 1st step= ',minval(Cdp),maxval(Cdp)
        !write(12,*)'Cd min/max at 1st step= ',minval(Cdp),maxval(Cdp)

!       Output warning for WBL if iteration didn't converge
#ifdef USE_WWM
        if(iwbl==1) then
           call mpi_reduce(ltmp1,ltmp2,1,MPI_LOGICAL,MPI_LOR,0,comm,ierr)
           if(myrank==0.and.ltmp2(1)) write(16,*)'WBL-GM did not converge'
           if(myrank==0) write(16,*)'Cumulative max. for GM iteration for rank 0= ',iwbl_itmax
!'
        endif !iwbl
#endif /*USE_WWM*/
      endif !nchi==1

      if(iabs(nchi)==1) then
        do i=1,nsa
          if(idry_s(i)==1) cycle
          Cd(i)=(Cdp(isidenode(1,i))+Cdp(isidenode(2,i)))/2
        enddo !i

!       Output Cd for first step
!        if(it==iths_main+1) then
!          fdb='Cd_0000'
!          lfdb=len_trim(fdb)
!          write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
!          open(32,file='outputs/'//trim(fdb),status='unknown')
!          !write(32,*)'Drag coefficents for nchi=1 or -1'
!          !write(32,*)nsa
!          do i=1,nsa
!            write(32,'(i6,2e14.6,1x,e11.3)')i,xcj(i),ycj(i),Cd(i)
!          enddo !i=1,ns
!          close(32)
!        endif
      endif !iabs(nchi)==1

!
!************************************************************************
!                                                                       *
!               Turbulence closure schemes                              *
!       Compute turbulence diffusivities dfv, dfh,                      *
!       and in MY-G, also dfq[1,2].                                     *
!                                                                       *
!************************************************************************
!

!...  Scheme 2: Pacanowski and Philander (1981)
      if(itur==2) then
        dfv=0; dfh=0 !for dry nodes
        do i=1,npa
          if(idry(i)==1) cycle
          if(prho(1,i)<-98) then
            write(errmsg,*)'Impossible dry 1'
            call parallel_abort(errmsg)
          endif

!         wet nodes
          if(dp(i)<=h1_pp) then
            vmax=vdmax_pp1
            vmin=vdmin_pp1
            tmin=tdmin_pp1
          else if(dp(i)<h2_pp) then
            vmax=vdmax_pp1+(vdmax_pp2-vdmax_pp1)*(dp(i)-h1_pp)/(h2_pp-h1_pp)
            vmin=vdmin_pp1+(vdmin_pp2-vdmin_pp1)*(dp(i)-h1_pp)/(h2_pp-h1_pp)
            tmin=tdmin_pp1+(tdmin_pp2-tdmin_pp1)*(dp(i)-h1_pp)/(h2_pp-h1_pp)
          else !dps >= h2
            vmax=vdmax_pp2
            vmin=vdmin_pp2
            tmin=tdmin_pp2
          endif

          do k=kbp(i),nvrt
            if(k==kbp(i).or.k==nvrt) then
              drhodz=0
            else
              drhodz=(prho(k+1,i)-prho(k-1,i))/(znl(k+1,i)-znl(k-1,i))
            endif
            bvf=-grav*(drhodz/rho0+grav/1.5e3**2)
            k2=min(k+1,nvrt)
            k1=max(k-1,kbp(i))
            dudz=(uu2(k2,i)-uu2(k1,i))/(znl(k2,i)-znl(k1,i))
            dvdz=(vv2(k2,i)-vv2(k1,i))/(znl(k2,i)-znl(k1,i))
            shear2=max(dudz**2+dvdz**2,1.0e-10_rkind) 
            rich=max(bvf/shear2,0._rkind)

!           vmax >= vmin
            dfv(k,i)=vmax/(1+5*rich)**2+vmin
            dfh(k,i)=dfv(k,i)/(1+5*rich)+tmin
          enddo !k      
        enddo !i=1,npa

        if(myrank==0) write(16,*) 'done turbulence closure (PP)...'
      endif !itur=2

!... Scheme 4: GOTM
!    In GOTM, all turbulence variables are defined at whole levels from bottom to F.S.
!    and mean flow variables at half levels. So the bottom is at level 0 (our kbp), 
!    F.S. is at level nlev (our nvrt).

      if(itur==4) then
#ifdef USE_GOTM
!        if(abs(cde-cmiu0**3)>1.e-4) then
!          write(11,*)'Mismatch in GOTM call:',cde,cmiu0**3
!          stop
!        endif
        if(myrank==0) write(16,*)'cde, cmiu0**3 = ',cde,cmiu0**3

        do j=1,npa
          if(idry(j)==1) then
            dfv(:,j)=diffmin(j)
            dfh(:,j)=diffmin(j)
            cycle
          endif
      
!         Friction velocity: [\niu*|du/dz|]^0.5 (m/s)
          u_taus=sqrt(sqrt(tau(1,j)**2+tau(2,j)**2))
          u_taub=sqrt(Cdp(j)*(uu2(kbp(j)+1,j)**2+vv2(kbp(j)+1,j)**2))
          nlev=nvrt-kbp(j)
          do k=0,nlev 
            klev=k+kbp(j) !kbp <= klev <= nvrt
            if(k/=0) h1d(k)=znl(klev,j)-znl(klev-1,j)
!           Shear frequency squared (1/s^2): (du/dz)^2+(dv/dz)^2 -add
!           vertical
!           Buoyancy frequency squared (1/s^2): -g/\rho0*(d\rho/dz))
            if(k==0.or.k==nlev) then
              if(dfv(klev,j)<=0) then
                write(errmsg,*)'Negative viscosity:',dfv(klev,j),iplg(j),klev
                call parallel_abort(errmsg)
              endif
              if(k==0) then
                SS1d(k)=(u_taub**2/dfv(klev,j))**2
              else
                SS1d(k)=(u_taus**2/dfv(klev,j))**2
              endif
              NN1d(k)=0
            else
              ztmp=znl(klev+1,j)-znl(klev-1,j)
              if(ztmp==0) then
                write(errmsg,*)'Zero layer:',iplg(j),klev
                call parallel_abort(errmsg)
              endif
              SS1d(k)=((uu2(klev+1,j)-uu2(klev-1,j))**2+(vv2(klev+1,j)-vv2(klev-1,j))**2)/ztmp**2
              NN1d(k)=-grav/rho0*(prho(klev+1,j)-prho(klev-1,j))/ztmp
            endif
            tke1d(k)=q2(klev,j)
            L1d(k)=xl(klev,j)
            if(tke1d(k)<0.or.L1d(k)<=0) then
              write(errmsg,*)'Negative tke,mixl:',tke1d(k),L1d(k),iplg(j),klev
              call parallel_abort(errmsg)
            endif
            eps1d(k)=cde*tke1d(k)**1.5/L1d(k) 
            num1d(k)=dfv(klev,j)
            nuh1d(k)=dfh(klev,j)

!           Debug11
!            if(myrank==2.and.iplg(j)==14178.and.it==3253) write(98,*)k,h1d(k),NN1d(k),SS1d(k),eps1d(k), &
!     &num1d(k),nuh1d(k),tke1d(k),L1d(k)

          enddo !k=0,nlev
!          h1d(0)=h1d(1)
          toth=eta2(j)+dp(j)
!         surface and bottom roughness length (m)
          z0s=min(0.1d0,toth/10)
          if(Cdp(j)==0) then
            z0b=0
          else
            z0b=(znl(kbp(j)+1,j)-znl(kbp(j),j))*exp(-0.4/sqrt(Cdp(j)))
          endif

!         Debug11
!          if(myrank==2.and.iplg(j)==14178.and.it==3253) then
!            write(99,*)j,'WOW1'
!            write(98,*)nlev,dt,toth,u_taus,u_taub,z0s,z0b
!          endif

          call do_turbulence(nlev,dt,toth,u_taus,u_taub,z0s,z0b,h1d,NN1d,SS1d)

#ifdef USE_TIMOR
          call flmud(j,dt,rough_p(j),SS1d,NN1d,tke1d,eps1d,L1d,num1d,nuh1d)
#endif /*USE_TIMOR*/

!         Debug11
!          if(myrank==2.and.iplg(j)==14178.and.it==3253) write(98,*)(k,h1d(k),NN1d(k),SS1d(k), &
!     &num1d(k),nuh1d(k),tke1d(k),L1d(k),k=0,nlev)

          q2(kbp(j):nvrt,j) = tke1d(0:nlev)
          xl(kbp(j):nvrt,j) = L1d(0:nlev)
!          eps(i,j,:) = eps1d
          do k=0,nlev
            klev=k+kbp(j)
!           Test if they are NaN or invalid numbers
            if(num1d(k)<0.or.nuh1d(k)<0.or.num1d(k)/=num1d(k).or.nuh1d(k)/=nuh1d(k)) then
              write(errmsg,*)'GOTM: problem with mixing:',num1d(k),nuh1d(k)
              call parallel_abort(errmsg)
            endif
 
#ifdef USE_TIMOR
            !Modify viscosity
            tmp=vts(klev,j)
            if(tmp/=tmp) call parallel_abort('GOTM: vts is NaN from TIMOR')
!'
            if(laddmud_v) num1d(k)=num1d(k)+tmp
#endif /*USE_TIMOR*/

            dfv(klev,j)=min(diffmax(j),num1d(k)+diffmin(j)) 
            dfh(klev,j)=min(diffmax(j),nuh1d(k)+diffmin(j))
          enddo !k
        enddo !j=1,npa
#endif /*USE_GOTM*/
      endif !itur==4
 
!... Scheme 3: Mellor-Yamada-Galperin & Umlauf-Burchard scheme
      if(itur==3) then
!------------------------------------------------------------
!     Debug
!      fdb='MY_0000'
!      lfdb=len_trim(fdb)
!      write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
!      open(32,file=trim(fdb),status='unknown')
!      rewind(32)

      do j=1,npa
        if(idry(j)==1) then
          do k=1,nvrt
            q2(k,j)=q2min; xl(k,j)=xlmin2(j)
            dfv(k,j)=0; dfh(k,j)=0; dfq1(k,j)=0; dfq2(k,j)=0
          enddo 
          cycle
        endif
        if(prho(1,j)<-98) call parallel_abort('Impossible dry 2')

!       Wet node; compute layer thickness etc.
!       Error: use ufg?
        do k=kbp(j)+1,nvrt
          dzz(k)=znl(k,j)-znl(k-1,j)
          dudz=(uu2(k,j)-uu2(k-1,j))/dzz(k)
          dvdz=(vv2(k,j)-vv2(k-1,j))/dzz(k)
          shearbt(k)=dudz**2+dvdz**2 !@ half levels
          rzbt(k)=grav/rho0*(prho(k,j)-prho(k-1,j))/dzz(k)
          q2ha(k)=(q2(k,j)+q2(k-1,j))/2
          xlha(k)=(xl(k,j)+xl(k-1,j))/2

!         Compute c_psi_3
          if(mid.eq.'MY') then
            cpsi3(k)=0.9
          else !GLS models
            if(rzbt(k)>0) then !unstable
              cpsi3(k)=1
            else !stable
              select case(mid)
                case('KL')
                  cpsi3(k)=2.53
                case('KE')
                  cpsi3(k)=-0.52
                case('KW')
                  cpsi3(k)=-0.58
                case('UB')
                  cpsi3(k)=0.1
                case default
                  write(errmsg,*)'Unknown closure model:',mid
                  call parallel_abort(errmsg)
              end select
            endif
          endif

!         Wall proximity function      
          if(mid.eq.'MY'.or.mid.eq.'KL') then
            zctr2=(znl(k,j)+znl(k-1,j))/2
            dists=eta2(j)-zctr2
            distb=zctr2+dp(j)
            if(dists==0.or.distb==0) then
              write(errmsg,*)'Zero in proximity function:',j,k
              call parallel_abort(errmsg)
            endif
            fwall=1+1.33*(xlha(k)/0.4/distb)**2+0.25*(xlha(k)/0.4/dists)**2
            cpsi2p(k)=fwall*cpsi2 !F_wall*cpsi2
          else !other GLS
            cpsi2p(k)=cpsi2
          endif
        enddo !k=kbp(j)+1,nvrt
!        rzbt(kbp(j))=0 !for Galperin's clipping

!        write(90,*)'WOW1',it,j

!	Compute upper bound for xl 
        do k=kbp(j)+1,nvrt
          dists=eta2(j)-znl(k,j)
          distb=znl(k,j)+dp(j)
!          if(k==kbp(j)) then
!            xlmax(k)=max(xlmin2(j),dzz(k+1)*0.4_rkind)
          if(k==nvrt) then
            xlmax(k)=max(xlmin2(j),dzz(k)*0.4_rkind)
          else !internal layers
            xlmax(k)=0.4*min(dists,distb)
          endif
!          xlmax(k)=max(0.4_rkind*min(dists,distb),xlmin2(j)) !can be very small
!          xlmax(k)=0.4*dists*distb/(dps(j)+etam)
!          xlmax(k)=0.4*min(dp(j)+eta2(j),xlmax00)
          if(xlmax(k)<=0) then
            write(errmsg,*)'Dist<0 in MY-G',j,k,eta2(j)+dp(j),dists,distb
            call parallel_abort(errmsg)
          endif
        enddo !k

!	b.c. (computed using values from previous time except wind)
        q2fs=16.6**(2.0/3)*sqrt(tau(1,j)**2+tau(2,j)**2)/2
        q2fs=max(q2fs,q2min)
        q2bot=16.6**(2.0/3)*Cdp(j)*(uu2(kbp(j)+1,j)**2+vv2(kbp(j)+1,j)**2)/2
        q2bot=max(q2bot,q2min)
        xlbot=max(xlmin2(j),min(2.5_rkind,xlsc0(j)*dzz(kbp(j)+1))*0.4_rkind) !"5" to prevent over-mixing

!modif AD :: modification of mixing layer as Delpey et al.
#ifdef USE_WWM
        tmp0=out_wwm(j,1) !Hs
        zsurf=0.2*tmp0
#else
        zsurf=dzz(nvrt)
#endif
        xlfs=max(xlmin2(j),xlsc0(j)*zsurf*0.4_rkind)

!       Debug
!        write(32,*)j,iplg(j),xlmin2(j),xlsc0(j),dzz(nvrt),xlfs
!        write(90,*)'WOW2',it,j

!	Matrix Q
        nqdim=nvrt-kbp(j) !+1
        do k=kbp(j)+1,nvrt
          kin=k-kbp(j) !+1 !row #
          alow(kin)=0
          bdia(kin)=0
          cupp(kin)=0
          rrhs(1,kin)=0
          if(k<nvrt) then
            tmp=(dfq1(k+1,j)+dfq1(k,j))/2*dt/dzz(k+1)
            bdia(kin)=bdia(kin)+dzz(k+1)/3+tmp
            cupp(kin)=cupp(kin)+dzz(k+1)/6-tmp
            rrhs(1,kin)=rrhs(1,kin)+dzz(k+1)/6*(2*q2(k,j)+q2(k+1,j))
            prod=(dfv(k+1,j)+dfv(k,j))/2*shearbt(k+1)
            buoy=(dfh(k+1,j)+dfh(k,j))/2*rzbt(k+1)
            if(prod+buoy>=0) then
              rrhs(1,kin)=rrhs(1,kin)+dt*dzz(k+1)/2*(prod+buoy)
            else
              tmp=dt*dzz(k+1)/6*(prod+buoy)/q2ha(k+1)
              bdia(kin)=bdia(kin)-2*tmp
              cupp(kin)=cupp(kin)-tmp
            endif
            diss=cmiu0**3*sqrt(q2ha(k+1))/xlha(k+1)*dzz(k+1)/6 !diss/k
            bdia(kin)=bdia(kin)+dt*diss*2
            cupp(kin)=cupp(kin)+dt*diss
          endif

          if(k>kbp(j)+1) then
            tmp=(dfq1(k,j)+dfq1(k-1,j))/2*dt/dzz(k)
            bdia(kin)=bdia(kin)+dzz(k)/3+tmp
            alow(kin)=alow(kin)+dzz(k)/6-tmp
            rrhs(1,kin)=rrhs(1,kin)+dzz(k)/6*(2*q2(k,j)+q2(k-1,j))
            prod=(dfv(k,j)+dfv(k-1,j))/2*shearbt(k)
            buoy=(dfh(k,j)+dfh(k-1,j))/2*rzbt(k)
            if(prod+buoy>=0) then
              rrhs(1,kin)=rrhs(1,kin)+dt*dzz(k)/2*(prod+buoy)
            else
              tmp=dt*dzz(k)/6*(prod+buoy)/q2ha(k)
              bdia(kin)=bdia(kin)-2*tmp
              alow(kin)=alow(kin)-tmp
            endif
            diss=cmiu0**3*sqrt(q2ha(k))/xlha(k)*dzz(k)/6 !diss/k
            bdia(kin)=bdia(kin)+dt*diss*2
            alow(kin)=alow(kin)+dt*diss
          endif
        enddo !k=kbp(j)+1,nvrt

!	Soln for q2 at new level
        call tridag(nvrt,100,nqdim,1,alow,bdia,cupp,rrhs,soln,gam)
        q2tmp(nvrt)=q2fs
        !Extrapolate to bottom mainly for diffusivities
        q2tmp(kbp(j):kbp(j)+1)=q2bot
        do k=kbp(j)+2,nvrt-1
          kin=k-kbp(j) !+1
!          if(k==nvrt) then
!            q2tmp(k)=q2fs
!          else if(k==kbp(j)+1) then
!            q2tmp(k)=q2bot
!          else
          q2tmp(k)=max(soln(1,kin),q2min)
!          endif
        enddo !k

!        write(90,*)'WOW4',it,j,(q2tmp(k),k=1,nvrt)
!        do k=1,nvrt
!          write(90,*)'Level ',k,alow(k),bdia(k),cupp(k)
!        enddo 

!	Matrix QL
        do k=kbp(j)+1,nvrt
          kin=k-kbp(j) !+1
          alow(kin)=0
          bdia(kin)=0
          cupp(kin)=0
          rrhs(1,kin)=0
          if(k<nvrt) then
            tmp=(dfq2(k+1,j)+dfq2(k,j))/2*dt/dzz(k+1)
            bdia(kin)=bdia(kin)+dzz(k+1)/3+tmp
            cupp(kin)=cupp(kin)+dzz(k+1)/6-tmp
            psi_n=cmiu0**rpub*q2(k,j)**rmub*xl(k,j)**rnub !psi^n_{j,k}
            psi_n1=cmiu0**rpub*q2(k+1,j)**rmub*xl(k+1,j)**rnub !psi^n_{j,k+1}
            rrhs(1,kin)=rrhs(1,kin)+dzz(k+1)/6*(2*psi_n+psi_n1)
            prod=cpsi1*(dfv(k+1,j)+dfv(k,j))/2*shearbt(k+1)
            buoy=cpsi3(k+1)*(dfh(k+1,j)+dfh(k,j))/2*rzbt(k+1)
            if(prod+buoy>=0) then
              rrhs(1,kin)=rrhs(1,kin)+dt*dzz(k+1)/2*(prod+buoy)*(psi_n+psi_n1)/2/q2ha(k+1)
            else
              tmp=dt*dzz(k+1)/6*(prod+buoy)/q2ha(k+1)
              bdia(kin)=bdia(kin)-2*tmp
              cupp(kin)=cupp(kin)-tmp
            endif
            diss=cpsi2p(k+1)*cmiu0**3*sqrt(q2ha(k+1))/xlha(k+1)*dzz(k+1)/6 !diss/k
            bdia(kin)=bdia(kin)+dt*diss*2
            cupp(kin)=cupp(kin)+dt*diss
          else !k=nvrt
            bdia(kin)=bdia(kin)+0.4*rnub*dt*dfq2(k,j)/xl(k,j)
          endif

          if(k>kbp(j)+1) then 
            tmp=(dfq2(k,j)+dfq2(k-1,j))/2*dt/dzz(k)
            bdia(kin)=bdia(kin)+dzz(k)/3+tmp
            alow(kin)=alow(kin)+dzz(k)/6-tmp
            psi_n=cmiu0**rpub*q2(k,j)**rmub*xl(k,j)**rnub !psi^n_{j,k}
            psi_n1=cmiu0**rpub*q2(k-1,j)**rmub*xl(k-1,j)**rnub !psi^n_{j,k-1}
            rrhs(1,kin)=rrhs(1,kin)+dzz(k)/6*(2*psi_n+psi_n1)
            prod=cpsi1*(dfv(k,j)+dfv(k-1,j))/2*shearbt(k)
            buoy=cpsi3(k)*(dfh(k,j)+dfh(k-1,j))/2*rzbt(k)
            if(prod+buoy>=0) then
              rrhs(1,kin)=rrhs(1,kin)+dt*dzz(k)/2*(prod+buoy)*(psi_n+psi_n1)/2/q2ha(k)
            else
              tmp=dt*dzz(k)/6*(prod+buoy)/q2ha(k)
              bdia(kin)=bdia(kin)-2*tmp
              alow(kin)=alow(kin)-tmp
            endif
            diss=cpsi2p(k)*cmiu0**3*sqrt(q2ha(k))/xlha(k)*dzz(k)/6 !diss/k
            bdia(kin)=bdia(kin)+dt*diss*2
            alow(kin)=alow(kin)+dt*diss
          else !k=kbp(j)+1
            bdia(kin)=bdia(kin)+0.4*rnub*dt*dfq2(k,j)/xl(k,j)
          endif
        enddo !k=kbp(j)+1,nvrt

!        write(90,*)'WOW5',it,j
!        do k=1,nvrt
!          write(90,*)'Level ',k,alow(k),bdia(k),cupp(k)
!        enddo 

!	Soln for q2l and xl at new level
        call tridag(nvrt,100,nqdim,1,alow,bdia,cupp,rrhs,soln,gam)

!        write(90,*)'WOW6',it,j

        do k=kbp(j)+1,nvrt
          kin=k-kbp(j) !+1
          q2l=max(soln(1,kin),psimin)
          if(k==nvrt) then
            xltmp(k)=xlfs
          else if(k==kbp(j)+1) then
            xltmp(k)=xlbot
          else
            xltmp(k)=(q2l*cmiu0**(-rpub)*q2tmp(k)**(-rmub))**(1/rnub)
          endif
!	  Galperin's clipping 
          if(rzbt(k)<0) then
            upper=sqrt(-0.56*q2tmp(k)/rzbt(k))
            xltmp(k)=min(xltmp(k),upper)
          endif
!	  Max. length based on dissipation; xlmin2 prevails
          xl_max=(cmiu0*sqrt(q2tmp(k)))**3/eps_min
          xltmp(k)=max(xlmin2(j),min(xl_max,xltmp(k)))
!	  Impose max. depth limit
          xltmp(k)=max(xlmin2(j),min(xltmp(k),xlmax(k)))

          q2(k,j)=q2tmp(k)
          xl(k,j)=xltmp(k)
          if(q2(k,j)<0) then
            write(errmsg,*)'Negative q2',q2(k,j),xl(k,j)
            call parallel_abort(errmsg)
          endif
        enddo !k=kbp(j)+1,nvrt

!       Extrapolate q2, xl to bottom mainly for diffusivities
        q2(kbp(j),j)=q2(kbp(j)+1,j)
        xl(kbp(j),j)=xl(kbp(j)+1,j)

!       Compute vertical diffusivities at new time
        do k=kbp(j),nvrt
          call asm(j,k,vd,td,qd1,qd2)
          dfv(k,j)=min(diffmax(j),max(diffmin(j),vd))
          dfh(k,j)=min(diffmax(j),max(diffmin(j),td))
          dfq1(k,j)=min(diffmax(j),max(diffmin(j),qd1))
          dfq2(k,j)=min(diffmax(j),max(diffmin(j),qd2))

!         Debug
!          write(90,*)'No. ',k,xl(k,j),dfh(k,j),dfv(k,j),dfq1(k,j),dfq2(k,j)
        enddo !k=kbp(j)+1,nvrt

!       Extend
        do k=1,kbp(j)-1
          q2(k,j)=q2(kbp(j),j)
          xl(k,j)=xl(kbp(j),j)
          dfv(k,j)=dfv(kbp(j),j)
          dfh(k,j)=dfh(kbp(j),j)
          dfq1(k,j)=dfq1(kbp(j),j)
          dfq2(k,j)=dfq2(kbp(j),j)
        enddo !k
      enddo !j=1,npa

!      if(it.eq.1739) write(90,*)'WOW7',it

      if(myrank==0) write(16,*)'done MYG-UB...'

!      close(32)
!------------------------------------------------------------
      endif !itur=3
   
#ifdef INCLUDE_TIMING
!     end turbulence
      wtmp2=mpi_wtime()
      wtimer(5,1)=wtimer(5,1)+wtmp2-wtmp1
!     start prepations
      wtmp1=wtmp2
#endif

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!   Backtracking
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!     Debug: test backtracking alone
      if(ibtrack_test==1) then
        !For first step, generate vertical profiles for T,S
        if(it==iths_main+1) then
          do i=1,npa
            do k=1,nvrt
              tnd(k,i)=tem0(1,i)-20*tanh(5*znl(k,i)/dp(i))
              snd(k,i)=sal0(1,i)-20*tanh(5*znl(k,i)/dp(i))
            enddo !k
          enddo !i
          do i=1,nsa
            t0=tsd(1,i); s0=ssd(1,i)
            do k=1,nvrt
              tsd(k,i)=t0-20*tanh(5*zs(k,i)/dps(i))
              ssd(k,i)=s0-20*tanh(5*zs(k,i)/dps(i))
            enddo !k
          enddo !i
        endif !it       

        eta1=0; eta2=0; we=0
        rot_per=3000 !period
        rot_f=2*pi/rot_per !angular freq.
!        xvel0=-1; yvel0=0.9
        do i=1,nsa
          do k=1,nvrt
            su2(k,i)=-ycj(i)*rot_f !xvel0
            sv2(k,i)=xcj(i)*rot_f
          enddo !k
        enddo !i
        do i=1,nea
          do k=1,nvrt
            do j=1,3
              nd=elnode(j,i)
              ufg(k,i,j)=-ynd(nd)*rot_f
              vfg(k,i,j)=xnd(nd)*rot_f
            enddo !j
          enddo !k
        enddo !i
        do i=1,npa
          do k=1,nvrt
            uu2(k,i)=-ynd(i)*rot_f
            vv2(k,i)=xnd(i)*rot_f
            ww2(k,i)=0 !-1.e-4*znl(k,i)*(50+znl(k,i))
          enddo !k
        enddo !i
      endif !ibtrack_test

!      fdb='btrack_0000'
!      lfdb=len_trim(fdb)
!      write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
!      open(10,file='outputs/'//fdb,status='unknown')

!     temp fix
!      if(ics==2) call zonal_flow

!     For ELM transport, pre-compute cubic spline coefficients
      if(iupwind_t==0) then
        do i=1,npa
          if(idry(i)==1) cycle
          call cubic_spline(nvrt-kbp(i)+1,znl(kbp(i):nvrt,i),tnd(kbp(i):nvrt,i), &
     &0._rkind,0._rkind,cspline_ypp_nd(kbp(i):nvrt,i,1))
          call cubic_spline(nvrt-kbp(i)+1,znl(kbp(i):nvrt,i),snd(kbp(i):nvrt,i), &
     &0._rkind,0._rkind,cspline_ypp_nd(kbp(i):nvrt,i,2))
        enddo !i
        do i=1,nsa
          if(idry_s(i)==1) cycle
          call cubic_spline(nvrt-kbs(i)+1,zs(kbs(i):nvrt,i),tsd(kbs(i):nvrt,i), &
     &0._rkind,0._rkind,cspline_ypp_sd(kbs(i):nvrt,i,1))
          call cubic_spline(nvrt-kbs(i)+1,zs(kbs(i):nvrt,i),ssd(kbs(i):nvrt,i), &
     &0._rkind,0._rkind,cspline_ypp_sd(kbs(i):nvrt,i,2))
        enddo !i
      endif !iupwind_t=0      

!...  From nodes and sidecenters, and whole levels
!...  ptbt, sdbt: interpolated values at whole levels
!...  webt: for vertical vel. at element (for non-hydrostatic model only)
!...  ptbt(:,:,:) for T,S only
!...  For ics=2, sdbt(1:2,:,:) are vel. vector in the sframe of the originating side,
!...  and all coordinates are expressed in the local frame (pframe or sframe) at originating
!...  node or side.

!     Pre-assign for dry and below-bottom nodes/sides.
      do i=1,np
        ptbt(1,:,i)=tnd(:,i)
        ptbt(2,:,i)=snd(:,i)
        ptbt(3,:,i)=dfv(:,i)
        ptbt(4,:,i)=dfh(:,i)
      enddo !i

      do i=1,ns
        sdbt(1,:,i)=su2(:,i)
        sdbt(2,:,i)=sv2(:,i)
        sdbt(3,:,i)=tsd(:,i)
        sdbt(4,:,i)=ssd(:,i)
      enddo !i

      webt=we

!     For vortex formulation, temporarily alter w-vel (will be restored after
!     btrack)
#ifdef USE_WWM
      if(RADFLAG.eq.'VOR') ww2=ww2+stokes_w_nd
#endif

!     Initialize inter-subdomain backtracking count
      nbtrk=0

!     Do for nodes or sides or elements; elements are only for non-hydro model
!      p_dis_max=-1 !max. error for node
!      s_dis_max=-1 !max. error for side
!      p_vdis_max=-1 !max. error for node (vel)
!      s_vdis_max=-1 !max. error for side (vel)
!      p_T_max=-1 !max. error for T
!      s_T_max=-1 !max. error for T
      do l=1,3
!       Resident nodes/sides only
        if(l==1) then !nodes
          limit=np
        else if(l==2) then
          limit=ns
        else
          limit=ne
        endif

        do i=1,limit
!         Bypass nodes if upwind scheme is used for both S,T
          if(l==1.and.iupwind_t/=0) cycle
!         Bypass elements for hydrostatic mode
          if(l==3.and.nonhydro==0) cycle

          if(l==1) then !nodes
            if(idry(i)==1) cycle
!           Wet node
            nd0=i
            jmin=kbp(nd0)
!   	    Find the 1st wet element (for consistency)
            ie0=0
            do m=1,nne(nd0)
              ie=indel(m,nd0)
              if(idry_e(ie)==0) then
                ie0=ie; exit
              endif
            enddo !m
            if(ie0==0) then
              write(errmsg,*)'MAIN: btrack finds no init. element (1):',iplg(nd0)
              call parallel_abort(errmsg)
            endif
            !swild to store global coord. of the starting node for ics=2 (not used for ics=1)
            swild(1)=xnd(nd0); swild(2)=ynd(nd0); swild(3)=znd(nd0)
            !swild10 to store pframe at starting pt for ics=2 (not used for ics=1)
            swild10(1:3,1:3)=pframe(:,:,nd0)
          else if(l==2) then !sides
            if(idry_s(i)==1) cycle
            isd0=i
            jmin=kbs(isd0)
            ie0=0
            do m=1,2
              ie=isdel(m,isd0)
              if(ie/=0) then; if(idry_e(ie)==0) then
                ie0=ie; exit
              endif; endif
            enddo !m
            if(ie0==0) then
              write(errmsg,*)'MAIN: btrack finds no init. element (2):',iplg(isidenode(1:2,isd0))
              call parallel_abort(errmsg)
            endif
            !swild to store global coord. of the starting side
            swild(1)=xcj(isd0); swild(2)=ycj(isd0); swild(3)=zcj(isd0)
            !swild10 to store sframe at starting pt for ics=2 (not used for ics=1)
            swild10(1:3,1:3)=sframe(:,:,isd0)
          else !elements; no lat/lon
            if(idry_e(i)==1) cycle
            ie0=i
            jmin=kbe(i)
            !swild to store global coord. of the starting element, although this is not used currently
            swild(1)=xctr(ie0); swild(2)=yctr(ie0); swild(3)=zctr(ie0)
            !swild10 to store eframe at starting pt, although this is not used currently
            swild10(1:3,1:3)=eframe(:,:,ie0)
          endif

          do j=jmin,nvrt 
!           Initialize (xt,yt,zt),nnel and vel.
!           For ics=2, the coord. are in local frames
!	    Caution! nnel must be initialized inside this loop as it is updated inside.
            if(l==1) then !nodes
              ipsgb=iplg(nd0)
              iadvf=iadv(nd0)
              if(ics==1) then
                xt=xnd(nd0)
                yt=ynd(nd0)
              else !lat/lon; in node frame
                xt=0 
                yt=0
                !centroid coord. for nudging
                call project_pt('g2l',xctr(ie0),yctr(ie0),zctr(ie0), &
     &(/xnd(nd0),ynd(nd0),znd(nd0)/),pframe(:,:,nd0),xctr2,yctr2,tmp)
              endif !ics
              zt=znl(j,nd0)
              uuint=uu2(j,nd0) !node frame for ics=2
              vvint=vv2(j,nd0)
              wwint=ww2(j,nd0)
              if(isbnd(1,nd0)/=0) then !on land or open bnd
                ifl_bnd=1
              else
                ifl_bnd=0
              endif
            else if(l==2) then !sides
              ipsgb=islg(isd0)
              n1=isidenode(1,isd0)
              n2=isidenode(2,isd0)
              iadvf=min(iadv(n1),iadv(n2))
              if(ics==1) then
                xt=xcj(isd0)
                yt=ycj(isd0)
              else !lat/lon; in side frame
                xt=0
                yt=0
                !centroid coord. for nudging
                call project_pt('g2l',xctr(ie0),yctr(ie0),zctr(ie0), &
     &(/xcj(isd0),ycj(isd0),zcj(isd0)/),sframe(:,:,isd0),xctr2,yctr2,tmp)
              endif !ics
              zt=zs(j,isd0)
              uuint=su2(j,isd0) !in side frame for ics=2
              vvint=sv2(j,isd0)
              wwint=(ww2(j,n1)+ww2(j,n2))/2 !in side frame for ics=2 (same vertical direction)
              if(isbs(isd0)/=0) then !on land or open bnd
                ifl_bnd=1
              else
                ifl_bnd=0
              endif
            else !elements; no lat/lon
              ipsgb=ielg(ie0)
              iadvf=min(iadv(elnode(1,i)),iadv(elnode(2,i)),iadv(elnode(3,i)))
              xt=xctr(ie0)
              yt=yctr(ie0)
              zt=ze(j,ie0)
              uuint=(su2(j,elside(1,ie0))+su2(j,elside(2,ie0))+su2(j,elside(3,ie0)))/3
              vvint=(sv2(j,elside(1,ie0))+sv2(j,elside(2,ie0))+sv2(j,elside(3,ie0)))/3
              wwint=we(j,ie0)
            endif !l
            !vmag=sqrt(uuint**2+vvint**2+wwint**2)
            vmag=sqrt(uuint**2+vvint**2)
            nnel=ie0
            jlev=j
!            jlev=min(j+1,nvrt) !make sure j>=2 for division()

!           vis_coe: blending factor between continuous and discontinuous vel
            if(indvel==1) then
              vis_coe=0
            else !indvel<=0
              if(l==1.or.l==3) then
                vis_coe=vis_coe2
              else !sides
                if(isdel(2,isd0)==0) then
                  vis_coe=vis_coe2
                else
                  vis_coe=vis_coe1
                endif
              endif
            endif

            if(vmag<=velmin_btrack) then !No activity 
              if(l==1) then
                ptbt(1,j,nd0)=tnd(j,nd0)
                ptbt(2,j,nd0)=snd(j,nd0)
                ptbt(3,j,nd0)=dfv(j,nd0)
                ptbt(4,j,nd0)=dfh(j,nd0)
              else if(l==2) then !sides
                sdbt(3,j,isd0)=tsd(j,isd0)
                sdbt(4,j,isd0)=ssd(j,isd0)
                sdbt(1,j,isd0)=su2(j,isd0)
                sdbt(2,j,isd0)=sv2(j,isd0)
              else !element
                webt(j,ie0)=we(j,ie0)
              endif
            else !do btrack
!              if(nadv>0) then
!                dtb_max=dtb_max1
!              else if(iadvf<=1) then !nadv=0
!                dtb_max=dtb_max1
!              else
!                dtb_max=dtb_max2
!              endif

!             Compute # of sub-division based on local gradients 
!              if(iadvf<=1) then !Euler
              ndelt_max=max(1.d0,dt/dtb_min)
              ndelt_min=max(1.d0,dt/dtb_max) 
              if(l==1) then !nodes
                suma=0
                do ii=1,nne(nd0)
                  ie=indel(ii,nd0)
                  id=iself(ii,nd0)
                  dudx=0; dudy=0; dvdx=0; dvdy=0
                  do jj=1,3
                    dudx=dudx+ufg(j,ie,jj)*dldxy(jj,1,ie) !not strictly along z; in element frame for ics=2
                    dudy=dudy+ufg(j,ie,jj)*dldxy(jj,2,ie) !not strictly along z
                    dvdx=dvdx+vfg(j,ie,jj)*dldxy(jj,1,ie) 
                    dvdy=dvdy+vfg(j,ie,jj)*dldxy(jj,2,ie) 
                  enddo !jj
                  suma=suma+dt*sqrt(dudx**2+dudy**2+dvdx**2+dvdy**2)/nne(nd0)
                enddo !ii=1,nne(nd0)
                ndelt=max0(ndelt_min,min0(ndelt_max,int(suma)*4)) !>=1
              else if(l==2) then !sides
                suma=0
                icount=0
                do ii=1,2
                  ie=isdel(ii,isd0)
                  if(ie==0) cycle
                  icount=icount+1

                  dudx=0; dudy=0; dvdx=0; dvdy=0
                  do jj=1,3
                    dudx=dudx+ufg(j,ie,jj)*dldxy(jj,1,ie) !not strictly along z; in element frame for ics=2
                    dudy=dudy+ufg(j,ie,jj)*dldxy(jj,2,ie) !not strictly along z
                    dvdx=dvdx+vfg(j,ie,jj)*dldxy(jj,1,ie)
                    dvdy=dvdy+vfg(j,ie,jj)*dldxy(jj,2,ie)
                  enddo !jj
                  suma=suma+dt*sqrt(dudx**2+dudy**2+dvdx**2+dvdy**2)
                enddo !ii=1,2
                if(icount==0) then
                  write(errmsg,*)'Impossible 77'
                  call parallel_abort(errmsg)
                endif
                ndelt=max0(ndelt_min,min0(ndelt_max,int(suma/icount)*4)) !>=1
              else !element; no lat/lon
                dudx=0; dudy=0; dvdx=0; dvdy=0
                do jj=1,3
                  dudx=dudx+ufg(j,ie0,jj)*dldxy(jj,1,ie0) !not strictly along z
                  dudy=dudy+ufg(j,ie0,jj)*dldxy(jj,2,ie0) !not strictly along z
                  dvdx=dvdx+vfg(j,ie0,jj)*dldxy(jj,1,ie0)
                  dvdy=dvdy+vfg(j,ie0,jj)*dldxy(jj,2,ie0)
                enddo !jj
                suma=dt*sqrt(dudx**2+dudy**2+dvdx**2+dvdy**2)
                ndelt=max0(ndelt_min,min0(ndelt_max,int(suma)*4)) !>=1
              endif
              dtbk=dt/ndelt !target btrack step; may be smaller sometimes
!              endif !iadvf<=1

!             Perturb starting pt to avoid underflow (except for wvel)
              if(l<=2) then
                eps=btrack_nudge
                if(ics==1) then
                  xt=(1-eps)*xt+eps*xctr(nnel)
                  yt=(1-eps)*yt+eps*yctr(nnel)
                else !lat/lon
                  xt=(1-eps)*xt+eps*xctr2
                  yt=(1-eps)*yt+eps*yctr2
                endif !ics
              endif

              time_rm=dt
              time_rm2=-99 !leftover from previous subdomain; init. as flag
              call btrack(l,ipsgb,ifl_bnd,j,iadvf,swild(1:3),swild10(1:3,1:3), &
     &dtbk,vis_coe,time_rm,time_rm2,uuint,vvint,wwint,nnel,jlev,xt,yt,zt,swild3,ltmp)

              if(ltmp) then !Backtracking exits augmented subdomain
                !Add trajectory to inter-subdomain backtracking list
                nbtrk=nbtrk+1
                if(nbtrk>mxnbt) call parallel_abort('MAIN: nbtrk > mxnbt')
!'
                btlist(nbtrk)%rank=myrank
                btlist(nbtrk)%l0=l
                btlist(nbtrk)%i0gb=ipsgb
                btlist(nbtrk)%isbndy=ifl_bnd
                btlist(nbtrk)%j0=j
                btlist(nbtrk)%adv=iadvf
!                btlist(nbtrk)%ndt=ndelt
                btlist(nbtrk)%dtbk=dtbk !dtb_max
                btlist(nbtrk)%vis=vis_coe
                btlist(nbtrk)%rt=time_rm
                btlist(nbtrk)%rt2=time_rm2
                btlist(nbtrk)%ut=uuint
                btlist(nbtrk)%vt=vvint
                btlist(nbtrk)%wt=wwint
                btlist(nbtrk)%iegb=ielg(nnel)
                btlist(nbtrk)%jvrt=jlev
                btlist(nbtrk)%xt=xt
                btlist(nbtrk)%yt=yt
                btlist(nbtrk)%zt=zt
                btlist(nbtrk)%gcor0=swild(1:3)
                btlist(nbtrk)%frame0=swild10(1:3,1:3)
              else !Backtracking completed within augmented subdomain
                if(l==1) then
                  ptbt(1:4,j,nd0)=swild3(1:4) !ttint
                  !ptbt(2,j,nd0)=ssint
                  !ptbt(3,j,nd0)=dfvint
                  !ptbt(4,j,nd0)=dfhint

                  !Check for ics=2 and zonal flow
!                  if(1==2.and.ics==2.and.j==nvrt) then
!                    call project_pt('l2g',xt,yt,0.d0,swild(1:3),swild10(1:3,1:3),xt4,yt4,zt4)
!                    !Exact
!                    !coorind. and lat/lon in the rotated frame
!                    hatx=xnd(nd0)*cos(alpha_zonal)+znd(nd0)*sin(alpha_zonal)
!                    hatz=-xnd(nd0)*sin(alpha_zonal)+znd(nd0)*cos(alpha_zonal)
!                    call compute_ll(hatx,ynd(nd0),hatz,rlon,rlat)
!                    rlam=rlon-omega_zonal*dt !btrack
!                    hatxex=rearth*cos(rlam)*cos(rlat)
!                    hatyex=rearth*sin(rlam)*cos(rlat)
!                    hatzex=rearth*sin(rlat)
!                    !coord. in the original frame
!                    xex=hatxex*cos(alpha_zonal)-hatzex*sin(alpha_zonal)
!                    yex=hatyex
!                    zex=hatxex*sin(alpha_zonal)+hatzex*cos(alpha_zonal)
!                    dis=sqrt((xt4-xex)**2+(yt4-yex)**2+(zt4-zex)**2)
!                    !rotated ll frame at foot; WARINING: assume nvrt>=3
!                    swild2(1,1)=-sin(rlam)*cos(alpha_zonal)
!                    swild2(2,1)=cos(rlam)
!                    swild2(3,1)=-sin(rlam)*sin(alpha_zonal)
!                    swild2(1,2)=-cos(rlam)*sin(rlat)*cos(alpha_zonal)-cos(rlat)*sin(alpha_zonal)
!                    swild2(2,2)=-sin(rlam)*sin(rlat)
!                    swild2(3,2)=-cos(rlam)*sin(rlat)*sin(alpha_zonal)+cos(rlat)*cos(alpha_zonal)
!                    call cross_product(swild2(1,1),swild2(2,1),swild2(3,1), &
!     &swild2(1,2),swild2(2,2),swild2(3,2),swild2(1,3),swild2(2,3),swild2(3,3))
!                    call project_hvec(uuint,vvint,swild10(1:3,1:3),swild2(1:3,1:3),u2,v2)
!                    uzonal=u00_zonal*cos(rlat)
!                    vdis=dsqrt((u2-uzonal)**2+v2*v2)
!                    write(12,*)'Node ',iplg(nd0),j,xt4,yt4,zt4,xex,yex,zex,dis,xnd(nd0),ynd(nd0),znd(nd0)
!                    write(12,*)'Node vel ',j,u2,v2,uzonal,0,vdis
!                    !Exact T (cosine_bell test)
!                    !ll of foot (original frame)
!                    call compute_ll(xex,yex,zex,rlon0,rlat0)
!                    rrr=rearth*acos(cos(rlat0)*cos(rlon0+pi/2))
!                    if(rrr<rearth/3) then
!                      tex=500*(1+cos(pi*rrr/rearth*3))
!                    else
!                      tex=0
!                    endif
!                    ter=swild3(1)-tex
!                    write(12,*)'Node T:',iplg(nd0),jlev,swild3(1),tex,ter,swild3(2)

                    !const. zonal flow
                    !xex=-uu2(nvrt,nd0)*dt !in pframe
                    !yex=-vv2(nvrt,nd0)*dt !in pframe
                    !dis=sqrt((xt-xex)**2+(yt-yex)**2)
                    !call project_hvec(uuint,vvint,pframe(:,:,nd0),pframe(:,:,nd0),utmp,vtmp)
                    !vdis=sqrt((utmp-uzonal)**2+(vtmp-vmer)**2)
                    !write(12,*)'Node level:',iplg(nd0),j,jlev !,xt,yt,xex,yex,dis,utmp,vtmp,vdis
 
!                    if(dis>p_dis_max) p_dis_max=dis
!                    if(vdis>p_vdis_max) p_vdis_max=vdis
!                    if(abs(ter)>p_T_max) p_T_max=abs(ter)
!                  endif !zonal flow

                else if(l==2) then !sides
                  sdbt(3,j,isd0)=swild3(1) !ttint
                  sdbt(4,j,isd0)=swild3(2) !ssint
                  if(iadvf==0) then
                    sdbt(1,j,isd0)=su2(j,isd0)
                    sdbt(2,j,isd0)=sv2(j,isd0)
                  else
                    sdbt(1,j,isd0)=uuint
                    sdbt(2,j,isd0)=vvint
                  endif

                  !Check for ics=2 and zonal flow
!                  if(1==2.and.ics==2.and.j==nvrt) then
!                    n1=isidenode(1,isd0)
!                    n2=isidenode(2,isd0)
!                    call project_pt('l2g',xt,yt,0.d0,swild(1:3),swild10(1:3,1:3),xt4,yt4,zt4)
!                    !Exact
!                    !coorind. and lat/lon in the rotated frame
!                    hatx=xcj(isd0)*cos(alpha_zonal)+zcj(isd0)*sin(alpha_zonal)
!                    hatz=-xcj(isd0)*sin(alpha_zonal)+zcj(isd0)*cos(alpha_zonal)
!                    call compute_ll(hatx,ycj(isd0),hatz,rlam,rlat)
!                    rlam=rlam-omega_zonal*dt
!                    !Have to reduce radius b/c initially sidecenter is not on earth surface
!                    rr0=sqrt(xcj(isd0)**2+ycj(isd0)**2+zcj(isd0)**2)
!                    hatxex=rr0*cos(rlam)*cos(rlat)
!                    hatyex=rr0*sin(rlam)*cos(rlat)
!                    hatzex=rr0*sin(rlat)
!                    !coord. in the original frame
!                    xex=hatxex*cos(alpha_zonal)-hatzex*sin(alpha_zonal)
!                    yex=hatyex
!                    zex=hatxex*sin(alpha_zonal)+hatzex*cos(alpha_zonal)
!                    dis=sqrt((xt4-xex)**2+(yt4-yex)**2+(zt4-zex)**2)
!                    !rotated ll frame at foot; WARINING: assume nvrt>=3
!                    swild2(1,1)=-sin(rlam)*cos(alpha_zonal)
!                    swild2(2,1)=cos(rlam)
!                    swild2(3,1)=-sin(rlam)*sin(alpha_zonal)
!                    swild2(1,2)=-cos(rlam)*sin(rlat)*cos(alpha_zonal)-cos(rlat)*sin(alpha_zonal)
!                    swild2(2,2)=-sin(rlam)*sin(rlat)
!                    swild2(3,2)=-cos(rlam)*sin(rlat)*sin(alpha_zonal)+cos(rlat)*cos(alpha_zonal)
!                    call cross_product(swild2(1,1),swild2(2,1),swild2(3,1), &
!     &swild2(1,2),swild2(2,2),swild2(3,2),swild2(1,3),swild2(2,3),swild2(3,3))
!                    call project_hvec(uuint,vvint,swild10(1:3,1:3),swild2(1:3,1:3),u2,v2)
!                    uzonal=u00_zonal*cos(rlat)
!                    vdis=dsqrt((u2-uzonal)**2+v2*v2)
!                    write(12,*)'Side ',iplg(isidenode(:,isd0)),j,xt4,yt4,zt4,xex,yex,zex,dis,&
!     &xcj(isd0),ycj(isd0),zcj(isd0),rr0
!                    write(12,*)'Side vel ',j,u2,v2,uzonal,0,vdis
!                    !Exact T
!                    !ll of foot (original frame)
!                    call compute_ll(xex,yex,zex,rlon0,rlat0)
!                    rrr=rearth*acos(cos(rlat0)*cos(rlon0+pi/2))
!                    if(rrr<rearth/3) then
!                      tex=500*(1+cos(pi*rrr/rearth*3))
!                    else
!                      tex=0
!                    endif
!                    ter=swild3(1)-tex
!                    write(12,*)'Side T:',jlev,swild3(1),tex,ter,swild3(2)
!                    if(dis>s_dis_max) s_dis_max=dis
!                    if(vdis>s_vdis_max) s_vdis_max=vdis
!                    if(abs(ter)>s_T_max) s_T_max=abs(ter)
!                  endif !zonal flow

                else !element
                  if(ics==2) call parallel_abort('MAIN: why am I here?')
                  if(iadvf==0) then
                    webt(j,ie0)=we(j,ie0)
                  else
                    webt(j,ie0)=wwint
                  endif
                endif 
              endif !ltmp
            endif !do backtrack

!           Debug
!            if(l<=3) then
!              xyzp(nd0,j,1)=xt; xyzp(nd0,j,2)=yt; xyzp(nd0,j,3)=zt;
!            else
!              xyzs(isd0,j,1)=xt; xyzs(isd0,j,2)=yt; xyzs(isd0,j,3)=zt;
!            endif

          enddo !j=jmin,nvrt

        enddo !i=1,limit
      enddo !l=1,3; nodes or sides or elements

!     Complete inter-subdomain backtracking (if necessary)
      if(nproc>1) then
        lbt(1)=(nbtrk/=0)
#ifdef INCLUDE_TIMING
        cwtmp=mpi_wtime()
#endif
        call mpi_allreduce(lbt,lbtgb,1,MPI_LOGICAL,MPI_LOR,comm,ierr)
#ifdef INCLUDE_TIMING
        wtimer(4,2)=wtimer(4,2)+mpi_wtime()-cwtmp
#endif
        if(ierr/=MPI_SUCCESS) call parallel_abort('MAIN: allreduce lbtgb',ierr)
!'

        if(lbtgb(1)) then
          if(myrank==0) write(16,*)'starting inter-subdomain btrack'
          call inter_btrack(it,nbtrk,btlist) !all ranks participate
          if(myrank==0) write(16,*)'done inter-subdomain btrack'
        endif

        if(lbt(1)) then !handle returned inter-subdomain trajectories
          do ibt=1,nbtrk
            if(btlist(ibt)%rank/=myrank) call parallel_abort('MAIN: not right rank')
!'
            l=btlist(ibt)%l0; iadvf=btlist(ibt)%adv
            j=btlist(ibt)%j0
            if(l==1) then !node
              if(ipgl(btlist(ibt)%i0gb)%rank/=myrank) then
                write(errmsg,*)'MAIN: not my node:',ipgl(btlist(ibt)%i0gb)%rank,l,btlist(ibt)%i0gb, &
                btlist(ibt)%j0,btlist(ibt)%adv,btlist(ibt)%iegb,btlist(ibt)%jvrt, &
                btlist(ibt)%vis,btlist(ibt)%rt,btlist(ibt)%ut,btlist(ibt)%vt,btlist(ibt)%wt,btlist(ibt)%sclr(1:4)
                call parallel_abort(errmsg)
              endif
              nd0=ipgl(btlist(ibt)%i0gb)%id
              ptbt(1:4,j,nd0)=btlist(ibt)%sclr(1:4) !tt
              !ptbt(2,j,nd0)=btlist(ibt)%st
              !ptbt(3,j,nd0)=btlist(ibt)%dfv
              !ptbt(4,j,nd0)=btlist(ibt)%dfh

            else if(l==2) then !sides
              if(isgl(btlist(ibt)%i0gb)%rank/=myrank) then
                write(errmsg,*)'MAIN: not my side:',isgl(btlist(ibt)%i0gb)%rank,l,btlist(ibt)%i0gb,&
                btlist(ibt)%j0,btlist(ibt)%adv,btlist(ibt)%iegb,btlist(ibt)%jvrt, &
                btlist(ibt)%vis,btlist(ibt)%rt,btlist(ibt)%ut,btlist(ibt)%vt,btlist(ibt)%wt,btlist(ibt)%sclr(1:4)
                call parallel_abort(errmsg)
              endif
              isd0=isgl(btlist(ibt)%i0gb)%id
              sdbt(3,j,isd0)=btlist(ibt)%sclr(1) !tt
              sdbt(4,j,isd0)=btlist(ibt)%sclr(2) !st
              if(iadvf==0) then
                sdbt(1,j,isd0)=su2(j,isd0)
                sdbt(2,j,isd0)=sv2(j,isd0)
              else
                sdbt(1,j,isd0)=btlist(ibt)%ut
                sdbt(2,j,isd0)=btlist(ibt)%vt
              endif

!              xyzs(isd0,j,1)=btlist(ibt)%xt; xyzs(isd0,j,2)=btlist(ibt)%yt; xyzs(isd0,j,3)=btlist(ibt)%zt;

            else if(l==3) then !element
              if(iegl(btlist(ibt)%i0gb)%rank/=myrank) then
                write(errmsg,*)'MAIN: not my element:',iegl(btlist(ibt)%i0gb)%rank,l,btlist(ibt)%i0gb,&
                btlist(ibt)%j0,btlist(ibt)%adv,btlist(ibt)%iegb,btlist(ibt)%jvrt, &
                btlist(ibt)%vis,btlist(ibt)%rt,btlist(ibt)%ut,btlist(ibt)%vt,btlist(ibt)%wt,btlist(ibt)%sclr(1:4)
                call parallel_abort(errmsg)
              endif
              ie0=iegl(btlist(ibt)%i0gb)%id
              if(iadvf==0) then
                webt(j,ie0)=we(j,ie0)
              else
                webt(j,ie0)=btlist(ibt)%wt
              endif
            else
              call parallel_abort('MAIN: interbtrack node/side/element index wrong')
!'
            endif !sides
          enddo !ibt
        endif !lbt(1)
      endif !nproc>1

!     Update ghost backtracked momentum
#ifdef INCLUDE_TIMING
      cwtmp=mpi_wtime()
#endif
      call exchange_p3d_4(ptbt)
      call exchange_s3d_4(sdbt)
      call exchange_e3dw(webt)

#ifdef INCLUDE_TIMING
      wtimer(4,2)=wtimer(4,2)+mpi_wtime()-cwtmp
#endif

!...  Update viscosity/diffusivity
!      do i=1,npa
!        dfv(:,i)=ptbt(3,:,i)
!        dfh(:,i)=ptbt(4,:,i)
!      enddo !i      

!     Debug
!      do i=1,np
!        th=pi/2+2*pi/3000*time
!        x0=1.8e3*cos(th)
!        y0=1.8e3*sin(th)
!        do k=1,nvrt
!          prho(k,i)=exp(-((xnd(i)-x0)**2+(ynd(i)-y0)**2)/2/600/600) !exact soln
!        enddo !k
!      enddo !i

!      write(12,*)'Max. node errors=',p_dis_max,p_vdis_max,p_T_max
!      write(12,*)'Max. side errors=',s_dis_max,s_vdis_max,s_T_max

!      etmax=-1 !max. T error
!      esmax=-1 !max. S error
!      do i=1,npa
!        !Exact soln
!        !in ll frame
!        xtmp=-uzonal*dt
!        ytmp=-vmer*dt
!        call project_pt('l2g',xtmp,ytmp,0.d0,(/xnd(i),ynd(i),znd(i)/),pframe(:,:,i),xg,yg,zg)
!        call compute_ll(xg,yg,zg,rlon,rlat)
!        rlon=rlon/pi*180
!        rlat=rlat/pi*180
!        tex=max(tempmin,min(tempmax,rlon+164+rlat-33))
!        sex=max(saltmin,min(saltmax,rlon+164-rlat+33))
!        ter=ptbt(1,1,i)-tex
!        ser=ptbt(2,1,i)-sex
!        write(12,*)'Node',i,iplg(i)
!        write(12,*)'T: ',ptbt(1,1,i),tex,ter
!        write(12,*)'S: ',ptbt(2,1,i),sex,ser
!        if(abs(ter)>etmax) etmax=abs(ter)
!        if(abs(ser)>esmax) esmax=abs(ser)
!      enddo !i
!      write(12,*)'Node max. T&S error:',etmax,esmax
!
!      etmax=-1 !max. T error
!      esmax=-1 !max. S error
!      do i=1,nsa
!        n1=isidenode(1,i)
!        n2=isidenode(2,i)
!        !Exact soln
!        !in ll frame
!        xtmp=-uzonal*dt
!        ytmp=-vmer*dt
!        swild10(1:3,1:3)=(pframe(:,:,n1)+pframe(:,:,n2))/2
!        call project_pt('l2g',xtmp,ytmp,0.d0,(/xcj(i),ycj(i),zcj(i)/),swild10(1:3,1:3),xg,yg,zg)
!        call compute_ll(xg,yg,zg,rlon,rlat)
!        rlon=rlon/pi*180
!        rlat=rlat/pi*180
!        tex=max(tempmin,min(tempmax,rlon+164+rlat-33))
!        sex=max(saltmin,min(saltmax,rlon+164-rlat+33))
!        ter=sdbt(3,1,i)-tex
!        ser=sdbt(4,1,i)-sex
!        write(12,*)'Side',i,iplg(n1),iplg(n2)
!        write(12,*)'T: ',sdbt(3,1,i),tex,ter
!        write(12,*)'S: ',sdbt(4,1,i),sex,ser
!        if(abs(ter)>etmax) etmax=abs(ter)
!        if(abs(ser)>esmax) esmax=abs(ser)
!      enddo !i
!      write(12,*)'Side max. T&S error:',etmax,esmax

!     Side
!      do i=1,ns
!        write(10,*)'Side',i,iplg(isidenode(1:2,i))
!
!       Exact soln
!!        r0=sqrt(xcj(i)**2+ycj(i)**2)
!!        if(r0==0) then
!!          x0=0; y0=0
!!        else
!!          th0=atan2(ycj(i),xcj(i))
!!          th=th0-2*pi/rot_per*time
!!          x0=r0*cos(th)
!!          y0=r0*sin(th)
!!        endif
!        x0=xcj(i)-xvel0*dt; y0=ycj(i)-yvel0*dt
!
!        do k=kbs(i),nvrt
!          if(abs(xyzs(i,k,1)-x0)+abs(xyzs(i,k,2)-y0)>difm) then
!            difm=abs(xyzs(i,k,1)-x0)+abs(xyzs(i,k,2)-y0) 
!            in1=i; in2=2
!          endif
!          write(10,*)k,xyzs(i,k,1:2),x0,y0
!        enddo !k
!!        if(abs(xyzs(i,nvrt,1)*10/3.e4+xyzs(i,nvrt,2)/6.e3-sdbt(3,nvrt,i))>0.1) &
!!        write(10,*)sdbt(3,nvrt,i),xyzs(i,nvrt,1)*10/3.e4+xyzs(i,nvrt,2)/6.e3
!      enddo !i
!      write(10,*)'Max diff=',difm,' at node/side ',in1,' which is a node/side',in2
!!'
!      close(10)
!
!      call parallel_finalize
!      stop
!     End debug


!...  bubt: total integrated value
      bubt=0
      do i=1,nea
        do j=1,3 !sides
          isd=elside(j,i)
          if(idry_s(isd)==0) then
            do k=kbs(isd)+1,nvrt !layer
              if(ics==1) then
                bubt(1,i)=bubt(1,i)+(sdbt(1,k,isd)+sdbt(1,k-1,isd))/2*(zs(k,isd)-zs(k-1,isd))*area(i)/3
                bubt(2,i)=bubt(2,i)+(sdbt(2,k,isd)+sdbt(2,k-1,isd))/2*(zs(k,isd)-zs(k-1,isd))*area(i)/3
              else 
                call project_hvec(sdbt(1,k,isd),sdbt(2,k,isd),sframe(:,:,isd),eframe(:,:,i),u2,v2)
                call project_hvec(sdbt(1,k-1,isd),sdbt(2,k-1,isd),sframe(:,:,isd),eframe(:,:,i),u1,v1)
                bubt(1,i)=bubt(1,i)+(u2+u1)/2*(zs(k,isd)-zs(k-1,isd))*area(i)/3 !in eframe
                bubt(2,i)=bubt(2,i)+(v2+v1)/2*(zs(k,isd)-zs(k-1,isd))*area(i)/3 !in eframe
              endif !ics
            enddo !k
          endif
        enddo !j
      enddo !i=1,nea

!     Restore w-vel
#ifdef USE_WWM
      if(RADFLAG.eq.'VOR') ww2=ww2-stokes_w_nd
#endif

      if(myrank==0) write(16,*)'done backtracking'

#ifdef INCLUDE_TIMING
!     End timing first backtracking section
      wtmp2=mpi_wtime()
      wtimer(4,1)=wtimer(4,1)+wtmp2-wtmp1
!     start turbulence timing
      wtmp1=wtmp2
#endif


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Wave continuity equation: preparation of matrix
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!...  compute elevation essential boundary conditions
!...  in case of border node on >1 bnd with imposed elevation, the bnd with largest segment # prevails.
      elbc=-9999 !flags
      do i=1,nope_global
        do j=1,nond_global(i)
          nd=iond_global(i,j) !global
          if(ipgl(nd)%rank==myrank) then
            ip=ipgl(nd)%id
            if(nramp_elev==1) then
              eic=etaic(ip)
            else
              eic=0
            endif

            if(iettype(i)==1.or.iettype(i)==2) then
              elbc(ip)=ramp*eth(1,i)+(1-ramp)*eic
            else if(iettype(i)==3) then
              elbc(ip)=(1-ramp)*eic !initialize
              do jfr=1,nbfr
                ncyc=int(amig(jfr)*time/2/pi)
                arg=amig(jfr)*time-ncyc*2*pi+face(jfr)-efa(i,j,jfr)
                elbc(ip)=elbc(ip)+ramp*ff(jfr)*emo(i,j,jfr)*cos(arg)
              enddo !jfr=1,nbfr
            else if(iettype(i)==4) then
              elbc(ip)=ramp*eth(j,i)+(1-ramp)*eic
            endif

            !Add inverse barometric effects
            if (inv_atm_bnd==1) elbc(ip)=elbc(ip)+ramp*(prmsl_ref-pr(ip))/grav/rho0

          endif !ipgl(nd)%rank==myrank
        enddo !j
      enddo !i=1,nope_global

!     Compute b.c. flag for all nodes for the matrix
      do i=1,npa
        if(elbc(i)>-9998) then
          lelbc(i)=.true.
        else
          lelbc(i)=.false.
        endif
      enddo !i

!...  Pre-compute some arrays: chi,hhat,bigu,ghat1
!...
      do i=1,nsa
        if(idry_s(i)==1) then
          chi(i)=0
          hhat(i)=0
          bigu(:,i)=0
          cycle
        endif

!	    Wet side
        n1=isidenode(1,i)
        n2=isidenode(2,i)
        if(idrag==1) then
          chi(i)=Cd(i)
        else
          chi(i)=Cd(i)*sqrt(sdbt(1,kbs(i)+1,i)**2+sdbt(2,kbs(i)+1,i)**2)
        endif

        if(lm2d) then
          hhat(i)=(eta2(n1)+eta2(n2))/2+dps(i)+chi(i)*dt
          if(hhat(i)<=0) then
            write(errmsg,*)'Impossible dry 53:',hhat(i),iplg(isidenode(1:2,i))
            call parallel_abort(errmsg)
          endif
        else
          hhat(i)=(eta2(n1)+eta2(n2))/2+dps(i)-chi(i)*dt
!	  Enforce positivity for 3D model
          if(ihhat==1) hhat(i)=max(0._rkind,hhat(i))
        endif

!	bigu1,2 (in sframe if ics=2)
        bigu(1,i)=0 !U^n_x
        bigu(2,i)=0 !U^n_y
        do k=kbs(i),nvrt-1
          bigu(1,i)=bigu(1,i)+(zs(k+1,i)-zs(k,i))*(su2(k,i)+su2(k+1,i))/2
          bigu(2,i)=bigu(2,i)+(zs(k+1,i)-zs(k,i))*(sv2(k,i)+sv2(k+1,i))/2
        enddo !k
      enddo !i=1,nsa

      if(myrank==0) write(16,*)'done 1st preparation'
        
!...  Baroclinic force at side and whole levels
      if(ibc==0) then
        bcc=0 
        if(iupwind_t==0) then
!==============================================================
!         ELM option
!==============================================================
          !Density mean profile (rho_mean) removed
          !dr_dxy in eframe if ics=2
          !Error: may not work for ics=2
          call hgrad_nodes(1,0,nvrt,npa,nsa,prho(:,1:npa)-rho_mean(:,1:npa),dr_dxy)  

!         bcc: -g/rho0* \int_z^\eta dr_dxy dz; trapzoidal rule (in sframe if ics=2)
!         ramp-up factor included
          do i=1,ns
            if(idry_s(i)==1) cycle
            bcc(1:2,nvrt,i)=0
            do k=nvrt-1,kbs(i),-1
              bcc(1:2,k,i)=bcc(1:2,k+1,i)-rampbc*grav/rho0*(zs(k+1,i)-zs(k,i))*(dr_dxy(1:2,k+1,i)+dr_dxy(1:2,k,i))/2
            enddo !k=kbs(i),nvrt
          enddo !i=1,ns

        else
!==============================================================
!         iupwind_t/=0; upwind or TVD
!==============================================================
!         Prepare cubic spline (2nd derivative stored in hp_int temporarily)
          hp_int=0 !temporary save of 2nd deriavtives
          do i=1,nea
            if(idry_e(i)==1) cycle

!           Density mean profile (rho_mean) removed
            swild(kbe(i)+1:nvrt)=(ze(kbe(i):nvrt-1,i)+ze(kbe(i)+1:nvrt,i))/2
            call cubic_spline(nvrt-kbe(i),swild(kbe(i)+1:nvrt),erho(kbe(i)+1:nvrt,i)-rho_mean(kbe(i)+1:nvrt,i), &
     &0._rkind,0._rkind,hp_int(kbe(i)+1:nvrt,i,1))
          enddo !i=1,nea

          dr_dxy=0 !@ half levels; in eframe if ics=2
          do i=1,ne !resident
            if(idry_e(i)==1) cycle

            swild(kbe(i)+1:nvrt)=(ze(kbe(i):nvrt-1,i)+ze(kbe(i)+1:nvrt,i))/2
!           Wet element; interpolate 3 neighbors
            do j=1,3 !neighbors
              ie=ic3(j,i)
              if(ie<0) then
                call parallel_abort('MAIN: bcc neighbor outside')
              else if(ie/=0) then; if(idry_e(ie)==0) then !internal and wet
                swild2(kbe(ie)+1:nvrt,4)=(ze(kbe(ie):nvrt-1,ie)+ze(kbe(ie)+1:nvrt,ie))/2
                eta_min=min(swild(nvrt),swild2(nvrt,4))
                zmax=max(swild(kbe(i)+1),swild2(kbe(ie)+1,4))
                if(-zmax>h_bcc1) then !deep sea
                  ibot_fl=0
                else !shallow
                  ibot_fl=1
                endif

                call eval_cubic_spline(nvrt-kbe(ie),swild2(kbe(ie)+1:nvrt,4), &
     &erho(kbe(ie)+1:nvrt,ie)-rho_mean(kbe(ie)+1:nvrt,ie), &
     &hp_int(kbe(ie)+1:nvrt,ie,1),nvrt-kbe(i),swild(kbe(i)+1:nvrt),ibot_fl,zmax,eta_min,swild2(kbe(i)+1:nvrt,j))
              endif; endif
            enddo !j=1,3

            do k=kbe(i)+1,nvrt
!             Maxtrix of 3 eqs.
              do j=1,3 !eqs
                ie=ic3(j,i)
                n1=elnode(nx(j,1),i)
                n2=elnode(nx(j,2),i)
                !max() added to avoid seg fault
                if(ie==0.or.ie/=0.and.idry_e(max(1,ie))==1) then
                  if(ics==1) then
                    xn1=xnd(n1)
                    yn1=ynd(n1)
                    xn2=xnd(n2)
                    yn2=ynd(n2)
                  else !to eframe
                    call project_pt('g2l',xnd(n1),ynd(n1),znd(n1),(/xctr(i),yctr(i),zctr(i)/), &
     &eframe(:,:,i),xn1,yn1,tmp)
                    call project_pt('g2l',xnd(n2),ynd(n2),znd(n2),(/xctr(i),yctr(i),zctr(i)/), &
     &eframe(:,:,i),xn2,yn2,tmp)
                  endif !ics
                  alow(j)=yn2-yn1 !ynd(n2)-ynd(n1)
                  bdia(j)=xn1-xn2 !xnd(n1)-xnd(n2)
                  cupp(j)=0
                else !internal and wet
                  if(ics==1) then
                    alow(j)=xctr(ie)-xctr(i)
                    bdia(j)=yctr(ie)-yctr(i)
                  else !to eframe
                    call project_pt('g2l',xctr(ie),yctr(ie),zctr(ie),(/xctr(i),yctr(i),zctr(i)/), &
     &eframe(:,:,i),xctr2,yctr2,tmp)
                    alow(j)=xctr2
                    bdia(j)=yctr2
                  endif !ics
                  cupp(j)=swild2(k,j)-(erho(k,i)-rho_mean(k,i))

#ifdef USE_TIMOR
                  !Limit density difference
                  cupp(j)=max(-80.d0,min(80.d0,cupp(j)))
#endif
                endif !ie/=0 etc
              enddo !j=1,3

!             Density gradient - average of 3
              icount=0
              do j=1,3 !pairs
                x10=alow(j); y10=bdia(j); bb1=cupp(j)
                x20=alow(nx(j,1)); y20=bdia(nx(j,1)); bb2=cupp(nx(j,1))
                rl10=sqrt(x10*x10+y10*y10)
                rl20=sqrt(x20*x20+y20*y20)
                delta=x10*y20-x20*y10
!                if(delta==0) then
!                  write(errmsg,*)'MAIN: baroc. failure (2):',ielg(i),j
!                  call parallel_abort(errmsg)
!                endif
                if(rl10==0.or.rl20==0) then
                  write(errmsg,*)'MAIN: baroc. failure (2):',ielg(i),j,k
                  call parallel_abort(errmsg)
                endif
                sintheta=abs(delta)/rl10/rl20
                if(sintheta>sin(pi/180)) then !use 1 degree as threshold
                  icount=icount+1
                  swild10(icount,1)=(y20*bb1-y10*bb2)/delta
                  swild10(icount,2)=(x10*bb2-x20*bb1)/delta
                endif
              enddo !j
              if(icount==0) then
                write(errmsg,*)'MAIN: baroc. failure (3):',ielg(i),k
                call parallel_abort(errmsg)
              endif
              dr_dxy(1,k,i)=sum(swild10(1:icount,1))/icount
              dr_dxy(2,k,i)=sum(swild10(1:icount,2))/icount
            enddo !k=kbe(i)+1,nvrt
          enddo !i=1,ne

#ifdef INCLUDE_TIMING
          cwtmp=mpi_wtime()
#endif
          call exchange_e3d_2(dr_dxy)
#ifdef INCLUDE_TIMING
          wtimer(6,2)=wtimer(6,2)+mpi_wtime()-cwtmp
#endif

!         Density gradient at sides and half levels
          do i=1,ns
            if(idry_s(i)==1) cycle
        
            swild2=0 !gradient at sidecenter and half level; sframe if ics=2
            do k=kbs(i)+1,nvrt
              icount=0
              do j=1,2
                ie=isdel(j,i) 
                if(ie==0.or.idry_e(max(1,ie))==1) cycle

!               Wet element
                icount=icount+1
                gam(kbe(ie)+1:nvrt)=(ze(kbe(ie):nvrt-1,ie)+ze(kbe(ie)+1:nvrt,ie))/2
                !dr_dxy at elements and half levels; eframe if ics=2
                rrhs(1,kbe(ie)+1:nvrt)=dr_dxy(1,kbe(ie)+1:nvrt,ie)
                rrhs(2,kbe(ie)+1:nvrt)=dr_dxy(2,kbe(ie)+1:nvrt,ie)
                call vinter(100,nvrt,2,(zs(k,i)+zs(k-1,i))/2,kbe(ie)+1,nvrt,k,gam,rrhs,swild,ibelow)
                if(ics==2) then !to sframe
                  call project_hvec(swild(1),swild(2),eframe(:,:,ie),sframe(:,:,i),tmp1,tmp2)
                  swild(1)=tmp1
                  swild(2)=tmp2
                endif !ics
                swild2(k,1:2)=swild2(k,1:2)+swild(1:2)
              enddo !j
              if(icount==0) call parallel_abort('MAIN: impossible 101')
              swild2(k,1:2)=swild2(k,1:2)/icount
            enddo !k=kbs(i)+1,nvrt

!           bcc (whole levels): -g/rho0* \int_z^\eta dr_dxy dz; trapzoidal rule
!           ramp-up factor included
!           In sframe if ics=2 
            bcc(1:2,nvrt,i)=0
            do k=nvrt-1,kbs(i),-1
              bcc(1:2,k,i)=bcc(1:2,k+1,i)-rampbc*grav/rho0*(zs(k+1,i)-zs(k,i))*swild2(k+1,1:2)
            enddo !k

          enddo !i=1,ns
!==============================================================
        endif !iupwind_t

!       Debug
!        if(myrank==0) then
!          do i=1,ne
!            if(idry_e(i)==1) cycle
!            write(98,*)'Element:',i
!            write(98,'(3(1x,e12.4))')((ze(k,i)+ze(k-1,i))/2,dr_dxy(1:2,k,i),k=kbe(i)+1,nvrt)
!          enddo !i
!        endif
#ifdef INCLUDE_TIMING
        cwtmp=mpi_wtime()
#endif
        call exchange_s3d_2(bcc)
#ifdef INCLUDE_TIMING
        wtimer(6,2)=wtimer(6,2)+mpi_wtime()-cwtmp
#endif
      endif !ibc==0

!     Non-hydrostatic pressure gradient
      if(nonhydro==1) then
        call hgrad_nodes(1,0,nvrt,npa,nsa,qnon,dqnon_dxy)
!       Exchange 
#ifdef INCLUDE_TIMING
        cwtmp=mpi_wtime()
#endif
        call exchange_s3d_2(dqnon_dxy)
#ifdef INCLUDE_TIMING
        wtimer(6,2)=wtimer(6,2)+mpi_wtime()-cwtmp
#endif
      endif !nonhydro==1

!     Debug
!      if(myrank==0) then
!        do i=1,ns
!          if(idry_s(i)==1) cycle
!          write(97,*)'Side:',i
!          write(97,'(3(1x,e12.4))')(zs(k,i),bcc(1:2,k,i),k=kbs(i),nvrt)
!        enddo !i
!      endif
!      call parallel_finalize
!      stop

!     ghat1 (in eframe if ics=2)
#ifdef USE_WWM
#ifdef DEBUG
      wafo = 0.d0
#endif
#endif
      do i=1,nea
        if(ihydraulics/=0.and.nhtblocks>0) then; if(isblock_el(i)>0) then !active block
          ghat1(1,i)=0
          ghat1(2,i)=0
          cycle
        endif; endif

        if(idry_e(i)==1) then
          ghat1(1,i)=0
          ghat1(2,i)=0
          cycle
        endif

!	    Wet elements
!	    Excluding hvis, baroclinc force first
!       Warning: \hat{G}_1 must include all: Coriolis, atmo. pressure, 
!                tidal potential, horizontal difusion, and baroclinic etc
!       Remember to update both f (botf) and F (bigf)
!       If ics=2, all *d{x,y} and ghat1 are in eframe
        tau_x=0
        tau_y=0
        detadx=0
        detady=0
        dprdx=0
        dprdy=0
        detpdx=0
        detpdy=0
        chigamma=0
        ubstar=0 !bottom advection * \chi
        vbstar=0
        hhat_bar=0
        h_bar=0
        bigf1=0 !all in F except baroclinic and hvis; in eframe
        bigf2=0
        botf1=0 !all in \chi*f_b except baroclinic and hvis; in eframe
        botf2=0 
        do j=1,3 !node or side
          nd=elnode(j,i)
!         idry_e(i) checked already
          detadx=detadx+eta2(nd)*dldxy(j,1,i) !in eframe if ics=2
          detady=detady+eta2(nd)*dldxy(j,2,i)
          dprdx=dprdx+pr(nd)*dldxy(j,1,i)
          dprdy=dprdy+pr(nd)*dldxy(j,2,i)
          if(dpe(i)>=tip_dp) then
            detpdx=detpdx+etp(nd)*dldxy(j,1,i)
            detpdy=detpdy+etp(nd)*dldxy(j,2,i)
          endif
          h_bar=h_bar+(dp(nd)+eta2(nd))/3

          isd=elside(j,i)
          chigamma=chigamma+chi(isd)/3
          hhat_bar=hhat_bar+hhat(isd)/3
          if(ics==1) then
            tau_x=tau_x+tau(1,nd)/3
            tau_y=tau_y+tau(2,nd)/3
            ubstar=ubstar+chi(isd)*sdbt(1,kbs(isd)+1,isd)/3 
            vbstar=vbstar+chi(isd)*sdbt(2,kbs(isd)+1,isd)/3
            bigf1=bigf1+cori(isd)*bigu(2,isd)/3
            bigf2=bigf2-cori(isd)*bigu(1,isd)/3
            botf1=botf1+chi(isd)*cori(isd)*sv2(kbs(isd)+1,isd)/3
            botf2=botf2-chi(isd)*cori(isd)*su2(kbs(isd)+1,isd)/3
          else !lat/lon; convert to eframe 
            call project_hvec(tau(1,nd),tau(2,nd),pframe(:,:,nd),eframe(:,:,i),taux2,tauy2)
            tau_x=tau_x+taux2/3 !in eframe
            tau_y=tau_y+tauy2/3
            call project_hvec(sdbt(1,kbs(isd)+1,isd),sdbt(2,kbs(isd)+1,isd),sframe(:,:,isd),eframe(:,:,i),ub2,vb2)
            call project_hvec(bigu(1,isd),bigu(2,isd),sframe(:,:,isd),eframe(:,:,i),bigu2,bigv2)
            call project_hvec(su2(kbs(isd)+1,isd),sv2(kbs(isd)+1,isd),sframe(:,:,isd),eframe(:,:,i),u2,v2)
            ubstar=ubstar+chi(isd)*ub2/3 !\bar{ki}*\bar{u}_b^\star
            vbstar=vbstar+chi(isd)*vb2/3
            bigf1=bigf1+cori(isd)*bigv2/3 !\bar{F}_r
            bigf2=bigf2-cori(isd)*bigu2/3
            botf1=botf1+chi(isd)*cori(isd)*v2/3 !\bar{ki}*\bar{f}_r
            botf2=botf2-chi(isd)*cori(isd)*u2/3
          endif !ics
        enddo !j=1,3
      
        bigf1=bigf1+h_bar*(0.69*grav*detpdx-dprdx/rho0)
        bigf2=bigf2+h_bar*(0.69*grav*detpdy-dprdy/rho0)
        botf1=botf1+chigamma*(0.69*grav*detpdx-dprdx/rho0)
        botf2=botf2+chigamma*(0.69*grav*detpdy-dprdy/rho0)

        if(.not.lm2d) then !3D
          ghat1(1,i)=bubt(1,i)+area(i)*dt*(bigf1+tau_x-ubstar-dt*botf1-grav*(1-thetai)*hhat_bar*detadx)
          ghat1(2,i)=bubt(2,i)+area(i)*dt*(bigf2+tau_y-vbstar-dt*botf2-grav*(1-thetai)*hhat_bar*detady)
        else !2D
          !Compute element averages (eframe)
          av_elem_x=0
          av_elem_y=0
          do j=1,3 !wet side
            isd=elside(j,i) 
            n1=isidenode(1,isd); n2=isidenode(2,isd)
!            etam=(eta2(n1)+eta2(n2))/2
            htot=zs(nvrt,isd)-zs(kbs(isd),isd) !etam+dps(isd)
            if(hhat(isd)<=0.or.htot<=0) then
              write(errmsg,*)'Impossible dry 51:',htot,hhat(isd),ielg(i),iplg(isidenode(1:2,isd))
              call parallel_abort(errmsg)
            endif
            !\hat{Gamma} (in eframe if ics=2)
            if(ics==1) then
              sdbtu=sdbt(1,1,isd)
              sdbtv=sdbt(2,1,isd)
              u2=su2(1,isd)
              v2=sv2(1,isd)
              taux2=(tau(1,n1)+tau(1,n2))/2
              tauy2=(tau(2,n1)+tau(2,n2))/2
            else !lat/lon; project to eframe
              call project_hvec(sdbt(1,1,isd),sdbt(2,1,isd),sframe(:,:,isd),eframe(:,:,i),sdbtu,sdbtv)
              call project_hvec(su2(1,isd),sv2(1,isd),sframe(:,:,isd),eframe(:,:,i),u2,v2)
              swild10(1:3,1:3)=(pframe(:,:,n1)+pframe(:,:,n2))/2
              call project_hvec((tau(1,n1)+tau(1,n2))/2,(tau(2,n1)+tau(2,n2))/2, &
     &swild10(1:3,1:3),eframe(:,:,i),taux2,tauy2)
            endif !ics
            hat_gam_x=sdbtu+dt*(-dprdx/rho0+0.69*grav*detpdx)+ &
     &(1-theta2)*cori(isd)*dt*v2-grav*(1-thetai)*dt*detadx+dt*taux2/htot
            hat_gam_y=sdbtv+dt*(-dprdy/rho0+0.69*grav*detpdy)- &
     &(1-theta2)*cori(isd)*dt*u2-grav*(1-thetai)*dt*detady+dt*tauy2/htot
#ifdef USE_WWM
            !wwave_force in eframe
            hat_gam_x=hat_gam_x+dt*wwave_force(1,1,isd)
            hat_gam_y=hat_gam_y+dt*wwave_force(2,1,isd)
#ifdef DEBUG
            wafo(:,isd,1)=wwave_force(1,1,isd)
            wafo(:,isd,2)=wwave_force(2,1,isd)
#endif
#endif /*USE_WWM*/
            del=hhat(isd)**2+(theta2*cori(isd)*dt*htot)**2 !delta>0
            !Gamma vector (eframe)
            gam_x=(hhat(isd)*htot*hat_gam_x+theta2*cori(isd)*htot*htot*dt*hat_gam_y)/del
            gam_y=(hhat(isd)*htot*hat_gam_y-theta2*cori(isd)*htot*htot*dt*hat_gam_x)/del
            av_elem_x=av_elem_x+htot*gam_x/3
            av_elem_y=av_elem_y+htot*gam_y/3
          enddo !j
          ghat1(1,i)=area(i)*av_elem_x
          ghat1(2,i)=area(i)*av_elem_y
        endif !lm2d

!       Horizontal viscosity
        horx=0 
        hory=0
        do j=1,3 !side
          isd=elside(j,i)
          do k=kbs(isd)+1,nvrt
            horx=horx+area(i)/3*(zs(k,isd)-zs(k-1,isd))*(d2uv(1,k,isd)+d2uv(1,k-1,isd))/2
            hory=hory+area(i)/3*(zs(k,isd)-zs(k-1,isd))*(d2uv(2,k,isd)+d2uv(2,k-1,isd))/2
          enddo !k
          horx=horx-dt*chigamma*area(i)/3*d2uv(1,kbs(isd)+1,isd)
          hory=hory-dt*chigamma*area(i)/3*d2uv(2,kbs(isd)+1,isd)
        enddo !j=1,3
        ghat1(1,i)=ghat1(1,i)+dt*horx
        ghat1(2,i)=ghat1(2,i)+dt*hory

!       Radiation stress (3D only; 2D part has been done above)
#ifdef USE_WWM
          if(.not.lm2d) then !3D
            rs1=0 
            rs2=0
            do j=1,3 !side
              isd=elside(j,i)
              do k=kbs(isd)+1,nvrt
                !wwave_force in eframe
                rs1=rs1+area(i)/3*(zs(k,isd)-zs(k-1,isd))* & !(swild(1)+swild(3))/2
     &(wwave_force(1,k,isd)+wwave_force(1,k-1,isd))/2
                rs2=rs2+area(i)/3*(zs(k,isd)-zs(k-1,isd))* & !(swild(2)+swild(4))/2
     &(wwave_force(2,k,isd)+wwave_force(2,k-1,isd))/2
              enddo !k
              rs1=rs1-dt*chigamma*area(i)/3*wwave_force(1,kbs(isd)+1,isd)
              rs2=rs2-dt*chigamma*area(i)/3*wwave_force(2,kbs(isd)+1,isd)
            enddo !j=1,3
            ghat1(1,i)=ghat1(1,i)+dt*rs1
            ghat1(2,i)=ghat1(2,i)+dt*rs2
          endif !lm2d
#endif /*USE_WWM*/

!       Baroclinic force
        if(ibc==0) then
          if(prho(1,elnode(1,i))<-98.or.prho(1,elnode(2,i))<-98.or.prho(1,elnode(3,i))<-98) then
            write(errmsg,*)'Impossible dry 5'
            call parallel_abort(errmsg)
          endif

          if(iupwind_t==0) then
!           ELM option
            swild(1:2)=0 !averaged bottom bcc (in eframe if ics=2)
            do j=1,3 !side
              isd=elside(j,i)
              !Rotate bcc to eframe if ics=2
              do k=kbs(isd),nvrt
                if(ics==1) then
                  swild10(k,1:2)=bcc(1:2,k,isd)
                else !lat/lon
                  call project_hvec(bcc(1,k,isd),bcc(2,k,isd),sframe(:,:,isd),eframe(:,:,i),swild10(k,1),swild10(k,2))
                endif !ics
              enddo !k
            
              bigfc1=0; bigfc2=0 !\int_{-h}^\eta f_c dz; (f_c=bcc); eframe if ics=2
              do k=kbs(isd)+1,nvrt
                bigfc1=bigfc1+(zs(k,isd)-zs(k-1,isd))*(swild10(k,1)+swild10(k-1,1))/2 !(bcc(1,k,isd)+bcc(1,k-1,isd))/2
                bigfc2=bigfc2+(zs(k,isd)-zs(k-1,isd))*(swild10(k,2)+swild10(k-1,2))/2 !(bcc(2,k,isd)+bcc(2,k-1,isd))/2
              enddo !k
              ghat1(1,i)=ghat1(1,i)+dt*area(i)*bigfc1/3
              ghat1(2,i)=ghat1(2,i)+dt*area(i)*bigfc2/3
              swild(1:2)=swild(1:2)+swild10(kbs(isd)+1,1:2)/3 !bcc(1:2,kbs(isd)+1,isd)/3
            enddo !side
            ghat1(1:2,i)=ghat1(1:2,i)-area(i)*chigamma*dt*dt*swild(1:2)

          else
!           Upwind/TVD
!           swild2(k,:) = \sum_{l=k}^N dr*dz; whole level (and eframe if ics=2)
            swild2(nvrt,1:2)=0
            do k=nvrt-1,kbe(i),-1
              !dr_dxy in eframe
              swild2(k,1:2)=swild2(k+1,1:2)+dr_dxy(1:2,k+1,i)*(ze(k+1,i)-ze(k,i))
            enddo !k

            do k=kbe(i)+1,nvrt
!             Ramp-up factor
              ghat1(1:2,i)=ghat1(1:2,i)-rampbc*dt*area(i)*grav/rho0*(ze(k,i)-ze(k-1,i))/2* &
     &(2*swild2(k,1:2)+dr_dxy(1:2,k,i)*(ze(k,i)-ze(k-1,i)))
            enddo !k
            ghat1(1:2,i)=ghat1(1:2,i)+rampbc*chigamma*dt*dt*area(i)*grav/rho0*swild2(kbe(i)+1,1:2)
          endif !iupwind_t
        endif !ibc==0

!       Non-hydrostatic press. gradient
        if(nonhydro==1) then
          swild(1:2)=0 !averaged bottom value
          do j=1,3 !side
            isd=elside(j,i)
            bigfc1=0; bigfc2=0 !\int_{-h}^\eta f dz; (f = -\nabla qnon)
            do k=kbs(isd)+1,nvrt
              bigfc1=bigfc1-(zs(k,isd)-zs(k-1,isd))*(dqnon_dxy(1,k,isd)+dqnon_dxy(1,k-1,isd))/2
              bigfc2=bigfc2-(zs(k,isd)-zs(k-1,isd))*(dqnon_dxy(2,k,isd)+dqnon_dxy(2,k-1,isd))/2
            enddo !k
            ghat1(1,i)=ghat1(1,i)+dt*area(i)*bigfc1/3
            ghat1(2,i)=ghat1(2,i)+dt*area(i)*bigfc2/3
            swild(1:2)=swild(1:2)-dqnon_dxy(1:2,kbs(isd)+1,isd)/3
          enddo !side
          ghat1(1:2,i)=ghat1(1:2,i)-area(i)*chigamma*dt*dt*swild(1:2)
        endif !nonhydro==1


!       Debug
!        if(myrank==irank0) write(96,*)i,ielg(i),ghat1(1:2,i)     

      enddo !i=1,nea

!      if(myrank==irank0) write(96,*)'=================================='

      if(myrank==0) write(16,*)'done 2nd preparation'

#ifdef INCLUDE_TIMING
! end preparations
      wtmp2=mpi_wtime()
      wtimer(6,1)=wtimer(6,1)+wtmp2-wtmp1
! start solver
      wtmp1=wtmp2
#endif

!...  setup coefficient matrix, sparsem, for the wave equation
!...  No elevation essential b.c. are imposed yet but other b.c. is imposed
      do i=1,np !resident only
        do j=0,nnp(i)
          sparsem(j,i)=0
        enddo !j
        qel(i)=0

!	Area integrals I_{1,4,7}
        do j=1,nne(i)
          ie=indel(j,i)
          id=iself(j,i)

          if(ihydraulics/=0.and.nhtblocks>0) then
            if(isblock_el(ie)>0) cycle !active block
          endif

!	  I_1
          !n2=elnode(nx(id,1),ie)
          !n3=elnode(nx(id,2),ie)
          !dot1=(xnd(n3)-xnd(n2))**2+(ynd(n3)-ynd(n2))**2
          !dot2=(xnd(n3)-xnd(n2))*(xnd(i)-xnd(n3))+(ynd(n3)-ynd(n2))*(ynd(i)-ynd(n3))
          id2=nx(id,1)
          id3=nx(id,2)
          dot1=(xel(id3,ie)-xel(id2,ie))**2+(yel(id3,ie)-yel(id2,ie))**2
          dot2=(xel(id3,ie)-xel(id2,ie))*(xel(id,ie)-xel(id3,ie))+ &
     &         (yel(id3,ie)-yel(id2,ie))*(yel(id,ie)-yel(id3,ie))
          dot3=-dot1-dot2
          if(lm2d) then
            hhatb=0
            avg2=0 !average for Coriolis part
            do jj=1,3 !side
              isd=elside(jj,ie)
              if(idry_s(isd)==0) then
                etam=(eta2(isidenode(1,isd))+eta2(isidenode(2,isd)))/2
                htot=etam+dps(isd)
                if(hhat(isd)<=0) then
                  write(errmsg,*)'Impossible dry 52:',hhat(isd),iplg(isidenode(1:2,isd))
                  call parallel_abort(errmsg)
                endif
                del=hhat(isd)*hhat(isd)+(theta2*cori(isd)*dt*htot)**2 !>0
                hhatb=hhatb+htot*htot*hhat(isd)/del/3
                avg2=avg2+cori(isd)*(htot*dt)**3/del/3
              endif !idry_s
            enddo !jj
          else !3D
            hhatb=(hhat(elside(1,ie))+hhat(elside(2,ie))+hhat(elside(3,ie)))/3
          endif 
          tmp0=area(ie)/6+grav*thetai**2*dt**2/4/area(ie)*hhatb*dot1
          tmpj=area(ie)/12+grav*thetai**2*dt**2/4/area(ie)*hhatb*dot2
          tmpj1=area(ie)/12+grav*thetai**2*dt**2/4/area(ie)*hhatb*dot3
          sparsem(0,i)=sparsem(0,i)+tmp0
          sparsem(j,i)=sparsem(j,i)+tmpj
          if(isbnd(1,i)==0.and.j==nne(i)) then
            indx=1
          else
            indx=j+1
          endif
          sparsem(indx,i)=sparsem(indx,i)+tmpj1

          if(lm2d) then !additional Coriolis part
            tmp=grav*theta2*thetai**2*avg2/2
            sparsem(j,i)=sparsem(j,i)+tmp
            sparsem(indx,i)=sparsem(indx,i)-tmp
          endif !lm2d

!	  Check dominance
          if(hhatb<0) then
            if(ihhat==0.and.ifort12(1)==0) then
              ifort12(1)=1
              write(12,*)'Modified depth < 0:',it,iplg(i),j,hhatb
            endif
            if(ihhat==1) then
              write(errmsg,*)'Impossible hhat:',hhatb
              call parallel_abort(errmsg)
            endif
          endif

!	  I_4
          isd1=elside(1,ie)
          isd2=elside(2,ie)
          isd3=elside(3,ie)
          if(ics==1) then
            dot1=dldxy(id,1,ie)*(bigu(1,isd1)+bigu(1,isd2)+bigu(1,isd3))/3+ &
     &dldxy(id,2,ie)*(bigu(2,isd1)+bigu(2,isd2)+bigu(2,isd3))/3
          else
            call project_hvec(bigu(1,isd1),bigu(2,isd1),sframe(:,:,isd1),eframe(:,:,ie),bigu1,bigv1)
            call project_hvec(bigu(1,isd2),bigu(2,isd2),sframe(:,:,isd2),eframe(:,:,ie),bigu2,bigv2)
            call project_hvec(bigu(1,isd3),bigu(2,isd3),sframe(:,:,isd3),eframe(:,:,ie),bigu3,bigv3)
            dot1=dldxy(id,1,ie)*(bigu1+bigu2+bigu3)/3+dldxy(id,2,ie)*(bigv1+bigv2+bigv3)/3
          endif !ics
          dot2=dldxy(id,1,ie)*ghat1(1,ie)+dldxy(id,2,ie)*ghat1(2,ie)
        
          qel(i)=qel(i)+(1-thetai)*dt*area(ie)*dot1+thetai*dt*dot2
#ifdef USE_WWM
          if(RADFLAG.eq.'VOR'.and.idry_e(ie)==0) then
            sum1=0; sum2=0 !in eframe
            do m=1,3 !wet sides
              isd=elside(m,ie)
              do k=kbs(isd),nvrt-1
                sum1=sum1+(zs(k+1,isd)-zs(k,isd))*(stokes_vel_sd(1,k+1,isd)+stokes_vel_sd(1,k,isd))/2/3
                sum2=sum2+(zs(k+1,isd)-zs(k,isd))*(stokes_vel_sd(2,k+1,isd)+stokes_vel_sd(2,k,isd))/2/3
              enddo !k
            enddo !m
            dot3=dldxy(id,1,ie)*sum1+dldxy(id,2,ie)*sum2
            qel(i)=qel(i)+dt*dot3
          endif
#endif

          do l=1,3
            if(id==l) then
              fac=2
            else
              fac=1
            endif
            nd=elnode(l,ie)
            if(imm==2) then
!Error: vel. not projected for ics=2
              call update_bdef(time,xctr(ie),yctr(ie),dep,swild)
              ubed=swild(1); vbed=swild(2); wbed=swild(3)
              dpdx=0; dpdy=0
              do m=1,3
                dpdx=dpdx+dp(elnode(m,ie))*dldxy(m,1,ie)
                dpdy=dpdy+dp(elnode(m,ie))*dldxy(m,2,ie)
              enddo !m   
              vnorm=ubed*dpdx+vbed*dpdy+wbed
              qel(i)=qel(i)+area(ie)/12*fac*(eta2(nd)+dt*vnorm)
            else
              qel(i)=qel(i)+area(ie)/12*fac*(eta2(nd)+bdef2(nd)-bdef1(nd))
            endif !imm
          enddo !l

!...      I_7: Impose Point Source volume
!Error: need to reconcile with ICM
!Error: did not add to 3D continuity eq.
          qel(i)=qel(i)+dt/3*vsource(ie)           

!...      Impose Point Source volume at the surface layer (including 3D continuity); added by YC
#ifdef USE_ICM
          if(iWQPS==2) then
            do l=1,3
              if(id==l) then
                fac=2
              else
                fac=1
              endif
              nd=elnode(l,ie)
!YC            raintmp=raintmp+area(ie)/12*fac*fluxprc(nd)/rho0*dt
!YC            qel(i)=qel(i)+area(ie)/12*fac*fluxprc(nd)/rho0*dt
              qel(i)=qel(i)-PSQ(ie)*dt/12*fac
            enddo !l
          endif
#endif /*USE_ICM*/
        enddo !j=1,nne(i)

!	bnd integrals I_{2,3,5,6}; they all vanish at land bnds 
!	I_2,6 are not needed if essential b.c. are enforced by elminating rows and columns
        if(isbnd(1,i)>0) then !open bnd node
!          ibnd=isbnd(1,i)
          do l=1,2 !two open bnd sides
            if(l==1) then
              ie=indel(1,i)
              id=iself(1,i)
              isd=elside(nx(id,2),ie)
              nj=elnode(nx(id,1),ie)
              ind=1
            else
              ie=indel(nne(i),i)
              id=iself(nne(i),i)
              isd=elside(nx(id,1),ie)
              nj=elnode(nx(id,2),ie)
              ind=nnp(i)
            endif

            nd=isidenode(1,isd)+isidenode(2,isd)-i
            if(nd/=nj) then
              write(errmsg,*)'Impossible 79'
              call parallel_abort(errmsg)
            endif

!	    I_3 
            if(isbs(isd)>0.and.ifltype(max(1,isbs(isd)))/=0) then !.and.(.not.lelbc(i))) then 
!             Natural or Flather 1 b.c.
!             Calculate I_3 even if i is on essential b.c. so as to check symmetry later
!             especially for Flather b.c.
              if(idry_s(isd)==1) then
                write(errmsg,*)'Dry flow bnd:',islg(isd),iplg(i),iplg(nd)
                call parallel_abort(errmsg)
              endif

              bigvn=0
              do k=kbs(isd),nvrt-1
                !uth, vth in lat/lon frame if ics=2
                if(ics==1) then
                  vn1=uth(k,isd)*sframe(1,1,isd)+vth(k,isd)*sframe(2,1,isd)
                  vn2=uth(k+1,isd)*sframe(1,1,isd)+vth(k+1,isd)*sframe(2,1,isd)
                else 
                  call project_hvec(uth(k,isd),vth(k,isd),pframe(:,:,i),sframe(:,:,isd),vn1,vtmp)
                  call project_hvec(uth(k+1,isd),vth(k+1,isd),pframe(:,:,i),sframe(:,:,isd),vn2,vtmp)
                endif !ics
                bigvn=bigvn+(zs(k+1,isd)-zs(k,isd))*(vn1+vn2)/2
              enddo !k
              ri3=distj(isd)*bigvn/2
              if(ifltype(isbs(isd))==-1) then !Flather 1
                if(eta_mean(i)<-98.or.eta_mean(nj)<-98) then
                  write(errmsg,*)'Mismatch 1'
                  call parallel_abort(errmsg)
                endif
                if(dps(isd)<=0) then
                  write(errmsg,*)'Negative depth at Flather bnd:',i,dps(isd)
                  call parallel_abort(errmsg)
                endif
                con0=distj(isd)/6*sqrt(grav*dps(isd)) !for coefficient matrix
                ri3=ri3-con0*(2*eta_mean(i)+eta_mean(nj))
                sparsem(0,i)=sparsem(0,i)+thetai*dt*con0*2
                sparsem(ind,i)=sparsem(ind,i)+thetai*dt*con0
              endif !Flather 1
              qel(i)=qel(i)-thetai*dt*ri3
            endif

!	    I_5
            if(isbs(isd)>0.and.idry_s(isd)==0) then
              if(ics==1) then
                Unbar=bigu(1,isd)*sframe(1,1,isd)+bigu(2,isd)*sframe(2,1,isd)
              else
                Unbar=bigu(1,isd)
              endif !ics
              tmp0=(1-thetai)*dt*distj(isd)*Unbar/2
              !Overwrite tmp0 for vortex formulation
#ifdef USE_WWM
              if(RADFLAG.eq.'VOR') then
                sum1=0 !integral; x-comp.
                sum2=0 !integral
                do k=kbs(isd),nvrt-1 !isd is wet
                  sum1=sum1+(zs(k+1,isd)-zs(k,isd))*(stokes_vel_sd(1,k+1,isd)+stokes_vel_sd(1,k,isd))/2
                  sum2=sum2+(zs(k+1,isd)-zs(k,isd))*(stokes_vel_sd(2,k+1,isd)+stokes_vel_sd(2,k,isd))/2
                enddo !k
                if(ics==1) then
                  Unbar=sum1*sframe(1,1,isd)+sum2*sframe(2,1,isd)
                else
                  call project_hvec(sum1,sum2,pframe(:,:,i),sframe(:,:,isd),Unbar,tmp)
                  !Unbar=bigu(1,isd)
                endif !ics
                tmp0=thetai*dt*distj(isd)*Unbar/2
              endif !RADFLAG
#endif/*USE_WWM*/

              qel(i)=qel(i)-tmp0 !(1-thetai)*dt*distj(isd)*Unbar/2
            endif
          enddo !l=1,2 sides
        endif !isbnd: bnd node i

        !Hydraulic blocks for I_3 and I_5
        if(ihydraulics/=0.and.nhtblocks>0) then; if(isblock_nd(1,i)>0) then
          do j=1,nne(i)
            ie=indel(j,i)
            id=iself(j,i)
            if(isbnd(1,i)/=0.and.j==1) then !bnd node
              lim=2 !1 extra side
            else
              lim=1
            endif
            do m=1,lim
              if(m==1) then
                isd=elside(nx(id,1),ie)
              else !bnd node; 1 extra side
                isd=elside(nx(id,2),ie)
              endif
              if(i/=isidenode(1,isd).and.i/=isidenode(2,isd)) call parallel_abort('MAIN: impossible 51')
              if(isblock_sd(1,isd)>0.and.isblock_sd(2,isd)>0) then !active block face side
                jblock=isblock_sd(1,isd)
                jface=isblock_sd(2,isd)
                !Compute if the local side normal is in/against block dir
                dot1=dot_product(dir_block(1:3,jblock),sframe(1:3,1,isd))
                ss=sign(1.d0,dot1)
                if(jface==1) then
                  !Out-normal for I_3,5 is along block dir
                else
                  !Out-normal for I_3,5 is against block dir
                  ss=-ss
                endif !jface

                !I_5
                if(ics==1) then
                  Unbar=bigu(1,isd)*sframe(1,1,isd)+bigu(2,isd)*sframe(2,1,isd)
                else
                  Unbar=bigu(1,isd)
                endif !ics
                Unbar=Unbar*ss
                qel(i)=qel(i)-(1-thetai)*dt*distj(isd)*Unbar/2

                !I_3
                if(idry_s(isd)==0) then
                  htot=zs(nvrt,isd)-zs(kbs(isd),isd)
                  if(htot<h0) call parallel_abort('MAIN: hydrau. dep<h0')
                  ri3=(1-block_nudge)*Unbar/2*distj(isd)+ &
     &block_nudge*vnth_block(jface,jblock)*(3-2*jface)*distj(isd)*htot/2 !sign added
                  qel(i)=qel(i)-thetai*dt*ri3
                endif !wet side

                !Check
                !write(12,*)'I_3,5:',iplg(i),j,jblock,jface,ri3,Unbar,vnth_block(jface,jblock)
              endif !isblock_sd
            enddo !m
          enddo !j=1,nne(i)
        endif; endif !ihydraulics

      enddo !i=1,np

!     Check symmetry 
!      fdb='spars_0000'
!      lfdb=len_trim(fdb)
!      write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
!      open(10,file=fdb,status='unknown')

!     2D has implicit Corioilis which destroys symmetry
!temp fix
!      if(1==2.and..not.lm2d) then
      if(.not.lm2d) then
        do i=1,np
          do j=1,nnp(i)
            nd=indnd(j,i)
            if(nd<=np) then !nd resident
              in1=0
              do l=1,nnp(nd)
                if(indnd(l,nd)==i) in1=l
              enddo !l
              if(in1==0) then
                write(errmsg,*)'Not resident:',iplg(i),iplg(nd)
                call parallel_abort(errmsg)
              endif
              if(abs(sparsem(j,i)-sparsem(in1,nd))>1.e-5) then
                write(errmsg,*)'Matrix not symmetric:',iplg(i),j,iplg(nd),sparsem(j,i),sparsem(in1,nd)
                call parallel_abort(errmsg)
              endif
              irank_s=myrank
            else !nd is ghost
              if(.not.associated(ipgl(iplg(nd))%next)) call parallel_abort('Wrong ghost')
              irank_s=ipgl(iplg(nd))%next%rank
            endif

!           Output
!           write(10,*)iplg(i),iplg(nd),irank_s,sparsem(j,i),sparsem(0,i)
          enddo !j
        enddo !i=1,np
      endif !.not.lm2d
!      close(10)

!     Save eta1
      eta1=eta2

!     Solve the wave equation for elevations at each element center in aug. subdomain
      call solve_jcg(it,moitn0,mxitn0,rtol0,sparsem,eta2,qel,elbc,lelbc)

!     Exchange eta2 to ensure consistency across processors
#ifdef INCLUDE_TIMING
      cwtmp=mpi_wtime()
#endif
      call exchange_p2d(eta2)
#ifdef INCLUDE_TIMING
      wtimer(7,2)=wtimer(7,2)+mpi_wtime()-cwtmp
#endif

      etatotl=0
      do i=1,np
        if(eta2(i)>elevmax(i)) elevmax(i)=eta2(i) !only for residents

        if(associated(ipgl(iplg(i))%next)) then !interface node
          if(ipgl(iplg(i))%next%rank<myrank) cycle !already in the sum so skip
        endif
        etatotl=etatotl+abs(eta2(i))
      enddo !i
#ifdef INCLUDE_TIMING
      cwtmp=mpi_wtime()
#endif
      call mpi_allreduce(etatotl,etatot,1,rtype,MPI_SUM,comm,ierr)
#ifdef INCLUDE_TIMING
      wtimer(7,2)=wtimer(7,2)+mpi_wtime()-cwtmp
#endif

      if(myrank==0) write(16,*)'done solver; ', 'etatot=',etatot

#ifdef INCLUDE_TIMING
!  end solver
      wtmp2=mpi_wtime()
      wtimer(7,1)=wtimer(7,1)+wtmp2-wtmp1
!  start momentum
      wtmp1=wtmp2
#endif

!
!************************************************************************
!									*
!		Momentum equations					*
!									*
!************************************************************************
!

!     Precompute elevation gradient, atmo. pressure and earth tidal potential
!     (in sframe if ics=2).
!     Initialize for dry sides and exchange
      deta2_dx=0; deta2_dy=0; deta1_dx=0; deta1_dy=0; dpr_dx=0; dpr_dy=0; detp_dx=0; detp_dy=0
      do j=1,ns !resident
        if(idry_s(j)==1) cycle

!       Wet side
        icount1=0 !for deta1 & dpr
        icount2=0 !for deta2
        icount3=0 !for detp
        do l=1,2 !elements
          ie=isdel(l,j)
          if(ie/=0) then
            itmp=0
            do m=1,3
             nd=elnode(m,ie)
             if(eta2(nd)+dp(nd)<=h0) itmp=1
            enddo !m
            if(itmp==0) then !wet
              icount2=icount2+1
              do m=1,3
                tmpx=eta2(elnode(m,ie))*dldxy(m,1,ie) !!eframe if ics=2
                tmpy=eta2(elnode(m,ie))*dldxy(m,2,ie)
                if(ics==2) then
                  call project_hvec(tmpx,tmpy,eframe(:,:,ie),sframe(:,:,j),tmpxs,tmpys)
                  tmpx=tmpxs
                  tmpy=tmpys
                endif !ics
                deta2_dx(j)=deta2_dx(j)+tmpx !sframe if ics=2
                deta2_dy(j)=deta2_dy(j)+tmpy
              enddo !m
            endif !wet at n+1
            if(idry_e(ie)==0) then
              icount1=icount1+1
              if(dpe(ie)>=tip_dp) icount3=icount3+1
              do m=1,3
                nd=elnode(m,ie)
                tmpx1=eta1(nd)*dldxy(m,1,ie) !eframe if ics=2
                tmpy1=eta1(nd)*dldxy(m,2,ie)
                tmpx2=pr(nd)*dldxy(m,1,ie)
                tmpy2=pr(nd)*dldxy(m,2,ie)
                tmpx3=etp(nd)*dldxy(m,1,ie)
                tmpy3=etp(nd)*dldxy(m,2,ie)
                if(ics==2) then
                  call project_hvec(tmpx1,tmpy1,eframe(:,:,ie),sframe(:,:,j),tmpx1s,tmpy1s)
                  call project_hvec(tmpx2,tmpy2,eframe(:,:,ie),sframe(:,:,j),tmpx2s,tmpy2s)
                  call project_hvec(tmpx3,tmpy3,eframe(:,:,ie),sframe(:,:,j),tmpx3s,tmpy3s)
                  tmpx1=tmpx1s; tmpy1=tmpy1s
                  tmpx2=tmpx2s; tmpy2=tmpy2s
                  tmpx3=tmpx3s; tmpy3=tmpy3s
                endif !ics
            
                deta1_dx(j)=deta1_dx(j)+tmpx1 !sframe if ics=2
                deta1_dy(j)=deta1_dy(j)+tmpy1
                dpr_dx(j)=dpr_dx(j)+tmpx2
                dpr_dy(j)=dpr_dy(j)+tmpy2
                if(dpe(ie)>=tip_dp) then
                  detp_dx(j)=detp_dx(j)+tmpx3
                  detp_dy(j)=detp_dy(j)+tmpy3
                endif
              enddo !m
            endif
          endif !ie/=0
        enddo !l=1,2
        if(icount1/=0) then
          deta1_dx(j)=deta1_dx(j)/icount1
          deta1_dy(j)=deta1_dy(j)/icount1
          dpr_dx(j)=dpr_dx(j)/icount1
          dpr_dy(j)=dpr_dy(j)/icount1
        endif
        if(icount3/=0) then
          detp_dx(j)=detp_dx(j)/icount3
          detp_dy(j)=detp_dy(j)/icount3
        endif
        if(icount2/=0) then
          deta2_dx(j)=deta2_dx(j)/icount2
          deta2_dy(j)=deta2_dy(j)/icount2
        endif
      enddo !j=1,ns

!     Compute bottom index for sides for zeroing out fluxes for Z layers
!     This is no longer done in the newer version and kbs_e is not
!     used now
      kbs_e=0 !larger of the 2 element bottom indices
      do j=1,ns
        if(idry_s(j)==1) cycle

!       Wet
        do l=1,2 !element
          ie=isdel(l,j)
          if(ie/=0.and.idry_e(max(1,ie))==0.and.kbe(max(1,ie))>kbs_e(j)) kbs_e(j)=kbe(ie)
        enddo !l
        if(kbs_e(j)==0) then
          write(errmsg,*)'Cannot find the higher bottom:',islg(j),(ielg(isdel(l,j)),l=1,2)
          call parallel_abort(errmsg)
        endif
!        if(kbs(j)>kbs_e(j)) then
!          write(errmsg,*)'Side index > elemnt:',kbs(j),kbs_e(j)
!          call parallel_abort(errmsg)
!        endif
        if(lm2d.and.kbs_e(j)/=1) then
          write(errmsg,*)'2D bottom index wrong:',kbs(j),kbs_e(j),iplg(isidenode(1:2,j))
          call parallel_abort(errmsg)
        endif
      enddo !j=1,ns

      allocate(swild99(9,nsa),stat=istat)
      if(istat/=0) call parallel_abort('MAIN: fail to allocate swild99')
!'
      swild99(1,:)=deta1_dx(:); swild99(2,:)=deta1_dy(:); swild99(3,:)=deta2_dx(:) 
      swild99(4,:)=deta2_dy(:); swild99(5,:)=dpr_dx(:); swild99(6,:)=dpr_dy(:)
      swild99(7,:)=detp_dx(:); swild99(8,:)=detp_dy(:); swild99(9,:)=kbs_e(:)
#ifdef INCLUDE_TIMING
      cwtmp=mpi_wtime()
#endif
      call exchange_s2d_9(swild99)
#ifdef INCLUDE_TIMING
      wtimer(8,2)=wtimer(8,2)+mpi_wtime()-cwtmp
#endif
      deta1_dx(:)=swild99(1,:); deta1_dy(:)=swild99(2,:); deta2_dx(:)=swild99(3,:)
      deta2_dy(:)=swild99(4,:); dpr_dx(:)=swild99(5,:); dpr_dy(:)=swild99(6,:)
      detp_dx(:)=swild99(7,:); detp_dy(:)=swild99(8,:); kbs_e(:)=nint(swild99(9,:))
      deallocate(swild99)

!     Save vel. at previous step (for hydraulics etc)
      allocate(swild98(2,nvrt,nsa))
      swild98(1,:,:)=su2(:,:)
      swild98(2,:,:)=sv2(:,:)

!...  Along each side
!     su2, sv2 in sframe if ics=2
#ifdef DEBUG
      bpgr = 0.d0
#endif

      do j=1,nsa !augumented
        if(idry_s(j)==1) then
          do k=1,nvrt
            su2(k,j)=0
            sv2(k,j)=0
          enddo !k
          cycle
        endif

!	Wet sides
        node1=isidenode(1,j)
        node2=isidenode(2,j)
!       ll frame at side
        swild10(1:3,1:3)=(pframe(:,:,node1)+pframe(:,:,node2))/2

        if(lm2d) then !2D
!-------------------------------------------------------------------------------------
          !Warning: don't use eta2 which is updated
          htot=zs(nvrt,j)-zs(kbs(j),j) !(eta2(node1)+eta2(node2))/2+dps(j)
          if(hhat(j)<=0) then
            write(errmsg,*)'Impossible dry 55:',hhat(j),iplg(isidenode(1:2,j))
            call parallel_abort(errmsg)
          endif
          del=hhat(j)*hhat(j)+(theta2*cori(j)*dt*htot)**2 !delta > 0
          !\hat{Gamma}
          !rotate wind to sframe
          taux2=(tau(1,node1)+tau(1,node2))/2
          tauy2=(tau(2,node1)+tau(2,node2))/2
          if(ics==2) then
            call project_hvec(taux2,tauy2,swild10(1:3,1:3),sframe(:,:,j),taux2s,tauy2s)
            taux2=taux2s
            tauy2=tauy2s
          endif !ics
          hat_gam_x=sdbt(1,1,j)+dt*(-dpr_dx(j)/rho0+0.69*grav*detp_dx(j))+ &
     &(1-theta2)*cori(j)*dt*sv2(1,j)-grav*(1-thetai)*dt*deta1_dx(j)+dt*taux2/htot
          hat_gam_y=sdbt(2,1,j)+dt*(-dpr_dy(j)/rho0+0.69*grav*detp_dy(j))- &
     &(1-theta2)*cori(j)*dt*su2(1,j)-grav*(1-thetai)*dt*deta1_dy(j)+dt*tauy2/htot
!         Radiation stress
#ifdef  USE_WWM
            !wwave_force in eframe
            if(ics==1) then
              tmp1=wwave_force(1,1,j)
              tmp2=wwave_force(2,1,j)
            else !use swild10 as approx.
              call project_hvec(wwave_force(1,1,j),wwave_force(2,1,j),swild10(1:3,1:3),sframe(:,:,j),tmp1,tmp2)
            endif !ics
            hat_gam_x=hat_gam_x+dt*tmp1
            hat_gam_y=hat_gam_y+dt*tmp2
#endif /*USE_WWM*/
          !Add deta2 to \hat{Gamma}
          hat_gam_x=hat_gam_x-grav*thetai*dt*deta2_dx(j)
          hat_gam_y=hat_gam_y-grav*thetai*dt*deta2_dy(j)
          su2(1,j)=(hhat(j)*htot*hat_gam_x+theta2*cori(j)*dt*htot*htot*hat_gam_y)/del
          sv2(1,j)=(hhat(j)*htot*hat_gam_y-theta2*cori(j)*dt*htot*htot*hat_gam_x)/del
          su2(1,j)=max(-rmaxvel,min(rmaxvel,su2(1,j)))
          sv2(1,j)=max(-rmaxvel,min(rmaxvel,sv2(1,j)))
#ifdef DEBUG
          bpgr(j,1) = -grav*(1-thetai)*deta1_dx(j)-grav*thetai*deta2_dx(j)
          bpgr(j,2) = -grav*(1-thetai)*deta1_dy(j)-grav*thetai*deta2_dy(j)
#endif

!-------------------------------------------------------------------------------------
        else !3D
!-------------------------------------------------------------------------------------
!       Define layer thickness & diffusivities
        do k=kbs(j)+1,nvrt
          dzz(k)=zs(k,j)-zs(k-1,j)
          if(dzz(k)<=0) call parallel_abort('MAIN: dzz=0 in momentum')
          dfz(k)=(ptbt(3,k,node1)+ptbt(3,k,node2)+ptbt(3,k-1,node1)+ptbt(3,k-1,node2))/4
        enddo !k

!       Define bottom level
        if(ibottom_bc==1) then !old @kbs(j)+1
          ibtm=kbs(j)+1
        else !new @ kbs
          ibtm=kbs(j)
        endif !ibottom_bc

!	Coefficient matrix 
        !ndim=nvrt-kbs(j)
        ndim=nvrt-ibtm+1
        do k=ibtm,nvrt !kbs(j)+1,nvrt
          kin=k-ibtm+1   !k-kbs(j) !eq. #
          alow(kin)=0 
          cupp(kin)=0
          bdia(kin)=0
          if(k<nvrt) then
            tmp=dt*dfz(k+1)/dzz(k+1)
            cupp(kin)=cupp(kin)+dzz(k+1)/6-tmp
            bdia(kin)=bdia(kin)+dzz(k+1)/3+tmp
          endif

          !if(k>kbs(j)+1) then
          if(k>ibtm) then
            tmp=dt*dfz(k)/dzz(k)
            alow(kin)=alow(kin)+dzz(k)/6-tmp
            bdia(kin)=bdia(kin)+dzz(k)/3+tmp
          else !b.c.
            bdia(kin)=bdia(kin)+dt*chi(j)
          endif
        enddo !k

!	RHS 
!	b.c. to be imposed at the end
        do k=ibtm,nvrt !kbs(j)+1,nvrt
          kin=k-ibtm+1  !k-kbs(j)
          rrhs(1,kin)=0
          rrhs(2,kin)=0
!	  Elevation gradient, atmo. pressure and tidal potential
          if(k<nvrt) then
            rrhs(1,kin)=rrhs(1,kin)-dzz(k+1)/2*dt*(grav*thetai*deta2_dx(j)+ &
                        grav*(1-thetai)*deta1_dx(j)+dpr_dx(j)/rho0-0.69*grav*detp_dx(j))
            rrhs(2,kin)=rrhs(2,kin)-dzz(k+1)/2*dt*(grav*thetai*deta2_dy(j)+ &
                        grav*(1-thetai)*deta1_dy(j)+dpr_dy(j)/rho0-0.69*grav*detp_dy(j))
          endif
          !if(k>kbs(j)+1) then 
          if(k>ibtm) then 
            rrhs(1,kin)=rrhs(1,kin)-dzz(k)/2*dt*(grav*thetai*deta2_dx(j)+ &
                        grav*(1-thetai)*deta1_dx(j)+dpr_dx(j)/rho0-0.69*grav*detp_dx(j))
            rrhs(2,kin)=rrhs(2,kin)-dzz(k)/2*dt*(grav*thetai*deta2_dy(j)+ &
                        grav*(1-thetai)*deta1_dy(j)+dpr_dy(j)/rho0-0.69*grav*detp_dy(j))
          endif

!	  Coriolis, advection, wind stress, and horizontal viscosity
          if(k<nvrt) then
            rrhs(1,kin)=rrhs(1,kin)+dzz(k+1)/6*(2*sdbt(1,k,j)+sdbt(1,k+1,j)+ &
     &dt*cori(j)*(2*sv2(k,j)+sv2(k+1,j))+dt*(2*d2uv(1,k,j)+d2uv(1,k+1,j)))
            rrhs(2,kin)=rrhs(2,kin)+dzz(k+1)/6*(2*sdbt(2,k,j)+sdbt(2,k+1,j)- &
     &dt*cori(j)*(2*su2(k,j)+su2(k+1,j))+dt*(2*d2uv(2,k,j)+d2uv(2,k+1,j)))
          else !k=nvrt
            taux2=(tau(1,node1)+tau(1,node2))/2
            tauy2=(tau(2,node1)+tau(2,node2))/2
            if(ics==2) then
              call project_hvec(taux2,tauy2,swild10(1:3,1:3),sframe(:,:,j),taux2s,tauy2s)
              taux2=taux2s
              tauy2=tauy2s
            endif !ics
            rrhs(1,kin)=rrhs(1,kin)+dt*taux2
            rrhs(2,kin)=rrhs(2,kin)+dt*tauy2
          endif

          !if(k>kbs(j)+1) then
          if(k>ibtm) then
            rrhs(1,kin)=rrhs(1,kin)+dzz(k)/6*(2*sdbt(1,k,j)+sdbt(1,k-1,j)+ &
     &dt*cori(j)*(2*sv2(k,j)+sv2(k-1,j))+dt*(2*d2uv(1,k,j)+d2uv(1,k-1,j)))
            rrhs(2,kin)=rrhs(2,kin)+dzz(k)/6*(2*sdbt(2,k,j)+sdbt(2,k-1,j)- &
     &dt*cori(j)*(2*su2(k,j)+su2(k-1,j))+dt*(2*d2uv(2,k,j)+d2uv(2,k-1,j)))
          endif 

!         Baroclinic
          if(ibc==0) then
            if(k<nvrt) then
               rrhs(1,kin)=rrhs(1,kin)+dzz(k+1)/6*dt*(2*bcc(1,k,j)+bcc(1,k+1,j))
               rrhs(2,kin)=rrhs(2,kin)+dzz(k+1)/6*dt*(2*bcc(2,k,j)+bcc(2,k+1,j))
            endif
            !if(k>kbs(j)+1) then
            if(k>ibtm) then
               rrhs(1,kin)=rrhs(1,kin)+dzz(k)/6*dt*(2*bcc(1,k,j)+bcc(1,k-1,j))
               rrhs(2,kin)=rrhs(2,kin)+dzz(k)/6*dt*(2*bcc(2,k,j)+bcc(2,k-1,j))
            endif
          endif !ibc==0

!         Non-hydrostatic
          if(nonhydro==1) then
            if(k<nvrt) then
              rrhs(1:2,kin)=rrhs(1:2,kin)-dzz(k+1)/6*dt*(2*dqnon_dxy(1:2,k,j)+dqnon_dxy(1:2,k+1,j))
            endif
            !if(k>kbs(j)+1) then
            if(k>ibtm) then
              rrhs(1:2,kin)=rrhs(1:2,kin)-dzz(k)/6*dt*(2*dqnon_dxy(1:2,k,j)+dqnon_dxy(1:2,k-1,j))
            endif
          endif !nonhydro==1

!         Radiation stress
#ifdef  USE_WWM
            if(ics==1) then
              if(k<nvrt) rrhs(1:2,kin)=rrhs(1:2,kin)+dzz(k+1)/6*dt* &
     &(2*wwave_force(1:2,k,j)+wwave_force(1:2,k+1,j))
              !if(k>kbs(j)+1) rrhs(1:2,kin)=rrhs(1:2,kin)+dzz(k)/6*dt* &
              if(k>ibtm) rrhs(1:2,kin)=rrhs(1:2,kin)+dzz(k)/6*dt* &
     &(2*wwave_force(1:2,k,j)+wwave_force(1:2,k-1,j))
            else !use swild10 as approx. to eframe
              call project_hvec(wwave_force(1,k,j),wwave_force(2,k,j), &
     &swild10(1:3,1:3),sframe(:,:,j),swild(1),swild(2))
              if(k<nvrt) then
                call project_hvec(wwave_force(1,k+1,j),wwave_force(2,k+1,j), &
     &swild10(1:3,1:3),sframe(:,:,j),swild(3),swild(4))
                rrhs(1:2,kin)=rrhs(1:2,kin)+dzz(k+1)/6*dt*(2*swild(1:2)+swild(3:4))
              endif !if(k<nvrt)
             
              !if(k>kbs(j)+1) then
              if(k>ibtm) then
                call project_hvec(wwave_force(1,k-1,j),wwave_force(2,k-1,j), &
     &swild10(1:3,1:3),sframe(:,:,j),swild(3),swild(4))
                rrhs(1:2,kin)=rrhs(1:2,kin)+dzz(k)/6*dt*(2*swild(1:2)+swild(3:4))
              endif !if(k>kbs(j)+1)
            endif !ics
#endif /*USE_WWM*/
        enddo !k=ibtm,nvrt

        call tridag(nvrt,100,ndim,2,alow,bdia,cupp,rrhs,soln,gam)
        !do k=kbs(j)+1,nvrt
        do k=ibtm,nvrt
          !kin=k-kbs(j)
          kin=k-ibtm+1
!         Impose limits
          su2(k,j)=max(-rmaxvel,min(rmaxvel,soln(1,kin)))
          sv2(k,j)=max(-rmaxvel,min(rmaxvel,soln(2,kin)))
        enddo !k
!-------------------------------------------------------------------------------------
        endif !lm2d

        if(imm==2) then !no slip
          call update_bdef(time,xcj(j),ycj(j),dep,swild)
          su2(kbs(j),j)=swild(1)
          sv2(kbs(j),j)=swild(2)
        else if(ibottom_bc==1) then
          if(.not.lm2d.and.Cd(j)==0) then
            su2(kbs(j),j)=su2(kbs(j)+1,j)
            sv2(kbs(j),j)=sv2(kbs(j)+1,j)
          else if(.not.lm2d) then !no slip bottom
            su2(kbs(j),j)=0
            sv2(kbs(j),j)=0
          endif
        endif

!       Extend
!        do k=1,kbs_e(j)-1
        do k=1,kbs(j)-1
          su2(k,j)=0 
          sv2(k,j)=0 
        enddo !k

!       Impose uniformity for 2D
        if(lm2d) then
          su2(2,j)=su2(1,j)
          sv2(2,j)=sv2(1,j)
        endif

!	    Impose b.c.
!        do k=kbs_e(j),nvrt
        do k=kbs(j),nvrt
          if(isbs(j)>0.and.ifltype(max(1,isbs(j)))/=0) then !open bnd side
            if(uth(k,j)<-98.or.vth(k,j)<-98) then
              write(errmsg,*)'Wrong vel. input:',uth(k,j),vth(k,j),node1,node2
              call parallel_abort(errmsg)
            endif
            !rotate uth, vth to sframe if ics=2; otherwise same
            uths=uth(k,j); vths=vth(k,j)
            if(ics==2) call project_hvec(uth(k,j),vth(k,j),swild10(1:3,1:3),sframe(:,:,j),uths,vths)

            if(ifltype(isbs(j))==-1) then !Flather 1
              if(eta_mean(node1)<-98.or.eta_mean(node2)<-98) then
                write(errmsg,*)'Flather bnd elevation not assigned:',isbs(j)
                call parallel_abort(errmsg)
              endif
              if(dps(j)<=0) then
                write(errmsg,*)'Flather bnd has negative depth:',isbs(j),dps(j)
                call parallel_abort(errmsg)
              endif

              vnorm=sqrt(grav/dps(j))*(eta2(node1)+eta2(node2)-eta_mean(node1)-eta_mean(node2))/2
              if(ics==1) then
                vnorm=vnorm+uth(k,j)*sframe(1,1,j)+vth(k,j)*sframe(2,1,j)
                su2(k,j)=vnorm*sframe(1,1,j)
                sv2(k,j)=vnorm*sframe(2,1,j)
              else
                vnorm=vnorm+uths
                su2(k,j)=vnorm
                sv2(k,j)=0
              endif !ics
            else if(ifltype(isbs(j))==-4) then !3D radiation
              if(ics==1) then
                vnorm=su2(k,j)*sframe(1,1,j)+sv2(k,j)*sframe(2,1,j)
              else
                vnorm=su2(k,j)
              endif !ics
              if(vnorm<=0) then !incoming
                su2(k,j)=(1-vobc1(isbs(j)))*su2(k,j)+vobc1(isbs(j))*uths 
                sv2(k,j)=(1-vobc1(isbs(j)))*sv2(k,j)+vobc1(isbs(j))*vths 
              else !outgoing
                su2(k,j)=(1-vobc2(isbs(j)))*su2(k,j)+vobc2(isbs(j))*uths 
                sv2(k,j)=(1-vobc2(isbs(j)))*sv2(k,j)+vobc2(isbs(j))*vths 
              endif
            else !not Flather or 3D radiation
              su2(k,j)=uths
              sv2(k,j)=vths
            endif !Flather or not
          endif !open bnd

          if(isbs(j)==-1) then !land bnd
            if(islip==0) then !free slip
              vnorm=0 !for most cases
              !Normal component from vortex formulation
#ifdef USE_WWM
              if(RADFLAG.eq.'VOR') then
                if(ics==1) then
                  vnorm=stokes_vel_sd(1,k,j)*sframe(1,1,j)+stokes_vel_sd(2,k,j)*sframe(2,1,j)
                else
                  call project_hvec(stokes_vel_sd(1,k,j),stokes_vel_sd(2,k,j), &
     &pframe(:,:,isidenode(1,j)),sframe(:,:,j),vnorm,vtmp)
                endif
              endif !RADFLAG
#endif               

              if(ics==1) then
                vtan=su2(k,j)*sframe(1,2,j)+sv2(k,j)*sframe(2,2,j)
                su2(k,j)=vtan*sframe(1,2,j)-vnorm*sframe(1,1,j)
                sv2(k,j)=vtan*sframe(2,2,j)-vnorm*sframe(2,1,j)
              else !lat/lon
                su2(k,j)=-vnorm
              endif !ics
            else !no slip
              su2(k,j)=0
              sv2(k,j)=0
            endif
          endif !land bnd

          !Hydraulic
          if(ihydraulics/=0.and.nhtblocks>0) then; if(isblock_sd(1,j)>0) then
            !Active block
            jblock=isblock_sd(1,j)
            if(isblock_sd(2,j)>0) then !face
              jface=isblock_sd(2,j)
              !Compute normal vel. in local sframe
              dot1=dot_product(dir_block(1:3,jblock),sframe(1:3,1,j))
              ss=sign(1.d0,dot1)
              vnorm=vnth_block(jface,jblock)*ss
              if(ics==1) then
                su2(k,j)=block_nudge*vnorm*sframe(1,1,j)+(1-block_nudge)*swild98(1,k,j) !su2(k,j)
                sv2(k,j)=block_nudge*vnorm*sframe(2,1,j)+(1-block_nudge)*swild98(2,k,j) !sv2(k,j)
              else !lat/lon
                su2(k,j)=block_nudge*vnorm+(1-block_nudge)*swild98(1,k,j) !su2(k,j)
                sv2(k,j)=0
              endif !ics
            else !internal side (for wet/dry) - use face 1 values
              if(ics==1) then
                tmp1=vnth_block(1,jblock)*dir_block(1,jblock)
                tmp2=vnth_block(1,jblock)*dir_block(2,jblock)
              else
                dot1=dot_product(dir_block(1:3,jblock),sframe(1:3,1,j))
                dot2=dot_product(dir_block(1:3,jblock),sframe(1:3,2,j))
                tmp1=vnth_block(1,jblock)*dot1
                tmp2=vnth_block(1,jblock)*dot2
              endif
              su2(k,j)=block_nudge*tmp1+(1-block_nudge)*swild98(1,k,j) !su2(k,j)
              sv2(k,j)=block_nudge*tmp2+(1-block_nudge)*swild98(2,k,j) !sv2(k,j)
            endif !face

            !Check
            !write(12,*)'Vel. b.c:',iplg(isidenode(1:2,j)),k,jblock,jface,real(su2(k,j)),real(sv2(k,j))
          endif; endif !ihydraulics
        enddo !k=kbs(j),nvrt
      enddo !j=1,nsa

!!!!!!!!!!!!!!! modif AD : post-processing
#ifdef USE_WWM
!if (it.eq.ntime) then
!
!  do i=1,npa
!    write(9000,'(I5,3F15.8)')i,xnd(i),ynd(i),eta2(i)
!  enddo
!
!  do i=1,nsa
!  do k=kbs(i),nvrt
!    write(9001,'(2I5,5F15.8)')i,k,xcj(i),ycj(i),zs(k,i),su2(k,i),sv2(k,i)
!!    write(9003,'(2I5,4F15.8)')i,k,xcj(i),ycj(i),zs(k,i),stokes_w_sd(k,i)
!    WRITE(9003,'(2I5,5F15.8)')i,k,xcj(i),ycj(i),zs(k,i),-dJ_dx(i),grav*deta2_dx(i)
!  enddo
!  enddo
!
!
!
!  do i=1,npa
!  do k=kbp(i),nvrt
!     write(9002,'(2I5,6F15.8)')i,k,xnd(i),ynd(i),znl(k,i),stokes_vel(1,k,i),stokes_vel(2,k,i),stokes_w_nd(k,i)
!  enddo
!  enddo
!endif

!JZ: below is only for debugging and should be removed afterwards. So
!I removed ntime
if (it.eq.1000) then
    write(12,*)'OUDELALI OUDELALA LABONNE AVENTURE',1000
    write(12,*)'npa =',np
    write(12,*)'i,xnd(i),ynd(i),eta2(i),dahv_x,dahv_y'
  do i=1,np
    write(12,'(I5,5F18.8)')i,xnd(i),ynd(i),eta2(i),dav(i,1),dav(i,2)
  enddo

    write(12,*)
    write(12,*)'nsa=',ns,'nvrt=',nvrt
    write(12,*)'i,k,xcj(i),ycj(i),zs(k,i),su2(k,i),sv2(k,i),-dJ_dx(i),grav*deta2_dx(i)'

  do i=1,ns
  do k=kbs(i),nvrt
    write(12,'(2I5,7F18.8)')i,k,xcj(i),ycj(i),zs(k,i),su2(k,i),sv2(k,i),-dJ_dx(i),grav*deta2_dx(i)
!    write(9003,'(2I5,4F15.8)')i,k,xcj(i),ycj(i),zs(k,i),stokes_w_sd(k,i)
!    WRITE(9003,'(2I5,5F15.8)')i,k,xcj(i),ycj(i),zs(k,i),-dJ_dx(i),grav*deta2_dx(i)
  enddo
  enddo


    write(12,*)
    write(12,*)'npa=',np,'nvrt=',nvrt
    write(12,*)'i,k,xnd(i),ynd(i),znl(k,i),stokes_vel(1,k,i),stokes_vel(2,k,i),stokes_w_nd(k,i)'
  do i=1,np
  do k=kbp(i),nvrt
     if ((stokes_vel(1,k,i).le.-100.).or.(stokes_vel(1,k,i).ge.100.).or.(stokes_vel(2,k,i).le.-100.).or.(stokes_vel(2,k,i).ge.100.)) then
     stokes_vel(1,k,i)=0.
     stokes_vel(2,k,i)=0.
     stokes_w_nd(k,i)=0.
     endif
     write(12,'(2I5,6F30.8)')i,k,xnd(i),ynd(i),znl(k,i),stokes_vel(1,k,i),stokes_vel(2,k,i),stokes_w_nd(k,i)
  enddo
  enddo
endif

#endif /*USE_WWM*/
!!!!!!!!!!!!!!! end modif AD

      deallocate(swild98)

!...  Shapiro filter (used only if indvel<=0)
!     use bcc as temporary variable (sframe)
      if(indvel<=0) then
        bcc=0
        do i=1,ns !residents only
          if(isdel(2,i)==0.or.idry_s(i)==1) cycle
          if(ihydraulics/=0.and.nhtblocks>0) then
            if(isblock_sd(1,i)/=0) cycle
          endif

!         Define bottom level
          if(ibottom_bc==1) then !old @kbs+1
            ibtm=kbs(i)+1
          else !new @ kbs
            ibtm=kbs(i)
          endif !ibottom_bc

!         Internal wet sides
          !do k=kbs(i)+1,nvrt
          do k=ibtm,nvrt
            suru=0
            surv=0
            do j=1,4
              id=isidenei2(j,i)
              if(idry_s(id)==1) then
                kin=k
              else
                kin=max(k,kbs(id)+1)
              endif
              if(ics==1) then
                utmp=su2(kin,id)
                vtmp=sv2(kin,id)
              else !lat/lon
                call project_hvec(su2(kin,id),sv2(kin,id),sframe(:,:,id),sframe(:,:,i),utmp,vtmp)
              endif !ics
              suru=suru+utmp
              surv=surv+vtmp
            enddo !j
            bcc(1,k,i)=su2(k,i)+shapiro/4*(suru-4*su2(k,i)) !sframe if ics=2
            bcc(2,k,i)=sv2(k,i)+shapiro/4*(surv-4*sv2(k,i))
          enddo !k
        enddo !i=1,ns

        do j=1,ns
          if(isdel(2,j)==0.or.idry_s(j)==1) cycle 
          if(ihydraulics/=0.and.nhtblocks>0) then
            if(isblock_sd(1,j)/=0) cycle
          endif

!         Define bottom level
          if(ibottom_bc==1) then !old @kbs(j)+1
            ibtm=kbs(j)+1
          else !new @ kbs
            ibtm=kbs(j)
          endif !ibottom_bc

          !do k=kbs(j)+1,nvrt
          do k=ibtm,nvrt
            su2(k,j)=bcc(1,k,j)
            sv2(k,j)=bcc(2,k,j)
          enddo !k

!         2D
          if(lm2d) then
            su2(1,j)=su2(nvrt,j)
            sv2(1,j)=sv2(nvrt,j)
          endif

!          do k=1,kbs_e(j)-1
          do k=1,kbs(j)-1
            su2(k,j)=0
            sv2(k,j)=0
          enddo !k
        enddo !j=1,ns

!       Exchange ghosts
        allocate(swild98(2,nvrt,nsa),stat=istat)
        if(istat/=0) call parallel_abort('MAIN: fail to allocate swild98')
!'
        swild98(1,:,:)=su2(:,:)
        swild98(2,:,:)=sv2(:,:)
#ifdef INCLUDE_TIMING
        cwtmp=mpi_wtime()
#endif
        call exchange_s3d_2(swild98)
#ifdef INCLUDE_TIMING
        wtimer(8,2)=wtimer(8,2)+mpi_wtime()-cwtmp
#endif
        su2(:,:)=swild98(1,:,:)
        sv2(:,:)=swild98(2,:,:)
        deallocate(swild98)
      endif !indvel<=0; Shapiro filter

      if(myrank==0) write(16,*)'done solving momentum eq...'

!**********************************************************************************
!     Non-hydrostatic part
!**********************************************************************************
      if(nonhydro==1) then
        if(ivcor/=2) call parallel_abort('MAIN: nonhydro cannot use other z-coor.')
        if(myrank==0) write(16,*)'start non-hydrostatic calculation...'
!'...    Solve for intermediate we
        hp_int(:,:,2)=we(:,:) !previous step temporarily saved as hp_int(:,:,2)
        we=0
        do i=1,nea
          if(idry_e(i)==1) cycle

!         Wet element
          n1=elnode(1,i); n2=elnode(2,i); n3=elnode(3,i)
!         Define layer thickness & diffusivities
          do k=kbe(i)+1,nvrt
            dzz(k)=ze(k,i)-ze(k-1,i)
            if(dzz(k)<=0) call parallel_abort('MAIN: dzz=0 in wvel')
            dfz(k)=(ptbt(3,k,n1)+ptbt(3,k,n2)+ptbt(3,k,n3)+ &
     &ptbt(3,k-1,n1)+ptbt(3,k-1,n2)+ptbt(3,k-1,n3))/6
          enddo !k

!         Coefficient matrix
          ndim=nvrt-kbe(i)+1
          alow(1)=0; bdia(1)=1; cupp(1)=0
          alow(ndim)=0; bdia(ndim)=1; cupp(ndim)=0
          do k=kbe(i)+1,nvrt-1
            kin=k-kbe(i)+1 !eq. #
            tmp1=dt*dfz(k+1)/dzz(k+1)
            tmp2=dt*dfz(k)/dzz(k)
            alow(kin)=dzz(k)/6-tmp2
            bdia(kin)=(dzz(k)+dzz(k+1))/3+tmp2+tmp1
            cupp(kin)=dzz(k+1)/6-tmp1
          enddo !k

!         RHS
!         b.c. first
          dhdx=0; dhdy=0; ubar1=0; vbar1=0; ubar2=0; vbar2=0; 
          eta1_bar=0; eta2_bar=0; swild(1:4)=0 !1:2 --> deta1_dxy; 3:4 -->deta2_dxy
          ifl=0 !flag for wet/dry for step n+1
          do j=1,3
            nd=elnode(j,i); isd=elside(j,i)
            if(eta2(nd)+dp(nd)<=h0) ifl=1
            dhdx=dhdx+dp(nd)*dldxy(j,1,i)
            dhdy=dhdy+dp(nd)*dldxy(j,2,i)
            ubar1=ubar1+su2(kbs(isd),isd)/3
            vbar1=vbar1+sv2(kbs(isd),isd)/3
            swild(1)=swild(1)+eta1(nd)*dldxy(j,1,i) !deta1_dx
            swild(2)=swild(2)+eta1(nd)*dldxy(j,2,i) !deta1_dy
            swild(3)=swild(3)+eta2(nd)*dldxy(j,1,i) !deta2_dx
            swild(4)=swild(4)+eta2(nd)*dldxy(j,2,i) !deta2_dy
            ubar2=ubar2+su2(nvrt,isd)/3
            vbar2=vbar2+sv2(nvrt,isd)/3
            eta1_bar=eta1_bar+eta1(nd)/3
            eta2_bar=eta2_bar+eta2(nd)/3
          enddo !j

          if(imm==2) then
            call update_bdef(time,xctr(i),yctr(i),dep,swild)
            ubed=swild(1); vbed=swild(2); wbed=swild(3)
            !vnorm=ubed*dhdx+vbed*dhdy+wbed
            rrhs(1,1)=wbed
          else !imm
            rrhs(1,1)=-dhdx*ubar1-dhdy*vbar1 !no deformation
          endif

          if(ifl==0) then
            rrhs(1,ndim)=(ubar2*(thetai*swild(3)+(1-thetai)*swild(1))+ &
     &vbar2*(thetai*swild(4)+(1-thetai)*swild(2))+(eta2_bar-eta1_bar)/dt- &
     &(1-thetai)*hp_int(nvrt,i,2))/thetai
          else
            rrhs(1,ndim)=0
          endif

!         Middle levels
          do k=kbe(i)+1,nvrt-1
            kin=k-kbe(i)+1
            qnon_e1=0; qnon_e2=0
            do j=1,3
              isd=elside(j,i)
              nd=elnode(j,i)
              qnon_e1=qnon_e1+qnon(k-1,nd)/3
              qnon_e2=qnon_e2+qnon(k+1,nd)/3
            enddo !j
            rrhs(1,kin)=dzz(k+1)/6*(2*webt(k,i)+webt(k+1,i))+dzz(k)/6*(2*webt(k,i)+webt(k-1,i))- &
     &dt/2*(qnon_e2-qnon_e1)
          enddo !k

          call tridag(nvrt,100,ndim,1,alow,bdia,cupp,rrhs,soln,gam)

!         Below bottom we=0 already assigned
          do k=kbe(i),nvrt
            kin=k-kbe(i)+1
            we(k,i)=soln(1,kin)
          enddo !k
        enddo !i=1,nea

        if(myrank==0) write(16,*)'done solving interim vertical vel...'

!...    Compute averaged divergence in each prism (temporarily saved as hp_int(nvrt,nea,1))
        hp_int(:,:,1)=0
        do i=1,nea
          if(idry_e(i)==1) cycle

!	      Wet elements with 3 wet nodes
!	      Compute upward normals and areas @ all levels
          n1=elnode(1,i)
          n2=elnode(2,i)
          n3=elnode(3,i)
          !av_bdef1=(bdef1(n1)+bdef1(n2)+bdef1(n3))/3 !average bed deformation
          !av_bdef2=(bdef2(n1)+bdef2(n2)+bdef2(n3))/3
          if(kbe(i)==0) then
            write(errmsg,*)'Impossible 95 (2)'
            call parallel_abort(errmsg)
          endif
          do l=kbe(i),nvrt
            xcon=(ynd(n2)-ynd(n1))*(znl(max(l,kbp(n3)),n3)-znl(max(l,kbp(n1)),n1))-(ynd(n3)-ynd(n1))* &
     &(znl(max(l,kbp(n2)),n2)-znl(max(l,kbp(n1)),n1))
            ycon=(xnd(n3)-xnd(n1))*(znl(max(l,kbp(n2)),n2)-znl(max(l,kbp(n1)),n1))-(xnd(n2)-xnd(n1))* &
     &(znl(max(l,kbp(n3)),n3)-znl(max(l,kbp(n1)),n1))
            zcon=area(i)*2
            area_e(l)=sqrt(xcon**2+ycon**2+zcon**2)/2
            if(area_e(l)==0) then
              write(errmsg,*)'Zero area (3):',i,l
              call parallel_abort(errmsg)
            endif
            sne(1,l)=xcon/area_e(l)/2
            sne(2,l)=ycon/area_e(l)/2
            sne(3,l)=zcon/area_e(l)/2 !>0
          enddo !l

          do l=kbe(i)+1,nvrt
            sum1=0 !all horizontal fluxes
            ubar=0 !av. vel. at level l
            vbar=0
            ubar1=0 !av. vel. at level l-1
            vbar1=0
            do j=1,3
              jsj=elside(j,i)
              vnor1=su2(l,jsj)*sframe(1,1,jsj)+sv2(l,jsj)*sframe(2,1,jsj)
              vnor2=su2(l-1,jsj)*sframe(1,1,jsj)+sv2(l-1,jsj)*sframe(2,1,jsj)
!              if(l-1<kbs(jsj).or.kbs(jsj)==0) then
!                write(errmsg,*)'Impossible 94 (2):',l,kbs(jsj),ielg(i)
!                call parallel_abort(errmsg)
!              endif
              sum1=sum1+ssign(j,i)*(zs(max(l,kbs(jsj)),jsj)-zs(max(l-1,kbs(jsj)),jsj))*distj(jsj)*(vnor1+vnor2)/2

              ubar=ubar+su2(l,jsj)/3    
              ubar1=ubar1+su2(l-1,jsj)/3    
              vbar=vbar+sv2(l,jsj)/3    
              vbar1=vbar1+sv2(l-1,jsj)/3    
            enddo !j=1,3

!           Impose b.c.
            if(l==kbe(i)+1) then
              if(imm==2) then
                call update_bdef(time,xctr(i),yctr(i),dep,swild)
                ubed=swild(1); vbed=swild(2); wbed=swild(3)
                bflux=ubed*sne(1,l-1)+vbed*sne(2,l-1)+wbed*sne(3,l-1)
              else
                bflux=0
              endif
            else
              bflux=ubar1*sne(1,l-1)+vbar1*sne(2,l-1)+we(l-1,i)*sne(3,l-1) !inner normal
            endif
            top=ubar*sne(1,l)+vbar*sne(2,l)+we(l,i)*sne(3,l)

            hp_int(l,i,1)=(sum1+top*area_e(l)-bflux*area_e(l-1))/area(i)/(ze(l,i)-ze(l-1,i))
          enddo !l=kbe(i)+1,nvrt
        enddo !i=1,nea

        if(myrank==0) write(16,*)'done computing divergence...'

!...    Solve for qhat - correction to non-hydrostatic pressure
!       Prepare derivatives for S prisms at each node (for dz_ds) and element (for \grad sig)
!       dz/dsigma is (temporarily) webt(nvrt,npa), and dsigma/d[x,y] are sdbt(1:2,nvrt,nea)
        sdbt(1:2,:,:)=-1.e25 !flag
        webt(:,:)=-1.e25 
        do i=1,npa
          if(idry(i)==1) cycle

!         Wet node
          if(kbp(i)>kz) call parallel_abort('MAIN: node level wrong')
          do k=kz,nvrt !no Z
            kin=k-kz+1
            if(hmod(i)<=h_c) then !traditional
              webt(k,i)=dp(i)+eta1(i) !use old levels
              if(webt(k,i)<=0) then
                write(errmsg,*)'MAIN: dzds<=0: ',iplg(i),webt(k,i)
                call parallel_abort(errmsg)
              endif
            else !S
              webt(k,i)=eta1(i)+h_c+(hmod(i)-h_c)*dcs(kin)
              if(webt(k,i)<=0) then
                write(errmsg,*)'MAIN: dzds<=0: ',iplg(i),webt(k,i)
                call parallel_abort(errmsg)
              endif
            endif
          enddo !k=kz,nvrt
        enddo !i=1,npa

        do i=1,nea
          if(idry_e(i)==1) cycle

!         Wet element
          if(kbe(i)>kz) call parallel_abort('MAIN: elem. level wrong')
          deta_dx=0; deta_dy=0; dhdx=0; dhdy=0
          hmin=1.e25
          do j=1,3
            nd=elnode(j,i)
            dhdx=dhdx+hmod(nd)*dldxy(j,1,i)
            dhdy=dhdy+hmod(nd)*dldxy(j,2,i)
            deta_dx=deta_dx+eta1(nd)*dldxy(j,1,i) !use old levels         
            deta_dy=deta_dy+eta1(nd)*dldxy(j,2,i) !use old levels         
            if(hmod(nd)<hmin) hmin=hmod(nd)
          enddo !j

          do k=kz,nvrt !no Z
            kin=k-kz+1
            dzds_av=0 !average dz_ds
            do j=1,3
              nd=elnode(j,i)
              if(webt(k,nd)<-1.e24) call parallel_abort('MAIN: impossible 102')
!'
              dzds_av=dzds_av+webt(k,nd)/3
            enddo !j
            if(dzds_av<=0) call parallel_abort('MAIN: impossible 103')

            if(hmin<=h_c) then !traditional
              css=sigma(kin)
            else !S
              css=cs(kin)
            endif

            sdbt(1,k,i)=-(deta_dx*(1+sigma(kin))+css*dhdx)/dzds_av
            sdbt(2,k,i)=-(deta_dy*(1+sigma(kin))+css*dhdy)/dzds_av
          enddo !k=kz,nvrt
        enddo !i=1,nea

        if(myrank==0) write(16,*)'done preparing derivatives...'

!       Construct matrix qmatr, RHS (qir)
!       Valid range for qmatr: (kbp_e(i):nvrt,-1:1,0:(mnei+1),i=1:np), and node i is wet; -1:1
!       represents levels k-1,k and k+1. So for each pt in 3D space (node and a level),
!       each eq. will have up to 3*(1+nnp(i)) unknowns.
        kbp_e=nvrt+1 !min. kbe from sourrounding _wet_ element for a _resident_ wet node
        qmatr=0; qir=0
        do i=1,np !resident
          if(idry(i)==1) cycle

!         Wet node
          do j=1,nne(i)
            ie=indel(j,i)
            if(idry_e(ie)==0.and.kbp_e(i)>kbe(ie)) kbp_e(i)=kbe(ie) 
          enddo !j
          if(kbp_e(i)==nvrt+1) call parallel_abort('MAIN: impossible 04')

          do k=kbp_e(i),nvrt
            icount=0 !# of wet elements
            do j=1,nne(i)
              ie=indel(j,i)
              if(idry_e(ie)==1) cycle

              !Wet element ie
              icount=icount+1
              id=iself(j,i)
              do l=0,1 !two prisms (ie,k+l)
                if(k+l>=kbe(ie)+1.and.k+l<=nvrt) then !prism exists
                  do m=1,3 !nodes
                    if(m==1) then
                      nd=i; ind=0 !index in the ball
                      ind2=id !index in the element
                    else
                      nd=elnode(nx(id,m-1),ie); ind=j+m-2
                      if(m==3.and.isbnd(1,i)==0.and.j==nne(i)) ind=1
                      ind2=nx(id,m-1)
                    endif

                    do mk=-1,0 !two levels k+l+mk; inside loop the index "l=(l1,l2)" (in the notes) is fixed
                      !\hat{I_0} & \hat{I_r}
                      dot1=dldxy(id,1,ie)*dldxy(ind2,1,ie)+dldxy(id,2,ie)*dldxy(ind2,2,ie) !dot product of gradient of shape function
                      if(k+l>=kz+1) then !S prism
                        dsigma=sigma(k+l-kz+1)-sigma(k+l-kz)
                        if(dsigma<=0) call parallel_abort('MAIN:dsig<=0')
                        dgam0=sign(1.,0.5-l)/dsigma !d{gamma} for \bar{l}_2
                        dgam1=sign(1.,0.5+mk)/dsigma !!d{gamma} for l_2
                        vol=area(ie)*dsigma

                        hat_i0=0
                        do mm=1,3 !nodes
                          nd_lam=elnode(mm,ie)
                          if(i==nd.and.i==nd_lam) then
                            iee=6
                          else if(i/=nd.and.i/=nd_lam.and.nd/=nd_lam) then
                            iee=1
                          else !2 equal
                            iee=2
                          endif
                          do mmk=-1,0 !two levels k+l+mmk for lambda (in the notes)
                            dzds=webt(k+l+mmk,nd_lam) !shorthand
                            dsdx=sdbt(1,k+l+mmk,ie)
                            dsdy=sdbt(2,k+l+mmk,ie)
                            if(min(dzds,dsdx,dsdy)<-1.e24.or.dzds<=0) then
                              write(errmsg,*)'MAIN: wrong derivatives:',iplg(i),ielg(ie),dzds,dsdx,dsdy
                              call parallel_abort(errmsg)
                            endif
                            dsig2=dsdx*dsdx+dsdy*dsdy !magnitude squared
                            if(k==k+l+mk.and.k==k+l+mmk) then
                              idel=1
                            else
                              idel=0
                            endif
                            hat_i0=hat_i0+dot1*dzds*(1+2*idel)/36+dgam0*dgam1*iee/120*(dsig2*dzds+1/dzds)+ &
     &dgam0*(dldxy(ind2,1,ie)*dsdx+dldxy(ind2,2,ie)*dsdy)*dzds/72*(1+kronecker(k+l+mk,k+l+mmk))*(1+kronecker(i,nd_lam))+&
     &dgam1*(dldxy(id,1,ie)*dsdx+dldxy(id,2,ie)*dsdy)*dzds/72*(1+kronecker(k,k+l+mmk))*(1+kronecker(nd,nd_lam))
                          enddo !mmk=-1,0; two levels
                        enddo !mm; 3 nodes
                        hat_i0=hat_i0*vol

                        !hat_ir is the term inside the summation (in the notes)
                        hat_ir=vol/72*(1+kronecker(i,nd))*(1+kronecker(k,k+l+mk))*webt(k+l+mk,nd)
                      else !Z prism
                        dz=ze(k+l,ie)-ze(k+l-1,ie)
                        if(dz<=0) call parallel_abort('MAIN: dz<=0')
                        dgam0=sign(1.,0.5-l)/dz
                        dgam1=sign(1.,0.5+mk)/dz
                        vol=area(ie)*dz
                        hat_i0=vol*((1.+kronecker(k,k+l+mk))/6*dot1+(1.+kronecker(i,nd))/12*dgam0*dgam1)
                        hat_ir=vol/36
                      endif !S or Z prism

                      if(ind==0.and.l+mk==0.and.hat_i0<=0) then !diagonal
                        write(errmsg,*)'MAIN: diagonal problem ',iplg(i),j,k,l,mk,m,k+l,kz+1,hat_i0
                        call parallel_abort(errmsg)
                      endif
   
                      qmatr(k,l+mk,ind,i)=qmatr(k,l+mk,ind,i)+hat_i0
                      qir(k,i)=qir(k,i)-hat_ir*hp_int(k+l,ie,1)/dt/thetai
                    enddo !mk=-1,0; two levels
                  enddo !m; 3 nodes
                endif !valid prism
              enddo !l=0,1 prisms
            enddo !j=1,nne(i)
            if(icount==0) call parallel_abort('MAIN: no wet element(5)')
          enddo !k=kbp_e(i),nvrt
        enddo !i=1,np

        if(myrank==0) write(16,*)'done preparing q-matrix...'

!...    Check diagonal, symmetry etc (the latter is time and memory consumimg)
        if(nproc==1) then
          allocate(swild99(nvrt*npa,nvrt*npa),stat=istat)
          if(istat/=0) call parallel_abort('MAIN: failed to allocate swild99 for symmetry check')
!'
          swild99=0
        endif

        tmp_max=0
        dia_min=1.e25 !min. of diagonal
        do i=1,np !resident
          if(idry(i)==1) cycle

          do k=kbp_e(i),nvrt
            sum1=0 !sum of all entries
            irow=(i-1)*nvrt+k !eq. #; note that some eqs. are 0=0 in global matrix (e.g., below bottom etc)
            do kk=-1,1
              if(k+kk<kbp_e(i).or.k+kk>nvrt) cycle

              do j=0,nnp(i)
                if(j==0) then
                  nd=i
                else
                  nd=indnd(j,i)
                endif
                icol=(nd-1)*nvrt+k+kk 
                if(nproc==1) swild99(irow,icol)=qmatr(k,kk,j,i)
                sum1=sum1+qmatr(k,kk,j,i)
              enddo !j=0,nnp(i)
              !write(12,*)i,k,kk,qmatr(k,kk,0:nnp(i),i)/qmatr(k,0,0,i)
            enddo !kk

            !check diagonal
            if(qmatr(k,0,0,i)<=0) then
              write(errmsg,*)'MAIN: qmatr diagonal<=0',iplg(i),k,qmatr(k,0,0,i),i,kbp_e(i)
              call parallel_abort(errmsg)
            endif
            tmp=sum1/qmatr(k,0,0,i)
            tmp_max=max(tmp_max,abs(tmp))
            if(tmp_max>1.e-5) then
              write(errmsg,*)'MAIN: sum/=0',iplg(i),k,qmatr(k,0,0,i),i,kbp_e(i),tmp_max
              call parallel_abort(errmsg)
            endif
            dia_min=min(dia_min,qmatr(k,0,0,i))
            !write(12,*)i,k,tmp,qmatr(k,0,0,i)
          enddo !k=kbp_e(i),nvrt
        enddo !i=1,np
        call mpi_reduce(tmp_max,tmp_max_gb,1,rtype,MPI_MAX,0,comm,ierr)
        call mpi_reduce(dia_min,dia_min_gb,1,rtype,MPI_MIN,0,comm,ierr)
        if(myrank==0) then
          write(16,*)'Max. sum/diagonal in qmatr= ',tmp_max_gb
          write(16,*)'Min. diagonal in qmatr= ',dia_min_gb
        endif

!...    To check symmetry please use nproc=1
        if(nproc==1) then
          df_max=-1
          do i=1,nvrt*np
            !if(swild99(i,i)<=0) call parallel_abort('MAIN: wrong diagnl')
            do j=1,nvrt*np
              tmp=abs(swild99(i,j)-swild99(j,i))
              if(tmp>1.e-4) then
                ieq=i/nvrt+1; k1=i-(ieq-1)*nvrt
                if(k1==0) then
                  ieq=ieq-1; k1=nvrt
                endif
                ij=j/nvrt+1; k2=j-(ij-1)*nvrt
                if(k2==0) then
                  ij=ij-1; k2=nvrt
                endif
                write(errmsg,*)'MAIN: qmatr not symmetric:',ieq,k1,ij,k2,swild99(i,j),swild99(j,i)
                call parallel_abort(errmsg)
              endif
              df_max=max(df_max,tmp)
            enddo !j
          enddo !i
          if(myrank==0) write(16,*)'Max. asymmetry=',df_max
          deallocate(swild99)
        endif !nproc==1

        if(myrank==0) write(16,*)'done checking q-matrix...'

!...    CG solver
!       No exchange for qhat afterwards so consistency may not be guarenteed
        call solve_jcg_qnon(it,moitn0,mxitn0,rtol0,nvrt,mnei,np,npa,ihydro,qmatr,qhat,qir)
        if(myrank==0) write(16,*)'done JCG solver...'
        
!...    Update qnon - non-hydrostatic pressure at n+1/2
        qnon=qnon+qhat

!...    Compute gradient of qhat
        call hgrad_nodes(1,0,nvrt,npa,nsa,qhat,dqnon_dxy)

!...    Solve for final u,v
        su2(:,1:ns)=su2(:,1:ns)-dt*thetai*dqnon_dxy(1,:,1:ns)
        sv2(:,1:ns)=sv2(:,1:ns)-dt*thetai*dqnon_dxy(2,:,1:ns)
        do j=1,ns
          if(idry_s(j)==1) cycle
          !do k=1,kbs_e(j)-1
          do k=1,kbs(j)-1
            su2(k,j)=0
            sv2(k,j)=0
          enddo !k
        enddo !j

!       Exchange ghosts
        allocate(swild98(2,nvrt,nsa),stat=istat)
        if(istat/=0) call parallel_abort('MAIN: fail to allocate swild98')
!'
        swild98(1,:,:)=su2(:,:)
        swild98(2,:,:)=sv2(:,:)
#ifdef INCLUDE_TIMING
        cwtmp=mpi_wtime()
#endif
        call exchange_s3d_2(swild98)
#ifdef INCLUDE_TIMING
        wtimer(8,2)=wtimer(8,2)+mpi_wtime()-cwtmp
#endif
        su2(:,:)=swild98(1,:,:)
        sv2(:,:)=swild98(2,:,:)
        deallocate(swild98)

!...    Update vertical velocity
        do i=1,nea
          if(idry_e(i)==1) cycle
          we(1:kbe(i)-1,i)=0

!         b.c. first
          dhdx=0; dhdy=0; ubar1=0; vbar1=0; ubar2=0; vbar2=0;
          eta1_bar=0; eta2_bar=0; swild(1:4)=0 !1:2 --> deta1_dxy; 3:4 -->deta2_dxy
          ifl=0 !flag for wet/dry for step n+1
          do j=1,3
            nd=elnode(j,i); isd=elside(j,i)
            if(eta2(nd)+dp(nd)<=h0) ifl=1
            dhdx=dhdx+dp(nd)*dldxy(j,1,i)
            dhdy=dhdy+dp(nd)*dldxy(j,2,i)
            ubar1=ubar1+su2(kbs(isd),isd)/3
            vbar1=vbar1+sv2(kbs(isd),isd)/3
            swild(1)=swild(1)+eta1(nd)*dldxy(j,1,i) !deta1_dx
            swild(2)=swild(2)+eta1(nd)*dldxy(j,2,i) !deta1_dy
            swild(3)=swild(3)+eta2(nd)*dldxy(j,1,i) !deta2_dx
            swild(4)=swild(4)+eta2(nd)*dldxy(j,2,i) !deta2_dy
            ubar2=ubar2+su2(nvrt,isd)/3
            vbar2=vbar2+sv2(nvrt,isd)/3
            eta1_bar=eta1_bar+eta1(nd)/3
            eta2_bar=eta2_bar+eta2(nd)/3
          enddo !j

          if(imm==2) then
            call update_bdef(time,xctr(i),yctr(i),dep,swild)
            ubed=swild(1); vbed=swild(2); wbed=swild(3)
            !vnorm=ubed*dhdx+vbed*dhdy+wbed
            we(kbe(i),i)=wbed
          else
            we(kbe(i),i)=-dhdx*ubar1-dhdy*vbar1 !no deformation
          endif

          if(ifl==0) then
            we(nvrt,i)=(ubar2*(thetai*swild(3)+(1-thetai)*swild(1))+ &
     &vbar2*(thetai*swild(4)+(1-thetai)*swild(2))+(eta2_bar-eta1_bar)/dt- &
     &(1-thetai)*hp_int(nvrt,i,2))/thetai
          else
            we(nvrt,i)=0
          endif

          do k=kbe(i)+1,nvrt-1
            qhat_e1=0; qhat_e2=0
            do j=1,3
              nd=elnode(j,i)
              qhat_e1=qhat_e1+qhat(k-1,nd)/3
              qhat_e2=qhat_e2+qhat(k+1,nd)/3
            enddo !j
            dqdz=(qhat_e2-qhat_e1)/(ze(k+1,i)-ze(k-1,i))
            we(k,i)=we(k,i)-dt*thetai*dqdz
          enddo !k

!          Cubic spline
!          swild=0
!          do k=kbe(i),nvrt
!            do j=1,3
!              nd=elnode(j,i)
!              swild(k)=swild(k)+qhat(k,nd)/3
!            enddo !j
!          enddo !k
!          call cubic_spline(nvrt-kbe(i)+1,ze(kbe(i):nvrt,i),swild(kbe(i):nvrt),0._rkind,0._rkind,swild2(kbe(i):nvrt,1))
!          do k=kbe(i)+1,nvrt-1
!            dqdz=(swild(k+1)-swild(k))/(ze(k+1,i)-ze(k,i))-(ze(k+1,i)-ze(k,i))/6*(2*swild2(k,1)+swild2(k+1,1))
!            we(k,i)=we(k,i)-dt*thetai*dqdz
!          enddo !k
        enddo !i=nea

        if(myrank==0) write(16,*)'finished non-hydrostatic calculation'
      endif !nonhydro==1
!**********************************************************************************
!     End non-hydrostatic part
!**********************************************************************************

!...  Sponge layer for elev. and vel.
      if(inu_elev==1) then
        do i=1,npa
          eta2(i)=eta2(i)*(1-elev_nudge(i)*dt)
        enddo !i
      endif !inu_elev

      if(inu_uv==1) then
        do i=1,nsa
          uvnu=(uv_nudge(isidenode(1,i))+uv_nudge(isidenode(2,i)))/2*dt
          su2(:,i)=su2(:,i)*(1-uvnu)
          sv2(:,i)=sv2(:,i)*(1-uvnu)
        enddo !i
      endif !inu_uv
           
!...  solve for vertical velocities using F.V.
!...  For hydrostatic model, this is the vertical vel; for non-hydrostatic
!...  model, this is only used in transport
!     swild98 for storing rotated hvel at 3 sides
      allocate(swild98(2,3,nvrt),stat=istat)
      we_fv=0 !for dry and below bottom levels; in eframe if ics=2
      do i=1,nea
        if(idry_e(i)==1) cycle

!	Wet elements with 3 wet nodes
!	Compute upward normals and areas @ all levels
        n1=elnode(1,i)
        n2=elnode(2,i)
        n3=elnode(3,i)
        av_bdef1=(bdef1(n1)+bdef1(n2)+bdef1(n3))/3 !average bed deformation
        av_bdef2=(bdef2(n1)+bdef2(n2)+bdef2(n3))/3
        if(kbe(i)==0) then
          write(errmsg,*)'Impossible 95'
          call parallel_abort(errmsg)
        endif
        do l=kbe(i),nvrt
          if(ics==1) then
            xcon=(ynd(n2)-ynd(n1))*(znl(max(l,kbp(n3)),n3)-znl(max(l,kbp(n1)),n1))-(ynd(n3)-ynd(n1))* &
     &(znl(max(l,kbp(n2)),n2)-znl(max(l,kbp(n1)),n1))
            ycon=(xnd(n3)-xnd(n1))*(znl(max(l,kbp(n2)),n2)-znl(max(l,kbp(n1)),n1))-(xnd(n2)-xnd(n1))* &
     &(znl(max(l,kbp(n3)),n3)-znl(max(l,kbp(n1)),n1))
            zcon=area(i)*2
          else !lat/lon
            !eframe
            call cross_product(xel(2,i)-xel(1,i),yel(2,i)-yel(1,i),znl(max(l,kbp(n2)),n2)-znl(max(l,kbp(n1)),n1), &
     &                         xel(3,i)-xel(1,i),yel(3,i)-yel(1,i),znl(max(l,kbp(n3)),n3)-znl(max(l,kbp(n1)),n1), &
     &                         xcon,ycon,zcon)
          endif !ics

          area_e(l)=sqrt(xcon**2+ycon**2+zcon**2)/2
          if(area_e(l)==0) then
            write(errmsg,*)'Zero area:',i,l
            call parallel_abort(errmsg)
          endif
          sne(1,l)=xcon/area_e(l)/2 !in eframe
          sne(2,l)=ycon/area_e(l)/2
          sne(3,l)=zcon/area_e(l)/2 !>0
        enddo !l

!       Rotate hvel. for 3 sides at all levels
        ubar=0; vbar=0 !average bottom hvel
        do m=1,3 !side
          isd=elside(m,i)
          do k=1,nvrt !cover all
            if(ics==1) then
              swild98(1,m,k)=su2(k,isd)
              swild98(2,m,k)=sv2(k,isd)
            else !to eframe
              call project_hvec(su2(k,isd),sv2(k,isd),sframe(:,:,isd), &
     &eframe(:,:,i),swild98(1,m,k),swild98(2,m,k))
            endif !ics
          enddo !k
          ubar=ubar+swild98(1,m,kbs(isd))/3
          vbar=vbar+swild98(2,m,kbs(isd))/3
        enddo !m=1,3

!       Bottom b.c.
        dhdx=dp(n1)*dldxy(1,1,i)+dp(n2)*dldxy(2,1,i)+dp(n3)*dldxy(3,1,i) !eframe
        dhdy=dp(n1)*dldxy(1,2,i)+dp(n2)*dldxy(2,2,i)+dp(n3)*dldxy(3,2,i)
!        ubar=(su2(kbs(elside(1,i)),elside(1,i))+su2(kbs(elside(2,i)),elside(2,i))+su2(kbs(elside(3,i)),elside(3,i)))/3
!        vbar=(sv2(kbs(elside(1,i)),elside(1,i))+sv2(kbs(elside(2,i)),elside(2,i))+sv2(kbs(elside(3,i)),elside(3,i)))/3
!       we_fv=0 unless Cd=0
!        if(nonhydro==1) then
!          we_fv(kbe(i),i)=(av_bdef2-av_bdef1)/dt-dhdx*ubar-dhdy*vbar
!        else
!          we_fv(kbe(i),i)=(av_bdef2-av_bdef1)/dt !-dhdx*ubar-dhdy*vbar
!        endif
        if(imm==2) then
          call update_bdef(time,xctr(i),yctr(i),dep,swild)
          ubed=swild(1); vbed=swild(2); wbed=swild(3)
          bflux0=ubed*sne(1,kbe(i))+vbed*sne(2,kbe(i))+wbed*sne(3,kbe(i)) !normal bed vel.
          we_fv(kbe(i),i)=wbed
        else
          we_fv(kbe(i),i)=(av_bdef2-av_bdef1)/dt-dhdx*ubar-dhdy*vbar
        endif

        do l=kbe(i),nvrt-1
          sum1=0
          ubar=0
          vbar=0
          ubar1=0
          vbar1=0
          do j=1,3
            jsj=elside(j,i)
            if(ics==1) then
              vnor1=su2(l,jsj)*sframe(1,1,jsj)+sv2(l,jsj)*sframe(2,1,jsj)
              vnor2=su2(l+1,jsj)*sframe(1,1,jsj)+sv2(l+1,jsj)*sframe(2,1,jsj)
            else !lat/lon
              vnor1=su2(l,jsj) !normal
              vnor2=su2(l+1,jsj)
            endif !ics
!            if(l<kbs(jsj).or.kbs(jsj)==0) then
!              write(errmsg,*)'Impossible 94:',l,kbs(jsj),ielg(i)
!              call parallel_abort(errmsg)
!            endif
            sum1=sum1+ssign(j,i)*(zs(max(l+1,kbs(jsj)),jsj)-zs(max(l,kbs(jsj)),jsj))*distj(jsj)*(vnor1+vnor2)/2

            !In eframe
            ubar=ubar+swild98(1,j,l)/3 !su2(l,jsj)/3    
            ubar1=ubar1+swild98(1,j,l+1)/3 !su2(l+1,jsj)/3    
            vbar=vbar+swild98(2,j,l)/3 !sv2(l,jsj)/3    
            vbar1=vbar1+swild98(2,j,l+1)/3 !sv2(l+1,jsj)/3    
          enddo !j=1,3

!         Impose bottom no-flux b.c.
          if(l==kbe(i)) then
            bflux=(av_bdef2-av_bdef1)/dt
            if(imm==2) bflux=bflux0
          else
            bflux=ubar*sne(1,l)+vbar*sne(2,l)+we_fv(l,i)*sne(3,l)
          endif

          we_fv(l+1,i)=(-sum1-(ubar1*sne(1,l+1)+vbar1*sne(2,l+1))*area_e(l+1) + &
     &bflux*area_e(l))/sne(3,l+1)/area_e(l+1)

!#ifdef USE_ICM
!          if(iWQPS==2) then
!            if(l==PSK(i)-1) we_fv(l+1,i)=we_fv(l+1,i)-(PSQ(i)*dt)/area(i)    !added by YC, need check
!          endif
!#endif /*USE_ICM*/

!         Debug
!          tmp1=sum1
!          tmp2=(ubar1*sne(1,l+1)+vbar1*sne(2,l+1)+we(l+1,i)*sne(3,l+1))*area_e(l+1)-bflux*area_e(l)
!          if(i==24044.and.it==2) write(97,*)l,tmp1,tmp2,tmp1+tmp2

        enddo !l=kbe(i),nvrt-1
      enddo !i=1,nea

      deallocate(swild98)
      if(nonhydro==0) we=we_fv

      if(myrank==0) write(16,*)'done solving w'

#ifdef INCLUDE_TIMING
!  end momentum
      wtmp2=mpi_wtime()
      wtimer(8,1)=wtimer(8,1)+wtmp2-wtmp1
!  start transport
      wtmp1=wtmp2
#endif

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Transport equation
!-------------------------------------------------------------------------------
!     Test backtracking alone with rotating Gausshill
      if(ibtrack_test==1) then
        eta1=0; eta2=0; we=0
        rot_per=3000 !period
        rot_f=2*pi/rot_per !freq.
!        xvel0=-1; yvel0=0.9
        do i=1,nsa
          do k=1,nvrt
            su2(k,i)=-ycj(i)*rot_f !xvel0
            sv2(k,i)=xcj(i)*rot_f
          enddo !k
        enddo !i
        do i=1,nea
          do k=1,nvrt
            do j=1,3
              nd=elnode(j,i)
              ufg(k,i,j)=-ynd(nd)*rot_f
              vfg(k,i,j)=xnd(nd)*rot_f
            enddo !j
          enddo !k
        enddo !i
        do i=1,npa
          do k=1,nvrt
            uu2(k,i)=-ynd(i)*rot_f
            vv2(k,i)=xnd(i)*rot_f
            ww2(k,i)=0 !-1.e-4*znl(k,i)*(50+znl(k,i))
          enddo !k
        enddo !i
      endif !ibtrack_test

      if(ibc==0.or.ibtp==1) then
!----------------------------------------------------------------------
!...  Initialize S,T as flags
!      tnd=-99; snd=-99; tsd=-99; ssd=-99 !flags

!*************************************************************************************
!        ELM option
!*************************************************************************************

      if(iupwind_t==0.or.iupwind_s==0) then
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!...  Along each node & side
!...      
      do l=1,2 !nodes&sides
        if(l==1) then
!         Resident nodes only; because radiation b.c. needs side vel. info, and the side may not be inside
          limit=np
        else
!         Aug. sides
          limit=nsa 
        endif
        do i=1,limit

          if(l==1) then; if(idry(i)==1) cycle; endif
          if(l==2) then; if(idry_s(i)==1) cycle; endif

!         Define nodes/sides, and layer thickness & diffusivities
          if(l==1) then !nodes
            nd0=i
            kbb=kbp(i)
            depth=dp(nd0)
          else !sides
            isd0=i
            kbb=kbs(i)
            node1=isidenode(1,isd0)
            node2=isidenode(2,isd0)
            depth=dps(isd0)
          endif
          do k=kbb+1,nvrt
            if(l==1) then !nodes
              dzz(k)=znl(k,nd0)-znl(k-1,nd0)
              dfz(k)=(ptbt(4,k,nd0)+ptbt(4,k-1,nd0))/2
            else !sides
              dzz(k)=zs(k,isd0)-zs(k-1,isd0)
              dfz(k)=(ptbt(4,k,node1)+ptbt(4,k,node2)+ptbt(4,k-1,node1)+ptbt(4,k-1,node2))/4
            endif
          enddo !k

!	  Coefficient matrix: mass lumping
          ndim=nvrt-kbb+1
          do k=kbb,nvrt
            kin=k-kbb+1
            alow(kin)=0
            cupp(kin)=0
            bdia(kin)=0
            if(k<nvrt) then
              tmp=dt*dfz(k+1)/dzz(k+1)
              cupp(kin)=cupp(kin)-tmp
              bdia(kin)=bdia(kin)+dzz(k+1)/2+tmp
            endif

            if(k>kbb) then
              tmp=dt*dfz(k)/dzz(k)
              alow(kin)=alow(kin)-tmp
              bdia(kin)=bdia(kin)+dzz(k)/2+tmp
            endif
          enddo !k

!	  RHS 
!	  b.c. to be imposed at the end
          do k=kbb,nvrt
            kin=k-kbb+1
            rrhs(1,kin)=0
            rrhs(2,kin)=0
            if(k<nvrt) then 
              if(l==1) then
                rrhs(1,kin)=rrhs(1,kin)+dzz(k+1)/2*ptbt(1,k,nd0)
                rrhs(2,kin)=rrhs(2,kin)+dzz(k+1)/2*ptbt(2,k,nd0)
              else
                rrhs(1,kin)=rrhs(1,kin)+dzz(k+1)/2*sdbt(3,k,isd0)
                rrhs(2,kin)=rrhs(2,kin)+dzz(k+1)/2*sdbt(4,k,isd0)
              endif
            else !surface fluxes
              if(ihconsv/=0) then !heat flux
                if(l==1) then
                  rrhs(1,kin)=rrhs(1,kin)+dt/rho0/shw*sflux(nd0)
                else
                  rrhs(1,kin)=rrhs(1,kin)+dt/rho0/shw*(sflux(node1)+sflux(node2))/2
                endif
              endif
              if(isconsv/=0) then !salt flux
                if(l==1) then
                  rrhs(2,kin)=rrhs(2,kin)+dt/rho0*snd(k,nd0)*(fluxevp(nd0)-fluxprc(nd0))
                else
                  rrhs(2,kin)=rrhs(2,kin)+dt/rho0*ssd(k,isd0)* &
     &(fluxevp(node1)+fluxevp(node2)-fluxprc(node1)-fluxprc(node2))/2
                endif
              endif
            endif

            if(k>kbb) then 
              if(l==1) then
                rrhs(1,kin)=rrhs(1,kin)+dzz(k)/2*ptbt(1,k,nd0)
                rrhs(2,kin)=rrhs(2,kin)+dzz(k)/2*ptbt(2,k,nd0)
              else
                rrhs(1,kin)=rrhs(1,kin)+dzz(k)/2*sdbt(3,k,isd0)
                rrhs(2,kin)=rrhs(2,kin)+dzz(k)/2*sdbt(4,k,isd0)
              endif
            endif

!           Compute fucntion F() for solar
!           solar flux= R*exp(z/d_1))+(1-R)*exp(z/d_2) (d_[1,2] are attentuation depths; smaller values for muddier water)
!           The values for R, d_1, d_2 are given below 
!           1: 0.58 0.35 23 (Jerlov type I)
!           2: 0.62 0.60 20 (Jerlov type IA)
!           3: 0.67 1.00 17 (Jerlov type IB)
!           4: 0.77 1.50 14 (Jerlov type II)
!           5: 0.78 1.40 7.9 (Jerlov type III)
!           6: 0.62 1.50 20 (Paulson and Simpson 1977; similar to type IA)
!           7: 0.80 0.90 2.1 (Mike Z.'s choice for estuary)

            if(ihconsv/=0) then !solar
              if(l==1) then
                zz1=max(-5.e2_rkind,znl(k,nd0)-znl(nvrt,nd0))
                rrhs(5,k)=srad(nd0) !F(); more to come
                itmp=iwater_type(nd0)
              else !sides
                zz1=max(-5.e2_rkind,zs(k,isd0)-zs(nvrt,isd0))
                rrhs(5,k)=(srad(node1)+srad(node2))/2
                itmp=max(iwater_type(node1),iwater_type(node2))
              endif
              if(zz1>0) then 
                write(errmsg,*)'Above f.s. (2):',l,i,k,zz1
                call parallel_abort(errmsg)
              endif
!             Water type
              select case(itmp)
                case(1)
                  rr=0.58; d_1=0.35; d_2=23
                case(2)
                  rr=0.62; d_1=0.60; d_2=20
                case(3)
                  rr=0.67; d_1=1.00; d_2=17
                case(4)
                  rr=0.77; d_1=1.50; d_2=14
                case(5)
                  rr=0.78; d_1=1.40; d_2=7.9
                case(6)
                  rr=0.62; d_1=1.50; d_2=20
                case(7)
                  rr=0.80; d_1=0.90; d_2=2.1
                case default
                  call parallel_abort('Unknown water type (2)')
              end select !itmp
              rrhs(5,k)=rrhs(5,k)*(rr*exp(zz1/d_1)+(1-rr)*exp(zz1/d_2)) !F()
              
              ! Zero out bottom
              !if(k==kbb) rrhs(5,k)=0
            endif !solar
          enddo !k=kbb,nvrt

!         Compute solar
          if(ihconsv/=0) then
            rrhs(1,nvrt-kbb+1)=rrhs(1,nvrt-kbb+1)+dt/rho0/shw/2*(rrhs(5,nvrt)-rrhs(5,nvrt-1))
            rrhs(1,1)=rrhs(1,1)+dt/rho0/shw/2*(rrhs(5,kbb+1)-rrhs(5,kbb))
            do k=kbb+1,nvrt-1
              kin=k-kbb+1
              rrhs(1,kin)=rrhs(1,kin)+dt/rho0/shw/2*(rrhs(5,k+1)-rrhs(5,k-1))
            enddo !k
          endif

          call tridag(nvrt,100,ndim,2,alow,bdia,cupp,rrhs,soln,gam)

!         Impose no flux condition at bottom B.L. for slipless bottom
          itmpf=0 !flag
          if(l==1) then; if(Cdp(nd0)/=0) itmpf=1; endif
          if(l==2) then; if(Cd(isd0)/=0) itmpf=1; endif
          if(itmpf==1) soln(1:2,1)=soln(1:2,2)

!         Correct overshoots for S,T
!         Debug
!          if(l==1.and.i==23) then
!            write(98,*)totalflux
!            do k=1,nvrt
!              write(98,*)k,soln(1,k),ptbt(1,k,nd0),rrhs(5,k),rrhs(4,k)*dt
!            enddo
!          endif

          tmin=100; tmax=-100
          smin=100; smax=-100
          do k=kbb,nvrt
            if(l==1) then
              if(tmin>ptbt(1,k,nd0)) tmin=ptbt(1,k,nd0)
              if(tmax<ptbt(1,k,nd0)) tmax=ptbt(1,k,nd0)
              if(smin>ptbt(2,k,nd0)) smin=ptbt(2,k,nd0)
              if(smax<ptbt(2,k,nd0)) smax=ptbt(2,k,nd0)
              rrhs(3,k)=ptbt(1,k,nd0) !for debugging only
            else
              if(tmin>sdbt(3,k,isd0)) tmin=sdbt(3,k,isd0)
              if(tmax<sdbt(3,k,isd0)) tmax=sdbt(3,k,isd0)
              if(smin>sdbt(4,k,isd0)) smin=sdbt(4,k,isd0)
              if(smax<sdbt(4,k,isd0)) smax=sdbt(4,k,isd0)
              rrhs(3,k)=sdbt(3,k,isd0)
            endif
          enddo !k
!         Reset extrema for heat exchange
!          if(ihconsv/=0) then
!            tmin=tempmin
!            tmax=tempmax
!          endif

!          if(tmin>tmax.or.tmin<0.or.smin>smax.or.smin<0) then
!            write(11,*)'Illegal min/max:',tmin,tmax,smin,smax,l,i
!            stop
!          endif

!         Store S,T in swild2 temporarily
          do k=kbb,nvrt
            kin=k-kbb+1
            swild2(k,1:2)=soln(1:2,kin)
            if(ihconsv/=0) swild2(k,1)=max(tempmin,min(tempmax,soln(1,kin)))
            if(isconsv/=0) swild2(k,2)=max(saltmin,min(saltmax,soln(2,kin)))
          enddo !k

!	  Impose b.c. on bnd nodes & sides
          if(l==1) then !nodes
            if(isbnd(1,nd0)>0) then 
              ibnd=isbnd(1,nd0) !impose b.c. according to 1st open bnd; global segment #
              ind=isbnd(-1,nd0) !global index

              isd=0 !flag
!             Two sides must be inside aug. domain as nd0 is resident
              do ll=1,2 !side
                if(ll==1) then
                  ie=indel(1,nd0)
                  id=iself(1,nd0)
                  isd3=elside(nx(id,2),ie)
                else
                  ie=indel(nne(nd0),nd0)
                  id=iself(nne(nd0),nd0)
                  isd3=elside(nx(id,1),ie)
                endif
                if(isbs(isd3)==ibnd) then
                  isd=isd3
                  exit
                endif
              enddo !ll
              if(isd==0) then
                write(errmsg,*)'Cannot find an open side:',nd0
                call parallel_abort(errmsg)
              endif
              do k=1,nvrt
                !Normal vel.
                if(ics==1) then
                  vnn=su2(k,isd)*sframe(1,1,isd)+sv2(k,isd)*sframe(2,1,isd)
                else
                  vnn=su2(k,isd)
                endif !ics
                if(vnn>0) cycle

                !Impose b.c. for inflow only
                if(itetype(ibnd)==1.or.itetype(ibnd)==2) then
                  swild2(k,1)=tobc(ibnd)*tth(1,1,ibnd)+(1-tobc(ibnd))*swild2(k,1)
                else if(itetype(ibnd)==3) then
                  swild2(k,1)=tobc(ibnd)*tem0(k,nd0)+(1-tobc(ibnd))*swild2(k,1)
                else if(itetype(ibnd)==4) then
                  swild2(k,1)=tobc(ibnd)*tth(k,ind,ibnd)+(1-tobc(ibnd))*swild2(k,1)
                endif

                if(isatype(ibnd)==1.or.isatype(ibnd)==2) then
                  swild2(k,2)=sobc(ibnd)*sth(1,1,ibnd)+(1-sobc(ibnd))*swild2(k,2)
                else if(isatype(ibnd)==3) then
                  swild2(k,2)=sobc(ibnd)*sal0(k,nd0)+(1-sobc(ibnd))*swild2(k,2)
                else if(isatype(ibnd)==4) then
                  swild2(k,2)=sobc(ibnd)*sth(k,ind,ibnd)+(1-sobc(ibnd))*swild2(k,2)
                endif
              enddo !k
            endif !isbnd>0
          else !sides
            if(isbs(isd0)>0) then
              ibnd=isbs(isd0)
              nwild(1:2)=0
              do ll=1,2 !nodes
                ndo=isidenode(ll,isd0)
                do lll=1,2 !2 possible bnds
                  if(isbnd(lll,ndo)==ibnd) then
                    nwild(ll)=isbnd(-lll,ndo) !global index
                    exit
                  endif
                enddo !lll
              enddo !ll
              in1=nwild(1); in2=nwild(2);
              if(in1<=0.or.in2<=0) then
                write(errmsg,*)'Impossible 102'
                call parallel_abort(errmsg)
              endif

              do k=1,nvrt
                !Normal vel.
                if(ics==1) then
                  vnn=su2(k,isd0)*sframe(1,1,isd0)+sv2(k,isd0)*sframe(2,1,isd0)
                else
                  vnn=su2(k,isd0)
                endif !ics
                if(vnn>0) cycle

                !Inflow
                if(itetype(ibnd)==1.or.itetype(ibnd)==2) then
                  swild2(k,1)=tobc(ibnd)*tth(1,1,ibnd)+(1-tobc(ibnd))*swild2(k,1)
                else if(itetype(ibnd)==3) then
                  swild2(k,1)=tobc(ibnd)*(tem0(k,node1)+tem0(k,node2))/2+(1-tobc(ibnd))*swild2(k,1)
                else if(itetype(ibnd)==4) then
                  swild2(k,1)=tobc(ibnd)*(tth(k,in1,ibnd)+tth(k,in2,ibnd))/2+(1-tobc(ibnd))*swild2(k,1)
                endif

                if(isatype(ibnd)==1.or.isatype(ibnd)==2) then
                  swild2(k,2)=sobc(ibnd)*sth(1,1,ibnd)+(1-sobc(ibnd))*swild2(k,2)
                else if(isatype(ibnd)==3) then
                  swild2(k,2)=sobc(ibnd)*(sal0(k,node1)+sal0(k,node2))/2+(1-sobc(ibnd))*swild2(k,2)
                else if(isatype(ibnd)==4) then
                  swild2(k,2)=sobc(ibnd)*(sth(k,in1,ibnd)+sth(k,in2,ibnd))/2+(1-sobc(ibnd))*swild2(k,2)
                endif
              enddo !k
            endif !isbs(isd0)>0
          endif

!         Nudging
          if(inu_st/=0) then
            do k=kbb,nvrt
              if(l==1) then !nodes
                if(znl(k,nd0)>=-vnh1) then
                  vnf=vnf1 !vertical nudging factor
                else if(znl(k,nd0)>=-vnh2) then
                  vnf=vnf1+(vnf2-vnf1)*(znl(k,nd0)+vnh1)/(-vnh2+vnh1) 
                else
                  vnf=vnf2
                endif
                tnu=t_nudge(nd0)*vnf*dt
                snu=s_nudge(nd0)*vnf*dt
                if(tnu<0.or.tnu>1.or.snu<0.or.snu>1.or.vnf<0.or.vnf>1) then
                  write(errmsg,*)'Nudging factor out of bound (1):',tnu,snu,vnf
                  call parallel_abort(errmsg)
                endif

                if(inu_st==1) then !to i.c.
                  swild2(k,1)=swild2(k,1)*(1-tnu)+tem0(k,nd0)*tnu
                  swild2(k,2)=swild2(k,2)*(1-snu)+sal0(k,nd0)*snu
                else if(inu_st==2) then
                  swild2(k,1)=swild2(k,1)*(1-tnu)+tnd_nu(k,nd0)*tnu
                  swild2(k,2)=swild2(k,2)*(1-snu)+snd_nu(k,nd0)*snu
                endif
              else !sides
                if(zs(k,isd0)>=-vnh1) then
                  vnf=vnf1 !vertical nudging factor
                else if(zs(k,isd0)>=-vnh2) then
                  vnf=vnf1+(vnf2-vnf1)*(zs(k,isd0)+vnh1)/(-vnh2+vnh1) 
                else
                  vnf=vnf2
                endif
                tnu=(t_nudge(node1)+t_nudge(node2))/2*vnf*dt
                snu=(s_nudge(node1)+s_nudge(node2))/2*vnf*dt
                if(tnu<0.or.tnu>1.or.snu<0.or.snu>1.or.vnf<0.or.vnf>1) then
                  write(errmsg,*)'Nudging factor out of bound (2):',tnu,snu,vnf
                  call parallel_abort(errmsg)
                endif

                if(inu_st==1) then !to i.c.
                  swild2(k,1)=swild2(k,1)*(1-tnu)+(tem0(k,node1)+tem0(k,node2))/2*tnu
                  swild2(k,2)=swild2(k,2)*(1-snu)+(sal0(k,node1)+sal0(k,node2))/2*snu
                else if(inu_st==2) then
                  swild2(k,1)=swild2(k,1)*(1-tnu)+(tnd_nu(k,node1)+tnd_nu(k,node2))/2*tnu
                  swild2(k,2)=swild2(k,2)*(1-snu)+(snd_nu(k,node1)+snd_nu(k,node2))/2*snu
                endif
              endif
            enddo !k
          endif !nudging

!         Extend
          do k=1,kbb-1
            swild2(k,1)=swild2(kbb,1)
            swild2(k,2)=swild2(kbb,2)
          enddo !k

!         Check bounds
          do k=1,nvrt
            if(swild2(k,1)<-98.or.swild2(k,2)<-98) then
              write(errmsg,*)'Werid ST (1):',l,i,k,swild2(k,1),swild2(k,2)
              call parallel_abort(errmsg)
            endif
          enddo !k

!         Output S,T
          do k=1,nvrt
            if(iupwind_t==0) then
              if(l==1) then
                tnd(k,nd0)=swild2(k,1)
              else
                tsd(k,isd0)=swild2(k,1)
              endif
            endif !iupwind_t

            if(iupwind_s==0) then
              if(l==1) then
                snd(k,nd0)=swild2(k,2)
              else
                ssd(k,isd0)=swild2(k,2)
              endif
            endif !iupwind_s
          enddo !k

        enddo !i=1,limit
      enddo !l=1,2

!     Exchange ghost nodes
#ifdef INCLUDE_TIMING
      cwtmp=mpi_wtime()
#endif
      if(iupwind_t==0) call exchange_p3dw(tnd)
      if(iupwind_s==0) call exchange_p3dw(snd)
#ifdef INCLUDE_TIMING
      wtimer(9,2)=wtimer(9,2)+mpi_wtime()-cwtmp
#endif

!     Compute tsel for tracer transport (temporary fix)
      do i=1,nea
        if(idry_e(i)==1) cycle

        do k=kbe(i)+1,nvrt
          tsel(1,k,i)=sum(tnd(k,elnode(1:3,i))+tnd(k-1,elnode(1:3,i)))/6
          tsel(2,k,i)=sum(snd(k,elnode(1:3,i))+snd(k-1,elnode(1:3,i)))/6
        enddo !k    
      enddo !i
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      endif
!     end of ELM option

!*************************************************************************************
!
!        Upwind and TVD option
!
!*************************************************************************************
      if(iupwind_t/=0) then
!       b.c. and body forces
        bdy_frc=0; flx_sf=0; flx_bt=0

!       Salt exchange
        if(isconsv/=0) then
          do i=1,nea
            if(idry_e(i)==1) cycle
       
            n1=elnode(1,i)
            n2=elnode(2,i)
            n3=elnode(3,i)
            evap=(fluxevp(n1)+fluxevp(n2)+fluxevp(n3))/3
            precip=(fluxprc(n1)+fluxprc(n2)+fluxprc(n3))/3
            flx_sf(2,i)=tsel(2,nvrt,i)*(evap-precip)/rho0
          enddo !i
        endif !isconsv/=0

!       Heat exchange
        if(ihconsv/=0) then
          do i=1,nea
            if(idry_e(i)==1) cycle

!           Wet element
            n1=elnode(1,i)
            n2=elnode(2,i)
            n3=elnode(3,i)

!           Surface flux
            sflux_e=(sflux(n1)+sflux(n2)+sflux(n3))/3
            flx_sf(1,i)=sflux_e/rho0/shw

!           Solar
!           Calculate water type
!           solar flux= R*exp(z/d_1))+(1-R)*exp(z/d_2) (d_[1,2] are attentuation depths; smaller values for muddier water)
!           The values for R, d_1, d_2 are given below
!           1: 0.58 0.35 23 (Jerlov type I)
!           2: 0.62 0.60 20 (Jerlov type IA)
!           3: 0.67 1.00 17 (Jerlov type IB)
!           4: 0.77 1.50 14 (Jerlov type II)
!           5: 0.78 1.40 7.9 (Jerlov type III)
!           6: 0.62 1.50 20 (Paulson and Simpson 1977; similar to type IA)
!           7: 0.80 0.90 2.1 (Mike Z.'s choice for estuary)
            itmp=max(iwater_type(n1),iwater_type(n2),iwater_type(n3))
            select case(itmp)
              case(1)
                rr=0.58; d_1=0.35; d_2=23
              case(2)
                rr=0.62; d_1=0.60; d_2=20
              case(3)
                rr=0.67; d_1=1.00; d_2=17
              case(4)
                rr=0.77; d_1=1.50; d_2=14
              case(5)
                rr=0.78; d_1=1.40; d_2=7.9
              case(6)
                rr=0.62; d_1=1.50; d_2=20
              case(7)
                rr=0.80; d_1=0.90; d_2=2.1
              case default
                call parallel_abort('Unknown water type (3)')
            end select !itmp
            
            srad_e(i)=(srad(n1)+srad(n2)+srad(n3))/3
            do k=kbe(i)+1,nvrt
!             Don't use eta2 as it has been updated but not znl()
              dp1=min(ze(nvrt,i)-ze(k-1,i),500._rkind) !to prevent underflow
              dp2=min(ze(nvrt,i)-ze(k,i),500._rkind) !to prevent underflow
              if(dp2<0.or.dp2>dp1) then
                write(errmsg,*)'Depth<0 in upwind transport:',i,k,dp1,dp2, &
     &ze(nvrt,i),(l,znl(l,elnode(1:3,i)),l=kbe(i),nvrt)
                call parallel_abort(errmsg)
              endif

!              if(k==kbe(i)+1) then
!                srad1=0
!              else
!              endif
              srad1=srad_e(i)*(rr*exp(-dp1/d_1)+(1-rr)*exp(-dp1/d_2))
              srad2=srad_e(i)*(rr*exp(-dp2/d_1)+(1-rr)*exp(-dp2/d_2))
              if(srad2<srad1.and.ifort12(19)==0) then
                ifort12(19)=1
                write(12,*)'Reset negative solar hearting:',ielg(i),k,srad2,srad1,srad2-srad1
              endif
              bdy_frc(1,k,i)=max(srad2-srad1,0._rkind)/rho0/shw/(ze(k,i)-ze(k-1,i)) !Q
            enddo !k=kbe(i)+1,nvrt
          enddo !i=1,nea
        endif !heat exchange

!       Point sources/sinks; at bottom layer
!Error: need to reconcile with ICM
!        if(if_source==1.and.mass_source==1) then
!          do i=1,nea
!            if(idry_e(i)==1.or.vsource(i)<=0) cycle
!
!            !Positive source only
!            bigv=area(i)*(ze(kbe(i)+1,i)-ze(kbe(i),i))
!            if(bigv==0) call parallel_abort('STEP: bigv==0')
!            do j=1,2 !T,S
!              bdy_frc(j,kbe(i)+1,i)=bdy_frc(j,kbe(i)+1,i)+ &
!     &(msource(j,i)-tsel(j,kbe(i)+1,i))*vsource(i)/bigv
!            enddo !j
!          enddo !i
!        endif !if_source

!       VIMS surface temperature; added by YC
#ifdef USE_ICM
        if(myrank==0) write(16,*)'start ICM adjust. surface T..'
        if(iSun==2) then
          do i=1,nea
            if(idry_e(i)==1) cycle
            n1=elnode(1,i)
            n2=elnode(2,i)
            n3=elnode(3,i)
            sflux_e=(surf_t(n1)+surf_t(n2)+surf_t(n3))/3
            tsel(1,nvrt,i)=sflux_e
          enddo !i
        if(myrank==0) write(16,*)'end ICM adjust. surface T..'
        endif !iSun=2
#endif /*USE_ICM*/

        up_tvd=iupwind_t>=2
        tr_el(1:2,:,:)=tsel(1:2,:,:)
        if(iupwind_t<=2) then !upwind or explicit TVD
          call do_transport_tvd(it,0,up_tvd,tvd_mid1,flimiter1,2,difnum_max_l,nvrt,npa,ptbt(4,:,:))
        else if(iupwind_t==3) then !implicit TVD
          call do_transport_tvd_imp(it,0,up_tvd,tvd_mid1,flimiter1,2,difnum_max_l,nvrt,npa,ptbt(4,:,:))
        endif !iupwind_t
        tsel(1:2,:,:)=tr_el(1:2,:,:)
        if(difnum_max_l>difnum_max_l2) difnum_max_l2=difnum_max_l

!       Debug
!        fdb='tsel_0000'
!        lfdb=len_trim(fdb)
!        write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
!        open(32,file=trim(fdb),status='unknown')        
!        do i=1,nea
!          do k=1,nvrt
!            write(32,*)i,k,tsel(1:2,k,i)
!          enddo !k
!        enddo !i
!        close(32)
!        call parallel_finalize
!        stop

!       Point sources/sinks using operator splitting (that guarentees max.
!       principle); at bottom layer
!       Do nothing for net sinks
!Error: need to reconcile with ICM
        if(if_source==1) then
          do i=1,nea
            if(idry_e(i)==1.or.vsource(i)<=0) cycle

            !Positive source only
            bigv=area(i)*(ze(kbe(i)+1,i)-ze(kbe(i),i))
            if(bigv<=0) call parallel_abort('STEP: bigv==0')
            rat=vsource(i)*dt/bigv !ratio of volumes (>0)
            do j=1,2 !T,S
              tsel(j,kbe(i)+1,i)=(tsel(j,kbe(i)+1,i)+rat*msource(j,i))/(1+rat)
            enddo !j
          enddo !i
        endif !if_source

#ifdef USE_ICM
      !ICM's method for imposing point source
      if(iWQPS==2) then 
        if(myrank==0) write(16,*)'impose ICM point source S..'
        !Impose point source bc
        do i=1,nea
          if(idry_e(i)==1) cycle
          if(PSQ(i).ne.0) then
            PSK(i)=max(kbe(i)+1,min(nvrt,PSK(i)))
            bigv=area(i)*(ze(PSK(i),i)-ze(PSK(i)-1,i)) !volume
            !YC bdy_frc(2,PSK(i),i)=area(i)*(tsel(2,PSK(i),i)*(PSQ(i))/area(i))/bigv
            ! tsel(2,PSK(i),i)=(tsel(2,PSK(i),i)*bigv)/(bigv-PSQ(i)*dt) !ZG
          endif
        enddo !i
        if(myrank==0) write(16,*)'end impose ICM point source S..'
      endif !iWQPS
#endif /*USE_ICM*/

!       Nudging
        do i=1,nea
          if(idry_e(i)==1) cycle

          n1=elnode(1,i)
          n2=elnode(2,i)
          n3=elnode(3,i)

          if(inu_st/=0) then
            do k=kbe(i)+1,nvrt
              if(ze(k,i)>=-vnh1) then
                vnf=vnf1 !vertical nudging factor
              else if(ze(k,i)>=-vnh2) then
                vnf=vnf1+(vnf2-vnf1)*(ze(k,i)+vnh1)/(-vnh2+vnh1)
              else
                vnf=vnf2
              endif
              tnu=(t_nudge(n1)+t_nudge(n2)+t_nudge(n3))/3*vnf*dt
              snu=(s_nudge(n1)+s_nudge(n2)+s_nudge(n3))/3*vnf*dt
              if(tnu<0.or.tnu>1.or.snu<0.or.snu>1.or.vnf<0.or.vnf>1) then
                write(errmsg,*)'Nudging factor out of bound (1):',tnu,snu,vnf
                call parallel_abort(errmsg)
              endif

              if(inu_st==1) then !to i.c.
                tsel01=(tem0(k,n1)+tem0(k,n2)+tem0(k,n3)+tem0(k-1,n1)+tem0(k-1,n2)+tem0(k-1,n3))/6
                tsel02=(sal0(k,n1)+sal0(k,n2)+sal0(k,n3)+sal0(k-1,n1)+sal0(k-1,n2)+sal0(k-1,n3))/6
                tsel(1,k,i)=tsel(1,k,i)*(1-tnu)+tsel01*tnu
                tsel(2,k,i)=tsel(2,k,i)*(1-snu)+tsel02*snu
              else if(inu_st==2) then
                tnd_nu_e=(tnd_nu(k,n1)+tnd_nu(k,n2)+tnd_nu(k,n3)+tnd_nu(k-1,n1)+tnd_nu(k-1,n2)+tnd_nu(k-1,n3))/6
                snd_nu_e=(snd_nu(k,n1)+snd_nu(k,n2)+snd_nu(k,n3)+snd_nu(k-1,n1)+snd_nu(k-1,n2)+snd_nu(k-1,n3))/6
                tsel(1,k,i)=tsel(1,k,i)*(1-tnu)+tnd_nu_e*tnu
                tsel(2,k,i)=tsel(2,k,i)*(1-snu)+snd_nu_e*snu
              endif
            enddo !k
          endif !inu_st/=0

!         Extend
          do k=1,kbe(i)
            tsel(1:2,k,i)=tsel(1:2,kbe(i)+1,i)
          enddo !k
        enddo !i=1,nea

!       Convert to S,T at nodes and sides and whole levels
!       Use hp_int to temporarily store values at elements and whole levels
        do i=1,nea
          if(idry_e(i)==1) cycle

          do k=kbe(i)+1,nvrt-1
            zrat=(ze(k+1,i)-ze(k,i))/(ze(k+1,i)-ze(k-1,i))
            if(zrat<=0.or.zrat>=1) then
              write(errmsg,*)'Ratio out of bound:',i,k,zrat
              call parallel_abort(errmsg)
            endif
            hp_int(k,i,1:2)=(1-zrat)*tsel(1:2,k+1,i)+zrat*tsel(1:2,k,i)
          enddo !k
          hp_int(nvrt,i,1:2)=tsel(1:2,nvrt,i)
          hp_int(kbe(i),i,1:2)=tsel(1:2,kbe(i)+1,i)
        enddo !i=1,nea

        do i=1,np
          if(idry(i)==1) cycle

          do k=1,nvrt
            tt1=0; ss1=0
            ta=0
            do j=1,nne(i)
              ie=indel(j,i)
              if(idry_e(ie)==0) then
                ta=ta+area(ie)
                kin=max0(k,kbe(ie))
                tt1=tt1+hp_int(kin,ie,1)*area(ie)
                ss1=ss1+hp_int(kin,ie,2)*area(ie)
              endif
            enddo !j
            if(ta==0) then !from levels(), a node is wet if and only if at least one surrounding element is wet
              write(errmsg,*)'Isolated wet node (9):',i
              call parallel_abort(errmsg)
            else
              if(iupwind_t/=0) tnd(k,i)=tt1/ta
              if(iupwind_s/=0) snd(k,i)=ss1/ta
            endif
          enddo !k
        enddo !i=1,np

        do i=1,ns
          if(idry_s(i)==1) cycle

          do k=1,nvrt
            tt1=0; ss1=0
            ta=0
            do j=1,2
              ie=isdel(j,i)
              if(ie/=0.and.idry_e(max(1,ie))==0) then
                ta=ta+area(ie)
                kin=max0(k,kbe(ie))
                tt1=tt1+hp_int(kin,ie,1)*area(ie)
                ss1=ss1+hp_int(kin,ie,2)*area(ie)
              endif
            enddo !j
            if(ta==0) then 
              write(errmsg,*)'Isolated wet side (9):',i,(isdel(j,i),j=1,2)
              call parallel_abort(errmsg)
            else
              if(iupwind_t/=0) tsd(k,i)=tt1/ta
              if(iupwind_s/=0) ssd(k,i)=ss1/ta
            endif
          enddo !k
        enddo !i=1,ns

!       Exchange
        if(iupwind_t/=0.and.iupwind_s/=0) then
!         More efficient exchange
          allocate(swild98(2,nvrt,nsa),stat=istat)
          if(istat/=0) call parallel_abort('MAIN: fail to allocate swild98')
!'
          swild98(1,:,1:npa)=tnd(:,:)
          swild98(2,:,1:npa)=snd(:,:)
#ifdef INCLUDE_TIMING
          cwtmp=mpi_wtime()
#endif
          call exchange_p3d_2(swild98)
#ifdef INCLUDE_TIMING
          wtimer(9,2)=wtimer(9,2)+mpi_wtime()-cwtmp
#endif
          tnd(:,:)=swild98(1,:,1:npa)
          snd(:,:)=swild98(2,:,1:npa)

          swild98(1,:,:)=tsd(:,:)
          swild98(2,:,:)=ssd(:,:)
#ifdef INCLUDE_TIMING
          cwtmp=mpi_wtime()
#endif
          call exchange_s3d_2(swild98)
#ifdef INCLUDE_TIMING
          wtimer(9,2)=wtimer(9,2)+mpi_wtime()-cwtmp
#endif
          tsd(:,:)=swild98(1,:,:)
          ssd(:,:)=swild98(2,:,:)
          deallocate(swild98)
        else !one of them uses ELM
#ifdef INCLUDE_TIMING
          cwtmp=mpi_wtime()
#endif
          if(iupwind_t/=0) then
            call exchange_p3dw(tnd)
            call exchange_s3dw(tsd)
          endif
          if(iupwind_s/=0) then
            call exchange_p3dw(snd)
            call exchange_s3dw(ssd)
          endif
#ifdef INCLUDE_TIMING
          wtimer(9,2)=wtimer(9,2)+mpi_wtime()-cwtmp
#endif
        endif

      endif !upwind
!----------------------------------------------------------------------
      endif !ibc.eq.0.or.ibtp.eq.1

!...  Tracer transport
      if(ntracers>0) then
!       Save previous results in tr_el (global array) to prepare for do_transport_* 
        tr_el(1:ntracers,:,:)=trel(1:ntracers,:,:)

        select case(flag_model)
          case(-1) !generic; modify arrays below
!           user-defined tracer part
!           define bdy_frc, flx_sf, flx_bt
!           bdy_frc(mntr,kbe(i)+1:nvrt,nea): body force at prism center Q_{i,k} (for all wet elements i);
!                                            has a dimension of [C]/s, where [C] is dimension of the tracer
!           flx_sf(mntr,nea): surface b.c. \kappa*dC/dz = flx_sf (at element center)
!           flx_bt(mntr,nea): bottom b.c.
            do i=1,nea
              if(idry_e(i)==1) cycle

              !Element wet
              do j=1,ntracers
                flx_sf(j,i)=0
                flx_bt(j,i)=0
                do k=kbe(i)+1,nvrt !all prisms along vertical
                  bdy_frc(j,k,i)= 0
                enddo !k
              enddo !j
            enddo !i
!           end user-defined tracer part

          case(0) !tracer age
            flx_bt=0
            flx_sf=0

            do i=1,nea
              if(idry_e(i)==1) cycle

              !Element wet
              do j=1,ntracers
                do k=kbe(i)+1,nvrt !all prisms along vertical
                  if(j<=ntracers/2) then
                    bdy_frc(j,k,i)=0
                  else
                    bdy_frc(j,k,i)=tr_el(j-ntracers/2,k,i)
                  endif
                enddo !k
              enddo !j
            enddo !i

            !Test pt source (imposed at bottom layer): use operator splitting as in
            !T,S
!            if(if_source==1.and.mass_source==1) then
!              do i=1,nea
!                if(idry_e(i)==1.or.vsource(i)<=0) cycle
!
!                !Positive source only
!                bigv=area(i)*(ze(kbe(i)+1,i)-ze(kbe(i),i))
!                if(bigv==0) call parallel_abort('STEP: bigv==0 (2)')
!                do j=1,ntracers
!                  bdy_frc(j,kbe(i)+1,i)=(msource(j+2,i)-trel(j,kbe(i)+1,i))*vsource(i)/bigv
!                enddo !j
!              enddo !i
!            endif !if_source

          case(1) !Sediment
#ifdef USE_SED
            if(myrank==0) write(16,*) 'Entering sediment model...'

!LLP
!----------------------------------------------------------------------
! Compute element depth averaged hvel for VRIJN bedload
!----------------------------------------------------------------------

            dav=0
            do i=1,npa
              if(idry(i)==1) cycle
              do k=kbp(i),nvrt-1
                dav(1,i)=dav(1,i)+(uu2(k+1,i)+uu2(k,i))/2*(znl(k+1,i)-znl(k,i))
                dav(2,i)=dav(2,i)+(vv2(k+1,i)+vv2(k,i))/2*(znl(k+1,i)-znl(k,i))
              enddo !k
              htot=eta2(i)+dp(i)
              if(htot<=h0) then
!                write(errmsg,*)'Impossible 24:',it,i,eta2(i),dp(i),htot,h0,iplg(i)
!                call parallel_abort(errmsg)
                !This is possible because level indices have not been updated
                dav(1:2,i)=0
              else
                dav(1:2,i)=dav(1:2,i)/htot
              endif
            enddo !i=1,npa

            dave=0.d0
            do i=1,nea
              if (idry_e(i)==1) cycle
              cff1=0.d0
              cff2=0.d0
              n1=elnode(1,i)
              n2=elnode(2,i)
              n3=elnode(3,i)
              cff1=(dav(1,n1)+dav(1,n2)+dav(1,n3))/3
              cff2=(dav(2,n1)+dav(2,n2)+dav(2,n3))/3
              dave(i)=sqrt(cff1*cff1+cff2*cff2)
            enddo !i

            bdy_frc = 0.d0
            flx_bt = 0.d0
            flx_sf = 0.d0

            call sediment(it,moitn0,mxitn0,rtol0,dave)
            !bdy_frc updated by the routine above; flx_* are not

#endif /*USE_SED*/
! LLP end
            if(myrank==0) write(16,*) 'done sediment model...'

          case(2) !EcoSim
#ifdef USE_ECO 
!...        Calculates spectral irradiance
!...        Gets hour and yday (day of th year)
            yday = yday + dt/86400
            hour = hour + dt/3600
            if (hour==24) hour = 0
            if (yday==366) yday = 1

            if(myrank==0) write(16,*) 'Calculating spectral irradiance (0)'
!'
            call spec_ir(Tair, Pair, Hair, cloud, Uwind, Vwind) !,SpecIr, avcos)

!...        Calculates sources and sinks terms of the ecological model
            if(myrank==0) write(16,*) 'Calculating ecological sources and sinks terms (0)'
!'
            bdy_frc = 0
            flx_bt = 0
            flx_sf = 0

            call ecosim(Uwind,Vwind)
            if(myrank==0) write(16,*) 'Done ecological sources and sinks terms (0)'
!'
#endif
          case(3) !Oil spill
          case(4) !NAPZD model: Spitz
#ifdef USE_NAPZD
            bdy_frc = 0.d0
            flx_bt = 0.d0
            flx_sf = 0.d0
            if(myrank==0) write(16,*) 'entering NAPZD model....'
            call napzd_spitz(nea,npa,nvrt,ntracers,ntracers2,srad_e)
!Debug
!            call parallel_barrier
            if(myrank==0) write(16,*) 'done NAPZD preparation....'
#endif /*USE_NAPZD*/
          case(5) !ICM
!            tr_el(1:ntracers,:,:)=trel(1:ntracers,:,:) !already init.
            bdy_frc = 0.d0                                               !added by YC
            flx_bt = 0.d0                                                !added by YC
            flx_sf = 0.d0
          case(6) !TIMOR
#ifdef USE_TIMOR
            !Treat settling vel. inside routine
            flx_bt=0
            flx_sf=0
            bdy_frc=0

            !write(12,*)'B4 do_trans'
            !do i=1,nea
            !  write(12,*)i,trel(1,:,i)
            !enddo !
#endif /*USE_TIMOR*/
          case default
            call parallel_abort('Unknown tracer model (9)')
        end select !flag_model

        up_tvd=itr_met>=2
        if(itr_met<=2) then !upwind or explicit TVD
          call do_transport_tvd(it,1,up_tvd,tvd_mid2,flimiter2,ntracers,difnum_max_l,nvrt,npa,ptbt(4,:,:))
        else if(itr_met==3) then !implicit TVD
          call do_transport_tvd_imp(it,1,up_tvd,tvd_mid2,flimiter2,ntracers,difnum_max_l,nvrt,npa,ptbt(4,:,:))
        endif !itr_met
        if(myrank==0) write(16,*)'done tracer transport...'

!#ifdef USE_TIMOR
        if(irouse_test==1) then
          tr_el(:,1:2,:)=1 
        endif
!#endif /*USE_TIMOR*/


        !Debug
        !do j=1,ntracers
        !  write(12,*)'After trc. trans.:',it,j,real(tr_el(j,:,8))
        !enddo !j
 

#ifdef USE_ICM
        if(myrank==0) write(16,*)'start ICM (5)..'
        if(amod(real(time),86400.)==0.) then
          call WQinput !(time)
          call wqm_out
        endif
        if(myrank==0) write(16,*) 'Calculating ecological sources and sinks terms...'     !added by YC
        !call ecosystem(iths_main,it,rnday,WSRP,WSLP,WSPB1,WSPB2,WSPB3,turb,WRea,PSQ)           !added by YC
        !Main routine of ICM
        call ecosystem(it)
        if(myrank==0) write(16,*) 'Done ecological sources and sinks terms...'            !added by YC
#endif /*USE_ICM*/

        trel(1:ntracers,:,:)=tr_el(1:ntracers,:,:)
        if(difnum_max_l>difnum_max_l2) difnum_max_l2=difnum_max_l

!       Point sources/sinks using operator splitting (that guarentees max.
!       principle); at bottom layer
!       Do nothing for net sinks
!Error: need to reconcile with ICM
        if(if_source==1) then
          do i=1,nea
            if(idry_e(i)==1.or.vsource(i)<=0) cycle

            !Positive source only
            bigv=area(i)*(ze(kbe(i)+1,i)-ze(kbe(i),i))
            if(bigv<=0) call parallel_abort('STEP: bigv==0 (3)')
            rat=vsource(i)*dt/bigv !ratio of volumes (>0)
            do j=1,ntracers
              trel(j,kbe(i)+1,i)=(trel(j,kbe(i)+1,i)+rat*msource(j+2,i))/(1+rat)
            enddo !j
          enddo !i
        endif !if_source

#ifdef USE_ICM
      if(iWQPS==2) then
        if(myrank==0) write(16,*)'adjust ICM pt source (trcr) ..'
        ! Impose point source bc added by YC
        do i=1,nea
          if(idry_e(i)==1) cycle
          if(PSQ(i).ne.0.) then
            bigv=area(i)*(ze(PSK(i),i)-ze(PSK(i)-1,i)) !volume
            do n=8,ntracers
              if(n.eq.8) total_loading=WWPRPOC(i)
              if(n.eq.9) total_loading=WWPLPOC(i)
              if(n.eq.10) total_loading=WWPDOCA(i)
              if(n.eq.11) total_loading=WWPRPON(i)
              if(n.eq.12) total_loading=WWPLPON(i)
              if(n.eq.13) total_loading=WWPDON(i)
              if(n.eq.14) total_loading=WWPNH4(i)
              if(n.eq.15) total_loading=WWPNO3(i)
              if(n.eq.16) total_loading=WWPRPOP(i)
              if(n.eq.17) total_loading=WWPLPOP(i)
              if(n.eq.18) total_loading=WWPDOP(i)
              if(n.eq.19) total_loading=WWPPO4t(i)
              if(n.eq.20) total_loading=WWPSU(i)
              if(n.eq.21) total_loading=WWPSAt(i)
              if(n.eq.22) total_loading=WWPCOD(i)
              if(n.eq.23) total_loading=WWPDO(i)
              ! trel(n,PSK(i),i)=(trel(n,PSK(i),i)*(bigv+PSQ(i)*dt)+(total_loading*dt/86400.))/bigv
            enddo
          endif
        enddo !i
        if(myrank==0) write(16,*)'done adjust ICM pt source..'
      endif
#endif /*USE_ICM*/

!       Nudging
        do i=1,nea
          if(idry_e(i)==1) cycle

          n1=elnode(1,i)
          n2=elnode(2,i)
          n3=elnode(3,i)

          if(inu_tr/=0) then
            do k=kbe(i)+1,nvrt
              !Horizontal relax. only
              trnu=(tr_nudge(n1)+tr_nudge(n2)+tr_nudge(n3))/3*dt
              if(trnu<0.or.trnu>1) then
                write(errmsg,*)'Nudging factor out of bound (2):',trnu
                call parallel_abort(errmsg)
              endif

              if(inu_tr==1) then !to i.c.
                trel(1:ntracers,k,i)=trel(1:ntracers,k,i)*(1-trnu)+trel0(1:ntracers,k,i)*trnu
              else if(inu_tr==2) then
                do j=1,ntracers
                  tmp=(trnd_nu(j,k,n1)+trnd_nu(j,k,n2)+trnd_nu(j,k,n3))/3
                  trel(j,k,i)=trel(j,k,i)*(1-trnu)+tmp*trnu
                enddo !j
              endif !inu_tr
            enddo !k
          endif !inu_tr/=0

!         Extend
          do k=1,kbe(i)
            trel(1:ntracers,k,i)=trel(1:ntracers,kbe(i)+1,i)
          enddo !k
        enddo !i=1,nea

!Debug
!        write(12,*)'stage 1'

!       Convert to nodes and whole levels
!       Use tr_el to temporarily store values at elements and whole levels
!        if(mod(it,nspool)==0) then
        do i=1,nea
          if(idry_e(i)==1) cycle

          do k=kbe(i)+1,nvrt-1
            zrat=(ze(k+1,i)-ze(k,i))/(ze(k+1,i)-ze(k-1,i))
            if(zrat<=0.or.zrat>=1) then
              write(errmsg,*)'Ratio out of bound (2):',i,k,zrat
              call parallel_abort(errmsg)
            endif
            tr_el(1:ntracers,k,i)=(1-zrat)*trel(1:ntracers,k+1,i)+zrat*trel(1:ntracers,k,i)
          enddo !k
          tr_el(1:ntracers,nvrt,i)=trel(1:ntracers,nvrt,i)
          tr_el(1:ntracers,kbe(i),i)=trel(1:ntracers,kbe(i)+1,i)

!LLP
#ifdef USE_SED
          !tr_el at surface and bottom. The total mass at centers equal to total mass at levels.
          !tr_tc - tracer vertical total mass at centers
          !tr_tl - tracer vertical total mass at levels

          tr_tc=0.d0
          tr_tl=0.d0

          do k=kbe(i)+1,nvrt
            vol=(ze(k,i)-ze(k-1,i))*area(i)
            tr_tc(1:ntracers,i)=tr_tc(1:ntracers,i)+vol*trel(1:ntracers,k,i)
          enddo !k
          do k=kbe(i)+1,nvrt-1
            vol=(ze(k+1,i)-ze(k-1,i))/2*area(i)
            tr_tl(1:ntracers,i)=tr_tl(1:ntracers,i)+vol*tr_el(1:ntracers,k,i)
          enddo !k

          n1=elnode(1,i)
          n2=elnode(2,i)
          n3=elnode(3,i)
!!...     difussivity of surface level (nvrt)
          av_df=(dfh(nvrt-1,n1)+dfh(nvrt-1,n2)+dfh(nvrt-1,n3))/3
          swild(1:ntracers)=av_df+Wsed(1:ntracers)*(ze(nvrt,i)-ze(nvrt-1,i))
          do j=1,ntracers
            if(swild(j)==0) call parallel_abort('MAIN: sed. div. by 0 (1)')
!'
          enddo !j
          tr_el(1:ntracers,nvrt,i)=(av_df*tr_el(1:ntracers,nvrt-1,i))/swild(1:ntracers)
!!... surface
          vol=((ze(nvrt,i)-ze(nvrt-1,i))/2)*area(i)
!!... bottom
          vol1=((ze(kbe(i)+1,i)-ze(kbe(i),i))/2)*area(i)
          tr_tl(1:ntracers,i)=tr_tl(1:ntracers,i)+vol*tr_el(1:ntracers,nvrt,i)
!          if(myrank==0)write(16,*)'vol',vol,tr_tl(1,i)
          if(vol1==0) call parallel_abort('MAIN: sed. div. by 0 (2)')
          tr_el(1:ntracers,kbe(i),i)=(tr_tc(1:ntracers,i)-tr_tl(1:ntracers,i))/vol1
!          if(myrank==0)write(16,*)'vol1',vol1,tr_tc(1,i),tr_tl(1,i)
#else
          tr_el(1:ntracers,nvrt,i)=trel(1:ntracers,nvrt,i)
          tr_el(1:ntracers,kbe(i),i)=trel(1:ntracers,kbe(i)+1,i)
#endif /*USE_SED*/
        enddo !i=1,nea

!       For rewetted nodes, use value at last wet step
!        tr_nd=-99 !for dry nodes
        do i=1,np
          if(idry(i)==1) cycle

          do k=1,nvrt
            swild(1:ntracers)=0
            ta=0
            do j=1,nne(i)
              ie=indel(j,i)
              if(idry_e(ie)==0) then
                ta=ta+area(ie)
                kin=max0(k,kbe(ie))
                swild(1:ntracers)=swild(1:ntracers)+tr_el(1:ntracers,kin,ie)*area(ie)
              endif
            enddo !j
            if(ta==0) then !from levels(), a node is wet if and only if at least one surrounding element is wet
              write(errmsg,*)'Isolated wet node (9):',i
              call parallel_abort(errmsg)
            else
              tr_nd(1:ntracers,k,i)=swild(1:ntracers)/ta
            endif
          enddo !k
        enddo !i=1,np

!Debug
!        write(12,*)'stage 2'

#ifdef USE_NAPZD
!CSD reuse tr_el array now to hold the Bio_bdef field interpolated to the cell faces.
        do i=1,nea
          if(idry_e(i)==1) cycle
          do k=kbe(i)+1,nvrt-1
            zrat=(ze(k+1,i)-ze(k,i))/(ze(k+1,i)-ze(k-1,i))
            if(zrat<=0.or.zrat>=1) then
              write(errmsg,*)'Ratio out of bound (2):',i,k,zrat
              call parallel_abort(errmsg)
            endif
            tr_el(1,k,i)=(1-zrat)*Bio_bdef(k+1,i)+zrat*Bio_bdef(k,i)
          enddo !k
          tr_el(1,nvrt,i)=Bio_bdef(nvrt,i)
          tr_el(1,kbe(i),i)=Bio_bdef(kbe(i)+1,i)
        enddo !i=1,nea

        Bio_bdefp=-99 !for dry nodes
        do i=1,np
          if(idry(i)==1) cycle

          do k=1,nvrt
            swild(1)=0
            ta=0
            do j=1,nne(i)
              ie=indel(j,i)
              if(idry_e(ie)==0) then
                ta=ta+area(ie)
                kin=max0(k,kbe(ie))
                swild(1)=swild(1)+tr_el(1,kin,ie)*area(ie)
              endif
            enddo !j
            if(ta==0) then !from levels(), a node is wet iff at least one surrounding element is wet
              write(errmsg,*)'Isolated wet node (9):',i
              call parallel_abort(errmsg)
            else
              Bio_bdefp(k,i)=swild(1)/ta
            endif
          enddo !k
        enddo !i=1,np
!Debug
!        write(12,*)'stage 3'

#endif /*USE_NAPZD*/

#ifdef INCLUDE_TIMING
          cwtmp=mpi_wtime()
#endif
          call exchange_p3d_tr(tr_nd)
#ifdef INCLUDE_TIMING
          wtimer(9,2)=wtimer(9,2)+mpi_wtime()-cwtmp
#endif
!        endif !mod(it,nspool)==0

      endif !ntracers>0
!...  End of tracer transport

      if(myrank==0) write(16,*)'done solving transport equation'

#ifdef USE_SED2D
#ifdef INCLUDE_TIMING
      cwtmp2=mpi_wtime() !start of timer
#endif

      if (nb_class==1) then
        call sed2d_main(it)
      elseif (nb_class>1) then
        call sed2d_main_mcml(it)
      endif

#ifdef INCLUDE_TIMING
      timer_ns(3)=timer_ns(3)+mpi_wtime()-cwtmp2 !end timing this section
#endif 
#endif /*USE_SED2D*/

#ifdef INCLUDE_TIMING
! end transport
      wtmp2=mpi_wtime()
      wtimer(9,1)=wtimer(9,1)+wtmp2-wtmp1
! start computing levels
      wtmp1=wtmp2
#endif

!...  Update bed deformation and depth info
      do i=1,npa
        bdef1(i)=bdef2(i)
        if(imm==1) then
          dp(i)=dp00(i)-bdef1(i)
        else if(imm==2) then
          call update_bdef(time,xnd(i),ynd(i),dep,swild)
          dp(i)=dep !min(1.,7-(xnd(i)+time))
        endif
        if(ivcor==2) hmod(i)=min(dp(i),h_s)
      enddo !i

      do i=1,nsa
        n1=isidenode(1,i)
        n2=isidenode(2,i)
        dps(i)=(dp(n1)+dp(n2))/2
      enddo !i
      do i=1,nea
        dpe(i)=1.e10
        do j=1,3
          if(dpe(i)>dp(elnode(j,i))) dpe(i)=dp(elnode(j,i))
        enddo !j
      enddo !i=1,nea

!...  Recompute vgrid and calculate rewetted pts
      if(inunfl==0) then
        call levels0(iths_main,it)
      else
        call levels1(iths_main,it)
      endif
      if(myrank==0) write(16,*) 'done recomputing levels...'

!...  Compute total mass of S,T (after levels are updated to 
!     approx. d/dt(total)
      if(1==1) then
        tot_heat=0
        tot_salt=0
        do i=1,ne
          if(idry_e(i)==1) cycle

          do k=kbe(i)+1,nvrt
            vol=(ze(k,i)-ze(k-1,i))*area(i)
            tot_heat=tot_heat+vol*tsel(1,k,i)
            tot_salt=tot_salt+vol*tsel(2,k,i)
          enddo !k
        enddo !i=1,ne

#ifdef INCLUDE_TIMING
        cwtmp=mpi_wtime()
#endif
        call mpi_allreduce(tot_heat,tot_heat_gb,1,rtype,MPI_SUM,comm,ierr)
        call mpi_allreduce(tot_salt,tot_salt_gb,1,rtype,MPI_SUM,comm,ierr)
#ifdef INCLUDE_TIMING
        wtimer(9,2)=wtimer(9,2)+mpi_wtime()-cwtmp
#endif
        if(myrank==0) write(91,*)time/86400,tot_heat_gb,tot_salt_gb

!...    Compute total tracers mass
!...
        if(ntracers>0 ) then
          swild(1:ntracers)=0
          do i=1,ne
            if(idry_e(i)==1) cycle
 
            do k=kbe(i)+1,nvrt
              vol=(ze(k,i)-ze(k-1,i))*area(i)
              swild(1:ntracers)=swild(1:ntracers)+vol*trel(1:ntracers,k,i)
            enddo !k
          enddo !i=1,ne
#ifdef INCLUDE_TIMING
          cwtmp=mpi_wtime()
#endif
          call mpi_allreduce(swild,swild3,ntracers,rtype,MPI_SUM,comm,ierr)
#ifdef INCLUDE_TIMING
          wtimer(9,2)=wtimer(9,2)+mpi_wtime()-cwtmp
#endif
          if(myrank==0) write(92,*)time/86400,swild3(1:ntracers)
        endif !ntracers>0
      endif !1==2

!...  Compute nodal vel. for output and next backtracking
      call nodalvel !(ifltype)

!...  Density (using new level indices)
      prho=-99
      do i=1,npa
        if(idry(i)==1) cycle
        do k=1,nvrt
          k2=max(k,kbp(i))

!new9: debug
#ifdef USE_TIMOR
          if(tr_nd(1,k2,i)<-98) then
            write(errmsg,*)'new9:',iplg(i),k,k2,tr_nd(1,k2,i),rhomud(1:ntracers,k2,i)
            call parallel_abort(errmsg)
          endif
#endif

          prho(k,i)=eqstate(1,tnd(k,i),snd(k,i)           &
#ifdef USE_SED
     &                      ,tr_nd(:,k,i),Srho(:)       &
#endif 
#ifdef USE_TIMOR
     &                      ,tr_nd(:,k2,i),rhomud(1:ntracers,k2,i),laddmud_d &
#endif 
     &                     )
        enddo !k
      enddo !i

      if(iupwind_t/=0) then
        erho=-99
        do i=1,nea
          if(idry_e(i)==1) cycle
          do k=1,nvrt
            k2=max(k,kbe(i))

#ifdef USE_TIMOR
            do m=1,ntracers
              swild(m)=sum(rhomud(m,k2,elnode(1:3,i)))/3
            enddo !m

!new9
            if(trel(1,k,i)<-98) then
              write(errmsg,*)'new9(2):',ielg(i),k,k2,swild(:),trel(1,k,i)
              call parallel_abort(errmsg)
            endif
#endif
            erho(k,i)=eqstate(2,tsel(1,k,i),tsel(2,k,i)   &
#ifdef USE_SED
     &                        ,trel(:,k,i),Srho(:)      &
#endif 
#ifdef USE_TIMOR
     &                        ,trel(:,k,i),swild(1:ntracers),laddmud_d &
#endif 
     &                       )          
          enddo !k
        enddo !i
      endif !iupwind_t

!...  Compute mean density profile at nodes or elements
      if(ibcc_mean==1.or.ihot==0.and.icst==2) then
        call mean_density
      else !other cases
        rho_mean=0
      endif

      if(myrank==0) write(16,*)'done density calculation...'

!...  Compute depth averaged h-vel.
!...  In pframe if ics=2
      dav=0
      do i=1,npa
        if(idry(i)==1) cycle
        do k=kbp(i),nvrt-1
          dav(1,i)=dav(1,i)+(uu2(k+1,i)+uu2(k,i))/2*(znl(k+1,i)-znl(k,i))
          dav(2,i)=dav(2,i)+(vv2(k+1,i)+vv2(k,i))/2*(znl(k+1,i)-znl(k,i))
        enddo !k
        htot=eta2(i)+dp(i)
        if(htot<=h0) then
          write(errmsg,*)'Impossible 24b:',it,i,eta2(i),dp(i),htot,h0,iplg(i)
          call parallel_abort(errmsg)
        endif
        dav(1:2,i)=dav(1:2,i)/htot

!       Max. dav (based on magnitude)
        dav_mag=sqrt(dav(1,i)**2+dav(2,i)**2)
        if(dav_mag>dav_maxmag(i)) then
          dav_maxmag(i)=dav_mag
          dav_max(1:2,i)=dav(1:2,i)
        endif
      enddo !i=1,npa

#ifdef INCLUDE_TIMING
! end computing levels
      wtmp2=mpi_wtime()
      wtimer(10,1)=wtimer(10,1)+wtmp2-wtmp1
! start flux compution
      wtmp1=wtmp2
#endif

!...  Optional computation of fluxes and total volume etc.
      if(iflux/=0) then
!--------------------------------------------------
!     Compute total mass etc.
      tvol=0 !total volume
      tmass=0 !total mass
      tpe=0 !total potential energy
      tkne=0 !total kinetic energy (quasi-2D only)
      enerf=0 !energy loss due to bottom friction; only correct for 2D model
      ener_ob=0 !total wave enery out of open bnds; only correct for 0 mean flows!
      do i=1,ne !residents only
        if(idry_e(i)==1) cycle

        n1=elnode(1,i)
        n2=elnode(2,i)
        n3=elnode(3,i)
        etam=(eta2(n1)+eta2(n2)+eta2(n3))/3
        tpe=tpe+0.5*rho0*grav*area(i)*etam**2
        av_dep=etam+(dp(n1)+dp(n2)+dp(n3))/3
        tvol=tvol+area(i)*av_dep
!        do k=kbe(i),nvrt-1
!          ah=(znl(k+1,n1)+znl(k+1,n2)+znl(k+1,n3)-znl(k,n1)-znl(k,n2)-znl(k,n3))/3
!        enddo !k

        do j=1,3 !node or side
          nd=elnode(j,i)
          do k=kbp(nd),nvrt-1
            tmass=tmass+area(i)*(prho(k,nd)+prho(k+1,nd))*(znl(k+1,nd)-znl(k,nd))/6
          enddo !k
          htot=eta2(nd)+dp(nd)
          if(htot<=h0) then
            write(errmsg,*)'Impossible dry (9):',ielg(i),j,iplg(nd),htot
            call parallel_abort(errmsg)
          endif

          isd=elside(j,i)
          do k=kbs(isd),nvrt-1
            vmag1=su2(k,isd)**2+sv2(k,isd)**2
            vmag2=su2(k+1,isd)**2+sv2(k+1,isd)**2
            tkne=tkne+rho0*area(i)/6*(zs(k+1,isd)-zs(k,isd))*(vmag1+vmag2)/2
          enddo !k

!         enerf only correct for quasi-2D model
          enerf=enerf+dt*area(i)/3*rho0*Cdp(nd)*sqrt(dav(1,nd)**2+dav(2,nd)**2)**3

!         ener_ob
          isd=elside(j,i)
          if(isbs(isd)>0) then !open bnd; no sharing between processes
            n1=isidenode(1,isd)
            n2=isidenode(2,isd)
            etam=(eta2(n1)+eta2(n2))/2
!Error: may not be accurate near poles
            vel_m1=(dav(1,n1)+dav(1,n2))/2 !both in ll frame
            vel_m2=(dav(2,n1)+dav(2,n2))/2
            ener_ob=ener_ob+rho0/2*sqrt(grav*dps(isd))*dt*(grav*etam**2+dps(isd)*(vel_m1**2+vel_m2**2))*distj(isd)
          endif
        enddo !j=1,3
      enddo !i=1,ne

      allocate(buf3(6)); buf3=0
      swild(1)=tvol; swild(2)=tmass; swild(3)=tpe; swild(4)=tkne; swild(5)=enerf; swild(6)=ener_ob
#ifdef INCLUDE_TIMING
      cwtmp=mpi_wtime()
#endif
      call mpi_reduce(swild,buf3,6,rtype,MPI_SUM,0,comm,ierr)
#ifdef INCLUDE_TIMING
      wtimer(11,2)=wtimer(11,2)+mpi_wtime()-cwtmp
#endif

      if(myrank==0) write(13,*)time/3600,buf3(1:4),buf3(3)+buf3(4),buf3(5:6)
      deallocate(buf3)

      !Fluxes
      !fluxes_vol(1:max_flreg): volume fluxes from region i to i-1, with i>=1
      !(i.e. excluding region -1)
      allocate(fluxes_vol(max_flreg),fluxes_vol_gb(max_flreg),stat=istat)
      if(istat/=0) call parallel_abort('STEP: fluxes_vol alloc')
      fluxes_vol=0 
      do i=1,ns
        if(idry_s(i)==1.or.isdel(2,i)==0) cycle

        !Wet internal side
        ie0=isdel(1,i); ie=isdel(2,i)
        if((iflux_e(ie0) .eq. -1) .or. (iflux_e(ie) .eq. -1)) cycle

        if(ie0<=0.or.ie<=0) call parallel_abort('STEP: isdel() out of bound') 
!'
        if(iflux_e(ie0) .ne. iflux_e(ie) .and. iabs(iflux_e(ie0) - iflux_e(ie)) .eq. 1) then
          if(associated(isgl(islg(i))%next)) then !interface side
            if(isgl(islg(i))%next%rank<myrank) cycle !already in the sum so skip
          endif

          itmp1=max(iflux_e(ie0),iflux_e(ie)) !'hi' region #
          itmp2=min(iflux_e(ie0),iflux_e(ie)) !'lo' region #
          if(itmp1<1) call parallel_abort('STEP: flux index <1')
          if(itmp1==iflux_e(ie0)) then
            fac=1
          else
            fac=-1
          endif

          do k=kbs(i),nvrt-1
            if(ics==1) then
              vnn=(su2(k+1,i)+su2(k,i))/2*sframe(1,1,i)+(sv2(k+1,i)+sv2(k,i))/2*sframe(2,1,i)
            else
              vnn=(su2(k+1,i)+su2(k,i))/2
            endif !ics
            ftmp=fac*distj(i)*(zs(k+1,i)-zs(k,i))*vnn
            fluxes_vol(itmp1)=fluxes_vol(itmp1)+ftmp
          enddo !k
        endif !side bording 2 regions
      enddo !i=1,ns

#ifdef INCLUDE_TIMING
      cwtmp=mpi_wtime()
#endif
      call mpi_reduce(fluxes_vol,fluxes_vol_gb,max_flreg,rtype,MPI_SUM,0,comm,ierr)
#ifdef INCLUDE_TIMING
      wtimer(11,2)=wtimer(11,2)+mpi_wtime()-cwtmp
#endif
      if(myrank==0) then
        write(9,'(f16.6,6000(1x,e14.4))')time/86400,fluxes_vol_gb(1:max_flreg)
        write(16,*)'done computing fluxes...'
      endif
      deallocate(fluxes_vol, fluxes_vol_gb)
!---------------------------------------------------------      
      endif !iflux ne 0
!...  end compute flux balance

#ifdef INCLUDE_TIMING
! end flux compution
      wtmp2=mpi_wtime()
      wtimer(11,1)=wtimer(11,1)+wtmp2-wtmp1
! Start timing global output section
      wtmp1=wtmp2
#endif

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Write global output data
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      do j=1,noutput
        if(iof(j)==1.and.mod(it,nspool)==0) then
          write(ichan(j)) real(time,out_rkind)
          write(ichan(j)) it
          write(ichan(j)) (real(eta2(i),out_rkind),i=1,np)
            if(j<=13) then
              if(j==1) then
                write(ichan(j)) (real(eta2(i),out_rkind),i=1,np)
              else if(j==2) then !.and.ihconsv.ne.0) then
                write(ichan(j)) (real(pr(i),out_rkind),i=1,np)
              else if(j==3.and.ihconsv/=0) then
                write(ichan(j)) (real(airt1(i),out_rkind),i=1,np)
              else if(j==4.and.ihconsv/=0) then
                write(ichan(j)) (real(shum1(i),out_rkind),i=1,np)
              else if(j==5.and.ihconsv/=0) then
                write(ichan(j)) (real(srad(i),out_rkind),i=1,np)
              else if(j==6.and.ihconsv/=0) then
                write(ichan(j)) (real(fluxsu(i),out_rkind),i=1,np)
              else if(j==7.and.ihconsv/=0) then
                write(ichan(j)) (real(fluxlu(i),out_rkind),i=1,np)
              else if(j==8.and.ihconsv/=0) then
                write(ichan(j)) (real(hradu(i),out_rkind),i=1,np)
              else if(j==9.and.ihconsv/=0) then
                write(ichan(j)) (real(hradd(i),out_rkind),i=1,np)
              else if(j==10.and.ihconsv/=0) then
                write(ichan(j)) (real(sflux(i),out_rkind),i=1,np)
              else if(j==11.and.isconsv/=0) then
                write(ichan(j)) (real(fluxevp(i),out_rkind),i=1,np)
              else if(j==12.and.isconsv/=0) then
                write(ichan(j)) (real(fluxprc(i),out_rkind),i=1,np)
              else if(j==13) then
                write(ichan(j)) (real(Cdp(i),out_rkind),i=1,np)
              else ! for undefined case as floatout=0 in old code
                write(ichan(j)) (0.0,i=1,np)
              endif
            else if(j<=16) then
              if(j==14) then
                if(nws==0) then
                  write(ichan(j))((0.0,0.0),i=1,np)
                else !in ll frame if ics=2
                  write(ichan(j)) (real(windx(i),out_rkind),real(windy(i),out_rkind),i=1,np)
                endif
              else if(j==15) then !in ll frame if ics=2
                write(ichan(j)) (real(tau(1,i),out_rkind),real(tau(2,i),out_rkind),i=1,np)
              else !j=16
                write(ichan(j)) (real(dav(1,i),out_rkind),real(dav(2,i),out_rkind),i=1,np)
              endif
            else if(j<27) then
                floatout=0 !for some undefined variables
                if(j==17) then
                  write(ichan(j)) ((real(ww2(k,i),out_rkind),k=max0(1,kbp00(i)),nvrt),i=1,np)
                else
                  if(j==18) then
                    do i=1,np
                      if(idry(i)==1) then
                        write(ichan(j)) (-99.,k=max0(1,kbp00(i)),nvrt)
                      else
                        write(ichan(j)) (real(tnd(k,i),out_rkind),k=max0(1,kbp00(i)),nvrt)
                      endif
                    enddo
                  else if(j==19) then
                    do i=1,np
                      if(idry(i)==1) then
                        write(ichan(j)) (-99.,k=max0(1,kbp00(i)),nvrt)
                      else
                        write(ichan(j)) (real(snd(k,i),out_rkind),k=max0(1,kbp00(i)),nvrt)
                      endif
                    enddo
                  else if(j==20) then
                    write(ichan(j)) ((real(prho(k,i),out_rkind),k=max0(1,kbp00(i)),nvrt),i=1,np)
                  else if(j==21) then
                    write(ichan(j)) ((real(dfh(k,i),out_rkind),k=max0(1,kbp00(i)),nvrt),i=1,np)
                  else if(j==22) then
                    write(ichan(j)) ((real(dfv(k,i),out_rkind),k=max0(1,kbp00(i)),nvrt),i=1,np)
                  else if(j==23) then
                    write(ichan(j)) ((real(q2(k,i),out_rkind),k=max0(1,kbp00(i)),nvrt),i=1,np)
                  else if(j==24) then
                    write(ichan(j)) ((real(xl(k,i),out_rkind),k=max0(1,kbp00(i)),nvrt),i=1,np)
                  else if(j==25) then
                    do i=1,np
                      if(idry(i)==1) then
                        write(ichan(j)) (0.,k=max0(1,kbp00(i)),nvrt)
                      else
                        write(ichan(j)) (real(znl(max0(k,kbp(i)),i),out_rkind),k=max0(1,kbp00(i)),nvrt)
                      endif
                    enddo
                  else if(j==26) then
                      write(ichan(j)) ((real(qnon(k,i),out_rkind),k=max0(1,kbp00(i)),nvrt),i=1,np)
                  endif
                endif
            else if(j==27) then
              write(ichan(j)) ((real(uu2(k,i),out_rkind),real(vv2(k,i),out_rkind),k=max0(1,kbp00(i)),nvrt),i=1,np)
            else if(j<=27+ntracers) then !tracers; implies ntracers>0
              write(ichan(j)) ((real(tr_nd(j-27,k,i),out_rkind),k=max0(1,kbp00(i)),nvrt),i=1,np)
            else !optional modules; MUST BE IN THE SAME ORDER AS BEFORE
#ifdef USE_SED
              if(flag_model/=1) call parallel_abort('MAIN: strange output (2)')
!'
              if(j<=indx_out(1,2)) then

                if(j==indx_out(1,1)) then
                  ! depth.61
                  write(ichan(j)) (real(dp(i)-dp00(i),out_rkind),i=1,np)
                else if(j<=indx_out(1,1)+ntracers) then
                  !qbdl_n.62
                  write(ichan(j)) (real(bedldu(i,j-indx_out(1,1)),out_rkind),real(bedldv(i,j-indx_out(1,1)),out_rkind),i=1,np)
                else if(j<=indx_out(1,1)+2*ntracers) then
                  !bfrac_n.61 (top layer only)
                  write(ichan(j)) (real(bed_fracn(i,j-indx_out(1,1)-ntracers),out_rkind),i=1,np)
                else if(j==indx_out(1,1)+1+2*ntracers) then
                  !bedd50.61
                  write(ichan(j)) (real(bed_d50n(i)*1000.d0,out_rkind),i=1,np) ! in mm
                else if (j==indx_out(1,1)+2+2*ntracers) then
                  !bstress.61
                  write(ichan(j)) (real(bed_taun(i)*rho0,out_rkind),i=1,np) ! in N.m-2
                else if (j==indx_out(1,1)+3+2*ntracers) then
                  !brough.61
                  write(ichan(j)) (real(bed_rough(i)*1000.d0,out_rkind),i=1,np) ! in mm
                endif
              endif !scope of SED model
#endif /*USE_SED*/

#ifdef USE_SED2D
              if((j>=indx_out(1,1)).and.(j<=indx_out(1,2))) then          
                if(j>=indx_out(1,1).and.(j<=indx_out(1,1)+3)) then !scalar
                  if(j==indx_out(1,1)) then
                     write(ichan(j)) (real(dp(i)-dp00(i),out_rkind),i=1,np)
                  else if(j==indx_out(1,1)+1) then
                     write(ichan(j)) (real(Cdsed(i),out_rkind),i=1,np)
                  else if(j==indx_out(1,1)+2) then
                     write(ichan(j)) (real(cflsed(i),out_rkind),i=1,np)
                  else if(j==indx_out(1,1)+3) then
                     write(ichan(j)) (real(d50moy(i,1),out_rkind),i=1,np)
                  endif
                else if(j>indx_out(1,1)+3) then !vector
                  if(j==indx_out(1,1)+4) then
                     write(ichan(j)) (real(qtot(i,1),out_rkind),real(qtot(i,2),out_rkind),i=1,np)
                  else if(j==indx_out(1,1)+5) then
                     write(ichan(j)) (real(qs(i,1),out_rkind),real(qs(i,2),out_rkind),i=1,np)
                  else if(j==indx_out(1,1)+6) then
                     write(ichan(j)) (real(qb(i,1),out_rkind),real(qb(i,2),out_rkind),i=1,np)
                  else if(j==indx_out(1,1)+7) then 
                     write(ichan(j)) (real(dpdxy(i,1),out_rkind),real(dpdxy(i,2),out_rkind),i=1,np)
                  else if(j==indx_out(1,1)+8) then
                     write(ichan(j)) (real(qav(i,1),out_rkind),real(qav(i,2),out_rkind),i=1,np)
                  endif
                endif
              endif !scope of SED2D model
#endif /*USE_SED2D*/

#ifdef USE_NAPZD
              if(j<=indx_out(2,2)) then
                if(j==indx_out(2,1)) then
                  write(ichan(j)) ((real(Bio_bdefp(k,i),out_rkind),k=max0(1,kbp00(i)),nvrt),i=1,np)
                else !total N
                  write(ichan(j)) ((real(sum(tr_nd(1:4,k,i)),out_rkind),k=max0(1,kbp00(i)),nvrt),i=1,np)
                endif
              endif !scope of NAPZD
#endif /*USE_NAPZD*/

#ifdef USE_WWM
              if((j>=indx_out(3,1)).and.(j<=indx_out(3,2))) then
                if(j<=indx_out(3,2)-2) then !scalar
                  itmp=j-indx_out(3,1)+1
                  if(itmp>24) call parallel_abort('MAIN: wwm_out over')
                  write(ichan(j)) (real(out_wwm(i,indx_wwm_out(itmp)),out_rkind),i=1,np)
                else !vectors
                  if (j==indx_out(3,2)-1) then
                    write(ichan(j)) (real(out_wwm(i,8),out_rkind),real(out_wwm(i,7),out_rkind),i=1,np)
                  else if (j==indx_out(3,2)) then
                    write(ichan(j)) (real(out_wwm(i,27),out_rkind),real(out_wwm(i,28),out_rkind),i=1,np)
                  endif
                endif !j
              endif !scope of WWM; j<=indx_out(3,2)
#endif /*USE_WWM*/

              if(flag_model==0) then; if(j>=indx_out(4,1).and.j<=indx_out(4,2)) then !age
                do i=1,np
                  do k=max0(1,kbp00(i)),nvrt
                    tmp1=max(1.d-5,tr_nd(j-indx_out(4,1)+1,k,i))
                    floatout=tr_nd(j-indx_out(4,1)+1+ntracers/2,k,i)/tmp1/86400
                    write(ichan(j)) real(floatout,out_rkind)
                  enddo !k
                enddo !i
              endif; endif !scope of age
            endif !j
          if(myrank==0) write(16,'(a48)')'done outputting '//variable_nm(j)
        endif !iof(j).eq.1.and.mod(it,nspool).eq.0
      enddo !j=1,noutput

      
!...  Non-standard outputs
      if(iof_ns(1)==1) then 
        call schism_output_custom(lwrite,6,2,201,'hvel',nvrt,nsa,su2,sv2)
        if(myrank==0.and.lwrite==1) write(16,*)'done outputting hvel.67'
      endif !iof_ns
      if(iof_ns(2)==1) then 
        call schism_output_custom(lwrite,8,1,202,'vert',nvrt,nea,we_fv)
        if(myrank==0.and.lwrite==1) write(16,*)'done outputting vert.68'
      endif !iof_ns
      if(iof_ns(3)==1) then 
        call schism_output_custom(lwrite,9,1,203,'temp',nvrt,nea,tsel(1,:,:))
        if(myrank==0.and.lwrite==1) write(16,*)'done outputting temp.70'
      endif !iof_ns
      if(iof_ns(4)==1) then 
        call schism_output_custom(lwrite,9,1,204,'salt',nvrt,nea,tsel(2,:,:))
        if(myrank==0.and.lwrite==1) write(16,*)'done outputting salt.70'
      endif !iof_ns
#ifdef USE_SED
      if(iof_ns(5)==1)then
        call schism_output_custom(lwrite,5,1,209,'z0st',1,nea,bottom(:,izbld))
        if(myrank==0.and.lwrite==1) write(16,*)'done outputting z0st.66'
      endif !iof_ns
      if(iof_ns(7)==1)then
        call schism_output_custom(lwrite,5,1,206,'z0cr',1,nea,bottom(:,izcr))
        if(myrank==0.and.lwrite==1) write(16,*)'done outputting z0cr.66'
      endif !iof_ns
      if(iof_ns(8)==1)then
        call schism_output_custom(lwrite,5,1,207,'z0sw',1,nea,bottom(:,izsw))
        if(myrank==0.and.lwrite==1) write(16,*)'done outputting z0sw.66'
      endif !iof_ns
      if(iof_ns(9)==1)then
        call schism_output_custom(lwrite,5,1,208,'z0wr',1,nea,bottom(:,izwr))
        if(myrank==0.and.lwrite==1) write(16,*)'done outputting z0wr.66'
      endif !iof_ns
      if(iof_ns(19)==1)then
        call schism_output_custom(lwrite,5,1,219,'bthk',1,nea,sum(bed(:,:,ithck),1))
        if(myrank==0.and.lwrite==1) write(16,*)'done outputting bthk.66'
      endif !iof_ns
      if(iof_ns(20)==1)then
        call schism_output_custom(lwrite,5,1,220,'bage',1,nea,sum(bed(:,:,iaged),1))
        if(myrank==0.and.lwrite==1) write(16,*)'done outputting bage.66'
      endif !iof_ns
#elif USE_SED2D
      if(iof_ns(6)==1)then
        call schism_output_custom(lwrite,5,1,206,'z0eq',1,nea,z0_e(:))
        if(myrank==0.and.lwrite==1) write(16,*)'done outputting z0eq.66'
      endif !iof_ns
      if(iof_ns(7)==1)then
        call schism_output_custom(lwrite,5,1,207,'z0cr',1,nea,z0cr_e(:))
        if(myrank==0.and.lwrite==1) write(16,*)'done outputting z0cr.66'
      endif !iof_ns
      if(iof_ns(8)==1)then
        call schism_output_custom(lwrite,5,1,208,'z0sw',1,nea,z0sw_e(:))
        if(myrank==0.and.lwrite==1) write(16,*)'done outputting z0sw.66'
      endif !iof_ns
      if(iof_ns(9)==1)then
        call schism_output_custom(lwrite,5,1,209,'z0wr',1,nea,z0wr_e(:))
        if(myrank==0.and.lwrite==1) write(16,*)'done outputting z0wr.66'
      endif !iof_ns
#endif

#ifdef DEBUG
      if(iof_ns(10)==1) then
        call schism_output_custom(lwrite,4,2,210,'bpgr',1,nsa,bpgr(:,1),bpgr(:,2))
        if(myrank==0.and.lwrite==1) write(16,*)'done outputting bpgr.65'
      endif !iof_ns
#ifdef USE_WWM
      if(iof_ns(11)==1) then
        call schism_output_custom(lwrite,6,2,211,'wafo',nvrt,nsa,wafo(:,:,1),wafo(:,:,2))
        if(myrank==0.and.lwrite==1) write(16,*)'done outputting wafo.67'
      endif !iof_ns
#endif /*USE_WWM*/
#endif /*DEBUG*/

#ifdef USE_ICM
      if(iof_ns(12)==1) then
        call schism_output_custom(lwrite,5,1,212,'bdoc',1,nea,BENDOC)
        if(myrank==0.and.lwrite==1) write(16,*)'done outputting bdoc.66'
      endif 
      if(iof_ns(13)==1) then
        call schism_output_custom(lwrite,5,1,213,'bnh4',1,nea,SED_BENNH4)
        if(myrank==0.and.lwrite==1) write(16,*)'done outputting bnh4.66'
      endif 
      if(iof_ns(14)==1) then
        call schism_output_custom(lwrite,5,1,214,'bno3',1,nea,SED_BENNO3)
        if(myrank==0.and.lwrite==1) write(16,*)'done outputting bno3.66'
      endif 
      if(iof_ns(15)==1) then
        call schism_output_custom(lwrite,5,1,215,'bpo4',1,nea,BENPO4)
        if(myrank==0.and.lwrite==1) write(16,*)'done outputting bpo4.66'
      endif 
      if(iof_ns(16)==1) then
        call schism_output_custom(lwrite,5,1,216,'bcod',1,nea,SED_BENCOD)
        if(myrank==0.and.lwrite==1) write(16,*)'done outputting bcod.66'
      endif 
      if(iof_ns(17)==1) then
        call schism_output_custom(lwrite,5,1,217,'sbdo',1,nea,sed_BENDO)
        if(myrank==0.and.lwrite==1) write(16,*)'done outputting sbdo.66'
      endif 
      if(iof_ns(18)==1) then
        call schism_output_custom(lwrite,5,1,218,'sbsa',1,nea,BENSA)
        if(myrank==0.and.lwrite==1) write(16,*)'done outputting sbsa.66'
      endif 
#endif /*USE_ICM*/

!     Test 
!      call schism_output_custom(lwrite,10,1,205,'elev',1,npa,eta2)
!      if(myrank==0.and.lwrite==1) write(16,*)'done outputting elev.71'

!     Open new global output files and write header data
      !if(it==ifile*ihfskip) then !.and.it/=ntime) then
      if(mod(it,ihfskip)==0) then
        ifile=ifile+1                   !output file #
        write(ifile_char,'(i12)') ifile !convert ifile to a string
        ifile_char=adjustl(ifile_char)  !place blanks at end
        ifile_len=len_trim(ifile_char)  !length without trailing blanks
        fgb=ifile_char(1:ifile_len)//'_0000'; lfgb=len_trim(fgb);
        write(fgb(lfgb-3:lfgb),'(i4.4)') myrank
        do i=1,noutput
          ichan(i)=100+i !output channel #
          if(iof(i)==1) then
            open(ichan(i),file='outputs/'//(fgb(1:lfgb)//'_'//outfile(i)),status='replace',form="unformatted",access="stream")
          endif
        enddo !i
      endif !it==ifile*ihfskip

!...  Station outputs
      if(iout_sta/=0) then
        do j=1,nvar_sta
          if(iof_sta(j)==0.or.mod(it,nspool_sta)/=0) cycle

          do i=1,nout_sta
            ie=iep_sta(i)
            if(ie==0) then !not parent
              iep_flag(i)=0 !for comm. later
              sta_out(i,j)=0
              sta_out3d(:,i,j)=0
              zta_out3d(:,i,j)=0
            else !is parent
              iep_flag(i)=1
              sta_out(i,j)=0 !initialize
              select case(j)
                case(1) !elev.
                  swild2(1,1:3)=eta2(elnode(1:3,ie))
                case(2) !air pressure
                  swild2(1,1:3)=pr(elnode(1:3,ie))
                case(3) !wind x
                  swild2(1,1:3)=windx(elnode(1:3,ie))
                case(4) !wind y
                  swild2(1,1:3)=windy(elnode(1:3,ie))
                case(5) !T
                  swild2(1:nvrt,1:3)=tnd(1:nvrt,elnode(1:3,ie))
                case(6) !S
                  swild2(1:nvrt,1:3)=snd(1:nvrt,elnode(1:3,ie))
                case(7) !u
!Error: may not be accurate near poles as pframe changes rapidly there
                  swild2(1:nvrt,1:3)=uu2(1:nvrt,elnode(1:3,ie))
                case(8) !v
                  swild2(1:nvrt,1:3)=vv2(1:nvrt,elnode(1:3,ie))
                case(9) !w
                  swild2(1:nvrt,1:3)=ww2(1:nvrt,elnode(1:3,ie))
                case default
                  call parallel_abort('MAIN: unknown sta. output')
              end select

              if(j<=4) then !2D var.
                sta_out(i,j)=sum(arco_sta(i,1:3)*swild2(1,1:3))
              else !3D var.
                if(idry_e(ie)==1) then !dry
                  sta_out(i,j)=-999.
                  sta_out3d(:,i,j)=-999.
                  zta_out3d(:,i,j)=-999.
                else !wet
                  do m=1,3 !wet nodes
                    nd=elnode(m,ie)
                    !Vertical interplation
                    if(zstal(i)<=znl(kbp(nd),nd)) then
                      k0=kbp(nd); zrat=0
                    else if(zstal(i)>=znl(nvrt,nd)) then
                      k0=nvrt-1; zrat=1
                    else
                      k0=0
                      do k=kbp(nd),nvrt-1
                        if(zstal(i)>=znl(k,nd).and.zstal(i)<=znl(k+1,nd)) then
                          k0=k
                          zrat=(zstal(i)-znl(k,nd))/(znl(k+1,nd)-znl(k,nd))
                          exit
                        endif
                      enddo !k
                      if(k0==0) call parallel_abort('MAIN: impossible (71)')
!'
                    endif !zstal
                    swild(m)=swild2(k0,m)*(1-zrat)+swild2(k0+1,m)*zrat
                  enddo !m=1,3

                  !Horizonal interplation
                  sta_out(i,j)=sum(arco_sta(i,1:3)*swild(1:3))

                  !Vertical profiles
                  do k=1,nvrt
                    do m=1,3
                      nd=elnode(m,ie)
                      if(k<kbp(nd)) then
                        swild4(1,m)=-9999 !zcor
                        swild4(2,m)=-9999 !var
                      else
                        swild4(1,m)=znl(k,nd)
                        swild4(2,m)=swild2(k,m)
                      endif
                    enddo !m
                    zta_out3d(k,i,j)=sum(arco_sta(i,1:3)*swild4(1,1:3))
                    sta_out3d(k,i,j)=sum(arco_sta(i,1:3)*swild4(2,1:3))
                  enddo !k
                endif !idry_e
              endif !j
            endif !ie
          enddo !i=1,nout_sta
!Debug
!          write(12,*)j,time,sta_out(:,j)
        enddo !j=1,nvar_sta

!       Output by rank 0
        call mpi_reduce(iep_flag,nwild2,nout_sta,itype,MPI_SUM,0,comm,ierr)
        call mpi_reduce(sta_out,sta_out_gb,nout_sta*nvar_sta,rtype,MPI_SUM,0,comm,ierr)
        call mpi_reduce(sta_out3d,sta_out3d_gb,nvrt*nout_sta*nvar_sta,rtype,MPI_SUM,0,comm,ierr)
        call mpi_reduce(zta_out3d,zta_out3d_gb,nvrt*nout_sta*nvar_sta,rtype,MPI_SUM,0,comm,ierr)

        if(myrank==0) then
!          write(290,*)nwild2(1:nout_sta)
          do i=1,nvar_sta
            if(iof_sta(i)==0.or.mod(it,nspool_sta)/=0) cycle
            do j=1,nout_sta
              if(nwild2(j)==0) then
                sta_out_gb(j,i)=-9999
                if(i>4) then !3D only
                  sta_out3d_gb(:,j,i)=-9999.
                  zta_out3d_gb(:,j,i)=-9999.
                endif
              else
                sta_out_gb(j,i)=sta_out_gb(j,i)/nwild2(j)
                if(i>4) then !3D only
                  sta_out3d_gb(:,j,i)=sta_out3d_gb(:,j,i)/nwild2(j) 
                  zta_out3d_gb(:,j,i)=zta_out3d_gb(:,j,i)/nwild2(j) 
                endif
              endif
            enddo !j
            write(250+i,'(e14.6,6000(1x,e14.6))')time,sta_out_gb(:,i)
!            if(i>4) write(250+i,'(e14.6,100000(1x,e14.6))')time,sta_out3d_gb(:,:,i),zta_out3d_gb(:,:,i)
          enddo !i
          write(16,*)'done station outputs...'
        endif !myrank
      endif !iout_sta/=0

#ifdef USE_HA
!...
!...  IF iharind=1 AND THE TIME STEP FALLS WITHIN THE SPECIFIED WINDOW
!...  AND ON THE SPECIFIED INCREMENT, USE MODEL RESULTS TO UPDATE
!...  HARMONIC ANALYSIS MATRIX AND LOAD VECTORS.  NOTE: AN 8 BYTE RECORD
!...  SHOULD BE USED THROUGHOUT THE HARMONIC ANALYSIS SUBROUTINES, EVEN
!...  ON 32 BIT WORKSTATIONS, SINCE IN THAT CASE THE HARMONIC ANALYSIS
!...  IS DONE IN DOUBLE PRECISION.
!...  Adapted from ADCIRC
      IF(iharind.EQ.1) THEN
         IF((it.GT.ITHAS).AND.(it.LE.ITHAF)) THEN
            IF(ICHA.EQ.NHAINC) ICHA=0
            ICHA=ICHA+1
            IF(ICHA.EQ.NHAINC) THEN
!...
!.....UPDATE THE LHS MATRIX
!...
               CALL LSQUPDLHS(time,it)
!... 
!.....IF DESIRED UPDATE GLOBAL ELEVATION LOAD VECTOR
!... 
               IF(NHAGE.EQ.1) CALL LSQUPDEG(ETA2,np)
!... 
!.....IF DESIRED UPDATE GLOBAL VELOCITY LOAD VECTOR
!               IF(NHAGV.EQ.1) CALL LSQUPDVG(UHA,VHA,np)

            ENDIF
         ENDIF

!...  LINES TO COMPUTE MEANS AND VARIANCES

         if (CHARMV) then
            IF(it.GT.ITMV) THEN
               NTSTEPS=NTSTEPS+1
               DO I=1,np
                  ELAV(I)=ELAV(I)+ETA2(I)
!                  XVELAV(I)=XVELAV(I)+UHA(I)
!                  YVELAV(I)=YVELAV(I)+VHA(I)
                  ELVA(I)=ELVA(I)+ETA2(I)*ETA2(I)
!                  XVELVA(I)=XVELVA(I)+UHA(I)*UHA(I)
!                  YVELVA(I)=YVELVA(I)+VHA(I)*VHA(I)
               END DO
            ENDIF
         endif                  !   charmv
      ENDIF
#endif /*USE_HA*/

#ifdef INCLUDE_TIMING
! End timing global output section
      wtmp2=mpi_wtime()
      wtimer(12,1)=wtimer(12,1)+wtmp2-wtmp1
! Start timing write hotstart section
      wtmp1=wtmp2
#endif

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Write hot start data
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      if(nhot==1.and.mod(it,nhot_write)==0) then
        !Flags for each module that needs hotstart outputs
        nwild=0 !init.
#ifdef USE_ICM
        nwild(1)=1
#endif
#ifdef USE_SED2D 
        nwild(2)=1
#endif
#ifdef USE_SED
        nwild(3)=1
#endif
#ifdef USE_HA
        nwild(4)=1
#endif

        write(it_char,'(i72)')it
        it_char=adjustl(it_char)
        lit=len_trim(it_char)
        it_char=it_char(1:lit)//'_0000'; lit=len_trim(it_char)
        write(it_char(lit-3:lit),'(i4.4)') myrank
        !Reserve 8 bytes for all integers as well
        ihot_len=8*(7+((3+2*ntracers)*nvrt+1)*ne+(4*nvrt+1)*ns+(2+11*nvrt)*np)
        open(36,file='outputs/'//it_char(1:lit)//'_hotstart', &
             access='direct',recl=ihot_len,status='replace') 

        write(36,rec=1)nwild(1:4),dble(time),it,ifile,(idry_e(i),(dble(we(j,i)),dble(tsel(1:2,j,i)), &
     &(dble(trel0(l,j,i)),dble(trel(l,j,i)),l=1,ntracers),j=1,nvrt),i=1,ne), &
     &(idry_s(i),(dble(su2(j,i)),dble(sv2(j,i)),dble(tsd(j,i)),dble(ssd(j,i)),j=1,nvrt),i=1,ns), &
     &(dble(eta2(i)),idry(i),(dble(tnd(j,i)),dble(snd(j,i)),dble(tem0(j,i)),dble(sal0(j,i)), &
     &dble(q2(j,i)),dble(xl(j,i)),dble(dfv(j,i)),dble(dfh(j,i)),dble(dfq1(j,i)),dble(dfq2(j,i)), &
     &dble(qnon(j,i)),j=1,nvrt),i=1,np) 
        close(36)

        !Save starting record # for other modules, assuming 8-byte per record
        IHOTSTP=ihot_len/8
#ifdef USE_ICM
        open(36,file='outputs/'//it_char(1:lit)//'_hotstart', &
             &access='direct',recl=8,status='old')
        do i=1,ne
          write(36,rec=IHOTSTP+1)sed_BENDO(i)
          write(36,rec=IHOTSTP+2)CTEMP(i)
          write(36,rec=IHOTSTP+3)BBM(i)
          write(36,rec=IHOTSTP+4)CPOS(i)
          write(36,rec=IHOTSTP+5)PO4T2TM1S(i)
          write(36,rec=IHOTSTP+6)NH4T2TM1S(i)
          write(36,rec=IHOTSTP+7)NO3T2TM1S(i)
          write(36,rec=IHOTSTP+8)HST2TM1S(i)
          write(36,rec=IHOTSTP+9)CH4T2TM1S(i)
          write(36,rec=IHOTSTP+10)CH41TM1S(i)
          write(36,rec=IHOTSTP+11)SO4T2TM1S(i)
          write(36,rec=IHOTSTP+12)SIT2TM1S(i)
          write(36,rec=IHOTSTP+13)BENSTR1S(i)
          write(36,rec=IHOTSTP+14)CPOP(i,1)
          write(36,rec=IHOTSTP+15)CPOP(i,2)
          write(36,rec=IHOTSTP+16)CPOP(i,3)
          write(36,rec=IHOTSTP+17)CPON(i,1)
          write(36,rec=IHOTSTP+18)CPON(i,2)
          write(36,rec=IHOTSTP+19)CPON(i,3)
          write(36,rec=IHOTSTP+20)CPOC(i,1)
          write(36,rec=IHOTSTP+21)CPOC(i,2)
          write(36,rec=IHOTSTP+22)CPOC(i,3)
          write(36,rec=IHOTSTP+23)NH41TM1S(i)
          write(36,rec=IHOTSTP+24)NO31TM1S(i)
          write(36,rec=IHOTSTP+25)HS1TM1S(i)
          write(36,rec=IHOTSTP+26)SI1TM1S(i)
          write(36,rec=IHOTSTP+27)PO41TM1S(i)
          write(36,rec=IHOTSTP+28)PON1TM1S(i)
          write(36,rec=IHOTSTP+29)PON2TM1S(i)
          write(36,rec=IHOTSTP+30)PON3TM1S(i)
          write(36,rec=IHOTSTP+31)POC1TM1S(i)
          write(36,rec=IHOTSTP+32)POC2TM1S(i)
          write(36,rec=IHOTSTP+33)POC3TM1S(i)
          write(36,rec=IHOTSTP+34)POP1TM1S(i)
          write(36,rec=IHOTSTP+35)POP2TM1S(i)
          write(36,rec=IHOTSTP+36)POP3TM1S(i)
          write(36,rec=IHOTSTP+37)PSITM1S(i)
          write(36,rec=IHOTSTP+38)BFORMAXS(i)
          write(36,rec=IHOTSTP+39)ISWBENS(i)
          write(36,rec=IHOTSTP+40)DFEEDM1S(i)
          IHOTSTP=IHOTSTP+40
        enddo !i=1,ne
        close(36)
#endif /*USE_ICM*/

        !write(12,*)'After hot trcr:',it,real(trel),real(trel0)

#ifdef USE_SED2D 
        open(36,file='outputs/'//it_char(1:lit)//'_hotstart', &
             &access='direct',recl=8,status='old')
        do i=1,np
          write(36,rec=IHOTSTP+1)dp(i)
          IHOTSTP=IHOTSTP+1
        enddo !i=1,np
        close(36)
#endif /*USE_SED2D*/

#ifdef USE_SED
        open(36,file='outputs/'//it_char(1:lit)//'_hotstart', &
             &access='direct',recl=8,status='old')
        write(36,rec=IHOTSTP+1)MBEDP
        write(36,rec=IHOTSTP+2)Nbed
        IHOTSTP=IHOTSTP+2
        do i=1,np
          write(36,rec=IHOTSTP+1)dp(i)
          write(36,rec=IHOTSTP+2)rough_p(i)
          IHOTSTP=IHOTSTP+2
        enddo !i=1,np

        do i=1,MBEDP
          do j=1,ne
            do k=1,Nbed
              write(36,rec=IHOTSTP+1)bed(k,j,i)
              IHOTSTP=IHOTSTP+1
            enddo !k
          enddo !j
        enddo !i

        do i=1,ntracers
          do k=1,ne
            do m=1,Nbed
              write(36,rec=IHOTSTP+1)bed_frac(m,k,i)
              IHOTSTP=IHOTSTP+1
            enddo !m
          enddo !k
        enddo !i

        close(36)
#endif /*USE_SED*/

#ifdef USE_HA
!...  not working properly yet
!...    IF APPROPRIATE ADD HARMONIC ANALYSIS INFORMATION TO HOT START FILE
!...    Adapted from ADCIRC
        open(36,file='outputs/'//it_char(1:lit)//'_hotstart', &
             &access='direct',recl=8,status='old') 
        !IHOTSTP = 3+((1+(3+2*ntracers)*nvrt)*ne)+((1+4*nvrt)*ns)+((2+9*nvrt)*np)
        IF((iharind.EQ.1).AND.(it.GT.ITHAS)) THEN
           WRITE(36,REC=IHOTSTP+1) ICHA
           IHOTSTP = IHOTSTP + 1
           CALL HAHOUT(np,0,0,0,0,NHAGE,NHAGV,36,IHOTSTP)

           IF(NHAGE.EQ.1) CALL HAHOUTEG(np,36,IHOTSTP)
           IF(NHAGV.EQ.1) CALL HAHOUTVG(np,36,IHOTSTP)
        ENDIF

        if (CHARMV) then
           IF((iharind.EQ.1).AND.(it.GT.ITMV)) THEN
              IHOTSTP=IHOTSTP+1
              WRITE(36,REC=IHOTSTP) NTSTEPS
              IF(NHAGE.EQ.1) THEN
                 DO I=1,np
                    WRITE(36,REC=IHOTSTP+1) ELAV(I)
                    WRITE(36,REC=IHOTSTP+2) ELVA(I)
                    IHOTSTP=IHOTSTP+2
                 END DO
              ENDIF
              IF(NHAGV.EQ.1) THEN
                 DO I=1,np
                    WRITE(36,REC=IHOTSTP+1) XVELAV(I)
                    WRITE(36,REC=IHOTSTP+2) YVELAV(I)
                    WRITE(36,REC=IHOTSTP+3) XVELVA(I)
                    WRITE(36,REC=IHOTSTP+4) YVELVA(I)
                    IHOTSTP=IHOTSTP+4
                 END DO
              ENDIF
           ENDIF
        endif               !  charmv
        close(36)
#endif /*USE_HA*/

        if(myrank==0) write(16,*) 'hot start written',it,time,ifile
      endif !nhot

#ifdef INCLUDE_TIMING
! End hotstart output section
      wtmp2=mpi_wtime()
      wtimer(13,1)=wtimer(13,1)+wtmp2-wtmp1
#endif

      if(myrank==0) then
        write(16,'(a,i12,a,f20.6)') 'TIME STEP= ',it,';  TIME= ',time
!'
        flush(16) !flush "mirror.out" for every time step
      endif
      call parallel_barrier !synchronize before starting next time step


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! End Time Stepping
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!      enddo !it

      first_call=.false.

!     Temp. fix for Stampede cluster problem
#ifdef STAMPEDE
      if(myrank==0) then
        open(32,file='die.stam',status='old')
        read(32,*)istat
        if(istat/=0) call parallel_abort('Aborting due to die.stam')
        close(32)
      endif
#endif

!     Deallocate temp. arrays to avoid memory leak
      if(if_source==1) deallocate(msource)
      if(nonhydro==1) deallocate(qhat,dqnon_dxy,qmatr,qir)
      deallocate(hp_int)

      if(ntracers>0.and.inu_tr==2) deallocate(swild9)

#ifdef DEBUG
      deallocate(bpgr,wafo)
#endif

#ifdef USE_NAPZD
      deallocate(Bio_bdefp)
#endif

#ifdef USE_SED
      deallocate(tr_tc,tr_tl)
#endif
 
#ifdef USE_WWM
      deallocate(stokes_w,stokes_w_nd,stokes_vel_sd)
#endif

      end subroutine schism_step

