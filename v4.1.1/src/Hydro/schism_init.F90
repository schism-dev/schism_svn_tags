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

!===============================================================================
!===============================================================================
! Initialize SCHISM 
!===============================================================================
!===============================================================================

      subroutine schism_init(iths,ntime)

!     Most mpi fortran compiler has mpif.h
!#ifdef USE_MPIMODULE
!      use mpi
!#endif
      use schism_glbl
      use schism_msgp

#ifdef USE_GOTM
      use turbulence, only: init_turbulence, cde, tke1d => tke, eps1d => eps, L1d => L, num1d => num, nuh1d => nuh
      use mtridiagonal, only: init_tridiagonal
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
                          &PSAt,PCOD,PDO     !added by YC
      USE icm_sed_mod, only: sed_BENDO,CTEMP,BBM,CPOS,PO4T2TM1S,NH4T2TM1S,NO3T2TM1S, &
                              &HST2TM1S,CH4T2TM1S,CH41TM1S,SO4T2TM1S,SIT2TM1S,BENSTR1S,CPOP,CPON,CPOC,  &
                              &NH41TM1S,NO31TM1S,HS1TM1S,SI1TM1S,PO41TM1S,PON1TM1S,PON2TM1S,PON3TM1S,POC1TM1S,POC2TM1S,&
                              &POC3TM1S,POP1TM1S,POP2TM1S,POP3TM1S,PSITM1S,BFORMAXS,ISWBENS,DFEEDM1S  !added YC 
#endif

#ifdef USE_NAPZD
      USE biology_napzd
#endif

#ifdef USE_SED
       USE sed_mod, only : Srho,Nbed,MBEDP,bed,bed_frac
#endif



#ifdef USE_OIL
#endif

#ifdef USE_HA
      USE harm
#endif

      USE hydraulic_structures

      implicit none
!#ifndef USE_MPIMODULE
      include 'mpif.h'
!#endif

      integer, intent(out) :: iths,ntime

!     Functions
      integer :: lindex_s
      real(rkind) :: covar,signa,eqstate

!     Output handles
      character(len=2) :: stringvalue
      character(len=8) :: date
      character(len=10) :: timestamp
!      character(len=72) :: it_char
      character(len=72) :: fgb  ! Processor specific global output file name
      integer :: lfgb       ! Length of processor specific global output file name
      real(4) :: floatout,floatout2
      real(rkind) :: double1 !for hotstart.in

!     Misc. arrays
      integer, allocatable :: ipiv(:)
      integer, allocatable :: nwild(:),nwild2(:),ibuf1(:,:),ibuf2(:,:)
      real(rkind), allocatable :: akr(:,:),akrp(:),work4(:),z_r2(:),xy_kr(:,:)
      real(rkind), allocatable :: swild(:),swild2(:,:)
      real(rkind), allocatable :: swild3(:),rwild(:,:)
      real(rkind), allocatable :: swild4(:,:),swild10(:,:) !double precision for hotstart.in
      real(rkind), allocatable :: swild99(:,:),swild98(:,:,:) !used for exchange etc (deallocate immediately afterwards)
!      real(rkind), allocatable :: buf1(:,:),buf2(:,:),buf3(:)
      real(4), allocatable :: swild8(:,:),swild9(:,:) !used in ST & tracer nudging

!     Local variables
      type(llist_type),pointer :: llp
      logical :: ltmp,ltmp1,ltmp2,lexist
      character(len=48) :: inputfile

      integer :: i,j,k,l,m,mm,im2d,itmp,itmp1,itmp2,izonal5,nadv,ncor, &
                 &iwindoff,istat,icount,indx2,lq,inter_mom,ic_elev, &
                 &ipgb,isgb,iegb,irr0,nn,ifl,nd,nd1,nd2,ii,nope1, &
                 &ntmp,nrecl_et,nrecl_fl,nrecl_te,nrecl_sa,nrecl_tr,nd_gb, &
                 &jblock,jface,isd,n1,n2,n3,ndgb,ndgb1,ndgb2,irank, &
                 &ihydro_region,iabort,ie,ie2,l0,id,id1,id2,iabort_gb,j1,j2, &
                 &ne_kr,nei_kr,npp,info,num,nz_r2,ip,IHABEG,il, &
                 &ninv,it,kl,noutput_ns,iside

      real(rkind) :: cwtmp2,wtmp1,tmp,slam0,sfea0,coricoef,dfv0,dfh0, &
                     &edge_min,edge_max,tmpmin,tmpmax,tmp1,tmp2,tmp3, &
                     &tmp4,xtmp,ytmp,zdev_max,dotp2,dpmax,eta_m0,rmag,x0,x1, &
                     &x2,x3,x4,y0,y1,y2,y3,y4,ar1,ar2,ar3,ar4,fc,beta,sphi, &
                     &ubl1,ubl2,ubl3,ubl4,ubl5,ubl6,ubl7,ubl8,xn1, &
                     &xn2,yn1,yn2,xstal,ystal,ae,THAS,THAF,err_max,rr,suma, &
                     &te,sa,wx1,wx2,wy1,wy2,aux1,aux2,time,ttt, &
                     &et,qq,tr,ft1,dep,sim_year,sim_month,sim_day,sim_hour, &
                     &sim_minute,sim_second

#ifdef USE_ICM
      real(rkind) :: yday
#endif

#ifdef USE_OIL
#endif


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Executable section of SCHISM
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

!==============================================================================
!     For loose coupling, bypass initialization if not first call
!==============================================================================
!      if(first_call) then
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

#ifdef INCLUDE_TIMING
!     Start timing total and init section
      wtimer=0d0 !Zero out timers
      wtmp1=mpi_wtime()
      wtimer(0,1)=wtmp1 !total
      wtimer(1,1)=wtmp1 !init section
#endif

!     All ranks open error message and other global output files
      call parallel_rrsync(1)
      if(myrank==0) then !open as replace
        open(11,file='fort.11',status='replace') !fatal errors
      else !open as old
        open(11,file='fort.11',status='old') !fatal errors
      endif
      call parallel_rrsync(-1)

      fdb='nonfatal_0000'
      lfdb=len_trim(fdb)
      write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
      open(12,file='outputs/'//fdb,status='replace') !non-fatal errors

!...  Rank 0 writes global & local volume, energy etc data
      if(myrank==0) then
        open(91,file='total_ST.dat',status='replace')
        open(92,file='total_TR.dat',status='replace')
        open(9,file='flux.dat',status='replace')
        open(13,file='total.dat',status='replace')
!        write(13,*)ntime
        write(13,'(a200)')'Time (hours), volume, mass, potential E, kinetic E, total E, friction loss (Joule), energy leak (Joule)'
!'
      endif

!     Rank 0 open mirror file
      if(myrank==0) open(16,file='mirror.out',status='replace')

!     Echo date and time
      call date_and_time(date,timestamp)
      if(myrank==0) write(16,'(4a)')'Run begins at ',date,', ',timestamp

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!     Read from param.in; some are needed in domain decomp.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      version='description' !not really used
!     2D flag (0: 3D model; other: 2D model)
      call get_param('param.in','im2d',1,im2d,tmp,stringvalue)
      if(im2d==0) then !3D
        lm2d=.false.
      else !2D
        lm2d=.true.
        !Implicit Coriolis
        call get_param('param.in','theta2',2,itmp,theta2,stringvalue)
      endif
!     Check modules for 2D model
      if(lm2d) then
#if defined USE_ECO || defined USE_ICM || defined USE_SED || defined PREC_EVAP || defined USE_GOTM || defined USE_NAPZD || defined USE_TIMOR
        write(errmsg,*)'2D model cannot have certain modules!'
        call parallel_abort(errmsg)
#endif
      else !3D
#ifdef USE_SED2D
        write(errmsg,*)'SED2D model can only runs in 2d mode!'
        call parallel_abort(errmsg)
#endif 
      endif !lm2d


!     Some param. from WWM for schism_msgp
      call get_param('param.in','msc2',1,msc2,tmp,stringvalue)
      call get_param('param.in','mdc2',1,mdc2,tmp,stringvalue)

      call get_param('param.in','ipre',1,ipre,tmp,stringvalue)
      if(nproc>1.and.ipre/=0) call parallel_abort('ipre/=0 is not enabled for nproc>1')

      call get_param('param.in','ntracers',1,ntracers,tmp,stringvalue)
      if(ntracers<0) call parallel_abort('check ntracers')

!     Lat/lon option
      call get_param('param.in','ics',1,ics,tmp,stringvalue)
      if(ics/=1.and.ics/=2) then
        write(errmsg,*) 'Unknown ics',ics
        call parallel_abort(errmsg)
      endif

!     Some modules are not available in lon/lat mode yet
#if defined USE_SED2D || defined USE_SED || defined USE_ICM || defined USE_TIMOR
      if(ics==2) then      
        write(errmsg,*)'Some models cannot be run on lon/lat!'
        call parallel_abort(errmsg)
      endif
#endif

      call get_param('param.in','nonhydro',1,nonhydro,tmp,stringvalue)
      if(nonhydro/=0.and.nonhydro/=1) call parallel_abort('INIT: check nonhydro')
      if(nonhydro==1) then
        if(ics==2) call parallel_abort('ics=2 and nonhydro==1')
        call get_param('param.in','ihydro_region',1,ihydro_region,tmp,stringvalue)
      endif

      call get_param('param.in','indvel',1,indvel,tmp,stringvalue)
      if(indvel<-1.or.indvel>1) then
        write(errmsg,*)'Illegal indvel:',indvel
        call parallel_abort(errmsg)
      endif

!     imm: 0: without bed deformation; 1: with bed deformation (e.g., tsunami);
!     2: 3D bed deformation model (needs user coding)
!     For imm=2, user needs to manually update bottom vel. etc in update_bdef()
!     (not working yet for ics=2)
      call get_param('param.in','imm',1,imm,tmp,stringvalue)
!     For moving bed, the output is from original bottom to nvrt
      if(imm<0.or.imm>2) then
        write(errmsg,*)'Unknown imm',imm
        call parallel_abort(errmsg)
      endif
      if(imm==2.and.ics==2) call parallel_abort('imm=ics=2')
      if(imm==1) then !read in deformation at all nodes
        call get_param('param.in','ibdef',1,ibdef,tmp,stringvalue)
      endif

!     hotstart option
      call get_param('param.in','ihot',1,ihot,tmp,stringvalue)
      if(ihot<0.or.ihot>2) then
        write(errmsg,*)'Unknown ihot',ihot
        call parallel_abort(errmsg)
      endif
!#ifdef USE_SWAN
!      if(ihot==2) then
!        write(errmsg,*)'ihot cannot be 2 for wave models',ihot
!        call parallel_abort(errmsg)
!      endif
!#endif

      call get_param('param.in','ihydraulics',1,ihydraulics,tmp,stringvalue)

!...  Option for Williamson test #5 (zonal flow over an isolated mount)
      call get_param('param.in','izonal5',1,izonal5,tmp,stringvalue)
      if(izonal5/=0.and.ics==1) call parallel_abort('ics=1 and izonal5/=0')

!...  Center of projection in degrees (used for beta-plane approx.)
      call get_param('param.in','cpp_lon',2,itmp,slam0,stringvalue) !This is not really used
      call get_param('param.in','cpp_lat',2,itmp,sfea0,stringvalue)

!...  Horizontal viscosity option
!     ihorcon =0 means all hvis=0 and no hvis.gr3 is needed
      call get_param('param.in','ihorcon',1,ihorcon,tmp,stringvalue)
      if(ihorcon/=0.and.ics==2) call parallel_abort('ics=2 and ihorcon/=0')
!     Land bnd friction coefficient, needed only if ihorcon/=0
      if(ihorcon/=0) then
        call get_param('param.in','cdh',2,itmp,cdh,stringvalue)
        if(cdh<0) call parallel_abort('MAIN: cdh<0')
      endif

!...  Horizontal diffusivity option; only works for upwind/TVD
!     ihdif=0 means all hdif=0 and no hdif.gr3 is needed
      call get_param('param.in','ihdif',1,ihdif,tmp,stringvalue)
!...  Implicitness factor
      call get_param('param.in','thetai',2,itmp,thetai,stringvalue)

!...  Baroclinic flags
      call get_param('param.in','ibcc',1,ibc,tmp,stringvalue)
      call get_param('param.in','itransport',1,ibtp,tmp,stringvalue)
      if(ibc/=0.and.ibc/=1) call parallel_abort('Unknown ibcc')
      if(ibtp/=0.and.ibtp/=1) call parallel_abort('Unknown itransport')

      if(ibc==0) then
        if(myrank==0) write(16,*)'You are using baroclinic model'
        call get_param('param.in','nrampbc',1,nrampbc,tmp,stringvalue)
        if(nrampbc/=0) call get_param('param.in','drampbc',2,itmp,drampbc,stringvalue)
      else !ibc=1
        if(ibtp==0) then
          if(myrank==0) write(16,*)'Barotropic model without ST calculation'
!'
        else !ibtp=1
          if(myrank==0) write(16,*)'Barotropic model with ST calculation'
!'
        endif
      endif

!     Run time in days
      call get_param('param.in','rnday',2,itmp,rnday,stringvalue)
!...  dramp not used if nramp=0
      call get_param('param.in','nramp',1,nramp,tmp,stringvalue)
      if(nramp/=0) call get_param('param.in','dramp',2,itmp,dramp,stringvalue)

      if(nramp/=0.and.nramp/=1) then
        write(errmsg,*)'Unknown nramp',nramp
        call parallel_abort(errmsg)
      endif

!     Time step in seconds
      call get_param('param.in','dt',2,itmp,dt,stringvalue)

!...  Advection flag for momentum eq.; 1-Euler; 2: R-K
      call get_param('param.in','nadv',1,nadv,tmp,stringvalue)
      if(nadv<0.or.nadv>2) then
        write(errmsg,*)'Unknown advection flag',nadv
        call parallel_abort(errmsg)
      endif

!...  Min/max. btracking step
      call get_param('param.in','dtb_min',2,itmp,dtb_min,stringvalue)
      call get_param('param.in','dtb_max',2,itmp,dtb_max,stringvalue)
      if(dtb_min>dtb_max.or.dtb_min<=0) call parallel_abort('dtb_min>dtb_max')
!'

!...  Minimum depth allowed
      call get_param('param.in','h0',2,itmp,h0,stringvalue)
      if(h0<=0) call parallel_abort('h0 must be positive')

!...  Bottom friction (-1: 2D; 0,1 for 3D)
      call get_param('param.in','bfric',1,nchi,tmp,stringvalue)
      if(iabs(nchi)>1) call parallel_abort('INIT: unknown nchi')
      if(lm2d.and.nchi/=-1) call parallel_abort('INIT: 2D requires nchi=-1')
      if(.not.lm2d.and.nchi==-1) call parallel_abort('INIT: 3D requires nchi/=-1')
      
      if(nchi==1) then
!       dzb_min: min. bottom boundary layer thickness [m]
        call get_param('param.in','dzb_min',2,itmp,dzb_min,stringvalue)
        call get_param('param.in','dzb_decay',2,itmp,dzb_decay,stringvalue)
        if(dzb_min<=0.or.dzb_decay>0) call parallel_abort('INIT: dzb_min<=0 or dzb_decay>0')
      endif
      if(nchi==-1) then
!       Min depth used in Manning formulation
        call get_param('param.in','hmin_man',2,itmp,hmin_man,stringvalue)
        if(hmin_man<=0) call parallel_abort('INIT: hmin wrong in Manning')
      endif

!     Coriolis options (must be 1 if ics=2)
      call get_param('param.in','ncor',1,ncor,tmp,stringvalue)
      if(iabs(ncor)>1.or.ics==2.and.ncor/=1) then
        write(errmsg,*)'Unknown ncor',ncor,ics
        call parallel_abort(errmsg)
      endif
      if(ncor==-1) then !latitude
        call get_param('param.in','latitude',2,itmp,tmp,stringvalue)
        coricoef=2*omega_e*sin(tmp/180*pi)
      else if(ncor==0) then
        call get_param('param.in','coriolis',2,itmp,coricoef,stringvalue)
      endif

!     Wind (nws=3: for conservation check; otherwise same as nws=2)
      call get_param('param.in','nws',1,nws,tmp,stringvalue)
      call get_param('param.in','wtiminc',2,itmp,wtiminc,stringvalue)
      if(nws<0.or.nws>4) then
        write(errmsg,*)'Unknown nws',nws
        call parallel_abort(errmsg)
      endif
      if(nws>0.and.dt>wtiminc) then
        write(errmsg,*)'wtiminc < dt'
        call parallel_abort(errmsg)
      endif

      iwind_form=0 !init.
      if(nws>=2.and.nws<=3) then !CORIE mode; read in hgrid.ll
        !If iwind_form=0, the stress is calculated from heat exchange
        !if iwind_form=-1, use old Pond formulation
        call get_param('param.in','iwind_form',1,iwind_form,tmp,stringvalue)
        if(iwind_form<-2.or.iwind_form>0) then
          write(errmsg,*)'Unknown iwind_form',iwind_form
          call parallel_abort(errmsg)
        endif
      endif !nws

      if(nws>0) then
        call get_param('param.in','nrampwind',1,nrampwind,tmp,stringvalue)
        call get_param('param.in','drampwind',2,itmp,drampwind,stringvalue)
        call get_param('param.in','iwindoff',1,iwindoff,tmp,stringvalue)
      endif !nws

!     Heat and salt conservation flags
      call get_param('param.in','ihconsv',1,ihconsv,tmp,stringvalue)
      call get_param('param.in','isconsv',1,isconsv,tmp,stringvalue)
      if(ihconsv<0.or.ihconsv>1.or.isconsv<0.or.isconsv>1) then
        write(errmsg,*)'Unknown ihconsv or isconsv',ihconsv,isconsv
        call parallel_abort(errmsg)
      endif
      if(isconsv/=0.and.ihconsv==0) call parallel_abort('Evap/precip model must be used with heat exchnage model')
!'
      if(ihconsv/=0.and.(nws<2.or.nws>3)) call parallel_abort('Heat budge model must have nws>=2')
!'
      if(isconsv/=0) then
#ifndef PREC_EVAP
        write(errmsg,*)'Pls enable PREC_EVAP:',isconsv
        call parallel_abort(errmsg)
!       USE_SFLUX and USE_NETCDF are definitely enabled in Makefile when
!       isconsv=1
#endif
      endif

!...  Turbulence closure options
      call get_param('param.in','itur',1,itur,tmp,stringvalue)
      if(itur<-2.or.itur>4) then
        write(errmsg,*)'Unknown turbulence closure model',itur
        call parallel_abort(errmsg)
      endif
      if(itur==0) then
        call get_param('param.in','dfv0',2,itmp,dfv0,stringvalue)
        call get_param('param.in','dfh0',2,itmp,dfh0,stringvalue)
!        dfv=dfv0; dfh=dfh0
      else if(itur==2) then !read in P&P coefficients
        call get_param('param.in','h1_pp',2,itmp,h1_pp,stringvalue)
        call get_param('param.in','h2_pp',2,itmp,h2_pp,stringvalue)
        call get_param('param.in','vdmax_pp1',2,itmp,vdmax_pp1,stringvalue)
        call get_param('param.in','vdmax_pp2',2,itmp,vdmax_pp2,stringvalue)
        call get_param('param.in','vdmin_pp1',2,itmp,vdmin_pp1,stringvalue)
        call get_param('param.in','vdmin_pp2',2,itmp,vdmin_pp2,stringvalue)
        call get_param('param.in','tdmin_pp1',2,itmp,tdmin_pp1,stringvalue)
        call get_param('param.in','tdmin_pp2',2,itmp,tdmin_pp2,stringvalue)
        if(h1_pp>=h2_pp) then
          write(errmsg,*)'h1_pp >= h2_pp in P&P'
          call parallel_abort(errmsg)
        endif
        if(vdmax_pp1<vdmin_pp1.or.vdmax_pp2<vdmin_pp2) then
          write(errmsg,*)'Wrong limits in P&P:',vdmax_pp1,vdmin_pp1,vdmax_pp2,vdmin_pp2
          call parallel_abort(errmsg)
        endif
      else if(itur==3) then
!         Closure name and stability function
          call get_param('param.in','turb_met',0,itmp,tmp,mid)
          call get_param('param.in','turb_stab',0,itmp,tmp,stab)
      endif !itur
      
!     i.c. for T,S
      call get_param('param.in','icst',1,icst,tmp,stringvalue)
      if(icst/=1.and.icst/=2) then
        write(errmsg,*)'Unknown i.c. flag',icst
        call parallel_abort(errmsg)
      endif

!     Mean T,S profile
!     If ibcc_mean=1, ts.ic is needed, which is the same input needed when
!     ihot=0 and icst=2.
      call get_param('param.in','ibcc_mean',1,ibcc_mean,tmp,stringvalue)
      if(ibcc_mean/=0.and.ibcc_mean/=1) then
        write(errmsg,*)'Unknown ibcc_mean flag',ibcc_mean
        call parallel_abort(errmsg)
      endif

!     Tracers 
      call get_param('param.in','flag_model',1,flag_model,tmp,stringvalue)
      call get_param('param.in','flag_ic',1,flag_ic,tmp,stringvalue)
      call get_param('param.in','sim_year',2,itmp,sim_year,stringvalue)
      call get_param('param.in','sim_month',2,itmp,sim_month,stringvalue)
      call get_param('param.in','sim_day',2,itmp,sim_day,stringvalue)
      call get_param('param.in','sim_hour',2,itmp,sim_hour,stringvalue)
      call get_param('param.in','sim_minute',1,itmp,sim_minute,stringvalue)
      call get_param('param.in','sim_second',2,itmp,sim_second,stringvalue)
!     Pass time info to EcoSim and ICM
#ifdef USE_ECO
      year=sim_year
      month=sim_month
      day=sim_day
      hour=sim_hour
      minutes=sim_minute
      seconds=sim_second
#endif
     
      call get_param('param.in','itr_met',1,itr_met,tmp,stringvalue)
      if(itr_met<1.or.itr_met>3) then
        write(errmsg,*)'Unknown tracer method',itr_met
        call parallel_abort(errmsg)
      endif
      if(itr_met>=2) then !TVD
        call get_param('param.in','tvd_mid2',0,itmp,tmp,tvd_mid2)
        call get_param('param.in','flimiter2',0,itmp,tmp,flimiter2)
        call get_param('param.in','h_tvd',2,itmp,h_tvd,stringvalue)
      endif

      call get_param('param.in','inu_tr',1,inu_tr,tmp,stringvalue)
      call get_param('param.in','step_nu_tr',2,itmp,step_nu_tr,stringvalue)
      if(inu_tr<0.or.inu_tr>2.or.step_nu_tr<dt) then
        write(errmsg,*)'Wrong inu_tr:',inu_tr,step_nu_tr
        call parallel_abort(errmsg)
      endif

!     Global output parameters (@nodes)
!     Make sure the order of optional modules is same during actual output
!     Array for storing local indices for optional modules
!     indx_out(i,j): i is model id (SED, NAPZD etc); j=1:2 is the start and end indices of each model
      allocate(indx_out(10,2),stat=istat)
      if(istat/=0) call parallel_abort('MAIN: indx_out failure')

      noutput=27+ntracers !all Hydro and generic tracers outputs
#ifdef USE_SED
      ! depth, d50, taub, z0, qbdl(ntracers), bedfrac(ntracers)
      noutput=noutput+4+2*ntracers
      indx_out(1,1)=noutput-(3+2*ntracers)
      indx_out(1,2)=noutput
#endif
#ifdef USE_SED2D
      noutput=noutput+9
      indx_out(1,1)=noutput-8
      indx_out(1,2)=noutput
#endif
#ifdef USE_NAPZD
      noutput=noutput+2
      indx_out(2,1)=noutput-1
      indx_out(2,2)=noutput
#endif
#ifdef USE_WWM
      noutput=noutput+26 
      indx_out(3,1)=noutput-25
      indx_out(3,2)=noutput
      !Index into out_wwm() for scalar outputs below
      allocate(indx_wwm_out(24),stat=istat)
      if(istat/=0) call parallel_abort('MAIN: indx_wwm_out failure')
!      indx_wwm_out=(/1,2,3,6,7,9,10,11,12,13,18,19,20,23,26,27,28/)
      icount=0
      do i=1,28
        if(i/=7.and.i/=8.and.i/=27.and.i/=28) then
          icount=icount+1
          if(icount>24) call parallel_abort('MAIN: indx_wwm_out')
          indx_wwm_out(icount)=i
        endif
      enddo !i
#endif

      if(flag_model==0) then !age
        noutput=noutput+ntracers/2
        indx_out(4,1)=noutput-ntracers/2+1
        indx_out(4,2)=noutput
      endif !flag_model

      if(noutput>mnout) then
        write(errmsg,*)'Increase mnout in schism_glbl to',noutput
        call parallel_abort(errmsg)
      endif

      outfile(1)='elev.61'
      outfile(2)='pres.61'
      outfile(3)='airt.61'
      outfile(4)='shum.61'
      outfile(5)='srad.61'
      outfile(6)='flsu.61'
      outfile(7)='fllu.61'
      outfile(8)='radu.61'
      outfile(9)='radd.61'
      outfile(10)='flux.61'
      outfile(11)='evap.61'
      outfile(12)='prcp.61'
      outfile(13)='bdrc.61'
      outfile(14)='wind.62'
      outfile(15)='wist.62'
      outfile(16)='dahv.62'
      outfile(17)='vert.63'
      outfile(18)='temp.63'
      outfile(19)='salt.63'
      outfile(20)='conc.63'
      outfile(21)='tdff.63'
      outfile(22)='vdff.63'
      outfile(23)='kine.63'
      outfile(24)='mixl.63'
      outfile(25)='zcor.63'
      outfile(26)='qnon.63'
      outfile(27)='hvel.64'
      variable_nm(1)='surface elevation'
      variable_nm(2)='atmopheric pressure'
      variable_nm(3)='air temperature'
      variable_nm(4)='specific humidity'
      variable_nm(5)='solar radiation'
      variable_nm(6)='fluxsu'
      variable_nm(7)='fluxlu'
      variable_nm(8)='hradu'
      variable_nm(9)='hradd'
      variable_nm(10)='total flux'
      variable_nm(11)='Evaporation rate (kg/m^2/s)'
      variable_nm(12)='Precipitation rate (kg/m^2/s)'
      variable_nm(13)='Bottom drag coefficient [-]'
      variable_nm(14)='wind speed'
      variable_nm(15)='wind stress (m^2/s^2)'
      variable_nm(16)='Depth averaged horizontal velocity'
      variable_nm(17)='vertical velocity'
      variable_nm(18)='temperature in C'
      variable_nm(19)='salinity in psu'
      variable_nm(20)='density anomaly in kg/m^3'
      variable_nm(21)='eddy diffusivity in m^2/s'
      variable_nm(22)='eddy viscosity in m^2/s'
      variable_nm(23)='turbulent kinetic energy'
      variable_nm(24)='turbulent mixing length'
      variable_nm(25)='z coordinates'
      variable_nm(26)='normalized non-hydrostatic pressure'
      variable_nm(27)='horizontal velocity'

      variable_dim(1:13)='2D scalar'
      variable_dim(14:16)='2D vector'
      variable_dim(17:26)='3D scalar'
      variable_dim(27)='3D vector'
    
      do i=1,ntracers
        write(ifile_char,'(i03)')i
        ifile_char=adjustl(ifile_char)  !place blanks at end
        ifile_len=len_trim(ifile_char)  !length without trailing blanks
        outfile(27+i)='trcr_'//ifile_char(1:ifile_len)//'.63' 
        variable_nm(27+i)='Tracer #'//trim(ifile_char)
        variable_dim(27+i)='3D scalar'
      enddo !i

      indx2=27+ntracers
#ifdef USE_SED
      outfile(indx2+1)='depth.61'
      variable_nm(indx2+1)='depth in m'
      variable_dim(indx2+1)='2D scalar'

      do i=1,ntracers
         write(ifile_char,'(i03)')i
         ifile_char=adjustl(ifile_char)
         ifile_len=len_trim(ifile_char)

         outfile(indx2+1+i)='qbdl_'//ifile_char(1:ifile_len)//'.62'
         variable_nm(indx2+1+i)='Bedload #'//trim(ifile_char)
         variable_dim(indx2+1+i)='2D vector'

         outfile(indx2+1+ntracers+i)='bfrac_'//ifile_char(1:ifile_len)//'.61'
         variable_nm(indx2+1+ntracers+i)='Bedfr (top layer) #'//trim(ifile_char)
         variable_dim(indx2+1+ntracers+i)='2D scalar'
      enddo !i

      outfile(indx2+2+2*ntracers)='bedd50.61'
      variable_nm(indx2+2+2*ntracers)='median grain size (mm)'
      variable_dim(indx2+2+2*ntracers)='2D scalar'

      outfile(indx2+3+2*ntracers)='bstress.61'
      variable_nm(indx2+3+2*ntracers)='bottom shear stress (N.m-2)'
      variable_dim(indx2+3+2*ntracers)='2D scalar'

      outfile(indx2+4+2*ntracers)='brough.61'
      variable_nm(indx2+4+2*ntracers)='bottom roughness lenght z0 (m)'
      variable_dim(indx2+4+2*ntracers)='2D scalar'

      indx2=indx2+4+2*ntracers
#endif /*USE_SED*/

#ifdef USE_SED2D
      outfile(indx2+1)='depth.61'
      variable_nm(indx2+1)='depth in m'
      variable_dim(indx2+1)='2D scalar'

      outfile(indx2+2)='cdsed.61'
      variable_nm(indx2+2)='drag coefficient in sed2d'
      variable_dim(indx2+2)='2D scalar'

      outfile(indx2+3)='cflsed.61'
      variable_nm(indx2+3)='Courant number in sed2d'
      variable_dim(indx2+3)='2D scalar'

      outfile(indx2+4)='d50moy.61'
      variable_nm(indx2+4)='surface layer mean d50'
      variable_dim(indx2+4)='2D scalar'

      outfile(indx2+5)='qtot.62'
      variable_nm(indx2+5)='total transport in m2.s-1'
      variable_dim(indx2+5)='2D vector'
         
      outfile(indx2+6)='qsus.62'
      variable_nm(indx2+6)='suspended transport in m2.s-1'
      variable_dim(indx2+6)='2D vector'
         
      outfile(indx2+7)='qbdl.62'
      variable_nm(indx2+7)='bed load transport in m2.s-1'
      variable_dim(indx2+7)='2D vector'
         
      outfile(indx2+8)='dpdxy.62'
      variable_nm(indx2+8)='Bottom slope in m.m-1'
      variable_dim(indx2+8)='2D vector'

      outfile(indx2+9)='qav.62'
      variable_nm(indx2+9)='averaged total transport in m2.s-1'
      variable_dim(indx2+9)='2D vector'

      indx2=indx2+9
#endif 

#ifdef USE_NAPZD
      outfile(indx2+1)='Bbdf.63'
      variable_nm(indx2+1)='Biological neglected loss'
      variable_dim(indx2+1)='3D scalar'
      outfile(indx2+2)='totN.63'
      variable_nm(indx2+2)='Total Nitrogyn'
      variable_dim(indx2+2)='3D scalar'
      indx2=indx2+2
#endif /*USE_NAPZD*/

#ifdef USE_WWM
      do i=1,26
        write(ifile_char,'(i03)')i
        ifile_char=adjustl(ifile_char)  !place blanks at end
        ifile_len=len_trim(ifile_char)  !length without trailing blanks
        variable_nm(indx2+i)='WWM #'//trim(ifile_char)
        if(i>24) then !vectors
          outfile(indx2+i)='wwm_'//ifile_char(1:ifile_len)//'.62'
          variable_dim(indx2+i)='2D vector'
        else
          outfile(indx2+i)='wwm_'//ifile_char(1:ifile_len)//'.61'
          variable_dim(indx2+i)='2D scalar'
        endif
      enddo !i
      indx2=indx2+26
#endif /*USE_WWM*/

      if(flag_model==0) then !age
        do i=1,ntracers/2
           write(ifile_char,'(i03)')i
           ifile_char=adjustl(ifile_char)
           ifile_len=len_trim(ifile_char)

           outfile(indx2+i)='age_'//ifile_char(1:ifile_len)//'.63'
           variable_nm(indx2+i)='Age #'//trim(ifile_char)
           variable_dim(indx2+i)='3D scalar'
        enddo !i
        indx2=indx2+ntracers/2
      endif !flag_model

!     Check
      if(indx2/=noutput) then
        write(errmsg,*)'MAIN: mismatch (1):',indx2,noutput
        call parallel_abort(errmsg)
      endif

!     nspool,ihfskip: output and file spools
      call get_param('param.in','nspool',1,nspool,tmp,stringvalue)
      call get_param('param.in','ihfskip',1,ihfskip,tmp,stringvalue)
      if(nspool==0.or.ihfskip==0) call parallel_abort('Zero nspool')
      if(mod(ihfskip,nspool)/=0) call parallel_abort('ihfskip/nspool /= integer')
!'

      do i=1,noutput
        call get_param('param.in',trim(adjustl(outfile(i))),1,iof(i),tmp,stringvalue)
        if(iof(i)/=0.and.iof(i)/=1) then
          write(errmsg,*)'Unknown output option',i,iof(i),outfile(i)
          call parallel_abort(errmsg)
        endif
      enddo !i=1,noutput
!      if(iof(24)==0) then
!        if(myrank==0) write(16,*)'Reset zcor output flag'
!        iof(24)=1
!      endif
!     For 2D model reset some flags
      if(lm2d) then
        iof(3:12)=0
        iof(17:26)=0
      endif !lm2d

!...  Non-standard outputs at sides, nodes and centroids and whole and half levels
!...  These outputs share same spool and stack size as standard ones
!...  Suffix convention:
!       65: 2D side 
!       66: 2D element
!       67: 3D side and whole level
!       68: 3D side and half level
!       69: 3D element and whole level
!       70: prism centers (centroid @ half levels)

      noutput_ns=20 !hvel.67,vert.69,temp.70,salt.70 etc
      allocate(outfile_ns(noutput_ns),iof_ns(noutput_ns),stat=istat)
      outfile_ns(1)='hvel.67' !must be consistent when calling the output routine later
      outfile_ns(2)='vert.69'
      outfile_ns(3)='temp.70'
      outfile_ns(4)='salt.70'
      outfile_ns(5)='z0st.66'
      outfile_ns(6)='z0eq.66'
      outfile_ns(7)='z0cr.66'
      outfile_ns(8)='z0sw.66'
      outfile_ns(9)='z0wr.66'
      outfile_ns(10)='bpgr.65'
      outfile_ns(11)='wafo.67'
      outfile_ns(12)='bdoc.66'
      outfile_ns(13)='bnh4.66'
      outfile_ns(14)='bno3.66'
      outfile_ns(15)='bpo4.66'
      outfile_ns(16)='bcod.66'
      outfile_ns(17)='sbdo.66'
      outfile_ns(18)='sbsa.66'
      outfile_ns(19)='bthk.66'
      outfile_ns(20)='bage.66'

      !varnm_ns is not used at the moment except for noting purpose
      !varnm_ns(1)='3D horizontal vel. at sides and whole levels'
      !varnm_ns(2)='Vertical vel. at centroids and whole levels'
      !varnm_ns(3)='temperature at prism centers'
      !varnm_ns(4)='salinity at prism centers'
      !varnm_ns(5)='Sediment transport z0 at elem. center (m)'
      !varnm_ns(6)='Roughness length z0 at elem. center(m)'
      !varnm_ns(7)='Current-ripples z0 at elem. center(m)'
      !varnm_ns(8)='Sand-waves z0 at elem. center (m)'
      !varnm_ns(9)='Wave-ripples z0 at elem. center (m)'
      !varnm_ns(10)='barotropic pressure gradient force at side centers'
      !varnm_ns(11)='wave force at side centers and whole levels'
      !varnm_ns(12)='ICM: BENDOC'
      !varnm_ns(13)='ICM: SED_BENNH4'
      !varnm_ns(14)='ICM: SED_BENNO3'
      !varnm_ns(15)='ICM: BENPO4'
      !varnm_ns(16)='ICM: SED_BENCOD'
      !varnm_ns(17)='ICM: sed_BENDO'
      !varnm_ns(18)='ICM: BENSA'
      !varnm_ns(19)='SED: total bed layer thickness (m)'
      !varnm_ns(20)='SED: total bed layer age (sec)'

      do i=1,noutput_ns
        call get_param('param.in',trim(adjustl(outfile_ns(i))),1,iof_ns(i),tmp,stringvalue)
        if(iof_ns(i)/=0.and.iof_ns(i)/=1) then
          write(errmsg,*)'Unknown output_ns option',i,iof_ns(i),outfile_ns(i)
          call parallel_abort(errmsg)
        endif
      enddo !i=1,noutput_ns

!...  input information about hot start output
!...
      call get_param('param.in','hotout',1,nhot,tmp,stringvalue)
      call get_param('param.in','hotout_write',1,nhot_write,tmp,stringvalue)
      if(nhot/=0.and.nhot/=1.or.nhot*mod(nhot_write,ihfskip)/=0) then
        write(errmsg,*)'Unknown hotout or hotout_write not multiple of ihfskip',nhot,ihfskip
!'
        call parallel_abort(errmsg)
      endif

!...  JCG solver parameters
!     moitn0: output interval; mxitn0: max iterations; rtol0: relative tolerance
      call get_param('param.in','slvr_output_spool',1,moitn0,tmp,stringvalue)
      call get_param('param.in','mxitn',1,mxitn0,tmp,stringvalue)
      call get_param('param.in','tolerance',2,itmp,rtol0,stringvalue)

!...  Compute flux flag
      call get_param('param.in','consv_check',1,iflux,tmp,stringvalue)

!...  Interpolation flag for S,T and vel. in ELM
!     Kriging in vel: no bnd nodes/sides vel. use Kriging as the filter is not
!     applied there
      call get_param('param.in','inter_st',1,lq,tmp,stringvalue)
      call get_param('param.in','inter_mom',1,inter_mom,tmp,stringvalue)
      if(lq<0.or.lq>2.or.inter_mom<-1.or.inter_mom>1) then
        write(errmsg,*)'Unknown interpolation flags inter_st or inter_mom:',lq,inter_mom
!'
        call parallel_abort(errmsg)
      endif

!...  Cut-off depth for option for hgrad_nodes() near bottom (like baroc. force)
      call get_param('param.in','depth_zsigma',2,itmp,h_bcc1,stringvalue)

!...  Sponge layer for elev. & vel. (relax. factor applied to 0 elev. or uv
!     similar to T,S)
      call get_param('param.in','inu_elev',1,inu_elev,tmp,stringvalue)
      call get_param('param.in','inu_uv',1,inu_uv,tmp,stringvalue)
      if(inu_elev<0.or.inu_elev>1.or.inu_uv<0.or.inu_uv>1) then
        write(errmsg,*)'Check sponge inputs:',inu_elev,inu_uv
        call parallel_abort(errmsg)
      endif

!...  Nudging options
      call get_param('param.in','inu_st',1,inu_st,tmp,stringvalue)
      call get_param('param.in','step_nu',2,itmp,step_nu,stringvalue)
      call get_param('param.in','vnh1',2,itmp,vnh1,stringvalue)
      call get_param('param.in','vnh2',2,itmp,vnh2,stringvalue)
      call get_param('param.in','vnf1',2,itmp,vnf1,stringvalue)
      call get_param('param.in','vnf2',2,itmp,vnf2,stringvalue)
      if(inu_st<0.or.inu_st>2.or.step_nu<dt) then
        write(errmsg,*)'Check nudging inputs:',inu_st,step_nu
        call parallel_abort(errmsg)
      endif

!...  Drag formulation
      call get_param('param.in','idrag',1,idrag,tmp,stringvalue)
      if(idrag/=1.and.idrag/=2) call parallel_abort('Unknown idrag')
      if(idrag==1.and.itur>0) call parallel_abort('Linear drag requires itur<=0')
!'
      if(idrag==1.and.nchi/=0) call parallel_abort('Linear drag requires nchi=0')
!'

!...  Option to limit \hat{H} to enhance stability for large friction in shallow
!area
      call get_param('param.in','ihhat',1,ihhat,tmp,stringvalue)
      if(ihhat/=0.and.ihhat/=1) then
        write(errmsg,*)'Unknown ihhat:',ihhat
        call parallel_abort(errmsg)
      endif
      if(lm2d) ihhat=0

!...  Transport options: ELM or upwind
!     iupwind_t: 0: ELM; 1: upwind; 2: TVD (explicit); 3: TVD (implicit vertical)
!     Fix iupwind_s=iupwind_t
      call get_param('param.in','iupwind_t',1,iupwind_t,tmp,stringvalue)
!     Reset for 2D model
      if(lm2d) iupwind_t=1

      iupwind_s=iupwind_t
      if(iupwind_t<0.or.iupwind_t>3) then
        write(errmsg,*)'Unknown iupwind:',iupwind_t,iupwind_s
        call parallel_abort(errmsg)
      endif

!     tvd_mid1: model AA (my own formulation); CC (Casulli's definition of upwind
!     ratio)
!     TVD scheme will be used if itvd_e=1 and min(total depth) >=h_tvd
      if(iupwind_t>=2) then
        call get_param('param.in','tvd_mid',0,itmp,tmp,tvd_mid1)
        call get_param('param.in','flimiter',0,itmp,tmp,flimiter1)
        call get_param('param.in','h_tvd',2,itmp,h_tvd,stringvalue)
      endif

!.... Blending factor for vel. in btrack (1 for internal sides; 2 for bnd sides
!or nodes)
      call get_param('param.in','blend_internal',2,itmp,vis_coe1,stringvalue)
      call get_param('param.in','blend_bnd',2,itmp,vis_coe2,stringvalue)
      if(vis_coe1<0.or.vis_coe1>1.or.vis_coe2<0.or.vis_coe2>1) then
        write(errmsg,*)'Illegal blend_internal or blend_bnd:',vis_coe1,vis_coe2
        call parallel_abort(errmsg)
      endif

!...  Shapiro filter (used if indvel<=0)
      call get_param('param.in','shapiro',2,itmp,shapiro,stringvalue)
      if(shapiro<0.or.shapiro>0.5) then
        write(errmsg,*)'Illegal shapiro:',shapiro
        call parallel_abort(errmsg)
      endif
      !ishapiro_violation=1: geo. check; 0: no geo. check
      call get_param('param.in','ishapiro_violation',1,ishapiro_violation,tmp,stringvalue)
      if(ishapiro_violation<0.or.ishapiro_violation>1) then
        write(errmsg,*)'Illegal ishapiro_violation:',ishapiro_violation
        call parallel_abort(errmsg)
      endif

!     Kriging option
!     Choice of generalized covariance fucntion
      call get_param('param.in','kr_co',1,kr_co,tmp,stringvalue)
      if(kr_co<=0.or.kr_co>4) then
        write(errmsg,*)'Wrong kr_co:',kr_co
        call parallel_abort(errmsg)
      endif

!...  Max. for vel. magnitude
      call get_param('param.in','rmaxvel',2,itmp,rmaxvel,stringvalue)
      if(rmaxvel<1) then
        write(errmsg,*)'Illegal rmaxvel:',rmaxvel
        call parallel_abort(errmsg)
      endif
      !Add noise for btrack
      rmaxvel=rmaxvel*1.011

!...  min. vel for invoking btrack and for abnormal exit in quicksearch
      call get_param('param.in','velmin_btrack',2,itmp,velmin_btrack,stringvalue)
      if(velmin_btrack<=0) then
        write(errmsg,*)'Illegal velmin_btrack:',velmin_btrack
        call parallel_abort(errmsg)
      endif

!...  Add more noise (nudge) in init. nudging in btrack
!     to avoid underflow. This should not need to be adjusted
!     normally; may need to lower it for some benchmark tests
!     Default: btrack_nudge=1.013e-3
      call get_param('param.in','btrack_nudge',2,itmp,btrack_nudge,stringvalue)
      if(btrack_nudge<=0.or.btrack_nudge>0.1) then
        write(errmsg,*)'Illegal btrack_nudge:',btrack_nudge
        call parallel_abort(errmsg)
      endif

!     Test btrack alone (1: rotating Gausshill)
      call get_param('param.in','ibtrack_test',1,ibtrack_test,tmp,stringvalue)
      if(ibtrack_test/=0.and.ibtrack_test/=1) then
        write(errmsg,*)'Illegal ibtrack_test:',ibtrack_test
        call parallel_abort(errmsg)
      endif

!     Rouse test
      call get_param('param.in','irouse_test',1,irouse_test,tmp,stringvalue)
      if(irouse_test/=0.and.irouse_test/=1) then
        write(errmsg,*)'Illegal irouse_test:',irouse_test
        call parallel_abort(errmsg)
      endif

      if(irouse_test==1) then
#if defined USE_TIMOR || defined USE_SED
#else
        call parallel_abort('Rouse test needs USE_TIMOR or USE_SED')
#endif
        if(ntracers/=1) call parallel_abort('Rouse test requires ntracers=1')
      endif

!...  Inundation algorithm flag (1: better algorithm for fine resolution)
      call get_param('param.in','inunfl',1,inunfl,tmp,stringvalue)
      if(inunfl/=0.and.inunfl/=1) then
        write(errmsg,*)'Illegal inunfl:',inunfl
        call parallel_abort(errmsg)
      endif

!     write mode; not used really
!      call get_param('param.in','iwrite',1,iwrite,tmp,stringvalue)

!     Elev. i.c. option (elev.ic)
      call get_param('param.in','ic_elev',1,ic_elev,tmp,stringvalue)
      if(ic_elev/=0.and.ic_elev/=1) then
        write(errmsg,*)'Illegal ic_elev:',ic_elev
        call parallel_abort(errmsg)
      endif

!     Elev. b.c. ramp option (=0: ramp up from eta=0; =1: from eta2 before
!     the time loop, after the hotstart loop)
      call get_param('param.in','nramp_elev',1,nramp_elev,tmp,stringvalue)
      if(nramp_elev/=0.and.nramp_elev/=1) then
        write(errmsg,*)'Illegal nramp_elev:',nramp_elev
        call parallel_abort(errmsg)
      endif

!     Inverse barometric effects on elev. b.c.
      call get_param('param.in','inv_atm_bnd',1,inv_atm_bnd,tmp,stringvalue)
      if(inv_atm_bnd/=0.and.inv_atm_bnd/=1) then
        write(errmsg,*)'Illegal inv_atm_bnd:',inv_atm_bnd
        call parallel_abort(errmsg)
      endif
      !Reference atmos. pressure 
      call get_param('param.in','prmsl_ref',2,itmp,prmsl_ref,stringvalue)

!     Scales for dimensioning in inter-subdomain btrack
!     mxnbt=s1_mxnbt*nmm*nvrt is the dimension of btlist (nmm is the max. of all
!     nsa);
!     mnbt=max(nbt)*s2_mxnbt is the dimension of btsendq,bttmp,btdone
!       (nbt is the initial # of inter-subdomain trajectories), and
!     mnbt*nnbr is the dimension of btrecvq() in routine inter_btrack (nnbr is #
!     of nbr processes).
      call get_param('param.in','s1_mxnbt',2,itmp,s1_mxnbt,stringvalue)
      call get_param('param.in','s2_mxnbt',2,itmp,s2_mxnbt,stringvalue)
      if(s1_mxnbt<=0.or.s2_mxnbt<=0) then
        write(errmsg,*)'Illegal s[12]_mxnbt:',s1_mxnbt,s2_mxnbt
        call parallel_abort(errmsg)
      endif

!     Station output option (/=0: need station.in)
!     If ics=2, the coord. in station.in must be in lat/lon (degrees)
      call get_param('param.in','iout_sta',1,iout_sta,tmp,stringvalue)

      if(iout_sta/=0) then
        call get_param('param.in','nspool_sta',1,nspool_sta,tmp,stringvalue) !output skip
        if(nspool_sta<=0) call parallel_abort('Wrong nspool_sta')
        if(mod(nhot_write,nspool_sta)/=0) call parallel_abort('mod(nhot_write,nspool_sta)/=0')
!'
      endif

!...  Read harmonic analysis information (Adapted from ADCIRC)
      call get_param('param.in','iharind',1,iharind,tmp,stringvalue)

!...  WWM 
!     Coupling flag
!     0: decoupled so 2 models will run independently;
!     1: full coupled (elevation, vel, and wind are all passed to WWM);
!     2: 1-way coupling: only R.S. from WWM feedback to SCHISM
      call get_param('param.in','icou_elfe_wwm',1,icou_elfe_wwm,tmp,stringvalue)
      if(icou_elfe_wwm<0.or.icou_elfe_wwm>7) then
        write(errmsg,*)'Wrong coupling flag:',icou_elfe_wwm
        call parallel_abort(errmsg)
      endif

!     Coupling interval (# of time steps)
      call get_param('param.in','nstep_wwm',1,nstep_wwm,tmp,stringvalue)
      if(nstep_wwm<1) then
        write(errmsg,*)'Wrong coupling interval:',nstep_wwm
        call parallel_abort(errmsg)
      endif

!     Min (total) depth in radiation stress calculation
      call get_param('param.in','hmin_radstress',2,itmp,hmin_radstress,stringvalue)      

!     Wave boundary layer option
      call get_param('param.in','iwbl',1,iwbl,tmp,stringvalue)
      if(iwbl<0.or.iwbl>1) then
        write(errmsg,*)'Wrong iwbl:',iwbl
        call parallel_abort(errmsg)
      endif
      if(iwbl/=0.and.(nchi/=1.or.icou_elfe_wwm==0)) then
        write(errmsg,*)'WBL requires nchi=1:',iwbl,nchi,icou_elfe_wwm
        call parallel_abort(errmsg)
      endif

!     Check compatability between 2D model and parameters
      if(lm2d) then
        if(ntracers/=0.or.nonhydro==1.or.ihdif/=0.or. &
     &ibc==0.or.ibtp==1.or.nchi/=-1.or.ihconsv==1.or.isconsv==1.or. &
     &itur/=0.or.ibcc_mean/=0.or.inu_st/=0.or.icst==2.or.nws==3) then
          write(errmsg,*)'Uncompatable params. for 2D model:',ntracers, &
     &nonhydro,ihdif,ibc,ibtp,nchi,ihconsv,isconsv,itur, &
     &ibcc_mean,inu_st,icst,nws
          call parallel_abort(errmsg)
        endif
      endif !lm2d

!     Volume and mass sources/sinks option. For mass, only upwind/TVD is
!     available
      call get_param('param.in','if_source',1,if_source,tmp,stringvalue)
      if(if_source/=0.and.if_source/=1) call parallel_abort('Wrong if_source')

!     Bottom b.c. option for momentum eq. (1: old; 2: new) 
      call get_param('param.in','ibottom_bc',1,ibottom_bc,tmp,stringvalue)
      if(ibottom_bc/=1.and.ibottom_bc/=2) call parallel_abort('Wrong ibottom_bc')

!...  Check parameter read in from param.in
      if(myrank==0) write(16,*)'done reading param.in; s2_mxnbt in param.in =',s2_mxnbt
!'
!-----------------------------------------------------------------
!...  End reading param.in

!     Setup cyclic node index (used in decomp.)
      do i=1,3
        do j=1,2
          nx(i,j)=i+j
          if(nx(i,j)>3) nx(i,j)=nx(i,j)-3
          if(nx(i,j)<1.or.nx(i,j)>3) then
            write(errmsg,*)'MAIN: nx wrong',i,j,nx(i,j)
            call parallel_abort(errmsg)
          endif
        enddo !j
      enddo !i

!     Aquire vertical grid
      call aquire_vgrid

!     Partition horizontal grid into subdomains
      call partition_hgrid

!     Aquire full horizontal grid based on partition
      call aquire_hgrid(.true.) 

!     Dump horizontal grid
      call dump_hgrid

#ifdef DEBUG
!     Test if ipgl and isgl are in ascending rank order for _residents_;
!     iegl has no problem as an element can be in no more than 2 processes
      do i=1,np
        ipgb=iplg(i)
        llp=>ipgl(ipgb)%next
        j=0
        do
          if(.not.associated(llp)) exit
          j=j+1
          if(j>1.and.llp%rank<=irr0) then
            write(errmsg,*)'Node not in order:',ipgb
            call parallel_abort(errmsg)
          endif
          irr0=llp%rank
          llp=>llp%next
        enddo
      enddo !i

      do i=1,ns
        isgb=islg(i)
        llp=>isgl(isgb)%next
        j=0
        do
          if(.not.associated(llp)) exit
          j=j+1
          if(j>1.and.llp%rank<=irr0) then
            write(errmsg,*)'Side not in order:',isgb
            call parallel_abort(errmsg)
          endif
          irr0=llp%rank
          llp=>llp%next
        enddo
      enddo !i
#endif

!     Construct parallel message-passing tables
      call msgp_tables

!     Initialize parallel message-passing datatypes
      call msgp_init

!     Synchronize
      call parallel_barrier

!     Debug
!      fdb='list_0000'
!      lfdb=len_trim(fdb)
!      write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
!      open(32,file=trim(fdb),status='unknown') 
!      do i=1,np_global
!        if(.not.associated(ipgl(i)%next)) write(32,*)i 
!      enddo !i
!      close(32)

!      call parallel_finalize
!      stop

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Allocate data arrays
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!     Check geometry
      if(npa>nea.or.nea>nsa) call parallel_abort('npa>nea.or.nea>nsa')
      if(np>ne.or.ne>ns) call parallel_abort('np>ne.or.ne>ns')

!     Allocate message passing arrays
      allocate(rrqst(0:nproc-1),stat=istat)
      if(istat/=0) call parallel_abort('main: rrqst allocation failure')
      allocate(rstat(MPI_STATUS_SIZE,0:nproc-1),stat=istat)
      if(istat/=0) call parallel_abort('main: rstat allocation failure')
      allocate(srqst(0:nproc-1),stat=istat)
      if(istat/=0) call parallel_abort('main: srqst allocation failure')
      allocate(sstat(MPI_STATUS_SIZE,0:nproc-1),stat=istat)
      if(istat/=0) call parallel_abort('main: sstat allocation failure')

!     Allocate the remaining grid geometry arrays held in schism_glbl
      allocate(kbe(nea),idry_e(nea),lqk(nea),ie_kr(nea), &
     &krvel(nea),itvd_e(nea),ze(nvrt,nea),dldxy(3,2,nea),dp00(npa),kfp(npa),kbp(npa), &
     &kbp00(npa),kbp_e(np),idry(npa),hmod(npa),znl(nvrt,npa), &
     &kbs(nsa),idry_s(nsa),isidenei2(4,ns),zs(nvrt,nsa), &
     &delj(ns),ibnd_ext_int(npa),pframe(3,3,npa),sigma_lcl(nvrt,npa),stat=istat)
      if(istat/=0) call parallel_abort('MAIN: grid geometry arrays allocation failure')
!'

!     Allocate the remaining arrays held in schism_glbl, except for Kriging related arrays 
!     ntracers,ntracers2 and mntr already read in from grid_subs.F90
      allocate(tsel(2,nvrt,nea),tem0(nvrt,npa),sal0(nvrt,npa),eta1(npa),eta2(npa), &
          & we(nvrt,nea),we_fv(nvrt,nea),su2(nvrt,nsa),sv2(nvrt,nsa),ufg(nvrt,nea,3),vfg(nvrt,nea,3), &
          & tsd(nvrt,nsa),ssd(nvrt,nsa),tnd(nvrt,npa),snd(nvrt,npa), &
          & prho(nvrt,npa),q2(nvrt,npa),xl(nvrt,npa),xlmin2(npa), &
          & uu2(nvrt,npa),vv2(nvrt,npa),ww2(nvrt,npa),bdef(npa),bdef1(npa),bdef2(npa),dfh(nvrt,npa), &
          & tr_el(mntr,nvrt,nea),bdy_frc(mntr,nvrt,nea),flx_sf(mntr,nea),flx_bt(mntr,nea), &
          & xlon_el(nea),ylat_el(nea),albedo(npa),cspline_ypp_nd(nvrt,npa,2), &
          & cspline_ypp_sd(nvrt,nsa,2),stat=istat)
      if(istat/=0) call parallel_abort('MAIN: dynamical arrays allocation failure')
!'

!     Allocate boundary forcings 
      allocate(iettype(max(1,nope_global)),ifltype(max(1,nope_global)), &
          & itetype(max(1,nope_global)),isatype(max(1,nope_global)), &
          & itrtype(max(1,nope_global)),trobc(nope_global),tobc(nope_global), &
          & sobc(nope_global),vobc1(nope_global),vobc2(nope_global), &
          & eth(mnond_global,nope_global),tth(nvrt,mnond_global,nope_global),sth(nvrt,mnond_global,nope_global), &
          & qthcon(nope_global),carea(nope_global), &
          & th_dt(max(1,ntracers),nthfiles),th_time(max(1,ntracers),2,nthfiles), &
          & uth(nvrt,nsa),vth(nvrt,nsa),ath(nope_global,max(1,ntracers),2,nthfiles), &
          & ath2(max(2,ntracers),nvrt,neta_global,2,nthfiles2), &
          & uthnd(nvrt,mnond_global,nope_global),vthnd(nvrt,mnond_global,nope_global), &
          & eta_mean(npa),trth(max(1,ntracers),nvrt,mnond_global,max(1,nope_global)),stat=istat)
!           iet1lg(nope),ifl1lg(nope),ite1lg(nope),isa1lg(nope),stat=istat)
      if(istat/=0) call parallel_abort('MAIN: 1st bnd forcings allocation failure')            
!'

!     All other arrays
      allocate(ptbt(4,nvrt,npa),sdbt(4,nvrt,nsa),webt(nvrt,nea), & !bubt(2,nea), & 
         &  windx1(npa),windy1(npa),windx2(npa),windy2(npa),windx(npa),windy(npa), &
         &  tau(2,npa),iadv(npa),windfactor(npa),pr1(npa),airt1(npa),shum1(npa), &
         &  pr2(npa),airt2(npa),shum2(npa),pr(npa),sflux(npa),srad(npa),tauxz(npa),tauyz(npa), &
         &  fluxsu(npa),fluxlu(npa),hradu(npa),hradd(npa),chi(nsa),cori(nsa),Cd(nsa), &
         &  Cdp(npa),rmanning(npa),rough_p(npa),dfv(nvrt,npa),elev_nudge(npa),uv_nudge(npa), &
         &  hdif(nvrt,npa),hvis(nvrt,nea),d2uv(2,nvrt,nsa),fluxprc(npa),fluxevp(npa), & 
         &  bcc(2,nvrt,nsa),sparsem(0:(mnei+1),np), & !sparsem for non-ghosts only
         &  t_nudge(npa),s_nudge(npa),tr_nudge(npa),dr_dxy(2,nvrt,nsa), &
         &  fun_lat(0:2,npa),dav(2,npa),elevmax(npa),dav_max(2,npa),dav_maxmag(npa), &
         &  tnd_nu1(nvrt,npa),snd_nu1(nvrt,npa),tnd_nu2(nvrt,npa),snd_nu2(nvrt,npa),tnd_nu(nvrt,npa),snd_nu(nvrt,npa), &
         &  diffmax(npa),diffmin(npa),dfq1(nvrt,npa),dfq2(nvrt,npa),xlsc0(npa), & 
!          Note: swild, swild2, swild10 will be re-dimensioned (larger dimension) later
         &  nwild(nea+12),nwild2(ne_global),swild(nsa+nvrt+12+ntracers),swild2(nvrt,12),swild10(max(3,nvrt),12), &
         &  swild3(50+ntracers),swild4(nvrt,3+2*ntracers),swild8(nvrt,2),&
         &  iwater_type(npa),rho_mean(nvrt,nea),erho(nvrt,nea),& 
         &  PSQ(nea),PSK(nea),surf_t1(npa),surf_t2(npa),surf_t(npa),etaic(npa),stat=istat)
      if(istat/=0) call parallel_abort('MAIN: other allocation failure')

!     Tracers
      if(ntracers==0) then !allocate some trivial arrays (for reference only)
        allocate(trel0(1,1,1),trel(1,1,1),tr_nd(1,1,1),stat=istat)
      else
        allocate(trel0(ntracers,nvrt,nea),trel(ntracers,nvrt,nea), &
                 tr_nd(ntracers,nvrt,npa),stat=istat)

        if(inu_tr==2) then
          allocate(swild9(ntracers,nvrt),trnd_nu1(ntracers,nvrt,npa),trnd_nu2(ntracers,nvrt,npa), &
     &trnd_nu(ntracers,nvrt,npa),stat=itmp)
          if(itmp/=0) call parallel_abort('INIT: alloc failed (56)')
        endif
      endif
      if(istat/=0) call parallel_abort('MAIN: other allocation failure')

#ifdef USE_ECO
      if(ntracers/=25) call parallel_abort('MAIN: ntracer/=25 in EcoSim')
      allocate(Pair(nea), Tair(nea), Hair(nea), Uwind(nea), Vwind(nea), cloud(nea), &
               SpecIr(nea,Nbands),avcos(nea,Nbands),stat=istat)
      if(istat/=0) call parallel_abort('MAIN: EcoSim allocation failure')
#endif

#ifdef USE_NAPZD
      allocate(Bio_bdef(nvrt,nea),stat=istat)
      if(istat/=0) call parallel_abort('MAIN: NAPZD allocation failure')
#endif

!     Wave model arrays
#ifdef  USE_WWM
      allocate(wwave_force(2,nvrt,nsa),out_wwm(npa,35),out_wwm_windpar(npa,10), &
     &stokes_vel(2,nvrt,npa),jpress(npa),sbr(2,npa),sbf(2,npa),stat=istat)
      if(istat/=0) call parallel_abort('MAIN: WWM allocation failure')
      wwave_force=0; out_wwm=0; out_wwm_windpar=0
      stokes_vel=0; jpress=0; sbr=0; sbf=0 
#endif

#ifdef USE_TIMOR
!     Allocate TIMOR arrays
#endif 

#ifdef USE_ICM
      call icm_init
#endif 

!     Non-hydrostatic arrays
      allocate(qnon(nvrt,npa),stat=istat)
      if(istat/=0) call parallel_abort('MAIN: Nonhydro allocation failure (1)')
!'
      qnon=0 !initialize

      if(nonhydro==1) then
        allocate(ihydro(npa),stat=istat) 
        if(istat/=0) call parallel_abort('MAIN: Nonhydro allocation failure')
!'
      endif

!     Alloc flux output arrays
      if(iflux/=0) then
        allocate(iflux_e(nea),stat=istat)
        if(istat/=0) call parallel_abort('MAIN: iflux_e alloc')

        open(32,file='fluxflag.prop',status='old')
        do i=1,ne_global
          read(32,*)j,tmp1
          if(tmp1<-1) call parallel_abort('MAIN: fluxflag.prop has <-1')
          if(iegl(i)%rank==myrank) iflux_e(iegl(i)%id)=tmp1
        enddo
        close(32)
!        do i=1,nea !must be aug.
!          iflux_e(i)=maxval(nwild2(elnode(1:3,i)))
!        enddo !i

        itmp1=maxval(iflux_e)
        call mpi_allreduce(itmp1,max_flreg,1,itype,MPI_MAX,comm,ierr)
        if(max_flreg<=0) call parallel_abort('INIT: fluxflag.prop flag wrong')
      endif !iflux

!     Test message passing here
!      call parallel_finalize
!      stop

!     Initialize some arrays and constants
      tempmin=0; tempmax=40; saltmin=0; saltmax=42
!     temp fix
!      tempmin=0; tempmax=1005; saltmin=0; saltmax=42
      pr1=0; pr2=0; pr=prmsl_ref !uniform pressure (the const. is unimportant)
      uthnd=-99; vthnd=-99; eta_mean=-99; uth=-99; vth=-99; !flags
      fluxsu00=0; srad00=0 !for nws/=3
      elevmax=-1.e34; dav_maxmag=-1; dav_max=0
      tsel=0; trel=0
      timer_ns=0

!     for output
      airt1=0; shum1=0;  airt2=0; shum2=0; srad=0; fluxsu=0; fluxlu=0
      hradu=0; hradd=0; sflux=0; windx=0; windy=0
      q2=0; xl=0 !for hotstart with itur/=3 only
      dfq1=0; dfq2=0 !for hotstart
      fluxevp=0; fluxprc=0
!     Fort.12 flags
      ifort12=0

!...  Test node, sidecenter and centroid lat/lon conversions
!      if(ics==2) then
!        errmax=-1 !max. distance
!        do i=1,npa
!          call compute_ll(xnd(i),ynd(i),znd(i),rlon,rlat)
!          x2=rearth*cos(rlat)*cos(rlon)
!          y2=rearth*cos(rlat)*sin(rlon)
!          z2=rearth*sin(rlat)
!          dis=sqrt((x2-xnd(i))**2+(y2-ynd(i))**2+(z2-znd(i))**2)
!          write(12,*)'Node ll:',iplg(i),rlon/pi*180,rlat/pi*180,xlon(i)/pi*180,ylat(i)/pi*180, &
!     &xnd(i),ynd(i),znd(i),x2,y2,z2,dis
!          if(dis>errmax) errmax=dis
!        enddo !i
!        write(12,*)'Node max. error=',errmax
!
!        errmax=-1 !max. distance
!        do i=1,nsa
!          n1=isidenode(1,i)
!          n2=isidenode(2,i)
!          call compute_ll(xcj(i),ycj(i),zcj(i),rlon,rlat)
!          rad=sqrt(xcj(i)**2+ycj(i)**2+zcj(i)**2)
!          x2=rad*cos(rlat)*cos(rlon)
!          y2=rad*cos(rlat)*sin(rlon)
!          z2=rad*sin(rlat)
!          rlon2=(xlon(n1)+xlon(n2))/2 !has problem around 180 deg. etc
!          rlat2=(ylat(n1)+ylat(n2))/2
!          dis=sqrt((x2-xcj(i))**2+(y2-ycj(i))**2+(z2-zcj(i))**2)
!          write(12,*)'Side ll:',iplg(isidenode(1:2,i)),rlon/pi*180,rlat/pi*180,rlon2/pi*180,rlat2/pi*180, &
!     &xcj(i),ycj(i),zcj(i),x2,y2,z2,dis
!          if(dis>errmax) errmax=dis
!        enddo !i
!        write(12,*)'Side max. error=',errmax
!
!        errmax=-1 !max. distance
!        do i=1,nea
!          call compute_ll(xctr(i),yctr(i),zctr(i),rlon,rlat)
!          rad=sqrt(xctr(i)**2+yctr(i)**2+zctr(i)**2)
!          x2=rad*cos(rlat)*cos(rlon)
!          y2=rad*cos(rlat)*sin(rlon)
!          z2=rad*sin(rlat)
!          rlon2=sum(xlon(elnode(1:3,i)))/3 !has problem around 180 deg.
!          rlat2=sum(ylat(elnode(1:3,i)))/3 !has problem around 180 deg.
!          dis=sqrt((x2-xctr(i))**2+(y2-yctr(i))**2+(z2-zctr(i))**2)
!          write(12,*)'Elem. ll:',iplg(elnode(:,i)),rlon/pi*180,rlat/pi*180,rlon2/pi*180,rlat2/pi*180, &
!     &xctr(i),yctr(i),zctr(i),x2,y2,z2,dis
!          if(dis>errmax) errmax=dis
!        enddo !i
!        write(12,*)'Elem. max. error=',errmax
!
!        call parallel_finalize
!        stop
!      endif !ics

!...  Finish off some remaining geometric calcualtions
!     Output side length distribution 
      edge_max=-1 !max. side length
      edge_min=1.e25
      do i=1,ns
        if(distj(i)>edge_max) edge_max=distj(i)
        if(distj(i)<edge_min) edge_min=distj(i)
      enddo !i
      call mpi_reduce(edge_min,tmpmin,1,rtype,MPI_MIN,0,comm,ierr)
      call mpi_reduce(edge_max,tmpmax,1,rtype,MPI_MAX,0,comm,ierr)
      if(myrank==0) write(16,*)'Max. & min. sidelength= ',tmpmax,tmpmin

!...  Compute transformation tensor for node frame pframe(i,j,ip) (in global frame) for ics=2 
!...  where j is the axis id, i is the component id, ip is the local node id
!...  pframe is along local ll direction: 2nd index indicates zonal or meridional or outward radial
!...  directions. Note that the vectors are strictly undefined at 2 poles, but can be calculated 
!...  as all we need is just a frame there.
!...  For ics=1 pframe is not needed
      pframe=0 !for ics=1
      if(ics==2) then
        do i=1,npa
          pframe(1,1,i)=-sin(xlon(i)) !zonal dir.
          pframe(2,1,i)=cos(xlon(i))
          pframe(3,1,i)=0
          pframe(1,2,i)=-cos(xlon(i))*sin(ylat(i)) !meri. dir.
          pframe(2,2,i)=-sin(xlon(i))*sin(ylat(i))
          pframe(3,2,i)=cos(ylat(i))
          call cross_product(pframe(1,1,i),pframe(2,1,i),pframe(3,1,i), &
                            &pframe(1,2,i),pframe(2,2,i),pframe(3,2,i), &
                            &pframe(1,3,i),pframe(2,3,i),pframe(3,3,i))
        enddo !i=1,npa

        !Check
        zdev_max=-1 !max. deviation of zp-axis and radial direction
        do i=1,npa
          tmp=sqrt(xnd(i)**2+ynd(i)**2+znd(i)**2)
          if(tmp==0) call parallel_abort('MAIN: node radial wrong')
          swild(1)=xnd(i)/tmp
          swild(2)=ynd(i)/tmp
          swild(3)=znd(i)/tmp
          dotp2=dot_product(swild(1:3),pframe(1:3,3,i))
          if(abs(dotp2-1)>zdev_max) zdev_max=abs(dotp2-1)
          !write(12,*)'pframe:',iplg(i),pframe(:,:,i) !dotp2
        enddo !i=1,npa

        call mpi_reduce(zdev_max,tmp,1,rtype,MPI_MAX,0,comm,ierr)
        if(myrank==0) then
          write(16,*)'Max. pframe dev. from radial= ',real(tmp) !zdev_max
          call flush(16) ! flush "mirror.out"
        endif
      endif !ics==2

!...  Modified depth
      dpmax=-1.e25 !max. depth
      do i=1,npa
        if(ivcor==2) hmod(i)=min(dp(i),h_s)
        if(dp(i)>dpmax) dpmax=dp(i)
      enddo !i=1,npa
!     Save intial depth for bed deformation case
      dp00=dp

      if(ivcor==2) then; if(ztot(1)>=-dpmax) then
        write(errmsg,*)'1st z-level must be below max. depth:',dpmax
        call parallel_abort(errmsg)
      endif; endif

!...  Read in sigma coord. and kbp from vgrid.in if ivcor=1
      if(ivcor==1) then
        open(19,file='vgrid.in',status='old')
        read(19,*); read(19,*) nvrt
        do i=1,np_global
          read(19,*)j,itmp,swild(itmp:nvrt)
          if(ipgl(i)%rank==myrank) then
            id1=ipgl(i)%id
            kbp(id1)=itmp
            sigma_lcl(itmp:nvrt,id1)=swild(itmp:nvrt)
          endif
        enddo !i
        close(19)
      endif !ivcor==1

!...  Derivatives of shape functions
!...  For ics=2, this is done inside element frame
      do i=1,nea
        do j=1,3
          dldxy(j,1,i)=(yel(nx(j,1),i)-yel(nx(j,2),i))/2/area(i) !dL_j/dx
          dldxy(j,2,i)=(xel(nx(j,2),i)-xel(nx(j,1),i))/2/area(i) !dL_j/dy
        enddo !j
      enddo !i=1,nea

!Debug
!      if(ics==2) then
!        do i=1,nea
!          do j=1,3
!            swild(j)=xel(j,i)+2*yel(j,i)
!          enddo !j
!          write(12,*)'db_dx-ana=',dot_product(swild,dldxy(1:3,1,i))-1
!          write(12,*)'db_dy-ana=',dot_product(swild,dldxy(1:3,2,i))-2
!        enddo !i
!      endif
!      call parallel_finalize
!      stop      

!...  Compute delj for internal resident sides only (used only in horizontal diffusion)
      do i=1,ns !resident only
        if(isdel(2,i)==0) cycle
        delj(i)=sqrt((xctr(isdel(2,i))-xctr(isdel(1,i)))**2+(yctr(isdel(2,i))-yctr(isdel(1,i)))**2+ &
     &(zctr(isdel(2,i))-zctr(isdel(1,i)))**2)    
        if(delj(i)==0) call parallel_abort('MAIN: Element distance =0')
      enddo !i

!...  Compute lat/lon at element center for EcoSim 
#ifdef USE_ECO 
      open(32,file='hgrid.ll',status='old')
      read(32,*)
      read(32,*) !ne,np
      do i=1,np_global
        read(32,*)j,xtmp,ytmp
        if(ipgl(i)%rank==myrank) then
          xlon(ipgl(i)%id)=xtmp*pi/180
          ylat(ipgl(i)%id)=ytmp*pi/180
        endif
      enddo !i
      close(32)
      lreadll=.true.

      do i=1,nea
        !Error: won't work near dateline!!!! Try to use compute_ll
        xlon_el(i)=(xlon(elnode(1,i)) &
                   +xlon(elnode(2,i)) &
                   +xlon(elnode(3,i)))/3*180/pi !in degrees
        ylat_el(i)=(ylat(elnode(1,i)) &
                  +ylat(elnode(2,i))+ &
                  ylat(elnode(3,i)))/3*180/pi
      enddo !i
#endif

!... Read lat/lon for spectral spatial interpolation  in WWM
#ifdef USE_WWM
      inquire(file='hgrid.ll',exist=lexist)
      if(lexist) then
        open(32,file='hgrid.ll',status='old')
        read(32,*)
        read(32,*) !ne,np
        do i=1,np_global
           read(32,*)j,xtmp,ytmp
           if(ipgl(i)%rank==myrank) then
             xlon(ipgl(i)%id)=xtmp*pi/180
             ylat(ipgl(i)%id)=ytmp*pi/180
           endif
        enddo !i
        close(32)
        lreadll=.true. 
      endif 
#endif

!...  Classify interior/exterior bnd node and calculate edge angles (for WWM only)
!...  WARNING: if WWM is used, the _land_ b.c. part of hgrid.gr3 must have flags for (exterior) land (0) and
!...           island (1) bnds, and no open bnd is allowed on islands
!...  ibnd_ext_int:
!       0: interior node; 
!       1: exterior bnd (including open and land bnds); 
!      -1: interior bnd (land bnd only)
#ifdef USE_WWM
        ibnd_ext_int=0 !not on bnd by default
        do i=1,npa
          if(isbnd(1,i)/=0) ibnd_ext_int(i)=1 !weed out island nodes later
        enddo!i
!       Identify island nodes
        open(14,file='hgrid.gr3',status='old')
        rewind(14)
        do i=1,2+np_global+ne_global; read(14,*); enddo;
        read(14,*); read(14,*);
        do k=1,nope_global
          read(14,*) nn
          do i=1,nn; read(14,*); enddo;
        enddo !k
        read(14,*) !nland_global
        read(14,*) !nvel_global
        do k=1,nland_global
          read(14,*) nn,ifl
          do i=1,nn 
            read(14,*)ipgb
            if(ifl/=0.and.ipgl(ipgb)%rank==myrank) then !island
              nd=ipgl(ipgb)%id
              if(isbnd(1,nd)>0) call parallel_abort('No open bnd on islands for WWM')
!'
              ibnd_ext_int(nd)=-1
            endif !ifl
          enddo !i
        enddo !k
        close(14)

#endif /*USE_WWM*/

!-------------------------------------------------------------------------------
! Read in boundary condition and tidal info
!-------------------------------------------------------------------------------
      open(31,file='bctides.in',status='old')
      read(31,'(a48)') start_time
!...  Earth tidal potential
      read(31,*) ntip,tip_dp !cut-off depth for applying tidal potential
      if(ntip>0) then
        allocate(tamp(ntip),tnf(ntip),tfreq(ntip),jspc(ntip),tear(ntip),stat=istat)
        if(istat/=0) call parallel_abort('MAIN: allocation failure for tamp etc')
!'
        open(32,file='hgrid.ll',status='old')
        read(32,*)
        read(32,*) !ne,np
        do i=1,np_global
          read(32,*)j,xtmp,ytmp
          if(ipgl(i)%rank==myrank) then
            ii=ipgl(i)%id
            xlon(ii)=xtmp*pi/180
            ylat(ii)=ytmp*pi/180
            !Pre-compute species function to save time
            fun_lat(0,ii)=3*sin(ylat(ii))**2-1
            fun_lat(1,ii)=sin(2*ylat(ii))
            fun_lat(2,ii)=cos(ylat(ii))**2
          endif
        enddo !i
        close(32)
        lreadll=.true.
      
        do i=1,ntip
          read(31,*) !tag
          read(31,*) jspc(i),tamp(i),tfreq(i),tnf(i),tear(i)
          if(jspc(i)<0.or.jspc(i)>2) then
            write(errmsg,*)'Illegal tidal species #',jspc(i)
            call parallel_abort(errmsg)
          endif
          tear(i)=tear(i)*pi/180
        enddo !i
      endif !ntip>0

!...  Boundary forcing freqs.
!     All b.c. arrays are global
      read(31,*) nbfr

      if(nbfr>0) then
        allocate(amig(nbfr),ff(nbfr),face(nbfr),emo(nope_global,mnond_global,nbfr), &
     &efa(nope_global,mnond_global,nbfr),vmo(nope_global,nbfr),vfa(nope_global,nbfr),stat=istat)
        if(istat/=0) call parallel_abort('MAIN: allocation failure for amig etc')
!'
        do i=1,nbfr
          read(31,*) !tag
          read(31,*) amig(i),ff(i),face(i) !freq., nodal factor and earth equil.
          face(i)=face(i)*pi/180
        enddo
      endif

      read(31,*) nope1
      if(nope1/=nope_global) then
        write(errmsg,*)'Inconsistent # of open bnds',nope1,nope_global
        call parallel_abort(errmsg)
      endif

      nettype=0 !total # of type I bnds; global variable
      nfltype=0
      ntetype=0
      nsatype=0
      nettype2=0 !total # of type IV bnds (3D input)
      nfltype2=0 
      ntetype2=0
      nsatype2=0
      nnode_et=0 !total # of open bnd nodes that require elev2D.th
      nnode_fl=0 !total # of open bnd nodes that require uv3D.th
      nnode_te=0
      nnode_sa=0
      lflbc=.false. !flag to indicate existence of ifltype/=0
      do k=1,nope_global
        read(31,*) ntmp,iettype(k),ifltype(k),itetype(k),isatype(k)
        lflbc= lflbc.or.ifltype(k)/=0
        if(ntmp/=nond_global(k)) then
          write(errmsg,*)'Inconsistent # of nodes at open boundary',k,ntmp,nond_global(k)
          call parallel_abort(errmsg)
        endif

        if(iettype(k)==1) then
          nettype=nettype+1
        else if(iettype(k)==2) then
          read(31,*) eth(1,k)
        else if(iettype(k)==3) then
          do i=1,nbfr
            read(31,*)  !freq. name
            do j=1,nond_global(k)
              read(31,*) emo(k,j,i),efa(k,j,i) !amp. and phase
              efa(k,j,i)=efa(k,j,i)*pi/180
            enddo
          enddo
        else if(iettype(k)==4) then
          nettype2=nettype2+1
          nnode_et=nnode_et+nond_global(k)
        else if(iettype(k)/=0) then
          call parallel_abort('Invalid iettype')
        endif

!       For ics=2, uthnd, vthnd, uth, vth are all in lat/lon frame 
!       (even at poles), with exception for uthnd, vthnd and Flather b.c. (see below)
!       For ics=1, they are in global frame
        if(ifltype(k)==1) then
          nfltype=nfltype+1
        else if(ifltype(k)==2) then
          read(31,*) qthcon(k)
        else if(ifltype(k)==3) then
          do i=1,nbfr
            read(31,*)
            read(31,*) vmo(k,i),vfa(k,i) !uniform amp. and phase along each segment
            vfa(k,i)=vfa(k,i)*pi/180
          enddo
        else if(iabs(ifltype(k))==4) then
!          if(ics==2) call parallel_abort('ics=2 and ifltype=4')
!         For radiation b.c. eta must be specified
          if(ifltype(k)==-4) then
            if(iettype(k)==0) then
              write(errmsg,*)'vel. obc needs elev. to be specified: ',k
              call parallel_abort(errmsg)
            endif
            read(31,*) vobc1(k),vobc2(k) !nudging factors for incoming and outgoing flow
          endif
          
          nfltype2=nfltype2+1
          nnode_fl=nnode_fl+nond_global(k)
        else if(ifltype(k)==-1) then !Flather 1
          if(iettype(k)/=0) then
            write(errmsg,*)'Flather obc requires iettype=0:',k
            call parallel_abort(errmsg)
          endif
          read(31,*) !'eta_mean'
          do j=1,nond_global(k)
            ipgb=iond_global(k,j)
            read(31,*) eta_m0
            if(ipgl(ipgb)%rank==myrank) eta_mean(ipgl(ipgb)%id)=eta_m0
          enddo !j
          read(31,*) !'vn_mean' - mean normal vel.
          do j=1,nond_global(k)
            read(31,*) uthnd(1:nvrt,j,k) !used to denote normal vel. (i.e. along xs axis)
          enddo !j
!         ifltype(k)=0: zero out vertical velocity for transport in the open bnd elements
        else if(ifltype(k)/=0) then
          write(errmsg,*) 'Invalid ifltype:',ifltype(k)
          call parallel_abort(errmsg)
        endif

        tobc(k)=0 !init. for checking below
        if(itetype(k)==1) then
          ntetype=ntetype+1
          read(31,*) tobc(k) !nudging factor for inflow (no b.c. for outflow)
        else if(itetype(k)==2) then
          read(31,*) tth(1,1,k)
          read(31,*) tobc(k) !nudging factor
        else if(itetype(k)==3) then
          read(31,*) tobc(k) !nudging factor
        else if(itetype(k)==4) then
          ntetype2=ntetype2+1
          nnode_te=nnode_te+nond_global(k)
          read(31,*) tobc(k) !nudging factor
        else if(itetype(k)/=0) then
          write(errmsg,*) 'INVALID VALUE FOR ITETYPE'
          call parallel_abort(errmsg)
        endif

        if(tobc(k)<0.or.tobc(k)>1) then
          write(errmsg,*)'Temp. obc nudging factor wrong:',tobc(k),k
          call parallel_abort(errmsg)
        endif

        sobc(k)=0
        if(isatype(k)==1) then
          nsatype=nsatype+1
          read(31,*) sobc(k) !nudging factor
        else if(isatype(k)==2) then
          read(31,*) sth(1,1,k)
          read(31,*) sobc(k) !nudging factor
        else if(isatype(k)==3) then
          read(31,*) sobc(k) !nudging factor
        else if(isatype(k)==4) then
          nsatype2=nsatype2+1
          nnode_sa=nnode_sa+nond_global(k)
          read(31,*) sobc(k) !nudging factor
        else if(isatype(k)/=0) then
          write(errmsg,*) 'INVALID VALUE FOR ISATYPE'
          call parallel_abort(errmsg)
        endif

        if(sobc(k)<0.or.sobc(k)>1) then
          write(errmsg,*)'Salt. obc nudging factor wrong:',sobc(k),k
          call parallel_abort(errmsg)
        endif
      enddo !k=1,nope_global

!...  Tracer transport
      ntrtype=0 !total # of type I bnds
      ntrtype2=0 !total # of type II bnds (tr3D.th)
      nnode_tr=0 !total # of open bnd nodes that require tr3D.th
      if(ntracers>0) then
!       b.c.
        read(31,*) !nope_global
        do k=1,nope_global
          read(31,*) itrtype(k)
          trobc(k)=0 !init.
          if(itrtype(k)==2) then
            read(31,*) trth(1:ntracers,1,1,k)
            read(31,*) trobc(k) !nudging factor
          else if(itrtype(k)==1) then
            read(31,*) trobc(k) !nudging factor
            ntrtype=ntrtype+1
            do m=1,ntracers
              write(ifile_char,'(i03)')m
              ifile_char=adjustl(ifile_char); ifile_len=len_trim(ifile_char)
              inputfile='tr_'//ifile_char(1:ifile_len)//'.th'
              open(300+m,file=inputfile,status='old')
            enddo 
          else if(itrtype(k)==3) then !nudge to i.c.
            read(31,*) trobc(k) !nudging factor
          else if(itrtype(k)==4) then 
            read(31,*) trobc(k) !nudging factor
            ntrtype2=ntrtype2+1
            nnode_tr=nnode_tr+nond_global(k)
          else if(itrtype(k)/=0) then
            write(errmsg,*)'Wrong itrtype:',k,itrtype(k)
            call parallel_abort(errmsg)
          endif

          if(trobc(k)<0.or.trobc(k)>1) then
            write(errmsg,*)'Tr. obc nudging factor wrong:',trobc(k),k
            call parallel_abort(errmsg)
          endif
        enddo !k
      endif !ntracers

!...  Done with bctides.in
      close(31)

!...  Read 1st 2 lines of ASCII .th
      if(nettype>0) then
        open(50,file='elev.th',status='old')
        read(50,*) tmp,ath(1:nettype,1,1,1)
        read(50,*) th_dt(1,1),ath(1:nettype,1,2,1)
        if(abs(tmp)>1.e-6.or.th_dt(1,1)<dt) call parallel_abort('SCHISM_INIT: elev.th start time wrong')
        th_time(1,1,1)=0
        th_time(1,2,1)=th_dt(1,1)
      endif !nettype

      if(nfltype>0) then
        open(51,file='flux.th',status='old')
        read(51,*) tmp,ath(1:nfltype,1,1,2)
        read(51,*) th_dt(1,2),ath(1:nfltype,1,2,2)
        if(abs(tmp)>1.e-6.or.th_dt(1,2)<dt) call parallel_abort('SCHISM_INIT: flux.th start time wrong')
        th_time(1,1,2)=0
        th_time(1,2,2)=th_dt(1,2)
      endif !nfltype

      if(ntetype>0) then
        open(52,file='temp.th',status='old')
        read(52,*) tmp,ath(1:ntetype,1,1,3)
        read(52,*) th_dt(1,3),ath(1:ntetype,1,2,3)
        if(abs(tmp)>1.e-6.or.th_dt(1,3)<dt) call parallel_abort('SCHISM_INIT: temp.th start time wrong')
        th_time(1,1,3)=0
        th_time(1,2,3)=th_dt(1,3)
      endif !ntetype

      if(nsatype>0) then
        open(53,file='salt.th',status='old')
        read(53,*) tmp,ath(1:nsatype,1,1,4)
        read(53,*) th_dt(1,4),ath(1:nsatype,1,2,4)
        if(abs(tmp)>1.e-6.or.th_dt(1,4)<dt) call parallel_abort('SCHISM_INIT: salt.th start time wrong')
        th_time(1,1,4)=0
        th_time(1,2,4)=th_dt(1,4)
      endif !nsatype

      if(ntrtype>0) then !type I
        do m=1,ntracers
          read(300+m,*)tmp,ath(1:ntrtype,m,1,5)
          read(300+m,*)th_dt(m,5),ath(1:ntrtype,m,2,5)
          if(abs(tmp)>1.e-6.or.th_dt(m,5)<dt) call parallel_abort('SCHISM_INIT: htr_.th start time wrong')
          th_time(m,1,5)=0
          th_time(m,2,5)=th_dt(m,5)
        enddo !m
      endif !ntrtype

!     Check dimension of ath2
      if(max(nnode_et,nnode_fl,nnode_te,nnode_sa,nnode_tr)>neta_global) then
        write(errmsg,*) 'MAIN: impossible! Dimension overflow for ath2:',nnode_et,nnode_fl,nnode_te,nnode_sa,nnode_tr
        call parallel_abort(errmsg)
      endif
!     Binary record length for *3D.th at each time step
      nrecl_et=nbyte*(1+nnode_et) !single precision
      nrecl_fl=nbyte*(1+nnode_fl*2*nvrt)
      nrecl_te=nbyte*(1+nnode_te*nvrt)
      nrecl_sa=nbyte*(1+nnode_sa*nvrt)
      nrecl_tr=nbyte*(1+nnode_tr*nvrt*ntracers)
      if(nettype2/=0) then
        open(54,file='elev2D.th',access='direct',recl=nrecl_et,status='old')
        read(54,rec=1) floatout,ath2(1,1,1:nnode_et,1,1)
        read(54,rec=2) floatout2,ath2(1,1,1:nnode_et,2,1)
        if(abs(floatout)>1.e-6.or.floatout2<dt) call parallel_abort('SCHISM_INIT: elev2D.th start wrong')
        th_dt2(1)=floatout2
        th_time2(1,1)=0
        th_time2(2,1)=th_dt2(1)
      endif !nettype2

      if(nfltype2/=0) then
        open(58,file='uv3D.th',access='direct',recl=nrecl_fl,status='old')
        read(58,rec=1) floatout,ath2(1:2,1:nvrt,1:nnode_fl,1,2)
        read(58,rec=2) floatout2,ath2(1:2,1:nvrt,1:nnode_fl,2,2)
        if(abs(floatout)>1.e-6.or.floatout2<dt) call parallel_abort('SCHISM_INIT: uv3D.th start wrong')
        th_dt2(2)=floatout2
        th_time2(1,2)=0
        th_time2(2,2)=th_dt2(2)
      endif !nfltype2

      if(ntetype2/=0) then
        open(56,file='temp3D.th',access='direct',recl=nrecl_te,status='old')
        read(56,rec=1) floatout,ath2(1,1:nvrt,1:nnode_te,1,3)
        read(56,rec=2) floatout2,ath2(1,1:nvrt,1:nnode_te,2,3)
        if(abs(floatout)>1.e-6.or.floatout2<dt) call parallel_abort('SCHISM_INIT: temp3D.th start wrong')
        th_dt2(3)=floatout2
        th_time2(1,3)=0
        th_time2(2,3)=th_dt2(3)
      endif !ntetype2

      if(nsatype2/=0) then
        open(57,file='salt3D.th',access='direct',recl=nrecl_sa,status='old')
        read(57,rec=1) floatout,ath2(1,1:nvrt,1:nnode_sa,1,4)
        read(57,rec=2) floatout2,ath2(1,1:nvrt,1:nnode_sa,2,4)
        if(abs(floatout)>1.e-6.or.floatout2<dt) call parallel_abort('SCHISM_INIT: salt3D.th start wrong')
        th_dt2(4)=floatout2
        th_time2(1,4)=0
        th_time2(2,4)=th_dt2(4)
      endif !nsatype2
  
      if(ntrtype2/=0) then
        open(59,file='tr3D.th',access='direct',recl=nrecl_tr,status='old')
        read(59,rec=1) floatout,ath2(1:ntracers,1:nvrt,1:nnode_tr,1,5)
        read(59,rec=2) floatout2,ath2(1:ntracers,1:nvrt,1:nnode_tr,2,5)
        if(abs(floatout)>1.e-6.or.floatout2<dt) call parallel_abort('SCHISM_INIT: tr3D.th start wrong')
        th_dt2(5)=floatout2
        th_time2(1,5)=0
        th_time2(2,5)=th_dt2(5)
      endif !ntrtype2

      !Update record #
      irec_th(1:5)=2

!...  Read in hydraulics.in
      if(ihydraulics/=0) then
        !Non-hydro model not working yet
        if(nonhydro/=0) call parallel_abort('MAIN: Non-hydro model cannot be used with hydraulics option')
!'
        !Specify blocks for hydraulic transfer structures (where fluxes are specified,
        !and tracers are conserved)
        call load_structures('hydraulics.in')
        if(block_nudge<0.or.block_nudge>=1) call parallel_abort('MAIN: wrong block_nudge')
!'
        if(nhtblocks>0) then
          allocate(isblock_nd(2,npa), &
                   dir_block(3,nhtblocks), &
                   isblock_sd(2,nsa), &
                   isblock_el(nea),q_block(nhtblocks),vnth_block(2,nhtblocks), &
                   q_block_lcl(nhtblocks),iq_block_lcl(nhtblocks), &
                   iq_block(nhtblocks),stat=istat)
          if(istat/=0) call parallel_abort('MAIN: Alloc failed (9)')

          isblock_nd=0 !init.

          ! Convert global nodes to local data structure
          do i=1,nhtblocks
            do j=1,structures(i)%npair
              do k=1,2
                nd_gb=structures(i)%node_pairs(k,j)
                if(ipgl(nd_gb)%rank==myrank) then
                  structures(i)%is_local = .true.  ! tell the global structure it is local
                  nd=ipgl(nd_gb)%id
                  isblock_nd(1,nd)=i
                  isblock_nd(2,nd)=k !face #
                  !Check
                  !write(12,*)'Block nodes:',i,k,nd_gb,nd
                endif
              enddo !k
            enddo !j
          enddo !i

          !Compute block elements and (mean) face directions for each block (assumed
          !to be outer normal of face "2" from the block element
          allocate(swild99(3,nhtblocks),ibuf1(nhtblocks,1),ibuf2(nhtblocks,1),stat=istat)
          if(istat/=0) call parallel_abort('MAIN: failed to alloc. (15)')
          swild99=0 !local copy of dir
          ibuf1=0 !local counter
          do i=1,ne
            jblock=minval(isblock_nd(1,elnode(1:3,i)))
            if(jblock>0) then
              !Look for face '2' side
              do j=1,3 !side
                isd=elside(j,i)
                n1=isidenode(1,isd); n2=isidenode(2,isd)
                if(isblock_nd(1,n1)/=isblock_nd(1,n2)) then
                  write(errmsg,*)'Check block nodes (0):',iplg(isidenode(1:2,isd))
                  call parallel_abort(errmsg)
                endif
                if(isblock_nd(2,n1)==2.and.isblock_nd(2,n2)==2) then
                  swild99(1:3,jblock)=swild99(1:3,jblock)+sframe(1:3,1,isd)*ssign(j,i)
                  ibuf1(jblock,1)=ibuf1(jblock,1)+1
                endif
              enddo !j
            endif !jblock>0
          enddo !i=1,ne

#ifdef INCLUDE_TIMING
          wtmp1=mpi_wtime()
#endif
          call mpi_allreduce(swild99,dir_block,3*nhtblocks,rtype,MPI_SUM,comm,ierr)
          call mpi_allreduce(ibuf1,ibuf2,nhtblocks,itype,MPI_SUM,comm,ierr)
#ifdef INCLUDE_TIMING
          wtimer(3,2)=wtimer(3,2)+mpi_wtime()-wtmp1
#endif

          do i=1,nhtblocks
            if(ibuf2(i,1)==0) then
              write(errmsg,*) 'MAIN: orphaned block face:',i
              call parallel_abort(errmsg)
            else
              dir_block(1:3,i)=dir_block(1:3,i)/ibuf2(i,1)
              !Re-normalize dir. vector
              rmag=sqrt(dir_block(1,i)**2+dir_block(2,i)**2+dir_block(3,i)**2)
              if(rmag==0) call parallel_abort('MAIN: 0 dir vector')
              dir_block(1:3,i)=dir_block(1:3,i)/rmag

              !Check
              !write(12,*)'Block face dir:',i,ibuf2(i,1),dir_block(1:3,i)
            endif
          enddo !i
          deallocate(swild99,ibuf1,ibuf2)

          !open(49,file='flux_blocks.th',status='old')
          !Build message passing arrays for 2 ref. nodes
          !The proc that has ref. node 1 will calculate the flux first before broadcasting
          allocate(nhtrecv1(0:nproc-1),nhtsend1(0:nproc-1), &
     &ihtrecv1_ndgb(nhtblocks,0:nproc-1),ihtsend1_ndgb(nhtblocks,0:nproc-1), &
     &ihtrecv1(nhtblocks,0:nproc-1),ihtsend1(nhtblocks,0:nproc-1), &
     &block_refnd2_eta(nhtblocks),recv_bsize(nhtblocks),send_bsize(nhtblocks), &
     &htrecv_type(0:nproc-1),htsend_type(0:nproc-1),stat=istat)
          if(istat/=0) call parallel_abort('MAIN: Failed to alloc (17)')
          nhtrecv1=0 !# of recvs
          do i=1,nhtblocks
            ndgb1=structures(i)%upnode
            if(ipgl(ndgb1)%rank==myrank) then; if(ipgl(ndgb1)%id<=np) then
              ndgb2=structures(i)%downnode
              irank=ipgl(ndgb2)%rank
              if(irank/=myrank) then
                itmp=nhtrecv1(irank)+1
                if(itmp>nhtblocks) call parallel_abort('MAIN: overflow (9)')
!'
                nhtrecv1(irank)=itmp
                ihtrecv1_ndgb(itmp,irank)=ndgb2 !global node #
                ihtrecv1(itmp,irank)=i-1 !displacement into recv arrays like block_refnd2_*
              endif !ref. node 2 is outisde myrank
            endif; endif !ipgl; ref. node 1 is in myrank and not ghost
          enddo !i=1,nhtblocks 

          !Comm. send counts
          call mpi_alltoall(nhtrecv1,1,itype,nhtsend1,1,itype,comm,ierr)
          if(ierr/=MPI_SUCCESS) call parallel_abort('MAIN: all2all(2)')

          !Get global node # for sends
          do i=0,nproc-1
            if(nhtrecv1(i)/=0) then
              if(i==myrank) call parallel_abort('MAIN: illegal comm.(1)')
              call mpi_isend(ihtrecv1_ndgb(1,i),nhtrecv1(i),itype,i,600,comm,srqst(i),ierr)
              if(ierr/=MPI_SUCCESS) call parallel_abort('MAIN: send error (0)')
!'
            else
              srqst(i)=MPI_REQUEST_NULL
            endif !nhtrecv1
          enddo !i

          do i=0,nproc-1
            if(nhtsend1(i)>nhtblocks) call parallel_abort('MAIN: nhtsend1(i)>nhtblocks')
!'
            if(nhtsend1(i)/=0) then
              if(i==myrank) call parallel_abort('MAIN: illegal comm.(2)')
              call mpi_irecv(ihtsend1_ndgb(1,i),nhtsend1(i),itype,i,600,comm,rrqst(i),ierr)
              if(ierr/=MPI_SUCCESS) call parallel_abort('MAIN: recv error (0)')
!'
            else
              rrqst(i)=MPI_REQUEST_NULL
            endif !nhtsend1
          enddo !i

          call mpi_waitall(nproc,rrqst,rstat,ierr)
          if(ierr/=MPI_SUCCESS) call parallel_abort('MAIN: mpi_waitall rrqst tag=600',ierr)
          call mpi_waitall(nproc,srqst,sstat,ierr)
          if(ierr/=MPI_SUCCESS) call parallel_abort('MAIN: mpi_waitall srqst tag=600',ierr)
!'

          !Build send list
          do i=0,nproc-1
            do j=1,nhtsend1(i) !nhtsend1>0
              ndgb=ihtsend1_ndgb(j,i)
              if(ipgl(ndgb)%rank/=myrank) call parallel_abort('MAIN: send node not mine')
!'              ihtsend1_nd(j,i)=ipgl(ndgb)%id !local index
              ihtsend1(j,i)=ipgl(ndgb)%id-1 !displacement (into elev. array etc) (local index)
            enddo !j
          enddo !i

          !Create data type for exchange
          !Send type
          send_bsize=1; recv_bsize=1
          do i=0,nproc-1
            if(nhtsend1(i)/=0) then
#if MPIVERSION==1
              call mpi_type_indexed(nhtsend1(i),send_bsize,ihtsend1(1,i),rtype, &
     &htsend_type(i),ierr)
#elif MPIVERSION==2
              call mpi_type_create_indexed_block(nhtsend1(i),1,ihtsend1(1,i),rtype, &
     &htsend_type(i),ierr)
#endif
              if(ierr/=MPI_SUCCESS) call parallel_abort('MAIN: create htsend_type',ierr)
              call mpi_type_commit(htsend_type(i),ierr)
              if(ierr/=MPI_SUCCESS) call parallel_abort('commit htsend_type',ierr)
!'
              !Debug
              !write(12,*)'htsend list:',i,nhtsend1(i),ihtsend1(1:nhtsend1(i),i)
            endif
          enddo !i

          !Recv type
          do i=0,nproc-1
            if(nhtrecv1(i)/=0) then
#if MPIVERSION==1
              call mpi_type_indexed(nhtrecv1(i),recv_bsize,ihtrecv1(1,i),rtype, &
     &htrecv_type(i),ierr)
#elif MPIVERSION==2
              call mpi_type_create_indexed_block(nhtrecv1(i),1,ihtrecv1(1,i),rtype, &
     &htrecv_type(i),ierr)
#endif
              if(ierr/=MPI_SUCCESS) call parallel_abort('MAIN: create htrecv_type',ierr)
              call mpi_type_commit(htrecv_type(i),ierr)
              if(ierr/=MPI_SUCCESS) call parallel_abort('commit htrecv_type',ierr)
!'
              !Debug
              !write(12,*)'htrecv list:',i,nhtrecv1(i),ihtrecv1(1:nhtrecv1(i),i)
            endif
          enddo !i

        endif !nhtblocks>0
        !Done with hydraulics.in
        close(31)

      endif !ihydraulics/=0

!     Read in source_sink.in and open t.h. files
      if(if_source==1) then
        open(31,file='source_sink.in',status='old')
        read(31,*)nsources
        allocate(ieg_source(max(1,nsources)),stat=istat)
        if(istat/=0) call parallel_abort('INIT: ieg_source failure')
        do i=1,nsources
          read(31,*)ieg_source(i) !global elem. #
        enddo !i

        read(31,*) !blank line
        read(31,*)nsinks
        allocate(ieg_sink(max(1,nsinks)),ath3(max(1,nsources,nsinks),2+ntracers,2,nthfiles3),stat=istat)
        if(istat/=0) call parallel_abort('INIT: ieg_sink failure')
        do i=1,nsinks
          read(31,*)ieg_sink(i)
        enddo !i
        close(31)

        if(nsources>0) then !read first 2 lines
          open(63,file='vsource.th',status='old') !values (>=0) in m^3/s
          read(63,*)tmp,ath3(1:nsources,1,1,1)
          read(63,*)th_dt3(1),ath3(1:nsources,1,2,1)
          if(abs(tmp)>1.e-6.or.th_dt3(1)<dt) call parallel_abort('SCHISM_INIT: vsource.th start time wrong')
          th_time3(1,1)=0
          th_time3(2,1)=th_dt3(1)

          !msource.th: values in concentration dimension (psu etc). At
          !each step, 2+ntracers lines are read (in order: T,S,tracers)
          open(65,file='msource.th',status='old') 
          !do j=1,2+ntracers
          read(65,*)tmp,ath3(1:nsources,1:2+ntracers,1,3)
          !enddo !j
          !do j=1,2+ntracers
          read(65,*)th_dt3(3),ath3(1:nsources,1:2+ntracers,2,3)
          !enddo !j
          if(abs(tmp)>1.e-6.or.th_dt3(3)<dt) call parallel_abort('SCHISM_INIT: msource.th start time wrong')
          th_time3(1,3)=0
          th_time3(2,3)=th_dt3(3)
        endif !nsources

        if(nsinks>0) then
          open(64,file='vsink.th',status='old') !values (<=0) in m^3/s
          read(64,*)tmp,ath3(1:nsinks,1,1,2)
          read(64,*)th_dt3(2),ath3(1:nsinks,1,2,2)
          if(abs(tmp)>1.e-6.or.th_dt3(2)<dt) call parallel_abort('SCHISM_INIT: vsink.th start time wrong')
          th_time3(1,2)=0
          th_time3(2,2)=th_dt3(2)
        endif !nsinks
      endif !if_source

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!  Initialize model
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!...  compute total number of time steps 
      ntime=rnday*86400.d0/dt+0.5
      nrec=min(ntime,ihfskip)/nspool

!...  Option for specifying hydrostatic region for non-hydrostatic model
      if(nonhydro==1) then
        if(ihydro_region==1) then
          open(32,file='hydro_region.gr3',status='old')
          read(32,*)
          read(32,*) itmp1,itmp2
          if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check hydro_region.gr3')
          do i=1,np_global
            read(32,*)j,xtmp,ytmp,tmp 
            if(ipgl(i)%rank==myrank) then
              ihydro(ipgl(i)%id)=nint(tmp)
              if(nint(tmp)/=0.and.nint(tmp)/=1) call parallel_abort('MAIN: check hydro_region.gr3')
!'
            endif
          enddo !i
          close(32)
        else
          ihydro=0 !0: non-hydro node; 1: hydrostatic node
        endif
      endif !nonhydro

!...  Compute neighborhood for internal sides for Shapiro filter
!...  isidenei2(4,ns): 4 neighboring sides of a _resident_ side
!...  Info for resident sides only!
      !Flag for checking violation
      ishapiro=0 
      if(indvel<=0) ishapiro=1
#ifdef USE_SED2D
      ishapiro=1
#endif

      if(ishapiro==1) then
        fdb='rogue_0000'
        lfdb=len_trim(fdb)
        write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
        open(10,file='outputs/'//fdb,status='replace')

        iabort=0 !abort flag
        loop18: do i=1,ns !resident sides only
          if(isdel(2,i)==0) cycle loop18

!         Internal sides
          do j=1,2
            ie=isdel(j,i)
            l0=lindex_s(i,ie)
            if(l0==0) then
              write(errmsg,*)'Cannot find a side'
              call parallel_abort(errmsg)
            endif
            nwild(2*j-1)=elside(nx(l0,1),ie)
            nwild(2*j)=elside(nx(l0,2),ie)
          enddo !j=1,2
          isidenei2(1:4,i)=nwild(1:4) !local index

          if(ishapiro_violation==0) cycle loop18

!         Check if pt "0" is inside
          if(ics==1) then
            x0=xcj(i); y0=ycj(i)
            x1=xcj(isidenei2(1,i)); y1=ycj(isidenei2(1,i))
            x2=xcj(isidenei2(2,i)); y2=ycj(isidenei2(2,i))
            x3=xcj(isidenei2(3,i)); y3=ycj(isidenei2(3,i))
            x4=xcj(isidenei2(4,i)); y4=ycj(isidenei2(4,i))
          else !ics=2
            x0=0; y0=0
            do j=1,4
              isd=isidenei2(j,i)
              call project_pt('g2l',xcj(isd),ycj(isd),zcj(isd),(/xcj(i),ycj(i),zcj(i)/),&
     &sframe(:,:,i),swild2(1,j),swild2(2,j),tmp)
            enddo !j
            x1=swild2(1,1); y1=swild2(2,1)
            x2=swild2(1,2); y2=swild2(2,2)
            x3=swild2(1,3); y3=swild2(2,3)
            x4=swild2(1,4); y4=swild2(2,4)
          endif !ics

          ar4=signa(x0,x4,x1,y0,y4,y1)
          ar3=signa(x0,x2,x3,y0,y2,y3)
          if(ar3<=0.and.ar4<=0) then
            write(errmsg,*)'Degenerate parallelogram'
            call parallel_abort(errmsg)
          endif

!         Enlarge stencil if pt 0 is outside
          if(ar3<=0.or.ar4<=0) then
            if(ar3<=0) then
              if(isdel(2,isidenei2(2,i))==0.or.isdel(2,isidenei2(3,i))==0) then
                iabort=1
                write(10,*)iplg(isidenode(1:2,i)),', bnd side (3)'
                cycle loop18
              endif

              nwild(1)=2; nwild(2)=3
              do k=1,2
                id=isidenei2(nwild(k),i)
                ie2=isdel(1,id)+isdel(2,id)-isdel(k,i) !inside aug. domain
                if(isdel(1,id)<=0.or.isdel(2,id)<0.or.ie2/=isdel(1,id).and.ie2/=isdel(2,id)) then
                  write(errmsg,*)'Filter sides out of aug. domain (1):',iplg(isidenode(1:2,i)),ie2,isdel(1:2,id)
                  call parallel_abort(errmsg)
                endif
                l0=lindex_s(id,ie2)
                if(l0==0) then
                  write(errmsg,*)'Cannot find a side (9):',k
                  call parallel_abort(errmsg)
                endif
                isidenei2(nwild(k),i)=elside(nx(l0,3-k),ie2)
              enddo !k
            endif !ar3
         
            if(ar4<=0) then
              if(isdel(2,isidenei2(1,i))==0.or.isdel(2,isidenei2(4,i))==0) then
                iabort=1
                write(10,*)iplg(isidenode(1:2,i)),', bnd side (4)'
                cycle loop18
              endif

              nwild(1)=1; nwild(2)=4
              do k=1,2
                id=isidenei2(nwild(k),i)
                ie2=isdel(1,id)+isdel(2,id)-isdel(k,i)
                if(isdel(1,id)<=0.or.isdel(2,id)<0.or.ie2/=isdel(1,id).and.ie2/=isdel(2,id)) then
                  write(errmsg,*)'Filter sides out of aug. domain (2):',iplg(isidenode(1:2,i)),ie2,isdel(1:2,id)
                  call parallel_abort(errmsg)
                endif
                l0=lindex_s(id,ie2)
                if(l0==0) then
                  write(errmsg,*)'Cannot find a side (8):',k
                  call parallel_abort(errmsg)
                endif
                isidenei2(nwild(k),i)=elside(nx(l0,k),ie2)
              enddo !k
            endif !ar4
         
!           Check convexity of quad 1-4
            if(ics==1) then
              x0=xcj(i); y0=ycj(i)
              x1=xcj(isidenei2(1,i)); y1=ycj(isidenei2(1,i))
              x2=xcj(isidenei2(2,i)); y2=ycj(isidenei2(2,i))
              x3=xcj(isidenei2(3,i)); y3=ycj(isidenei2(3,i))
              x4=xcj(isidenei2(4,i)); y4=ycj(isidenei2(4,i))
            else !ics=2
              x0=0; y0=0
              do j=1,4
                isd=isidenei2(j,i)
                call project_pt('g2l',xcj(isd),ycj(isd),zcj(isd),(/xcj(i),ycj(i),zcj(i)/),&
     &sframe(:,:,i),swild2(1,j),swild2(2,j),tmp)
              enddo !j
              x1=swild2(1,1); y1=swild2(2,1)
              x2=swild2(1,2); y2=swild2(2,2)
              x3=swild2(1,3); y3=swild2(2,3)
              x4=swild2(1,4); y4=swild2(2,4)
            endif !ics

            ar1=signa(x1,x2,x3,y1,y2,y3)
            ar2=signa(x1,x3,x4,y1,y3,y4)
            ar3=signa(x1,x2,x4,y1,y2,y4)
            ar4=signa(x2,x3,x4,y2,y3,y4)
            if(ar1<=0.or.ar2<=0.or.ar3<=0.or.ar4<=0) then
              iabort=1
              write(10,*)iplg(isidenode(1:2,i)),'  Concave quad '
!              write(10,*)((isidenode(mm,isidenei2(m,i)),mm=1,2),m=1,4)
              write(10,*)ar1,ar2,ar3,ar4
              write(10,*)'--------------------------------------------'
              cycle loop18
            endif

!           Check if pt "0" is inside
            ar1=signa(x1,x2,x0,y1,y2,y0)
            ar2=signa(x2,x3,x0,y2,y3,y0)
            ar3=signa(x3,x4,x0,y3,y4,y0)
            ar4=signa(x4,x1,x0,y4,y1,y0)
            if(ar1<=0.or.ar2<=0.or.ar3<=0.or.ar4<=0) then
              iabort=1
              write(10,*)iplg(isidenode(1:2,i)),'  pt outside quad '
              write(10,*)ar1,ar2,ar3,ar4
!              write(10,*)((isidenode(mm,isidenei2(m,i)),mm=1,2),m=1,4)
              write(10,*)'----------------------------------------'
              cycle loop18
            endif
          endif !pt 0 outside
        end do loop18 !i=1,ns

        close(10)
        call mpi_allreduce(iabort,iabort_gb,1,itype,MPI_SUM,comm,ierr)
        if(iabort_gb>0) then
          write(errmsg,*)'Check rogue_* for problems in side neighborhood'
!'
          call parallel_abort(errmsg)
        endif
      endif !ishapiro

!     End of pre-processing
      if(ipre/=0) then !nproc=1
        write(*,*)'Pre-processing completed successfully!'
        call parallel_finalize
        stop
      endif

!     imm: 0: without bed deformation; 1: with bed deformation (e.g., tsunami);
!     2: 3D bed deformation model (needs user coding)
!     For imm=2, user needs to manually update bottom vel. etc in update_bdef()
!     (not working yet for ics=2)

!     Initialize variables used in tsunami model (but bdef[1,2] and ibdef are available for all models)
      bdef=0 !total deformation
      ibdef=1 !# of time steps for deformation (deformation rate=0 when it>ibdef)
      if(imm==1) then !read in deformation at all nodes
        open(32,file='bdef.gr3',status='old') !connectivity part not used
        read(32,*)
        read(32,*) itmp1,itmp2
        if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check bdef.gr3')
        do i=1,np_global
          read(32,*)j,xtmp,ytmp,tmp !total deformation
          if(ipgl(i)%rank==myrank) bdef(ipgl(i)%id)=tmp
        enddo !i
        close(32)
      endif

!...  Center of projection in degrees (used for beta-plane approx.)
      slam0=slam0*pi/180
      sfea0=sfea0*pi/180

!...  Horizontal viscosity option
!     ihorcon =0 means all hvis=0 and no hvis.gr3 is needed
      if(ihorcon/=0) then
!        if(ics==2) call parallel_abort('ics=2 and ihorcon/=0')
        open(32,file='hvis.gr3',status='old')
        read(32,*)
        read(32,*) itmp1,itmp2
        if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check hvis.gr3')
        do i=1,np_global
          read(32,*)j,xtmp,ytmp,tmp 
          if(ipgl(i)%rank==myrank) swild(ipgl(i)%id)=tmp
        enddo !i
        close(32)
        do i=1,nea
          hvis(:,i)=(swild(elnode(1,i))&
                    +swild(elnode(2,i))&
                    +swild(elnode(3,i)))/3
        enddo !i
      endif !ihorcon/=0
      
!...  Horizontal diffusivity option; only works for upwind/TVD
!     ihdif=0 means all hdif=0 and no hdif.gr3 is needed
      if(ihdif==0) then
        hdif=0
      else
        open(32,file='hdif.gr3',status='old')
        read(32,*)
        read(32,*) itmp1,itmp2
        if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check hdif.gr3')
        do i=1,np_global
          read(32,*)j,xtmp,ytmp,tmp 
          if(tmp<0) then
            write(errmsg,*)'hdif out of bound:',tmp,i
            call parallel_abort(errmsg)
          endif
          if(ipgl(i)%rank==myrank) hdif(:,ipgl(i)%id)=tmp
        enddo !i
        close(32)
      endif !ihdif/=0
      
!     Advection flags
      if(nadv==0) then
        open(10,file='adv.gr3',status='old')
        read(10,*)
        read(10,*) itmp1,itmp2
        if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check adv.gr3')
        do i=1,np_global
          read(10,*)j,xtmp,ytmp,tmp
          if(int(tmp)<0.or.int(tmp)>2) then
            write(errmsg,*)'Unknown iadv',i
            call parallel_abort(errmsg)
          endif
          if(ipgl(i)%rank==myrank) iadv(ipgl(i)%id)=int(tmp)
        enddo
        close(10)
      else !nadv/=0
        iadv=nadv
      endif

!...  Bottom friction
      if(nchi==-1) then !read in Manning's n for 2D model
        open(32,file='manning.gr3',status='old')
        read(32,*)
        read(32,*) itmp1,itmp2
        if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check manning.gr3')
        do i=1,np_global
          read(32,*)j,xtmp,ytmp,tmp
          if(tmp<0) then
            write(errmsg,*)'Negative Manning',tmp
            call parallel_abort(errmsg)
          endif
          if(ipgl(i)%rank==myrank) rmanning(ipgl(i)%id)=tmp
        enddo
        close(32)
      else if(nchi==0) then !read in drag coefficients for 3D model
        open(32,file='drag.gr3',status='old')
        read(32,*)
        read(32,*) itmp1,itmp2
        if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check drag.gr3')
        do i=1,np_global
          read(32,*)j,xtmp,ytmp,tmp
          if(tmp<0) then
            write(errmsg,*)'Negative bottom drag',tmp
            call parallel_abort(errmsg)
          endif
          if(ipgl(i)%rank==myrank) Cdp(ipgl(i)%id)=tmp
        enddo
        do i=1,nsa
          n1=isidenode(1,i)
          n2=isidenode(2,i)
          Cd(i)=(Cdp(n1)+Cdp(n2))/2
!         Debug
!          if(myrank==0) write(99,*)i,iplg(n1),iplg(n2),Cd(i)
        enddo
        close(32)
      else if(nchi==1) then !read in roughness in meters (3D)
        open(32,file='rough.gr3',status='old')
        read(32,*)
        read(32,*) itmp1,itmp2
        if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check rough.gr3')
        do i=1,np_global
          read(32,*)j,xtmp,ytmp,tmp
          if(ipgl(i)%rank==myrank) rough_p(ipgl(i)%id)=tmp
        enddo !i
        close(32)
      else
        write(errmsg,*)'Unknown bfric', nchi
        call parallel_abort(errmsg)
      endif !nchi

!     Coriolis options (must be 1 if ics=2)
      if(ncor==-1) then !latitude
        cori=coricoef
      else if(ncor==0) then
        cori=coricoef
      else !ncor=1
        if(ics==1.and.myrank==0) write(16,*)'Check slam0 and sfea0 as variable Coriolis is used'
!'
        open(32,file='hgrid.ll',status='old')
        read(32,*)
        read(32,*) !ne,np
        do i=1,np_global
          read(32,*)j,xtmp,ytmp
          if(ipgl(i)%rank==myrank) then
            xlon(ipgl(i)%id)=xtmp*pi/180
            ylat(ipgl(i)%id)=ytmp*pi/180
          endif
        enddo !i
        close(32)
        lreadll=.true.

        fc=2*omega_e*sin(sfea0)
        beta=2*omega_e*cos(sfea0)
        if(myrank==0) open(31,file='coriolis.out',status='replace')
        do i=1,nsa
          id1=isidenode(1,i)
          id2=isidenode(2,i)
          sphi=(ylat(id1)+ylat(id2))/2
          if(ics==1) then
            cori(i)=fc+beta*(sphi-sfea0)
          else !ics=2
            cori(i)=2*omega_e*sin(sphi)
          endif !ics
          if(myrank==0) write(31,*)i,xlon(id1)/pi*180,ylat(id1)/pi*180,cori(i)
        enddo !i=1,nsa
        if(myrank==0) close(31)
      endif !ncor

!     Wind 
      if(nws>=2.and.nws<=3) then !CORIE mode; read in hgrid.ll and open debug outputs
        open(32,file='hgrid.ll',status='old')
        read(32,*)
        read(32,*) !ne,np
        do i=1,np_global
          read(32,*)j,tmp1,tmp2
          if(ipgl(i)%rank==myrank) then
            xlon(ipgl(i)%id)=tmp1*pi/180
            ylat(ipgl(i)%id)=tmp2*pi/180
          endif
        enddo !i
        close(32)
        lreadll=.true.

#ifdef DEBUG
        fdb='sflux_0000'
        lfdb=len_trim(fdb)
        write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
        open(38,file='outputs/'//fdb,status='replace')
#endif
      endif

      windfactor=1 !intialize for default
      if(nws>0) then
        if(iwindoff/=0) then
          open(32,file='windfactor.gr3',status='old')
          read(32,*)
          read(32,*) itmp1,itmp2
          if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check windfactor.gr3')
          do i=1,np_global
            read(32,*)j,xtmp,ytmp,tmp
            if(tmp<0) then
              write(errmsg,*)'Wind scaling factor must be positive:',i,tmp
              call parallel_abort(errmsg)
            endif
            if(ipgl(i)%rank==myrank) windfactor(ipgl(i)%id)=tmp
          enddo !i
          close(32)
        endif
      endif !nws>0

!     Alloc. the large array for nws=4 option (may consider changing to unformatted binary read)
      if(nws==4) then
         allocate(rwild(np_global,3),stat=istat)
         if(istat/=0) call parallel_abort('MAIN: failed to alloc. (71)')
      endif !nws=4

!     Heat and salt conservation flags
      if(ihconsv/=0) then
        if(myrank==0) then
          write(16,*)'Warning: you have chosen a heat conservation model'
          write(16,*)'which assumes start time specified in sflux_inputs.txt!'
        endif

!       Read in albedo
        open(32,file='albedo.gr3',status='old')
        read(32,*)
        read(32,*) itmp1,itmp2
        if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check albedo.gr3')
        do i=1,np_global
          read(32,*)j,xtmp,ytmp,tmp
          if(tmp<0.or.tmp>1) then
            write(errmsg,*)'Albedo out of bound:',i,tmp
            call parallel_abort(errmsg)
          endif
          if(ipgl(i)%rank==myrank) albedo(ipgl(i)%id)=tmp
        enddo !i
        close(32)

!       Read in water type; the values for R, d_1, d_2 are given below 
!       solar flux= R*exp(z/d_1))+(1-R)*exp(z/d_2) (d_[1,2] are attentuation depths; smaller values for muddier water)
!       1: 0.58 0.35 23 (Jerlov type I)
!       2: 0.62 0.60 20 (Jerlov type IA)
!       3: 0.67 1.00 17 (Jerlov type IB)
!       4: 0.77 1.50 14 (Jerlov type II)
!       5: 0.78 1.40 7.9 (Jerlov type III)
!       6: 0.62 1.50 20 (Paulson and Simpson 1977; similar to type IA)
!       7: 0.80 0.90 2.1 (Mike Z.'s choice for estuary)
        open(32,file='watertype.gr3',status='old')
        read(32,*)
        read(32,*) itmp1,itmp2
        if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check watertype.gr3')
        do i=1,np_global
          read(32,*)j,xtmp,ytmp,tmp
          if(int(tmp)<1.or.int(tmp)>7) then
            write(errmsg,*)'Unknown water type:',i,tmp
            call parallel_abort(errmsg)
          endif
          if(ipgl(i)%rank==myrank) iwater_type(ipgl(i)%id)=tmp
        enddo !i
        close(32)
      endif
   
!...  Turbulence closure options
      if(itur==0) then
        dfv=dfv0; dfh=dfh0
      else if(itur==-1) then !VVD
        open(10,file='vvd.in',status='old')
        read(10,*) !nvrt
        do j=1,nvrt
          read(10,*)k,dfv0,dfh0
          dfv(j,:)=dfv0
          dfh(j,:)=dfh0
        enddo !j
        close(10)
      else if(itur==-2) then !HVD
        open(10,file='hvd.mom',status='old')
        open(32,file='hvd.tran',status='old')
        read(10,*)
        read(10,*) !np
        read(32,*)
        read(32,*) !np
        do i=1,np_global
          read(10,*)k,xtmp,ytmp,dfv0
          read(32,*)k,xtmp,ytmp,dfh0
          if(ipgl(i)%rank==myrank) then
            dfv(:,ipgl(i)%id)=dfv0
            dfh(:,ipgl(i)%id)=dfh0
          endif
        enddo !i
        close(10)
        close(32)
      else if(itur==2) then !read in P&P coefficients
      else if(itur==3.or.itur==4) then !read in const. (cf. Umlauf and Burchard 2003)
!       Common variables for both models
        cmiu0=sqrt(0.3d0)
!       read in mixing limits
        open(31,file='diffmin.gr3',status='old')
        open(32,file='diffmax.gr3',status='old')
        read(31,*)
        read(31,*) itmp1,itmp2
        if(itmp1/=ne_global.or.itmp2/=np_global) & 
     &call parallel_abort('Check diffmin.gr3')
        read(32,*)
        read(32,*) itmp1,itmp2
        if(itmp1/=ne_global.or.itmp2/=np_global) & 
     &call parallel_abort('Check diffmax.gr3')
        do i=1,np_global
          read(31,*)j,xtmp,ytmp,tmp1
          read(32,*)j,xtmp,ytmp,tmp2
          if(tmp2<tmp1) then
            write(errmsg,*)'diffmin > diffmax:',i
            call parallel_abort(errmsg)
          endif
          if(ipgl(i)%rank==myrank) then
            diffmin(ipgl(i)%id)=tmp1
            diffmax(ipgl(i)%id)=tmp2
          endif
        enddo !i
        close(31)
        close(32)

        if(itur==3) then
!	  Constants used in GLS; cpsi3 later
          a2_cm03=2/cmiu0**3
          eps_min=1.e-12

          select case(mid)
            case('MY') 
              rpub=0; rmub=1; rnub=1; cpsi1=0.9; cpsi2=0.5
              q2min=5.e-6; psimin=1.e-8
              if(stab.ne.'GA') then
                write(errmsg,*)'MY must use Galperins ASM:',stab
                call parallel_abort(errmsg)
              endif
            case('KL')
              rpub=0; rmub=1; rnub=1; schk=2.44; schpsi=2.44; cpsi1=0.9; cpsi2=0.5
              q2min=5.e-6; psimin=1.e-8
            case('KE')
              rpub=3; rmub=1.5; rnub=-1; schk=1; schpsi=1.3; cpsi1=1.44; cpsi2=1.92
              !q2min=1.0e-9; psimin=1.e-8
              q2min=7.6e-6; psimin=1.e-12 !Warner et al., OM, 2005, pp. 87
            case('KW')
              rpub=-1; rmub=0.5; rnub=-1; schk=2; schpsi=2; cpsi1=0.555; cpsi2=0.833
              !q2min=1.0e-9; psimin=1.e-8 
              q2min=7.6e-6; psimin=1.e-12 !Warner et al., OM, 2005, pp. 87
            case('UB')
              rpub=2; rmub=1; rnub=-0.67; schk=0.8; schpsi=1.07; cpsi1=1; cpsi2=1.22
              !q2min=1.0e-9; psimin=1.e-8 
              q2min=7.6e-6; psimin=1.e-12 !Warner et al., OM, 2005, pp. 87
            case default
              write(errmsg,*)'Unknown turb_met:',mid
              call parallel_abort(errmsg)
          end select
          if(rnub==0) then
            write(errmsg,*)'Wrong input for rnub:',rnub
            call parallel_abort(errmsg)
          endif

          if(stab.ne.'GA'.and.stab.ne.'KC') then
            write(errmsg,*)'Unknown turb_stab:',stab
            call parallel_abort(errmsg)
          endif

!	  Consts. used in Canuto's ASM (Model A)
          ubl1=0.1070
          ubl2=0.0032
          ubl3=0.0864
          ubl4=0.12
          ubl5=11.9
          ubl6=0.4
          ubl7=0
          ubl8=0.48
          ubs0=1.5*ubl1*ubl5**2
          ubs1=-ubl4*(ubl6+ubl7)+2*ubl4*ubl5*(ubl1-ubl2/3-ubl3)+1.5*ubl1*ubl5*ubl8
          ubs2=-0.375*ubl1*(ubl6**2-ubl7**2)
          ubs4=2*ubl5
          ubs5=2*ubl4
          ubs6=2*ubl5/3*(3*ubl3**2-ubl2**2)-0.5*ubl5*ubl1*(3*ubl3-ubl2)+0.75*ubl1*(ubl6-ubl7)
          ubd0=3*ubl5**2
          ubd1=ubl5*(7*ubl4+3*ubl8)
          ubd2=ubl5**2*(3*ubl3**2-ubl2**2)-0.75*(ubl6**2-ubl7**2)
          ubd3=ubl4*(4*ubl4+3*ubl8)
          ubd4=ubl4*(ubl2*ubl6-3*ubl3*ubl7-ubl5*(ubl2**2-ubl3**2))+ubl5*ubl8*(3*ubl3**2-ubl2**2)
          ubd5=0.25*(ubl2**2-3*ubl3**2)*(ubl6**2-ubl7**2)
!  	  print*, 'ubd2=',ubd2,',ubd4=',ubd4,',ubd2/ubd4=',ubd2/ubd4

!         Initialize k and l
          do i=1,npa
            xlmin2(i)=2*q2min*0.1*max(h0,dp(i)) !min. xl for non-surface layers
            q2(:,i)=q2min
            xl(:,i)=xlmin2(i)
          enddo !i
          dfv=0; dfh=0; dfq1=0; dfq2=0 !initialize for closure eqs.

!         Read in xlsc0
          open(32,file='xlsc.gr3',status='old')
          read(32,*)
          read(32,*) itmp1,itmp2
          if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check xlsc.gr3')
          do i=1,np_global
            read(32,*)j,xtmp,ytmp,tmp
            if(tmp<0.or.tmp>1) then
              write(errmsg,*)'Wrong xlsc0:',i,tmp
              call parallel_abort(errmsg)
            endif
            if(ipgl(i)%rank==myrank) xlsc0(ipgl(i)%id)=tmp
          enddo !i
          close(32)

        else !itur=4
#ifndef USE_GOTM
          write(errmsg,*)'Compile with GOTM:',itur
          call parallel_abort(errmsg)
#endif

        endif !itur=3 or 4
      endif !itur

!...  Interpolation flag for S,T and vel. in ELM
!     Kriging in vel: no bnd nodes/sides vel. use Kriging as the filter is not applied there
      if(lq/=0) then
        lqk=lq
      else !lq==0
        open(32,file='lqk.gr3',status='old')
        read(32,*)
        read(32,*) itmp1,itmp2
        if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check lqk.gr3')
        do i=1,np_global
          read(32,*)j,xtmp,ytmp,tmp
          if(tmp<1.or.tmp>2) then
            write(errmsg,*)'Unknown interpolation flag in lqk.gr3:',i
            call parallel_abort(errmsg)
          endif
          if(ipgl(i)%rank==myrank) swild(ipgl(i)%id)=tmp
        enddo !i
        close(32)
        do i=1,nea
          n1=elnode(1,i); n2=elnode(2,i); n3=elnode(3,i)
          lqk(i)=min(swild(n1),swild(n2),swild(n3))
        enddo !i
      endif
      
      if(inter_mom/=-1) then
        krvel=inter_mom
      else !-1
        open(32,file='krvel.gr3',status='old')
        read(32,*)
        read(32,*) itmp1,itmp2
        if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check krvel.gr3')
        do i=1,np_global
          read(32,*)j,xtmp,ytmp,tmp
          if(tmp<0.or.tmp>1) then
            write(errmsg,*)'Unknown interpolation flag in krvel.gr3:',i
            call parallel_abort(errmsg)
          endif
          if(ipgl(i)%rank==myrank) swild(ipgl(i)%id)=tmp
        enddo !i
        close(32)
        do i=1,nea
          n1=elnode(1,i); n2=elnode(2,i); n3=elnode(3,i)
          krvel(i)=min(swild(n1),swild(n2),swild(n3))
        enddo !i
      endif

!...  Land b.c. option (inactive)
!      read(15,*) !islip !0: free slip; otherwise no slip
      islip=0
!      if(islip/=0.and.islip/=1) then
!        write(errmsg,*)'Unknow islip:',islip
!        call parallel_abort(errmsg)
!      endif
!      if(islip==1) read(15,*) hdrag0

!...  Sponge layer for elev. & vel. (relax. factor applied to 0 elev. or uv -similar to T,S)
      if(inu_elev==1) then
        open(10,file='elev_nudge.gr3',status='old')
        read(10,*)
        read(10,*) itmp1,itmp2
        if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check elev_nudge.gr3')
        do i=1,np_global
          read(10,*)j,xtmp,ytmp,tmp1
          if(tmp1<0.or.tmp1*dt>1) then
            write(errmsg,*)'Wrong nudging factor at node (1):',i,tmp1
            call parallel_abort(errmsg)
          endif
          if(ipgl(i)%rank==myrank) elev_nudge(ipgl(i)%id)=tmp1 !Dimension: sec^-1
        enddo !i
        close(10)
      endif !inu_elev

      if(inu_uv==1) then
        open(10,file='uv_nudge.gr3',status='old')
        read(10,*)
        read(10,*) itmp1,itmp2
        if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check uv_nudge.gr3')
        do i=1,np_global
          read(10,*)j,xtmp,ytmp,tmp1
          if(tmp1<0.or.tmp1*dt>1) then
            write(errmsg,*)'Wrong nudging factor at node (2):',i,tmp1
            call parallel_abort(errmsg)
          endif
          if(ipgl(i)%rank==myrank) uv_nudge(ipgl(i)%id)=tmp1 !Dimension: sec^-1
        enddo !i
        close(10)
      endif !inu_uv

!...  Nudging options for T,S
      if(inu_st/=0) then
        if(vnh1>=vnh2.or.vnf1<0.or.vnf1>1.or.vnf2<0.or.vnf2>1) then
          write(errmsg,*)'Check vertical nudging limits:',vnh1,vnf1,vnh2,vnf2
          call parallel_abort(errmsg)
        endif

        open(10,file='t_nudge.gr3',status='old')
        open(32,file='s_nudge.gr3',status='old')
        read(10,*)
        read(10,*) itmp1,itmp2
        if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check t_nudge.gr3')
        read(32,*)
        read(32,*) itmp1,itmp2
        if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check s_nudge.gr3')
        do i=1,np_global
          read(10,*)j,xtmp,ytmp,tmp1
          read(32,*)j,xtmp,ytmp,tmp2
          if(tmp1<0.or.tmp1*dt>1.or.tmp2<0.or.tmp2*dt>1) then
            write(errmsg,*)'Wrong nudging factor at node:',i,tmp1,tmp2
            call parallel_abort(errmsg)
          endif
          if(ipgl(i)%rank==myrank) then
            !Dimension: sec^-1
            t_nudge(ipgl(i)%id)=tmp1
            s_nudge(ipgl(i)%id)=tmp2
          endif
        enddo !i
        close(10)
        close(32)

        if(inu_st==2) then
!          nrec_nu=nbyte*(1+np*nvrt) !single precision
          open(37,file='temp_nu.in',form='unformatted',status='old')
          open(35,file='salt_nu.in',form='unformatted',status='old')
        endif
      endif !inu_st

!...  Nudging for tracers
      if(inu_tr/=0) then
        open(10,file='tracer_nudge.gr3',status='old')
        read(10,*)
        read(10,*) itmp1,itmp2
        if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check tracer_nudge.gr3')
        do i=1,np_global
          read(10,*)j,xtmp,ytmp,tmp1
          if(tmp1<0.or.tmp1>1) then
            write(errmsg,*)'Wrong nudging factor at node (1):',i,tmp1
            call parallel_abort(errmsg)
          endif
          if(ipgl(i)%rank==myrank) tr_nudge(ipgl(i)%id)=tmp1
        enddo !i
        close(10)

        if(inu_tr==2) open(45,file='tr_nu.in',form='unformatted',status='old') !single precision
      endif !inu_tr/=0

!...  Surface min. mixing length for f.s. and max. for all; inactive 
!      read(15,*) !xlmax00

!     TVD scheme will be used if itvd_e=1 and min(total depth @ 3 nodes) >=h_tvd. itvd_e and h_tvd are shared 
!     between T,S and all tracers
      itvd_e=0 !init. for upwind
      if(iupwind_t>=2.or.itr_met>=2) then
        open(32,file='tvd.prop',status='old')
        do i=1,ne_global
          read(32,*)j,tmp
          itmp=nint(tmp)
          if(itmp/=0.and.itmp/=1) then
            write(errmsg,*)'Unknown TVD flag:',i,tmp
            call parallel_abort(errmsg)
          endif
          if(iegl(i)%rank==myrank) itvd_e(iegl(i)%id)=itmp
        enddo !i
        close(32)
      endif !iupwind_t

!     Station output option (/=0: need station.in)
!     If ics=2, the coord. in station.in must be in lat/lon (degrees)
      if(iout_sta/=0) then
        nvar_sta=9 !# of output variables
        allocate(iof_sta(nvar_sta),stat=istat)
        if(istat/=0) call parallel_abort('Sta. allocation failure (1)')
        open(32,file='station.in',status='old')
!       Output variables in order: elev, air pressure, windx, windy, T, S, u, v, w
        read(32,*)iof_sta(1:nvar_sta) !on-off flag for each variable
        read(32,*)nout_sta
!       Following is needed for dimension of nwild2
        if(nout_sta>ne_global) call parallel_abort('MAIN: too many stations')
!'

!       Allocate: zstal is vertical up; xsta, ysta, zsta are global coord. if ics=2
        allocate(xsta(nout_sta),ysta(nout_sta),zstal(nout_sta),zsta(nout_sta),iep_sta(nout_sta),iep_flag(nout_sta), &
     &arco_sta(nout_sta,3),sta_out(nout_sta,nvar_sta),sta_out_gb(nout_sta,nvar_sta), &
     &sta_out3d(nvrt,nout_sta,nvar_sta),sta_out3d_gb(nvrt,nout_sta,nvar_sta), &
     &zta_out3d(nvrt,nout_sta,nvar_sta),zta_out3d_gb(nvrt,nout_sta,nvar_sta),stat=istat)
        if(istat/=0) call parallel_abort('MAIN: sta. allocation failure')

        do i=1,nout_sta
          read(32,*)j,xsta(i),ysta(i),zstal(i) !z not used for 2D variables; xsta, ysta in lat/lon if ics=2
          if(ics==2) then
            xtmp=xsta(i)/180*pi
            ytmp=ysta(i)/180*pi
            xsta(i)=rearth*cos(ytmp)*cos(xtmp)
            ysta(i)=rearth*cos(ytmp)*sin(xtmp)
            zsta(i)=rearth*sin(ytmp)
          endif !ics
        enddo !i
        close(32)

!       Find parent elements and initialize outputs
        iep_sta=0 !flag for no-parent
        do i=1,ne
          do l=1,nout_sta
            if(iep_sta(l)/=0) cycle
            do j=1,3
              n1=elnode(nx(j,1),i)
              n2=elnode(nx(j,2),i)
              if(ics==1) then
                xn1=xnd(n1)
                yn1=ynd(n1)
                xn2=xnd(n2)
                yn2=ynd(n2)
                xstal=xsta(l)
                ystal=ysta(l)
              else !to eframe
                call project_pt('g2l',xnd(n1),ynd(n1),znd(n1), &
     &(/xctr(i),yctr(i),zctr(i)/),eframe(:,:,i),xn1,yn1,tmp)
                call project_pt('g2l',xnd(n2),ynd(n2),znd(n2), &
     &(/xctr(i),yctr(i),zctr(i)/),eframe(:,:,i),xn2,yn2,tmp)
                call project_pt('g2l',xsta(l),ysta(l),zsta(l), &
     &(/xctr(i),yctr(i),zctr(i)/),eframe(:,:,i),xstal,ystal,tmp)
              endif !ics
              swild(j)=signa(xn1,xn2,xstal,yn1,yn2,ystal) !temporary storage
            enddo !j
            ae=minval(swild(1:3)/area(i)) !abs(sum(abs(swild(1:3)))-area(i))/area(i)
            !if(ae<=small2) then
            if(ae>-small1) then
              iep_sta(l)=i
              arco_sta(l,1:3)=swild(1:3)/area(i)
              arco_sta(l,1)=max(0.d0,min(1.d0,arco_sta(l,1)))
              arco_sta(l,2)=max(0.d0,min(1.d0,arco_sta(l,2)))
              if(arco_sta(l,1)+arco_sta(l,2)>1) then 
                arco_sta(l,3)=0
                arco_sta(l,2)=1-arco_sta(l,1)
              else
                arco_sta(l,3)=1-arco_sta(l,1)-arco_sta(l,2)
              endif
            endif !ae
          enddo !l; build pts

          ifl=0 !flag
          do l=1,nout_sta
            if(iep_sta(l)==0) then
              ifl=1
              exit
            endif
          end do !l
          if(ifl==0) exit
        enddo !i=1,ne

!       See if any pt is outside
        call mpi_reduce(iep_sta,nwild2,nout_sta,itype,MPI_SUM,0,comm,ierr)
        if(myrank==0) then
          do i=1,nout_sta
            if(nwild2(i)==0) write(16,*)'Station pts outside domain:',i
          enddo !i
        endif

!       Open output file from rank 0
        if(myrank==0) then
          do i=1,nvar_sta
            write(ifile_char,'(i03)')i
            ifile_char=adjustl(ifile_char)  !place blanks at end
            ifile_len=len_trim(ifile_char)
            if(ihot==2) then
              open(250+i,file='outputs/staout_'//ifile_char(1:ifile_len),status='old')
            else
              open(250+i,file='outputs/staout_'//ifile_char(1:ifile_len),status='replace')
            endif
          enddo !i
          write(16,*)'done preparing station outputs'
          call flush(16) ! flush "mirror.out"
        endif
      endif !iout_sta

#ifdef USE_HA
!...  Read harmonic analysis information (Adapted from ADCIRC)
      if(iharind/=0) then
        open(31,file='harm.in',status='old')
!...
!...  READ AND CHECK INFORMATION ABOUT HARMONIC ANALYSIS OF MODEL RESULTS
!...  
        READ(31,*) NFREQ 
        if (myrank==0) then
          WRITE(16,99392) NFREQ  
99392     FORMAT(////,1X,'HARMONIC ANALYSIS INFORMATION OUTPUT : ',//,5X,'HARMONIC ANALYSIS PERFORMED FOR ',I4,' CONSTITUENTS',/)
        endif
        MNHARF = NFREQ
        IF (NFREQ.EQ.0) MNHARF = 1

!...  allocate harmonic analysis arrays

        IF (NFREQ.GT.0) THEN
          CALL ALLOC_HA(np)
          ALLOCATE (XVELAV(np),YVELAV(np),XVELVA(np),YVELVA(np))
          ALLOCATE (ELAV(np),ELVA(np))
        ENDIF
      
        IF(NFREQ.LT.0) THEN
          WRITE(errmsg,99391)
99391     FORMAT(////,1X,'!!!!!!!!!!  WARNING - FATAL ERROR !!!!!!!!!',//,1X,'YOUR SELECTION OF NHARFR (A harm.in INPUT PARAMETER) IS NOT AN ALLOWABLE VALUE',/,1X,'PLEASE CHECK YOUR INPUT',//,1X,'!!!!!! EXECUTION WILL NOW BE TERMINATED !!!!!!',//)
          call parallel_abort(errmsg)
        ENDIF
        IF(NFREQ.GT.0 .AND. myrank.EQ. 0) WRITE(16,2330)
 2330   FORMAT(/,7X,'FREQUENCY',4X,'NODAL FACTOR',6X,'EQU.ARG(DEG)',1X,'CONSTITUENT',/)
!'
        DO I=1,NFREQ  
           READ(31,'(A10)') NAMEFR(I)
           READ(31,*) HAFREQ(I),HAFF(I),HAFACE(I)
           if (myrank==0) WRITE(16,2331) HAFREQ(I),HAFF(I),HAFACE(I),NAMEFR(I)
 2331      FORMAT(4X,F15.12,2X,F10.7,5X,F10.3,7X,A10)
        enddo

!...  read in interval information for harmonic analysis
!...  compute thas and thaf in terms of the number of time steps
        READ(31,*) THAS,THAF,NHAINC,FMV
        ITHAS=INT(THAS*(86400.D0/dt) + 0.5d0)
        THAS=ITHAS*dt/86400.D0
        ITHAF=INT(THAF*(86400.D0/dt) + 0.5d0)
        THAF=ITHAF*dt/86400.D0
        ITMV = ITHAF - (ITHAF-ITHAS)*FMV
        if (myrank==0) then
          IF(NFREQ.GT.0) THEN
            WRITE(16,34634) THAS,ITHAS,THAF,ITHAF,NHAINC
34634       FORMAT(/,5X,'HARMONIC ANALYSIS WILL START AFTER THAS =',F8.3,' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',I9,' TIME STEPS INTO THE SIMULATION',//,5X,'HARMONIC ANALYSIS WILL STOP AFTER THAF =',F8.3,' DAY(S) RELATIVE',/,9X,'TO THE STARTING TIME OR',I9,' TIME STEPS INTO THE SIMULATION',//,5X,'INFORMATION WILL BE ANALYZED EVERY ','NHAINC =',I8,' TIME STEPS.')
            WRITE(16,34639) FMV*100.,ITMV
34639       FORMAT(/,5X,'MEANS AND VARIANCES WILL BE COMPUTED FOR THE FINAL ',F10.5,' %',/9X,'OF THE HARMONIC ANALYSIS PERIOD OR AFTER ',I9,' TIME STEPS INTO THE SIMULATION.',/9X,' RESULTS ARE WRITTEN TO UNIT 55.')
!'
         
          ELSE
            WRITE(16,34645)
34645       FORMAT(///,5X,'NO HARMONIC ANALYSIS WILL BE DONE')
          ENDIF
        endif
      
        IF ((FMV.GT.0.).AND.(NFREQ.GT.0)) CHARMV = .TRUE.
      
!...  read in and write out information on where harmonic analysis will
!...  be done

        READ(31,*) NHAGE,NHAGV
        IF((NHAGE.LT.0).OR.(NHAGE.GT.1)) THEN
           WRITE(errmsg,99663)
99663      FORMAT(////,1X,'!!!!!!!!!!  WARNING - INPUT ERROR  !!!!!!!!!',//,1X,'YOUR SELECTION OF NHAGE (A harm.in INPUT PARAMETER) IS NOT AN ALLOWABLE VALUE',/,1X,'PLEASE CHECK YOUR INPUT')
           call parallel_abort(errmsg)
        ENDIF
        IF(NHAGE.EQ.1) THEN
           if (myrank==0) WRITE(16,34643)
34643      FORMAT(///,5X,'GLOBAL ELEVATION HARMONIC ANAL WILL BE WRITTEN TO FILE harme.53')
!'
        ENDIF
        IF((NHAGV.LT.0).OR.(NHAGV.GT.1)) THEN
          WRITE(errmsg,99664)
99664     FORMAT(////,1X,'!!!!!!!!!!  WARNING - INPUT ERROR  !!!!!!!!!',//,1X,'YOUR SELECTION OF NHAGV (A harm.ina INPUT PARAMETER) IS NOT AN ALLOWABLE VALUE',/,1X,'PLEASE CHECK YOUR INPUT')
          call parallel_abort(errmsg)
        ENDIF
        IF(NHAGV.EQ.1) THEN
          if (myrank==0) WRITE(16,34644)
34644     FORMAT(///,5X,'GLOBAL VELOCITY HARMONIC ANAL WILL BE WRITTEN TO FILE harmv.54')
!'
        ENDIF
        
!...  compute flag telling whether any harmonic analysis will be done
        iharind=NFREQ*(NHAGE+NHAGV)
        if (iharind.GT.0) iharind=1
        close(31)
      endif !iharind/=0
#endif /*USE_HA*/

!     Inter-subdomain backtracking
      call init_inter_btrack !setup datatypes and mxnbt and _r
!      allocate(btlist(mxnbt),stat=istat)
!      if(istat/=0) call parallel_abort('MAIN: btlist allocation failure')

!...  Compute neighborhood for 2-tier Kriging and invert matrix for resident elements only
!     Compute ne_kr for dimensioning
!      if(ics==2.and.inter_mom/=0) call parallel_abort('ics=2 and inter_mom/=0')
!'
      ie_kr=0 !non-zero value points to local nth Kriging elements
      ne_kr=0 !total # of elements in Kriging zone
      do i=1,ne !resident
        if(krvel(i)==1) then
          ne_kr=ne_kr+1
          ie_kr(i)=ne_kr
        endif
      enddo !i

!     Compute mnei_kr (max. # of Kriging pts) for dimensioning
      mnei_kr=3
      do i=1,ne !resident
        if(krvel(i)/=1) cycle
  
        nei_kr=3 !# of Kriging nodes for i
        nwild(1:3)=elnode(1:3,i) !temporarily save Kriging nodes
        do j=1,3 !resident nodes
          nd=elnode(j,i)
          loop14: do m=1,nnp(nd)
            nd2=indnd(m,nd)
            if(nd2<=0) then
              write(errmsg,*)'MAIN: node outside:',ielg(i),iplg(nd)
              call parallel_abort(errmsg)
            endif
            ! Check if present
            do l=1,nei_kr
              if(nwild(l)==nd2) cycle loop14
            enddo !l
            ! New node
            nei_kr=nei_kr+1
            nwild(nei_kr)=nd2
          enddo loop14 !m=1,nnp(nd)
        enddo !j=1,3 nodes
        if(nei_kr>mnei_kr) mnei_kr=nei_kr
      enddo !i=1,ne
      write(12,*)'Max. # of Kriging points = ',mnei_kr

!     Allocate arrays
      allocate(itier_nd(0:mnei_kr,ne_kr),akrmat_nd(mnei_kr+3,mnei_kr+3,ne_kr), &
              &akr(mnei_kr+3,mnei_kr+3),akrp((mnei_kr+3)*(mnei_kr+4)/2), &
              &xy_kr(2,mnei_kr),ipiv(mnei_kr+3),work4(mnei_kr+3))

!     Compute Kriging neighborhood
      do i=1,ne !resident
        if(krvel(i)/=1) cycle
 
        ie=ie_kr(i)
        nei_kr=3 !# of Kriging nodes for i
        itier_nd(1:3,ie)=elnode(1:3,i) !temporarily save Kriging nodes
        do j=1,3 !resident nodes
          nd=elnode(j,i)
          loop15: do m=1,nnp(nd)
            nd2=indnd(m,nd)
            ! Check if present
            do l=1,nei_kr
              if(itier_nd(l,ie)==nd2) cycle loop15
            enddo !l
            ! New node
            nei_kr=nei_kr+1
            itier_nd(nei_kr,ie)=nd2
          enddo loop15 !m=1,nnp(nd)
        enddo !j=1,3 nodes
        itier_nd(0,ie)=nei_kr

!       Debug
!        if(myrank==3) write(99,*)ielg(i),itier_nd(0,ie),iplg(itier_nd(1:nei_kr,ie))
      enddo !i=1,ne
      
!...  Invert Kriging matrices
      akrmat_nd=-1.e34 !initialization for debugging
      err_max=0 !max. error in computing the inverse matices
      do k=1,ne !resident
        if(ie_kr(k)==0) cycle

        ie=ie_kr(k) !local index
        npp=itier_nd(0,ie)
        do i=1,npp
          n1=itier_nd(i,ie)
          if(ics==2) then
            call project_pt('g2l',xnd(n1),ynd(n1),znd(n1), &
     &(/xctr(k),yctr(k),zctr(k)/),eframe(:,:,k),xy_kr(1,i),xy_kr(2,i),tmp)
          endif !ics

          do j=1,npp
            n2=itier_nd(j,ie)
            if(ics==1) then
              rr=sqrt((xnd(n1)-xnd(n2))**2+(ynd(n1)-ynd(n2))**2+(znd(n1)-znd(n2))**2)
            else
              call project_pt('g2l',xnd(n2),ynd(n2),znd(n2), &
     &(/xctr(k),yctr(k),zctr(k)/),eframe(:,:,k),xn2,yn2,tmp)
              rr=sqrt((xy_kr(1,i)-xn2)**2+(xy_kr(2,i)-yn2)**2)
            endif !ics
            akr(i,j)=covar(kr_co,rr)
          enddo !j

          akr(i,npp+1)=1
          if(ics==1) then
            akr(i,npp+2)=xnd(n1)
            akr(i,npp+3)=ynd(n1)
          else
            akr(i,npp+2)=xy_kr(1,i)
            akr(i,npp+3)=xy_kr(2,i)
          endif !ics
        enddo !i=1,npp

        akr(npp+1,1:npp)=1
        if(ics==1) then
          akr(npp+2,1:npp)=xnd(itier_nd(1:npp,ie))
          akr(npp+3,1:npp)=ynd(itier_nd(1:npp,ie))
        else
          akr(npp+2,1:npp)=xy_kr(1,1:npp)
          akr(npp+3,1:npp)=xy_kr(2,1:npp)
        endif !ics
        akr((npp+1):(npp+3),(npp+1):(npp+3))=0
!        bkr(1:(npp+3),1)=0 !does not matter

!       Debug
        akrmat_nd(1:(npp+3),1:(npp+3),ie)=akr(1:(npp+3),1:(npp+3))

!        call gaussj(akr,npp+3,mnei_kr+3,bkr,1,1)

!       LAPACK routines for positive definite symmetric matrix below did not work
!       Note: in LAPACK, the matrix dimension is (LDA,*) so the dimensions will match
!        call dpotrf('U',npp+3,akr,mnei_kr+3,info)
!        if(info/=0) then
!          write(11,*)'Failed dpotrf:',info
!          stop
!        endif
!        call dpotri('U',npp+3,akr,mnei_kr+3,info)
!        if(info/=0) then
!          write(11,*)'Failed dpotri:',info
!          stop
!        endif
!        do i=1,npp+3
!          do j=i+1,npp+3
!            akr(j,i)=akr(i,j)
!          enddo !j
!        enddo !i

!       Pack symmetric matrix
        do j=1,npp+3
          do i=1,j
            akrp(i+j*(j-1)/2)=akr(i,j)
          enddo !i
        enddo !j
        call dsptrf('U',npp+3,akrp,ipiv,info)
        if(info/=0) then
          write(errmsg,*)'MAIN: Failed dsptrf:',info,ielg(k),(i,(j,akr(i,j),j=1,npp+3),i=1,npp+3)
          call parallel_abort(errmsg) 
        endif
        call dsptri('U',npp+3,akrp,ipiv,work4,info)
        if(info/=0) then
          write(errmsg,*)'Failed dsptri:',info,ielg(k)
          call parallel_abort(errmsg)
        endif
!       Unpack
        do j=1,npp+3
          do i=1,j
            akr(i,j)=akrp(i+j*(j-1)/2)
          enddo !i
        enddo !j

        do i=1,npp+3
          do j=i+1,npp+3
            akr(j,i)=akr(i,j)
          enddo !j
        enddo !i

!       Check
        do i=1,npp+3
          do j=1,npp+3
            suma=0
            do l=1,npp+3
              suma=suma+akrmat_nd(i,l,ie)*akr(l,j)
            enddo !l
            if(i==j) suma=suma-1

!            if(k==22910) then
!              write(96,*)i,j,akrmat_nd(i,j,ie),akr(i,j),suma
!            endif

!            if(abs(suma)>1.e-8) write(98,*)k,i,j,suma
            if(abs(suma)>err_max) err_max=abs(suma)
          enddo !j
        enddo !i

        akrmat_nd(1:(npp+3),1:(npp+3),ie)=akr(1:(npp+3),1:(npp+3))
      enddo !k=1,ne

      write(12,*)'Max. error in inverting Kriging maxtrice= ',err_max

      deallocate(ipiv,akr,akrp,work4,xy_kr)
!...  End Kriging preparation

      if(myrank==0) then
        write(16,*)'done init (1)...'
        call flush(16)  ! flush "mirror.out"
      endif 

!	
!*******************************************************************
!
!	Initialization flow for both cold and hot start 
!
!*******************************************************************
!
!...  initialize elevations and vel.; may be overwritten by hotstart later 
!...  Outside ihot==0 loop for initializing levels0()
!     Read in elev.ic
      if(ic_elev==0) then
        eta2=0
      else    
        open(32,file='elev.ic',status='old')
        read(32,*)
        read(32,*)
        do i=1,np_global
          read(32,*)j,xtmp,ytmp,tmp
          if(ipgl(i)%rank==myrank) eta2(ipgl(i)%id)=tmp
        enddo !i
        close(32)
      endif

!     For ics=1, (su2,sv2) are defined in the _global_ Cartesian frame
!     For ics=2, they are defined in the _side_ frame
      su2=0; sv2=0
      we=0 !in element frame

!     Williamson test #5 - zonal flow over an isolated mount
      if(izonal5/=0) call zonal_flow

!...  Initialize levels
      if(inunfl==0) then
        call levels0(0,0)
      else
        call levels1(0,0)
      endif

!...  Bottom roughness length (m) for nchi==0
      if(nchi==0) then
        do j=1,npa
          if(idry(j)==1.or.Cdp(j)==0) then
            rough_p(j)=0
          else
            rough_p(j)=(znl(kbp(j)+1,j)-znl(kbp(j),j))*exp(-0.4/sqrt(Cdp(j)))
          endif
        enddo !j
      endif !nchi

!...  kbp00 will be used for output only
      kbp00=kbp

!...  Calculate mean density profile, and for ihot==0 & icst=2, initialize T,S 
!...  at nodes, sides and elements as well (which will be over-written for other cases)
      if(ibcc_mean==1.or.ihot==0.and.icst==2) then
!       Read in intial mean S,T
        open(32,file='ts.ic',status='old')
        read(32,*)nz_r
        if(nz_r<2) then
          write(errmsg,*)'Change nz_r:',nz_r
          call parallel_abort(errmsg)
        endif
        allocate(z_r(nz_r),tem1(nz_r),sal1(nz_r),cspline_ypp(nz_r,2),stat=istat)
        deallocate(swild,stat=istat)
        allocate(swild(max(nsa+nvrt+12+ntracers,nz_r)),stat=istat)
        do k=1,nz_r
          !z_r in local frame if ics=2
          read(32,*)j,z_r(k),tem1(k),sal1(k)
          if(tem1(k)<tempmin.or.tem1(k)>tempmax.or.sal1(k)<saltmin.or.sal1(k)>saltmax) then
            write(errmsg,*)'Initial invalid S,T at',k,tem1(k),sal1(k)
            call parallel_abort(errmsg)
          endif
          if(k>=2) then; if(z_r(k)<=z_r(k-1)) then
            write(errmsg,*)'Inverted z-level (0):',k
            call parallel_abort(errmsg)
          endif; endif
        enddo !k
        close(32)

!       Cubic spline coefficients (save for interpolation later)
        call cubic_spline(nz_r,z_r,tem1,0._rkind,0._rkind,swild)
        cspline_ypp(1:nz_r,1)=swild(1:nz_r)
        call cubic_spline(nz_r,z_r,sal1,0._rkind,0._rkind,swild)
        cspline_ypp(1:nz_r,2)=swild(1:nz_r)

!       T,S @ nodes
        do i=1,npa
          if(idry(i)==1) then
            tem0(:,i)=tem1(nz_r)
            sal0(:,i)=sal1(nz_r)
            cycle
          endif

!         Wet nodes
          if(znl(kbp(i),i)<z_r(1)) then !.or.znl(nvrt,i)>z_r(nz_r)) then
            call parallel_abort('MAIN: node depth too big for ts.ic')
          endif 
          call eval_cubic_spline(nz_r,z_r,tem1,cspline_ypp(1:nz_r,1),nvrt-kbp(i)+1,znl(kbp(i):nvrt,i), &
     &0,z_r(1),z_r(nz_r),tem0(kbp(i):nvrt,i))
          call eval_cubic_spline(nz_r,z_r,sal1,cspline_ypp(1:nz_r,2),nvrt-kbp(i)+1,znl(kbp(i):nvrt,i), &
     &0,z_r(1),z_r(nz_r),sal0(kbp(i):nvrt,i))

!         Impose no slip b.c. to be consistent with ELM transport
          if(Cdp(i)/=0.or.rough_p(i)/=0) then
            tem0(kbp(i),i)=tem0(kbp(i)+1,i)
            sal0(kbp(i),i)=sal0(kbp(i)+1,i)
          endif

!         Extend
          do k=1,kbp(i)-1
            tem0(k,i)=tem0(kbp(i),i)
            sal0(k,i)=sal0(kbp(i),i)
          enddo !k
        enddo !i=1,npa

!       T,S @ sides 
        do i=1,nsa
          if(idry_s(i)==1) then
            tsd(:,i)=tem1(nz_r)
            ssd(:,i)=sal1(nz_r)
            cycle
          else !wet side
            if(zs(kbs(i),i)<z_r(1)) then !.or.zs(nvrt,i)>z_r(nz_r)) then
              call parallel_abort('MAIN: side depth too big for ts.ic')
            endif 
            call eval_cubic_spline(nz_r,z_r,tem1,cspline_ypp(1:nz_r,1),nvrt-kbs(i)+1,zs(kbs(i):nvrt,i), &
     &0,z_r(1),z_r(nz_r),tsd(kbs(i):nvrt,i))
            call eval_cubic_spline(nz_r,z_r,sal1,cspline_ypp(1:nz_r,2),nvrt-kbs(i)+1,zs(kbs(i):nvrt,i), &
     &0,z_r(1),z_r(nz_r),ssd(kbs(i):nvrt,i))
          endif !wet or dry side

!         Extend
          do k=1,kbs(i)-1
            tsd(k,i)=tsd(kbs(i),i)
            ssd(k,i)=ssd(kbs(i),i)
          enddo !k
        enddo !i=1,nsa

!       T,S @ elements
        do i=1,nea
          if(idry_e(i)==1) then
            tsel(1,:,i)=tem1(nz_r)
            tsel(2,:,i)=sal1(nz_r)
            cycle
          else !wet element
            if(ze(kbe(i),i)<z_r(1)) then !.or.ze(nvrt,i)>z_r(nz_r)) then
              call parallel_abort('MAIN: ele. depth too big for ts.ic')
            endif 

            do k=kbe(i)+1,nvrt
              swild(k)=(ze(k,i)+ze(k-1,i))/2
            enddo !k
            call eval_cubic_spline(nz_r,z_r,tem1,cspline_ypp(1:nz_r,1),nvrt-kbe(i),swild(kbe(i)+1:nvrt), &
     &0,z_r(1),z_r(nz_r),tsel(1,kbe(i)+1:nvrt,i))
            call eval_cubic_spline(nz_r,z_r,sal1,cspline_ypp(1:nz_r,2),nvrt-kbe(i),swild(kbe(i)+1:nvrt), &
     &0,z_r(1),z_r(nz_r),tsel(2,kbe(i)+1:nvrt,i))
          endif !wet or dry elements

!         Extend
          do k=1,kbe(i)
            tsel(1:2,k,i)=tsel(1:2,kbe(i)+1,i)
          enddo !k
        enddo !i=1,nea
      endif !ibcc_mean==1.or.ihot==0.and.icst==2

!								   
!*******************************************************************
!								   
!	Initialization for cold start alone
!								   
!*******************************************************************
!								   
      if(ihot==0) then
!------------------------------------------------------------------
!...  read the initial salinity and temperature values from 
!...  salt.ic and temp.ic files. Initial S,T fields may vary
!...  either horizontally (and vertically homogeneous) or vertically 
!...  (horizontally homogeneous). For more general 3D case, use hot start.
!...
      if(ibc==1.and.ibtp==0) then
!	Reset icst
        icst=1
        tem0=10; sal0=0; tsd=10; ssd=0; tsel(1,:,:)=10; tsel(2,:,:)=0
      else !read in S,T
        if(icst==1) then
          open(31,file='temp.ic',status='old')
          open(32,file='salt.ic',status='old')
          read(31,*) 
          read(31,*) !np
          do i=1,np_global
            read(31,*) itmp,xtmp,ytmp,te
            if(te<tempmin.or.te>tempmax) then
              write(errmsg,*)'Initial invalid T at',i,te
              call parallel_abort(errmsg)
            endif
            if(ipgl(i)%rank==myrank) tem0(:,ipgl(i)%id)=te
          enddo !i

          read(32,*) 
          read(32,*) !np
          do i=1,np_global
            read(32,*) itmp,xtmp,ytmp,sa
            if(sa<saltmin.or.sa>saltmax) then
              write(errmsg,*)'Initial invalid S at',i,sa
              call parallel_abort(errmsg)
            endif
            if(ipgl(i)%rank==myrank) sal0(:,ipgl(i)%id)=sa
          enddo
          close(31)
          close(32)

!         T,S @ sides and elements
          do i=1,nsa
            n1=isidenode(1,i)
            n2=isidenode(2,i)
            do k=1,nvrt
              tsd(k,i)=(tem0(k,n1)+tem0(k,n2))/2
              ssd(k,i)=(sal0(k,n1)+sal0(k,n2))/2
            enddo !k
          enddo !i

          do i=1,nea
            n1=elnode(1,i)
            n2=elnode(2,i)
            n3=elnode(3,i)
            do k=2,nvrt
              tsel(1,k,i)=(tem0(k,n1)+tem0(k,n2)+tem0(k,n3)+tem0(k-1,n1)+tem0(k-1,n2)+tem0(k-1,n3))/6
              tsel(2,k,i)=(sal0(k,n1)+sal0(k,n2)+sal0(k,n3)+sal0(k-1,n1)+sal0(k-1,n2)+sal0(k-1,n3))/6
            enddo !k
            tsel(1,1,i)=tsel(1,2,i) !mainly for hotstart format
            tsel(2,1,i)=tsel(2,2,i)
          enddo !i

        else !icst=2 
!         Already initialized

        endif !icst
      endif !ibc.eq.1.and.ibtp.eq.0

!...  initialize S,T @ nodes
      tnd=tem0; snd=sal0

!     Debug
!      if(myrank==0) then
!        itmp=7898
!        write(98,'(3(1x,f10.3))')(znl(k,itmp),tnd(k,itmp),snd(k,itmp),k=kbp(itmp),nvrt)
!        write(98,*)
!        write(98,'(3(1x,f10.3))')(zs(k,itmp),tsd(k,itmp),ssd(k,itmp),k=kbs(itmp),nvrt)
!        write(98,*)
!        write(98,'(3(1x,f10.3))')((ze(k,itmp)+ze(k-1,itmp))/2,tsel(1:2,k,itmp),k=kbe(itmp)+1,nvrt)
!      endif
!      call parallel_finalize
!      stop

!...  initialize wind for nws=1,2 (first two lines)
!...  Wind vector always in lat/lon frame and so will have problem at poles
      if(nws==1) then
        open(22,file='wind.th',status='old')
        read(22,*)tmp1,wx1,wy1
        read(22,*)tmp2,wx2,wy2
        if(abs(tmp1)>1.e-4.or.abs(tmp2-wtiminc)>1.e-4) &
     &call parallel_abort('check time stamp in wind.th')
        do i=1,npa
          windx1(i)=wx1
          windy1(i)=wy1
          windx2(i)=wx2
          windy2(i)=wy2
        enddo
        wtime1=0
        wtime2=wtiminc 
      endif

      if(nws==4) then
        open(22,file='wind.th',status='old')
        read(22,*)tmp1,rwild(:,:)
        do i=1,np_global
          if(ipgl(i)%rank==myrank) then
            nd=ipgl(i)%id
            windx1(nd)=rwild(i,1)
            windy1(nd)=rwild(i,2)
            pr1(nd)=rwild(i,3)
          endif
        enddo !i

        read(22,*)tmp2,rwild(:,:)
        do i=1,np_global
          if(ipgl(i)%rank==myrank) then
            nd=ipgl(i)%id
            windx2(nd)=rwild(i,1)
            windy2(nd)=rwild(i,2)
            pr2(nd)=rwild(i,3)
          endif
        enddo !i
        if(abs(tmp1)>1.e-4.or.abs(tmp2-wtiminc)>1.e-4) &
     &call parallel_abort('check time stamp in wind.th (4)')

        wtime1=0
        wtime2=wtiminc
      endif !nws=4

!	CORIE mode
      if(nws>=2.and.nws<=3) then
        wtime1=0
        wtime2=wtiminc 
!       wind speed upon output is rotated to the map projection
!       For ics=2, make sure windrot* =0 (i.e. true east/north direction)
        call get_wind(wtime1,windx1,windy1,pr1,airt1,shum1)
        call get_wind(wtime2,windx2,windy2,pr2,airt2,shum2)

!       Uncomment the following to overwrite wind (similary airt etc)
!       Search for "Overwrite wind with wind.th"; 
!       WARNING: wind.th has time step of dt not wtiminc, and 
!       starts from t=0,dt,2*dt ... (format: windu,windv; similar to nws=1).
!       Overwrite wind with wind.th

!        open(22,file='wind.th',status='old')
!        read(22,*)tmp,wx1,wy1
!        windx1=wx1; windx2=wx1
!        windy1=wy1; windy2=wy1
!       End

        if(nws==3) then
!         Open sflux.th (time interval is dt, not wtiminc!)
          open(23,file='sflux.th',status='old')
          read(23,*) !t=0
!         To heat up water, fluxsu00<0, srad00>0
          read(23,*) tmp,fluxsu00,srad00 !time, total surface flux, solar radiation
        endif !nws==3
      endif !nws>=2

!     VIMS Point source loading added by YC
#ifdef USE_ICM 
      if(myrank==0) write(16,*)'start reading ICM point source...'
      call WQCO1(dt,rnday,NDTWQ)
      if(iWQPS==2) then
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
!
        x1 = 1.0E3 !conversion from kg to g
        open(61,file='ps.in',status='old')
        read(61,*)
        read(61,*)
        do i=1,nps
          read(61,*) iegb,xPSK,xPSQ,PRPOC,PLPOC,PDOCA,PRPON,PLPON,PDON, &
          &             PNH4,PNO3,PRPOP,PLPOP,PDOP,PPO4t,PSU,PSAt,PCOD,PDO
          if(iegl(iegb)%rank==myrank) then
            PSQ(iegl(iegb)%id)     = xPSQ !m^3/s - usually negative
            PSK(iegl(iegb)%id)     = xPSK !vertical layer where the source is applied
            WWPRPOC(iegl(iegb)%id) = PRPOC * x1  ! kg/d * 10^3 = g per day
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
        npstime=0
        npstiminc=86400.  !added by YC, users can change by themselves
!org yc        npstime1=0
!org yc        npstime2=npstiminc 
        if(myrank==0) write(16,*)'end reading ICM point source...'
      endif ! iWQPS=2

!     VIMS surface temperature mode added by YC
!     May want to use heat exchange instead
      if(iSun==2) then
        if(myrank==0) write(16,*)'start reading ICM surface T...'
        open(62,file='surface_t.th',status='old')
        read(62,*)
        read(62,*)
        read(62,*)
        do i=1,np_global
          read(62,*) ipgb,xtmp
          if(ipgl(ipgb)%rank==myrank) then
           surf_t1(ipgl(ipgb)%id)=xtmp
          endif
        enddo
        read(62,*)
        do i=1,np_global
          read(62,*) ipgb,xtmp
          if(ipgl(ipgb)%rank==myrank) then
           surf_t2(ipgl(ipgb)%id)=xtmp
          endif
        enddo
        surf_time1=0.
        surf_time2=86400. 

        if(myrank==0) write(16,*)'end reading ICM surface T...'
      endif ! iSUN=2   !added by YC
#endif /*USE_ICM*/

!...  Read initial nudging S,T
      if(inu_st==2) then
        read(37)floatout
        read(35)floatout
        do i=1,np_global
          read(37)(swild8(j,1),j=1,nvrt)
          read(35)(swild8(j,2),j=1,nvrt)
          if(ipgl(i)%rank==myrank) then
            tnd_nu1(:,ipgl(i)%id)=swild8(1:nvrt,1)
            snd_nu1(:,ipgl(i)%id)=swild8(1:nvrt,2)
          endif
        enddo !i
        read(37)floatout
        read(35)floatout
        do i=1,np_global
          read(37)(swild8(j,1),j=1,nvrt)
          read(35)(swild8(j,2),j=1,nvrt)
          if(ipgl(i)%rank==myrank) then
            tnd_nu2(:,ipgl(i)%id)=swild8(1:nvrt,1)
            snd_nu2(:,ipgl(i)%id)=swild8(1:nvrt,2)
          endif
        enddo !i
        irec_nu=2
        time_nu=step_nu
      endif !inu_st

!...  Read initial nudging for tracers
      if(ntracers>0.and.inu_tr==2) then
        read(45)floatout
        do i=1,np_global
          read(45)swild9
          if(ipgl(i)%rank==myrank) then
            trnd_nu1(:,:,ipgl(i)%id)=swild9
          endif
        enddo !i
        read(45)floatout
        do i=1,np_global
          read(45)swild9
          if(ipgl(i)%rank==myrank) then
            trnd_nu2(:,:,ipgl(i)%id)=swild9
          endif
        enddo !i
        irec_nu_tr=2
        time_nu_tr=step_nu_tr
      endif !inu_tr

#ifdef USE_HA
!...
!....INITIALIZE HARMONIC ANALYSIS MATRICES, MEAN AND SQUARE VECTORS
!... Adapted from ADCIRC
      IF (iharind.EQ.1) THEN
        ICHA=0
        CALL HACOLDS(HAFREQ)
        IF(NHAGE.EQ.1) CALL HACOLDSEG(np)
        IF(NHAGV.EQ.1) CALL HACOLDSVG(np)
        IF (CHARMV) THEN
          ELAV  =0.D0
          XVELAV=0.D0
          YVELAV=0.D0
          ELVA  =0.D0
          XVELVA=0.D0
          YVELVA=0.D0
!          DO I=1,np
!            ELAV(I)=0.D0
!            XVELAV(I)=0.D0
!            YVELAV(I)=0.D0
!            ELVA(I)=0.D0
!            XVELVA(I)=0.D0
!            YVELVA(I)=0.D0
!          ENDDO
        ENDIF
      ENDIF
#endif /*USE_HA*/

!------------------------------------------------------------------
      endif !ihot=0

!...  Finish off init. for both cold and hotstart
!...  Tracers; user-defined tracer part
!     This part needs T,S i.c. (tsel)
      if(ntracers>0) then
        tr_nd=0
!        trel0(1,:,:)=1 
!        trel0(2,:,:)=0 
!        trel=trel0

        select case(flag_model)
          case(-1) ! for generic use by users
            if(myrank==0) write(16,*)'Generic tracer transport model'
          case(0) ! Tracer age
            !Method: all i.c. =0; conc=1 at relevant bnd(s), and itrtype=0 at ocean bnd 
            if(myrank==0) write(16,*)'tracer age calculation'
            if(mod(ntracers,2)/=0) call parallel_abort('STEP: ntracers must be even')
          case(1) ! Sediment model

!LLP
#ifndef USE_SED
            call parallel_abort('MAIN: need to turn on USE_SED')
#else
            if(ntracers<=0) call parallel_abort('MAIN: ntracers must be > 0 for sediment model')
            if(ics==2) call parallel_abort('MAIN: Sediment model cannot be used with lat/long coordinates (ics=2)')
            if(imm/=0) call parallel_abort('MAIN: imm and sediment model cannot be used at same time')
!'
! * FG. - Moving most of sediment initializations within sed_init.F90
!           Reads sediment model inputs (sediment.in file)
            CALL read_sed_input
!           Allocation of sediment arrays
            CALL sed_alloc
#endif /*USE_SED*/
!LLP end

          case(2) ! EcoSim
            if(myrank==0) write(16,*) 'Ecological model invoked'

#ifndef USE_ECO
              call parallel_abort('MAIN: need to turn on USE_ECO')
#else

!             Initialize tracers indices and some scalars
              call initialize_param
              call initialize_scalars

              if(myrank==0) write(16,*)'Numbert of Biological Tracers (NBIT)=', NBIT
!'

!             and converts to Julian day 
              if(month==1)then
                yday = day
              else if(month==2)then
                yday = day + 31
              else if(month==3)then
                yday = day + 59
              else if(month==4)then
                yday = day + 90
              else if(month==5)then
                yday = day + 120
              else if(month==6)then
                yday = day + 151
              else if(month==7)then
                yday = day + 181
              else if(month==8)then
                yday = day + 212
              else if(month==9)then
                yday = day + 243
              else if(month==10)then
                yday = day + 273
              else if(month==11)then
                yday = day + 304
              else
                yday = day + 334
              endif

!...          Writes tracers indices for output
              if(myrank==0) then
                open(31, file='ecosim.out', status='replace')
                write(31,*) 'Ecological model output'
                write(31,*) 'Output identifiers'
                write(31,*) 'itemp', itemp
                write(31,*) 'isalt', isalt
                write(31,*) 'Nutrients and DIC'
                write(31,*) '  iDIC_', iDIC_
                if(IRON==1) write(31,*) '  iFeO_', iFeO_
                write(31,*) '  iNH4_', iNH4_
                write(31,*) '  iNO3_', iNO3_
                write(31,*) '  iPO4_', iPO4_
                write(31,*) '  iSiO_', iSiO_
                write(31,*) 'Bacterioplankton group'
                do i=1,Nbac
                  write(31,*) 'Nbac ', i
                  write(31,*) '  iBacC', iBacC(i)
                  if(IRON==1) write(31,*) '  iBacF', iBacF(i)
                  write(31,*) '  iBacN', iBacN(i)
                  write(31,*) '  iBacP', iBacP(i)
                enddo
                write(31,*) 'Dissolved organic matter'
                write(31,*)'Ndom ', 1
                if(CDOC==1) write(31,*)'  iCDMC', iCDMC(1)
                write(31,*)'  iDOMC', iDOMC(1)
                write(31,*)'  iDOMN', iDOMN(1)
                write(31,*)'  iDOMP', iDOMP(1)
                if(Ndom==2) then
                  write(31,*)'Ndom ', 2
                  if(CDOC==1)  write(31,*)'  iCDMC', iCDMC(2)
                  write(31,*)'  iDOMC', iDOMC(2)
                  write(31,*)'  iDOMN', iDOMN(2)
                endif
                write(31,*) 'Fecal particulate matter'
                do i=1,Nfec
                  write(31,*)'Nfec ', i
                  write(31,*)'  iFecC', iFecC(i)
                  if(IRON==1) write(31,*)'  iFecF', iFecF(i)
                  write(31,*)'  iFecN', iFecN(i)
                  write(31,*)'  iFecP', iFecP(i)
                  write(31,*)'  iFecS', iFecS(i)
                enddo
                write(31,*) 'Phytoplankton group'
                do i=1,Nphy
                  write(31,*)'Nphy ', i
                  write(31,*)'  iPhyC',iPhyC(i)
                  if(IRON==1) write(31,*)'  iPhyF',iPhyF(i)
                  write(31,*)'  iPhyN',iPhyN(i)
                  write(31,*)'  iPhyP',iPhyP(i)
                  if(PHY(i)<=2) write(31,*)'  iPhyS', iPhyS(i)
                  do j=1,Npig
                    if(PIG(PHY(i),j)==1) write(31,*) '  iPigs',j, iPigs(i,j)
                  enddo
                enddo
                write(31,*)'Zooplankton group'
                do i=1,Nzoo
                  write(31,*)'Nzoo ', i
                  write(31,*)'  iZooC',izooC(i)
                  write(31,*)'  iZooN',izooN(i)
                  write(31,*)'  iZooP', izooP(i)
                enddo
                write(31,*)'Oxygen'
                write(31,*)'  iDO_', iDO_
                write(31,*)'  iCOD_', iCOD_
                close(31)
              endif !myrank==0

!...          Reads ecological model inputs (ecosim.in file)
              if(myrank==0) write(16,*) 'Reading ecological model parameters inputs...'
!'
              call read_ecoin
              call initialize_biology

!...          Reads constant atmospheric parameters (!MFR - to be used when nws=0... to clean later...)
              if(nws/=0) then
                call parallel_abort('EcoSim must use nws=0 currently')
              else !nws=0
                open(31,file='atmos.in', status='old')
                if(myrank==0) write(16,*) 'Reading atmospheric parameters from atmos.in...'
!'       
                read(31,*)(swild(i),i=1,6) !Uwind(1),Vwind(1),Pair(1),Hair(1),Tair(1),cloud(1)
                Uwind=swild(1)
                Vwind=swild(2)
                Pair=swild(3)
                Hair=swild(4)
                Tair=swild(5)
                cloud=swild(6)
                close(31)
              endif    
#endif /*USE_ECO*/
          case(3) ! Oil spill
            call parallel_abort('Oil spill model under consruction')
          case(4) ! NAPZD Spitz
#ifndef USE_NAPZD
              call parallel_abort('Need to turn on USE_NAPZD')
#else
              if(myrank==0) write(16,*)'reading inputs from NAPZD model'
              call read_napzd_input
#endif /*USE_NAPZD*/
          case(5) ! ICM
            if(myrank==0) write(16,*) 'ICM model invoked'

#ifndef USE_ICM
              call parallel_abort('MAIN: need to turn on USE_ICM')
#else

!             Initialize tracers indices and some scalars
              if(sim_month==1)then
                yday = sim_day
              else if(sim_month==2)then
                yday = sim_day + 31
              else if(sim_month==3)then
                yday = sim_day + 59
              else if(sim_month==4)then
                yday = sim_day + 90
              else if(sim_month==5)then
                yday = sim_day + 120
              else if(sim_month==6)then
                yday = sim_day + 151
              else if(sim_month==7)then
                yday = sim_day + 181
              else if(sim_month==8)then
                yday = sim_day + 212
              else if(sim_month==9)then
                yday = sim_day + 243
              else if(sim_month==10)then
                yday = sim_day + 273
              else if(sim_month==11)then
                yday = sim_day + 304
              else
                yday = sim_day + 334
              endif

!...          Reads ecological model inputs (ecosim.in file)
              if(myrank==0) write(16,*)'Reading ICM parameters inputs'

!!YC02062013  call WQCO1(dt,rnday,NDTWQ)
              allocate(WSRP(nea),WSLP(nea),WSPB1(nea),WSPB2(nea),WSPB3(nea),turb(nea),WRea(nea),stat=istat)  !added by YC
              if(istat/=0) call parallel_abort('Failed to allocate (11)')
              call WQCO2(WSRP,WSLP,WSPB1,WSPB2,WSPB3,turb,WRea) !added by YC
              call WQinput !(time) !added by YC, still need debuging
              call wqm_out 
              if(myrank==0) write(16,*)'done reading ICM parameters'
#endif /*USE_ICM*/
          case(6) ! TIMOR
#ifndef USE_TIMOR
              call parallel_abort('Need to turn on USE_TIMOR')
#else
!             Init. TIMOR (tr_nd)
              call init_flmud
#endif /*USE_TIMOR*/

          case default
            call parallel_abort('Unknown tracer model')
        end select !flag_model

!       i.c.
        select case(flag_ic)
          case(1)                
!	    Horizontally varying
            do m=1,ntracers
              write(ifile_char,'(i03)')m
              ifile_char=adjustl(ifile_char); ifile_len=len_trim(ifile_char)
              inputfile='htr_'//ifile_char(1:ifile_len)//'.ic'
              open(10,file=inputfile,status='old')
              read(10,*)
              read(10,*) !np
              do j=1,np_global
                read(10,*) num,xtmp,ytmp,tr_tmp1
                if(tr_tmp1<0) then
                  write(errmsg,*)'Initial invalid tracer at:',j,tr_tmp1
                  call parallel_abort(errmsg)
                endif
                if(ipgl(j)%rank==myrank) tr_nd(m,:,ipgl(j)%id)=tr_tmp1
              enddo !j
              close(10)
            enddo !m
          case(2)
!	    Vertically varying
            do m=1,ntracers
              write(ifile_char,'(i03)')m
              ifile_char=adjustl(ifile_char) 
              ifile_len=len_trim(ifile_char)
              inputfile='vtr_'//ifile_char(1:ifile_len)//'.ic'
              open(10,file=inputfile,status='old')
              read(10,*)nz_r2
              if(nz_r2<2) then
                write(errmsg,*)'Change nz_r2:',nz_r2
                call parallel_abort(errmsg)
              endif
              deallocate(swild10,swild,swild2,stat=istat)
              allocate(z_r2(nz_r2),stat=istat)
              allocate(swild(max(nsa+nvrt+12+ntracers,nz_r2)),swild2(max(nvrt,nz_r2),12), &
     &swild10(max(3,nvrt,nz_r2),12),stat=istat)
              do k=1,nz_r2
                read(10,*)j,z_r2(k),swild10(k,1)
                if(swild10(k,1)<0) then
                  write(errmsg,*)'Initial invalid Tr at',k,swild10(k,1)
                  call parallel_abort(errmsg)
                endif
                if(k>=2) then; if(z_r2(k)<=z_r2(k-1)) then
                  write(errmsg,*)'Inverted z-level (10):',k
                  call parallel_abort(errmsg)
                endif; endif
              enddo !k
              close(10)

!             Cubic spline coefficients
              call cubic_spline(nz_r2,z_r2,swild10(1:nz_r2,1),0._rkind,0._rkind,swild)
              swild2(1:nz_r2,1)=swild(1:nz_r2)

              do i=1,npa
                if(kbp(i)==0) cycle
!               Wet node
                !do k=kbp(i)+1,nvrt
                !  swild(k)=(ze(k,i)+ze(k-1,i))/2
                !enddo !k
                call eval_cubic_spline(nz_r2,z_r2,swild10(1:nz_r2,1),swild2(1:nz_r2,1), &
     &nvrt-kbp(i)+1,znl(kbp(i):nvrt,i),0,z_r2(1),z_r2(nz_r2),tr_nd(m,kbp(i):nvrt,i))

! 	        Extend
                do k=1,kbp(i)-1
                  tr_nd(m,k,i)=tr_nd(m,kbp(i),i)
                enddo !k
              enddo !i=1,npa

              deallocate(z_r2)
            enddo !m=1,ntracers

          case(0)
!	        Model sets own i.c.
#ifdef USE_ECO
            call bio_init(trel0)
#endif

#ifdef USE_TIMOR
            !Already init'ed in init_flmud
            !tr_el(1:ntracers,:,1:npa)=tr_nd
#endif
          case default
            call parallel_abort('INIT: Check  flag_ic!!!')
        end select !flag_ic

!       Initialize Tracer
        ltmp=.false.
        if(flag_ic/=0) ltmp=.true.
#ifdef USE_TIMOR
        ltmp=.true.
#endif
        !if(flag_ic/=0)then
        if(ltmp)then
          do m=1,ntracers
            do i=1,nea
              n1=elnode(1,i)
              n2=elnode(2,i)
              n3=elnode(3,i)
              do k=2,nvrt
                trel0(m,k,i)=(tr_nd(m,k,n1)+tr_nd(m,k,n2)+tr_nd(m,k,n3)+ &
     & tr_nd(m,k-1,n1)+tr_nd(m,k-1,n2)+tr_nd(m,k-1,n3))/6
              enddo !k
              trel0(m,1,i)=trel0(m,2,i) !mainly for hotstart format
            enddo !i
          enddo !m
        endif !ltmp

        if(irouse_test==1) then
          trel0=0
          trel0(:,1:2,:)=1
        endif

        trel=trel0

        if(myrank==0) write(16,*)'done init. tracers..'
      endif !ntracers>0
!     end user-defined tracer part

!...  Initialize GOTM for both cold and hot starts (for cde etc).
!...  For real hot start, q2, xl, dfv and dfh will use the values in hotstart.in;
!...  otherwise they will be assigned values below.
      if(itur==4) then
#ifdef USE_GOTM
          call init_turbulence(8,'gotmturb.inp',nvrt-1) !GOTM starts from level 0
          call init_tridiagonal(nvrt-1)
#endif
      endif

#ifdef USE_SED2D
#ifdef INCLUDE_TIMING
      cwtmp2=mpi_wtime() !start of timer
#endif

      call sed2d_init

#ifdef INCLUDE_TIMING
      timer_ns(2)=timer_ns(2)+mpi_wtime()-cwtmp2 !end timing this section
#endif 
#endif /*USE_SED2D*/

      if(myrank==0) write(16,*)'done initializing cold start'
      
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Hot start section
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      time=0.d0
      if(ihot/=0) then
        open(36,file='hotstart.in',form='unformatted',status='old')
        read(36) double1,iths,ifile
        time=double1

        ! Element data
        do i=1,ne_global
          read(36) iegb,itmp1,((swild4(j,l),l=1,3+2*ntracers),j=1,nvrt)
          if(iegl(iegb)%rank==myrank) then
            ie=iegl(iegb)%id
            idry_e(ie)=itmp1
            we(:,ie)=swild4(:,1)
            tsel(1,:,ie)=swild4(:,2)
            tsel(2,:,ie)=swild4(:,3)
            do l=1,ntracers
!              if(flag_model==0) cycle !use i.c. for testing model
#ifndef USE_NAPZD
              trel0(l,:,ie)=swild4(:,3+2*l-1)
              trel(l,:,ie)=swild4(:,3+2*l)
#endif
            enddo !l
          endif
        enddo !i=1,ne_global

        ! Side data
        do i=1,ns_global
          read(36) isgb,itmp1,((swild10(j,l),l=1,4),j=1,nvrt)
          if(isgl(isgb)%rank==myrank) then
            iside=isgl(isgb)%id
            idry_s(iside)=itmp1
            su2(:,iside)=swild10(:,1)
            sv2(:,iside)=swild10(:,2)
            tsd(:,iside)=swild10(:,3)
            ssd(:,iside)=swild10(:,4)
          endif
        enddo !i=1,ns_global

        ! Node data (include non-hydrostatic pressure even for hydrostatic model for convenience)
        do i=1,np_global
          read(36) ipgb,double1,itmp,((swild10(j,l),l=1,11),j=1,nvrt)
          if(ipgl(ipgb)%rank==myrank) then
            ip=ipgl(ipgb)%id
            eta2(ip)=double1
            idry(ip)=itmp
            tnd(:,ip)=swild10(:,1)
            snd(:,ip)=swild10(:,2)
            tem0(:,ip)=swild10(:,3)
            sal0(:,ip)=swild10(:,4)
            q2(:,ip)=swild10(:,5)
            xl(:,ip)=swild10(:,6)
            dfv(:,ip)=swild10(:,7)
            dfh(:,ip)=swild10(:,8)
            dfq1(:,ip)=swild10(:,9)
            dfq2(:,ip)=swild10(:,10)
            qnon(:,ip)=swild10(:,11)
          endif
        enddo !i=1,np_global

#ifdef USE_ICM
        do i=1,ne_global       
          read(36) iegb,swild3(1:40)
          if(iegl(iegb)%rank==myrank) then
            ie=iegl(iegb)%id
            sed_BENDO(ie)=swild3(1)
            CTEMP(ie)=swild3(2)
            BBM(ie)=swild3(3)
            CPOS(ie)=swild3(4)
            PO4T2TM1S(ie)=swild3(5)
            NH4T2TM1S(ie)=swild3(6)
            NO3T2TM1S(ie)=swild3(7)
            HST2TM1S(ie)=swild3(8)
            CH4T2TM1S(ie)=swild3(9)
            CH41TM1S(ie)=swild3(10)
            SO4T2TM1S(ie)=swild3(11)
            SIT2TM1S(ie)=swild3(12)
            BENSTR1S(ie)=swild3(13)
            CPOP(ie,1)=swild3(14)
            CPOP(ie,2)=swild3(15)
            CPOP(ie,3)=swild3(16)
            CPON(ie,1)=swild3(17)
            CPON(ie,2)=swild3(18)
            CPON(ie,3)=swild3(19)
            CPOC(ie,1)=swild3(20)
            CPOC(ie,2)=swild3(21)
            CPOC(ie,3)=swild3(22)
            NH41TM1S(ie)=swild3(23)
            NO31TM1S(ie)=swild3(24)
            HS1TM1S(ie)=swild3(25)
            SI1TM1S(ie)=swild3(26)
            PO41TM1S(ie)=swild3(27)
            PON1TM1S(ie)=swild3(28)
            PON2TM1S(ie)=swild3(29)
            PON3TM1S(ie)=swild3(30)
            POC1TM1S(ie)=swild3(31)
            POC2TM1S(ie)=swild3(32)
            POC3TM1S(ie)=swild3(33)
            POP1TM1S(ie)=swild3(34)
            POP2TM1S(ie)=swild3(35)
            POP3TM1S(ie)=swild3(36)
            PSITM1S(ie)=swild3(37)
            BFORMAXS(ie)=swild3(38)
            ISWBENS(ie)=swild3(39)
            DFEEDM1S(ie)=swild3(40)

            !write(12,*)'ICM:',iegb,swild3(1:40)
          endif
        enddo !i
#endif

#ifdef USE_SED2D 
        do i=1,np_global
          read(36) ipgb,double1
          if(ipgl(ipgb)%rank==myrank) then
            ip=ipgl(ipgb)%id
            dp(ip)=double1
          endif
        enddo !i
#endif

#ifdef USE_SED
        do i=1,np_global
          read(36) ipgb,double1
          read(36) ipgb,tmp1
          if(ipgl(ipgb)%rank==myrank) then
            ip=ipgl(ipgb)%id
            dp(ip)=double1
            rough_p(ip)=tmp1
          endif
        enddo !i

        deallocate(swild4)
        allocate(swild4(Nbed,ntracers),swild99(Nbed,MBEDP),stat=istat)
        if(istat/=0) call parallel_abort('MAIN: alloc failure (55)')
        do i=1,ne_global
          read(36) iegb,swild4,swild99
          if(iegl(iegb)%rank==myrank) then
            ie=iegl(iegb)%id
            bed_frac(:,ie,:)=swild4
            bed(:,ie,:)=swild99
          endif
        enddo !i
        deallocate(swild99)

#endif /*USE_SED*/

#ifdef USE_HA
!...
!......HOT START INFORMATION FOR HARMONIC ANALYSIS
!......Adapted from ADCIRC
!...
        IF(iharind.EQ.1) THEN
          IHABEG=ITHAS+NHAINC
!...
!.......IF HARMONIC ANALYSIS HAS NOT BEGUN, COLD START THE HARMONIC ANALYSIS
!...
          IF(iths.LT.IHABEG) THEN
            ICHA=0
            CALL HACOLDS(HAFREQ)
            IF(NHAGE.EQ.1) CALL HACOLDSEG(np)
            IF(NHAGV.EQ.1) CALL HACOLDSVG(np)
            IF (CHARMV) THEN
              ELAV  =0.D0
              XVELAV=0.D0
              YVELAV=0.D0
              ELVA  =0.D0
              XVELVA=0.D0
              YVELVA=0.D0
!              DO I=1,np
!                ELAV(I)=0.D0
!                XVELAV(I)=0.D0
!                YVELAV(I)=0.D0
!                ELVA(I)=0.D0
!                XVELVA(I)=0.D0
!                YVELVA(I)=0.D0
!              END DO
            ENDIF   !   charmv
          ENDIF

!...
!........IF HARMONIC ANALYSIS HAS ALREADY BEGUN, READ IN HOT START
!........HARMONIC ANALYSIS, MEAN AND SQUARE INFO
!...
          IF(iths.GT.ITHAS) READ(36) ICHA
          IF(iths.GE.IHABEG) THEN
            CALL HAHOTS(0,0,np_global,0,0,NHAGE,NHAGV,36,myrank)
            IF(NHAGE.EQ.1) CALL HAHOTSEG(np_global,36,myrank)
            IF(NHAGV.EQ.1) CALL HAHOTSVG(np_global,36,myrank)
          ENDIF

!..Read in Means and Squares

          IF(CHARMV) THEN
            IF((FMV.NE.0.).AND.(iths.GT.ITMV)) THEN
              READ(36) NTSTEPS
              IF(NHAGE.EQ.1) THEN
                DO I=1,np_global
                  READ(36) tmp1
                  READ(36) tmp2
                  if (ipgl(i)%rank==myrank) then	!for augmented region
                    il = ipgl(i)%id
                    ELAV(il)=tmp1
                    ELVA(il)=tmp2
                  endif
                ENDDO
              ENDIF
              IF(NHAGV.EQ.1) THEN
                DO I=1,np_global
                  READ(36) tmp1
                  READ(36) tmp2
                  READ(36) tmp3
                  READ(36) tmp4
                  if (ipgl(i)%rank==myrank) then	!for augmented region
                    il = ipgl(i)%id
                    XVELAV(il) = tmp1
                    YVELAV(il) = tmp2
                    XVELVA(il) = tmp3
                    YVELVA(il) = tmp4
                  endif
                ENDDO
              ENDIF
            ENDIF
          ENDIF   !  charmv
        ENDIF    !HARIND
#endif /*USE_HA*/

        ! Close hotstart file
        close(36)

        if(itur==3) then
          do i=1,npa
            do j=1,nvrt
              q2(j,i)=max(q2min,q2(j,i))
              xl(j,i)=max(xlmin2(i),xl(j,i))
            enddo
          enddo
        endif

!...    change time and iteration for forecast mode
!...    Causion: this affects all t.h. files (fort.5[0-3]) and wind files
        if(ihot==1) then
          time=0
          iths=0
        endif

        if(myrank==0) write(16,*)'hot start at time=',time,iths,' ; stack #=',ifile

!'
!...  find position in the wind input file for nws=1,2, and read in wind[x,y][1,2]
!...  Wind vector always in lat/lon frame
        if(nws==1) then
          open(22,file='wind.th',status='old')
          rewind(22)
          ninv=time/wtiminc
          wtime1=ninv*wtiminc 
          wtime2=(ninv+1)*wtiminc 
          do it=0,ninv
            read(22,*)tmp,wx1,wy1
          enddo
          read(22,*)tmp,wx2,wy2
          windx1=wx1
          windy1=wy1
          windx2=wx2
          windy2=wy2
        endif

        if(nws==4) then
          open(22,file='wind.th',status='old')
          rewind(22)
          ninv=time/wtiminc
          wtime1=ninv*wtiminc
          wtime2=(ninv+1)*wtiminc
          do it=0,ninv
            read(22,*)tmp,rwild(:,:)
          enddo
          do i=1,np_global
            if(ipgl(i)%rank==myrank) then
              nd=ipgl(i)%id
              windx1(nd)=rwild(i,1)
              windy1(nd)=rwild(i,2)
              pr1(nd)=rwild(i,3)
            endif
          enddo !i
          read(22,*)tmp,rwild(:,:)
          do i=1,np_global
            if(ipgl(i)%rank==myrank) then
              nd=ipgl(i)%id
              windx2(nd)=rwild(i,1)
              windy2(nd)=rwild(i,2)
              pr2(nd)=rwild(i,3)
            endif
          enddo !i
        endif !nws=4

        if(nws>=2.and.nws<=3) then
          ninv=time/wtiminc
          wtime1=ninv*wtiminc 
          wtime2=(ninv+1)*wtiminc 
          call get_wind(wtime1,windx1,windy1,pr1,airt1,shum1)
          call get_wind(wtime2,windx2,windy2,pr2,airt2,shum2)

!         Overwrite wind with wind.th
!          open(22,file='wind.th',status='old')
!          rewind(22)
!          do it=0,iths
!            read(22,*)tmp,wx1,wy1
!          enddo !it
!          windx1=wx1; windx2=wx1
!          windy1=wy1; windy2=wy1
!         End

          if(nws==3) then
!           Read sflux.th
            open(23,file='sflux.th',status='old')
            rewind(23)
            do it=0,iths
              read(23,*)
            enddo
            read(23,*)tmp,fluxsu00,srad00
          endif !nws==3
        endif !nws

!       VIMS Point source loading added by YC
#ifdef USE_ICM 
        if(myrank==0) write(16,*)'hotstart ICM point source..'
        call WQCO1(dt,rnday,NDTWQ)                                         !added by YC
        allocate(WSRP(nea),WSLP(nea),WSPB1(nea),WSPB2(nea),WSPB3(nea),turb(nea),WRea(nea),stat=istat)  !added by YC
        call WQCO2(WSRP,WSLP,WSPB1,WSPB2,WSPB3,turb,WRea)                  !added by YC
        if(iWQPS==2) then
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
!
          x1 = 1.0E3 !from kg to g
          npstiminc=86400.
          open(61,file='ps.in',status='old')
          rewind(61)
          ninv=time/npstiminc
          npstime=ninv*npstiminc
!org yc        npstime1=ninv*npstiminc
!org yc        npstime2=(ninv+1)*npstiminc
          read(61,*)      !title
          do it=0,ninv
            read(61,*)    !time
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
              if(myrank==0)write(16,*)'ICM: xPSK= ',xPSK
            enddo !i
          enddo !ninv
        endif ! iWQPS=2

        npstiminc=86400.
        ninv=time/npstiminc
        npstime=ninv*npstiminc
 
        do it=0,ninv
          call WQinput !(time)
        enddo
        call wqm_out

        if(myrank==0) write(16,*)'end hotstart ICM point source..'
        
!    VIMS surface temperature mode added by YC02062013
        if(iSun==2) then
          if(myrank==0) write(16,*)'doing ICM hotstart surface T..'
          open(62,file='surface_t.th',status='old')
          rewind(62)
          read(62,*)
          read(62,*)
          do it=0,ninv
            read(62,*)
            do i=1,np_global
              read(62,*) ipgb,xtmp
              if(ipgl(ipgb)%rank==myrank) then
                surf_t1(ipgl(ipgb)%id)=xtmp
              endif
            enddo
          enddo
     
          read(62,*)
          do i=1,np_global
            read(62,*) ipgb,xtmp
            if(ipgl(ipgb)%rank==myrank) then
             surf_t2(ipgl(ipgb)%id)=xtmp
            endif
          enddo

          surf_time1=int(time/npstiminc)*npstiminc  !added by wangzg
          surf_time2=surf_time1+86400.
        
        endif ! iSun=2
        if(myrank==0) write(16,*)'done ICM hotstart surface T..'

#endif /*USE_ICM*/
!End of VIMS Point source loading added by YC

!...    Nudging 
        if(inu_st==2) then
          irec_nu=time/step_nu+1
          time_nu=irec_nu*step_nu
       
          do it=1,irec_nu+1
            read(37)floatout
            read(35)floatout
            do i=1,np_global
              read(37)(swild8(j,1),j=1,nvrt)
              read(35)(swild8(j,2),j=1,nvrt)
              if(it==irec_nu.and.ipgl(i)%rank==myrank) then
                tnd_nu1(:,ipgl(i)%id)=swild8(1:nvrt,1)
                snd_nu1(:,ipgl(i)%id)=swild8(1:nvrt,2)
              endif
              if(it==irec_nu+1.and.ipgl(i)%rank==myrank) then
                tnd_nu2(:,ipgl(i)%id)=swild8(1:nvrt,1)
                snd_nu2(:,ipgl(i)%id)=swild8(1:nvrt,2)
              endif
            enddo !i
          enddo !it
          irec_nu=irec_nu+1
        endif !inu_st

        if(ntracers>0.and.inu_tr==2) then
          irec_nu_tr=time/step_nu_tr+1
          time_nu_tr=irec_nu_tr*step_nu_tr
       
          do it=1,irec_nu_tr+1
            read(45)floatout
            do i=1,np_global
              read(45)swild9
              if(it==irec_nu_tr.and.ipgl(i)%rank==myrank) then
                trnd_nu1(:,:,ipgl(i)%id)=swild9
              endif
              if(it==irec_nu_tr+1.and.ipgl(i)%rank==myrank) then
                trnd_nu2(:,:,ipgl(i)%id)=swild9
              endif
            enddo !i
          enddo !it
          irec_nu_tr=irec_nu_tr+1
        endif !inu_tr

!...    Find positions in t.h. files 
        if(nettype>0) then
          rewind(50)
          ninv=time/th_dt(1,1)
          do it=0,ninv
            read(50,*) ttt,ath(1:nettype,1,1,1)
          enddo
          th_time(1,1,1)=ttt
          read(50,*) ttt,ath(1:nettype,1,2,1)
          th_time(1,2,1)=ttt
        endif !nettype

        if(nfltype>0) then 
          rewind(51)
          ninv=time/th_dt(1,2)
          do it=0,ninv
            read(51,*) ttt,ath(1:nfltype,1,1,2)
          enddo 
          th_time(1,1,2)=ttt
          read(51,*) ttt,ath(1:nfltype,1,2,2)
          th_time(1,2,2)=ttt
        endif !nfltype

        if(ntetype>0) then 
          rewind(52)
          ninv=time/th_dt(1,3)
          do it=0,ninv
            read(52,*) ttt,ath(1:ntetype,1,1,3)
          enddo 
          th_time(1,1,3)=ttt
          read(52,*) ttt,ath(1:ntetype,1,2,3)
          th_time(1,2,3)=ttt
        endif

        if(nsatype>0) then 
          rewind(53)
          ninv=time/th_dt(1,4)
          do it=0,ninv
            read(53,*) ttt,ath(1:nsatype,1,1,4)
          enddo 
          th_time(1,1,4)=ttt
          read(53,*) ttt,ath(1:nsatype,1,2,4)
          th_time(1,2,4)=ttt
        endif !nsatype

        if(ntrtype>0) then
          do m=1,ntracers
            rewind(300+m)
            ninv=time/th_dt(m,5)
            do it=0,ninv
              read(300+m,*) ttt,ath(1:ntrtype,m,1,5)
            enddo
            th_time(m,1,5)=ttt
            read(300+m,*) ttt,ath(1:ntrtype,m,2,5)
            th_time(m,2,5)=ttt
          enddo
        endif 

        if(nettype2>0) then
          ninv=time/th_dt2(1)
          th_time2(1,1)=ninv*th_dt2(1)
          th_time2(2,1)=th_time2(1,1)+th_dt2(1)
          read(54,rec=ninv+1)floatout,ath2(1,1,1:nnode_et,1,1)
          read(54,rec=ninv+2)floatout,ath2(1,1,1:nnode_et,2,1)
          irec_th(1)=ninv+2
        endif !nettype2

        if(nfltype2>0) then
          ninv=time/th_dt2(2)
          th_time2(1,2)=ninv*th_dt2(2)
          th_time2(2,2)=th_time2(1,2)+th_dt2(2)
          read(58,rec=ninv+1)floatout,ath2(1:2,1:nvrt,1:nnode_fl,1,2)
          read(58,rec=ninv+2)floatout,ath2(1:2,1:nvrt,1:nnode_fl,2,2)
          irec_th(2)=ninv+2
        endif !nfltype2

        if(ntetype2>0) then
          ninv=time/th_dt2(3)
          th_time2(1,3)=ninv*th_dt2(3)
          th_time2(2,3)=th_time2(1,3)+th_dt2(3)
          read(56,rec=ninv+1)floatout,ath2(1,1:nvrt,1:nnode_te,1,3)
          read(56,rec=ninv+2)floatout,ath2(1,1:nvrt,1:nnode_te,2,3)
          irec_th(3)=ninv+2
        endif !ntetype2

        if(nsatype2>0) then
          ninv=time/th_dt2(4)
          th_time2(1,4)=ninv*th_dt2(4)
          th_time2(2,4)=th_time2(1,4)+th_dt2(4)
          read(57,rec=ninv+1)floatout,ath2(1,1:nvrt,1:nnode_sa,1,4)
          read(57,rec=ninv+2)floatout,ath2(1,1:nvrt,1:nnode_sa,2,4)
          irec_th(4)=ninv+2
        endif !nsatype2

        if(ntrtype2>0) then
          ninv=time/th_dt2(5)
          th_time2(1,5)=ninv*th_dt2(5)
          th_time2(2,5)=th_time2(1,5)+th_dt2(5)
          read(59,rec=ninv+1)floatout,ath2(1:ntracers,1:nvrt,1:nnode_tr,1,5)
          read(59,rec=ninv+2)floatout,ath2(1:ntracers,1:nvrt,1:nnode_tr,2,5)
          irec_th(5)=ninv+2
        endif !ntrtype2

        !if(ihydraulics/=0.and.nhtblocks>0) then; do it=1,iths; read(49,*) ttt,tmp; enddo; endif;

        if(if_source==1) then
          if(nsources>0) then
            ninv=time/th_dt3(1)
            rewind(63)
            do it=0,ninv
              read(63,*)tmp,ath3(1:nsources,1,1,1)
            enddo !it
            th_time3(1,1)=tmp
            read(63,*)tmp,ath3(1:nsources,1,2,1)
            th_time3(2,1)=tmp

            ninv=time/th_dt3(3)
            rewind(65)
            do it=0,ninv
              !do j=1,2+ntracers
              read(65,*)tmp,ath3(1:nsources,1:2+ntracers,1,3)
              !enddo !j
            enddo !it
            th_time3(1,3)=tmp
            !do j=1,2+ntracers
            read(65,*)tmp,ath3(1:nsources,1:2+ntracers,2,3)
            !enddo !j
            th_time3(2,3)=tmp
          endif !nsources
     
          if(nsinks>0) then
            ninv=time/th_dt3(2)
            rewind(64)
            do it=0,ninv
              read(64,*)tmp,ath3(1:nsinks,1,1,2)
            enddo !it
            th_time3(1,2)=tmp
            read(64,*)tmp,ath3(1:nsinks,1,2,2)
            th_time3(2,2)=tmp
          endif !nsinks
        endif !if_source

!       Station output
        if(iout_sta/=0.and.myrank==0) then
          do i=1,nvar_sta
            rewind(250+i)    
            do it=1,iths
              if(iof_sta(i)==1.and.mod(it,nspool_sta)==0) then
                read(250+i,*)
!                if(i>4) read(250+i,*) !vertical profile for 3D var
              endif
            enddo !it
          enddo !i
        endif !myrank

!...  end hot start section
      endif !ihot/=0

#ifdef USE_SED
!...    Sediment model initialization
        call sed_init
#endif /*USE_SED*/

!       Initialize time series for hydraulic structures that use them, including 
!       opening files and "fast forwarding" to the restart time
        
      if(ihydraulics/=0.and.nhtblocks>0) then
         call init_struct_time_series(time)
      endif

!...  init. eta1 (for some routines like WWM) and i.c. (for ramp function)
      eta1=eta2 
      etaic=eta2

      if(myrank==0) write(16,'(a)')'Done initializing variables'

!     NZPZD model: overwrite NO3
#ifdef USE_NAPZD
!CSD HARDWIRE some code here to set up an initial NO3 distribution for a cold start initialization.
!CSD Let's initialize the NO3 (1) field using an analytical function with high nitrate concentration
!CSD uniformly with depth in the estuary...melding into our oceanic initial profile for NO3 offshore.
!CSD  call mpi_barrier(comm,ierr)
        do i=1,nea
          if(xctr(i)/1000>=334 .and.                           &
     &      yctr(i)/1000<=7.5/5.0*xctr(i)/1000-206.5 .and.   &
     &      yctr(i)/1000>=292) then
!CSD temporarily move the nitrate to x coordinate 392.
            ft1=0.5+0.5*tanh((xctr(i)/1000-392)/4)
            trel0(1,:,i)=20.0*ft1+trel0(1,:,i)*(1.0-ft1);
          else if(xctr(i)/1000>=340.3 .and. yctr(i)/1000<=7.5/5.*xctr(i)/1000-206.5        &
     &            .and. yctr(i)/1000>=-13.5/19*xctr(i)/1000+528.3) then
!           ft1=0.5+0.5*tanh((xnd(i)/1000-338)/4)
!CSD temporarily move the nitrate to x coordinate 392.
            ft1=0.5+0.5*tanh((xctr(i)/1000-338)/4);
            trel0(1,:,i)=20.0*ft1+trel0(1,:,i)*(1.0-ft1);
          endif
          trel(1,:,i)=trel0(1,:,i)
        enddo !i=1,nea
#endif /*USE_NAPZD*/

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Open output files
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

!     Open global output files and write header data
      if(ihot<=1) ifile=1 !reset output file #
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

!      if(iwrite==0) then !one global output
!      call write_header0(iths,iths)
!  elseif(iwrite==1) then !each processor output a separate file
!    call write_ggtbls
!      call write_header1
!  elseif(iwrite==2) then !separate files with sequential unformatted
!    call write_ggtbls
!    call write_header2(nt,iths,iths)
!  endif

!     Write local to global mapping and header info for combining scripts
      fdb='local_to_global_0000'
      lfdb=len_trim(fdb)
      write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
      open(10,file='outputs/'//fdb,status='replace')

!     header info (except ivs, i23d, variable_nm, and variable_dim; the 1st two can be inferred from output names and 
!     last two are not important)
      write(10,*)ne_global,np_global,nvrt,nproc,ntracers !global info
!     local to global mapping      
      write(10,*)'local to global mapping:'
      write(10,*)ne
      do ie=1,ne
        write(10,*)ie,ielg(ie)
      enddo
      write(10,*)np 
      do ip=1,np
        write(10,*)ip,iplg(ip)
      enddo
      write(10,*)ns
      do isd=1,ns
        write(10,*)isd,islg(isd)
      enddo

      write(10,*)'Header:'
      write(10,'(a)')data_format,version,start_time
      write(10,*)nrec,real(dt*nspool),nspool,nvrt,kz,real(h0),real(h_s),real(h_c),real(theta_b),real(theta_f)
      write(10,*)(real(ztot(k)),k=1,kz-1),(real(sigma(k)),k=1,nvrt-kz+1)
      if(ics==1) then
        write(10,*)np,ne,(real(xnd(m)),real(ynd(m)),real(dp00(m)),kbp00(m),m=1,np), &
     &             (3,(elnode(mm,m),mm=1,3),m=1,ne)
      else !lat/lon
        write(10,*)np,ne,(real(xlon(m)/pi*180),real(ylat(m)/pi*180),real(dp00(m)),kbp00(m),m=1,np), &
     &             (3,(elnode(mm,m),mm=1,3),m=1,ne)
      endif !ics

      close(10)
      
      if(myrank==0) write(16,'(a)')'Done initializing outputs'

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Preparation for time stepping
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      if(ihot==0) iths=0
!...  Compute initial bed deformation and update depths info
      do i=1,npa
        bdef1(i)=bdef(i)/ibdef*min0(iths,ibdef)
        if(imm==1) then
          !Add conditional to avoid conflict with sediment morph model
          dp(i)=dp00(i)-bdef1(i)
        else if(imm==2) then
          call update_bdef(iths*dt,xnd(i),ynd(i),dep,swild)
          dp(i)=dep !min(1.,7-(xnd(i)+iths*dt))
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
      enddo !i=1,ne

!...  Compute initial vgrid
      if(inunfl==0) then
        call levels0(iths,iths)
      else
        call levels1(iths,iths)
      endif

      if(myrank==0) write(16,*)'done computing initial vgrid...'

!...  Compute nodal vel. 
      call nodalvel !(ifltype)
      if(myrank==0) write(16,*)'done computing initial nodal vel...'

!...  Compute initial density at nodes or elements
      prho=-99
      do i=1,npa
        if(idry(i)==1) cycle
        do k=1,nvrt
          kl=max(k,kbp(i))
          prho(k,i)=eqstate(1,tnd(k,i),snd(k,i)              &
#ifdef USE_SED
     &                      ,tr_nd(1:ntracers,k,i),Srho(:) &
#endif /*USE_SED*/
#ifdef USE_TIMOR
     &                      ,tr_nd(1:ntracers,kl,i),rhomud(1:ntracers,kl,i),laddmud_d &
#endif
     &                      )
        enddo !k
      enddo !i

      if(iupwind_t/=0) then
        erho=-99
        do i=1,nea
          if(idry_e(i)==1) cycle

          do k=1,nvrt
            kl=max(k,kbe(i))
#ifdef USE_TIMOR
            do m=1,ntracers
              swild(m)=sum(rhomud(m,kl,elnode(1:3,i)))/3
            enddo !m
#endif
            erho(k,i)=eqstate(2,tsel(1,k,i),tsel(2,k,i)      &
!LLP
#ifdef USE_SED
     &                        ,trel(:,k,i),Srho(:)         &
#endif /*USE_SED*/
#ifdef USE_TIMOR
     &                        ,trel(:,k,i),swild(1:ntracers),laddmud_d &
#endif
!LLP end
     &                       )
          enddo !k
        enddo !i
      endif !iupwind_t

!...  Compute mean density profile at nodes or elements (using current z-coord.)
!Error: did not account for sediment
      if(ibcc_mean==1.or.ihot==0.and.icst==2) then
        call mean_density
      else !other cases
        rho_mean=0
      endif

      if(myrank==0) write(16,*)'done computing initial density...'

!...  Initialize heat budget model
      if(ihconsv/=0) then
!#ifdef USE_SFLUX
        call surf_fluxes(wtime1,windx1,windy1,pr1,airt1,shum1,srad,fluxsu,fluxlu,hradu,hradd,tauxz,tauyz, &
#ifdef PREC_EVAP
     &                   fluxprc,fluxevp, &
#endif
     &                   nws,fluxsu00,srad00)
!#endif
!       fluxsu: the turbulent flux of sensible heat (upwelling) (W/m^2)
!       fluxlu: the turbulent flux of latent heat (upwelling) (W/m^2)
!       hradu: upwelling infrared (longwave) radiative fluxes at surface (W/m^2)
!       hradd: downwelling infrared (longwave) radiative fluxes at surface (W/m^2)
!       srad: solar radiation (W/m^2)
!       tauxz,tauyz: wind stress (in true E-N direction if ics=2)
        do i=1,npa
          sflux(i)=-fluxsu(i)-fluxlu(i)-(hradu(i)-hradd(i))
        enddo
        if(myrank==0) write(16,*)'heat budge model completes...'
      endif !nws>=2

!...  Assign variables in GOTM for cold starts
      if(itur==4.and.(ihot==0.or.ihot==1.and.nramp==1)) then
#ifdef USE_GOTM
!          call init_turbulence(8,'gotmturb.inp',nvrt-1) !GOTM starts from level 0
!          call init_tridiagonal(nvrt-1)
!         Debug
!          do k=0,nvrt-1
!            write(99,*)k,tke1d(k),L1d(k),num1d(k),nuh1d(k)
!          enddo !i
!          stop

          do j=1,npa
            q2(:,j) = tke1d(0:(nvrt-1))
            xl(:,j) = L1d(0:(nvrt-1))
            dfv(:,j) = min(diffmax(j),max(diffmin(j),num1d(0:(nvrt-1))))
            dfh(:,j) = min(diffmax(j),max(diffmin(j),nuh1d(0:(nvrt-1))))
          enddo !j
#endif
      endif !itur==4 etc

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Begin time stepping
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
#ifdef INCLUDE_TIMING
! End timing init section & begin timing time-stepping section
      wtmp1=mpi_wtime()
      wtimer(1,1)=wtmp1-wtimer(1,1)
      wtimer(2,1)=wtmp1 !time-stepping section
#endif

      if(myrank==0) then
        write(16,*)'time stepping begins...',iths+1,ntime
        call flush(16) ! flush "mirror.out"
      endif

      difnum_max_l2=0 !max. horizontal diffusion number reached by each process (check stability)
      iwbl_itmax=0 !cumulative max. of iterations for WBL (Grant-Madsen formulation) for a rank 

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!      endif !first_call


!-------------------------------------------------------------------------------
!   Initialize wind wave model (WWM)
!-------------------------------------------------------------------------------
#ifdef USE_WWM
        !Init. windx,y for WWM using windx,y1
        windx=windx1
        windy=windy1
        CALL INITIALIZE_WWM
#endif      

!     Broadcast to global module
      iths_main=iths
 
!     Deallocate temp. arrays to release memory
      deallocate(nwild,nwild2,swild,swild2,swild3,swild4,swild8,swild10)
      if(nws==4) deallocate(rwild)
      if(ntracers>0.and.inu_tr==2) deallocate(swild9)

      end subroutine schism_init
