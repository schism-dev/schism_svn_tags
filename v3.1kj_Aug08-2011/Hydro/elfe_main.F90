!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                       !
!                                SELFE:                                                 !
!         A Three-Dimensional Baroclinic Model for Unstructured Grids                   !
!                       MPI Version 3.1kj (June 28, 2011) 
!                                                                                       !
!                 Center for Coastal and Land-Margin Research                           !
!             Department of Environmental Science and Engineering                       !
!                   OGI School of Science and Engineering,                              !
!                     Oregon Health & Science University                                !
!                       Beaverton, Oregon 97006, USA                                    !
!                                                                                       !
!                   Code development:                                                   !
!                            Lead: Joseph Zhang (OHSU)
!                            Air-sea exchnage: Mike Zulauf (OHSU)
!                            Ecology: Marta Rodrigues/Anabela Oliveira (LNEC)
!                            Sediment: Ligia Pinto/Andre Fortunato (LNEC)
!                            Oil Spill: Alberto Azvedo/Anabela Oliveira (LNEC)
!                            Waves: Aron Roland (Zanke & Partner), Tai-Wen Hsu (NCKU),
!                                   Ulrich Zanke (Zanke & Partner) 
!                            Water quality: Harry Wang... (VIMS)
!
!                   Scientific direction: Antonio Baptista                              !
!                                                                                       !
!               Copyright 2003-2004 Oregon Health and Science University                !
!                              All Rights Reserved                                      !
!                                                                                       !
!       The heat exchange module makes use of the bulk aerodynamic surface flux         !
!       algorithm introduced by Zeng et al (1998), and the polynomial fits to           !
!       saturation vapor pressure of Flatau et al (1992):                               ! 
!       Zeng, X., M. Zhao, and R. E. Dickinson, 1998:  Intercomparison of bulk          !
!       aerodynamic algorithms for the computation of sea surface fluxes using          !
!       TOGA COARE and TAO data.  J. Clim., 11, 2628-2644.                              !
!       Flatau, P. J., R. L. Walko and W. R. Cotton, 1992:  Polynomial fits to          !
!       saturation vapor pressure.  J. Appl. Meteor., 31, 1507-1513.                    !
!                                                                                       !
!       Attenuation of solar radiation (and solar heating) within the water column      !
!       is based upon the expression given by Paulson and Simpson (1977), for the       !
!       water types defined by Jerlov (1968):                                           !
!       Jerlov, N. G., Optical Oceanography, Elsevier, 1968.                            !
!       Paulson, C. A., and J. J. Simpson, Irradiance measurements in the upper         !
!       ocean, J. Phys. Oceanogr., 7, 952-956, 1977.                                    !
!                                                                                       !
!       In addition, the module must be linked with netcdf library.
!
!       The GOTM option was taken from gotm.net.
!
!       A very special thanks to Dr. Tim Campbell, the author of MPI ELCIRC. We have
!       largely followed his style in MPI SELFE. We are indebted to his generous
!       help thru'out the process of parallelizing SELFE.
!
!                                                                                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                       !
!===============================================================================
!===============================================================================
! SELFE MAIN PROGRAM
!===============================================================================
!===============================================================================

#ifdef USE_SWAN
     subroutine elfe(comm1,istart_elfe,iend_elfe,np_gb,eta_gb,dav_gb) 
#else
      program elfe
#endif

#ifdef USE_MPIMODULE
      use mpi
#endif
      use elfe_glbl
      use elfe_msgp

#ifdef USE_GOTM
      use turbulence, only: init_turbulence, do_turbulence, cde, tke1d => tke, eps1d => eps, L1d => L, num1d => num, nuh1d => nuh
      use mtridiagonal, only: init_tridiagonal

#endif

#ifdef USE_ECO
      USE bio_param
      USE biology
      USE eclight
#endif

#ifdef USE_ICM
      USE icm_param, only: iSun,iWQPS,nps,DTD,WWPRPOC,WWPLPOC, &
                          &WWPDOCA,WWPRPON,WWPLPON,WWPDON,WWPNH4,WWPNO3, &
                          &WWPRPOP,WWPLPOP,WWPDOP,WWPPO4t,WWPSU,WWPSAt, &
                          &WWPCOD,WWPDO,xPSQ,xPSK,PRPOC,PLPOC,PDOCA,PRPON, &
                          &PLPON,PDON,PNH4,PNO3,PRPOP,PLPOP,PDOP,PPO4t,PSU, &
                          &PSAt,PCOD,PDO     !added by YC
#endif

#ifdef USE_NAPZD
      USE biology_napzd
#endif

#ifdef USE_SED
!      USE sed_param
      USE sed_mod, only : Wsed,Srho
      USE ocean_mod, only : bedldu,bedldv,Zob
#if defined BEDLOAD_VR || defined BEDLOAD_MPM && defined SED_MORPH
        USE ocean_mod, only : mcoefd
#endif
#endif USE_SED

#ifdef USE_OIL
#endif

#ifdef USE_HA
      USE harm
#endif

      implicit real(rkind)(a-h,o-z),integer(i-n)
#ifndef USE_MPIMODULE
      include 'mpif.h'
#endif

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Data specification section of MAIN
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

#ifdef USE_SWAN
      integer, intent(in) :: comm1,istart_elfe,iend_elfe,np_gb
      real (8), intent(out) :: eta_gb(np_gb),dav_gb(2,np_gb)
      save
#endif

! Message passing arrays
      integer,allocatable :: srqst(:),sstat(:,:)
      integer,allocatable :: rrqst(:),rstat(:,:)
      integer,allocatable :: imapping(:,:),iegath(:)

! Output handles
      character(len=8) :: date
      character(len=10) :: timestamp
      character(len=72) :: it_char
      character(72) :: fgb  ! Processor specific global output file name
      integer :: lfgb       ! Length of processor specific global output file name
      real(4) :: floatout,floatout2
      real(8) :: double1 !for hotstart.in

! Input handles added by YC
      character(len=72) :: windfile,apssfile         !added by YC
! arrays for waterquality model added by YC
      real(8),save,allocatable :: WSRP(:),WSLP(:),WSPB1(:),WSPB2(:), &
     &WSPB3(:),turb(:),WRea(:)    
      real(8),allocatable :: PSQ(:)
      integer,allocatable :: PSK(:)

! Inter-subdomain backtracking
      logical :: lbt,lbtgb
      integer :: nbtrk
      type(bt_type),allocatable :: btlist(:)

! Geometry
      allocatable :: sigmap(:,:),sigma_prod(:,:,:)
      allocatable :: icolor1(:),icolor2(:),ifront(:),ifront2(:)
      allocatable :: akr(:,:)
      allocatable :: akrp(:),ipiv(:),work4(:)

! Boundary forcings
      logical :: lbc   !local dummy
      logical :: lerbc !flag to indicate a partitioned radiation boundary segment
      logical :: lflbc !flag to indicate existence of ifltype/=0
      logical :: ltobc !flag to indicate a patitioned temperature open boundary segment
      logical :: lsobc !flag to indicate a patitioned salinity open boundary segment
      logical,allocatable :: lelbc(:)
      allocatable :: iettype(:),ifltype(:),itetype(:),isatype(:),itrtype(:),trobc(:), &
     &tobc(:),sobc(:),vobc1(:),vobc2(:)
      allocatable :: tamp(:),tnf(:),tfreq(:),jspc(:),tear(:)
      allocatable :: amig(:),ff(:),face(:)
      allocatable :: emo(:,:,:),efa(:,:,:),vmo(:,:),vfa(:,:)
      allocatable :: eth(:,:),qthcon(:),ath(:)
      real(4), allocatable :: a2th(:,:,:)
      allocatable :: uth(:,:),vth(:,:),uthnd(:,:,:),vthnd(:,:,:)
      allocatable :: carea(:),z_r2(:),eta_mean(:)
!      allocatable :: iet1lg(:),ifl1lg(:),ite1lg(:),isa1lg(:)

! Flow arrays
      character(len=2) :: tvd_mid,flimiter,tvd_mid2,flimiter2,stringvalue
      logical :: up_tvd
      allocatable :: ptbt(:,:,:),sdbt(:,:,:),webt(:,:),bubt(:,:)
      allocatable :: out3(:,:),out2(:)
      allocatable :: windx1(:),windy1(:),windx2(:),windy2(:),surf_t1(:),surf_t2(:) !YC
      allocatable :: tau(:,:),iadv(:),nsubd(:),windfactor(:),surf_t(:) !YC
      allocatable :: pr1(:),airt1(:),shum1(:),pr2(:),airt2(:),shum2(:),pr(:)
      allocatable :: sflux(:),srad(:),srad_e(:),tauxz(:),tauyz(:),fluxsu(:),fluxlu(:),hradu(:),hradd(:)
      allocatable :: chi(:),cori(:),Cd(:),rough_p(:)
      allocatable :: dfz(:),dzz(:),erho(:,:)
      allocatable :: hvis(:,:),d2uv(:,:,:) 
      allocatable :: icoef(:),jcoef(:),e2coef(:),qel(:),sparsem(:,:),elbc(:) !,imap(:),qel2(:),eta3(:)
      allocatable :: hhat(:),bigu(:,:),ghat1(:,:),sne(:,:),area_e(:)
      allocatable :: bcc(:,:,:),hp_int(:,:,:),ctmp(:) !hp_int indices reversed
      allocatable :: ibt_p(:),ibt_s(:),t_nudge(:),s_nudge(:),tr_nudge(:),dr_dxy(:,:,:),fun_lat(:,:),etp(:)
      allocatable :: elev_nudge(:),uv_nudge(:)
      allocatable :: fluxprc(:),fluxevp(:),dav(:,:),elevmax(:),dav_max(:,:),dav_maxmag(:) !max. elev. & dahv at nodes for all steps for tsunami
      allocatable :: tr_nd(:,:,:)
      allocatable :: deta2_dx(:),deta2_dy(:),deta1_dx(:),deta1_dy(:),dpr_dx(:),dpr_dy(:),detp_dx(:),detp_dy(:)
      allocatable :: kbs_e(:),iwater_type(:)
      real(rkind), allocatable :: h1d(:),SS1d(:),NN1d(:)
      real(4), dimension(:,:), allocatable :: tnd_nu1,snd_nu1,tnd_nu2,snd_nu2,tnd_nu,snd_nu
!     Turbulence closure arrays
      allocatable :: diffmax(:),diffmin(:),dfq1(:,:),dfq2(:,:),q2tmp(:),xltmp(:)
      allocatable :: rzbt(:),shearbt(:),xlmax(:),cpsi3(:),cpsi2p(:),q2ha(:),xlha(:),xlsc0(:)

!     Wild-card arrays
      allocatable :: nwild(:),nwild2(:),swild(:),swild2(:,:) !swild2 dimension must match that in vinter()
      allocatable :: swild3(:),rwild(:,:)
      real(8), allocatable :: swild4(:,:),swild10(:,:) !double precision for hotstart.in
      allocatable :: swild99(:,:),swild98(:,:,:) !used for exchange (deallocate immediately afterwards)
      allocatable :: buf1(:,:),buf2(:,:),buf3(:)
      real(4), allocatable :: swild8(:,:) !used in ST nudging
      real(4) :: swild_sng(8)
      type(llist_type),pointer :: llp
      logical :: ltmp1,ltmp2

!     Non-hydrostatic arrays
      allocatable :: qnon(:,:),qhat(:,:),dqnon_dxy(:,:,:),qmatr(:,:,:,:),qir(:,:),ihydro(:)

!     Solver arrays for TRIDAG
      allocatable :: alow(:),bdia(:),cupp(:),rrhs(:,:),soln(:,:),gam(:)

!     Some constants
      logical :: first_call=.true.
!      real(rkind),parameter :: omega_e=7.292e-5 !angular freq. of earth rotation 
!      real(rkind),parameter :: grav=9.81
      real(rkind),parameter :: shw=4184.  !specific heat of pure water

!     Tracers
      character(len=48) :: inputfile
      integer :: flag_model,flag_ic

#ifdef USE_ECO 
!...  MFR - other variables to atmospheric parameters (when nws=0 ... probably to clean later...)
      real(rkind), allocatable :: Pair(:), Tair(:), Hair(:), Uwind(:), Vwind(:), cloud(:)

!...  MFR - Tracer models
      real(rkind) :: tr_tmp1
!      allocatable :: SpecIr(:,:), avcos(:,:)
#endif

#ifdef USE_NAPZD
      allocatable :: Bio_bdefp(:,:)
#endif

#ifdef USE_SED
      allocatable :: tr_tc(:,:),tr_tl(:,:)     
#if defined BEDLOAD_VR
        real(8),allocatable :: dave(:),dahv(:,:)
#endif

!#   ifdef SED_MORPH
!      allocatable :: dhnd(:),bed_thick(:,:)
!#   endif
#endif USE_SED

#ifdef USE_OIL
#endif

!     Wind wave model arrays
!# ifdef  USE_WWM
!      allocatable :: out_wwm(:,:),wwave_force(:,:,:)
!# endif

#ifdef USE_HA
      INTEGER NTSTEPS,ITMV
      REAL(8) TIMEBEG
      REAL(rkind) FMV
      allocatable :: XVELAV(:),YVELAV(:),XVELVA(:),YVELVA(:),ELAV(:),ELVA(:)
#endif

!     Station and other output arrays
      allocatable :: xsta(:),ysta(:),zstal(:),zsta(:),iep_sta(:),iep_flag(:),arco_sta(:,:), &
     &iof_sta(:),sta_out(:,:),sta_out_gb(:,:),indx_out(:,:),indx_wwm_out(:)

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Executable section of SELFE
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!     Initialize parallel environment
!      call parallel_init
!     Initialize MPI
#ifdef USE_SWAN
      comm=comm1
#else
      call mpi_init(ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort(error=ierr)
  
! Duplicate communication space to isolate ELFE communication
      call mpi_comm_dup(MPI_COMM_WORLD,comm,ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort(error=ierr)
#endif

! Get number of processors
      call mpi_comm_size(comm,nproc,ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort(error=ierr)

! Get rank
      call mpi_comm_rank(comm,myrank,ierr)
      if(ierr/=MPI_SUCCESS) call parallel_abort(error=ierr)

!==============================================================================
!     For loose coupling, bypass initialization if not forst call
!==============================================================================
      if(first_call) then
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

      if(myrank==0) then
        open(91,file='total_ST.dat',status='replace')
!        open(92,file='total_TR.dat',status='replace')
      endif

!     Debug
!      open(10,file='rank.in',status='old')
!      read(10,*) irank0 
!      close(10)

!     Rank 0 open mirror file
      if(myrank==0) open(16,file='mirror.out',status='replace')

!     Echo date and time
      call date_and_time(date,timestamp)
      if(myrank==0) write(16,'(4a)')'Run begins at ',date,', ',timestamp

!     Setup cyclic node index
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

!     Read from param.in the 2D flag (0: 3D model; other: 2D model)
      call get_param('im2d',1,im2d,realvalue,stringvalue)
      if(im2d==0) then !3D
        lm2d=.false.
      else !2D
        lm2d=.true.
        !Implicit Coriolis
        call get_param('theta2',2,itmp,theta2,stringvalue)
      endif

!     Check modules for 2D model
      if(lm2d) then
#if defined USE_ECO || defined USE_ICM || defined USE_SED || defined PREC_EVAP || defined USE_GOTM || defined USE_NAPZD
        write(errmsg,*)'2D model cannot have certain modules!'
        call parallel_abort(errmsg)
#endif
      endif !lm2d

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
      allocate(rrqst(nproc-1),stat=istat)
      if(istat/=0) call parallel_abort('main: rrqst allocation failure')
      allocate(rstat(MPI_STATUS_SIZE,nproc-1),stat=istat)
      if(istat/=0) call parallel_abort('main: rstat allocation failure')
      allocate(srqst(nproc-1),stat=istat)
      if(istat/=0) call parallel_abort('main: srqst allocation failure')
      allocate(sstat(MPI_STATUS_SIZE,nproc-1),stat=istat)
      if(istat/=0) call parallel_abort('main: sstat allocation failure')

!     Allocate the remaining grid geometry arrays held in elfe_glbl
      allocate(kbe(nea),idry_e(nea),idry_e0(nea),lqk(nea),ie_kr(nea),krvel(nea),ze(nvrt,nea),dl(nea,3,2), &
           dp00(npa),kfp(npa),kbp(npa),kbp00(npa),kbp_e(np),idry(npa),iback(npa),hmod(npa),znl(nvrt,npa), &
           kbs(nsa),idry_s(nsa),isidenei2(ns,4),zs(nvrt,nsa),side_ac(nsa,2,2),side_x(nsa,2),delj(ns), &
           ibnd_ext_int(npa),edge_angle(2,npa),pframe(3,3,npa),stat=istat)
      if(istat/=0) call parallel_abort('MAIN: grid geometry arrays allocation failure')
!'

!     Allocate the remaining arrays held in elfe_glbl, except for Kriging related arrays
!     ntracers,ntracers2 and mntr already read in from grid_subs.F90
      allocate(tsel(2,nvrt,nea),tem0(nvrt,npa),sal0(nvrt,npa),eta1(npa),eta2(npa), &
           we(nvrt,nea),we_fv(nvrt,nea),su2(nvrt,nsa),sv2(nvrt,nsa),ufg(nvrt,nea,3),vfg(nvrt,nea,3), &
           tsd(nvrt,nsa),ssd(nvrt,nsa),tnd(nvrt,npa),snd(nvrt,npa), &
           prho(nvrt,npa),q2(npa,nvrt),xl(npa,nvrt),xlmin2(npa), &
           uu2(nvrt,npa),vv2(nvrt,npa),ww2(nvrt,npa),bdef(npa),bdef1(npa),bdef2(npa),dfh(npa,nvrt), &
           tr_el(mntr,nvrt,nea),bdy_frc(mntr,nvrt,nea),flx_sf(mntr,nea),flx_bt(mntr,nea), &
           xlon_el(nea),ylat_el(nea),albedo(npa),cspline_ypp_nd(2,nvrt,npa), &
           cspline_ypp_sd(2,nvrt,nsa),stat=istat)
      if(istat/=0) call parallel_abort('MAIN: dynamical arrays allocation failure')
!'

!     Allocate boundary forcings (remaining are allocated after dimensioning parameters are read)
      allocate(lelbc(npa),iettype(max(1,nope_global)),ifltype(max(1,nope_global)), &
           itetype(max(1,nope_global)),isatype(max(1,nope_global)), &
           itrtype(max(1,nope_global)),trobc(nope_global),tobc(nope_global), &
           sobc(nope_global),vobc1(nope_global),vobc2(nope_global), &
           eth(nope_global,mnond_global),tth(nope_global,mnond_global,nvrt),sth(nope_global,mnond_global,nvrt), &
           qthcon(nope_global),ath(nope_global),a2th(max(2,ntracers),nvrt,neta_global),carea(nope_global), &
           uth(nsa,nvrt),vth(nsa,nvrt), &
           uthnd(nope_global,mnond_global,nvrt),vthnd(nope_global,mnond_global,nvrt), &
           eta_mean(npa),trth(max(1,ntracers),nvrt,mnond_global,max(1,nope_global)),stat=istat)
!           iet1lg(nope),ifl1lg(nope),ite1lg(nope),isa1lg(nope),stat=istat)
      if(istat/=0) call parallel_abort('MAIN: 1st bnd forcings allocation failure')            
!'

!     All other arrays
      allocate(sigmap(nvrt,10),sigma_prod(nvrt,nvrt,-4:4),icolor1(npa),icolor2(npa), &
           ptbt(4,nvrt,npa),sdbt(4,nvrt,nsa),webt(nvrt,nea),bubt(nea,2),out3(nvrt,3),out2(12), &
           windx1(npa),windy1(npa),windx2(npa),windy2(npa),windx(npa),windy(npa), &
           tau(npa,2),iadv(npa),nsubd(npa),windfactor(npa),pr1(npa),airt1(npa),shum1(npa), &
           pr2(npa),airt2(npa),shum2(npa),pr(npa),sflux(npa),srad(npa),srad_e(nea),tauxz(npa),tauyz(npa), &
           fluxsu(npa),fluxlu(npa),hradu(npa),hradd(npa),chi(nsa),cori(nsa),Cd(nsa), &
           Cdp(npa),rough_p(npa),dfv(npa,nvrt),dfz(2:nvrt),dzz(2:nvrt), &
           hdif(nvrt,npa),hvis(nvrt,nea),d2uv(2,nvrt,nsa), & !horcon(nsa), &
           icoef(npa+1),jcoef(npa*(mnei+1)),e2coef(npa*(mnei+1)),qel(npa),sparsem(np,0:(mnei+1)), & !sparsem for non-ghosts only
           elbc(npa),hhat(nsa),bigu(nsa,2),ghat1(nea,2), &
           sne(nvrt,3),area_e(nvrt),bcc(2,nvrt,nsa),hp_int(nvrt,nea,2),ctmp(0:nvrt),&
           ibt_p(npa),ibt_s(nsa),t_nudge(npa),s_nudge(npa),tr_nudge(npa),dr_dxy(2,nvrt,nsa), &
           elev_nudge(npa),uv_nudge(npa), &
           fun_lat(npa,0:2),etp(npa),dav(npa,2),elevmax(npa),dav_max(npa,2),dav_maxmag(npa), &
           fluxprc(npa),fluxevp(npa),h1d(0:nvrt),SS1d(0:nvrt),NN1d(0:nvrt), &
           tnd_nu1(npa,nvrt),snd_nu1(npa,nvrt),tnd_nu2(npa,nvrt),snd_nu2(npa,nvrt),tnd_nu(npa,nvrt),snd_nu(npa,nvrt), &
           diffmax(npa),diffmin(npa),dfq1(npa,nvrt),dfq2(npa,nvrt),q2tmp(nvrt),xltmp(nvrt), &
           rzbt(nvrt),shearbt(2:nvrt),xlmax(nvrt),cpsi3(2:nvrt),cpsi2p(2:nvrt),q2ha(2:nvrt),xlha(2:nvrt),xlsc0(npa), &
!          Note: swild, swild2, swild10 will be re-dimensioned (larger dimension) later
           nwild(nea+12),nwild2(ne_global),swild(nsa+nvrt+12+ntracers),swild2(nvrt,12),swild10(max(3,nvrt),12), &
           swild3(20+ntracers),swild4(nvrt,3+2*ntracers),swild8(nvrt,2),&
           alow(max(3,nvrt)),bdia(max(3,nvrt)),cupp(max(3,nvrt)),rrhs(nvrt,100),soln(nvrt,100),gam(nvrt), &
           deta1_dx(nsa),deta1_dy(nsa),deta2_dx(nsa),deta2_dy(nsa),dpr_dx(nsa),dpr_dy(nsa),detp_dx(nsa),detp_dy(nsa), &
           kbs_e(nsa),iwater_type(npa),rho_mean(nvrt,nea),erho(nvrt,nea),& 
           PSQ(nea),PSK(nea),surf_t1(npa),surf_t2(npa),surf_t(npa),rwild(np_global,3),stat=istat)
      if(istat/=0) call parallel_abort('MAIN: other allocation failure')

!     Tracers
      if(ntracers==0) then !allocate some trivial arrays (for reference only)
        allocate(trel0(1,1,1),trel(1,1,1),tr_nd(1,1,1),stat=istat)
      else
        allocate(trel0(ntracers,nvrt,nea),trel(ntracers,nvrt,nea), &
                 tr_nd(ntracers,nvrt,np),stat=istat)
      endif
      if(istat/=0) call parallel_abort('MAIN: other allocation failure')

#ifdef USE_ECO
      if(ntracers/=25) call parallel_abort('MAIN: ntracer/=25 in EcoSim')
      allocate(Pair(nea), Tair(nea), Hair(nea), Uwind(nea), Vwind(nea), cloud(nea), &
               SpecIr(nea,Nbands),avcos(nea,Nbands),stat=istat)
      if(istat/=0) call parallel_abort('MAIN: EcoSim allocation failure')
#endif

#ifdef USE_NAPZD
      allocate(Bio_bdef(nvrt,nea), Bio_bdefp(nvrt,np), stat=istat)
      if(istat/=0) call parallel_abort('MAIN: NAPZD allocation failure')
#endif

#ifdef USE_SED
      allocate(tr_tc(ntracers,nea),tr_tl(ntracers,nea),Zob(nea),stat=istat)
      if (istat/=0) call parallel_abort('Main: sed. allocation failure1')
#ifdef BEDLOAD_VR
        allocate(dave(nea),dahv(npa,2),stat=istat)
        if (istat/=0) call parallel_abort('Main:sed. allocation failure')
#endif
#if defined BEDLOAD_MPM || defined BEDLOAD_VR && defined SED_MORPH
        allocate(mcoefd(np,0:(mnei+1)),stat=istat)
        if(istat/=0) call parallel_abort('MAIN: Sed. allocation failure')
#endif
#endif USE_SED

!     Wave model arrays
#ifdef  USE_WWM
      allocate(wwave_force(nvrt,nsa,2),out_wwm(npa,30),stat=istat)
      if(istat/=0) call parallel_abort('MAIN: WWM allocation failure')
#endif

!     Non-hydrostatic arrays; flag nonhydro already read in from aquire_hgrid
      allocate(qnon(nvrt,npa),stat=istat)
      if(istat/=0) call parallel_abort('MAIN: Nonhydro allocation failure (1)')
!'
      qnon=0 !initialize

      if(nonhydro==1) then
        allocate(qhat(nvrt,npa),dqnon_dxy(2,nvrt,nsa),qmatr(nvrt,-1:1,0:(mnei+1),np), &
     &qir(nvrt,np),ihydro(npa),stat=istat) 
        if(istat/=0) call parallel_abort('MAIN: Nonhydro allocation failure')
!'
      endif

!     Test message passing here
!      call parallel_finalize
!      stop

!     Initialize some arrays and constants
      tempmin=0; tempmax=40; saltmin=0; saltmax=42
!     temp fix
!      tempmin=0; tempmax=1005; saltmin=0; saltmax=42
      pr1=0; pr2=0; pr=0 !uniform pressure (the const. is unimportant)
      uthnd=-99; vthnd=-99; eta_mean=-99; uth=-99; vth=-99; !flags
      fluxsu00=0; srad00=0 !for nws/=3
      elevmax=-1.e34; dav_maxmag=-1
      tsel=0; trel=0

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
!          n1=isidenode(i,1)
!          n2=isidenode(i,2)
!          call compute_ll(xcj(i),ycj(i),zcj(i),rlon,rlat)
!          rad=sqrt(xcj(i)**2+ycj(i)**2+zcj(i)**2)
!          x2=rad*cos(rlat)*cos(rlon)
!          y2=rad*cos(rlat)*sin(rlon)
!          z2=rad*sin(rlat)
!          rlon2=(xlon(n1)+xlon(n2))/2 !has problem around 180 deg. etc
!          rlat2=(ylat(n1)+ylat(n2))/2
!          dis=sqrt((x2-xcj(i))**2+(y2-ycj(i))**2+(z2-zcj(i))**2)
!          write(12,*)'Side ll:',iplg(isidenode(i,1:2)),rlon/pi*180,rlat/pi*180,rlon2/pi*180,rlat2/pi*180, &
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
!          rlon2=sum(xlon(nm(i,1:3)))/3 !has problem around 180 deg.
!          rlat2=sum(ylat(nm(i,1:3)))/3 !has problem around 180 deg.
!          dis=sqrt((x2-xctr(i))**2+(y2-yctr(i))**2+(z2-zctr(i))**2)
!          write(12,*)'Elem. ll:',iplg(nm(i,:)),rlon/pi*180,rlat/pi*180,rlon2/pi*180,rlat2/pi*180, &
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
        if(myrank==0) write(16,*)'Max. pframe dev. from radial= ',real(tmp) !zdev_max
      endif !ics==2

!...  Modified depth
      dpmax=-1.e25 !max. depth
      do i=1,npa
        hmod(i)=min(dp(i),h_s)
        if(dp(i)>dpmax) dpmax=dp(i)
      enddo !i=1,npa
!     Save intial depth for bed deformation case
      dp00=dp

      if(ztot(1)>=-dpmax) then
        write(errmsg,*)'1st z-level must be below max. depth:',dpmax
        call parallel_abort(errmsg)
      endif

!...  Derivatives of shape functions
!...  For ics=2, this is done inside element frame
      do i=1,nea
        do j=1,3
          dl(i,j,1)=(yel(nx(j,1),i)-yel(nx(j,2),i))/2/area(i) !dL_1/dx
          dl(i,j,2)=(xel(nx(j,2),i)-xel(nx(j,1),i))/2/area(i) !dL_1/dy
        enddo !j

!        n1=nm(i,1)
!        n2=nm(i,2)
!        n3=nm(i,3)
!        dl(i,1,1)=(ynd(n2)-ynd(n3))/2/area(i) !dL_1/dx
!        dl(i,2,1)=(ynd(n3)-ynd(n1))/2/area(i) !dL_2/dx
!        dl(i,3,1)=(ynd(n1)-ynd(n2))/2/area(i)
!        dl(i,1,2)=(xnd(n3)-xnd(n2))/2/area(i) !dL_1/dy
!        dl(i,2,2)=(xnd(n1)-xnd(n3))/2/area(i)
!        dl(i,3,2)=(xnd(n2)-xnd(n1))/2/area(i)
      enddo !i=1,nea

!Debug
!      if(ics==2) then
!        do i=1,nea
!          do j=1,3
!            swild(j)=xel(j,i)+2*yel(j,i)
!          enddo !j
!          write(12,*)'db_dx-ana=',dot_product(swild,dl(i,1:3,1))-1
!          write(12,*)'db_dy-ana=',dot_product(swild,dl(i,1:3,2))-2
!        enddo !i
!      endif
!      call parallel_finalize
!      stop      

!...  Compute delj for internal resident sides only (used only in horizontal diffusion)
      do i=1,ns !resident only
        if(is(i,2)==0) cycle
        delj(i)=sqrt((xctr(is(i,2))-xctr(is(i,1)))**2+(yctr(is(i,2))-yctr(is(i,1)))**2+ &
     &(zctr(is(i,2))-zctr(is(i,1)))**2)    
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
      do i=1,nea
        !Error: won't work near dateline!!!! Try to use compute_ll
        xlon_el(i)=(xlon(nm(i,1))+xlon(nm(i,2))+xlon(nm(i,3)))/3*180/pi !in degrees
        ylat_el(i)=(ylat(nm(i,1))+ylat(nm(i,2))+ylat(nm(i,3)))/3*180/pi
      enddo !i
#endif

!...  Classify interior/exterior bnd node and calculate edge angles (for WWM only)
!...  WARNING: if WWM is used, the _land_ b.c. part of hgrid.gr3 must have flags for (exterior) land (0) and
!...           island (1) bnds, and no open bnd is allowed on islands
!...  ibnd_ext_int:
!       0: interior node; 
!       1: exterior bnd (including open and land bnds); 
!      -1: interior bnd (land bnd only)
#ifdef USE_WWM
!Error: in lat/lon, angles do not make sense
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

!       Compute 2 orientation angles of 2 bnd sides of a bnd node
        edge_angle=-99 !flag for non-bnd nodes
        do i=1,np
          if(isbnd(1,i)==0) cycle
!         Bnd node; find 2 neighbors
          do ll=1,2 !side
            if(ll==1) then
              ie=ine(i,1)
              id=iself(i,1)
              nd2=nm(ie,nx(id,1))
            else
              ie=ine(i,nne(i))
              id=iself(i,nne(i))
              nd3=nm(ie,nx(id,2))
            endif
          enddo !ll
          if(isbnd(1,nd2)==0.or.isbnd(1,nd3)==0) then
            write(errmsg,*)'Wrong bnd nodes:',iplg(i),iplg(nd2),iplg(nd3)
            call parallel_abort(errmsg)
          endif

          if(ics==1) then
            tmp1=atan2(ynd(nd2)-ynd(i),xnd(nd2)-xnd(i))
            tmp2=atan2(ynd(nd3)-ynd(i),xnd(nd3)-xnd(i))
          else !lat/lon
            call project_pt('g2l',xnd(nd2),ynd(nd2),znd(nd2),(/xnd(i),ynd(i),znd(i)/), &
     &pframe(:,:,i),xtmp2,ytmp2,ztmp2)
            call project_pt('g2l',xnd(nd3),ynd(nd3),znd(nd3),(/xnd(i),ynd(i),znd(i)/), &
     &pframe(:,:,i),xtmp3,ytmp3,ztmp3)
            tmp1=atan2(ytmp2,xtmp2)
            tmp2=atan2(ytmp3,xtmp3)
          endif !ics

          if(tmp1<0) tmp1=tmp1+2*pi
          if(tmp2<0) tmp2=tmp2+2*pi
          edge_angle(1,i)=tmp1; edge_angle(2,i)=tmp2
        enddo !i=1,np

        swild(1:npa)=edge_angle(1,1:npa)
        call exchange_p2d(swild)
        edge_angle(1,1:npa)=swild(1:npa)
        swild(1:npa)=edge_angle(2,1:npa)
        call exchange_p2d(swild)
        edge_angle(2,1:npa)=swild(1:npa)
        !Debug
        !do i=1,npa
        !  if(edge_angle(i)>=0) write(12,*)'After:',iplg(i),edge_angle(i)
        !enddo !i
#endif USE_WWM

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
            fun_lat(ii,0)=3*sin(ylat(ii))**2-1
            fun_lat(ii,1)=sin(2*ylat(ii))
            fun_lat(ii,2)=cos(ylat(ii))**2
          endif
        enddo !i
        close(32)
      
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
        write(errmsg,*)'Inconsistent # of open bnds',nope1,nope
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
      nnode_et=0 !total # of open bnd nodes that require elev3D.th
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
!	  Mock reading
          open(50,file='elev.th',status='old')
!          do j=1,ntime
!            read(50,*) ttt,et
!          enddo !j
!          rewind(50)
        else if(iettype(k)==2) then
          read(31,*) eth(k,1)
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
!          open(54,file='elev3D.th',status='old')
        else if(iettype(k)/=0) then
          call parallel_abort('Invalid iettype')
        endif

!       For ics=2, uthnd, vthnd, uth, vth are all in lat/lon frame 
!       (even at poles), with exception for uthnd, vthnd and Flather b.c. (see below)
!       For ics=1, they are in global frame
        if(ifltype(k)==1) then
          nfltype=nfltype+1
          open(51,file='flux.th',status='old')
!          do j=1,ntime
!            read(51,*) ttt,qq
!          enddo
!          rewind(51)
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
!          open(58,file='uv3D.th',status='old')
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
            read(31,*) uthnd(k,j,1:nvrt) !used to denote normal vel. (i.e. along xs axis)
          enddo !j
!       ifltype(k)=0: zero out vertical velocity for transport in the open bnd elements
        else if(ifltype(k)/=0) then
          write(errmsg,*) 'Invalid ifltype:',ifltype(k)
          call parallel_abort(errmsg)
        endif

        tobc(k)=0 !init. for checking below
        if(itetype(k)==1) then
          ntetype=ntetype+1
          read(31,*) tobc(k) !nudging factor for inflow (no b.c. for outflow)
          open(52,file='temp.th',status='old')
!          do j=1,ntime
!            read(52,*) ttt,temp
!          enddo
!          rewind(52)
        else if(itetype(k)==2) then
          read(31,*) tth(k,1,1)
          read(31,*) tobc(k) !nudging factor
        else if(itetype(k)==3) then
          read(31,*) tobc(k) !nudging factor
        else if(itetype(k)==4) then
          ntetype2=ntetype2+1
          nnode_te=nnode_te+nond_global(k)
!          open(56,file='temp3D.th',status='old')
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
          open(53,file='salt.th',status='old')
!          do j=1,ntime
!            read(53,*) ttt,sal
!          enddo
!          rewind(53) 
        else if(isatype(k)==2) then
          read(31,*) sth(k,1,1)
          read(31,*) sobc(k) !nudging factor
        else if(isatype(k)==3) then
          read(31,*) sobc(k) !nudging factor
        else if(isatype(k)==4) then
          nsatype2=nsatype2+1
          nnode_sa=nnode_sa+nond_global(k)
!          open(57,file='salt3D.th',status='old')
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
        read(31,*) itr_met !=1: upwind; 2: TVD
        if(itr_met/=1.and.itr_met/=2) then
          write(errmsg,*)'Unknown tracer method',itr_met
          call parallel_abort(errmsg)
        endif
        if(itr_met==2) read(31,*) tvd_mid2,flimiter2

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
              inputfile='htr_'//ifile_char(1:ifile_len)//'.th'
!              inputfile='htr_'//trim(ifile_char)//'.th'
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

!       Nudging
        read(31,*) inu_tr
        if(inu_tr/=0.and.inu_tr/=1) then
          write(errmsg,*)'Wrong inu_tr:',inu_tr
          call parallel_abort(errmsg)
        endif
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
        endif !inu_tr/=0
      endif !ntracers

      close(31)

!     Check dimension of a2th
      if(max(nnode_et,nnode_fl,nnode_te,nnode_sa,nnode_tr)>neta_global) then
        write(errmsg,*) 'MAIN: impossible! Dimension overflow for a2th:',nnode_et,nnode_fl,nnode_te,nnode_sa,nnode_tr
        call parallel_abort(errmsg)
      endif
!     Binary record length for *3D.th at each time step
      nrecl_et=nbyte*(1+nnode_et) !single precision
      nrecl_fl=nbyte*(1+nnode_fl*2*nvrt)
      nrecl_te=nbyte*(1+nnode_te*nvrt)
      nrecl_sa=nbyte*(1+nnode_sa*nvrt)
      nrecl_tr=nbyte*(1+nnode_tr*nvrt*ntracers)
      if(nettype2/=0) open(54,file='elev3D.th',access='direct',recl=nrecl_et,status='old')
      if(nfltype2/=0) open(58,file='uv3D.th',access='direct',recl=nrecl_fl,status='old')
      if(ntetype2/=0) open(56,file='temp3D.th',access='direct',recl=nrecl_te,status='old')
      if(nsatype2/=0) open(57,file='salt3D.th',access='direct',recl=nrecl_sa,status='old')
      if(ntrtype2/=0) open(59,file='tr3D.th',access='direct',recl=nrecl_tr,status='old')

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Read input parameters from unit 15
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!      open(15,file='param.in',status='old')
!      read(15,'(a48)') version
!      read(15,'(a48)') start_time
!      read(15,*) ipre !pre-processing flag (to output obs.out)
!      read(15,*) ntracers
!      call get_param('version',0,itmp,tmp,version)
!      call get_param('start_time',0,itmp,tmp,start_time)
      version='description' !not really used
      call get_param('ipre',1,ipre,tmp,stringvalue) !pre-processing flag (to output obs.out)
!      call get_param('ntracers',1,ntracers,tmp,stringvalue)

!...  Option for Williamson test #5 (zonal flow over an isolated mount)
      call get_param('izonal5',1,izonal5,tmp,stringvalue)
      if(izonal5/=0.and.ics==1) call parallel_abort('ics=1 and izonal5/=0')
!'

!...  Option for specifying hydrostatic region for non-hydrostatic model
      if(nonhydro==1) then
        if(ics==2) call parallel_abort('ics=2 and nonhydro==1')
        call get_param('ihydro_region',1,ihydro_region,tmp,stringvalue)
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

!...  Option for calculating nodal vel.
!     -1: vel. interpolated (in btrack) from P_1_NC shape function (using side vel.); side vel. filtered with Shapiro filter; 
!     0: vel. interpolated (in btrack) from P_1 shape function; side vel. filtered with Shapiro filter; 
!...  1: averaging w/o Shapiro filter: more diffusion)
      call get_param('indvel',1,indvel,tmp,stringvalue)
      if(indvel<-1.or.indvel>1) then
        write(errmsg,*)'Illegal indvel:',indvel
        call parallel_abort(errmsg)
      endif

!...  Compute neighborhood for internal sides for Shapiro filter
!...  isidenei2(ns,4): 4 neighboring sides of a _resident_ side
!...  Info for resident sides only!
      if(indvel<=0) then
        fdb='rogue_0000'
        lfdb=len_trim(fdb)
        write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
        open(10,file='outputs/'//fdb,status='replace')

        iabort=0 !abort flag
        loop18: do i=1,ns !resident sides only
          if(is(i,2)==0) cycle loop18

!         Internal sides
          do j=1,2
            ie=is(i,j)
            l0=lindex_s(i,ie)
            if(l0==0) then
              write(errmsg,*)'Cannot find a side'
              call parallel_abort(errmsg)
            endif
            nwild(2*j-1)=js(ie,nx(l0,1))
            nwild(2*j)=js(ie,nx(l0,2))
          enddo !j=1,2
          isidenei2(i,1:4)=nwild(1:4) !local index

!         Check if pt "0" is inside
          if(ics==1) then
            x0=xcj(i); y0=ycj(i)
            x1=xcj(isidenei2(i,1)); y1=ycj(isidenei2(i,1))
            x2=xcj(isidenei2(i,2)); y2=ycj(isidenei2(i,2))
            x3=xcj(isidenei2(i,3)); y3=ycj(isidenei2(i,3))
            x4=xcj(isidenei2(i,4)); y4=ycj(isidenei2(i,4))
          else !ics=2
            x0=0; y0=0
            do j=1,4
              isd=isidenei2(i,j)
!              swild2(1,j)=(xcj(isd)-xcj(i))*sframe(1,1,i)+(ycj(isd)-ycj(i))*sframe(2,1,i)+ &
!     &(zcj(isd)-zcj(i))*sframe(3,1,i)
!              swild2(2,j)=(xcj(isd)-xcj(i))*sframe(1,2,i)+(ycj(isd)-ycj(i))*sframe(2,2,i)+ &
!     &(zcj(isd)-zcj(i))*sframe(3,2,i)
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
              if(is(isidenei2(i,2),2)==0.or.is(isidenei2(i,3),2)==0) then
                iabort=1
                write(10,*)iplg(isidenode(i,1:2)),', bnd side (3)'
                cycle loop18
              endif

              nwild(1)=2; nwild(2)=3
              do k=1,2
                id=isidenei2(i,nwild(k))
                ie2=is(id,1)+is(id,2)-is(i,k) !inside aug. domain
                if(is(id,1)<=0.or.is(id,2)<0.or.ie2/=is(id,1).and.ie2/=is(id,2)) then
                  write(errmsg,*)'Filter sides out of aug. domain (1):',iplg(isidenode(i,1:2)),ie2,is(id,1:2)
                  call parallel_abort(errmsg)
                endif
                l0=lindex_s(id,ie2)
                if(l0==0) then
                  write(errmsg,*)'Cannot find a side (9):',k
                  call parallel_abort(errmsg)
                endif
                isidenei2(i,nwild(k))=js(ie2,nx(l0,3-k))
              enddo !k
            endif !ar3
         
            if(ar4<=0) then
              if(is(isidenei2(i,1),2)==0.or.is(isidenei2(i,4),2)==0) then
                iabort=1
                write(10,*)iplg(isidenode(i,1:2)),', bnd side (4)'
                cycle loop18
              endif

              nwild(1)=1; nwild(2)=4
              do k=1,2
                id=isidenei2(i,nwild(k))
                ie2=is(id,1)+is(id,2)-is(i,k)
                if(is(id,1)<=0.or.is(id,2)<0.or.ie2/=is(id,1).and.ie2/=is(id,2)) then
                  write(errmsg,*)'Filter sides out of aug. domain (2):',iplg(isidenode(i,1:2)),ie2,is(id,1:2)
                  call parallel_abort(errmsg)
                endif
                l0=lindex_s(id,ie2)
                if(l0==0) then
                  write(errmsg,*)'Cannot find a side (8):',k
                  call parallel_abort(errmsg)
                endif
                isidenei2(i,nwild(k))=js(ie2,nx(l0,k))
              enddo !k
            endif !ar4
         
!           Check convexity of quad 1-4
            if(ics==1) then
              x0=xcj(i); y0=ycj(i)
              x1=xcj(isidenei2(i,1)); y1=ycj(isidenei2(i,1))
              x2=xcj(isidenei2(i,2)); y2=ycj(isidenei2(i,2))
              x3=xcj(isidenei2(i,3)); y3=ycj(isidenei2(i,3))
              x4=xcj(isidenei2(i,4)); y4=ycj(isidenei2(i,4))
            else !ics=2
              x0=0; y0=0
              do j=1,4
                isd=isidenei2(i,j)
!                swild2(1,j)=(xcj(isd)-xcj(i))*sframe(1,1,i)+(ycj(isd)-ycj(i))*sframe(2,1,i)+ &
!     &(zcj(isd)-zcj(i))*sframe(3,1,i)
!                swild2(2,j)=(xcj(isd)-xcj(i))*sframe(1,2,i)+(ycj(isd)-ycj(i))*sframe(2,2,i)+ &
!     &(zcj(isd)-zcj(i))*sframe(3,2,i)
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
              write(10,*)iplg(isidenode(i,1:2)),'  Concave quad '
!              write(10,*)((isidenode(isidenei2(i,m),mm),mm=1,2),m=1,4)
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
              write(10,*)iplg(isidenode(i,1:2)),'  pt outside quad '
              write(10,*)ar1,ar2,ar3,ar4
!              write(10,*)((isidenode(isidenei2(i,m),mm),mm=1,2),m=1,4)
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
      endif !indvel<=0

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
      call get_param('imm',1,imm,tmp,stringvalue)
!     For moving bed, the output is from original bottom to nvrt
      if(imm<0.or.imm>2) then
        write(errmsg,*)'Unknown imm',imm
        call parallel_abort(errmsg)
      endif
      if(imm==2.and.ics==2) call parallel_abort('imm=ics=2')

!     Initialize variables used in tsunami model (but bdef[1,2] and ibdef are available for all models)
      bdef=0 !total deformation
      ibdef=1 !# of time steps for deformation (deformation rate=0 when it>ibdef)
      if(imm==1) then !read in deformation at all nodes
        call get_param('ibdef',1,ibdef,tmp,stringvalue)
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

!     hotstart option
      call get_param('ihot',1,ihot,tmp,stringvalue)
      if(ihot<0.or.ihot>2) then
        write(errmsg,*)'Unknown ihot',ihot
        call parallel_abort(errmsg)
      endif
#ifdef USE_SWAN
      if(ihot==2) then
        write(errmsg,*)'ihot cannot be 2 for wave models',ihot
        call parallel_abort(errmsg)
      endif
#endif

!...  Center of projection in degrees (used for f-plane approx.)
      call get_param('cpp_lon',2,itmp,slam0,stringvalue)
      call get_param('cpp_lat',2,itmp,sfea0,stringvalue)
      slam0=slam0*pi/180
      sfea0=sfea0*pi/180

!...  Horizontal viscosity option
!     ihorcon =0 means all hvis=0 and no hvis.gr3 is needed
      call get_param('ihorcon',1,ihorcon,tmp,stringvalue)
      if(ihorcon/=0) then
        if(ics==2) call parallel_abort('ics=2 and ihorcon/=0')
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
          hvis(:,i)=(swild(nm(i,1))+swild(nm(i,2))+swild(nm(i,3)))/3
        enddo !i
      endif !ihorcon/=0
      
!...  Horizontal diffusivity option; only works for upwind/TVD
!     ihdif=0 means all hdif=0 and no hdif.gr3 is needed
      call get_param('ihdif',1,ihdif,tmp,stringvalue)
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
      
!...  Implicitness factor
      call get_param('thetai',2,itmp,thetai,stringvalue)

!...  Baroclinic flags
      call get_param('ibcc',1,ibc,tmp,stringvalue)
      call get_param('itransport',1,ibtp,tmp,stringvalue)
      if(ibc/=0.and.ibc/=1) call parallel_abort('Unknown ibcc')
      if(ibtp/=0.and.ibtp/=1) call parallel_abort('Unknown itransport')

      if(ibc==0) then
        if(myrank==0) write(16,*)'You are using baroclinic model'
        call get_param('nrampbc',1,nrampbc,tmp,stringvalue)
        if(nrampbc/=0) call get_param('drampbc',2,itmp,drampbc,stringvalue)
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
      call get_param('rnday',2,itmp,rnday,stringvalue)

!...  dramp not used if nramp=0
      call get_param('nramp',1,nramp,tmp,stringvalue)
      if(nramp/=0) call get_param('dramp',2,itmp,dramp,stringvalue)

      if(nramp/=0.and.nramp/=1) then
        write(errmsg,*)'Unknown nramp',nramp
        call parallel_abort(errmsg)
      endif

!     Time step in seconds
      call get_param('dt',2,itmp,dt,stringvalue)

!...  compute total number of time steps 
      ntime=rnday*86400.d0/dt+0.5

!...  Rank 0 writes global & local volume, energy etc data
      if(myrank==0) then
        open(9,file='flux.dat',status='replace')
        open(13,file='total.dat',status='replace')
        write(13,*)ntime
        write(13,'(a200)')'Time (hours), volume, mass, potential E, kinetic E, total E, friction loss (Joule), energy leak (Joule)'
!'
      endif

!...  Advection flag for momentum eq.; 1-Euler; 2: R-K
      call get_param('nadv',1,nadv,tmp,stringvalue)
      if(nadv<0.or.nadv>2) then
        write(errmsg,*)'Unknown advection flag',nadv
        call parallel_abort(errmsg)
      endif

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

!...  Min/max. btracking step
      call get_param('dtb_min',2,itmp,dtb_min,stringvalue)
      call get_param('dtb_max',2,itmp,dtb_max,stringvalue)
      if(dtb_min>dtb_max.or.dtb_min<=0) call parallel_abort('dtb_min>dtb_max')
!'

!...  Minimum depth allowed
      call get_param('h0',2,itmp,h0,stringvalue)
      if(h0<=0) call parallel_abort('h0 must be positive')

!...  Bottom friction
      call get_param('bfric',1,nchi,tmp,stringvalue)
      if(nchi==0) then !read in drag coefficients
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
          n1=isidenode(i,1)
          n2=isidenode(i,2)
          Cd(i)=(Cdp(n1)+Cdp(n2))/2

!         Debug
!          if(myrank==0) write(99,*)i,iplg(n1),iplg(n2),Cd(i)
        enddo
        close(32)
      else if(nchi==1) then !read in roughness in meters
!       Cdmax: max. Cd
        call get_param('Cdmax',2,itmp,Cdmax,stringvalue)
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
!        do i=1,nsa
!          n1=isidenode(i,1)
!          n2=isidenode(i,2)
!          sm=min(rough_p(n1),rough_p(n2))
!          if(sm<0) then
!            rough(i)=sm !<0
!          else !both non-negative
!            rough(i)=(rough_p(n1)+rough_p(n2))/2 !>=0
!          endif
!        enddo !i
      else
        write(errmsg,*)'Unknown bfric', nchi
        call parallel_abort(errmsg)
      endif !nchi

!     Coriolis options (must be 1 if ics=2)
      call get_param('ncor',1,ncor,tmp,stringvalue)
      if(iabs(ncor)>1.or.ics==2.and.ncor/=1) then
        write(errmsg,*)'Unknown ncor',ncor,ics
        call parallel_abort(errmsg)
      endif
      if(ncor==-1) then !lattitude
        call get_param('lattitude',2,itmp,tmp,stringvalue)
        coricoef=2*omega_e*sin(tmp/180*pi)
        cori=coricoef
      else if(ncor==0) then
        call get_param('coriolis',2,itmp,coricoef,stringvalue)
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

        fc=2*omega_e*sin(sfea0)
        beta=2*omega_e*cos(sfea0)
        if(myrank==0) open(31,file='coriolis.out',status='replace')
        do i=1,nsa
          id1=isidenode(i,1)
          id2=isidenode(i,2)
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

!     Wind (nws=3: for conservation check; otherwise same as nws=2)
      call get_param('nws',1,nws,tmp,stringvalue)
      call get_param('wtiminc',2,itmp,wtiminc,stringvalue)
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
        call get_param('iwind_form',1,iwind_form,tmp,stringvalue)
        if(iwind_form/=0.and.iwind_form/=-1) then
          write(errmsg,*)'Unknown iwind_form',iwind_form
          call parallel_abort(errmsg)
        endif

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
      endif

      windfactor=1 !intialize for default
      if(nws>0) then
        call get_param('nrampwind',1,nrampwind,tmp,stringvalue)
        call get_param('iwindoff',1,iwindoff,tmp,stringvalue)
        call get_param('drampwind',2,itmp,drampwind,stringvalue)
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
      endif

!     Heat and salt conservation flags
      call get_param('ihconsv',1,ihconsv,tmp,stringvalue)      
      call get_param('isconsv',1,isconsv,tmp,stringvalue)      
      if(ihconsv<0.or.ihconsv>1.or.isconsv<0.or.isconsv>1) then
        write(errmsg,*)'Unknown ihconsv or isconsv',ihconsv,isconsv
        call parallel_abort(errmsg)
      endif
      if(isconsv/=0.and.ihconsv==0) call parallel_abort('Evap/precip model must be used with heat exchnage model')
!'
      if(ihconsv/=0.and.(nws<2.or.nws>3)) call parallel_abort('Heat budge model must have nws>=2')
!'
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
   
      if(isconsv/=0) then 
#ifndef PREC_EVAP
        write(errmsg,*)'Pls enable PREC_EVAP:',isconsv
        call parallel_abort(errmsg)
!       USE_SFLUX and USE_NETCDF are definitely enabled in Makefile when isconsv=1
#endif
      endif

!...  Turbulence closure options
      call get_param('itur',1,itur,tmp,stringvalue)
      if(itur<-2.or.itur>4) then
        write(errmsg,*)'Unknown turbulence closure model',itur
        call parallel_abort(errmsg)
      endif

      if(itur==0) then
        call get_param('dfv0',2,itmp,dfv0,stringvalue)
        call get_param('dfh0',2,itmp,dfh0,stringvalue)
        dfv=dfv0; dfh=dfh0
      else if(itur==-1) then !VVD
        open(10,file='vvd.dat',status='old')
        read(10,*) !nvrt
        do j=1,nvrt
          read(10,*)k,dfv0,dfh0
          dfv(:,j)=dfv0
          dfh(:,j)=dfh0
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
            dfv(ipgl(i)%id,:)=dfv0
            dfh(ipgl(i)%id,:)=dfh0
          endif
        enddo !i
        close(10)
        close(32)
      else if(itur==2) then !read in P&P coefficients
        call get_param('h1_pp',2,itmp,h1_pp,stringvalue)
        call get_param('h2_pp',2,itmp,h2_pp,stringvalue)
        call get_param('vdmax_pp1',2,itmp,vdmax_pp1,stringvalue)
        call get_param('vdmax_pp2',2,itmp,vdmax_pp2,stringvalue)
        call get_param('vdmin_pp1',2,itmp,vdmin_pp1,stringvalue)
        call get_param('vdmin_pp2',2,itmp,vdmin_pp2,stringvalue)
        call get_param('tdmin_pp1',2,itmp,tdmin_pp1,stringvalue)
        call get_param('tdmin_pp2',2,itmp,tdmin_pp2,stringvalue)
        if(h1_pp>=h2_pp) then
          write(errmsg,*)'h1_pp >= h2_pp in P&P'
          call parallel_abort(errmsg)
        endif
        if(vdmax_pp1<vdmin_pp1.or.vdmax_pp2<vdmin_pp2) then
          write(errmsg,*)'Wrong limits in P&P:',vdmax_pp1,vdmin_pp1,vdmax_pp2,vdmin_pp2
          call parallel_abort(errmsg)
        endif
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
!         Closure name and stability function
          call get_param('turb_met',0,itmp,tmp,mid)
          call get_param('turb_stab',0,itmp,tmp,stab)

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
            q2(i,:)=q2min
            xl(i,:)=xlmin2(i)
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

!     i.c. for T,S
      call get_param('icst',1,icst,tmp,stringvalue)
      if(icst/=1.and.icst/=2) then
        write(errmsg,*)'Unknown i.c. flag',icst
        call parallel_abort(errmsg)
      endif

!     Mean T,S profile
!     If ibcc_mean=1, ts.ic is needed, which is the same input needed when ihot=0 and icst=2.
      call get_param('ibcc_mean',1,ibcc_mean,tmp,stringvalue)
      if(ibcc_mean/=0.and.ibcc_mean/=1) then
        write(errmsg,*)'Unknown ibcc_mean flag',ibcc_mean
        call parallel_abort(errmsg)
      endif

!     Global output parameters
!     Make sure the order of optional modules is same during actual output
!     Array for storing local indices for optional modules
!     indx_out(i,j): i is model id (SED, NAPZD etc); j=1:2 is the start and end indices of each model
      allocate(indx_out(10,2),stat=istat)
      if(istat/=0) call parallel_abort('MAIN: indx_out failure')

      noutput=26+ntracers
#ifdef USE_SED 
      noutput=noutput+1+2*ntracers !27+3*ntracers
      indx_out(1,1)=27+ntracers
      indx_out(1,2)=noutput
#endif
#ifdef USE_NAPZD
      noutput=noutput+2
      indx_out(2,1)=noutput-1
      indx_out(2,2)=noutput
#endif
#ifdef USE_WWM
      noutput=noutput+25 
      indx_out(3,1)=noutput-24
      indx_out(3,2)=noutput
      !Index into out_wwm() for scalar outputs below
      allocate(indx_wwm_out(23),stat=istat)
      if(istat/=0) call parallel_abort('MAIN: indx_wwm_out failure')
!      indx_wwm_out=(/1,2,3,6,7,9,10,11,12,13,18,19,20,23,26,27,28/)
      icount=0
      do i=1,27
        if(i/=6.and.i/=7.and.i/=26.and.i/=27) then
          icount=icount+1
          if(icount>23) call parallel_abort('MAIN: indx_wwm_out')
          indx_wwm_out(icount)=i
        endif
      enddo !i
#endif

      if(noutput>mnout) then
        write(errmsg,*)'Increase mnout in the header to',noutput
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
      outfile(13)='wind.62'
      outfile(14)='wist.62'
      outfile(15)='dahv.62'
      outfile(16)='vert.63'
      outfile(17)='temp.63'
      outfile(18)='salt.63'
      outfile(19)='conc.63'
      outfile(20)='tdff.63'
      outfile(21)='vdff.63'
      outfile(22)='kine.63'
      outfile(23)='mixl.63'
      outfile(24)='zcor.63'
      outfile(25)='qnon.63'
      outfile(26)='hvel.64'
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
      variable_nm(13)='wind speed'
      variable_nm(14)='wind stress (m^2/s^2)'
      variable_nm(15)='Depth averaged horizontal velocity'
      variable_nm(16)='vertical velocity'
      variable_nm(17)='temperature in C'
      variable_nm(18)='salinity in psu'
      variable_nm(19)='density anomaly in kg/m^3'
      variable_nm(20)='eddy diffusivity in m^2/s'
      variable_nm(21)='eddy viscosity in m^2/s'
      variable_nm(22)='turbulent kinetic energy'
      variable_nm(23)='turbulent mixing length'
      variable_nm(24)='z coordinates'
      variable_nm(25)='normalized non-hydrostatic pressure'
      variable_nm(26)='horizontal velocity'

      variable_dim(1:12)='2D scalar'
      variable_dim(13:15)='2D vector'
      variable_dim(16:25)='3D scalar'
      variable_dim(26)='3D vector'
    
      do i=1,ntracers
        write(ifile_char,'(i03)')i
        ifile_char=adjustl(ifile_char)  !place blanks at end
        ifile_len=len_trim(ifile_char)  !length without trailing blanks
        outfile(26+i)='trcr_'//ifile_char(1:ifile_len)//'.63' 
        variable_nm(26+i)='Tracer #'//trim(ifile_char)
        variable_dim(26+i)='3D scalar'
      enddo !i

      indx2=26+ntracers
#ifdef USE_SED
      outfile(indx2+1)='depth.61'
      variable_nm(indx2+1)='depth in m'
      variable_dim(indx2+1)='2D scalar'
      do i=1,ntracers
         write(ifile_char,'(i03)')i
         ifile_char=adjustl(ifile_char)
         ifile_len=len_trim(ifile_char)
         outfile(indx2+1+i)='bedlu_'//ifile_char(1:ifile_len)//'.61'
         variable_nm(indx2+1+i)='Bedlu #'//trim(ifile_char)
         variable_dim(indx2+1+i)='2D scalar'
         outfile(indx2+1+ntracers+i)='bedlv_'//ifile_char(1:ifile_len)//'.61'
!'
         variable_nm(indx2+1+ntracers+i)='Bedlv #'//trim(ifile_char)
         variable_dim(indx2+1+ntracers+i)='2D scalar'
      enddo !i
      indx2=indx2+1+2*ntracers
#endif USE_SED

#ifdef USE_NAPZD
      outfile(indx2+1)='Bbdf.63'
      variable_nm(indx2+1)='Biological neglected loss'
      variable_dim(indx2+1)='3D scalar'
      outfile(indx2+2)='totN.63'
      variable_nm(indx2+2)='Total Nitrogyn'
      variable_dim(indx2+2)='3D scalar'
      indx2=indx2+2
#endif USE_NAPZD

#ifdef USE_WWM
      do i=1,25
        write(ifile_char,'(i03)')i
        ifile_char=adjustl(ifile_char)  !place blanks at end
        ifile_len=len_trim(ifile_char)  !length without trailing blanks
        variable_nm(indx2+i)='WWM #'//trim(ifile_char)
        if(i>23) then !vectors
          outfile(indx2+i)='wwm_'//ifile_char(1:ifile_len)//'.62'
          variable_dim(indx2+i)='2D vector'
        else
          outfile(indx2+i)='wwm_'//ifile_char(1:ifile_len)//'.61'
          variable_dim(indx2+i)='2D scalar'
        endif
      enddo !i
      indx2=indx2+25
#endif USE_WWM

!     Check
      if(indx2/=noutput) then
        write(errmsg,*)'MAIN: mismatch (1):',indx2,noutput
        call parallel_abort(errmsg)
      endif

!     nspool,ihfskip: output and file spools
      call get_param('nspool',1,nspool,tmp,stringvalue)
      call get_param('ihfskip',1,ihfskip,tmp,stringvalue)
      if(nspool==0.or.ihfskip==0) call parallel_abort('Zero nspool')
      if(mod(ihfskip,nspool)/=0) call parallel_abort('ihfskip/nspool /= integer')
!'
      nrec=min(ntime,ihfskip)/nspool

      do i=1,noutput
        call get_param(trim(adjustl(outfile(i))),1,iof(i),tmp,stringvalue)
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
        iof(16:25)=0
      endif !lm2d

!...  Test output parameters
      call get_param('testout',1,noutgm,tmp,stringvalue)
      if(noutgm/=0) call parallel_abort('Test output unavailable')
      
      if(noutgm/=1.and.noutgm/=0) then
        write(errmsg,*)'Unknown testout',noutgm
        call parallel_abort(errmsg)
      endif
      
!...  input information about hot start output
!...
      call get_param('hotout',1,nhot,tmp,stringvalue)
      call get_param('hotout_write',1,nhot_write,tmp,stringvalue)
      if(nhot/=0.and.nhot/=1.or.mod(nhot_write,ihfskip)/=0) then
        write(errmsg,*)'Unknown hotout or hotout_write not multiple of ihfskip',nhot,ihfskip
!'
        call parallel_abort(errmsg)
      endif

!...  JCG solver parameters
!     moitn: output interval; mxitn: max iterations; rtol: relative tolerance
      call get_param('slvr_output_spool',1,moitn,tmp,stringvalue)
      call get_param('mxitn',1,mxitn,tmp,stringvalue)
      call get_param('tolerance',2,itmp,rtol,stringvalue)

!...  Compute flux flag
      call get_param('consv_check',1,iflux,tmp,stringvalue)

!...  Interpolation flag for S,T and vel. in ELM
!     Kriging in vel: no bnd nodes/sides vel. use Kriging as the filter is not applied there
      call get_param('inter_st',1,lq,tmp,stringvalue)
      call get_param('inter_mom',1,inter_mom,tmp,stringvalue)
      if(lq<0.or.lq>2.or.inter_mom<-1.or.inter_mom>1) then
        write(errmsg,*)'Unknown interpolation flags inter_st or inter_mom:',lq,inter_mom
!'
        call parallel_abort(errmsg)
      endif
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
          n1=nm(i,1); n2=nm(i,2); n3=nm(i,3)
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
          n1=nm(i,1); n2=nm(i,2); n3=nm(i,3)
          krvel(i)=min(swild(n1),swild(n2),swild(n3))
        enddo !i
      endif

!...  Interpolation mode (1: along Z; 2: along S) - removed
!      if(lm2d) then !2D
!        interpol=2
!      else !3D
!        open(32,file='interpol.gr3',status='old')
!        read(32,*)
!        read(32,*) itmp1,itmp2
!        if(itmp1/=ne_global.or.itmp2/=np_global) &
!     &call parallel_abort('Check interpol.gr3')
!        do i=1,np_global
!          read(32,*)j,xtmp,ytmp,tmp
!          if(tmp/=1.and.tmp/=2) then
!            write(errmsg,*)'Unknown interpolation flag in interpol.gr3:',i
!            call parallel_abort(errmsg)
!          endif
!          if(ipgl(i)%rank==myrank) swild(ipgl(i)%id)=tmp
!        enddo !i
!        close(32)
!        do i=1,nea
!          n1=nm(i,1); n2=nm(i,2); n3=nm(i,3)
!          interpol(i)=min(swild(n1),swild(n2),swild(n3))
!        enddo !i
!      endif !lm2d
!
!     Make sure lqk=2 & interpol=2 are in pure S region
!      do i=1,nea
!        !if(lqk(i)==2.or.interpol(i)==2) then
!        if(interpol(i)==2) then
!          if(dp(nm(i,1))>h_s.or.dp(nm(i,2))>h_s.or.dp(nm(i,3))>h_s) then
!            write(errmsg,*)'interpol=2 must be inside pure S region:',ielg(i),interpol(i)
!!'
!            call parallel_abort(errmsg)
!          endif
!        endif
!      enddo !i

!...  Cut-off depth for option for hgrad_nodes() near bottom (like baroc. force)
      call get_param('depth_zsigma',2,itmp,h_bcc1,stringvalue)

!...  Land b.c. option (inactive)
!      read(15,*) !islip !0: free slip; otherwise no slip
      islip=0
!      if(islip/=0.and.islip/=1) then
!        write(errmsg,*)'Unknow islip:',islip
!        call parallel_abort(errmsg)
!      endif
!      if(islip==1) read(15,*) hdrag0

!...  Sponge layer for elev. & vel. (relax. factor applied to 0 elev. or uv -similar to T,S)
      call get_param('inu_elev',1,inu_elev,tmp,stringvalue)
      call get_param('inu_uv',1,inu_uv,tmp,stringvalue)
      if(inu_elev<0.or.inu_elev>1.or.inu_uv<0.or.inu_uv>1) then
        write(errmsg,*)'Check sponge inputs:',inu_elev,inu_uv
        call parallel_abort(errmsg)
      endif

      if(inu_elev==1) then
        open(10,file='elev_nudge.gr3',status='old')
        read(10,*)
        read(10,*) itmp1,itmp2
        if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check elev_nudge.gr3')
        do i=1,np_global
          read(10,*)j,xtmp,ytmp,tmp1
          if(tmp1<0.or.tmp1>1) then
            write(errmsg,*)'Wrong nudging factor at node (1):',i,tmp1
            call parallel_abort(errmsg)
          endif
          if(ipgl(i)%rank==myrank) elev_nudge(ipgl(i)%id)=tmp1
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
          if(tmp1<0.or.tmp1>1) then
            write(errmsg,*)'Wrong nudging factor at node (2):',i,tmp1
            call parallel_abort(errmsg)
          endif
          if(ipgl(i)%rank==myrank) uv_nudge(ipgl(i)%id)=tmp1
        enddo !i
        close(10)
      endif !inu_uv

!...  Nudging options
      call get_param('inu_st',1,inu_st,tmp,stringvalue)
      call get_param('step_nu',2,itmp,step_nu,stringvalue)
      call get_param('vnh1',2,itmp,vnh1,stringvalue)
      call get_param('vnh2',2,itmp,vnh2,stringvalue)
      call get_param('vnf1',2,itmp,vnf1,stringvalue)
      call get_param('vnf2',2,itmp,vnf2,stringvalue)
      if(inu_st<0.or.inu_st>2.or.step_nu<dt) then
        write(errmsg,*)'Check nudging inputs:',inu_st,step_nu
        call parallel_abort(errmsg)
      endif
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
          if(tmp1<0.or.tmp1>1.or.tmp2<0.or.tmp2>1) then
            write(errmsg,*)'Wrong nudging factor at node:',i,tmp1,tmp2
            call parallel_abort(errmsg)
          endif
          if(ipgl(i)%rank==myrank) then
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
      endif

!...  Surface min. mixing length for f.s. and max. for all; inactive 
!      read(15,*) !xlmax00

!...  Order of (vertical) integration for baroclinicity
!      call get_param('bcc_order',1,mmm,tmp,stringvalue)
      mmm=0
      if(mmm<0) call parallel_abort('bcc_order<0')
!     Pre-compute sigmap & sigma_prod for rint_lag()
      if(mmm>0) then
        do k=1,nsig
          do j=1,2*mmm+1
            if(j==1) then
              sigmap(k,j)=sigma(k)
            else
              sigmap(k,j)=sigmap(k,j-1)*sigma(k)
            endif
          enddo !j 

          if(k<nsig) then
            do l=1,k
              j1=max0(l,k-mmm)
              j2=min0(nsig,k+mmm)
              if(j1>=j2) then
                write(errmsg,*)'Weird indices:',j1,j2,k,l
                call parallel_abort(errmsg)
              endif

              do i=j1,j2
                if(abs(i-k)>4) call parallel_abort('sigma_prod index out of bound')
!'
                sigma_prod(l,k,i-k)=1
                do j=j1,j2
                  if(j/=i) sigma_prod(l,k,i-k)=sigma_prod(l,k,i-k)*(sigma(i)-sigma(j))
                enddo !j
                if(sigma_prod(l,k,i-k)==0) call parallel_abort('Impossible in sigma_prod')
!'
              enddo !i
            enddo !l
          endif !k<nsig
        enddo !k=1,nsig
      endif !mmm>0

!...  Drag formulation
      call get_param('idrag',1,idrag,tmp,stringvalue)
      if(idrag/=1.and.idrag/=2) call parallel_abort('Unknown idrag')
      if(idrag==1.and.itur>0) call parallel_abort('Linear drag requires itur<=0')
!'
      if(idrag==1.and.nchi/=0) call parallel_abort('Linear drag requires nchi=0')
!'

!...  ELAD correction option for heat exchange (inactive)
!      read(15,*) !ielad 

!...  Option to limit \hat{H} to enhance stability for large friction in shallow area
      call get_param('ihhat',1,ihhat,tmp,stringvalue)
      if(ihhat/=0.and.ihhat/=1) then
        write(errmsg,*)'Unknown ihhat:',ihhat
        call parallel_abort(errmsg)
      endif
      if(lm2d) ihhat=0

!...  Transport options: ELM or upwind
!     iupwind_t: 0: ELM; 1: upwind; 2: TVD
!     Fix iupwind_s=iupwind_t
      call get_param('iupwind_t',1,iupwind_t,tmp,stringvalue)
!      call get_param('iupwind_s',1,iupwind_s,tmp,stringvalue)
!     Reset for 2D model
      if(lm2d) iupwind_t=1

      iupwind_s=iupwind_t
      if(iupwind_t<0.or.iupwind_t>2) then
        write(errmsg,*)'Unknown iupwind:',iupwind_t,iupwind_s
        call parallel_abort(errmsg)
      endif
!      if(iupwind_t+iupwind_s==3) then
!        write(errmsg,*)'TVD cannot be combined with upwind:',iupwind_t,iupwind_s
!        call parallel_abort(errmsg)
!      endif

!     tvd_mid: model AA (my own formulation); CC (Casulli's definition of upwind ratio)
      if(iupwind_t==2) then
        call get_param('tvd_mid',0,itmp,tmp,tvd_mid)
        call get_param('flimiter',0,itmp,tmp,flimiter)
        call get_param('h_tvd',2,itmp,h_tvd,stringvalue)
      endif

!.... Blending factor for vel. in btrack (1 for internal sides; 2 for bnd sides or nodes)
      call get_param('blend_internal',2,itmp,vis_coe1,stringvalue)
      call get_param('blend_bnd',2,itmp,vis_coe2,stringvalue)
      if(vis_coe1<0.or.vis_coe1>1.or.vis_coe2<0.or.vis_coe2>1) then
        write(errmsg,*)'Illegal blend_internal or blend_bnd:',vis_coe1,vis_coe2
        call parallel_abort(errmsg)
      endif

      call get_param('shapiro',2,itmp,shapiro,stringvalue)
      if(shapiro<0.or.shapiro>0.5) then
        write(errmsg,*)'Illegal shapiro:',shapiro
        call parallel_abort(errmsg)
      endif

!     Kriging option
!     Choice of generalized covariance fucntion
      call get_param('kr_co',1,kr_co,tmp,stringvalue)
      if(kr_co<=0.or.kr_co>4) then
        write(errmsg,*)'Wrong kr_co:',kr_co
        call parallel_abort(errmsg)
      endif

!...  Max. for vel. magnitude
      call get_param('rmaxvel',2,itmp,rmaxvel,stringvalue)
      if(rmaxvel<5) then
        write(errmsg,*)'Illegal rmaxvel:',rmaxvel
        call parallel_abort(errmsg)
      endif
      !Add noise for btrack
      rmaxvel=rmaxvel*1.011

!...  min. vel for invoking btrack and for abnormal exit in quicksearch
      call get_param('velmin_btrack',2,itmp,velmin_btrack,stringvalue)
      if(velmin_btrack<=0) then
        write(errmsg,*)'Illegal velmin_btrack:',velmin_btrack
        call parallel_abort(errmsg)
      endif

!...  Add more noise (nudge) in init. nudging in btrack 
!     to avoid underflow. This should not need to be adjusted
!     normally; may need to lower it for some benchmark tests
!     Default: btrack_nudge=1.013e-3
      call get_param('btrack_nudge',2,itmp,btrack_nudge,stringvalue)
      if(btrack_nudge<=0.or.btrack_nudge>0.1) then
        write(errmsg,*)'Illegal btrack_nudge:',btrack_nudge
        call parallel_abort(errmsg)
      endif

!     Test btrack alone (1: rotating Gausshill)
      call get_param('ibtrack_test',1,ibtrack_test,tmp,stringvalue)
      if(ibtrack_test/=0.and.ibtrack_test/=1) then
        write(errmsg,*)'Illegal ibtrack_test:',ibtrack_test
        call parallel_abort(errmsg)
      endif

!...  Inundation algorithm flag (1: better algorithm for fine resolution) 
      call get_param('inunfl',1,inunfl,tmp,stringvalue)
      if(inunfl/=0.and.inunfl/=1) then
        write(errmsg,*)'Illegal inunfl:',inunfl
        call parallel_abort(errmsg)
      endif

!     write mode; not used really
      call get_param('iwrite',1,iwrite,tmp,stringvalue)

!     Elev. i.c. option (elev.ic)
      call get_param('ic_elev',1,ic_elev,tmp,stringvalue)
      if(ic_elev/=0.and.ic_elev/=1) then
        write(errmsg,*)'Illegal ic_elev:',ic_elev
        call parallel_abort(errmsg)
      endif

!     Scales for dimensioning in inter-subdomain btrack
!     mxnbt=s1_mxnbt*nmm*nvrt is the dimension of btlist (nmm is the max. of all nsa);
!     mnbt=max(nbt)*s2_mxnbt is the dimension of btsendq,bttmp,btdone 
!       (nbt is the initial # of inter-subdomain trajectories), and
!     mnbt*nnbr is the dimension of btrecvq() in routine inter_btrack (nnbr is # of nbr processes).
      call get_param('s1_mxnbt',2,itmp,s1_mxnbt,stringvalue)
      call get_param('s2_mxnbt',2,itmp,s2_mxnbt,stringvalue)
      if(s1_mxnbt<=0.or.s2_mxnbt<=0) then
        write(errmsg,*)'Illegal s[12]_mxnbt:',s1_mxnbt,s2_mxnbt
        call parallel_abort(errmsg)
      endif

!     Station output option (/=0: need station.in)
!     If ics=2, the coord. in station.in must be in lat/lon (degrees)
      call get_param('iout_sta',1,iout_sta,tmp,stringvalue)

      if(iout_sta/=0) then
        call get_param('nspool_sta',1,nspool_sta,tmp,stringvalue) !output skip
        if(nspool_sta<=0) call parallel_abort('Wrong nspool_sta')
        if(mod(nhot_write,nspool_sta)/=0) call parallel_abort('mod(nhot_write,nspool_sta)/=0')
!'

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
     &arco_sta(nout_sta,3),sta_out(nout_sta,nvar_sta),sta_out_gb(nout_sta,nvar_sta),stat=istat)
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
              n1=nm(i,nx(j,1))
              n2=nm(i,nx(j,2))
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
              arco_sta(l,1)=max(0.,min(1.,arco_sta(l,1)))
              arco_sta(l,2)=max(0.,min(1.,arco_sta(l,2)))
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
        endif
      endif !iout_sta

#ifdef USE_HA
!...  Read harmonic analysis information (Adapted from ADCIRC)
      call get_param('iharind',1,iharind,tmp,stringvalue)
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
#endif USE_HA

!...  Check parameter read in from param.in
      if(myrank==0) write(16,*)'s2_mxnbt in param.in =',s2_mxnbt
!     Almost done reading param.in

!     Inter-subdomain backtracking
      call init_inter_btrack !setup datatypes and mxnbt and _r
      allocate(btlist(mxnbt),stat=istat)
      if(istat/=0) call parallel_abort('MAIN: btlist allocation failure')

!...  Compute neighborhood for 2-tier Kriging and invert matrix for resident elements only
!     Compute ne_kr for dimensioning
      if(ics==2.and.inter_mom/=0) call parallel_abort('ics=2 and inter_mom/=0')
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
        nwild(1:3)=nm(i,1:3) !temporarily save Kriging nodes
        do j=1,3 !resident nodes
          nd=nm(i,j)
          loop14: do m=1,nnp(nd)
            nd2=inp(nd,m)
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
      allocate(itier_nd(ne_kr,0:mnei_kr),akrmat_nd(ne_kr,mnei_kr+3,mnei_kr+3), &
              &akr(mnei_kr+3,mnei_kr+3),akrp((mnei_kr+3)*(mnei_kr+4)/2),ipiv(mnei_kr+3),work4(mnei_kr+3))

!     Compute Kriging neighborhood
      do i=1,ne !resident
        if(krvel(i)/=1) cycle
 
        ie=ie_kr(i)
        nei_kr=3 !# of Kriging nodes for i
        itier_nd(ie,1:3)=nm(i,1:3) !temporarily save Kriging nodes
        do j=1,3 !resident nodes
          nd=nm(i,j)
          loop15: do m=1,nnp(nd)
            nd2=inp(nd,m)
            ! Check if present
            do l=1,nei_kr
              if(itier_nd(ie,l)==nd2) cycle loop15
            enddo !l
            ! New node
            nei_kr=nei_kr+1
            itier_nd(ie,nei_kr)=nd2
          enddo loop15 !m=1,nnp(nd)
        enddo !j=1,3 nodes
        itier_nd(ie,0)=nei_kr

!       Debug
!        if(myrank==3) write(99,*)ielg(i),itier_nd(ie,0),iplg(itier_nd(ie,1:nei_kr))
      enddo !i=1,ne
      
!...  Invert Kriging matrices
!TODO: for ics=2, use eframe
      akrmat_nd=-1.e34 !initialization for debugging
      err_max=0 !max. error in computing the inverse matices
      do k=1,ne !resident
        if(ie_kr(k)==0) cycle

        ie=ie_kr(k) !local index
        npp=itier_nd(ie,0)
        do i=1,npp
          n1=itier_nd(ie,i)
          do j=1,npp
            n2=itier_nd(ie,j)
            rr=sqrt((xnd(n1)-xnd(n2))**2+(ynd(n1)-ynd(n2))**2+(znd(n1)-znd(n2))**2)
            akr(i,j)=covar(kr_co,rr)
          enddo !j
          akr(i,npp+1)=1
          akr(i,npp+2)=xnd(n1)
          akr(i,npp+3)=ynd(n1)
        enddo !i=1,npp

        akr(npp+1,1:npp)=1
        akr(npp+2,1:npp)=xnd(itier_nd(ie,1:npp))
        akr(npp+3,1:npp)=ynd(itier_nd(ie,1:npp))
        akr((npp+1):(npp+3),(npp+1):(npp+3))=0
!        bkr(1:(npp+3),1)=0 !does not matter

!       Debug
        akrmat_nd(ie,1:(npp+3),1:(npp+3))=akr(1:(npp+3),1:(npp+3))

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
              suma=suma+akrmat_nd(ie,i,l)*akr(l,j)
            enddo !l
            if(i==j) suma=suma-1

!            if(k==22910) then
!              write(96,*)i,j,akrmat_nd(ie,i,j),akr(i,j),suma
!            endif

!            if(abs(suma)>1.e-8) write(98,*)k,i,j,suma
            if(abs(suma)>err_max) err_max=abs(suma)
          enddo !j
        enddo !i

        akrmat_nd(ie,1:(npp+3),1:(npp+3))=akr(1:(npp+3),1:(npp+3))
      enddo !k=1,ne

      write(12,*)'Max. error in inverting Kriging maxtrice= ',err_max
!...  End Kriging preparation

      if(myrank==0) write(16,*)'done reading inputs...'

!								   *
!*******************************************************************
!								   *
!	Initialization for both cold and hot start 
!								   *
!*******************************************************************
!								   *
!-------------------------------------------------------------------------------
!   Wind wave model (WWM)
!-------------------------------------------------------------------------------
      iwbl=0 !init.
#ifdef USE_WWM
!       Coupling flag
!       0: decoupled so 2 models will run independently; 
!       1: full coupled (elevation, vel, and wind are all passed to WWM); 
!       2: 1-way coupling: only R.S. from WWM feedback to SELFE
        call get_param('icou_elfe_wwm',1,icou_elfe_wwm,tmp,stringvalue)
        if(icou_elfe_wwm<0.or.icou_elfe_wwm>2) then
          write(errmsg,*)'Wrong coupling flag:',icou_elfe_wwm
          call parallel_abort(errmsg)
        endif

!       Coupling interval (# of time steps)
        call get_param('nstep_wwm',1,nstep_wwm,tmp,stringvalue)
        if(nstep_wwm<1) then
          write(errmsg,*)'Wrong coupling interval:',nstep_wwm
          call parallel_abort(errmsg)
        endif

!       Wave boundary layer option
        call get_param('iwbl',1,iwbl,tmp,stringvalue)
        if(iwbl<0.or.iwbl>1) then
          write(errmsg,*)'Wrong iwbl:',iwbl
          call parallel_abort(errmsg)
        endif
        if(iwbl/=0.and.(nchi/=1.or.icou_elfe_wwm==0)) then
          write(errmsg,*)'WBL requires nchi=1:',iwbl,nchi,icou_elfe_wwm
          call parallel_abort(errmsg)
        endif
#endif USE_WWM

!     Check compatability between 2D model and parameters
      if(lm2d) then
        if(ntracers/=0.or.nonhydro==1.or.ihdif/=0.or. &
     &ibc==0.or.ibtp==1.or.nchi==1.or.ihconsv==1.or.isconsv==1.or. &
     &itur/=0.or.ibcc_mean/=0.or.inu_st/=0.or.icst==2.or.nws==3) then
          write(errmsg,*)'Uncompatable params. for 2D model:',ntracers, &
     &nonhydro,ihdif,ibc,ibtp,nchi,ihconsv,isconsv,itur, &
     &ibcc_mean,inu_st,icst,nws
          call parallel_abort(errmsg)
        endif
      endif !lm2d

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

!								   *
!*******************************************************************
!								   *
!	Initialization for cold start alone			   *
!								   *
!*******************************************************************
!								   *
      if(ihot==0) then
!------------------------------------------------------------------
!     Parabolic bowl problem 
!      bA=(1.25**2-1)/(1.25**2+1)
!      do i=1,np
!        eta2(i)=sqrt(1-bA**2)/(1-bA)-1-(xnd(i)**2+ynd(i)**2)/1500/1500*((1-bA**2)/(1-bA)**2-1)
!      enddo


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
            read(31,*) num,xtmp,ytmp,te
            if(te<tempmin.or.te>tempmax) then
              write(errmsg,*)'Initial invalid T at',i,te
              call parallel_abort(errmsg)
            endif
            if(ipgl(i)%rank==myrank) tem0(:,ipgl(i)%id)=te
          enddo !i

          read(32,*) 
          read(32,*) !np
          do i=1,np_global
            read(32,*) num,xtmp,ytmp,sa
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
            n1=isidenode(i,1)
            n2=isidenode(i,2)
            do k=1,nvrt
              tsd(k,i)=(tem0(k,n1)+tem0(k,n2))/2
              ssd(k,i)=(sal0(k,n1)+sal0(k,n2))/2
            enddo !k
          enddo !i

          do i=1,nea
            n1=nm(i,1)
            n2=nm(i,2)
            n3=nm(i,3)
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
        read(22,*) wx1,wy1
        read(22,*) wx2,wy2
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
        read(22,*)rwild(:,:)
        do i=1,np_global
          if(ipgl(i)%rank==myrank) then
            nd=ipgl(i)%id
            windx1(nd)=rwild(i,1)
            windy1(nd)=rwild(i,2)
            pr1(nd)=rwild(i,3)
          endif
        enddo !i

        read(22,*)rwild(:,:)
        do i=1,np_global
          if(ipgl(i)%rank==myrank) then
            nd=ipgl(i)%id
            windx2(nd)=rwild(i,1)
            windy2(nd)=rwild(i,2)
            pr2(nd)=rwild(i,3)
          endif
        enddo !i
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
!        read(22,*) wx1,wy1
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
        x1 = 1.0E3 !* DTD
        open(61,file='ps.in',status='old')
        read(61,*)
        read(61,*)
        do i=1,nps
          read(61,*) iegb,xPSK,xPSQ,PRPOC,PLPOC,PDOCA,PRPON,PLPON,PDON, &
          &             PNH4,PNO3,PRPOP,PLPOP,PDOP,PPO4t,PSU,PSAt,PCOD,PDO
          if(iegl(iegb)%rank==myrank) then
            PSQ(iegl(iegb)%id)     = xPSQ
            PSK(iegl(iegb)%id)     = xPSK
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
#endif USE_ICM

!...  Read initial nudging S,T
      if(inu_st==2) then
        read(37)floatout
        read(35)floatout
        do i=1,np_global
          read(37)(swild8(j,1),j=1,nvrt)
          read(35)(swild8(j,2),j=1,nvrt)
          if(ipgl(i)%rank==myrank) then
            tnd_nu1(ipgl(i)%id,1:nvrt)=swild8(1:nvrt,1)
            snd_nu1(ipgl(i)%id,1:nvrt)=swild8(1:nvrt,2)
          endif
        enddo !i
        read(37)floatout
        read(35)floatout
        do i=1,np_global
          read(37)(swild8(j,1),j=1,nvrt)
          read(35)(swild8(j,2),j=1,nvrt)
          if(ipgl(i)%rank==myrank) then
            tnd_nu2(ipgl(i)%id,1:nvrt)=swild8(1:nvrt,1)
            snd_nu2(ipgl(i)%id,1:nvrt)=swild8(1:nvrt,2)
          endif
        enddo !i
        irec_nu=2
        time_nu=step_nu
      endif !inu_st

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
#endif USE_HA

!------------------------------------------------------------------
      endif !ihot=0

!...  Tracers; user-defined tracer part
!     This part needs T,S i.c. (tsel)
      if(ntracers>0) then
!        trel0(1,:,:)=1 
!        trel0(2,:,:)=0 
!        trel=trel0

        !open(32,file='tracer_param.in',status='old')
        !read(32,*)
        !read(32,*) flag_model
        call get_param('flag_model',1,flag_model,tmp,stringvalue)
        call get_param('flag_ic',1,flag_ic,tmp,stringvalue)

        select case(flag_model)
          case(-1) ! for generic use by users
            if(myrank==0) write(16,*)'Generic tracer transport model'
          case(0) ! for test
            if(myrank==0) write(16,*)'testing tracer transport model'
          case(1) ! Sediment model

!LLP
#ifndef USE_SED
              call parallel_abort('MAIN: need to turn on USE_SED')
#else
              if(ntracers<=0) call parallel_abort('MAIN: ntracers must be > 0')
!'
              if(ics==2) call parallel_abort('MAIN: sed. model & ics=2')

#if defined SED_MORPH
                if(imm/=0) call parallel_abort('MAIN: imm and SED_MORPH cannot be used same time')
!'
#endif SED_MORPH

!...          Initialize tracers indices and some scalars
              call initialize_scalars
              call initialize_ocean
#if defined BEDLOAD_VR 
!Ligia: uu2,vv2=0 at this pt as nodalvel has not been called
!...            Initialized htot1
!                htot1=0.d0
!                do i=1,npa
!                  if(idry(i)==1)cycle
!                  htot1(i)=dp(i)+eta2(i)
!                enddo
!----------------------------------------------------------------
! Compute depth averaged hvel for VRIJN bedload
!----------------------------------------------------------------------
                dav=0.d0
                dave=0.d0

!                do i=1,npa
!                  if(idry(i)==1) cycle
!                  do k=kbp(i),nvrt-1
!                    dav(i,1)=dav(i,1)+(uu2(k+1,i)+uu2(k,i))/2*(z(k+1,i)-z(k,i))
!                    dav(i,2)=dav(i,2)+(vv2(k+1,i)+vv2(k,i))/2*(z(k+1,i)-z(k,i))
!                  enddo !k
!                  if(htot1(i)<=h0) then
!                    write(errmsg,*)'Impossible 24n',i,htot1(i),eta2(i),dp(i),idry(i),h0
!                    call parallel_abort(errmsg)
!                  endif
!                  dav(i,1)=dav(i,1)/htot1(i)
!                  dav(i,2)=dav(i,2)/htot1(i)
!                enddo !i=1,npa
#endif BEDLOAD_VR

#if defined BEDLOAD_MPM || defined BEDLOAD_VR && defined SED_MORPH
!               Compute mass matrix
                mcoefd=0
                aux1=22./108
                aux2=7./108
                do i=1,np !residents
                  do j=1,nne(i)
                    ie=ine(i,j)
                    mcoefd(i,0)=mcoefd(i,0)+area(ie)*aux1 !diagonal
                    mcoefd(i,j)=mcoefd(i,j)+area(ie)*aux2 
                    if(isbnd(1,i)==0.and.j==nne(i)) then !internal ball
                      mcoefd(i,1)=mcoefd(i,1)+area(ie)*aux2
                    else
                      mcoefd(i,j+1)=mcoefd(i,j+1)+area(ie)*aux2
                    endif
                  enddo !j

                  !Debug
                  !tmp=sum(mcoefd(i,1:nnp(i)))
                  !write(12,*)i,isbnd(1,i),tmp/mcoefd(i,0)
                enddo !i=1,np
#endif 

!             reads sediment model inputs (sedim.in file)
              if(myrank==0) write(16,*)'Reading sedi. model parameters'
              call read_sed_input
              if(myrank==0) write(16,*)'end reading sediment parameters'
#endif USE_SED
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

!             Reads date and time (same as param.in but in different format)
              !read(32,*) day, month, year, hour, minutes, seconds
              call get_param('sim_day',2,itmp,day,stringvalue)     
              call get_param('sim_month',1,month,tmp,stringvalue)
              call get_param('sim_year',2,itmp,year,stringvalue)     
              call get_param('sim_hour',2,itmp,hour,stringvalue)     
              call get_param('sim_minute',1,minutes,tmp,stringvalue)
              call get_param('sim_second',2,itmp,seconds,stringvalue)
 
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
#endif USE_ECO
          case(3) ! Oil spill
            call parallel_abort('Oil spill model under consruction')
          case(4) ! NAPZD Spitz
#ifndef USE_NAPZD
              call parallel_abort('Need to turn on USE_NAPZD')
#else
!             Reads date and time (same as param.in but in different format)
              call get_param('sim_day',2,itmp,day,stringvalue)
              call get_param('sim_month',1,month,tmp,stringvalue)
              call get_param('sim_year',2,itmp,year,stringvalue)
              call get_param('sim_hour',2,itmp,hour,stringvalue)
              call get_param('sim_minute',1,minutes,tmp,stringvalue)
              call get_param('sim_second',2,itmp,seconds,stringvalue)
              if(myrank==0) write(16,*)'reading inputs from NAPZD model'
              call read_napzd_input
#endif USE_NAPZD
          case(5) ! ICM
            if(myrank==0) write(16,*) 'ICM model invoked'

#ifndef USE_ICM
              call parallel_abort('MAIN: need to turn on USE_ICM')
#else

!             Initialize tracers indices and some scalars
!             Reads date and time (same as param.in but in different format)
              !read(32,*) day, month, year, hour, minutes, seconds
              call get_param('sim_day',2,itmp,day,stringvalue)     
              call get_param('sim_month',1,month,tmp,stringvalue)
              call get_param('sim_year',2,itmp,year,stringvalue)     
              call get_param('sim_hour',2,itmp,hour,stringvalue)     
              call get_param('sim_minute',1,minutes,tmp,stringvalue)
              call get_param('sim_second',2,itmp,seconds,stringvalue)
 
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

!...          Reads ecological model inputs (ecosim.in file)
              if(myrank==0) write(16,*)'Reading ICM parameters inputs'

              call WQCO1(dt,rnday,NDTWQ)
              allocate(WSRP(nea),WSLP(nea),WSPB1(nea),WSPB2(nea),WSPB3(nea),turb(nea),WRea(nea),stat=istat)  !added by YC
              if(istat/=0) call parallel_abort('Failed to allocate (11)')
              call WQCO2(WSRP,WSLP,WSPB1,WSPB2,WSPB3,turb,WRea) !added by YC
              call WQinput !(time) !added by YC, still need debuging
              call wqm_out 
              if(myrank==0) write(16,*)'done reading ICM parameters'
#endif USE_ICM
          case default
            call parallel_abort('Unknown tracer model')
        end select !flag_model

!       Initial conditions (may be overwritten when hotstart)
!       This is for convenience of some models that use i.c. even for hotstart
!        trel0(1,:,:)=1 
!        trel0(2,:,:)=0 
!        trel=trel0

!       Use tr_el temporaily to store i.c. at _nodes_
        tr_el=0

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
                if(ipgl(j)%rank==myrank) tr_el(m,:,ipgl(j)%id)=tr_tmp1
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
     &nvrt-kbp(i)+1,znl(kbp(i):nvrt,i),0,z_r2(1),z_r2(nz_r2),tr_el(m,kbp(i):nvrt,i))

! 	        Extend
                do k=1,kbp(i)-1
                  tr_el(m,k,i)=tr_el(m,kbp(i),i)
                enddo !k
              enddo !i=1,npa
            enddo !m=1,ntracers

          case(3)
!	    Analytical
#ifdef USE_ECO
              call bio_init(trel0)
#endif
          case default
            call parallel_abort('MAIN: Check  flag_ic!!!')
        end select !flag_ic
!#       endif

!       Initialize Tracer
        if(flag_ic/=3)then
          do m=1,ntracers
            do i=1,nea
              n1=nm(i,1)
              n2=nm(i,2)
              n3=nm(i,3)
              do k=2,nvrt
                trel0(m,k,i)=(tr_el(m,k,n1)+tr_el(m,k,n2)+tr_el(m,k,n3)+ &
     & tr_el(m,k-1,n1)+tr_el(m,k-1,n2)+tr_el(m,k-1,n3))/6
              enddo !k
              trel0(m,1,i)=trel0(m,2,i) !mainly for hotstart format
            enddo !i
          enddo !m
        endif !flag_ic/=3
        trel=trel0

        if(myrank==0) write(16,*)'done init. tracers..'
      endif !ntracers>0
!     end user-defined tracer part

!     Finally finished param.in!
!...  Initialize GOTM for both cold and hot starts (for cde etc).
!...  For real hot start, q2, xl, dfv and dfh will use the values in hotstart.in;
!...  otherwise they will be assigned values below.
      if(itur==4) then
#ifdef USE_GOTM
          call init_turbulence(8,'gotmturb.inp',nvrt-1) !GOTM starts from level 0
          call init_tridiagonal(nvrt-1)
#endif
      endif

      if(myrank==0) write(16,*)'done initializing cold start'
      
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Hot start section
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
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
              if(flag_model==0) cycle !use i.c. for testing model
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
            q2(ip,:)=swild10(:,5)
            xl(ip,:)=swild10(:,6)
            dfv(ip,:)=swild10(:,7)
            dfh(ip,:)=swild10(:,8)
            dfq1(ip,:)=swild10(:,9)
            dfq2(ip,:)=swild10(:,10)
            qnon(:,ip)=swild10(:,11)
          endif
        enddo !i=1,np_global

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
#endif USE_HA

        ! Close hotstart file
        close(36)

        if(itur==3) then
          do i=1,npa
            do j=1,nvrt
              q2(i,j)=max(q2min,q2(i,j))
              xl(i,j)=max(xlmin2(i),xl(i,j))
            enddo
          enddo
        endif

!...    change time and iteration for forecast mode
!...    Causion: this affects all t.h. files (fort.5[0-3]) and wind files
        if(ihot==1) then
          time=0
          iths=0
        endif

        if(myrank==0) then
!          write(*,*)'hot start at time=',time,iths,' ; stack #=',ifile
          write(16,*)'hot start at time=',time,iths,' ; stack #=',ifile
        endif

!...  find position in the wind input file for nws=1,2, and read in wind[x,y][1,2]
!...  Wind vector always in lat/lon frame
        if(nws==1) then
          open(22,file='wind.th',status='old')
          rewind(22)
          ninv=time/wtiminc
          wtime1=ninv*wtiminc 
          wtime2=(ninv+1)*wtiminc 
          do it=0,ninv
            read(22,*)wx1,wy1
          enddo
          read(22,*)wx2,wy2
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
            read(22,*)rwild(:,:)
          enddo
          do i=1,np_global
            if(ipgl(i)%rank==myrank) then
              nd=ipgl(i)%id
              windx1(nd)=rwild(i,1)
              windy1(nd)=rwild(i,2)
              pr1(nd)=rwild(i,3)
            endif
          enddo !i
          read(22,*)rwild(:,:)
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
!            read(22,*)wx1,wy1
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
          x1 = 1.0E3 !* DTD
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
        do it=0,ninv
          call WQinput !(time)
        enddo
        call wqm_out

        if(myrank==0) write(16,*)'end hotstart ICM point source..'
#endif USE_ICM
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
                tnd_nu1(ipgl(i)%id,1:nvrt)=swild8(1:nvrt,1)
                snd_nu1(ipgl(i)%id,1:nvrt)=swild8(1:nvrt,2)
              endif
              if(it==irec_nu+1.and.ipgl(i)%rank==myrank) then
                tnd_nu2(ipgl(i)%id,1:nvrt)=swild8(1:nvrt,1)
                snd_nu2(ipgl(i)%id,1:nvrt)=swild8(1:nvrt,2)
              endif
            enddo !i
          enddo !it
          irec_nu=irec_nu+1
        endif !inu_st

!...    Find positions in t.h. files 
        if(nettype>0) then; do it=1,iths; read(50,*) ttt,et; enddo; endif;
        if(nfltype>0) then; do it=1,iths; read(51,*) ttt,qq; enddo; endif;
        if(ntetype>0) then; do it=1,iths; read(52,*) ttt,te; enddo; endif;
        if(nsatype>0) then; do it=1,iths; read(53,*) ttt,sal; enddo; endif;
        if(ntrtype>0) then  !added by YC
          do it=1,iths
            do m=1,ntracers
              read(300+m,*) ttt,tr
            enddo
          enddo
        endif 

!       Station output
        if(iout_sta/=0.and.myrank==0) then
          do i=1,nvar_sta
            rewind(250+i)    
            do it=1,iths
              if(iof_sta(i)==1.and.mod(it,nspool_sta)==0) read(250+i,*)
            enddo !it
          enddo !i
        endif !myrank

!...  end hot start section
      endif !ihot/=0

!...  init. eta1 (for some routines like WWM)
      eta1=eta2 
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
#endif USE_NAPZD

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
#ifdef USE_OPEN64
            !openMPI has trouble with no adv. write
            open(ichan(i),file='outputs/'//(fgb(1:lfgb)//'_'//outfile(i)),status='replace',form='BINARY')
#else
            open(ichan(i),file='outputs/'//(fgb(1:lfgb)//'_'//outfile(i)),status='replace')
#endif
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
     &             (3,(nm(m,mm),mm=1,3),m=1,ne)
      else !lat/lon
        write(10,*)np,ne,(real(xlon(m)/pi*180),real(ylat(m)/pi*180),real(dp00(m)),kbp00(m),m=1,np), &
     &             (3,(nm(m,mm),mm=1,3),m=1,ne)
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
        hmod(i)=min(dp(i),h_s)
      enddo !i
      do i=1,nsa
        n1=isidenode(i,1)
        n2=isidenode(i,2)
        dps(i)=(dp(n1)+dp(n2))/2
      enddo !i
      do i=1,nea
        dpe(i)=1.e10
        do j=1,3
          if(dpe(i)>dp(nm(i,j))) dpe(i)=dp(nm(i,j))
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
      call nodalvel(ifltype)
      if(myrank==0) write(16,*)'done computing initial nodal vel...'

!...  Compute initial density at nodes or elements
      prho=-99
      swild=0    !sed. conc.
      do i=1,npa
        if(idry(i)==1) cycle
        do k=1,nvrt
          prho(k,i)=eqstate(tnd(k,i),snd(k,i) &
! LLP
#ifdef DENSED
     &                      ,swild(1:ntracers),Srho(:)  &
#endif DENSED
! LLP end
     &                      )
        enddo !k
      enddo !i

      if(iupwind_t/=0) then
        erho=-99
        do i=1,nea
          if(idry_e(i)==1) cycle
          do k=1,nvrt
            erho(k,i)=eqstate(tsel(1,k,i),tsel(2,k,i) &
!LLP
#ifdef DENSED
         &                     ,trel(:,k,i),Srho(:)        &
#endif DENSED
!LLP end
         &                    )
          enddo !k
        enddo !i
      endif !iupwind_t

!...  Compute mean density profile at nodes or elements (using current z-coord.)
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
            q2(j,1:nvrt) = tke1d(0:(nvrt-1))
            xl(j,1:nvrt) = L1d(0:(nvrt-1))
            dfv(j,1:nvrt) = min(diffmax(j),max(diffmin(j),num1d(0:(nvrt-1))))
            dfh(j,1:nvrt) = min(diffmax(j),max(diffmin(j),nuh1d(0:(nvrt-1))))
          enddo !j
#endif
      endif !itur==4 etc

#ifdef USE_SED
!...    initialize sediment variables
!...    Bottom roughness length (m)
        if(myrank==0) write(16,*)'initialize sediment...'
        if(nchi==0) then
          do j=1,npa
            if(idry(j)==1.or.Cdp(j)==0) then
              rough_p(j)=0
            else
              rough_p(j)=(znl(kbp(j)+1,j)-znl(kbp(j),j))*exp(-0.4/sqrt(Cdp(j)))
            endif
          enddo !j
        endif !nchi

        Zob=0.d0
        do j=1,nea
          n1=nm(j,1); n2=nm(j,2); n3=nm(j,3)
          Zob(j)=(rough_p(n1)+rough_p(n2)+rough_p(n3))/3
        enddo

        if(myrank==0) write(16,*)'call sed_init'
        call sed_init
        if(myrank==0) write(16,*)'end sed_init' 

#endif USE_SED

!temp fix
!Output T at start
!      fdb='temp_0000'
!      lfdb=len_trim(fdb)
!      write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
!      open(10,file='outputs/'//fdb,status='replace')
!      write(10,*)np,nproc
!      do i=1,np
!        write(10,'(i10,2(1x,e20.14),1x,e9.3)')iplg(i),xlon(i),ylat(i),tnd(1,i)
!      enddo !i
!      close(10)
      

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

      if(myrank==0) write(16,*)'time stepping begins...',iths+1,ntime

      difnum_max_l2=0 !max. horizontal diffusion number reached by each process (check stability)
      iwbl_itmax=0 !cumulative max. of iterations for WBL (Grant-Madsen formulation) for a rank 

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      endif !first_call


!-------------------------------------------------------------------------------
!   Initialize wind wave model (WWM)
!-------------------------------------------------------------------------------
#ifdef USE_WWM
        CALL INITIALIZE_WWM()
#endif      

#ifdef USE_SWAN
      do it=istart_elfe,iend_elfe
#else
      do it=iths+1,ntime
#endif

#ifdef INCLUDE_TIMING
      wtmp1=mpi_wtime() !Forcing preparation section
#endif

      time=it*dt 

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
      else
        ramp=1
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
              nd=nm(i,j)
              dzdx=dzdx+znl(k,nd)*dl(i,j,1)
              dzdy=dzdy+znl(k,nd)*dl(i,j,2)
              dzz1=znl(k2,nd)-znl(k1,nd)
              if(k1==k2.or.dzz1<=0) then
                write(errmsg,*)'Impossible 91:',k1,k2,dzz1
                call parallel_abort(errmsg)
              endif
              dudz=dudz+(uu2(k2,nd)-uu2(k1,nd))/dzz1/3
              dvdz=dvdz+(vv2(k2,nd)-vv2(k1,nd))/dzz1/3
              dudx2=dudx2+uu2(k,nd)*dl(i,j,1) !on S-plane for k>=kz
              dudy2=dudy2+uu2(k,nd)*dl(i,j,2)
              dvdx2=dvdx2+vv2(k,nd)*dl(i,j,1)
              dvdy2=dvdy2+vv2(k,nd)*dl(i,j,2)
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
              ie=is(i,j)
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
            if(isbs(i)==-1) swild98(:,k,i)=0 !free slip land bnd
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

!       Compute horizontal viscosity
        do j=1,ns !residents only
          do k=kbs(j)+1,nvrt !viscosity = 0 at bottom
            ta=0
            do l=1,2 !element
              ie=is(j,l)
              if(ie/=0) then
                ta=ta+area(ie)
                do i=1,3 !sides
                  jsj=js(ie,i)
                  if(is(j,2)==0.or.jsj/=j) then
                    d2uv(1:2,k,j)=d2uv(1:2,k,j)+ssign(ie,i)*distj(jsj)*hvis(k,ie)*swild98(1:2,k,jsj)
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
        do jf=1,ntip
          ncyc=int(tfreq(jf)*time/2/pi)
          arg=tfreq(jf)*time-ncyc*2*pi+jspc(jf)*xlon(i)+tear(jf)
          etp(i)=etp(i)+ramp*tamp(jf)*tnf(jf)*fun_lat(i,jspc(jf))*cos(arg)
        enddo !jf
      enddo !i

!...  process new wind info 
!...  Wind vectors always in lat/lon frame 
      if(nws==1) then
        if(time>=wtime2) then
          wtime1=wtime2
          wtime2=wtime2+wtiminc
          read(22,*) wx2,wy2
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
          read(22,*)rwild(:,:) 
          do i=1,np_global
            if(ipgl(i)%rank==myrank) then
              nd=ipgl(i)%id
              windx2(nd)=rwild(i,1)
              windy2(nd)=rwild(i,2)
              pr2(nd)=rwild(i,3)
            endif
          enddo !i
        endif

        wtratio=(time-wtime1)/wtiminc
        do i=1,npa
          windx(i)=windx1(i)+wtratio*(windx2(i)-windx1(i))
          windy(i)=windy1(i)+wtratio*(windy2(i)-windy1(i))
          pr(i)=pr1(i)+wtratio*(pr2(i)-pr1(i))
        enddo !i
      endif !nws=4

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
!#ifdef USE_SFLUX
          call get_wind(wtime2,windx2,windy2,pr2,airt2,shum2)
!#endif
        endif !time>=wtime2

        wtratio=(time-wtime1)/wtiminc
        do i=1,npa
          windx(i)=windx1(i)+wtratio*(windx2(i)-windx1(i))
          windy(i)=windy1(i)+wtratio*(windy2(i)-windy1(i))
          pr(i)=pr1(i)+wtratio*(pr2(i)-pr1(i))
        enddo !i

!       Overwrite wind with wind.th
!        read(22,*)wx2,wy2
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

!...  compute wind stress components (in lat/lon frame if ics=2; in map projection E-N direction if ics=1)
      dragcmin=1.0d-3*(0.61+0.063*6)
      dragcmax=1.0d-3*(0.61+0.063*50)
      do i=1,npa
        if(nws==0) then
          tau(i,1)=0
          tau(i,2)=0
        else if(nws==1.or.nws==4.or.nws>=2.and.ihconsv==0.or.iwind_form==-1) then
          wmag=sqrt(windx(i)**2+windy(i)**2)
          dragcoef=1.0d-3*(0.61+0.063*wmag)
          dragcoef=min(max(dragcoef,dragcmin),dragcmax)
          tau(i,1)=dragcoef*0.001293*wmag*windx(i)*rampwind
          tau(i,2)=dragcoef*0.001293*wmag*windy(i)*rampwind
        else !nws>=2 and ihconsv !=0 and iwind_form=0; tauxz and tauyz defined
          if(idry(i)==1) then
            tau(i,1)=0
            tau(i,2)=0
          else !rescale as well
            tau(i,1)=-tauxz(i)/rho0*rampwind*windfactor(i)**2 !sign and scale difference between stresses tauxz and tau
            tau(i,2)=-tauyz(i)/rho0*rampwind*windfactor(i)**2
          endif
        endif !nws
      enddo !i=1,npa

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
            read(62,*) ipgb,xtmp
            if(ipgl(ipgb)%rank==myrank) then
              surf_t2(ipgl(ipgb)%id)=xtmp
            endif
          enddo
        endif !time>=wtime2+idwindrv*86400.

        stratio=(time-surf_time1)/86400.

        do i=1,npa
          surf_t(i)=surf_t1(i)+stratio*(surf_t2(i)-surf_t1(i))
        enddo !i
      endif ! iSun=2
      if(myrank==0) write(16,*)'done ICM surface T..'
#endif USE_ICM
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
              tnd_nu2(ipgl(i)%id,1:nvrt)=swild8(1:nvrt,1)
              snd_nu2(ipgl(i)%id,1:nvrt)=swild8(1:nvrt,2)
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

!...  Get new t.h. values *.th
!...
      if(nettype>0) then
        read(50,*) ttt,(ath(i),i=1,nettype)
        if(it==iths+1.and.abs(ttt-time)>1.e-4) then
          write(errmsg,*)'Starting time wrong for eta',it,ttt
          call parallel_abort(errmsg)
        endif
      
        icount=0
        do k=1,nope_global
          if(iettype(k)==1) then
            icount=icount+1
            if(icount>nettype) call parallel_abort('Wrong counting 1')
            eth(k,1)=ath(icount)
          endif
        enddo 
      endif

      if(nfltype>0) then
        read(51,*) ttt,(ath(i),i=1,nfltype)
        if(it==iths+1.and.abs(ttt-time)>1.e-4) then
          write(errmsg,*)'Starting time wrong for flux',it,ttt,time
          call parallel_abort(errmsg)
        endif

        icount=0
        do k=1,nope_global
          if(ifltype(k)==1) then
            icount=icount+1
            if(icount>nfltype) call parallel_abort('Wrong counting 2')
            qthcon(k)=ath(icount)
          endif
        enddo !k
      endif

      if(ntetype>0) then
        read(52,*) ttt,(ath(i),i=1,ntetype)
        if(it==iths+1.and.abs(ttt-time)>1.e-4) then
          write(errmsg,*)'Starting time wrong for temp',it,ttt
          call parallel_abort(errmsg)
        endif

        icount=0
        do k=1,nope_global
          if(itetype(k)==1) then
            icount=icount+1
            if(icount>ntetype) call parallel_abort('Wrong counting 3')
            tth(k,1,1)=ath(icount)
          endif
        enddo !k
      endif

      if(nsatype>0) then
        read(53,*) ttt,(ath(i),i=1,nsatype)
        if(it==iths+1.and.abs(ttt-time)>1.e-4) then
          write(errmsg,*)'Starting time wrong for salt',it,ttt
          call parallel_abort(errmsg)
        endif

        icount=0
        do k=1,nope_global
          if(isatype(k)==1) then
            icount=icount+1
            if(icount>nsatype) call parallel_abort('Wrong counting 4')
            sth(k,1,1)=ath(icount)
          endif
        enddo !k
      endif

      if(ntrtype>0) then !type I
        do m=1,ntracers
          read(300+m,*) ttt,(ath(i),i=1,ntrtype)
          if(it==iths+1.and.abs(ttt-time)>1.e-4) then
            write(errmsg,*)'Starting time wrong for tracer',it,ttt
            call parallel_abort(errmsg)
          endif

          icount=0
          do k=1,nope_global
            if(itrtype(k)==1) then
              icount=icount+1
              if(icount>ntrtype) call parallel_abort('Wrong counting 5')
              trth(m,1,1,k)=ath(icount)
            endif
          enddo !k
        enddo !# of tracers
      endif

      if(nettype2>0) then
        read(54,rec=it) floatout,(a2th(1,1,i),i=1,nnode_et)
        if(it==iths+1.and.abs(floatout-time)>1.e-4) then
          write(errmsg,*)'Starting time wrong for eta 2',it,floatout
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
              eth(k,j)=a2th(1,1,icount2)
            enddo !j
          endif
        enddo !k
      endif

      if(nfltype2>0) then
        read(58,rec=it) floatout,((a2th(1:2,l,i),l=1,nvrt),i=1,nnode_fl)
        if(it==iths+1.and.abs(floatout-time)>1.e-4) then
          write(errmsg,*)'Starting time wrong for flux 2',it,floatout
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
              uthnd(k,j,1:nvrt)=a2th(1,1:nvrt,icount2) !ll frame if ics=2
              vthnd(k,j,1:nvrt)=a2th(2,1:nvrt,icount2)
!              read(58,*)(uthnd(k,j,l),vthnd(k,j,l),l=1,nvrt)
!              if(nd/=nd2) then
!                write(11,*)'Wrong node # in uv.th',nd,nd2
!                stop
!              endif
            enddo !j
          endif
        enddo !k
      endif

      if(ntetype2>0) then
        read(56,rec=it) floatout,((a2th(1,l,i),l=1,nvrt),i=1,nnode_te) 
        if(it==iths+1.and.abs(floatout-time)>1.e-4) then
          write(errmsg,*)'Starting time wrong for temp. 2',it,floatout
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
              tth(k,j,1:nvrt)=a2th(1,1:nvrt,icount2)
!              read(56,*)nd2,(tth(k,j,l),l=1,nvrt)
!              if(nd/=nd2) then
!                write(11,*)'Wrong node # in temp3D.th',nd,nd2
!                stop
!              endif
            enddo !j
          endif
        enddo !k
      endif

      if(nsatype2>0) then
        read(57,rec=it) floatout,((a2th(1,l,i),l=1,nvrt),i=1,nnode_sa)
        if(it==iths+1.and.abs(floatout-time)>1.e-4) then
          write(errmsg,*)'Starting time wrong for salt 2',it,floatout
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
              sth(k,j,1:nvrt)=a2th(1,1:nvrt,icount2)
!              read(57,*)nd2,(sth(k,j,l),l=1,nvrt)
!              if(nd/=nd2) then
!                write(11,*)'Wrong node # in salt3D.th',nd,nd2
!                stop
!              endif
            enddo !j
          endif
        enddo !k
      endif

!     Tracers
      if(ntrtype2>0) then
        read(59,rec=it) floatout,((a2th(1:ntracers,l,i),l=1,nvrt),i=1,nnode_te) 
        if(it==iths+1.and.abs(floatout-time)>1.e-4) then
          write(errmsg,*)'Starting time wrong for tracers 2',it,floatout
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
              trth(1:ntracers,1:nvrt,j,k)=a2th(1:ntracers,1:nvrt,icount2)
            enddo !j
          endif !itetype
        enddo !k
      endif !ntrtype2

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
                ie=ine(n1,j)
                if(ie>0) then
                  do l=1,3
                    isd=js(ie,l)
                    if((isidenode(isd,1)==n1.or.isidenode(isd,2)==n1).and. &
                       (isidenode(isd,1)==n2.or.isidenode(isd,2)==n2)) then
                       isd0=isd; exit loop01
                    endif
                  enddo !l
                endif !ie>0
              end do loop01 !j=1,nne(n1)

              if(isd0==0.or.isd0>ns) cycle !skip ghost to avoid duplication

              htot=dps(isd0)+(eta2(n1)+eta2(n2))/2
              if(htot<=h0) then
                write(errmsg,*)'Dry bnd side:',htot,kk,i,iplg(n1)
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
        n1=isidenode(i,1)
        n2=isidenode(i,2)
        ibnd=isbs(i) !global bnd #
        if(ibnd<=0) cycle

!       Open bnds
!       ll frame at side
        swild10(1:3,1:3)=(pframe(:,:,n1)+pframe(:,:,n2))/2

!       Find bnd node indices for n1,n2
        nwild(1:2)=0
        do j=1,2
          do jj=1,2
            if(isbnd(jj,isidenode(i,j))==ibnd) then
              nwild(j)=isbnd(-jj,isidenode(i,j)) !global index
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
            uth(i,:)=vnth0*sframe(1,1,i)
            vth(i,:)=vnth0*sframe(2,1,i)
          else !lat/lon
            call project_hvec(vnth0,0.d0,sframe(:,:,i),swild10(1:3,1:3),utmp,vtmp)
            uth(i,:)=utmp
            vth(i,:)=vtmp
          endif !ics
          !do k=1,nvrt
          !enddo !k

          !write(12,*)it,i,iplg(n1),iplg(n2),uth(i,1),vth(i,1),sqrt(uth(i,1)**2+vth(i,1)**2)

        else if(ifltype(ibnd)==-1) then !Flather 1
!         uthnd is the normal vel.; no ramp up
          do k=1,nvrt
            if(uthnd(ibnd,nwild(1),k)<-98.or.uthnd(ibnd,nwild(2),k)<-98) then
              write(errmsg,*)'MAIN: Problem with Flather:',iplg(n1),iplg(n2)
              call parallel_abort(errmsg)
            endif
            tmp=(uthnd(ibnd,nwild(1),k)+uthnd(ibnd,nwild(2),k))/2
            if(ics==1) then
              uth(i,k)=tmp*sframe(1,1,i)
              vth(i,k)=tmp*sframe(2,1,i) 
            else !lat/lon
              call project_hvec(tmp,0.d0,sframe(:,:,i),swild10(1:3,1:3),utmp,vtmp)
              uth(i,k)=utmp
              vth(i,k)=vtmp
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
            uth(i,:)=vnth0*sframe(1,1,i)
            vth(i,:)=vnth0*sframe(2,1,i)
          else !lat/lon
            call project_hvec(vnth0,0.d0,sframe(:,:,i),swild10(1:3,1:3),utmp,vtmp)
            uth(i,:)=utmp
            vth(i,:)=vtmp
          endif !ics
        else if(iabs(ifltype(ibnd))==4) then
          do k=1,nvrt
            if(uthnd(ibnd,nwild(1),k)<-98.or.uthnd(ibnd,nwild(2),k)<-98.or. &
               vthnd(ibnd,nwild(1),k)<-98.or.vthnd(ibnd,nwild(2),k)<-98) then
              write(errmsg,*)'Wrong time series of vel.'
              call parallel_abort(errmsg)
            endif
            uth(i,k)=ramp*(uthnd(ibnd,nwild(1),k)+uthnd(ibnd,nwild(2),k))/2
            vth(i,k)=ramp*(vthnd(ibnd,nwild(1),k)+vthnd(ibnd,nwild(2),k))/2
          enddo !k

!          if(myrank==4) write(99,*)i,iplg(n1),iplg(n2),uth(i,1),vth(i,1),uth(i,nvrt),vth(i,nvrt)

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


!-------------------------------------------------------------------------------
!   Wind wave model (WWM)
!-------------------------------------------------------------------------------
#ifdef USE_WWM
!     Error: 2D not checked
      if(mod(it,nstep_wwm)==0) then
        wtmp1=mpi_wtime()
        if(myrank==0) write(16,*)'starting WWM'
        call WWM_II(it,icou_elfe_wwm,dt,nstep_wwm) !,wwave_force,out_wwm) !,windx,windy

       !&windx,windy,sdbt(1,:,1:npa),sdbt(3,:,1:npa),sdbt(2,:,1:npa),out_wwm)

!     out_wwm(npa,30): output variables from WWM (all 2D); see names in NVARS() in the routine
!                      BASIC_PARAMETER() in wwm_initio.F90 for details; below is a snapshot from there:
         !OUTPAR(1) = HS        ! Significant wave height
         !OUTPAR(2) = TM01      ! Mean average period
         !OUTPAR(3) = TM02      ! Zero down crossing period for comparison with buoy.
         !OUTPAR(4) = KLM       ! Mean wave number
         !OUTPAR(5) = WLM       ! Mean wave length
         !OUTPAR(6)  = ETOTS    ! Etot energy in y-direction
         !OUTPAR(7)  = ETOTC    ! Etot energy in x-direction
         !OUTPAR(8)  = DM       ! Mean average energy transport direction
         !OUTPAR(9)  = DSPR     ! Mean directional spreading
         !OUTPAR(10)  = FPP     ! Peak frequency (Hz)
         !OUTPAR(11)  = TPP     ! Peak period (Tp) (sec)
         !OUTPAR(12)  = CPP     ! Peak phase vel. (m/s)
         !OUTPAR(13)  = WNPP    ! Peak n-factor
         !OUTPAR(14)  = CGPP     ! Peak group vel.
         !OUTPAR(15)  = KPP      ! Peak wave number
         !OUTPAR(16)  = LPP      ! Peak wave length.
         !OUTPAR(17)  = PEAKD    ! Peak (dominant) direction (degr)
         !OUTPAR(18)  = PEAKDSPR ! Peak directional spreading
         !OUTPAR(19)  = DPEAK    ! Discrete peak direction
         !OUTPAR(20)  = UBOT     ! Orbital vel. (m/s)
         !OUTPAR(21)  = ORBITAL  ! RMS Orbital vel. (m/s)
         !OUTPAR(22)  = BOTEXPER ! Bottom excursion period.
         !OUTPAR(23)  = TMBOT    ! bottom wave period (sec)
         !OUTPAR(24) = URSELL    ! Uresell number based on peak period ...
         !OUTPAR(25)  = USTOKES   ! Bottom stokes drift
         !OUTPAR(26)  = USTOKES_X ! X-Component
         !OUTPAR(27)  = USTOKES_Y ! Y-Component.

!!     wwave_force(:,1:nsa,1:2) = Rsx, Rsy in my notes (the terms in momen. eq.)
!!     and has a dimension of m/s/s
!      wwave_force=0
!      call hgrad_nodes(0,nvrt,npa,nsa,sdbt(1,:,1:npa),dr_dxy)
!      call exchange_s3d_2(dr_dxy)
!      do i=1,nsa
!        if(idry_s(i)==0) then
!          etam=(eta2(isidenode(i,1))+eta2(isidenode(i,2)))/2
!          if(dps(i)+etam<=0) call parallel_abort('MAIN: (999)')
!          wwave_force(:,i,1)=wwave_force(:,i,1)-dr_dxy(1,:,i)/(dps(i)+etam)
!        endif
!      enddo !i
!
!      call hgrad_nodes(0,nvrt,npa,nsa,sdbt(3,:,1:npa),dr_dxy)
!      call exchange_s3d_2(dr_dxy)
!      do i=1,nsa
!        if(idry_s(i)==0) then
!          etam=(eta2(isidenode(i,1))+eta2(isidenode(i,2)))/2
!          if(dps(i)+etam<=0) call parallel_abort('MAIN: (998)')
!          wwave_force(:,i,2)=wwave_force(:,i,2)-dr_dxy(2,:,i)/(dps(i)+etam)
!        endif
!      enddo !i
!
!      call hgrad_nodes(0,nvrt,npa,nsa,sdbt(2,:,1:npa),dr_dxy)
!      call exchange_s3d_2(dr_dxy)
!
!!Debug
!!      write(12,*)'Checking R.S.',it,time
!      do i=1,nsa
!        if(idry_s(i)==0) then
!          etam=(eta2(isidenode(i,1))+eta2(isidenode(i,2)))/2
!          if(dps(i)+etam<=0) call parallel_abort('MAIN: (997)')
!          wwave_force(:,i,1)=wwave_force(:,i,1)-dr_dxy(2,:,i)/(dps(i)+etam)
!          wwave_force(:,i,2)=wwave_force(:,i,2)-dr_dxy(1,:,i)/(dps(i)+etam)
!
!!Debug: check force
!!          if(i<=ns.and.(abs(xcj(i)-1960)<1.e-4.or.abs(xcj(i)-5960)<1.e-4)) then
!!            write(12,*)i,iplg(isidenode(i,1:2)),wwave_force(:,i,1)
!!            write(12,*)i,iplg(isidenode(i,1:2)),wwave_force(:,i,2)
!!          endif
!        endif
!      enddo !i

!     Reset wave force if decoupled
!      if(icou_elfe_wwm==0) wwave_force=0

        ttime=mpi_wtime()-wtmp1
        if(myrank==0) write(16,*)'WWM-RS part took (sec) ',ttime
      endif !mod()

!Debug
!      if(it>=1856) then
!        write(12,*)'VARS in WWM:',it,out_wwm
!        write(12,*)'R.S. in WWM:',it,wwave_force
!        write(12,*)'Max. & min. R.S. from WWM (m/s/s):',it,maxval(wwave_force),minval(wwave_force)
!      endif
      
#endif USE_WWM

!...  Bottom drag coefficients for nchi=1; Cd and Cdp for nchi=0 already read in
      if(nchi==1) then !idrag=2 
        Cdp=0; Cd=0 !for dry pts
!       Drag at nodes
        ltmp1=.false. !for WBL iteration
        do i=1,npa
          if(idry(i)==1) cycle
!         Wet node
          htot=dp(i)+eta2(i)
          if(rough_p(i)<=0) then !time-independent Cd
            Cdp(i)=abs(rough_p(i))
          else !roughness >0
            bthick=znl(kbp(i)+1,i)-znl(kbp(i),i) !thickness of bottom bnd layer
           if(bthick<=rough_p(i)) then
              if(ifort12(5)==0) then
                ifort12(5)=1
                write(12,*)'BL too fine (2):',i,bthick,rough_p(i),htot
              endif
              Cdp(i)=Cdmax
            else
              Cdp(i)=1/(2.5*log(bthick/rough_p(i)))**2 

              !WBL
#ifdef USE_WWM
              if(iwbl==1) then
                vmag=sqrt(uu2(kbp(i)+1,i)**2+vv2(kbp(i)+1,i)**2)
                taubx=Cdp(i)*vmag*uu2(kbp(i)+1,i)
                tauby=Cdp(i)*vmag*vv2(kbp(i)+1,i)
                ubm=out_wwm(i,21) !orbital vel.
                wfr=2*pi/max(0.1d0,out_wwm(i,11)) !angular freq.
                wdir=out_wwm(i,17) !wave direction
                call wbl_GM(taubx,tauby,rough_p(i),ubm,wfr,wdir,z0b,fw,delta_wc,iter,ifl)
                if(ifl==2) ltmp1=.true.                
                if(iter>iwbl_itmax) iwbl_itmax=iter
                !Debug
!                if(it==1000) write(12,*)'WBL:',iplg(i),dp(i),ubm,out_wwm(i,4),wdir, &
!     &rough_p(i),z0b,iter,ifl,delta_wc

                if(bthick<=z0b) then
                  Cdp(i)=Cdmax
                else
                  Cdp(i)=1/(2.5*log(bthick/z0b))**2
                endif
              endif !iwbl             
#endif USE_WWM              

              Cdp(i)=min(Cdp(i),Cdmax)
            endif
          endif
        enddo !i=1,npa

!       Output warning for WBL if iteration didn't converge
#ifdef USE_WWM
        if(iwbl==1) then
           call mpi_reduce(ltmp1,ltmp2,1,MPI_LOGICAL,MPI_LOR,0,comm,ierr)
           if(myrank==0.and.ltmp2) write(16,*)'WBL-GM did not converge'
           if(myrank==0) write(16,*)'Cumulative max. for GM iteration for rank 0= ',iwbl_itmax
!'
        endif !iwbl
#endif USE_WWM

        do i=1,nsa
          if(idry_s(i)==1) cycle
          Cd(i)=(Cdp(isidenode(i,1))+Cdp(isidenode(i,2)))/2
        enddo !i

!       Output Cd for first step
        if(it==iths+1) then
          fdb='Cd_0000'
          lfdb=len_trim(fdb)
          write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
          open(32,file='outputs/'//trim(fdb),status='unknown')
          write(32,*)'Drag coefficents for nchi=1'
          write(32,*)nsa
          do i=1,nsa
            write(32,'(i6,2e14.6,1x,e9.3)')i,xcj(i),ycj(i),Cd(i)
          enddo !i=1,ns
          close(32)
        endif
      endif !nchi==1

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
            dfv(i,k)=vmax/(1+5*rich)**2+vmin
            dfh(i,k)=dfv(i,k)/(1+5*rich)+tmin
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
            dfv(j,:)=diffmin(j)
            dfh(j,:)=diffmin(j)
            cycle
          endif
      
!         Friction velocity: [\niu*|du/dz|]^0.5 (m/s)
          u_taus=sqrt(sqrt(tau(j,1)**2+tau(j,2)**2))
          u_taub=sqrt(Cdp(j)*(uu2(kbp(j)+1,j)**2+vv2(kbp(j)+1,j)**2))
          nlev=nvrt-kbp(j)
          do k=0,nlev 
            klev=k+kbp(j) !kbp <= klev <= nvrt
            if(k/=0) h1d(k)=znl(klev,j)-znl(klev-1,j)
!           Shear frequency squared (1/s^2): (du/dz)^2+(dv/dz)^2
!           Buoyancy frequency squared (1/s^2): -g/\rho0*(d\rho/dz))
            if(k==0.or.k==nlev) then
              if(dfv(j,klev)<=0) then
                write(errmsg,*)'Negative viscosity:',dfv(j,klev),iplg(j),klev
                call parallel_abort(errmsg)
              endif
              if(k==0) then
                SS1d(k)=u_taub**2/dfv(j,klev)
              else
                SS1d(k)=u_taus**2/dfv(j,klev)
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
            tke1d(k)=q2(j,klev)
            L1d(k)=xl(j,klev)
            if(tke1d(k)<=0.or.L1d(k)<=0) then
              write(errmsg,*)'Negative tke,mixl:',tke1d(k),L1d(k),iplg(j),klev
              call parallel_abort(errmsg)
            endif
            eps1d(k)=cde*tke1d(k)**1.5/L1d(k) 
            num1d(k)=dfv(j,klev)
            nuh1d(k)=dfh(j,klev)

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

!         Debug11
!          if(myrank==2.and.iplg(j)==14178.and.it==3253) write(98,*)(k,h1d(k),NN1d(k),SS1d(k), &
!     &num1d(k),nuh1d(k),tke1d(k),L1d(k),k=0,nlev)

          q2(j,kbp(j):nvrt) = tke1d(0:nlev)
          xl(j,kbp(j):nvrt) = L1d(0:nlev)
!          eps(i,j,:) = eps1d
          do k=0,nlev
            klev=k+kbp(j)
!           Test if they are NaN or invalid numbers
            if(num1d(k)<0.or.nuh1d(k)<0) then
              write(errmsg,*)'GOTM: problem with mixing:',num1d(k),nuh1d(k)
              call parallel_abort(errmsg)
            endif
            dfv(j,klev)=min(diffmax(j),num1d(k)+diffmin(j)) 
            dfh(j,klev)=min(diffmax(j),nuh1d(k)+diffmin(j))
          enddo !k
        enddo !j=1,npa
#endif
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
            q2(j,k)=q2min; xl(j,k)=xlmin2(j)
            dfv(j,k)=0; dfh(j,k)=0; dfq1(j,k)=0; dfq2(j,k)=0
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
          q2ha(k)=(q2(j,k)+q2(j,k-1))/2
          xlha(k)=(xl(j,k)+xl(j,k-1))/2

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
        rzbt(kbp(j))=0 !for Galperin's clipping

!        write(90,*)'WOW1',it,j

!	Compute upper bound for xl 
        do k=kbp(j),nvrt
          dists=eta2(j)-znl(k,j)
          distb=znl(k,j)+dp(j)
          if(k==kbp(j)) then
            xlmax(k)=max(xlmin2(j),dzz(k+1)*0.4_rkind)
          else if(k==nvrt) then
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
        q2fs=16.6**(2.0/3)*sqrt(tau(j,1)**2+tau(j,2)**2)/2
        q2fs=max(q2fs,q2min)
        q2bot=16.6**(2.0/3)*Cdp(j)*(uu2(kbp(j)+1,j)**2+vv2(kbp(j)+1,j)**2)/2
        q2bot=max(q2bot,q2min)
        xlfs=max(xlmin2(j),xlsc0(j)*dzz(nvrt)*0.4_rkind) 
        xlbot=max(xlmin2(j),min(2.5_rkind,xlsc0(j)*dzz(kbp(j)+1))*0.4_rkind) !"5" to prevent over-mixing

!       Debug
!        write(32,*)j,iplg(j),xlmin2(j),xlsc0(j),dzz(nvrt),xlfs
!        write(90,*)'WOW2',it,j

!	Matrix Q
        nqdim=nvrt-kbp(j)+1
        do k=kbp(j),nvrt
          kin=k-kbp(j)+1 !row #
          alow(kin)=0
          bdia(kin)=0
          cupp(kin)=0
          rrhs(kin,1)=0
          if(k<nvrt) then
            tmp=(dfq1(j,k+1)+dfq1(j,k))/2*dt/dzz(k+1)
            bdia(kin)=bdia(kin)+dzz(k+1)/3+tmp
            cupp(kin)=cupp(kin)+dzz(k+1)/6-tmp
            rrhs(kin,1)=rrhs(kin,1)+dzz(k+1)/6*(2*q2(j,k)+q2(j,k+1))
            prod=(dfv(j,k+1)+dfv(j,k))/2*shearbt(k+1)
            buoy=(dfh(j,k+1)+dfh(j,k))/2*rzbt(k+1)
            if(prod+buoy>=0) then
              rrhs(kin,1)=rrhs(kin,1)+dt*dzz(k+1)/2*(prod+buoy)
            else
              tmp=dt*dzz(k+1)/6*(prod+buoy)/q2ha(k+1)
              bdia(kin)=bdia(kin)-2*tmp
              cupp(kin)=cupp(kin)-tmp
            endif
            diss=cmiu0**3*sqrt(q2ha(k+1))/xlha(k+1)*dzz(k+1)/6 !diss/k
            bdia(kin)=bdia(kin)+dt*diss*2
            cupp(kin)=cupp(kin)+dt*diss
          endif

          if(k>kbp(j)) then
            tmp=(dfq1(j,k)+dfq1(j,k-1))/2*dt/dzz(k)
            bdia(kin)=bdia(kin)+dzz(k)/3+tmp
            alow(kin)=alow(kin)+dzz(k)/6-tmp
            rrhs(kin,1)=rrhs(kin,1)+dzz(k)/6*(2*q2(j,k)+q2(j,k-1))
            prod=(dfv(j,k)+dfv(j,k-1))/2*shearbt(k)
            buoy=(dfh(j,k)+dfh(j,k-1))/2*rzbt(k)
            if(prod+buoy>=0) then
              rrhs(kin,1)=rrhs(kin,1)+dt*dzz(k)/2*(prod+buoy)
            else
              tmp=dt*dzz(k)/6*(prod+buoy)/q2ha(k)
              bdia(kin)=bdia(kin)-2*tmp
              alow(kin)=alow(kin)-tmp
            endif
            diss=cmiu0**3*sqrt(q2ha(k))/xlha(k)*dzz(k)/6 !diss/k
            bdia(kin)=bdia(kin)+dt*diss*2
            alow(kin)=alow(kin)+dt*diss
          endif
        enddo !k=kbp(j),nvrt

!	Soln for q2 at new level
        call tridag(nvrt,100,nqdim,1,alow,bdia,cupp,rrhs,soln,gam)
        do k=kbp(j),nvrt
          kin=k-kbp(j)+1
          if(k==nvrt) then
            q2tmp(k)=q2fs
          else if(k==kbp(j)) then
            q2tmp(k)=q2bot
          else
            q2tmp(k)=max(soln(kin,1),q2min)
          endif
        enddo !k

!        write(90,*)'WOW4',it,j,(q2tmp(k),k=1,nvrt)
!        do k=1,nvrt
!          write(90,*)'Level ',k,alow(k),bdia(k),cupp(k)
!        enddo 

!	Matrix QL
        do k=kbp(j),nvrt
          kin=k-kbp(j)+1
          alow(kin)=0
          bdia(kin)=0
          cupp(kin)=0
          rrhs(kin,1)=0
          if(k<nvrt) then
            tmp=(dfq2(j,k+1)+dfq2(j,k))/2*dt/dzz(k+1)
            bdia(kin)=bdia(kin)+dzz(k+1)/3+tmp
            cupp(kin)=cupp(kin)+dzz(k+1)/6-tmp
            psi_n=cmiu0**rpub*q2(j,k)**rmub*xl(j,k)**rnub !psi^n_{j,k}
            psi_n1=cmiu0**rpub*q2(j,k+1)**rmub*xl(j,k+1)**rnub !psi^n_{j,k+1}
            rrhs(kin,1)=rrhs(kin,1)+dzz(k+1)/6*(2*psi_n+psi_n1)
            prod=cpsi1*(dfv(j,k+1)+dfv(j,k))/2*shearbt(k+1)
            buoy=cpsi3(k+1)*(dfh(j,k+1)+dfh(j,k))/2*rzbt(k+1)
            if(prod+buoy>=0) then
              rrhs(kin,1)=rrhs(kin,1)+dt*dzz(k+1)/2*(prod+buoy)*(psi_n+psi_n1)/2/q2ha(k+1)
            else
              tmp=dt*dzz(k+1)/6*(prod+buoy)/q2ha(k+1)
              bdia(kin)=bdia(kin)-2*tmp
              cupp(kin)=cupp(kin)-tmp
            endif
            diss=cpsi2p(k+1)*cmiu0**3*sqrt(q2ha(k+1))/xlha(k+1)*dzz(k+1)/6 !diss/k
            bdia(kin)=bdia(kin)+dt*diss*2
            cupp(kin)=cupp(kin)+dt*diss
          else !k=nvrt
            bdia(kin)=bdia(kin)+0.4*rnub*dt*dfq2(j,k)/xl(j,k)
          endif

          if(k>kbp(j)) then 
            tmp=(dfq2(j,k)+dfq2(j,k-1))/2*dt/dzz(k)
            bdia(kin)=bdia(kin)+dzz(k)/3+tmp
            alow(kin)=alow(kin)+dzz(k)/6-tmp
            psi_n=cmiu0**rpub*q2(j,k)**rmub*xl(j,k)**rnub !psi^n_{j,k}
            psi_n1=cmiu0**rpub*q2(j,k-1)**rmub*xl(j,k-1)**rnub !psi^n_{j,k-1}
            rrhs(kin,1)=rrhs(kin,1)+dzz(k)/6*(2*psi_n+psi_n1)
            prod=cpsi1*(dfv(j,k)+dfv(j,k-1))/2*shearbt(k)
            buoy=cpsi3(k)*(dfh(j,k)+dfh(j,k-1))/2*rzbt(k)
            if(prod+buoy>=0) then
              rrhs(kin,1)=rrhs(kin,1)+dt*dzz(k)/2*(prod+buoy)*(psi_n+psi_n1)/2/q2ha(k)
            else
              tmp=dt*dzz(k)/6*(prod+buoy)/q2ha(k)
              bdia(kin)=bdia(kin)-2*tmp
              alow(kin)=alow(kin)-tmp
            endif
            diss=cpsi2p(k)*cmiu0**3*sqrt(q2ha(k))/xlha(k)*dzz(k)/6 !diss/k
            bdia(kin)=bdia(kin)+dt*diss*2
            alow(kin)=alow(kin)+dt*diss
          else !k=kbp(j)
            bdia(kin)=bdia(kin)+0.4*rnub*dt*dfq2(j,k)/xl(j,k)
          endif
        enddo !k=kbp(j),nvrt

!        write(90,*)'WOW5',it,j
!        do k=1,nvrt
!          write(90,*)'Level ',k,alow(k),bdia(k),cupp(k)
!        enddo 

!	Soln for q2l and xl at new level
        call tridag(nvrt,100,nqdim,1,alow,bdia,cupp,rrhs,soln,gam)

!        write(90,*)'WOW6',it,j

        do k=kbp(j),nvrt
          kin=k-kbp(j)+1
          q2l=max(soln(kin,1),psimin)
          if(k==nvrt) then
            xltmp(k)=xlfs
          else if(k==kbp(j)) then
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

          q2(j,k)=q2tmp(k)
          xl(j,k)=xltmp(k)
          if(q2(j,k)<0) then
            write(errmsg,*)'Negative q2',q2(j,k),xl(j,k)
            call parallel_abort(errmsg)
          endif

!         Compute vertical diffusivities at new time
          call asm(j,k,vd,td,qd1,qd2)
          dfv(j,k)=min(diffmax(j),max(diffmin(j),vd))
          dfh(j,k)=min(diffmax(j),max(diffmin(j),td))
          dfq1(j,k)=min(diffmax(j),max(diffmin(j),qd1))
          dfq2(j,k)=min(diffmax(j),max(diffmin(j),qd2))

!         Debug
!          write(90,*)'No. ',k,xl(j,k),dfh(j,k),dfv(j,k),dfq1(j,k),dfq2(j,k)
        enddo !k=kbp(j),nvrt

!       Extend
        do k=1,kbp(j)-1
          q2(j,k)=q2(j,kbp(j))
          xl(j,k)=xl(j,kbp(j))
          dfv(j,k)=dfv(j,kbp(j))
          dfh(j,k)=dfh(j,kbp(j))
          dfq1(j,k)=dfq1(j,kbp(j))
          dfq2(j,k)=dfq2(j,kbp(j))
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
        if(it==iths+1) then
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
              nd=nm(i,j)
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
     &0._rkind,0._rkind,cspline_ypp_nd(1,kbp(i):nvrt,i))
          call cubic_spline(nvrt-kbp(i)+1,znl(kbp(i):nvrt,i),snd(kbp(i):nvrt,i), &
     &0._rkind,0._rkind,cspline_ypp_nd(2,kbp(i):nvrt,i))
        enddo !i
        do i=1,nsa
          if(idry_s(i)==1) cycle
          call cubic_spline(nvrt-kbs(i)+1,zs(kbs(i):nvrt,i),tsd(kbs(i):nvrt,i), &
     &0._rkind,0._rkind,cspline_ypp_sd(1,kbs(i):nvrt,i))
          call cubic_spline(nvrt-kbs(i)+1,zs(kbs(i):nvrt,i),ssd(kbs(i):nvrt,i), &
     &0._rkind,0._rkind,cspline_ypp_sd(2,kbs(i):nvrt,i))
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
        ptbt(3,:,i)=dfv(i,:)
        ptbt(4,:,i)=dfh(i,:)
      enddo !i

      do i=1,ns
        sdbt(1,:,i)=su2(:,i)
        sdbt(2,:,i)=sv2(:,i)
        sdbt(3,:,i)=tsd(:,i)
        sdbt(4,:,i)=ssd(:,i)
      enddo !i

      webt=we

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
              ie=ine(nd0,m)
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
              ie=is(isd0,m)
              if(ie/=0) then; if(idry_e(ie)==0) then
                ie0=ie; exit
              endif; endif
            enddo !m
            if(ie0==0) then
              write(errmsg,*)'MAIN: btrack finds no init. element (2):',iplg(isidenode(isd0,1:2))
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
              n1=isidenode(isd0,1)
              n2=isidenode(isd0,2)
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
              iadvf=min(iadv(nm(i,1)),iadv(nm(i,2)),iadv(nm(i,3)))
              xt=xctr(ie0)
              yt=yctr(ie0)
              zt=ze(j,ie0)
              uuint=(su2(j,js(ie0,1))+su2(j,js(ie0,2))+su2(j,js(ie0,3)))/3
              vvint=(sv2(j,js(ie0,1))+sv2(j,js(ie0,2))+sv2(j,js(ie0,3)))/3
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
                if(is(isd0,2)==0) then
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
                ptbt(3,j,nd0)=dfv(nd0,j)
                ptbt(4,j,nd0)=dfh(nd0,j)
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
                  ie=ine(nd0,ii)
                  id=iself(nd0,ii)
                  dudx=0; dudy=0; dvdx=0; dvdy=0
                  do jj=1,3
                    dudx=dudx+ufg(j,ie,jj)*dl(ie,jj,1) !not strictly along z; in element frame for ics=2
                    dudy=dudy+ufg(j,ie,jj)*dl(ie,jj,2) !not strictly along z
                    dvdx=dvdx+vfg(j,ie,jj)*dl(ie,jj,1) 
                    dvdy=dvdy+vfg(j,ie,jj)*dl(ie,jj,2) 
                  enddo !jj
                  suma=suma+dt*sqrt(dudx**2+dudy**2+dvdx**2+dvdy**2)/nne(nd0)
                enddo !ii=1,nne(nd0)
                ndelt=max0(ndelt_min,min0(ndelt_max,int(suma)*4)) !>=1
              else if(l==2) then !sides
                suma=0
                icount=0
                do ii=1,2
                  ie=is(isd0,ii)
                  if(ie==0) cycle
                  icount=icount+1

                  dudx=0; dudy=0; dvdx=0; dvdy=0
                  do jj=1,3
                    dudx=dudx+ufg(j,ie,jj)*dl(ie,jj,1) !not strictly along z; in element frame for ics=2
                    dudy=dudy+ufg(j,ie,jj)*dl(ie,jj,2) !not strictly along z
                    dvdx=dvdx+vfg(j,ie,jj)*dl(ie,jj,1)
                    dvdy=dvdy+vfg(j,ie,jj)*dl(ie,jj,2)
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
                  dudx=dudx+ufg(j,ie0,jj)*dl(ie0,jj,1) !not strictly along z
                  dudy=dudy+ufg(j,ie0,jj)*dl(ie0,jj,2) !not strictly along z
                  dvdx=dvdx+vfg(j,ie0,jj)*dl(ie0,jj,1)
                  dvdy=dvdy+vfg(j,ie0,jj)*dl(ie0,jj,2)
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
     &dtbk,vis_coe,time_rm,time_rm2,uuint,vvint,wwint,nnel,jlev,xt,yt,zt,swild3,lbt)

              if(lbt) then !Backtracking exits augmented subdomain
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
                  if(1==2.and.ics==2.and.j==nvrt) then
                    call project_pt('l2g',xt,yt,0.d0,swild(1:3),swild10(1:3,1:3),xt4,yt4,zt4)
                    !Exact
                    !coorind. and lat/lon in the rotated frame
                    hatx=xnd(nd0)*cos(alpha_zonal)+znd(nd0)*sin(alpha_zonal)
                    hatz=-xnd(nd0)*sin(alpha_zonal)+znd(nd0)*cos(alpha_zonal)
                    call compute_ll(hatx,ynd(nd0),hatz,rlon,rlat)
                    rlam=rlon-omega_zonal*dt !btrack
                    hatxex=rearth*cos(rlam)*cos(rlat)
                    hatyex=rearth*sin(rlam)*cos(rlat)
                    hatzex=rearth*sin(rlat)
                    !coord. in the original frame
                    xex=hatxex*cos(alpha_zonal)-hatzex*sin(alpha_zonal)
                    yex=hatyex
                    zex=hatxex*sin(alpha_zonal)+hatzex*cos(alpha_zonal)
                    dis=sqrt((xt4-xex)**2+(yt4-yex)**2+(zt4-zex)**2)
                    !rotated ll frame at foot; WARINING: assume nvrt>=3
                    swild2(1,1)=-sin(rlam)*cos(alpha_zonal)
                    swild2(2,1)=cos(rlam)
                    swild2(3,1)=-sin(rlam)*sin(alpha_zonal)
                    swild2(1,2)=-cos(rlam)*sin(rlat)*cos(alpha_zonal)-cos(rlat)*sin(alpha_zonal)
                    swild2(2,2)=-sin(rlam)*sin(rlat)
                    swild2(3,2)=-cos(rlam)*sin(rlat)*sin(alpha_zonal)+cos(rlat)*cos(alpha_zonal)
                    call cross_product(swild2(1,1),swild2(2,1),swild2(3,1), &
     &swild2(1,2),swild2(2,2),swild2(3,2),swild2(1,3),swild2(2,3),swild2(3,3))
                    call project_hvec(uuint,vvint,swild10(1:3,1:3),swild2(1:3,1:3),u2,v2)
                    uzonal=u00_zonal*cos(rlat)
                    vdis=dsqrt((u2-uzonal)**2+v2*v2)
                    write(12,*)'Node ',iplg(nd0),j,xt4,yt4,zt4,xex,yex,zex,dis,xnd(nd0),ynd(nd0),znd(nd0)
                    write(12,*)'Node vel ',j,u2,v2,uzonal,0,vdis
                    !Exact T (cosine_bell test)
                    !ll of foot (original frame)
                    call compute_ll(xex,yex,zex,rlon0,rlat0)
                    rrr=rearth*acos(cos(rlat0)*cos(rlon0+pi/2))
                    if(rrr<rearth/3) then
                      tex=500*(1+cos(pi*rrr/rearth*3))
                    else
                      tex=0
                    endif
                    ter=swild3(1)-tex
                    write(12,*)'Node T:',iplg(nd0),jlev,swild3(1),tex,ter,swild3(2)

                    !const. zonal flow
                    !xex=-uu2(nvrt,nd0)*dt !in pframe
                    !yex=-vv2(nvrt,nd0)*dt !in pframe
                    !dis=sqrt((xt-xex)**2+(yt-yex)**2)
                    !call project_hvec(uuint,vvint,pframe(:,:,nd0),pframe(:,:,nd0),utmp,vtmp)
                    !vdis=sqrt((utmp-uzonal)**2+(vtmp-vmer)**2)
                    !write(12,*)'Node level:',iplg(nd0),j,jlev !,xt,yt,xex,yex,dis,utmp,vtmp,vdis
 
                    if(dis>p_dis_max) p_dis_max=dis
                    if(vdis>p_vdis_max) p_vdis_max=vdis
                    if(abs(ter)>p_T_max) p_T_max=abs(ter)
                  endif !zonal flow

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
                  if(1==2.and.ics==2.and.j==nvrt) then
                    n1=isidenode(isd0,1)
                    n2=isidenode(isd0,2)
                    call project_pt('l2g',xt,yt,0.d0,swild(1:3),swild10(1:3,1:3),xt4,yt4,zt4)
                    !Exact
                    !coorind. and lat/lon in the rotated frame
                    hatx=xcj(isd0)*cos(alpha_zonal)+zcj(isd0)*sin(alpha_zonal)
                    hatz=-xcj(isd0)*sin(alpha_zonal)+zcj(isd0)*cos(alpha_zonal)
                    call compute_ll(hatx,ycj(isd0),hatz,rlam,rlat)
                    rlam=rlam-omega_zonal*dt
                    !Have to reduce radius b/c initially sidecenter is not on earth surface
                    rr0=sqrt(xcj(isd0)**2+ycj(isd0)**2+zcj(isd0)**2)
                    hatxex=rr0*cos(rlam)*cos(rlat)
                    hatyex=rr0*sin(rlam)*cos(rlat)
                    hatzex=rr0*sin(rlat)
                    !coord. in the original frame
                    xex=hatxex*cos(alpha_zonal)-hatzex*sin(alpha_zonal)
                    yex=hatyex
                    zex=hatxex*sin(alpha_zonal)+hatzex*cos(alpha_zonal)
                    dis=sqrt((xt4-xex)**2+(yt4-yex)**2+(zt4-zex)**2)
                    !rotated ll frame at foot; WARINING: assume nvrt>=3
                    swild2(1,1)=-sin(rlam)*cos(alpha_zonal)
                    swild2(2,1)=cos(rlam)
                    swild2(3,1)=-sin(rlam)*sin(alpha_zonal)
                    swild2(1,2)=-cos(rlam)*sin(rlat)*cos(alpha_zonal)-cos(rlat)*sin(alpha_zonal)
                    swild2(2,2)=-sin(rlam)*sin(rlat)
                    swild2(3,2)=-cos(rlam)*sin(rlat)*sin(alpha_zonal)+cos(rlat)*cos(alpha_zonal)
                    call cross_product(swild2(1,1),swild2(2,1),swild2(3,1), &
     &swild2(1,2),swild2(2,2),swild2(3,2),swild2(1,3),swild2(2,3),swild2(3,3))
                    call project_hvec(uuint,vvint,swild10(1:3,1:3),swild2(1:3,1:3),u2,v2)
                    uzonal=u00_zonal*cos(rlat)
                    vdis=dsqrt((u2-uzonal)**2+v2*v2)
                    write(12,*)'Side ',iplg(isidenode(isd0,:)),j,xt4,yt4,zt4,xex,yex,zex,dis,&
     &xcj(isd0),ycj(isd0),zcj(isd0),rr0
                    write(12,*)'Side vel ',j,u2,v2,uzonal,0,vdis
                    !Exact T
                    !ll of foot (original frame)
                    call compute_ll(xex,yex,zex,rlon0,rlat0)
                    rrr=rearth*acos(cos(rlat0)*cos(rlon0+pi/2))
                    if(rrr<rearth/3) then
                      tex=500*(1+cos(pi*rrr/rearth*3))
                    else
                      tex=0
                    endif
                    ter=swild3(1)-tex
                    write(12,*)'Side T:',jlev,swild3(1),tex,ter,swild3(2)
                    if(dis>s_dis_max) s_dis_max=dis
                    if(vdis>s_vdis_max) s_vdis_max=vdis
                    if(abs(ter)>s_T_max) s_T_max=abs(ter)
                  endif !zonal flow

                else !element
                  if(ics==2) call parallel_abort('MAIN: why am I here?')
                  if(iadvf==0) then
                    webt(j,ie0)=we(j,ie0)
                  else
                    webt(j,ie0)=wwint
                  endif
                endif 
              endif !lbt
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
        lbt=(nbtrk/=0)
#ifdef INCLUDE_TIMING
        cwtmp=mpi_wtime()
#endif
        call mpi_allreduce(lbt,lbtgb,1,MPI_LOGICAL,MPI_LOR,comm,ierr)
#ifdef INCLUDE_TIMING
        wtimer(4,2)=wtimer(4,2)+mpi_wtime()-cwtmp
#endif
        if(ierr/=MPI_SUCCESS) call parallel_abort('MAIN: allreduce lbtgb',ierr)
!'

        if(lbtgb) then
          if(myrank==0) write(16,*)'starting inter-subdomain btrack'
          call inter_btrack(it,nbtrk,btlist) !all ranks participate
          if(myrank==0) write(16,*)'done inter-subdomain btrack'
        endif

        if(lbt) then !handle returned inter-subdomain trajectories
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
        endif !lbt
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
!        dfv(i,:)=ptbt(3,:,i)
!        dfh(i,:)=ptbt(4,:,i)
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
!        n1=isidenode(i,1)
!        n2=isidenode(i,2)
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
!        write(10,*)'Side',i,iplg(isidenode(i,1:2))
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
          isd=js(i,j)
          if(idry_s(isd)==0) then
            do k=kbs(isd)+1,nvrt !layer
              if(ics==1) then
                bubt(i,1)=bubt(i,1)+(sdbt(1,k,isd)+sdbt(1,k-1,isd))/2*(zs(k,isd)-zs(k-1,isd))*area(i)/3
                bubt(i,2)=bubt(i,2)+(sdbt(2,k,isd)+sdbt(2,k-1,isd))/2*(zs(k,isd)-zs(k-1,isd))*area(i)/3
              else 
                call project_hvec(sdbt(1,k,isd),sdbt(2,k,isd),sframe(:,:,isd),eframe(:,:,i),u2,v2)
                call project_hvec(sdbt(1,k-1,isd),sdbt(2,k-1,isd),sframe(:,:,isd),eframe(:,:,i),u1,v1)
                bubt(i,1)=bubt(i,1)+(u2+u1)/2*(zs(k,isd)-zs(k-1,isd))*area(i)/3 !in eframe
                bubt(i,2)=bubt(i,2)+(v2+v1)/2*(zs(k,isd)-zs(k-1,isd))*area(i)/3 !in eframe
              endif !ics
            enddo !k
          endif
        enddo !j
      enddo !i=1,nea

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
            if(iettype(i)==1.or.iettype(i)==2) then
              elbc(ip)=ramp*eth(i,1)
            else if(iettype(i)==3) then
              elbc(ip)=0 !initialize
              do jfr=1,nbfr
                ncyc=int(amig(jfr)*time/2/pi)
                arg=amig(jfr)*time-ncyc*2*pi+face(jfr)-efa(i,j,jfr)
                elbc(ip)=elbc(ip)+ramp*ff(jfr)*emo(i,j,jfr)*cos(arg)
              enddo !jfr=1,nbfr
            else if(iettype(i)==4) then
              elbc(ip)=ramp*eth(i,j)
            endif
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
          bigu(i,1)=0
          bigu(i,2)=0
          cycle
        endif

!	Wet side
        n1=isidenode(i,1)
        n2=isidenode(i,2)
        if(idrag==1) then
          chi(i)=Cd(i)
        else
          chi(i)=Cd(i)*sqrt(sdbt(1,kbs(i)+1,i)**2+sdbt(2,kbs(i)+1,i)**2)
        endif

        if(lm2d) then
          hhat(i)=(eta2(n1)+eta2(n2))/2+dps(i)+chi(i)*dt
          if(hhat(i)<=0) then
            write(errmsg,*)'Impossible dry 53:',hhat(i),iplg(isidenode(i,1:2))
            call parallel_abort(errmsg)
          endif
        else
          hhat(i)=(eta2(n1)+eta2(n2))/2+dps(i)-chi(i)*dt
!	  Enforce positivity for 3D model
          if(ihhat==1) hhat(i)=max(0._rkind,hhat(i))
        endif

!	bigu1,2 (in sframe if ics=2)
        bigu(i,1)=0 !U^n_x
        bigu(i,2)=0 !U^n_y
        do k=kbs(i),nvrt-1
          bigu(i,1)=bigu(i,1)+(zs(k+1,i)-zs(k,i))*(su2(k,i)+su2(k+1,i))/2
          bigu(i,2)=bigu(i,2)+(zs(k+1,i)-zs(k,i))*(sv2(k,i)+sv2(k+1,i))/2
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
          !dr_dxy in sframe if ics=2
          call hgrad_nodes(0,nvrt,npa,nsa,prho(:,1:npa)-rho_mean(:,1:npa),dr_dxy)  

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
              ie=ic3(i,j)
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
                ie=ic3(i,j)
                n1=nm(i,nx(j,1))
                n2=nm(i,nx(j,2))
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
                ie=is(i,j) 
                if(ie==0.or.idry_e(max(1,ie))==1) cycle

!               Wet element
                icount=icount+1
                gam(kbe(ie)+1:nvrt)=(ze(kbe(ie):nvrt-1,ie)+ze(kbe(ie)+1:nvrt,ie))/2
                !dr_dxy at elements and half levels; eframe if ics=2
                rrhs(kbe(ie)+1:nvrt,1)=dr_dxy(1,kbe(ie)+1:nvrt,ie)
                rrhs(kbe(ie)+1:nvrt,2)=dr_dxy(2,kbe(ie)+1:nvrt,ie)
                call vinter(nvrt,100,2,(zs(k,i)+zs(k-1,i))/2,kbe(ie)+1,nvrt,k,gam,rrhs,swild,ibelow)
                if(ics==2) then !to sframe
                  call project_hvec(swild(1),swild(2),eframe(:,:,ie),sframe(:,:,i),drdx,drdy)
                  swild(1)=drdx
                  swild(2)=drdy
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
        call hgrad_nodes(0,nvrt,npa,nsa,qnon,dqnon_dxy)
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
      do i=1,nea
        if(idry_e(i)==1) then
          ghat1(i,1)=0
          ghat1(i,2)=0
          cycle
        endif

!	Wet elements
!	Excluding hvis, baroclinc force first
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
          nd=nm(i,j)
!         idry_e(i) checked already
          detadx=detadx+eta2(nd)*dl(i,j,1) !in eframe if ics=2
          detady=detady+eta2(nd)*dl(i,j,2)
          dprdx=dprdx+pr(nd)*dl(i,j,1)
          dprdy=dprdy+pr(nd)*dl(i,j,2)
          if(dpe(i)>=tip_dp) then
            detpdx=detpdx+etp(nd)*dl(i,j,1)
            detpdy=detpdy+etp(nd)*dl(i,j,2)
          endif
          h_bar=h_bar+(dp(nd)+eta2(nd))/3

          isd=js(i,j)
          chigamma=chigamma+chi(isd)/3
          hhat_bar=hhat_bar+hhat(isd)/3
          if(ics==1) then
            tau_x=tau_x+tau(nd,1)/3
            tau_y=tau_y+tau(nd,2)/3
            ubstar=ubstar+chi(isd)*sdbt(1,kbs(isd)+1,isd)/3 
            vbstar=vbstar+chi(isd)*sdbt(2,kbs(isd)+1,isd)/3
            bigf1=bigf1+cori(isd)*bigu(isd,2)/3
            bigf2=bigf2-cori(isd)*bigu(isd,1)/3
            botf1=botf1+chi(isd)*cori(isd)*sv2(kbs(isd)+1,isd)/3
            botf2=botf2-chi(isd)*cori(isd)*su2(kbs(isd)+1,isd)/3
          else !lat/lon; convert to eframe 
            call project_hvec(tau(nd,1),tau(nd,2),pframe(:,:,nd),eframe(:,:,i),taux2,tauy2)
            tau_x=tau_x+taux2/3 !in eframe
            tau_y=tau_y+tauy2/3
            call project_hvec(sdbt(1,kbs(isd)+1,isd),sdbt(2,kbs(isd)+1,isd),sframe(:,:,isd),eframe(:,:,i),ub2,vb2)
            call project_hvec(bigu(isd,1),bigu(isd,2),sframe(:,:,isd),eframe(:,:,i),bigu2,bigv2)
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
          ghat1(i,1)=bubt(i,1)+area(i)*dt*(bigf1+tau_x-ubstar-dt*botf1-grav*(1-thetai)*hhat_bar*detadx)
          ghat1(i,2)=bubt(i,2)+area(i)*dt*(bigf2+tau_y-vbstar-dt*botf2-grav*(1-thetai)*hhat_bar*detady)
        else !2D
          !Compute element averages (eframe)
          av_elem_x=0
          av_elem_y=0
          do j=1,3 !wet side
            isd=js(i,j) 
            n1=isidenode(isd,1); n2=isidenode(isd,2)
!            etam=(eta2(n1)+eta2(n2))/2
            htot=zs(nvrt,isd)-zs(kbs(isd),isd) !etam+dps(isd)
            if(hhat(isd)<=0) then
              write(errmsg,*)'Impossible dry 51:',hhat(isd),ielg(i),iplg(isidenode(isd,1:2))
              call parallel_abort(errmsg)
            endif
            !\hat{Gamma} (in eframe if ics=2)
            if(ics==1) then
              sdbtu=sdbt(1,1,isd)
              sdbtv=sdbt(2,1,isd)
              u2=su2(1,isd)
              v2=sv2(1,isd)
              taux2=(tau(n1,1)+tau(n2,1))/2
              tauy2=(tau(n1,2)+tau(n2,2))/2
            else !lat/lon; project to eframe
              call project_hvec(sdbt(1,1,isd),sdbt(2,1,isd),sframe(:,:,isd),eframe(:,:,i),sdbtu,sdbtv)
              call project_hvec(su2(1,isd),sv2(1,isd),sframe(:,:,isd),eframe(:,:,i),u2,v2)
              swild10(1:3,1:3)=(pframe(:,:,n1)+pframe(:,:,n2))/2
              call project_hvec((tau(n1,1)+tau(n2,1))/2,(tau(n1,2)+tau(n2,2))/2, &
     &swild10(1:3,1:3),eframe(:,:,i),taux2,tauy2)
            endif !ics
            hat_gam_x=sdbtu+dt*(-dprdx/rho0+0.69*grav*detpdx)+ &
     &(1-theta2)*cori(isd)*dt*v2-grav*(1-thetai)*dt*detadx+dt*taux2/htot
            hat_gam_y=sdbtv+dt*(-dprdy/rho0+0.69*grav*detpdy)- &
     &(1-theta2)*cori(isd)*dt*u2-grav*(1-thetai)*dt*detady+dt*tauy2/htot
#ifdef USE_WWM
            hat_gam_x=hat_gam_x+dt*wwave_force(1,isd,1)
            hat_gam_y=hat_gam_y+dt*wwave_force(1,isd,2)
#endif USE_WWM
            del=hhat(isd)**2+(theta2*cori(isd)*dt*htot)**2 !delta>0
            !Gamma vector (eframe)
            gam_x=(hhat(isd)*htot*hat_gam_x+theta2*cori(isd)*htot*htot*dt*hat_gam_y)/del
            gam_y=(hhat(isd)*htot*hat_gam_y-theta2*cori(isd)*htot*htot*dt*hat_gam_x)/del
            av_elem_x=av_elem_x+htot*gam_x/3
            av_elem_y=av_elem_y+htot*gam_y/3
          enddo !j
          ghat1(i,1)=area(i)*av_elem_x
          ghat1(i,2)=area(i)*av_elem_y
        endif !lm2d

!       Horizontal viscosity
        horx=0 
        hory=0
        do j=1,3 !side
          isd=js(i,j)
          do k=kbs(isd)+1,nvrt
            horx=horx+area(i)/3*(zs(k,isd)-zs(k-1,isd))*(d2uv(1,k,isd)+d2uv(1,k-1,isd))/2
            hory=hory+area(i)/3*(zs(k,isd)-zs(k-1,isd))*(d2uv(2,k,isd)+d2uv(2,k-1,isd))/2
          enddo !k
          horx=horx-dt*chigamma*area(i)/3*d2uv(1,kbs(isd)+1,isd)
          hory=hory-dt*chigamma*area(i)/3*d2uv(2,kbs(isd)+1,isd)
        enddo !j=1,3
        ghat1(i,1)=ghat1(i,1)+dt*horx
        ghat1(i,2)=ghat1(i,2)+dt*hory

!       Radiation stress (3D only; 2D part has been done above)
#ifdef USE_WWM
          if(.not.lm2d) then !3D
            rs1=0 
            rs2=0
            do j=1,3 !side
              isd=js(i,j)
              do k=kbs(isd)+1,nvrt
                rs1=rs1+area(i)/3*(zs(k,isd)-zs(k-1,isd))* &
     &(wwave_force(k,isd,1)+wwave_force(k-1,isd,1))/2
                rs2=rs2+area(i)/3*(zs(k,isd)-zs(k-1,isd))* &
     &(wwave_force(k,isd,2)+wwave_force(k-1,isd,2))/2
              enddo !k
              rs1=rs1-dt*chigamma*area(i)/3*wwave_force(kbs(isd)+1,isd,1)
              rs2=rs2-dt*chigamma*area(i)/3*wwave_force(kbs(isd)+1,isd,2)
            enddo !j=1,3
            ghat1(i,1)=ghat1(i,1)+dt*rs1
            ghat1(i,2)=ghat1(i,2)+dt*rs2
          endif !lm2d
#endif USE_WWM

!       Baroclinic force
        if(ibc==0) then
          if(prho(1,nm(i,1))<-98.or.prho(1,nm(i,2))<-98.or.prho(1,nm(i,3))<-98) then
            write(errmsg,*)'Impossible dry 5'
            call parallel_abort(errmsg)
          endif

          if(iupwind_t==0) then
!           ELM option
            swild(1:2)=0 !averaged bottom bcc (in eframe if ics=2)
            do j=1,3 !side
              isd=js(i,j)
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
              ghat1(i,1)=ghat1(i,1)+dt*area(i)*bigfc1/3
              ghat1(i,2)=ghat1(i,2)+dt*area(i)*bigfc2/3
              swild(1:2)=swild(1:2)+swild10(kbs(isd)+1,1:2)/3 !bcc(1:2,kbs(isd)+1,isd)/3
            enddo !side
            ghat1(i,1:2)=ghat1(i,1:2)-area(i)*chigamma*dt*dt*swild(1:2)

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
              ghat1(i,1:2)=ghat1(i,1:2)-rampbc*dt*area(i)*grav/rho0*(ze(k,i)-ze(k-1,i))/2* &
     &(2*swild2(k,1:2)+dr_dxy(1:2,k,i)*(ze(k,i)-ze(k-1,i)))
            enddo !k
            ghat1(i,1:2)=ghat1(i,1:2)+rampbc*chigamma*dt*dt*area(i)*grav/rho0*swild2(kbe(i)+1,1:2)
          endif !iupwind_t
        endif !ibc==0

!       Non-hydrostatic press. gradient
        if(nonhydro==1) then
          swild(1:2)=0 !averaged bottom value
          do j=1,3 !side
            isd=js(i,j)
            bigfc1=0; bigfc2=0 !\int_{-h}^\eta f dz; (f = -\nabla qnon)
            do k=kbs(isd)+1,nvrt
              bigfc1=bigfc1-(zs(k,isd)-zs(k-1,isd))*(dqnon_dxy(1,k,isd)+dqnon_dxy(1,k-1,isd))/2
              bigfc2=bigfc2-(zs(k,isd)-zs(k-1,isd))*(dqnon_dxy(2,k,isd)+dqnon_dxy(2,k-1,isd))/2
            enddo !k
            ghat1(i,1)=ghat1(i,1)+dt*area(i)*bigfc1/3
            ghat1(i,2)=ghat1(i,2)+dt*area(i)*bigfc2/3
            swild(1:2)=swild(1:2)-dqnon_dxy(1:2,kbs(isd)+1,isd)/3
          enddo !side
          ghat1(i,1:2)=ghat1(i,1:2)-area(i)*chigamma*dt*dt*swild(1:2)
        endif !nonhydro==1


!       Debug
!        if(myrank==irank0) write(96,*)i,ielg(i),ghat1(i,1:2)     

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
          sparsem(i,j)=0
        enddo !j
        qel(i)=0

!	Area integrals I_{1,4}
        do j=1,nne(i)
          ie=ine(i,j)
          id=iself(i,j)

!	  I_1
          !n2=nm(ie,nx(id,1))
          !n3=nm(ie,nx(id,2))
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
              isd=js(ie,jj)
              if(idry_s(isd)==0) then
                etam=(eta2(isidenode(isd,1))+eta2(isidenode(isd,2)))/2
                htot=etam+dps(isd)
                if(hhat(isd)<=0) then
                  write(errmsg,*)'Impossible dry 52:',hhat(isd),iplg(isidenode(isd,1:2))
                  call parallel_abort(errmsg)
                endif
                del=hhat(isd)*hhat(isd)+(theta2*cori(isd)*dt*htot)**2 !>0
                hhatb=hhatb+htot*htot*hhat(isd)/del/3
                avg2=avg2+cori(isd)*(htot*dt)**3/del/3
              endif !idry_s
            enddo !jj
          else !3D
            hhatb=(hhat(js(ie,1))+hhat(js(ie,2))+hhat(js(ie,3)))/3
          endif 
          tmp0=area(ie)/6+grav*thetai**2*dt**2/4/area(ie)*hhatb*dot1
          tmpj=area(ie)/12+grav*thetai**2*dt**2/4/area(ie)*hhatb*dot2
          tmpj1=area(ie)/12+grav*thetai**2*dt**2/4/area(ie)*hhatb*dot3
          sparsem(i,0)=sparsem(i,0)+tmp0
          sparsem(i,j)=sparsem(i,j)+tmpj
          if(isbnd(1,i)==0.and.j==nne(i)) then
            indx=1
          else
            indx=j+1
          endif
          sparsem(i,indx)=sparsem(i,indx)+tmpj1

          if(lm2d) then !additional Coriolis part
            tmp=grav*theta2*thetai**2*avg2/2
            sparsem(i,j)=sparsem(i,j)+tmp
            sparsem(i,indx)=sparsem(i,indx)-tmp
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
          isd1=js(ie,1)
          isd2=js(ie,2)
          isd3=js(ie,3)
          if(ics==1) then
            dot1=dl(ie,id,1)*(bigu(isd1,1)+bigu(isd2,1)+bigu(isd3,1))/3+ &
     &dl(ie,id,2)*(bigu(isd1,2)+bigu(isd2,2)+bigu(isd3,2))/3
          else
            call project_hvec(bigu(isd1,1),bigu(isd1,2),sframe(:,:,isd1),eframe(:,:,ie),bigu1,bigv1)
            call project_hvec(bigu(isd2,1),bigu(isd2,2),sframe(:,:,isd2),eframe(:,:,ie),bigu2,bigv2)
            call project_hvec(bigu(isd3,1),bigu(isd3,2),sframe(:,:,isd3),eframe(:,:,ie),bigu3,bigv3)
            dot1=dl(ie,id,1)*(bigu1+bigu2+bigu3)/3+dl(ie,id,2)*(bigv1+bigv2+bigv3)/3
          endif !ics
          dot2=dl(ie,id,1)*ghat1(ie,1)+dl(ie,id,2)*ghat1(ie,2)
        
          qel(i)=qel(i)+(1-thetai)*dt*area(ie)*dot1+thetai*dt*dot2
          do l=1,3
            if(id==l) then
              fac=2
            else
              fac=1
            endif
            nd=nm(ie,l)
            if(imm==2) then
!Error: vel. not projected for ics=2
              call update_bdef(time,xctr(ie),yctr(ie),dep,swild)
              ubed=swild(1); vbed=swild(2); wbed=swild(3)
              dpdx=0; dpdy=0
              do m=1,3
                dpdx=dpdx+dp(nm(ie,m))*dl(ie,m,1)
                dpdy=dpdy+dp(nm(ie,m))*dl(ie,m,2)
              enddo !m   
              vnorm=ubed*dpdx+vbed*dpdy+wbed
              qel(i)=qel(i)+area(ie)/12*fac*(eta2(nd)+dt*vnorm)
            else
              qel(i)=qel(i)+area(ie)/12*fac*(eta2(nd)+bdef2(nd)-bdef1(nd))
            endif !imm
          enddo !l

!...      Impose Point Source volume at the surface layer; added by YC
#ifdef USE_ICM
         if(iWQPS==2) then
           do l=1,3
             if(id==l) then
               fac=2
             else
              fac=1
             endif
             nd=nm(ie,l)
!YC            raintmp=raintmp+area(ie)/12*fac*fluxprc(nd)/rho0*dt
!YC            qel(i)=qel(i)+area(ie)/12*fac*fluxprc(nd)/rho0*dt
             qel(i)=qel(i)-PSQ(ie)*dt/12*fac
           enddo !l
         endif
#endif USE_ICM
        enddo !j=1,nne(i)

!	bnd integrals I_{2,3,5,6}; they all vanish at land bnds 
!	I_2,6 are not needed if essential b.c. are enforced by elminating rows and columns
        if(isbnd(1,i)>0) then !open bnd node
!          ibnd=isbnd(1,i)
          do l=1,2 !two open bnd sides
            if(l==1) then
              ie=ine(i,1)
              id=iself(i,1)
              isd=js(ie,nx(id,2))
              nj=nm(ie,nx(id,1))
              ind=1
            else
              ie=ine(i,nne(i))
              id=iself(i,nne(i))
              isd=js(ie,nx(id,1))
              nj=nm(ie,nx(id,2))
              ind=nnp(i)
            endif

            nd=isidenode(isd,1)+isidenode(isd,2)-i
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
                  vn1=uth(isd,k)*sframe(1,1,isd)+vth(isd,k)*sframe(2,1,isd)
                  vn2=uth(isd,k+1)*sframe(1,1,isd)+vth(isd,k+1)*sframe(2,1,isd)
                else 
                  call project_hvec(uth(isd,k),vth(isd,k),pframe(:,:,i),sframe(:,:,isd),vn1,vtmp)
                  call project_hvec(uth(isd,k+1),vth(isd,k+1),pframe(:,:,i),sframe(:,:,isd),vn2,vtmp)
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
                sparsem(i,0)=sparsem(i,0)+thetai*dt*con0*2
                sparsem(i,ind)=sparsem(i,ind)+thetai*dt*con0
              endif !Flather 1
              qel(i)=qel(i)-thetai*dt*ri3
            endif

!	    I_5
            if(isbs(isd)>0.and.idry_s(isd)==0) then
              if(ics==1) then
                Unbar=bigu(isd,1)*sframe(1,1,isd)+bigu(isd,2)*sframe(2,1,isd)
              else
                Unbar=bigu(isd,1)
              endif !ics
              qel(i)=qel(i)-(1-thetai)*dt*distj(isd)*Unbar/2
            endif
          enddo !l=1,2 sides
        endif !isbnd: bnd node i
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
            nd=inp(i,j)
            if(nd<=np) then !nd resident
              in1=0
              do l=1,nnp(nd)
                if(inp(nd,l)==i) in1=l
              enddo !l
              if(in1==0) then
                write(errmsg,*)'Not resident:',iplg(i),iplg(nd)
                call parallel_abort(errmsg)
              endif
              if(abs(sparsem(i,j)-sparsem(nd,in1))>1.e-5) then
                write(errmsg,*)'Matrix not symmetric:',iplg(i),j,iplg(nd),sparsem(i,j),sparsem(nd,in1)
                call parallel_abort(errmsg)
              endif
              irank_s=myrank
            else !nd is ghost
              if(.not.associated(ipgl(iplg(nd))%next)) call parallel_abort('Wrong ghost')
              irank_s=ipgl(iplg(nd))%next%rank
            endif

!           Output
!           write(10,*)iplg(i),iplg(nd),irank_s,sparsem(i,j),sparsem(i,0)
          enddo !j
        enddo !i=1,np
      endif !.not.lm2d
!      close(10)

!     Save eta1
      eta1=eta2

!     Solve the wave equation for elevations at each element center in aug. subdomain
      call solve_jcg(it,moitn,mxitn,rtol,sparsem,eta2,qel,elbc,lelbc)

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
          ie=is(j,l)
          if(ie/=0) then
            itmp=0
            do m=1,3
             nd=nm(ie,m)
             if(eta2(nd)+dp(nd)<=h0) itmp=1
            enddo !m
            if(itmp==0) then !wet
              icount2=icount2+1
              do m=1,3
                tmpx=eta2(nm(ie,m))*dl(ie,m,1) !!eframe if ics=2
                tmpy=eta2(nm(ie,m))*dl(ie,m,2)
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
                nd=nm(ie,m)
                tmpx1=eta1(nd)*dl(ie,m,1) !eframe if ics=2
                tmpy1=eta1(nd)*dl(ie,m,2)
                tmpx2=pr(nd)*dl(ie,m,1)
                tmpy2=pr(nd)*dl(ie,m,2)
                tmpx3=etp(nd)*dl(ie,m,1)
                tmpy3=etp(nd)*dl(ie,m,2)
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
      kbs_e=0 !larger of the 2 element bottom indices
      do j=1,ns
        if(idry_s(j)==1) cycle

!       Wet
        do l=1,2 !element
          ie=is(j,l)
          if(ie/=0.and.idry_e(max(1,ie))==0.and.kbe(max(1,ie))>kbs_e(j)) kbs_e(j)=kbe(ie)
        enddo !l
        if(kbs_e(j)==0) then
          write(errmsg,*)'Cannot find the higher bottom:',islg(j),(ielg(is(j,l)),l=1,2)
          call parallel_abort(errmsg)
        endif
        if(kbs(j)>kbs_e(j)) then
          write(errmsg,*)'Side index > elemnt:',kbs(j),kbs_e(j)
          call parallel_abort(errmsg)
        endif
        if(lm2d.and.kbs_e(j)/=1) then
          write(errmsg,*)'2D bottom index wrong:',kbs(j),kbs_e(j),iplg(isidenode(j,1:2))
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

!...  Along each side
!     su2, sv2 in sframe if ics=2
      do j=1,nsa !augumented
        if(idry_s(j)==1) then
          do k=1,nvrt
            su2(k,j)=0
            sv2(k,j)=0
          enddo !k
          cycle
        endif

!	Wet sides
        node1=isidenode(j,1)
        node2=isidenode(j,2)
!       ll frame at side
        swild10(1:3,1:3)=(pframe(:,:,node1)+pframe(:,:,node2))/2

        if(lm2d) then !2D
!-------------------------------------------------------------------------------------
          !Warning: don't use eta2 which is updated
          htot=zs(nvrt,j)-zs(kbs(j),j) !(eta2(node1)+eta2(node2))/2+dps(j)
          if(hhat(j)<=0) then
            write(errmsg,*)'Impossible dry 55:',hhat(j),iplg(isidenode(j,1:2))
            call parallel_abort(errmsg)
          endif
          del=hhat(j)*hhat(j)+(theta2*cori(j)*dt*htot)**2 !delta > 0
          !\hat{Gamma}
          !rotate wind to sframe
          taux2=(tau(node1,1)+tau(node2,1))/2
          tauy2=(tau(node1,2)+tau(node2,2))/2
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
            hat_gam_x=hat_gam_x+dt*wwave_force(1,j,1)
            hat_gam_y=hat_gam_y+dt*wwave_force(1,j,2)
#endif USE_WWM
          !Add deta2 to \hat{Gamma}
          hat_gam_x=hat_gam_x-grav*thetai*dt*deta2_dx(j)
          hat_gam_y=hat_gam_y-grav*thetai*dt*deta2_dy(j)
          su2(1,j)=(hhat(j)*htot*hat_gam_x+theta2*cori(j)*dt*htot*htot*hat_gam_y)/del
          sv2(1,j)=(hhat(j)*htot*hat_gam_y-theta2*cori(j)*dt*htot*htot*hat_gam_x)/del
          su2(1,j)=max(-rmaxvel,min(rmaxvel,su2(1,j)))
          sv2(1,j)=max(-rmaxvel,min(rmaxvel,sv2(1,j)))

!-------------------------------------------------------------------------------------
        else !3D
!-------------------------------------------------------------------------------------
!       Define layer thickness & diffusivities
        do k=kbs(j)+1,nvrt
          dzz(k)=zs(k,j)-zs(k-1,j)
          if(dzz(k)<=0) call parallel_abort('MAIN: dzz=0 in momentum')
          !dfz(k)=(dfv(node1,k)+dfv(node2,k)+dfv(node1,k-1)+dfv(node2,k-1))/4
          dfz(k)=(ptbt(3,k,node1)+ptbt(3,k,node2)+ptbt(3,k-1,node1)+ptbt(3,k-1,node2))/4
        enddo !k

!	Coefficient matrix 
        ndim=nvrt-kbs(j)
        do k=kbs(j)+1,nvrt
          kin=k-kbs(j) !eq. #
          alow(kin)=0 
          cupp(kin)=0
          bdia(kin)=0
          if(k<nvrt) then
            tmp=dt*dfz(k+1)/dzz(k+1)
            cupp(kin)=cupp(kin)+dzz(k+1)/6-tmp
            bdia(kin)=bdia(kin)+dzz(k+1)/3+tmp
          endif

          if(k>kbs(j)+1) then
            tmp=dt*dfz(k)/dzz(k)
            alow(kin)=alow(kin)+dzz(k)/6-tmp
            bdia(kin)=bdia(kin)+dzz(k)/3+tmp
          else !b.c.
            bdia(kin)=bdia(kin)+dt*chi(j)
          endif
        enddo !k

!	RHS 
!	b.c. to be imposed at the end
        do k=kbs(j)+1,nvrt
          kin=k-kbs(j)
          rrhs(kin,1)=0
          rrhs(kin,2)=0
!	  Elevation gradient, atmo. pressure and tidal potential
          if(k<nvrt) then
            rrhs(kin,1)=rrhs(kin,1)-dzz(k+1)/2*dt*(grav*thetai*deta2_dx(j)+ &
                        grav*(1-thetai)*deta1_dx(j)+dpr_dx(j)/rho0-0.69*grav*detp_dx(j))
            rrhs(kin,2)=rrhs(kin,2)-dzz(k+1)/2*dt*(grav*thetai*deta2_dy(j)+ &
                        grav*(1-thetai)*deta1_dy(j)+dpr_dy(j)/rho0-0.69*grav*detp_dy(j))
          endif
          if(k>kbs(j)+1) then 
            rrhs(kin,1)=rrhs(kin,1)-dzz(k)/2*dt*(grav*thetai*deta2_dx(j)+ &
                        grav*(1-thetai)*deta1_dx(j)+dpr_dx(j)/rho0-0.69*grav*detp_dx(j))
            rrhs(kin,2)=rrhs(kin,2)-dzz(k)/2*dt*(grav*thetai*deta2_dy(j)+ &
                        grav*(1-thetai)*deta1_dy(j)+dpr_dy(j)/rho0-0.69*grav*detp_dy(j))
          endif

!	  Coriolis, advection, wind stress, and horizontal viscosity
          if(k<nvrt) then
            rrhs(kin,1)=rrhs(kin,1)+dzz(k+1)/6*(2*sdbt(1,k,j)+sdbt(1,k+1,j)+ &
     &dt*cori(j)*(2*sv2(k,j)+sv2(k+1,j))+dt*(2*d2uv(1,k,j)+d2uv(1,k+1,j)))
            rrhs(kin,2)=rrhs(kin,2)+dzz(k+1)/6*(2*sdbt(2,k,j)+sdbt(2,k+1,j)- &
     &dt*cori(j)*(2*su2(k,j)+su2(k+1,j))+dt*(2*d2uv(2,k,j)+d2uv(2,k+1,j)))
          else !k=nvrt
            taux2=(tau(node1,1)+tau(node2,1))/2
            tauy2=(tau(node1,2)+tau(node2,2))/2
            if(ics==2) then
              call project_hvec(taux2,tauy2,swild10(1:3,1:3),sframe(:,:,j),taux2s,tauy2s)
              taux2=taux2s
              tauy2=tauy2s
            endif !ics
            rrhs(kin,1)=rrhs(kin,1)+dt*taux2
            rrhs(kin,2)=rrhs(kin,2)+dt*tauy2
          endif

          if(k>kbs(j)+1) then
            rrhs(kin,1)=rrhs(kin,1)+dzz(k)/6*(2*sdbt(1,k,j)+sdbt(1,k-1,j)+ &
     &dt*cori(j)*(2*sv2(k,j)+sv2(k-1,j))+dt*(2*d2uv(1,k,j)+d2uv(1,k-1,j)))
            rrhs(kin,2)=rrhs(kin,2)+dzz(k)/6*(2*sdbt(2,k,j)+sdbt(2,k-1,j)- &
     &dt*cori(j)*(2*su2(k,j)+su2(k-1,j))+dt*(2*d2uv(2,k,j)+d2uv(2,k-1,j)))
          endif 

!         Baroclinic
          if(ibc==0) then
            if(k<nvrt) then
               rrhs(kin,1)=rrhs(kin,1)+dzz(k+1)/6*dt*(2*bcc(1,k,j)+bcc(1,k+1,j))
               rrhs(kin,2)=rrhs(kin,2)+dzz(k+1)/6*dt*(2*bcc(2,k,j)+bcc(2,k+1,j))
            endif
            if(k>kbs(j)+1) then
               rrhs(kin,1)=rrhs(kin,1)+dzz(k)/6*dt*(2*bcc(1,k,j)+bcc(1,k-1,j))
               rrhs(kin,2)=rrhs(kin,2)+dzz(k)/6*dt*(2*bcc(2,k,j)+bcc(2,k-1,j))
            endif
          endif !ibc==0

!         Non-hydrostatic
          if(nonhydro==1) then
            if(k<nvrt) then
              rrhs(kin,1:2)=rrhs(kin,1:2)-dzz(k+1)/6*dt*(2*dqnon_dxy(1:2,k,j)+dqnon_dxy(1:2,k+1,j))
            endif
            if(k>kbs(j)+1) then
              rrhs(kin,1:2)=rrhs(kin,1:2)-dzz(k)/6*dt*(2*dqnon_dxy(1:2,k,j)+dqnon_dxy(1:2,k-1,j))
            endif
          endif !nonhydro==1

!         Radiation stress
#ifdef  USE_WWM
            if(k<nvrt) rrhs(kin,1:2)=rrhs(kin,1:2)+dzz(k+1)/6*dt* &
     &(2*wwave_force(k,j,1:2)+wwave_force(k+1,j,1:2))
            if(k>kbs(j)+1) rrhs(kin,1:2)=rrhs(kin,1:2)+dzz(k)/6*dt* &
     &(2*wwave_force(k,j,1:2)+wwave_force(k-1,j,1:2))
#endif USE_WWM
        enddo !k=kbs(j)+1,nvrt

        call tridag(nvrt,100,ndim,2,alow,bdia,cupp,rrhs,soln,gam)
        do k=kbs(j)+1,nvrt
          kin=k-kbs(j)
!         Impose limits
          su2(k,j)=max(-rmaxvel,min(rmaxvel,soln(kin,1)))
          sv2(k,j)=max(-rmaxvel,min(rmaxvel,soln(kin,2)))
        enddo !k
!-------------------------------------------------------------------------------------
        endif !lm2d

        if(imm==2) then !no slip
          call update_bdef(time,xcj(j),ycj(j),dep,swild)
          su2(kbs(j),j)=swild(1)
          sv2(kbs(j),j)=swild(2)
        else if(.not.lm2d.and.Cd(j)==0) then
          su2(kbs(j),j)=su2(kbs(j)+1,j)
          sv2(kbs(j),j)=sv2(kbs(j)+1,j)
        else if(.not.lm2d) then !no slip bottom
          su2(kbs(j),j)=0
          sv2(kbs(j),j)=0
        endif

!       Extend
        do k=1,kbs_e(j)-1
          su2(k,j)=0 !su2(kbs(j),j)
          sv2(k,j)=0 !sv2(kbs(j),j)
        enddo !k

!       Impose uniformity for 2D
        if(lm2d) then
          su2(2,j)=su2(1,j)
          sv2(2,j)=sv2(1,j)
        endif

!	Impose b.c.
        do k=kbs_e(j),nvrt
          if(isbs(j)>0.and.ifltype(max(1,isbs(j)))/=0) then !open bnd side
            if(uth(j,k)<-98.or.vth(j,k)<-98) then
              write(errmsg,*)'Wrong vel. input:',uth(j,k),vth(j,k),node1,node2
              call parallel_abort(errmsg)
            endif
            !rotate uth, vth to sframe if ics=2; otherwise same
            uths=uth(j,k); vths=vth(j,k)
            if(ics==2) call project_hvec(uth(j,k),vth(j,k),swild10(1:3,1:3),sframe(:,:,j),uths,vths)

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
                vnorm=vnorm+uth(j,k)*sframe(1,1,j)+vth(j,k)*sframe(2,1,j)
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
                su2(k,j)=(1-vobc1(isbs(j)))*su2(k,j)+vobc1(isbs(j))*uths !uth(j,k)
                sv2(k,j)=(1-vobc1(isbs(j)))*sv2(k,j)+vobc1(isbs(j))*vths !vth(j,k)
              else !outgoing
                su2(k,j)=(1-vobc2(isbs(j)))*su2(k,j)+vobc2(isbs(j))*uths !uth(j,k)
                sv2(k,j)=(1-vobc2(isbs(j)))*sv2(k,j)+vobc2(isbs(j))*vths !vth(j,k)
              endif
            else !not Flather or 3D radiation
              su2(k,j)=uths !uth(j,k)
              sv2(k,j)=vths !vth(j,k)
            endif !Flather or not
          endif !open bnd

          if(isbs(j)==-1) then !land bnd
            if(islip==0) then !free slip
              if(ics==1) then
                vtan=su2(k,j)*sframe(1,2,j)+sv2(k,j)*sframe(2,2,j)
                su2(k,j)=vtan*sframe(1,2,j)
                sv2(k,j)=vtan*sframe(2,2,j)
              else !lat/lon
                su2(k,j)=0
              endif !ics
            else !no slip
              su2(k,j)=0
              sv2(k,j)=0
            endif
          endif
        enddo !k
      enddo !j=1,nsa

!...  Shapiro filter (used only if indvel<=0)
!     use bcc as temporary variable (sframe)
      if(indvel<=0) then
        bcc=0
        do i=1,ns !residents only
          if(is(i,2)==0.or.idry_s(i)==1) cycle

!         Internal wet sides
          do k=kbs(i)+1,nvrt
            suru=0
            surv=0
            do j=1,4
              id=isidenei2(i,j)
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
          if(is(j,2)==0.or.idry_s(j)==1) cycle

          do k=kbs(j)+1,nvrt
            su2(k,j)=bcc(1,k,j)
            sv2(k,j)=bcc(2,k,j)
          enddo !k

!         2D
          if(lm2d) then
            su2(1,j)=su2(nvrt,j)
            sv2(1,j)=sv2(nvrt,j)
          endif

          do k=1,kbs_e(j)-1
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
        if(myrank==0) write(16,*)'start non-hydrostatic calculation...'
!...    Solve for intermediate we
        hp_int(:,:,2)=we(:,:) !previous step temporarily saved as hp_int(:,:,2)
        we=0
        do i=1,nea
          if(idry_e(i)==1) cycle

!         Wet element
          n1=nm(i,1); n2=nm(i,2); n3=nm(i,3)
!         Define layer thickness & diffusivities
          do k=kbe(i)+1,nvrt
            dzz(k)=ze(k,i)-ze(k-1,i)
            if(dzz(k)<=0) call parallel_abort('MAIN: dzz=0 in wvel')
            !dfz(k)=(dfv(n1,k)+dfv(n2,k)+dfv(n3,k)+dfv(n1,k-1)+dfv(n2,k-1)+dfv(n3,k-1))/6
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
            nd=nm(i,j); isd=js(i,j)
            if(eta2(nd)+dp(nd)<=h0) ifl=1
            dhdx=dhdx+dp(nd)*dl(i,j,1)
            dhdy=dhdy+dp(nd)*dl(i,j,2)
            ubar1=ubar1+su2(kbs(isd),isd)/3
            vbar1=vbar1+sv2(kbs(isd),isd)/3
            swild(1)=swild(1)+eta1(nd)*dl(i,j,1) !deta1_dx
            swild(2)=swild(2)+eta1(nd)*dl(i,j,2) !deta1_dy
            swild(3)=swild(3)+eta2(nd)*dl(i,j,1) !deta2_dx
            swild(4)=swild(4)+eta2(nd)*dl(i,j,2) !deta2_dy
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
            rrhs(ndim,1)=(ubar2*(thetai*swild(3)+(1-thetai)*swild(1))+ &
     &vbar2*(thetai*swild(4)+(1-thetai)*swild(2))+(eta2_bar-eta1_bar)/dt- &
     &(1-thetai)*hp_int(nvrt,i,2))/thetai
          else
            rrhs(ndim,1)=0
          endif

!         Middle levels
          do k=kbe(i)+1,nvrt-1
            kin=k-kbe(i)+1
            qnon_e1=0; qnon_e2=0
            do j=1,3
              isd=js(i,j)
              nd=nm(i,j)
              qnon_e1=qnon_e1+qnon(k-1,nd)/3
              qnon_e2=qnon_e2+qnon(k+1,nd)/3
            enddo !j
            rrhs(kin,1)=dzz(k+1)/6*(2*webt(k,i)+webt(k+1,i))+dzz(k)/6*(2*webt(k,i)+webt(k-1,i))- &
     &dt/2*(qnon_e2-qnon_e1)
          enddo !k

          call tridag(nvrt,100,ndim,1,alow,bdia,cupp,rrhs,soln,gam)

!         Below bottom we=0 already assigned
          do k=kbe(i),nvrt
            kin=k-kbe(i)+1
            we(k,i)=soln(kin,1)
          enddo !k
        enddo !i=1,nea

        if(myrank==0) write(16,*)'done solving interim vertical vel...'

!...    Compute averaged divergence in each prism (temporarily saved as hp_int(nvrt,nea,1))
        hp_int(:,:,1)=0
        do i=1,nea
          if(idry_e(i)==1) cycle

!	  Wet elements with 3 wet nodes
!	  Compute upward normals and areas @ all levels
          n1=nm(i,1)
          n2=nm(i,2)
          n3=nm(i,3)
          !av_bdef1=(bdef1(n1)+bdef1(n2)+bdef1(n3))/3 !average bed deformation
          !av_bdef2=(bdef2(n1)+bdef2(n2)+bdef2(n3))/3
          if(kbe(i)==0) then
            write(errmsg,*)'Impossible 95 (2)'
            call parallel_abort(errmsg)
          endif
          do l=kbe(i),nvrt
            xcon=(ynd(n2)-ynd(n1))*(znl(l,n3)-znl(l,n1))-(ynd(n3)-ynd(n1))*(znl(l,n2)-znl(l,n1))
            ycon=(xnd(n3)-xnd(n1))*(znl(l,n2)-znl(l,n1))-(xnd(n2)-xnd(n1))*(znl(l,n3)-znl(l,n1))
            zcon=area(i)*2
            area_e(l)=sqrt(xcon**2+ycon**2+zcon**2)/2
            if(area_e(l)==0) then
              write(errmsg,*)'Zero area (3):',i,l
              call parallel_abort(errmsg)
            endif
            sne(l,1)=xcon/area_e(l)/2
            sne(l,2)=ycon/area_e(l)/2
            sne(l,3)=zcon/area_e(l)/2 !>0
          enddo !l

          do l=kbe(i)+1,nvrt
            sum1=0 !all horizontal fluxes
            ubar=0 !av. vel. at level l
            vbar=0
            ubar1=0 !av. vel. at level l-1
            vbar1=0
            do j=1,3
              jsj=js(i,j)
              vnor1=su2(l,jsj)*sframe(1,1,jsj)+sv2(l,jsj)*sframe(2,1,jsj)
              vnor2=su2(l-1,jsj)*sframe(1,1,jsj)+sv2(l-1,jsj)*sframe(2,1,jsj)
              if(l-1<kbs(jsj).or.kbs(jsj)==0) then
                write(errmsg,*)'Impossible 94 (2):',l,kbs(jsj),ielg(i)
                call parallel_abort(errmsg)
              endif
              sum1=sum1+ssign(i,j)*(zs(l,jsj)-zs(l-1,jsj))*distj(jsj)*(vnor1+vnor2)/2

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
                bflux=ubed*sne(l-1,1)+vbed*sne(l-1,2)+wbed*sne(l-1,3)
              else
                bflux=0
              endif
            else
              bflux=ubar1*sne(l-1,1)+vbar1*sne(l-1,2)+we(l-1,i)*sne(l-1,3) !inner normal
            endif
            top=ubar*sne(l,1)+vbar*sne(l,2)+we(l,i)*sne(l,3)

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
            nd=nm(i,j)
            dhdx=dhdx+hmod(nd)*dl(i,j,1)
            dhdy=dhdy+hmod(nd)*dl(i,j,2)
            deta_dx=deta_dx+eta1(nd)*dl(i,j,1) !use old levels         
            deta_dy=deta_dy+eta1(nd)*dl(i,j,2) !use old levels         
            if(hmod(nd)<hmin) hmin=hmod(nd)
          enddo !j

          do k=kz,nvrt !no Z
            kin=k-kz+1
            dzds_av=0 !average dz_ds
            do j=1,3
              nd=nm(i,j)
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
            ie=ine(i,j)
            if(idry_e(ie)==0.and.kbp_e(i)>kbe(ie)) kbp_e(i)=kbe(ie) 
          enddo !j
          if(kbp_e(i)==nvrt+1) call parallel_abort('MAIN: impossible 04')

          do k=kbp_e(i),nvrt
            icount=0 !# of wet elements
            do j=1,nne(i)
              ie=ine(i,j)
              if(idry_e(ie)==1) cycle

              !Wet element ie
              icount=icount+1
              id=iself(i,j)
              do l=0,1 !two prisms (ie,k+l)
                if(k+l>=kbe(ie)+1.and.k+l<=nvrt) then !prism exists
                  do m=1,3 !nodes
                    if(m==1) then
                      nd=i; ind=0 !index in the ball
                      ind2=id !index in the element
                    else
                      nd=nm(ie,nx(id,m-1)); ind=j+m-2
                      if(m==3.and.isbnd(1,i)==0.and.j==nne(i)) ind=1
                      ind2=nx(id,m-1)
                    endif

                    do mk=-1,0 !two levels k+l+mk; inside loop the index "l=(l1,l2)" (in the notes) is fixed
                      !\hat{I_0} & \hat{I_r}
                      dot1=dl(ie,id,1)*dl(ie,ind2,1)+dl(ie,id,2)*dl(ie,ind2,2) !dot product of gradient of shape function
                      if(k+l>=kz+1) then !S prism
                        dsigma=sigma(k+l-kz+1)-sigma(k+l-kz)
                        if(dsigma<=0) call parallel_abort('MAIN:dsig<=0')
                        dgam0=sign(1.,0.5-l)/dsigma !d{gamma} for \bar{l}_2
                        dgam1=sign(1.,0.5+mk)/dsigma !!d{gamma} for l_2
                        vol=area(ie)*dsigma

                        hat_i0=0
                        do mm=1,3 !nodes
                          nd_lam=nm(ie,mm)
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
     &dgam0*(dl(ie,ind2,1)*dsdx+dl(ie,ind2,2)*dsdy)*dzds/72*(1+kronecker(k+l+mk,k+l+mmk))*(1+kronecker(i,nd_lam))+&
     &dgam1*(dl(ie,id,1)*dsdx+dl(ie,id,2)*dsdy)*dzds/72*(1+kronecker(k,k+l+mmk))*(1+kronecker(nd,nd_lam))
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
                  nd=inp(i,j)
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
        call solve_jcg_qnon(it,moitn,mxitn,rtol,nvrt,mnei,np,npa,ihydro,qmatr,qhat,qir)
        if(myrank==0) write(16,*)'done JCG solver...'
        
!...    Update qnon - non-hydrostatic pressure at n+1/2
        qnon=qnon+qhat

!...    Compute gradient of qhat
        call hgrad_nodes(0,nvrt,npa,nsa,qhat,dqnon_dxy)

!...    Solve for final u,v
        su2(:,1:ns)=su2(:,1:ns)-dt*thetai*dqnon_dxy(1,:,1:ns)
        sv2(:,1:ns)=sv2(:,1:ns)-dt*thetai*dqnon_dxy(2,:,1:ns)
        do j=1,ns
          if(idry_s(j)==1) cycle
          do k=1,kbs_e(j)-1
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
            nd=nm(i,j); isd=js(i,j)
            if(eta2(nd)+dp(nd)<=h0) ifl=1
            dhdx=dhdx+dp(nd)*dl(i,j,1)
            dhdy=dhdy+dp(nd)*dl(i,j,2)
            ubar1=ubar1+su2(kbs(isd),isd)/3
            vbar1=vbar1+sv2(kbs(isd),isd)/3
            swild(1)=swild(1)+eta1(nd)*dl(i,j,1) !deta1_dx
            swild(2)=swild(2)+eta1(nd)*dl(i,j,2) !deta1_dy
            swild(3)=swild(3)+eta2(nd)*dl(i,j,1) !deta2_dx
            swild(4)=swild(4)+eta2(nd)*dl(i,j,2) !deta2_dy
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
              nd=nm(i,j)
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
!              nd=nm(i,j)
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
          eta2(i)=eta2(i)*(1-elev_nudge(i))
        enddo !i
      endif !inu_elev

      if(inu_uv==1) then
        do i=1,nsa
          uvnu=(uv_nudge(isidenode(i,1))+uv_nudge(isidenode(i,2)))/2
          su2(:,i)=su2(:,i)*(1-uvnu)
          sv2(:,i)=sv2(:,i)*(1-uvnu)
        enddo !i
      endif !inu_uv

!...  solve for vertical velocities using F.V.
!...  For hydrostatic model, this is the vertical vel; for non-hydrostatic
!...  model, this is only used in upwind transport
!     swild98 for storing rotated hvel at 3 sides
      allocate(swild98(2,3,nvrt),stat=istat)
      we_fv=0 !for dry and below bottom levels; in eframe if ics=2
      do i=1,nea
        if(idry_e(i)==1) cycle

!	Wet elements with 3 wet nodes
!	Compute upward normals and areas @ all levels
        n1=nm(i,1)
        n2=nm(i,2)
        n3=nm(i,3)
        av_bdef1=(bdef1(n1)+bdef1(n2)+bdef1(n3))/3 !average bed deformation
        av_bdef2=(bdef2(n1)+bdef2(n2)+bdef2(n3))/3
        if(kbe(i)==0) then
          write(errmsg,*)'Impossible 95'
          call parallel_abort(errmsg)
        endif
        do l=kbe(i),nvrt
          if(ics==1) then
            xcon=(ynd(n2)-ynd(n1))*(znl(l,n3)-znl(l,n1))-(ynd(n3)-ynd(n1))*(znl(l,n2)-znl(l,n1))
            ycon=(xnd(n3)-xnd(n1))*(znl(l,n2)-znl(l,n1))-(xnd(n2)-xnd(n1))*(znl(l,n3)-znl(l,n1))
            zcon=area(i)*2
          else !lat/lon
            !eframe
            call cross_product(xel(2,i)-xel(1,i),yel(2,i)-yel(1,i),znl(l,n2)-znl(l,n1), &
     &                         xel(3,i)-xel(1,i),yel(3,i)-yel(1,i),znl(l,n3)-znl(l,n1), &
     &                         xcon,ycon,zcon)
          endif !ics

          area_e(l)=sqrt(xcon**2+ycon**2+zcon**2)/2
          if(area_e(l)==0) then
            write(errmsg,*)'Zero area:',i,l
            call parallel_abort(errmsg)
          endif
          sne(l,1)=xcon/area_e(l)/2 !in eframe
          sne(l,2)=ycon/area_e(l)/2
          sne(l,3)=zcon/area_e(l)/2 !>0
        enddo !l

!       Rotate hvel. for 3 sides at all levels
        ubar=0; vbar=0 !average bottom hvel
        do m=1,3 !side
          isd=js(i,m)
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
        dhdx=dp(n1)*dl(i,1,1)+dp(n2)*dl(i,2,1)+dp(n3)*dl(i,3,1) !eframe
        dhdy=dp(n1)*dl(i,1,2)+dp(n2)*dl(i,2,2)+dp(n3)*dl(i,3,2)
!        ubar=(su2(kbs(js(i,1)),js(i,1))+su2(kbs(js(i,2)),js(i,2))+su2(kbs(js(i,3)),js(i,3)))/3
!        vbar=(sv2(kbs(js(i,1)),js(i,1))+sv2(kbs(js(i,2)),js(i,2))+sv2(kbs(js(i,3)),js(i,3)))/3
!       we_fv=0 unless Cd=0
!        if(nonhydro==1) then
!          we_fv(kbe(i),i)=(av_bdef2-av_bdef1)/dt-dhdx*ubar-dhdy*vbar
!        else
!          we_fv(kbe(i),i)=(av_bdef2-av_bdef1)/dt !-dhdx*ubar-dhdy*vbar
!        endif
        if(imm==2) then
          call update_bdef(time,xctr(i),yctr(i),dep,swild)
          ubed=swild(1); vbed=swild(2); wbed=swild(3)
          bflux0=ubed*sne(kbe(i),1)+vbed*sne(kbe(i),2)+wbed*sne(kbe(i),3) !normal bed vel.
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
            jsj=js(i,j)
            if(ics==1) then
              vnor1=su2(l,jsj)*sframe(1,1,jsj)+sv2(l,jsj)*sframe(2,1,jsj)
              vnor2=su2(l+1,jsj)*sframe(1,1,jsj)+sv2(l+1,jsj)*sframe(2,1,jsj)
            else !lat/lon
              vnor1=su2(l,jsj) !normal
              vnor2=su2(l+1,jsj)
            endif !ics
            if(l<kbs(jsj).or.kbs(jsj)==0) then
              write(errmsg,*)'Impossible 94:',l,kbs(jsj),ielg(i)
              call parallel_abort(errmsg)
            endif
            sum1=sum1+ssign(i,j)*(zs(l+1,jsj)-zs(l,jsj))*distj(jsj)*(vnor1+vnor2)/2

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
            bflux=ubar*sne(l,1)+vbar*sne(l,2)+we_fv(l,i)*sne(l,3)
          endif

          we_fv(l+1,i)=(-sum1-(ubar1*sne(l+1,1)+vbar1*sne(l+1,2))*area_e(l+1) + &
     &bflux*area_e(l))/sne(l+1,3)/area_e(l+1)

#ifdef USE_ICM
          if(iWQPS==2) then
            if(l==PSK(i)-1) we_fv(l+1,i)=we_fv(l+1,i)-(PSQ(i)*dt)/area(i)    !added by YC, need check
          endif
#endif USE_ICM

!         Debug
!          tmp1=sum1
!          tmp2=(ubar1*sne(l+1,1)+vbar1*sne(l+1,2)+we(l+1,i)*sne(l+1,3))*area_e(l+1)-bflux*area_e(l)
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
              nd=nm(i,j)
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
            node1=isidenode(isd0,1)
            node2=isidenode(isd0,2)
            depth=dps(isd0)
          endif
          do k=kbb+1,nvrt
            if(l==1) then !nodes
              dzz(k)=znl(k,nd0)-znl(k-1,nd0)
              !dfz(k)=(dfh(nd0,k)+dfh(nd0,k-1))/2
              dfz(k)=(ptbt(4,k,nd0)+ptbt(4,k-1,nd0))/2
            else !sides
              dzz(k)=zs(k,isd0)-zs(k-1,isd0)
              !dfz(k)=(dfh(node1,k)+dfh(node2,k)+dfh(node1,k-1)+dfh(node2,k-1))/4
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
            rrhs(kin,1)=0
            rrhs(kin,2)=0
            if(k<nvrt) then 
              if(l==1) then
                rrhs(kin,1)=rrhs(kin,1)+dzz(k+1)/2*ptbt(1,k,nd0)
                rrhs(kin,2)=rrhs(kin,2)+dzz(k+1)/2*ptbt(2,k,nd0)
              else
                rrhs(kin,1)=rrhs(kin,1)+dzz(k+1)/2*sdbt(3,k,isd0)
                rrhs(kin,2)=rrhs(kin,2)+dzz(k+1)/2*sdbt(4,k,isd0)
              endif
            else !surface fluxes
              if(ihconsv/=0) then !heat flux
                if(l==1) then
                  rrhs(kin,1)=rrhs(kin,1)+dt/rho0/shw*sflux(nd0)
                else
                  rrhs(kin,1)=rrhs(kin,1)+dt/rho0/shw*(sflux(node1)+sflux(node2))/2
                endif
              endif
              if(isconsv/=0) then !salt flux
                if(l==1) then
                  rrhs(kin,2)=rrhs(kin,2)+dt/rho0*snd(k,nd0)*(fluxevp(nd0)-fluxprc(nd0))
                else
                  rrhs(kin,2)=rrhs(kin,2)+dt/rho0*ssd(k,isd0)* &
     &(fluxevp(node1)+fluxevp(node2)-fluxprc(node1)-fluxprc(node2))/2
                endif
              endif
            endif

            if(k>kbb) then 
              if(l==1) then
                rrhs(kin,1)=rrhs(kin,1)+dzz(k)/2*ptbt(1,k,nd0)
                rrhs(kin,2)=rrhs(kin,2)+dzz(k)/2*ptbt(2,k,nd0)
              else
                rrhs(kin,1)=rrhs(kin,1)+dzz(k)/2*sdbt(3,k,isd0)
                rrhs(kin,2)=rrhs(kin,2)+dzz(k)/2*sdbt(4,k,isd0)
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
                rrhs(k,5)=srad(nd0) !F(); more to come
                itmp=iwater_type(nd0)
              else !sides
                zz1=max(-5.e2_rkind,zs(k,isd0)-zs(nvrt,isd0))
                rrhs(k,5)=(srad(node1)+srad(node2))/2
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
              rrhs(k,5)=rrhs(k,5)*(rr*exp(zz1/d_1)+(1-rr)*exp(zz1/d_2)) !F()
              
              ! Zero out bottom
              !if(k==kbb) rrhs(k,5)=0
            endif !solar
          enddo !k=kbb,nvrt

!         Compute solar
          if(ihconsv/=0) then
            rrhs(nvrt-kbb+1,1)=rrhs(nvrt-kbb+1,1)+dt/rho0/shw/2*(rrhs(nvrt,5)-rrhs(nvrt-1,5))
            rrhs(1,1)=rrhs(1,1)+dt/rho0/shw/2*(rrhs(kbb+1,5)-rrhs(kbb,5))
            do k=kbb+1,nvrt-1
              kin=k-kbb+1
              rrhs(kin,1)=rrhs(kin,1)+dt/rho0/shw/2*(rrhs(k+1,5)-rrhs(k-1,5))
            enddo !k
          endif

          call tridag(nvrt,100,ndim,2,alow,bdia,cupp,rrhs,soln,gam)

!         Impose no flux condition at bottom B.L. for slipless bottom
!          if(l==1.and.Cdp(nd0)/=0.or.l==2.and.Cd(isd0)/=0) soln(1,1:2)=soln(2,1:2)
          itmpf=0 !flag
          if(l==1) then; if(Cdp(nd0)/=0) itmpf=1; endif
          if(l==2) then; if(Cd(isd0)/=0) itmpf=1; endif
          if(itmpf==1) soln(1,1:2)=soln(2,1:2)

!         Correct overshoots for S,T
!         Debug
!          if(l==1.and.i==23) then
!            write(98,*)totalflux
!            do k=1,nvrt
!              write(98,*)k,soln(k,1),ptbt(1,k,nd0),rrhs(k,5),rrhs(k,4)*dt
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
              rrhs(k,3)=ptbt(1,k,nd0) !for debugging only
            else
              if(tmin>sdbt(3,k,isd0)) tmin=sdbt(3,k,isd0)
              if(tmax<sdbt(3,k,isd0)) tmax=sdbt(3,k,isd0)
              if(smin>sdbt(4,k,isd0)) smin=sdbt(4,k,isd0)
              if(smax<sdbt(4,k,isd0)) smax=sdbt(4,k,isd0)
              rrhs(k,3)=sdbt(3,k,isd0)
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
            swild2(k,1:2)=soln(kin,1:2)
            if(ihconsv/=0) swild2(k,1)=max(tempmin,min(tempmax,soln(kin,1)))
            if(isconsv/=0) swild2(k,2)=max(saltmin,min(saltmax,soln(kin,2)))
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
                  ie=ine(nd0,1)
                  id=iself(nd0,1)
                  isd3=js(ie,nx(id,2))
                else
                  ie=ine(nd0,nne(nd0))
                  id=iself(nd0,nne(nd0))
                  isd3=js(ie,nx(id,1))
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
                  swild2(k,1)=tobc(ibnd)*tth(ibnd,1,1)+(1-tobc(ibnd))*swild2(k,1)
                else if(itetype(ibnd)==3) then
                  swild2(k,1)=tobc(ibnd)*tem0(k,nd0)+(1-tobc(ibnd))*swild2(k,1)
                else if(itetype(ibnd)==4) then
                  swild2(k,1)=tobc(ibnd)*tth(ibnd,ind,k)+(1-tobc(ibnd))*swild2(k,1)
                !else if(itetype(ibnd)==-4) then
                  !if(ics==1) then
                  !  vnn=su2(k,isd)*sframe(1,1,isd)+sv2(k,isd)*sframe(2,1,isd)
                  !else
                  !  vnn=su2(k,isd)
                  !endif !ics
                  !if(vnn<0) swild2(k,1)=tobc(ibnd)*tth(ibnd,ind,k)+(1-tobc(ibnd))*swild2(k,1)
                !else if(itetype(ibnd)==-1) then
                !  if(ics==1) then
                !    vnn=su2(k,isd)*sframe(1,1,isd)+sv2(k,isd)*sframe(2,1,isd)
                !  else
                !    vnn=su2(k,isd)
                !  endif !ics
                !  if(vnn<0) swild2(k,1)=tobc(ibnd)*tem0(k,nd0)+(1-tobc(ibnd))*swild2(k,1)
                endif

                if(isatype(ibnd)==1.or.isatype(ibnd)==2) then
                  swild2(k,2)=sobc(ibnd)*sth(ibnd,1,1)+(1-sobc(ibnd))*swild2(k,2)
                else if(isatype(ibnd)==3) then
                  swild2(k,2)=sobc(ibnd)*sal0(k,nd0)+(1-sobc(ibnd))*swild2(k,2)
                else if(isatype(ibnd)==4) then
                  swild2(k,2)=sobc(ibnd)*sth(ibnd,ind,k)+(1-sobc(ibnd))*swild2(k,2)
!                else if(isatype(ibnd)==-4) then
!                  if(ics==1) then
!                    vnn=su2(k,isd)*sframe(1,1,isd)+sv2(k,isd)*sframe(2,1,isd)
!                  else
!                    vnn=su2(k,isd)
!                  endif !ics
!                  if(vnn<0) swild2(k,2)=sobc(ibnd)*sth(ibnd,ind,k)+(1-sobc(ibnd))*swild2(k,2)
!                else if(isatype(ibnd)==-1) then
!                  if(ics==1) then
!                    vnn=su2(k,isd)*sframe(1,1,isd)+sv2(k,isd)*sframe(2,1,isd)
!                  else
!                    vnn=su2(k,isd)
!                  endif !ics
!                  if(vnn<0) swild2(k,2)=sobc(ibnd)*sal0(k,nd0)+(1-sobc(ibnd))*swild2(k,2)
                endif
              enddo !k
            endif !isbnd>0
          else !sides
            if(isbs(isd0)>0) then
              ibnd=isbs(isd0)
              nwild(1:2)=0
              do ll=1,2 !nodes
                ndo=isidenode(isd0,ll)
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
                  swild2(k,1)=tobc(ibnd)*tth(ibnd,1,1)+(1-tobc(ibnd))*swild2(k,1)
                else if(itetype(ibnd)==3) then
                  swild2(k,1)=tobc(ibnd)*(tem0(k,node1)+tem0(k,node2))/2+(1-tobc(ibnd))*swild2(k,1)
                else if(itetype(ibnd)==4) then
                  swild2(k,1)=tobc(ibnd)*(tth(ibnd,in1,k)+tth(ibnd,in2,k))/2+(1-tobc(ibnd))*swild2(k,1)
!                else if(itetype(ibnd)==-4) then
!                  if(ics==1) then
!                    vnn=su2(k,isd0)*sframe(1,1,isd0)+sv2(k,isd0)*sframe(2,1,isd0)
!                  else
!                    vnn=su2(k,isd0)
!                  endif !ics
!                  if(vnn<0) swild2(k,1)=tobc(ibnd)*(tth(ibnd,in1,k)+tth(ibnd,in2,k))/2+(1-tobc(ibnd))*swild2(k,1)
!                else if(itetype(ibnd)==-1) then
!                  if(ics==1) then
!                    vnn=su2(k,isd0)*sframe(1,1,isd0)+sv2(k,isd0)*sframe(2,1,isd0)
!                  else
!                    vnn=su2(k,isd0)
!                  endif !ics
!                  if(vnn<0) swild2(k,1)=tobc(ibnd)*(tem0(k,node1)+tem0(k,node2))/2+(1-tobc(ibnd))*swild2(k,1)
                endif

                if(isatype(ibnd)==1.or.isatype(ibnd)==2) then
                  swild2(k,2)=sobc(ibnd)*sth(ibnd,1,1)+(1-sobc(ibnd))*swild2(k,2)
                else if(isatype(ibnd)==3) then
                  swild2(k,2)=sobc(ibnd)*(sal0(k,node1)+sal0(k,node2))/2+(1-sobc(ibnd))*swild2(k,2)
                else if(isatype(ibnd)==4) then
                  swild2(k,2)=sobc(ibnd)*(sth(ibnd,in1,k)+sth(ibnd,in2,k))/2+(1-sobc(ibnd))*swild2(k,2)
!                else if(isatype(ibnd)==-4) then
!                  if(ics==1) then
!                    vnn=su2(k,isd0)*sframe(1,1,isd0)+sv2(k,isd0)*sframe(2,1,isd0)
!                  else
!                    vnn=su2(k,isd0)
!                  endif !ics
!                  if(vnn<0) swild2(k,2)=sobc(ibnd)*(sth(ibnd,in1,k)+sth(ibnd,in2,k))/2+(1-sobc(ibnd))*swild2(k,2)
!                else if(isatype(ibnd)==-1) then
!                  if(ics==1) then
!                    vnn=su2(k,isd0)*sframe(1,1,isd0)+sv2(k,isd0)*sframe(2,1,isd0)
!                  else
!                    vnn=su2(k,isd0)
!                  endif !ics
!                  if(vnn<0) swild2(k,2)=sobc(ibnd)*(sal0(k,node1)+sal0(k,node2))/2+(1-sobc(ibnd))*swild2(k,2)
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
                tnu=t_nudge(nd0)*vnf
                snu=s_nudge(nd0)*vnf
                if(tnu<0.or.tnu>1.or.snu<0.or.snu>1.or.vnf<0.or.vnf>1) then
                  write(errmsg,*)'Nudging factor out of bound (1):',tnu,snu,vnf
                  call parallel_abort(errmsg)
                endif

                if(inu_st==1) then !to i.c.
                  swild2(k,1)=swild2(k,1)*(1-tnu)+tem0(k,nd0)*tnu
                  swild2(k,2)=swild2(k,2)*(1-snu)+sal0(k,nd0)*snu
                else if(inu_st==2) then
                  swild2(k,1)=swild2(k,1)*(1-tnu)+tnd_nu(nd0,k)*tnu
                  swild2(k,2)=swild2(k,2)*(1-snu)+snd_nu(nd0,k)*snu
                endif
              else !sides
                if(zs(k,isd0)>=-vnh1) then
                  vnf=vnf1 !vertical nudging factor
                else if(zs(k,isd0)>=-vnh2) then
                  vnf=vnf1+(vnf2-vnf1)*(zs(k,isd0)+vnh1)/(-vnh2+vnh1) 
                else
                  vnf=vnf2
                endif
                tnu=(t_nudge(node1)+t_nudge(node2))/2*vnf
                snu=(s_nudge(node1)+s_nudge(node2))/2*vnf
                if(tnu<0.or.tnu>1.or.snu<0.or.snu>1.or.vnf<0.or.vnf>1) then
                  write(errmsg,*)'Nudging factor out of bound (2):',tnu,snu,vnf
                  call parallel_abort(errmsg)
                endif

                if(inu_st==1) then !to i.c.
                  swild2(k,1)=swild2(k,1)*(1-tnu)+(tem0(k,node1)+tem0(k,node2))/2*tnu
                  swild2(k,2)=swild2(k,2)*(1-snu)+(sal0(k,node1)+sal0(k,node2))/2*snu
                else if(inu_st==2) then
                  swild2(k,1)=swild2(k,1)*(1-tnu)+(tnd_nu(node1,k)+tnd_nu(node2,k))/2*tnu
                  swild2(k,2)=swild2(k,2)*(1-snu)+(snd_nu(node1,k)+snd_nu(node2,k))/2*snu
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
          tsel(1,k,i)=sum(tnd(k,nm(i,1:3))+tnd(k-1,nm(i,1:3)))/6
          tsel(2,k,i)=sum(snd(k,nm(i,1:3))+snd(k-1,nm(i,1:3)))/6
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
      if(iupwind_t/=0.or.iupwind_s/=0) then
!       b.c. and body forces
        bdy_frc=0; flx_sf=0; flx_bt=0

!       Salt exchange
        if(isconsv/=0) then
          do i=1,nea
            if(idry_e(i)==1) cycle
       
            n1=nm(i,1)
            n2=nm(i,2)
            n3=nm(i,3)
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
            n1=nm(i,1)
            n2=nm(i,2)
            n3=nm(i,3)

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
     &ze(nvrt,i),(l,znl(l,nm(i,1:3)),l=kbe(i),nvrt)
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

!       VIMS surface temperature; added by YC
#ifdef USE_ICM
        if(myrank==0) write(16,*)'start ICM adjust. surface T..'
        if(iSun==2) then
          do i=1,nea
            if(idry_e(i)==1) cycle
            n1=nm(i,1)
            n2=nm(i,2)
            n3=nm(i,3)
            sflux_e=(surf_t(n1)+surf_t(n2)+surf_t(n3))/3
            tsel(1,nvrt,i)=sflux_e
          enddo !i
        if(myrank==0) write(16,*)'end ICM adjust. surface T..'
        endif !iSun=2
#endif USE_ICM

        up_tvd=(iupwind_t==2.or.iupwind_s==2)
        tr_el(1:2,:,:)=tsel(1:2,:,:)
        call do_transport_tvd(it,0,up_tvd,tvd_mid,flimiter,2,ifltype, &
     &itetype,isatype,itrtype,tobc,sobc,trobc,difnum_max_l,nvrt,npa,ptbt(4,:,:))
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


#ifdef USE_ICM
      if(myrank==0) write(16,*)'impose ICM point source S..'
      if(iWQPS==2) then 
!       Impose point source bc
        do i=1,nea
          if(idry_e(i)==1) cycle
          if(PSQ(i).ne.0.) then
            bigv=area(i)*(ze(PSK(i),i)-ze(PSK(i)-1,i)) !volume
!YC            bdy_frc(2,PSK(i),i)=area(i)*(tsel(2,PSK(i),i)*(PSQ(i))/area(i))/bigv
!!YC            tsel(2,PSK(i),i)=0.
            tsel(2,PSK(i),i)=(tsel(2,PSK(i),i)*(bigv+PSQ(i)*dt)+(10.01*(-PSQ(i))*dt))/bigv   !zhujz
            do k=PSK(i)-1,PSK(i)-1
              tsel(2,k,i)=tsel(2,PSK(i),i)
            enddo
            if(tsel(2,PSK(i),i).lt.0.) tsel(2,PSK(i),i)=0.
          endif
        enddo !i
      endif !iWQPS
      if(myrank==0) write(16,*)'end impose ICM point source S..'
#endif USE_ICM

!       Nudging
        do i=1,nea
          if(idry_e(i)==1) cycle

          n1=nm(i,1)
          n2=nm(i,2)
          n3=nm(i,3)

          if(inu_st/=0) then
            do k=kbe(i)+1,nvrt
              if(ze(k,i)>=-vnh1) then
                vnf=vnf1 !vertical nudging factor
              else if(ze(k,i)>=-vnh2) then
                vnf=vnf1+(vnf2-vnf1)*(ze(k,i)+vnh1)/(-vnh2+vnh1)
              else
                vnf=vnf2
              endif
              tnu=(t_nudge(n1)+t_nudge(n2)+t_nudge(n3))/3*vnf
              snu=(s_nudge(n1)+s_nudge(n2)+s_nudge(n3))/3*vnf
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
                tnd_nu_e=(tnd_nu(n1,k)+tnd_nu(n2,k)+tnd_nu(n3,k)+tnd_nu(n1,k-1)+tnd_nu(n2,k-1)+tnd_nu(n3,k-1))/6
                snd_nu_e=(snd_nu(n1,k)+snd_nu(n2,k)+snd_nu(n3,k)+snd_nu(n1,k-1)+snd_nu(n2,k-1)+snd_nu(n3,k-1))/6
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
              ie=ine(i,j)
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
              ie=is(i,j)
              if(ie/=0.and.idry_e(max(1,ie))==0) then
                ta=ta+area(ie)
                kin=max0(k,kbe(ie))
                tt1=tt1+hp_int(kin,ie,1)*area(ie)
                ss1=ss1+hp_int(kin,ie,2)*area(ie)
              endif
            enddo !j
            if(ta==0) then 
              write(errmsg,*)'Isolated wet side (9):',i,(is(i,j),j=1,2)
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
!       Save previous results in tr_el (global array) to prepare for do_transport_* and for ecosim
        tr_el(1:ntracers,:,:)=trel(1:ntracers,:,:)

        select case(flag_model)
          case(-1) !generic; modify arrays below
!           user-defined tracer part
!           define bdy_frc, flx_sf, flx_bt
!           bdy_frc(mntr,kbe(i)+1:nvrt,nea): body force at prism center Q_{i,k} (for all wet elements i)
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

          case(0) !for testing 
            bdy_frc = 0
            flx_bt = 0
            flx_sf = 0
          case(1) !Sediment
#ifdef USE_SED
              if(myrank==0) write(16,*) 'Entering sediment model...'

!LLP
#if defined BEDLOAD_VR
!----------------------------------------------------------------------
! Compute element depth averaged hvel for VRIJN bedload
!----------------------------------------------------------------------
                dave=0.d0
                do i=1,nea
                  if (idry_e(i)==1) cycle
                  cff1=0.d0
                  cff2=0.d0
                  n1=nm(i,1)
                  n2=nm(i,2)
                  n3=nm(i,3)
                  cff1=(dav(n1,1)+dav(n2,1)+dav(n3,1))/3
                  cff2=(dav(n1,2)+dav(n2,2)+dav(n3,2))/3
                  dave(i)=sqrt(cff1*cff1+cff2*cff2)
                enddo !i
#endif

              bdy_frc = 0.d0
              flx_bt = 0.d0
              flx_sf = 0.d0

              call sediment(it   &
#if defined BEDLOAD_MPM || defined BEDLOAD_VR && defined SED_MORPH
     &                      ,moitn,mxitn,rtol &
#endif
#if defined BEDLOAD_VR
     &                     ,dave             &
#endif
     &                      )

#endif USE_SED
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

            if(myrank==0) write(16,*) 'Calculating spectral irradiance...'
!'
            call spec_ir(Tair, Pair, Hair, cloud, Uwind, Vwind) !,SpecIr, avcos)

!...        Calculates sources and sinks terms of the ecological model
            if(myrank==0) write(16,*) 'Calculating ecological sources and sinks terms...'
!'
            bdy_frc = 0
            flx_bt = 0
            flx_sf = 0

            call ecosim(Uwind,Vwind)
            if(myrank==0) write(16,*) 'Done ecological sources and sinks terms...'
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
#endif USE_NAPZD
          case(5) !ICM
!            tr_el(1:ntracers,:,:)=trel(1:ntracers,:,:) !already init.
            bdy_frc = 0.d0                                               !added by YC
            flx_bt = 0.d0                                                !added by YC
            flx_sf = 0.d0
          case default
            call parallel_abort('Unknown tracer model (9)')
        end select !flag_model

        up_tvd=itr_met==2
        call do_transport_tvd(it,1,up_tvd,tvd_mid2,flimiter2,ntracers,ifltype, &
     &itetype,isatype,itrtype,tobc,sobc,trobc,difnum_max_l,nvrt,npa,ptbt(4,:,:))
        if(myrank==0) write(16,*)'done tracer transport...'

#ifdef USE_ICM
        if(myrank==0) write(16,*)'start ICM (5)..'
        if(amod(real(time),86400.)==0.) then
          call WQinput !(time)
          call wqm_out
        endif
        if(myrank==0) write(16,*) 'Calculating ecological sources and sinks terms...'     !added by YC
        call ecosystem(iths,it,rnday,WSRP,WSLP,WSPB1,WSPB2,WSPB3,turb,WRea,PSQ)           !added by YC
        if(myrank==0) write(16,*) 'Done ecological sources and sinks terms...'            !added by YC
#endif USE_ICM

        trel(1:ntracers,:,:)=tr_el(1:ntracers,:,:)
        if(difnum_max_l>difnum_max_l2) difnum_max_l2=difnum_max_l

#ifdef USE_ICM
      if(myrank==0) write(16,*)'adjust ICM pt source..'
      if(iWQPS==2) then
!       Impose point source bc added by YC
        do i=1,nea
          if(idry_e(i)==1) cycle
          if(PSQ(i).ne.0.) then
!YC	    if(it.eq.1) then
!YC	      TotV_old = 0.
!YC	      TotM_old = 0.
!YC              do k=kbe(i)+1,nvrt
!YC	        e_vol=(ze(k,i)-ze(k-1,i))*area(i)
!YC		TotM_old=  TotM_old + e_vol*trel(23,k,i)
!YC                TotV_old = TotV_old + e_vol
!YC              enddo
!YC              init_23_tr= TotM_old/TotV_old
!YC	    endif
!YC	    TotV_new = 0.
!YC	    TotM_new = 0.
!YC            do k=kbe(i)+1,nvrt
!YC	      e_vol=(ze(k,i)-ze(k-1,i))*area(i)
!YC	      TotM_new=  TotM_new + e_vol*trel(23,k,i)
!YC              TotV_new = TotV_new + e_vol
!YC            enddo
!YC            if (init_23_tr.le.WWPDO(i)/(86400*abs(PSQ(i)))) then
!YC	      do k=kbe(i)+1,nvrt
!YC                trel(23,k,i)=(TotM_old+(WWPDO(i)*dt/86400.))/(TotV_new)                  !tested by YC
!YC                if (trel(23,k,i).gt.WWPDO(i)/(86400*abs(PSQ(i)))) then
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
!!YC            trel(n,PSK(i),i)=total_loading/(86400*abs(PSQ(i)))
            trel(n,PSK(i),i)=(trel(n,PSK(i),i)*(bigv+PSQ(i)*dt)+(total_loading*dt/86400.))/bigv
!!YC            trel(n,PSK(i)-1:PSK(i)-4,i)=trel(n,PSK(i),i)
            do k=PSK(i)-1,PSK(i)-6
              trel(n,k,i)=trel(n,PSK(i),i)
            enddo
!YC            trel(n,nvrt,i)=WWPDO(i)/(86400*abs(PSQ(i)))
          enddo
!YC                endif
!YC              enddo
!YC            endif    
!YC            if (init_23_tr.gt.WWPDO(i)/(86400*abs(PSQ(i)))) then
!YC              do k=kbe(i)+1,nvrt
!YC                trel(23,nvrt,i)=WWPDO(i)/(86400*abs(PSQ(i)))
!YC                trel(23,k,i)=(TotM_old+(WWPDO(i)*dt/86400.))/(TotV_new)                  !tested by YC
!YC                if (trel(23,k,i).lt.WWPDO(i)/(86400*abs(PSQ(i)))) then
!YC                  trel(23,k,i)=WWPDO(i)/(86400*abs(PSQ(i)))
!YC                endif
!YC              enddo
!YC	    endif
!YC            if(myrank==0)write(995,'(24f18.4)')time/86400.,PSQ(i),WWPDO(i),trel(23,nvrt,i)
!YC 	    TotV_old=TotV_new
!YC         TotM_old=TotM_new
          endif
        enddo !i
      endif
      if(myrank==0) write(16,*)'done adjust ICM pt source..'
#endif USE_ICM

!       Nudging
        do i=1,nea
          if(idry_e(i)==1) cycle

          n1=nm(i,1)
          n2=nm(i,2)
          n3=nm(i,3)

          if(inu_tr/=0) then
            do k=kbe(i)+1,nvrt
              trnu=(tr_nudge(n1)+tr_nudge(n2)+tr_nudge(n3))/3
              if(trnu<0.or.trnu>1) then
                write(errmsg,*)'Nudging factor out of bound (2):',trnu
                call parallel_abort(errmsg)
              endif

              if(inu_tr==1) then !to i.c.
                trel(1:ntracers,k,i)=trel(1:ntracers,k,i)*(1-trnu)+trel0(1:ntracers,k,i)*trnu
              !else if(inu_st==2) then
              endif
            enddo !k
          endif !inu_st/=0

!         Extend
          do k=1,kbe(i)
            trel(1:ntracers,k,i)=trel(1:ntracers,kbe(i)+1,i)
          enddo !k
        enddo !i=1,nea

!Debug
!        write(12,*)'stage 1'

!       Convert to nodes and whole levels only for output; defined at resident nodes only
!       Use tr_el to temporarily store values at elements and whole levels
        if(mod(it,nspool)==0) then
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

            n1=nm(i,1)
            n2=nm(i,2)
            n3=nm(i,3)
!!...      difussivity of surface level (nvrt)
            av_df=(dfh(n1,nvrt-1)+dfh(n2,nvrt-1)+dfh(n3,nvrt-1))/3
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
!            if(myrank==0)write(16,*)'vol',vol,tr_tl(1,i)
            if(vol1==0) call parallel_abort('MAIN: sed. div. by 0 (2)')
            tr_el(1:ntracers,kbe(i),i)=(tr_tc(1:ntracers,i)-tr_tl(1:ntracers,i))/vol1
!            if(myrank==0)write(16,*)'vol1',vol1,tr_tc(1,i),tr_tl(1,i)
#else
            tr_el(1:ntracers,nvrt,i)=trel(1:ntracers,nvrt,i)
            tr_el(1:ntracers,kbe(i),i)=trel(1:ntracers,kbe(i)+1,i)
#endif USE_SED
          enddo !i=1,nea

          tr_nd=-99 !for dry nodes
          do i=1,np
            if(idry(i)==1) cycle

            do k=1,nvrt
              swild(1:ntracers)=0
              ta=0
              do j=1,nne(i)
                ie=ine(i,j)
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
                ie=ine(i,j)
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

#endif USE_NAPZD

!#ifdef INCLUDE_TIMING
!          cwtmp=mpi_wtime()
!#endif
!          call exchange_p3d_tr(tr_nd)
!#ifdef INCLUDE_TIMING
!          wtimer(9,2)=wtimer(9,2)+mpi_wtime()-cwtmp
!#endif
        endif !mod(it,nspool)==0

!...    Compute total mass (please move this part after levels are updated)
!...
!        swild(1:ntracers)=0
!        do i=1,ne
!          if(idry_e(i)==1) cycle
!
!          do k=kbe(i)+1,nvrt
!            vol=(ze(k,i)-ze(k-1,i))*area(i)
!            swild(1:ntracers)=swild(1:ntracers)+vol*trel(1:ntracers,k,i)
!          enddo !k
!        enddo !i=1,ne
!#ifdef INCLUDE_TIMING
!        cwtmp=mpi_wtime()
!#endif
!        call mpi_allreduce(swild,swild3,ntracers,rtype,MPI_SUM,comm,ierr)
!#ifdef INCLUDE_TIMING
!        wtimer(9,2)=wtimer(9,2)+mpi_wtime()-cwtmp
!#endif
!        if(myrank==0) write(92,*)time/86400,swild3(1:ntracers)
      endif !ntracers>0
!...  End of tracer transport

      if(myrank==0) write(16,*)'done solving transport equation'

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
        hmod(i)=min(dp(i),h_s)
      enddo !i

      do i=1,nsa
        n1=isidenode(i,1)
        n2=isidenode(i,2)
        dps(i)=(dp(n1)+dp(n2))/2
      enddo !i
      do i=1,nea
        dpe(i)=1.e10
        do j=1,3
          if(dpe(i)>dp(nm(i,j))) dpe(i)=dp(nm(i,j))
        enddo !j
      enddo !i=1,nea

!...  Recompute vgrid and calculate rewetted pts
      if(inunfl==0) then
        call levels0(iths,it)
      else
        call levels1(iths,it)
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
      endif !1==2

!...  Compute nodal vel. for output and next backtracking
      call nodalvel(ifltype)

!...  Density (using new level indices)
      prho=-99
      do i=1,npa
        if(idry(i)==1) cycle
        do k=1,nvrt
          prho(k,i)=eqstate(tnd(k,i),snd(k,i))
        enddo !k
      enddo !i

      if(iupwind_t/=0) then
        erho=-99
        do i=1,nea
          if(idry_e(i)==1) cycle
          do k=1,nvrt
            erho(k,i)=eqstate(tsel(1,k,i),tsel(2,k,i))
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
          dav(i,1)=dav(i,1)+(uu2(k+1,i)+uu2(k,i))/2*(znl(k+1,i)-znl(k,i))
          dav(i,2)=dav(i,2)+(vv2(k+1,i)+vv2(k,i))/2*(znl(k+1,i)-znl(k,i))
        enddo !k
        htot=eta2(i)+dp(i)
        if(htot<=h0) then
          write(errmsg,*)'Impossible 24'
          call parallel_abort(errmsg)
        endif
        dav(i,1)=dav(i,1)/htot
        dav(i,2)=dav(i,2)/htot

!       Max. dav (based on magnitude)
        dav_mag=sqrt(dav(i,1)**2+dav(i,2)**2)
        if(dav_mag>dav_maxmag(i)) then
          dav_maxmag(i)=dav_mag
          dav_max(i,1:2)=dav(i,1:2)
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

        n1=nm(i,1)
        n2=nm(i,2)
        n3=nm(i,3)
        etam=(eta2(n1)+eta2(n2)+eta2(n3))/3
        tpe=tpe+0.5*rho0*grav*area(i)*etam**2
        av_dep=etam+(dp(n1)+dp(n2)+dp(n3))/3
        tvol=tvol+area(i)*av_dep
!        do k=kbe(i),nvrt-1
!          ah=(znl(k+1,n1)+znl(k+1,n2)+znl(k+1,n3)-znl(k,n1)-znl(k,n2)-znl(k,n3))/3
!        enddo !k

        do j=1,3 !node or side
          nd=nm(i,j)
          do k=kbp(nd),nvrt-1
            tmass=tmass+area(i)*(prho(k,nd)+prho(k+1,nd))*(znl(k+1,nd)-znl(k,nd))/6
          enddo !k
          htot=eta2(nd)+dp(nd)
          if(htot<=h0) then
            write(errmsg,*)'Impossible dry (9):',ielg(i),j,iplg(nd),htot
            call parallel_abort(errmsg)
          endif

          isd=js(i,j)
          do k=kbs(isd),nvrt-1
            vmag1=su2(k,isd)**2+sv2(k,isd)**2
            vmag2=su2(k+1,isd)**2+sv2(k+1,isd)**2
            tkne=tkne+rho0*area(i)/6*(zs(k+1,isd)-zs(k,isd))*(vmag1+vmag2)/2
          enddo !k

!         enerf only correct for quasi-2D model
          enerf=enerf+dt*area(i)/3*rho0*Cdp(nd)*sqrt(dav(nd,1)**2+dav(nd,2)**2)**3

!         ener_ob
          isd=js(i,j)
          if(isbs(isd)>0) then !open bnd; no sharing between processes
            n1=isidenode(isd,1)
            n2=isidenode(isd,2)
            eta_m=(eta2(n1)+eta2(n2))/2
!Error: may not be accurate near poles
            vel_m1=(dav(n1,1)+dav(n2,1))/2 !both in ll frame
            vel_m2=(dav(n1,2)+dav(n2,2))/2
            ener_ob=ener_ob+rho0/2*sqrt(grav*dps(isd))*dt*(grav*eta_m**2+dps(isd)*(vel_m1**2+vel_m2**2))*distj(isd)
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

!     Fluxes
      open(32,file='fluxflag.gr3',status='old')
      read(32,*)
      read(32,*) itmp1,itmp2
      if(itmp1/=ne_global.or.itmp2/=np_global) &
     &call parallel_abort('Check fluxflag.gr3')
      do i=1,np_global
        read(32,*)j,xtmp,ytmp,tmp
        if(ipgl(i)%rank==myrank) nwild2(ipgl(i)%id)=tmp
      enddo
      close(32)
      do i=1,nea !must be aug.
        nwild(i)=0
        do j=1,3
          nd=nm(i,j)
          if(nwild(i)<nwild2(nd)) nwild(i)=nwild2(nd)
        enddo !j
      enddo !i

      tvol12=0 !total volume inside rgns 1 and 2
      fluxbnd=0 !total flux across natural open bnds (positive outward)
      fluxchan=0  !total flux across unnatural bnds (fluxchan1+fluxchan2)
      fluxchan1=0  !flux at bnd 1 (positive outward)
      fluxchan2=0  !flux at bnd 2 (positive outward)
 
      tot_s=0 !total salt inside rgns 1 and 2; [PSU*m^3]
      flux_s=0 !flux out of  rgns 1 and 2 (positive out) [PSU*m^3/s]
      do i=1,ne
        if(nwild(i)==0) cycle

        do j=1,3 !nodes or sides
          if(idry_e(i)==0) then
            nd=nm(i,j)
            do k=kbp(nd),nvrt-1
              ah=znl(k+1,nd)-znl(k,nd)
              tvol12=tvol12+area(i)*ah/3
              tot_s=tot_s+(snd(k+1,nd)+snd(k,nd))/2*ah*area(i)/3
            enddo !k
          endif !wet element

          isd=js(i,j)
          if(idry_s(isd)==0) then
            do k=kbs(isd),nvrt-1
              if(ic3(i,j)==0) then !bnd side
                if(ics==1) then
                  vnn=(su2(k+1,isd)+su2(k,isd))/2*sframe(1,1,isd)+(sv2(k+1,isd)+sv2(k,isd))/2*sframe(2,1,isd) 
                else
                  vnn=(su2(k+1,isd)+su2(k,isd))/2
                endif !ics
                ftmp=ssign(i,j)*distj(isd)*(zs(k+1,isd)-zs(k,isd))*vnn
                fluxbnd=fluxbnd+ftmp
              else if(ic3(i,j)>0.and.nwild(max(1,ic3(i,j)))==0) then !channel side
                if(ics==1) then
                  vnn=(su2(k+1,isd)+su2(k,isd))/2*sframe(1,1,isd)+(sv2(k+1,isd)+sv2(k,isd))/2*sframe(2,1,isd)
                else
                  vnn=(su2(k+1,isd)+su2(k,isd))/2
                endif !ics
                ftmp=ssign(i,j)*distj(isd)*(zs(k+1,isd)-zs(k,isd))*vnn
                fluxchan=fluxchan+ftmp
                if(nwild(i)==1) then
                  fluxchan1=fluxchan1+ftmp
                else !nwild(i)=2
                  fluxchan2=fluxchan2+ftmp
                endif

                flux_s=flux_s+ftmp*(ssd(k+1,isd)+ssd(k,isd))/2 !*dt=total salt lost
              endif
            enddo !k
          endif !wet side
        enddo !j=1,3
      enddo !i=1,ne

      buf3=0
      swild(1)=tvol12; swild(2)=tot_s; swild(3)=fluxbnd; swild(4)=fluxchan; swild(5)=fluxchan1; swild(6)=fluxchan2
#ifdef INCLUDE_TIMING
      cwtmp=mpi_wtime()
#endif
      call mpi_reduce(swild,buf3,6,rtype,MPI_SUM,0,comm,ierr)
#ifdef INCLUDE_TIMING
      wtimer(11,2)=wtimer(11,2)+mpi_wtime()-cwtmp
#endif
      if(myrank==0) then
        write(9,*)time/86400,buf3(5),-buf3(6),buf3(3),buf3(1),buf3(4),buf3(2),0.
        write(16,*)'done computing fluxes...'
      endif
      deallocate(buf3)
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
      ! Open existing processor specific file
      !fgb=ifile_char(1:ifile_len)//'_0000'; lfgb=len_trim(fgb)
      !write(fgb(lfgb-3:lfgb),'(i4.4)') myrank

      do j=1,noutput
        if(iof(j)==1.and.mod(it,nspool)==0) then
          a_4 = transfer(source=real(time),mold=a_4)
#ifdef USE_OPEN64
          !openMPI has trouble with no adv. write
          write(ichan(j)) a_4
#else
          write(ichan(j),"(a4)",advance="no") a_4
#endif
          a_4 = transfer(source=it,mold=a_4)
#ifdef USE_OPEN64
          write(ichan(j)) a_4
#else
          write(ichan(j),"(a4)",advance="no") a_4
#endif

          do i=1,np !residents only
            a_4 = transfer(source=real(eta2(i)),mold=a_4)
#ifdef USE_OPEN64
            write(ichan(j)) a_4
#else
            write(ichan(j),"(a4)",advance="no") a_4
#endif
          enddo !i

          do i=1,np !residents only
            if(j<=12) then
              floatout=0 !in case some are not defined
              if(j.eq.1) then
                if(idry(i)==1) then
                  floatout=-9999
                else
                  floatout=eta2(i)
                endif
              else if(j.eq.2) then !.and.ihconsv.ne.0) then
                floatout=pr(i)
              else if(j.eq.3.and.ihconsv.ne.0) then
                floatout=airt1(i)
              else if(j.eq.4.and.ihconsv.ne.0) then
                floatout=shum1(i)
              else if(j.eq.5.and.ihconsv.ne.0) then
                floatout=srad(i)
              else if(j.eq.6.and.ihconsv.ne.0) then
                floatout=fluxsu(i)
              else if(j.eq.7.and.ihconsv.ne.0) then
                floatout=fluxlu(i)
              else if(j.eq.8.and.ihconsv.ne.0) then
                floatout=hradu(i)
              else if(j.eq.9.and.ihconsv.ne.0) then
                floatout=hradd(i)
              else if(j.eq.10.and.ihconsv.ne.0) then
                floatout=sflux(i)
              else if(j.eq.11.and.isconsv.ne.0) then
                floatout=fluxevp(i)
              else if(j.eq.12.and.isconsv.ne.0) then
                floatout=fluxprc(i)
              endif

              a_4 = transfer(source=floatout,mold=a_4)
#ifdef USE_OPEN64
              write(ichan(j)) a_4
#else
              write(ichan(j),"(a4)",advance="no") a_4
#endif
            else if(j<=15) then
              if(j==13) then
                if(nws==0) then
                  floatout=0
                  floatout2=0
                else
                  floatout=windx(i) !in ll frame if ics=2
                  floatout2=windy(i)
                endif
              else if(j==14) then
                floatout=tau(i,1) !in ll frame if ics=2
                floatout2=tau(i,2)
              else !j=15
                floatout=dav(i,1)
                floatout2=dav(i,2)
              endif
              a_4 = transfer(source=floatout,mold=a_4)
#ifdef USE_OPEN64
              write(ichan(j)) a_4
#else
              write(ichan(j),"(a4)",advance="no") a_4
#endif
              a_4 = transfer(source=floatout2,mold=a_4)
#ifdef USE_OPEN64
              write(ichan(j)) a_4
#else
              write(ichan(j),"(a4)",advance="no") a_4
#endif
            else if(j<26) then
              do k=max0(1,kbp00(i)),nvrt
                floatout=0 !for some undefined variables
                if(j.eq.16) then
                  floatout=ww2(k,i)
                else
                  if(j.eq.17) then
                    if(idry(i)==1) then
                      floatout=-99
                    else
                      floatout=tnd(k,i)
                    endif
                  else if(j.eq.18) then
                    if(idry(i)==1) then
                      floatout=-99
                    else
                      floatout=snd(k,i)
                    endif
                  else if(j.eq.19) then
                    floatout=prho(k,i)
                  else if(j.eq.20) then
                    floatout=dfh(i,k)
                  else if(j.eq.21) then
                    floatout=dfv(i,k)
                  else if(j.eq.22) then
                    floatout=q2(i,k)
                  else if(j.eq.23) then
                    floatout=xl(i,k)
                  else if(j==24) then
                    if(idry(i)==1) then
                      floatout=0
                    else
                      floatout=znl(max0(k,kbp(i)),i)
                    endif
                  else if(j==25) then
                    floatout=qnon(k,i)
                  endif
                endif

                a_4 = transfer(source=floatout,mold=a_4)
#ifdef USE_OPEN64
                write(ichan(j)) a_4
#else
                write(ichan(j),"(a4)",advance="no") a_4
#endif
              enddo !k
            else if(j==26) then
              do k=max0(1,kbp00(i)),nvrt
                a_4 = transfer(source=real(uu2(k,i)),mold=a_4)
#ifdef USE_OPEN64
                write(ichan(j)) a_4
#else
                write(ichan(j),"(a4)",advance="no") a_4
#endif
               a_4 = transfer(source=real(vv2(k,i)),mold=a_4)
#ifdef USE_OPEN64
                write(ichan(j)) a_4
#else
                write(ichan(j),"(a4)",advance="no") a_4
#endif
              enddo !k

            else if(j<=26+ntracers) then !tracers; implies ntracers>0
              do k=max0(1,kbp00(i)),nvrt
                floatout=tr_nd(j-26,k,i) !warning: '26' may need to change
                a_4 = transfer(source=floatout,mold=a_4)
#ifdef USE_OPEN64
                write(ichan(j)) a_4
#else
                write(ichan(j),"(a4)",advance="no") a_4
#endif
              enddo !k
            else !optional modules; MUST BE IN THE SAME ORDER AS BEFORE
#ifdef USE_SED
              if(flag_model/=1) call parallel_abort('MAIN: strange output (2)')
!'
              if(j<=indx_out(1,2)) then
                if(j==indx_out(1,1)) then
                  floatout=dp(i)
                else if(j<=indx_out(1,1)+ntracers) then
                  floatout=bedldu(i,j-indx_out(1,1))
                else
                  floatout=bedldv(i,j-indx_out(1,1)-ntracers)
                endif
                a_4 = transfer(source=floatout,mold=a_4)
#ifdef USE_OPEN64
                write(ichan(j)) a_4
#else
                write(ichan(j),"(a4)",advance="no") a_4
#endif
              endif !scope of SED model
#endif USE_SED

#ifdef USE_NAPZD
              if(j<=indx_out(2,2)) then
                do k=max0(1,kbp00(i)),nvrt
                  if(j==indx_out(2,1)) then
                    floatout=Bio_bdefp(k,i)
                  else !total N
                    floatout=sum(tr_nd(1:4,k,i))
                  endif
                  a_4 = transfer(source=floatout,mold=a_4)
#ifdef USE_OPEN64
                  write(ichan(j)) a_4
#else
                  write(ichan(j),"(a4)",advance="no") a_4
#endif
                enddo !k
              endif !scope of NAPZD
#endif USE_NAPZD

#ifdef USE_WWM
              if(j<=indx_out(3,2)) then
                if(j<=indx_out(3,2)-2) then !scalar
                  itmp=j-indx_out(3,1)+1
                  if(itmp>23) call parallel_abort('MAIN: wwm_out over')
                  floatout=out_wwm(i,indx_wwm_out(itmp))
                else !vectors
                  if (j==indx_out(3,2)-1) then
                    floatout=out_wwm(i,7); floatout2=out_wwm(i,6);
                  else if (j==indx_out(3,2)) then
                    floatout=out_wwm(i,26); floatout2=out_wwm(i,27);
                  endif
                endif !j

                a_4 = transfer(source=floatout,mold=a_4)
#ifdef USE_OPEN64
                write(ichan(j)) a_4
#else
                write(ichan(j),"(a4)",advance="no") a_4
#endif

                if(j>indx_out(3,2)-2) then !vectors
                  a_4 = transfer(source=floatout2,mold=a_4)
#ifdef USE_OPEN64
                  write(ichan(j)) a_4
#else
                  write(ichan(j),"(a4)",advance="no") a_4
#endif
                endif !vectors
              endif !scope of WWM; j<=indx_out(3,2)
#endif USE_WWM
            endif !j
          enddo !i=1,np

          if(myrank==0) write(16,'(a48)')'done outputting '//variable_nm(j)
        endif !iof(j).eq.1.and.mod(it,nspool).eq.0
      enddo !j=1,noutput

!...  Test output 
! TODO

!     Open new global output files and write header data
      !if(it==ifile*ihfskip) then !.and.it/=ntime) then
      if(mod(it,ihfskip)==0) then
        ifile=ifile+1                   !output file #
        write(ifile_char,'(i12)') ifile !convert ifile to a string
        ifile_char=adjustl(ifile_char)  !place blanks at end
        ifile_len=len_trim(ifile_char)  !length without trailing blanks
!        call write_header1
        fgb=ifile_char(1:ifile_len)//'_0000'; lfgb=len_trim(fgb);
        write(fgb(lfgb-3:lfgb),'(i4.4)') myrank
        do i=1,noutput
          ichan(i)=100+i !output channel #
!          if(i>=13.and.i<=15.or.i==25) then
!            ivs=2
!          else
!            ivs=1
!          endif
!          if(i<=15) then
!            i23d=2 !2 or 3D
!          else
!            i23d=3
!          endif

          if(iof(i)==1) then
#ifdef USE_OPEN64
            open(ichan(i),file='outputs/'//(fgb(1:lfgb)//'_'//outfile(i)),status='replace', form='BINARY')
#else
            open(ichan(i),file='outputs/'//(fgb(1:lfgb)//'_'//outfile(i)),status='replace')
#endif
          endif
        enddo !i
      endif !it==ifile*ihfskip

!...  Station outputs
      if(iout_sta/=0) then
        do j=1,nvar_sta
          if(iof_sta(j)==0.or.mod(it,nspool_sta)/=0) cycle

          do i=1,nout_sta
            ie=iep_sta(i)
            if(ie==0) then !no parent
              iep_flag(i)=0 !for comm. later
              sta_out(i,j)=0
            else !has parent
              iep_flag(i)=1
              sta_out(i,j)=0 !initialize
              select case(j)
                case(1) !elev.
                  swild2(1,1:3)=eta2(nm(ie,1:3))
                case(2) !air pressure
                  swild2(1,1:3)=pr(nm(ie,1:3))
                case(3) !wind x
                  swild2(1,1:3)=windx(nm(ie,1:3))
                case(4) !wind y
                  swild2(1,1:3)=windy(nm(ie,1:3))
                case(5) !T
                  swild2(1:nvrt,1:3)=tnd(1:nvrt,nm(ie,1:3))
                case(6) !S
                  swild2(1:nvrt,1:3)=snd(1:nvrt,nm(ie,1:3))
                case(7) !u
!Error: may not be accurate near poles as pframe changes rapidly there
                  swild2(1:nvrt,1:3)=uu2(1:nvrt,nm(ie,1:3))
                case(8) !v
                  swild2(1:nvrt,1:3)=vv2(1:nvrt,nm(ie,1:3))
                case(9) !w
                  swild2(1:nvrt,1:3)=ww2(1:nvrt,nm(ie,1:3))
                case default
                  call parallel_abort('MAIN: unknown sta. output')
              end select

              if(j<=4) then !2D var.
                sta_out(i,j)=sum(arco_sta(i,1:3)*swild2(1,1:3))
              else !3D var.
                if(idry_e(ie)==1) then !dry
                  sta_out(i,j)=-999.
                else !wet
                  do m=1,3 !wet nodes
                    nd=nm(ie,m)
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

        if(myrank==0) then
!          write(290,*)nwild2(1:nout_sta)
          do i=1,nvar_sta
            if(iof_sta(i)==0.or.mod(it,nspool_sta)/=0) cycle
            do j=1,nout_sta
              if(nwild2(j)==0) then
                sta_out_gb(j,i)=-9999
              else
                sta_out_gb(j,i)=sta_out_gb(j,i)/nwild2(j)
              endif
            enddo !j
            write(250+i,'(e14.6,6000(1x,e14.6))')time,sta_out_gb(:,i)
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
#endif USE_HA

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
        write(it_char,'(i72)')it
        it_char=adjustl(it_char)
        lit=len_trim(it_char)
        it_char=it_char(1:lit)//'_0000'; lit=len_trim(it_char)
        write(it_char(lit-3:lit),'(i4.4)') myrank
        ihot_len=nbyte*(4+((6+4*ntracers)*nvrt+1)*ne+(8*nvrt+1)*ns+(3+22*nvrt)*np)
        open(36,file='outputs/'//it_char(1:lit)//'_hotstart', &
             access='direct',recl=ihot_len,status='replace') 

        write(36,rec=1)dble(time),it,ifile,(idry_e(i),(dble(we(j,i)),dble(tsel(1:2,j,i)), &
     &(dble(trel0(l,j,i)),dble(trel(l,j,i)),l=1,ntracers),j=1,nvrt),i=1,ne), &
     &(idry_s(i),(dble(su2(j,i)),dble(sv2(j,i)),dble(tsd(j,i)),dble(ssd(j,i)),j=1,nvrt),i=1,ns), &
     &(dble(eta2(i)),idry(i),(dble(tnd(j,i)),dble(snd(j,i)),dble(tem0(j,i)),dble(sal0(j,i)), &
     &dble(q2(i,j)),dble(xl(i,j)),dble(dfv(i,j)),dble(dfh(i,j)),dble(dfq1(i,j)),dble(dfq2(i,j)), &
     &dble(qnon(j,i)),j=1,nvrt),i=1,np)
        close(36)

#ifdef USE_HA
!... 
!...    IF APPROPRIATE ADD HARMONIC ANALYSIS INFORMATION TO HOT START FILE
!...    Adapted from ADCIRC
        open(36,file='outputs/'//it_char(1:lit)//'_hotstart', &
             access='direct',recl=8,status='old') 
        IHOTSTP = 3+((1+(3+2*ntracers)*nvrt)*ne)+((1+4*nvrt)*ns)+((2+9*nvrt)*np)
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
#endif USE_HA

        if(myrank==0) write(16,*) 'hot start written',it,time,ifile
      endif !nhot

#ifdef INCLUDE_TIMING
! End hotstart output section
      wtmp2=mpi_wtime()
      wtimer(13,1)=wtimer(13,1)+wtmp2-wtmp1
#endif

      if(myrank==0) write(16,'(a,i12,a,f20.6)') 'TIME STEP= ',it,';  TIME= ',time
!'
      call parallel_barrier !synchronize before starting next time step


!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! End Time Stepping
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
      enddo !it=iths+1,ntime

!...  Output max. diffusion number reached
      if(difnum_max_l2>0.5) write(12,*)'MAIN: max. diffusion # exceeds 0.5:',difnum_max_l2
!'

!...  Output max. # of iterations for all ranks for WBL (Grant-Madsen formulation)
#ifdef USE_WWM
        if(iwbl==1) then
           call mpi_reduce(iwbl_itmax,iwbl_itmax_gb,1,itype,MPI_MAX,0,comm,ierr)
           if(myrank==0) write(16,*)'Max. iteration for Grant-Madsen = ',iwbl_itmax_gb
        endif !iwbl
#endif USE_WWM


!...  Output max. elevations & dahv
      fdb='maxelev_0000'
      lfdb=len_trim(fdb)
      write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
      open(10,file='outputs/'//fdb,status='replace')
      write(10,*)np,nproc
      do i=1,np
        write(10,*)iplg(i),real(xnd(i)),real(ynd(i)),real(elevmax(i))
      enddo !i
      close(10)

      fdb='maxdahv_0000'
      lfdb=len_trim(fdb)
      write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
      open(10,file='outputs/'//fdb,status='replace')
      write(10,*)np,nproc
      do i=1,np
        write(10,*)iplg(i),real(xnd(i)),real(ynd(i)),real(dav_maxmag(i)),real(dav_max(i,1:2))
      enddo !i
      close(10)

#ifdef USE_HA
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! SOLVE THE HARMONIC ANALYSIS PROBLEM (Adapted from ADCIRC)
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!     Use 12 to write from each process
      write(12,*)'HA',iharind,it,ITHAS,CHARMV
      IF ((iharind.EQ.1).AND.(it.GT.ITHAS)) THEN
!...Compute means and variances for checking the harmonic analysis results
!...Accumulate mean and variance at each node.
        if (CHARMV) then
          IF (FMV.NE.0.) THEN
            DO I=1,np
              ELAV(I)   = ELAV(I)/NTSTEPS
              XVELAV(I) = XVELAV(I)/NTSTEPS
              YVELAV(I) = YVELAV(I)/NTSTEPS
              ELVA(I)   = ELVA(I)/NTSTEPS   - ELAV(I)*ELAV(I)
              XVELVA(I) = XVELVA(I)/NTSTEPS - XVELAV(I)*XVELAV(I)
              YVELVA(I) = YVELVA(I)/NTSTEPS - YVELAV(I)*YVELAV(I)
            END DO
            TIMEBEG=ITMV*dt
            write(it_char(1:4),'(i4.4)') myrank
            open(55,file='outputs/harme.55'//it_char(1:4))
            WRITE(55,*) np
          ENDIF
        endif

!......Fill out and decompose the LHS harmonic analaysis matrix

        CALL FULSOL(0)

!......Solve the harmonic analysis problem and write the output

        write(12,*)'myrank=',myrank,TIMEBEG,ITMV,dt
        IF(NHAGE.EQ.1) CALL LSQSOLEG(np,ELAV,ELVA,CHARMV,myrank,TIMEBEG,DT,FMV,NTSTEPS,ITMV)
        IF(NHAGV.EQ.1) CALL LSQSOLVG(np,XVELAV,YVELAV,XVELVA,YVELVA,CHARMV,myrank,TIMEBEG,DT,FMV,NTSTEPS,ITMV)
      ENDIF
#endif USE_HA

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Finalize
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

#ifdef INCLUDE_TIMING
      wtmp2=mpi_wtime()
      wtimer(2,1)=wtmp2-wtimer(2,1) !time-stepping section
      wtimer(0,1)=wtmp2-wtimer(0,1) !total
!     Report timing
      call report_timers
#endif

#ifdef USE_SWAN
      first_call=.false.
!     Gather outputs for coupler
      
      end subroutine elfe
#else
      call date_and_time(date,timestamp)
      if(myrank==0) write(16,'(/4a)') 'Run completed successfully at ',date,', ',timestamp

! Finalize parallel environment
      call parallel_finalize
      end program elfe
#endif
!===============================================================================
!===============================================================================
! END SELFE MAIN PROGRAM
!===============================================================================
!===============================================================================
