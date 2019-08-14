!===============================================================================
! Read in binary outputs (rank-specific) from parallel code and combine them into
! one global output in v5.0 format or netcdf format. Works for partial outputs.
! Global-local mappings are read in from separate files.
! Run this program inside the directory outputs/, where some of the input files below
! can be found.

! Inputs:
!        rank-specific binary files (from SELFE outputs/); 
!        local_to_global_* (from SELFE outputs/);
!        combine_output.in (1st line: elev.61 etc; 2nd line: start and end file indices;
!                          3rd line: inc (1 for netcdf option)); 
!        additional inputs for non-standard outputs (hvel.67 etc): 
!             sidecenters.gr3, centers.gr3 (in the same dir as hgrid.gr3)
! Output: combined binary file (for nc file: e.g. *_salt.nc; netcdf not working for 
!         non-standard outputs!).
!
!  Compile on tsunami:
!  ifort -Bstatic -O2 -assume byterecl -o combine_output7 combine_output7.f90 selfe_geometry.f90 -Vaxlib -I/share/apps/netcdf/include/ -L/share/apps/netcdf/lib/ -lnetcdf

!  On TYPHOON:
!  pgf90 -O2 -mcmodel=medium  -Bstatic -o combine_output7 combine_output7.f90 ~yinglong/SELFE/svn/trunk/src/Utility/UtilLib/selfe_geometry.f90 -I/usr/local/amd64a/opteron/pgi/netcdf-3.6.2/include -L/usr/local/amd64a/opteron/pgi/netcdf-3.6.2/lib -lnetcdf

!  History: 
!           2013-12-09 drf: NF_FLOATs, add SZ metadata, netcdf error checking towards ugrid v0.9
!              Add sigma attributes and sigma_theta_f, sigma_theta_b, sigma_h_c' variables for CF 
!              Add mesh variable for ugrid compliance
!              Parse binary file time (from bctides.in header) into udunits2 compatible units  
!   
!           changed non-standard outputs which can be viz'ed with vis6 now (vertical
!           coordinates not dependable).

! TODO: mark dry elev. according to SURA standard
!===============================================================================


program read_iwrite1
!-------------------------------------------------------------------------------
!  use typeSizes
!  use netcdf

  implicit real(4)(a-h,o-z),integer(i-n)
  include 'netcdf.inc'

  parameter(nbyte=4)
  character*30 file63
  character*12 it_char
  character*48 start_time,version,variable_nm,variable_dim
  character*48 data_format
  character(72) :: fgb,fgb2,fdb  ! Processor specific global output file name
  integer :: lfgb,lfdb       ! Length of processor specific global output file name
  character(len= 4) :: a_4
  integer,allocatable :: elnode(:,:)
  integer,allocatable :: elside(:,:)
  integer,allocatable :: isdel(:,:)
  allocatable ne(:),np(:),ns(:),ihot_len(:)
  allocatable ztot(:),sigma(:),cs(:),outb(:,:,:),eta2(:),outeb(:,:,:),outsb(:,:,:)
  allocatable i34(:),nm2(:,:),xctr(:),yctr(:),dpe(:)
  allocatable x(:),y(:),dp(:),kbp00(:),iplg(:,:),ielg(:,:),kbp01(:,:)
  allocatable islg(:,:),kbs(:),kbe(:),xcj(:),ycj(:),dps(:),ic3(:,:),isidenode(:,:)
  allocatable nm_e(:,:),nm_s(:,:),eta2s(:),eta2e(:)

  !netcdf variables
  character(len=50) fname
  character(len=256) cbuffer, stdname,units,longname
  character(len=8) varname,vname,uname
  integer :: time_dims(1),ele_dims(2),x_dims(1),y_dims(1),z_dims(1),sigma_dims(1), &
            &var1d_dims(1),var2d_dims(2),var3d_dims(3), &
            &data_start_1d(1),data_start_2d(2),data_start_3d(3), &
            &data_count_1d(1),data_count_2d(2),data_count_3d(3), &
            & int_buffer(4)
  real*4  :: hc_array(1),hs_array(1),thetab_array(1),thetaf_array(1),real_buffer(4) 
 
 

  
      
!-------------------------------------------------------------------------------
! Aquire user inputs
!-------------------------------------------------------------------------------

  open(10,file='combine_output.in',status='old')
  read(10,'(a30)') file63 !e.g. 'hvel.64'
  read(10,*) ibgn,iend
  read(10,*) inc !inc=1: netcdf option
  close(10)

! Read local_to_global_0000 for global info
  open(10,file='local_to_global_0000',status='old')
  read(10,*)ne_global,np_global,nvrt,nproc,ntracers
  close(10)

  allocate(x(np_global),y(np_global),dp(np_global),kbp00(np_global),kbe(ne_global), &
           np(0:nproc-1),ns(0:nproc-1),ne(0:nproc-1),elnode(3,ne_global),nm2(ne_global,3),eta2(np_global), &
           ztot(nvrt),sigma(nvrt),cs(nvrt),outb(np_global,nvrt,2),ihot_len(0:nproc-1), &
           outeb(ne_global,nvrt,2),dpe(ne_global),xctr(ne_global),yctr(ne_global), &
           eta2e(ne_global),stat=istat)
  if(istat/=0) stop 'Allocation error: x,y'

! Initialize outb for ivalid pts (below bottom etc)
  outb=-9999.
  outeb=-9999.

!-------------------------------------------------------------------------------
! Read rank-specific local_to_global*
!-------------------------------------------------------------------------------

  ! Compute ivs and i23d; i23d_out for vis
  file63=adjustl(file63)
  lfile63=len_trim(file63)
  if(file63((lfile63-1):lfile63).eq.'61'.or.file63((lfile63-1):lfile63).eq.'63') then
    ivs=1
  else if(file63((lfile63-1):lfile63).eq.'62'.or.file63((lfile63-1):lfile63).eq.'64') then
    ivs=2
  else 
    write(it_char,'(i12)')ibgn
    it_char=adjustl(it_char)  !place blanks at end
    it_len=len_trim(it_char)  !length without trailing blanks
    fgb=it_char(1:it_len)//'_0000'; lfgb=len_trim(fgb)
    open(63,file=fgb(1:lfgb)//'_'//file63,access='direct',recl=nbyte,status='old')
    read(63,rec=3)ivs
    close(63)
    print*, 'ivs=',ivs
  endif
  if(file63((lfile63-1):lfile63).eq.'61'.or.file63((lfile63-1):lfile63).eq.'62') then
    i23d=2; i23d_out=2
  else if(file63((lfile63-1):lfile63).eq.'63'.or.file63((lfile63-1):lfile63).eq.'64') then
    i23d=3; i23d_out=3
  else if(file63((lfile63-1):lfile63).eq.'65') then
    i23d=4; i23d_out=2 !2D side
  else if(file63((lfile63-1):lfile63).eq.'66') then
    i23d=5; i23d_out=2 !2D element
  else if(file63((lfile63-1):lfile63).eq.'67') then
    i23d=6; i23d_out=3 !3D side and whole level
  else if(file63((lfile63-1):lfile63).eq.'68') then
    i23d=7; i23d_out=3 !3D side and half level
  else if(file63((lfile63-1):lfile63).eq.'69') then
    i23d=8; i23d_out=3 !3D elem. and whole level
  else if(file63((lfile63-1):lfile63).eq.'70') then
    i23d=9; i23d_out=3 !3D elem. and half level
  else if(file63((lfile63-1):lfile63).eq.'71') then
    i23d=10; i23d_out=2 !2D node
  else if(file63((lfile63-1):lfile63).eq.'72') then
    i23d=11; i23d_out=3 !3D node and whole level
  else if(file63((lfile63-1):lfile63).eq.'73') then
    i23d=12; i23d_out=3 !3D node and half level
  else 
    stop 'Unknown suffix'
  endif
   
  !Read sidecenters.gr3 or centers.gr3 for non-standard outputs to get conectivity tables
  if(i23d==4.or.i23d==6.or.i23d==7) then
    open(14,file='../sidecenters.gr3',status='old')
    read(14,*); read(14,*)ns_e,ns1
    allocate(nm_s(ns_e,3),stat=istat)
    if(istat/=0) stop 'Allocation error: side grid'
    do i=1,ns1; read(14,*); enddo
    do i=1,ns_e
      read(14,*)j,k,nm_s(i,1:3)
    enddo !i
    close(14)
  endif !i23d

  if(i23d==5.or.i23d==8.or.i23d==9) then
    open(14,file='../centers.gr3',status='old')
    read(14,*); read(14,*)ne_e,ne1
    if(ne1/=ne_global) then
      print*, 'centers.gr3 not consistent with grid:',ne1,ne_global
      stop
    endif
    allocate(nm_e(ne_e,3),stat=istat)
    if(istat/=0) stop 'Allocation error: center grid'
    do i=1,ne1; read(14,*); enddo
    do i=1,ne_e
      read(14,*)j,k,nm_e(i,1:3)
    enddo !i
    close(14)
  endif !i23d

  ! Read in local-global mappings from all ranks
  fdb='local_to_global_0000'
  lfdb=len_trim(fdb)

  !Find max. for dimensioning
  do irank=0,nproc-1
    write(fdb(lfdb-3:lfdb),'(i4.4)') irank
    open(10,file=fdb,status='old')
    read(10,*) !global info
    read(10,*) !info
    read(10,*)ne(irank)
    do i=1,ne(irank)
      read(10,*)!j,ielg(irank,i)
    enddo !i
    read(10,*)np(irank)
    do i=1,np(irank)
      read(10,*)
    enddo !i
    read(10,*)ns(irank)
    close(10)
  enddo !irank
  np_max=maxval(np(:))
  ns_max=maxval(ns(:))
  ne_max=maxval(ne(:))

  allocate(iplg(0:nproc-1,np_max),kbp01(0:nproc-1,np_max), &
     &ielg(0:nproc-1,ne_max),islg(0:nproc-1,ns_max),stat=istat)
  if(istat/=0) stop 'Allocation error (2)'


  !Re-read
  ns_global=0
  do irank=0,nproc-1
    write(fdb(lfdb-3:lfdb),'(i4.4)') irank
    open(10,file=fdb,status='old')

    read(10,*) !global info

    read(10,*) !info
    read(10,*)ne(irank)
    do i=1,ne(irank)
      read(10,*)j,ielg(irank,i)
    enddo !i
    read(10,*)np(irank)
    do i=1,np(irank)
      read(10,*)j,iplg(irank,i)
    enddo
    read(10,*)ns(irank) !sides
    do i=1,ns(irank)
      read(10,*)j,islg(irank,i)
      if(ns_global<islg(irank,i)) ns_global=islg(irank,i)
    enddo

    read(10,*) !'Header:'
    read(10,'(a)')data_format,version,start_time
    read(10,*)nrec,dtout,nspool,nvrt,kz,h0,h_s,h_c,theta_b,theta_f
    read(10,*)(ztot(k),k=1,kz-1),(sigma(k),k=1,nvrt-kz+1)
    read(10,*)np(irank),ne(irank),(x(iplg(irank,m)),y(iplg(irank,m)),dp(iplg(irank,m)),kbp01(irank,m),m=1,np(irank)), &
    &         (ntmp,(nm2(m,mm),mm=1,3),m=1,ne(irank))

    close(10)

!   Compute C(s) for output
    do klev=kz,nvrt
      k=klev-kz+1
      cs(k)=(1-theta_b)*sinh(theta_f*sigma(k))/sinh(theta_f)+ &
     &theta_b*(tanh(theta_f*(sigma(k)+0.5))-tanh(theta_f*0.5))/2/tanh(theta_f*0.5)
    enddo !klev

!   Compute kbp00 (to avoid mismatch of indices) - larger rank prevails
    do m=1,np(irank)
      ipgb=iplg(irank,m)
      kbp00(ipgb)=kbp01(irank,m)
    enddo !m
 
!   Reconstruct connectivity table
    do m=1,ne(irank)
      iegb=ielg(irank,m)
      if(iegb>ne_global) stop 'Overflow!'
      do mm=1,3
        itmp=nm2(m,mm)
        if(itmp>np(irank).or.itmp<=0) then
          write(*,*)'Overflow:',m,mm,itmp
          stop
        endif
        elnode(mm,iegb)=iplg(irank,itmp)
      enddo !mm
    enddo !m
  enddo !irank=0,nproc-1

! Compute geometry
  call compute_nside(np_global,ne_global,elnode,ns2)
  allocate(ic3(3,ne_global),elside(3,ne_global),isdel(2,ns2),isidenode(2,ns2),xcj(ns2),ycj(ns2),stat=istat)
  if(istat/=0) stop 'Allocation error: side(0)'
  call selfe_geometry(np_global,ne_global,ns2,x,y,elnode,ic3,elside,isdel,isidenode,xcj,ycj)
  if(ns2/=ns_global) then
    write(*,*)'Mismatch in side:',ns2,ns_global
    stop
  endif

  if(i23d==4.or.i23d==6.or.i23d==7) then
    if(ns1/=ns2) then
      write(*,*)'Mismatch in side:',ns1,ns2
      stop
    endif
  endif

! Allocate side arrays
  allocate(dps(ns_global),kbs(ns_global),outsb(ns_global,nvrt,2),eta2s(ns_global),stat=istat)
  if(istat/=0) stop 'Allocation error: side'
  outsb=-9999.

! Compute side/element bottom index
  do i=1,ne_global
    kbe(i)=maxval(kbp00(elnode(1:3,i)))
    dpe(i)=sum(dp(elnode(1:3,i)))/3.d0
    xctr(i)=sum(x(elnode(1:3,i)))/3.d0
    yctr(i)=sum(y(elnode(1:3,i)))/3.d0
  enddo !i
  do i=1,ns_global
    kbs(i)=maxval(kbp00(isidenode(1:2,i)))
    dps(i)=sum(dp(isidenode(1:2,i)))/2
  enddo !i
! Compute record length to read from each rank-specific binary output per time step
  do irank=0,nproc-1
    ihot_len(irank)=nbyte*(2+np(irank)) !time,it,eta
    if(i23d>3) ihot_len(irank)=ihot_len(irank)+nbyte !ivs

    if(i23d==2.or.i23d==10) then
      ihot_len(irank)=ihot_len(irank)+nbyte*ivs*np(irank)
    else if(i23d==3) then
      do i=1,np(irank)
        do k=max0(1,kbp01(irank,i)),nvrt
          do m=1,ivs
            ihot_len(irank)=ihot_len(irank)+nbyte
          enddo !m
        enddo !k
      enddo !i
    else if(i23d==4) then !2D side
      ihot_len(irank)=ihot_len(irank)+nbyte*ns(irank)*ivs
    else if(i23d==5) then !2D elem
      ihot_len(irank)=ihot_len(irank)+nbyte*ne(irank)*ivs
    else if(i23d==6.or.i23d==7) then !3D side and whole/half level
      ihot_len(irank)=ihot_len(irank)+nbyte*ns(irank)*nvrt*ivs
    else if(i23d==8.or.i23d==9) then !3D element and whole/half level 
      ihot_len(irank)=ihot_len(irank)+nbyte*ne(irank)*nvrt*ivs
    else if(i23d==11.or.i23d==12) then !3D node and whole/half level 
      ihot_len(irank)=ihot_len(irank)+nbyte*np(irank)*nvrt*ivs
    endif
  enddo !irank

!-------------------------------------------------------------------------------
! Time iteration -- select "node" data
!-------------------------------------------------------------------------------
 
 
! Loop over input files
  do iinput=ibgn,iend
    write(it_char,'(i12)')iinput
    it_char=adjustl(it_char)  !place blanks at end
    it_len=len_trim(it_char)  !length without trailing blanks
    fgb=it_char(1:it_len)//'_0000'; lfgb=len_trim(fgb);

    ! Read actual number of spools in this file
!    open(63,file=fgb(1:lfgb)//'_'//file63,access='direct',recl=nbyte,status='old')
!    read(63,rec=irec_nrec) nrec
!    close(63)

    !Write header to output files 
    if(inc==0) then
!     open(65,file=it_char(1:it_len)//'_'//file63,status='replace')
      open(65,file=it_char(1:it_len)//'_'//file63,status='replace',form="unformatted",access="stream") !,buffered="yes")
      data_format='DataFormat v5.0'
      variable_nm=file63 !not important
      variable_dim=file63

      write(65) data_format
      write(65) version
      write(65) start_time
      write(65) variable_nm
      write(65) variable_dim

      
      write(65) nrec
      write(65) dtout
      write(65) nspool
      write(65) ivs
      write(65) i23d_out

      !Vertical grid
      write(65) nvrt
      write(65) kz
      write(65) h0
      write(65) h_s
      write(65) h_c
      write(65) theta_b
      write(65) theta_f


      write(65) (ztot(k),k=1,kz-1)
      write(65)  (sigma(kin),kin=1,nvrt-kz+1)

      !Horizontal grid
      if(i23d<=3.or.i23d>=10) then !standard
        write(65) np_global
        write(65) ne_global
        write(65) ((x(m),y(m),dp(m),kbp00(m)),m=1,np_global)
        !write(65) ((3,(elnode(mm,m),mm=1,3)),m=1,ne_global) 
        !the line above doing the same job,but not easy to understand
        do m=1,ne_global
          write(65) 3, (elnode(mm,m),mm=1,3)
        enddo !m
      else if(i23d==4.or.i23d==6.or.i23d==7) then !side based
        write(65)  ns_global
        write(65) ns_e
        write(65) ((xcj(m),ycj(m),dps(m),kbs(m)),m=1,ns_global)
        do m=1,ns_e
          write(65) 3,(nm_s(m,mm),mm=1,3)
        enddo !m
      else !elem. based
        write(65) ne_global
        write(65) ne_e
        write(65) ((xctr(m),yctr(m),dpe(m),kbe(m)),m=1,ne_global)
        do m=1,ne_e
          write(65) 3,(nm_e(m,mm),mm=1,3)
        enddo !m
      endif !i23d

    else !netcdf
!     enter define mode
      fname=it_char(1:it_len)//'_'//file63(1:lfile63-3)//'.nc'
      iret = nf_create(trim(fname), NF_CLOBBER, ncid)
!     define dimensions
      iret = nf_def_dim(ncid, 'node',np_global, node_dim)
      iret = nf_def_dim(ncid, 'nele',ne_global, nele_dim)
      iret = nf_def_dim(ncid, 'nface',3, nface_dim)
      iret = nf_def_dim(ncid, 'nv',nvrt, nv_dim)
      iret = nf_def_dim(ncid, 'one',1, none_dim)
      iret = nf_def_dim(ncid, 'sigma',nvrt-kz+1, nsigma_dim)
      if(kz/=1) iret = nf_def_dim(ncid, 'nz',kz-1, nz_dim)
      iret = nf_def_dim(ncid, 'time', NF_UNLIMITED, ntime_dim)

!     define variables
      time_dims(1) = ntime_dim
      iret=nf_def_var(ncid,'time',NF_FLOAT,1,time_dims,itime_id)
      if(iret.ne.NF_NOERR) then
         print*, nf_strerror(iret); stop
      endif
      read(start_time(7:10),*) int_buffer(1)  ! Year
      read(start_time(4:6), *) int_buffer(2)  ! Month
      read(start_time(1:2), *) int_buffer(3)  ! Day
      read(start_time(12:13), *) int_buffer(4)  ! Hour
      write(cbuffer,20) int_buffer(1),int_buffer(2),int_buffer(3),start_time(11:len(start_time))
20    FORMAT('seconds since ',I0.4,'-',I0.2,'-',I0.2,' ',A38)
      iret=nf_put_att_text(ncid,itime_id,'long_name',4,'Time')
      iret=nf_put_att_text(ncid,itime_id,'units',len(trim(cbuffer)),cbuffer)
      iret=nf_put_att_text(ncid,itime_id,'base_date',len(start_time),start_time)
      iret=nf_put_att_text(ncid,itime_id,'standard_name',4,'time')
!     write time stamps later

      ele_dims(1)=nface_dim; ele_dims(2)=nele_dim
      iret=nf_def_var(ncid,'ele',NF_INT,2,ele_dims,iele_id)
      iret=nf_put_att_text(ncid,iele_id,'long_name',35,'Horizontal Triangular Element Table')
      iret=nf_put_att_text(ncid,iele_id,'units',15,'non-dimensional')
      iret=nf_put_att_text(ncid,iele_id,'cf_role',22,'face_node_connectivity')
      int_buffer(1) = 1  ! fortran indexing starts at 1
      iret=nf_put_att_int(ncid,iele_id,'start_index',NF_INT,1,int_buffer)
      x_dims(1)=node_dim

      iret=nf_def_var(ncid,'x',NF_FLOAT,1,x_dims,ix_id)
      iret=nf_put_att_text(ncid,ix_id,'long_name',13,'x-coordinates')
      iret=nf_put_att_text(ncid,ix_id,'units',6,'meters')
      iret=nf_put_att_text(ncid,ix_id,'mesh',4,'mesh')
      iret=nf_put_att_text(ncid,ix_id,'location',4,'node')

      iret=nf_def_var(ncid,'y',NF_FLOAT,1,x_dims,iy_id)
      iret=nf_put_att_text(ncid,iy_id,'long_name',13,'y-coordinates')
      iret=nf_put_att_text(ncid,iy_id,'units',6,'meters')
      iret=nf_put_att_text(ncid,iy_id,'mesh',4,'mesh')
      iret=nf_put_att_text(ncid,iy_id,'location',4,'node')

      iret=nf_def_var(ncid,'depth',NF_FLOAT,1,x_dims,idepth_id)
      iret=nf_put_att_text(ncid,idepth_id,'long_name',10,'Bathymetry')
      iret=nf_put_att_text(ncid,idepth_id,'units',6,'meters')
      iret=nf_put_att_text(ncid,idepth_id,'positive',4,'down')
      iret=nf_put_att_text(ncid,idepth_id,'mesh',4,'mesh')
      iret=nf_put_att_text(ncid,idepth_id,'location',4,'node')

      sigma_dims(1)=nsigma_dim
      ! See http://cf-pcmdi.llnl.gov/documents/cf-conventions/1.7-draft1/apd.html 
      ! section "Ocean s-coordinate"  for the CF calculation
      ! See trunk/src/Utility/Post-Processing-Fortran/compute_zcor.f90 for the SELFE calc.
      ! CF:SELFE corresponence: depth_c:h_c s:sigma C(k):cs(k) a:theta_f b:theta_b 
      hs_array(1)=h_s
      hc_array(1)=h_c
      thetab_array(1)=theta_b
      thetaf_array(1)=theta_f
      iret=nf_def_var(ncid,'sigma',NF_FLOAT,1,sigma_dims,isigma_id)
      iret=nf_put_att_text(ncid,isigma_id,'long_name',29,'S coordinates at whole levels')
      iret=nf_put_att_text(ncid,isigma_id,'units',1,'1')
      cbuffer='ocean_s_coordinate'
      iret=nf_put_att_text(ncid,isigma_id,'standard_name',len(trim(cbuffer)),cbuffer)
      iret=nf_put_att_text(ncid,isigma_id,'positive',2,'up')
      print *, 'h_s, h_c, theta_b theta_f:', hs_array(1), hc_array(1),thetab_array(1),thetaf_array(1)
      iret=nf_put_att_real(ncid,isigma_id,'h_s',NF_FLOAT,1,hs_array)
      iret=nf_put_att_real(ncid,isigma_id,'h_c',NF_FLOAT,1,hc_array)
      iret=nf_put_att_real(ncid,isigma_id,'theta_b',NF_FLOAT,1,thetab_array)
      iret=nf_put_att_real(ncid,isigma_id,'theta_f',NF_FLOAT,1,thetaf_array)
      cbuffer='s: sigma eta: elev depth: depth a: sigma_theta_f b: sigma_theta_b depth_c: sigma_h_c'
      iret=nf_put_att_text(ncid,isigma_id,'formula_terms',len(trim(cbuffer)),cbuffer)

      var1d_dims(1)=none_dim
      iret=nf_def_var(ncid,'sigma_h_c',NF_FLOAT,1,var1d_dims,ihc_id)
      if(iret.ne.NF_NOERR) then
         print*, nf_strerror(iret); stop
      endif
      cbuffer='ocean_s_coordinate h_c constant'
      iret=nf_put_att_text(ncid,ihc_id,'long_name',len(trim(cbuffer)),cbuffer)
      iret=nf_put_att_text(ncid,idepth_id,'units',6,'meters')
      iret=nf_put_att_text(ncid,idepth_id,'positive',4,'down')

      var1d_dims(1)=none_dim
      iret=nf_def_var(ncid,'sigma_theta_b',NF_FLOAT,1,var1d_dims,itheta_b_id)
      cbuffer='ocean_s_coordinate theta_b constant'
      iret=nf_put_att_text(ncid,itheta_a_id,'long_name',len(trim(cbuffer)),cbuffer)

      var1d_dims(1)=none_dim
      iret=nf_def_var(ncid,'sigma_theta_f',NF_FLOAT,1,var1d_dims,itheta_f_id)
      cbuffer='ocean_s_coordinate theta_f constant'
      iret=nf_put_att_text(ncid,itheta_f_id,'long_name',len(trim(cbuffer)),cbuffer)

      var1d_dims(1)=none_dim
      iret=nf_def_var(ncid,'sigma_maxdepth',NF_FLOAT,1,var1d_dims,ihs_id)
      cbuffer='ocean_s_coordinate maximum depth cutoff (mixed s over z boundary'
      iret=nf_put_att_text(ncid,ihs_id,'long_name',len(trim(cbuffer)),cbuffer)
      iret=nf_put_att_text(ncid,idepth_id,'units',6,'meters')
      iret=nf_put_att_text(ncid,idepth_id,'positive',4,'down')

       
      iret=nf_def_var(ncid,'Cs',NF_FLOAT,1,sigma_dims,ics_id)
      iret=nf_put_att_text(ncid,ics_id,'long_name',29,'Function C(s) at whole levels')
      iret=nf_put_att_text(ncid,ics_id,'units',15,'non-dimensional')
      iret=nf_put_att_text(ncid,ics_id,'positive',2,'up')

      if(kz/=1) then
        z_dims(1)=nz_dim
        iret=nf_def_var(ncid,'z',NF_FLOAT,1,z_dims,iz_id)
        iret=nf_put_att_text(ncid,iz_id,'long_name',29,'Z coordinates at whole levels')
        iret=nf_put_att_text(ncid,iz_id,'units',6,'meters')
        iret=nf_put_att_text(ncid,iz_id,'positive',2,'up')
      endif

      if(i23d==2) then
        var2d_dims(1)=node_dim; var2d_dims(2)=ntime_dim
        if(file63((lfile63-1):lfile63).eq.'62') then
          iret=nf_def_var(ncid,'dahv_u',NF_FLOAT,2,var2d_dims,iu_id)
          iret=nf_put_att_text(ncid,iu_id,'long_name',27,'Eastward Velocity or stress')
          iret=nf_put_att_text(ncid,iu_id,'missing_value',5,'-999.')
          iret=nf_put_att_text(ncid,iu_id,'mesh',4,'mesh')
          iret=nf_put_att_text(ncid,iu_id,'location',4,'node')
          iret=nf_def_var(ncid,'dahv_v',NF_FLOAT,2,var2d_dims,iv_id)
          iret=nf_put_att_text(ncid,iv_id,'long_name',28,'Northward Velocity or stress')
          iret=nf_put_att_text(ncid,iv_id,'missing_value',5,'-999.')
          iret=nf_put_att_text(ncid,iv_id,'mesh',4,'mesh')
          iret=nf_put_att_text(ncid,iv_id,'location',4,'node')
        else
          iret=nf_def_var(ncid,file63(1:lfile63-3),NF_FLOAT,2,var2d_dims,ivar2d_id)
          iret=nf_put_att_text(ncid,ivar2d_id,'long_name',lfile63,file63(1:lfile63))
          iret=nf_put_att_text(ncid,ivar2d_id,'mesh',4,'mesh')
          iret=nf_put_att_text(ncid,ivar2d_id,'location',4,'node')
        endif !file63
      else !3D
        var3d_dims(1)=node_dim; var3d_dims(2)=nv_dim; var3d_dims(3)=ntime_dim
        if(file63((lfile63-1):lfile63).eq.'64') then
          iret=nf_def_var(ncid,'u',NF_FLOAT,3,var3d_dims,iu_id)
          iret=nf_put_att_text(ncid,iu_id,'long_name',23,'Eastward Water Velocity')
          iret=nf_put_att_text(ncid,iu_id,'missing_value',5,'-999.')
          iret=nf_put_att_text(ncid,iu_id,'mesh',4,'mesh')
          iret=nf_put_att_text(ncid,iu_id,'location',4,'node')
          iret=nf_def_var(ncid,'v',NF_FLOAT,3,var3d_dims,iv_id)
          iret=nf_put_att_text(ncid,iv_id,'long_name',24,'Northward Water Velocity')
          iret=nf_put_att_text(ncid,iv_id,'missing_value',5,'-999.')
          iret=nf_put_att_text(ncid,iv_id,'mesh',4,'mesh')
          iret=nf_put_att_text(ncid,iv_id,'location',4,'node')
        else
          iret=nf_def_var(ncid,file63(1:lfile63-3),NF_FLOAT,3,var3d_dims,ivar3d_id)
          iret=nf_put_att_text(ncid,ivar3d_id,'long_name',lfile63,file63(1:lfile63))
          iret=nf_put_att_text(ncid,ivar3d_id,'missing_value',5,'-999.')
          iret=nf_put_att_text(ncid,ivar3d_id,'mesh',4,'mesh')
          iret=nf_put_att_text(ncid,ivar3d_id,'location',4,'node')
        endif
      endif

      iref=nf_def_var(ncid,'mesh',NF_INT,0,ele_dims,imesh_id)
      iret=nf_put_att_text(ncid,imesh_id,'cf_role',13,'mesh_topology')
      iret=nf_put_att_text(ncid,imesh_id,'long_name',37,'Topology data of 2d unstructured mesh')
      int_buffer(1) = 2  ! 2d mesh dimension
      iret=nf_put_att_int(ncid,imesh_id,'topology_dimension',NF_INT,1,int_buffer)
      cbuffer='x y'
      iret=nf_put_att_text(ncid,imesh_id,'node_coordinates',len(trim(cbuffer)),cbuffer)
      cbuffer='ele'
      iret=nf_put_att_text(ncid,imesh_id,'face_node_connectivity',len(trim(cbuffer)),cbuffer)

      iret = nf_put_att_text(ncid, NF_GLOBAL, 'Conventions', 6, 'CF-1.0')
!     leave define mode
      iret = nf_enddef(ncid)

!     Write mode (header part only)
      data_start_2d(1:2)=1
      data_count_2d(1)=3; data_count_2d(2)=ne_global
      iret=nf_put_vara_int(ncid,iele_id,data_start_2d,data_count_2d,elnode)
      if(iret.ne.NF_NOERR) then
         print*, nf_strerror(iret); stop
      endif
      iret=nf_put_vara_real(ncid,ix_id,1,np_global,x)
      iret=nf_put_vara_real(ncid,iy_id,1,np_global,y)
      iret=nf_put_vara_real(ncid,idepth_id,1,np_global,dp)
      iret=nf_put_vara_real(ncid,isigma_id,1,nvrt-kz+1,sigma)
      iret=nf_put_vara_real(ncid,ics_id,1,nvrt-kz+1,cs)
      iret=nf_put_vara_real(ncid,ihc_id,1,1,h_c)
      iret=nf_put_vara_real(ncid,itheta_b_id,1,1,theta_b)
      iret=nf_put_vara_real(ncid,itheta_f_id,1,1,theta_f)
      iret=nf_put_vara_real(ncid,ihs_id,1,1,h_s)
      if(kz/=1) iret=nf_put_vara_real(ncid,iz_id,1,kz-1,ztot)
    endif !inc
 
    !print*, 'Last element:',elnode(1:3,ne_global)
    !end output header
   
    ! Loop over output spools in file
    do ispool=1,nrec
       
      !Gather all ranks
      do irank=0,nproc-1
        !Open input file
        fgb2=fgb
        write(fgb2(lfgb-3:lfgb),'(i4.4)') irank

        open(63,file=fgb2(1:lfgb)//'_'//file63,access='direct',recl=ihot_len(irank),status='old')
 
        if(i23d==2) then
          read(63,rec=ispool)time,it,(eta2(iplg(irank,i)),i=1,np(irank)),((outb(iplg(irank,i),1,m),m=1,ivs),i=1,np(irank))
        else if(i23d==3) then !3D
          read(63,rec=ispool)time,it,(eta2(iplg(irank,i)),i=1,np(irank)), &
     &                       (((outb(iplg(irank,i),k,m),m=1,ivs),k=max0(1,kbp01(irank,i)),nvrt),i=1,np(irank))
        else if(i23d==4) then !2D side 
          read(63,rec=ispool)time,it,itmp,(eta2(iplg(irank,i)),i=1,np(irank)), &
     &                       (outsb(islg(irank,i),1,1:ivs),i=1,ns(irank))
        else if(i23d==5) then !2D elem
          read(63,rec=ispool)time,it,itmp,(eta2(iplg(irank,i)),i=1,np(irank)), &
     &                       (outeb(ielg(irank,i),1,1:ivs),i=1,ne(irank))
        else if(i23d==6.or.i23d==7) then !3D side and whole/half level
          read(63,rec=ispool)time,it,itmp,(eta2(iplg(irank,i)),i=1,np(irank)), &
     &                       (((outsb(islg(irank,i),k,m),m=1,ivs),k=1,nvrt),i=1,ns(irank))
        else if(i23d==8.or.i23d==9) then !3D element at whole/half levels
          read(63,rec=ispool)time,it,itmp,(eta2(iplg(irank,i)),i=1,np(irank)), &
     &                       (((outeb(ielg(irank,i),k,m),m=1,ivs),k=1,nvrt),i=1,ne(irank))
        else if(i23d==10) then !2D node (n.s.)
          read(63,rec=ispool)time,it,itmp,(eta2(iplg(irank,i)),i=1,np(irank)), &
     &((outb(iplg(irank,i),1,m),m=1,ivs),i=1,np(irank))
        else !3D node
          read(63,rec=ispool)time,it,itmp,(eta2(iplg(irank,i)),i=1,np(irank)), &
     &                       (((outb(iplg(irank,i),k,m),m=1,ivs),k=1,nvrt),i=1,np(irank))
        endif
        ! Close input file
        close(63)
      enddo !irank

      
      !t_loop_over_output=t_loop_over_output+t4-t3

      !Output
      !print*, 'it=',it,' time=',time/86400.0
     
      if(inc==0) then !binary
        !Compute eta2 at sides/elem.
        do i=1,ns_global
          n1=isidenode(1,i); n2=isidenode(2,i)
          eta2s(i)=(eta2(n1)+eta2(n2))/2.d0
        enddo !i
        do i=1,ne_global
          eta2e(i)=sum(eta2(elnode(1:3,i)))/3.d0
        enddo !i

       
        write(65) time
        write(65) it
      
        if(i23d<4.or.i23d>=10) then !node based
         
          write(65) (eta2(i),i=1,np_global)
          if(i23d==2.or.i23d==10) then
             write(65) ((outb(i,1,m),m=1,ivs),i=1,np_global)
          else !if(i23d==3) then
             write(65) (((outb(i,k,m),m=1,ivs),k=max0(1,kbp00(i)),nvrt),i=1,np_global)
          endif !i23d
         
        else if(i23d==4) then !2D side
          write(65) (eta2s(i),i=1,ns_global)
          write(65) ((outsb(i,1,m),m=1,ivs),i=1,ns_global)
        else if(i23d==5) then !2D elem
          write(65) (eta2e(i),i=1,ne_global)
          write(65) ((outeb(i,1,m),m=1,ivs),i=1,ne_global)
        else if(i23d==6.or.i23d==7) then !3D side and whole/half level
          write(65) (eta2s(i),i=1,ns_global)
          write(65) (((outsb(i,k,m),m=1,ivs),k=max0(1,kbs(i)),nvrt),i=1,ns_global)
            !Debug
            !if(iinput==iend.and.ispool==nrec) then
            !  write(89,*)xcj(i),ycj(i),outsb(i,nvrt,1:ivs)
            !endif
        else !3D element and whole/half level
          write(65) (eta2e(i),i=1,ne_global)
          write(65) (((outeb(i,k,m),m=1,ivs),k=max0(1,kbe(i)),nvrt),i=1,ne_global)
        endif !i23d
      else !netcdf
        iret=nf_put_vara_real(ncid,itime_id,ispool,1,time)
        if(i23d==2) then
          data_start_2d(1)=1; data_start_2d(2)=ispool
          data_count_2d(1)=np_global; data_count_2d(2)=1
          if(file63((lfile63-1):lfile63).eq.'62') then
            iret=nf_put_vara_real(ncid,iu_id,data_start_2d,data_count_2d,outb(:,1,1))
            iret=nf_put_vara_real(ncid,iv_id,data_start_2d,data_count_2d,outb(:,1,2))
          else
            iret=nf_put_vara_real(ncid,ivar2d_id,data_start_2d,data_count_2d,outb(:,1,1))
          endif
        else !3D
          data_start_3d(1:2)=1; data_start_3d(3)=ispool
          data_count_3d(1)=np_global; data_count_3d(2)=nvrt; data_count_3d(3)=1
          if(file63((lfile63-1):lfile63).eq.'64') then
            iret=nf_put_vara_real(ncid,iu_id,data_start_3d,data_count_3d,outb(:,:,1))
            iret=nf_put_vara_real(ncid,iv_id,data_start_3d,data_count_3d,outb(:,:,2))
          else
            iret=nf_put_vara_real(ncid,ivar3d_id,data_start_3d,data_count_3d,outb(:,:,1))
          endif
        endif !i23d
        print*, 'done record # ',ispool
      endif !inc
       
    
    enddo !ispool=1,nrec
    ! Close output file
    if (inc==0)  then !binary
       close(65)
    else
       iret = nf_close(ncid)
    endif
  enddo ! do iinput=ibgn,iend

  
 


  deallocate(x,y,dp,kbp00,kbe,np,ns,ne,elnode,nm2,eta2,ztot,sigma,cs,outb,ihot_len, &
     &outeb,dpe,xctr,yctr,eta2e,iplg,kbp01,ielg,islg,ic3,elside,isdel,isidenode, &
     &xcj,ycj,dps,kbs,outsb,eta2s)
  if(i23d==4.or.i23d==6.or.i23d==7) deallocate(nm_s)
  if(i23d==5.or.i23d==8.or.i23d==9) deallocate(nm_e)

 
  stop
end program read_iwrite1
