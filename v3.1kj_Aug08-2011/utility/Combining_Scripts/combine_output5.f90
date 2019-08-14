!===============================================================================
! Read in binary outputs (rank-specific) from parallel code and combine them into
! one global output in v5.0 format or netcdf format. Works for partial outputs.
! Gobal-local mappings are read in from separate files.
! Run this program inside the directory outputs/, where some of the input files below
! can be found.

! Inputs:
!        rank-specific binary files (from SELFE outputs); 
!        local_to_global_* (from SELFE outputs);
!        combine_output.in (1st line: elev.61 etc; 2nd line: start and end file indices;
!                          3rd line: inc (1 for netcdf option)); 
! Output: combined binary file (for nc file: e.g. *_salt.nc).
!
!  Compile on canopus01:
!  ifort -Bstatic -O3 -assume byterecl -o combine_output5_canopus combine_output5.f90 -Vaxlib -I/usr/local/include/ -L/usr/local/lib -lnetcdf

!  Compile on amb6402:
!  ifort -Bstatic -O3 -assume byterecl -o combine_output5_canopus combine_output5.f90 -Vaxlib -I/usr/local/netcdf/include/ -L/usr/local/netcdf/lib -lnetcdf
!
!  On Ranger:
!  pgf90 -O2 -mcmodel=medium  -Bstatic -o combine_output5_canopus combine_output5.f90 .....

!  History: (1) added netcdf option.
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
  allocatable ne(:),np(:),nsproc(:),ihot_len(:)
  allocatable ztot(:),sigma(:),cs(:),outb(:,:,:),eta2(:)
  allocatable i34(:),nm(:,:),nm2(:,:),js(:,:),xctr(:),yctr(:),dpe(:)
  allocatable x(:),y(:),dp(:),kbp00(:),iplg(:,:),ielg(:,:),kbp01(:,:)
  allocatable xcj(:),ycj(:),dps(:)

  !netcdf variables
  character(len=50) fname
  integer :: time_dims(1),ele_dims(2),x_dims(1),sigma_dims(1),var2d_dims(2),var3d_dims(3), &
            &data_start_2d(2),data_start_3d(3),data_count_2d(2),data_count_3d(3),z_dims(1)
      
!-------------------------------------------------------------------------------
! Aquire user inputs
!-------------------------------------------------------------------------------

  open(10,file='combine_output.in',status='old')
  read(10,'(a30)') file63
  read(10,*) ibgn,iend
  read(10,*) inc !inc=1: netcdf option
  close(10)

! Read local_to_global_0000 for global info
  open(10,file='local_to_global_0000',status='old')
  read(10,*)ne_global,np_global,nvrt,nproc,ntracers
  close(10)

  allocate(x(np_global),y(np_global),dp(np_global),kbp00(np_global), &
           np(0:nproc-1),ne(0:nproc-1),nm(ne_global,3),nm2(ne_global,3),eta2(np_global), &
           ztot(nvrt),sigma(nvrt),cs(nvrt),outb(np_global,nvrt,2),ihot_len(0:nproc-1),stat=istat)
  if(istat/=0) stop 'Allocation error: x,y'

! Initialize outb for ivalid pts (below bottom etc)
  outb=-9999.

!-------------------------------------------------------------------------------
! Read rank-specific local_to_global*
!-------------------------------------------------------------------------------

  ! Compute ivs and i23d
  file63=adjustl(file63)
  lfile63=len_trim(file63)
  if(file63((lfile63-1):lfile63).eq.'61'.or.file63((lfile63-1):lfile63).eq.'63') then
    ivs=1
  else if(file63((lfile63-1):lfile63).eq.'62'.or.file63((lfile63-1):lfile63).eq.'64') then
    ivs=2
  else
    print*, 'Unknown file type:',file63
  endif
  if(file63((lfile63-1):lfile63).eq.'61'.or.file63((lfile63-1):lfile63).eq.'62') then
    i23d=2
  else
    i23d=3
  endif
   
  ! Read in local-global mappings from all ranks
  fdb='local_to_global_0000'
  lfdb=len_trim(fdb)

  !Find max. for dimensioning
  np_max=0; ne_max=0
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
    close(10)
    np_max=max(np_max,np(irank))
    ne_max=max(ne_max,ne(irank))
  enddo !irank

  allocate(iplg(0:nproc-1,np_max),kbp01(0:nproc-1,np_max), &
     &ielg(0:nproc-1,ne_max),stat=istat)
  if(istat/=0) stop 'Allocation error (2)'

  !Re-read
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
    read(10,*)itmp !sides
    do i=1,itmp
      read(10,*)
    enddo
!    print*, 'Mapping:',irank,ielg(irank,ne(irank))

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

!   Compute kbp00 (to avoid mismatch of indices)
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
        nm(iegb,mm)=iplg(irank,itmp)
      enddo !mm
    enddo !m
  enddo !irank=0,nproc-1

! Compute record length for each rank-specific binary output per time step
  do irank=0,nproc-1
    ihot_len(irank)=nbyte*(2+np(irank))
    if(i23d==2) then
      ihot_len(irank)=ihot_len(irank)+nbyte*ivs*np(irank)
    else
      do i=1,np(irank)
        do k=max0(1,kbp01(irank,i)),nvrt
          do m=1,ivs
            ihot_len(irank)=ihot_len(irank)+nbyte
          enddo !m
        enddo !k
      enddo !i
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
      open(65,file=it_char(1:it_len)//'_'//file63,status='replace')
      data_format='DataFormat v5.0'
      variable_nm=file63 !not important
      variable_dim=file63

      write(65,'(a48)',advance="no") data_format
      write(65,'(a48)',advance="no") version
      write(65,'(a48)',advance="no") start_time
      write(65,'(a48)',advance="no") variable_nm
      write(65,'(a48)',advance="no") variable_dim

      a_4 = transfer(source=nrec,mold=a_4)
      write(65,"(a4)",advance="no") a_4
      a_4 = transfer(source=dtout,mold=a_4)
      write(65,"(a4)",advance="no") a_4
      a_4 = transfer(source=nspool,mold=a_4)
      write(65,"(a4)",advance="no") a_4
      a_4 = transfer(source=ivs,mold=a_4)
      write(65,"(a4)",advance="no") a_4
      a_4 = transfer(source=i23d,mold=a_4)
      write(65,"(a4)",advance="no") a_4

      !Vertical grid
      a_4 = transfer(source=nvrt,mold=a_4)
      write(65,'(a4)',advance="no") a_4
      a_4 = transfer(source=kz,mold=a_4)
      write(65,'(a4)',advance="no") a_4
      a_4 = transfer(source=h0,mold=a_4)
      write(65,'(a4)',advance="no") a_4
      a_4 = transfer(source=h_s,mold=a_4)
      write(65,'(a4)',advance="no") a_4
      a_4 = transfer(source=h_c,mold=a_4)
      write(65,'(a4)',advance="no") a_4
      a_4 = transfer(source=theta_b,mold=a_4)
      write(65,'(a4)',advance="no") a_4
      a_4 = transfer(source=theta_f,mold=a_4)
      write(65,'(a4)',advance="no") a_4

      do k=1,kz-1
        a_4 = transfer(source=ztot(k),mold=a_4)
        write(65,'(a4)',advance="no") a_4
      enddo
      do k=kz,nvrt
        kin=k-kz+1
        a_4 = transfer(source=sigma(kin),mold=a_4)
        write(65,'(a4)',advance="no") a_4
      enddo !k

      !Horizontal grid
      a_4 = transfer(source=np_global,mold=a_4)
      write(65,'(a4)',advance="no") a_4
      a_4 = transfer(source=ne_global,mold=a_4)
      write(65,'(a4)',advance="no") a_4

      do m=1,np_global
        a_4 = transfer(source=x(m),mold=a_4)
        write(65,'(a4)',advance="no") a_4
        a_4 = transfer(source=y(m),mold=a_4)
        write(65,'(a4)',advance="no") a_4
        a_4 = transfer(source=dp(m),mold=a_4)
        write(65,'(a4)',advance="no") a_4
        a_4 = transfer(source=kbp00(m),mold=a_4)
        write(65,'(a4)',advance="no") a_4
      enddo !m=1,np
      do m=1,ne_global
        a_4 = transfer(source=3,mold=a_4)
        write(65,'(a4)',advance="no") a_4
        do mm=1,3
          a_4 = transfer(source=nm(m,mm),mold=a_4)
          write(65,'(a4)',advance="no") a_4
        enddo !mm
      enddo !m

    else !netcdf
!     enter define mode
      fname=it_char(1:it_len)//'_'//file63(1:lfile63-3)//'.nc'
      iret = nf_create(trim(fname), NF_CLOBBER, ncid)
!     define dimensions
      iret = nf_def_dim(ncid, 'node',np_global, node_dim)
      iret = nf_def_dim(ncid, 'nele',ne_global, nele_dim)
      iret = nf_def_dim(ncid, 'nface',3, nface_dim)
      iret = nf_def_dim(ncid, 'nv',nvrt, nv_dim)
      iret = nf_def_dim(ncid, 'sigma',nvrt-kz+1, nsigma_dim)
      if(kz/=1) iret = nf_def_dim(ncid, 'nz',kz-1, nz_dim)
      iret = nf_def_dim(ncid, 'time', NF_UNLIMITED, ntime_dim)

!     define variables
      time_dims(1) = ntime_dim
      iret=nf_def_var(ncid,'time',NF_REAL,1,time_dims,itime_id)
      if(iret.ne.NF_NOERR) then
         print*, nf_strerror(iret); stop
      endif
      iret=nf_put_att_text(ncid,itime_id,'long_name',4,'Time')
      iret=nf_put_att_text(ncid,itime_id,'units',7,'seconds')
      iret=nf_put_att_text(ncid,itime_id,'base_date',len(start_time),start_time)
      iret=nf_put_att_text(ncid,itime_id,'standard_name',4,'time')
!     write time stamps later

      ele_dims(1)=nele_dim; ele_dims(2)=nface_dim
      iret=nf_def_var(ncid,'ele',NF_INT,2,ele_dims,iele_id)
      iret=nf_put_att_text(ncid,iele_id,'long_name',35,'Horizontal Triangular Element Table')
      iret=nf_put_att_text(ncid,iele_id,'units',15,'non-dimensional')

      x_dims(1)=node_dim
      iret=nf_def_var(ncid,'x',NF_REAL,1,x_dims,ix_id)
      iret=nf_put_att_text(ncid,ix_id,'long_name',13,'x-coordinates')
      iret=nf_put_att_text(ncid,ix_id,'units',6,'meters')

      iret=nf_def_var(ncid,'y',NF_REAL,1,x_dims,iy_id)
      iret=nf_put_att_text(ncid,iy_id,'long_name',13,'y-coordinates')
      iret=nf_put_att_text(ncid,iy_id,'units',6,'meters')

      iret=nf_def_var(ncid,'depth',NF_REAL,1,x_dims,idepth_id)
      iret=nf_put_att_text(ncid,idepth_id,'long_name',10,'Bathymetry')
      iret=nf_put_att_text(ncid,idepth_id,'units',6,'meters')
      iret=nf_put_att_text(ncid,idepth_id,'positive',6,'down')

      sigma_dims(1)=nsigma_dim
      iret=nf_def_var(ncid,'sigma',NF_REAL,1,sigma_dims,isigma_id)
      iret=nf_put_att_text(ncid,isigma_id,'long_name',29,'S coordinates at whole levels')
      iret=nf_put_att_text(ncid,isigma_id,'units',15,'non-dimensional')
      iret=nf_put_att_text(ncid,isigma_id,'positive',2,'up')

      iret=nf_def_var(ncid,'Cs',NF_REAL,1,sigma_dims,ics_id)
      iret=nf_put_att_text(ncid,ics_id,'long_name',29,'Function C(s) at whole levels')
      iret=nf_put_att_text(ncid,ics_id,'units',15,'non-dimensional')
      iret=nf_put_att_text(ncid,ics_id,'positive',2,'up')

      if(kz/=1) then
        z_dims(1)=nz_dim
        iret=nf_def_var(ncid,'z',NF_REAL,1,z_dims,iz_id)
        iret=nf_put_att_text(ncid,iz_id,'long_name',29,'Z coordinates at whole levels')
        iret=nf_put_att_text(ncid,iz_id,'units',6,'meters')
        iret=nf_put_att_text(ncid,iz_id,'positive',2,'up')
      endif

      if(i23d==2) then
        var2d_dims(1)=node_dim; var2d_dims(2)=ntime_dim
        iret=nf_def_var(ncid,file63(1:lfile63-3),NF_REAL,2,var2d_dims,ivar2d_id)
        iret=nf_put_att_text(ncid,ivar2d_id,'long_name',lfile63,file63(1:lfile63))
      else !3D
        var3d_dims(1)=node_dim; var3d_dims(2)=nv_dim; var3d_dims(3)=ntime_dim
        if(file63((lfile63-1):lfile63).eq.'64') then
          iret=nf_def_var(ncid,'u',NF_REAL,3,var3d_dims,iu_id)
          iret=nf_put_att_text(ncid,iu_id,'long_name',23,'Eastward Water Velocity')
          iret=nf_put_att_text(ncid,iu_id,'missing_value',6,'-9999.')
          iret=nf_def_var(ncid,'v',NF_REAL,3,var3d_dims,iv_id)
          iret=nf_put_att_text(ncid,iv_id,'long_name',24,'Northward Water Velocity')
          iret=nf_put_att_text(ncid,iv_id,'missing_value',6,'-9999.')
        else
          iret=nf_def_var(ncid,file63(1:lfile63-3),NF_REAL,3,var3d_dims,ivar3d_id)
          iret=nf_put_att_text(ncid,ivar3d_id,'long_name',lfile63,file63(1:lfile63))
          iret=nf_put_att_text(ncid,ivar3d_id,'missing_value',6,'-9999.')
        endif
      endif

      iret = nf_put_att_text(ncid, NF_GLOBAL, 'Conventions', 6, 'CF-1.0')
!     leave define mode
      iret = nf_enddef(ncid)

!     Write mode (header part only)
      data_start_2d(1:2)=1
      data_count_2d(1)=ne_global; data_count_2d(2)=3
      iret=nf_put_vara_int(ncid,iele_id,data_start_2d,data_count_2d,nm)
      iret=nf_put_vara_real(ncid,ix_id,1,np_global,x)
      iret=nf_put_vara_real(ncid,iy_id,1,np_global,y)
      iret=nf_put_vara_real(ncid,idepth_id,1,np_global,dp)
      iret=nf_put_vara_real(ncid,isigma_id,1,nvrt-kz+1,sigma)
      iret=nf_put_vara_real(ncid,ics_id,1,nvrt-kz+1,cs)
      if(kz/=1) iret=nf_put_vara_real(ncid,iz_id,1,kz-1,ztot)
    endif !inc
 
    !print*, 'Last element:',nm(ne_global,1:3)
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
        else !3D
          read(63,rec=ispool)time,it,(eta2(iplg(irank,i)),i=1,np(irank)), &
     &                       (((outb(iplg(irank,i),k,m),m=1,ivs),k=max0(1,kbp01(irank,i)),nvrt),i=1,np(irank))
        endif
        ! Close input file
        close(63)
      enddo !irank

      !Output
      !print*, 'it=',it,' time=',time/86400.0

      if(inc==0) then !binary
        a_4 = transfer(source=time,mold=a_4)
        write(65,"(a4)",advance="no") a_4
        a_4 = transfer(source=it,mold=a_4)
        write(65,"(a4)",advance="no") a_4

        do i=1,np_global
          a_4 = transfer(source=eta2(i),mold=a_4)
          write(65,"(a4)",advance="no") a_4
        enddo !i
      
        do i=1,np_global
          if(i23d==2) then
            do m=1,ivs
              a_4 = transfer(source=outb(i,1,m),mold=a_4)
              write(65,"(a4)",advance="no") a_4
            enddo !m
          else !i23d=3 
            do k=max0(1,kbp00(i)),nvrt
              do m=1,ivs
                a_4 = transfer(source=outb(i,k,m),mold=a_4)
                write(65,"(a4)",advance="no") a_4
              enddo !m
            enddo !k
          endif !i23d
        enddo !i
      else !netcdf
        iret=nf_put_vara_real(ncid,itime_id,ispool,1,time)
        if(i23d==2) then
          data_start_2d(1)=1; data_start_2d(2)=ispool
          data_count_2d(1)=np_global; data_count_2d(2)=1
          iret=nf_put_vara_real(ncid,ivar2d_id,data_start_2d,data_count_2d,outb(:,1,1))
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
  enddo !iinput=1,ninput_files

! Close output file
  close(65)

  iret = nf_close(ncid)

  stop
end program read_iwrite1
