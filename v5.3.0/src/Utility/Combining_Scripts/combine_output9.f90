!===============================================================================
! Read in binary outputs (rank-specific) from parallel code and combine them into
! one global output in v5.0 format or netcdf format. Works for partial outputs.
! Global-local mappings are read in from separate files.
! Run this program inside the directory outputs/, where some of the input files below
! can be found.

! Inputs:
!        rank-specific binary files (from SCHISM outputs/); 
!        local_to_global_* (from SCHISM outputs/);
!        combine_output.in (1st line: elev.61 etc; 2nd line: start and end file indices;
!                          3rd line: inc (1 for netcdf option)); 
!        additional inputs for non-standard outputs (hvel.67 etc): 
!             sidecenters.gr3, centers.gr3 (in the same dir as hgrid.gr3)
! Output: combined binary file (for nc file: e.g. *_salt.nc; netcdf not working for 
!         non-standard outputs!).
!

!  ifort -O2 -Bstatic -assume byterecl -o combine_output9 combine_output9.f90 netcdf_var_names.f90 ../UtilLib/schism_geometry.f90 ../UtilLib/argparse.F90 -L/nasa/netcdf/3.6.0/intel/lib -lnetcdf -I/nasa/netcdf/3.6.0/intel/include
!  On TYPHOON:
!  pgf90 -O2 -mcmodel=medium -Bstatic -Mpreprocess -o combine_output9 combine_output9.f90 netcdf_var_names.f90 ../UtilLib/schism_geometry.f90 ../UtilLib/argparse.F90 -I/usr/local/amd64a/opteron/pgi/netcdf-3.6.2/include -L/usr/local/amd64a/opteron/pgi/netcdf-3.6.2/lib -lnetcdf -I../UtilLib
!  With gcc:
!  gfortran -O2 -mcmodel=medium -Bstatic -Mpreprocess -o combine_output9 combine_output9.f90 netcdf_var_names.f90 ../UtilLib/schism_geometry.f90 ../UtilLib/argparse.F90 -Mbounds -I/usr/include -I/opt/netcdf-fortran/include -L/usr/lib64 -lnetcdf -L/opt/netcdf-fortran/lib -lnetcdff  -I../UtilLib 
!  History: 
!          2015-12-04  changed mesh, dim and var names, added support of LSC2 and mixed triangle and quad mesh
!                      used standard and long name as CF convention

!           2013-12-09 drf: NF_FLOATs, add SZ metadata, netcdf error checking towards ugrid v0.9
!              Add sigma attributes and sigma_theta_f, sigma_theta_b, sigma_h_c' variables for CF 
!              Add mesh variable for ugrid compliance
!              Parse binary file time (from bctides.in header) into udunits2 compatible units  
!   
!           changed non-standard outputs which can be viz'ed with vis6 now (vertical
!           coordinates not dependable).

!           From combine_output8: do not split quads.
! TODO: mark dry elev. according to SURA standard
subroutine combine_output9(files,nfile,ibgn,iend,inc)
!-------------------------------------------------------------------------------
!  use typeSizes
!  use netcdf
  use netcdf_var_names
  implicit real(4)(a-h,o-z),integer(i-n)
  include 'netcdf.inc'
  
  integer                              :: nfile  !< number of file base names
  character(len=30),dimension(nfile)   :: files  !< base names (e.g. elev.61) of files to combine
  integer                              :: ibgn   !< first output file index to process
  integer                              :: iend   !< last output file index to process
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

  integer invalid_index
  
  !netcdf variables
  integer :: coords_type ! lat/long or project coords
  character*48 lat_coord_standard_name
  character*48 lon_coord_standard_name
  integer ::  lat_str_len,lon_str_len ! length of coord name
  character(len=50) fname
  character(len=256) cbuffer, stdname,units,longname
  character(len=8) varname,vname,uname
  integer :: time_dims(1),ele_dims(2),x_dims(1),y_dims(1),z_dims(1),sigma_dims(1), &
            &var1d_dims(1),var2d_dims(2),var3d_dims(3), &
            &data_start_1d(1),data_start_2d(2),data_start_3d(3), &
            &data_count_1d(1),data_count_2d(2),data_count_3d(3), &
            & int_buffer(4),dummy_dim(1),ihgrid_id, tempint_array(1),&
            & chunks(3)
  real*4  :: hc_array(1),hs_array(1),thetab_array(1),thetaf_array(1),real_buffer(4) 
  
  character(len=long_name_len) netcdf_var_long_name(2) 
  character(len=long_name_len) netcdf_var_standard_name(2) 
  character(len=dataset_name_len) netcdf_var_name(2)  
  character(len=4) netcdf_out_location
  character(len=5) netcdf_level_location
  integer netcdf_var_dim
  logical found_netcdf_var
      


! Read local_to_global_0000 for global info
  open(10,file='local_to_global_0000',status='old')
  read(10,*)ne_global,np_global,nvrt,nproc !,ntracers
  close(10)
  invalid_index = -99999
  do ifile = 1,nfile
      file63=files(ifile)
      allocate(x(np_global),y(np_global),dp(np_global),kbp00(np_global),kbe(ne_global),i34(ne_global), &
              &np(0:nproc-1),ns(0:nproc-1),ne(0:nproc-1),elnode(4,ne_global),nm2(ne_global,4),eta2(np_global), &
              &ztot(nvrt),sigma(nvrt),cs(nvrt),outb(np_global,nvrt,2),ihot_len(0:nproc-1), &
              &outeb(ne_global,nvrt,2),dpe(ne_global),xctr(ne_global),yctr(ne_global), &
              &eta2e(ne_global),stat=istat)
      if(istat/=0) stop 'Allocation error: x,y'

    ! Initialize outb for ivalid pts (below bottom etc)
      outb=-9999.
      outeb=-9999.

      elnode = invalid_index
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
        &(i34(ielg(irank,m)),(nm2(m,mm),mm=1,i34(ielg(irank,m))),m=1,ne(irank))

        close(10)

        !Debug
        !write(98,*)irank,(i34(ielg(irank,m)),m=1,ne(irank))

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
          do mm=1,i34(iegb)
            itmp=nm2(m,mm)
            if(itmp>np(irank).or.itmp<=0) then
              write(*,*)'Overflow:',m,mm,itmp
              stop
            endif
            elnode(mm,iegb)=iplg(irank,itmp)
          enddo !mm

          !Debug
          !write(99,*)iegb,i34(iegb),elnode(1:i34(iegb),iegb)

        enddo !m
      enddo !irank=0,nproc-1

    ! Compute geometry
      call compute_nside(np_global,ne_global,i34,elnode(1:4,1:ne_global),ns2)
      allocate(ic3(4,ne_global),elside(4,ne_global),isdel(2,ns2),isidenode(2,ns2),xcj(ns2),ycj(ns2),stat=istat)
      if(istat/=0) stop 'Allocation error: side(0)'
      call schism_geometry(np_global,ne_global,ns2,x,y,i34,elnode(1:4,1:ne_global),ic3(1:4,1:ne_global), &
         &elside(1:4,1:ne_global),isdel,isidenode,xcj,ycj)

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
        kbe(i)=minval(kbp00(elnode(1:i34(i),i))) !maxval(kbp00(elnode(1:i34(i),i)))
        dpe(i)=sum(dp(elnode(1:i34(i),i)))/i34(i)
        xctr(i)=sum(x(elnode(1:i34(i),i)))/i34(i)
        yctr(i)=sum(y(elnode(1:i34(i),i)))/i34(i)
      enddo !i
      do i=1,ns_global
        kbs(i)=minval(kbp00(isidenode(1:2,i))) !maxval(kbp00(isidenode(1:2,i)))
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
          open(65,file=it_char(1:it_len)//'_'//file63,status='replace',form="unformatted",access="stream")
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
            write(65) (x(m),y(m),dp(m),kbp00(m),m=1,np_global)
            do m=1,ne_global
              write(65) i34(m), (elnode(mm,m),mm=1,i34(m))
            enddo !m
          else if(i23d==4.or.i23d==6.or.i23d==7) then !side based
            write(65)  ns_global
            write(65) ns_e
            write(65) (xcj(m),ycj(m),dps(m),kbs(m),m=1,ns_global)
            do m=1,ns_e
              write(65) 3,(nm_s(m,mm),mm=1,3)
            enddo !m
          else !elem. based
            write(65) ne_global
            write(65) ne_e
            write(65) (xctr(m),yctr(m),dpe(m),kbe(m),m=1,ne_global)
            do m=1,ne_e
              write(65) 3,(nm_e(m,mm),mm=1,3)
            enddo !m
          endif !i23d

        else !netcdf
        
          !PRINT *, "netcf version is ",nf_inq_libvers()
    !     enter define mode
          fname=it_char(1:it_len)//'_'//file63(1:lfile63-3)//'.nc'
#ifdef   NETCDF_4
          iret = nf_create(trim(fname), OR(NF_CLOBBER,NF_NETCDF4), ncid)
#else
          iret = nf_create(trim(fname), NF_CLOBBER, ncid)
#endif
    !     define dimensions
          iret = nf_def_dim(ncid, 'nSCHISM_hgrid_node',np_global, node_dim)
          iret = nf_def_dim(ncid, 'nSCHISM_hgrid_face',ne_global, nele_dim)
          iret = nf_def_dim(ncid, 'nSCHISM_hgrid_edge',ns_global, nedge_dim)
          iret = nf_def_dim(ncid, 'nMaxSCHISM_hgrid_face_nodes',4, nfour_dim)
          iret = nf_def_dim(ncid, 'nSCHISM_vgrid_layers',nvrt, nv_dim)
          iret = nf_def_dim(ncid, 'one',1, none_dim)
          iret = nf_def_dim(ncid, 'two',2, ntwo_dim)
          iret = nf_def_dim(ncid, 'sigma',nvrt-kz+1, nsigma_dim)
          if(kz/=1) iret = nf_def_dim(ncid, 'nz',kz-1, nz_dim)
          iret = nf_def_dim(ncid, 'time', NF_UNLIMITED, ntime_dim)
      
      
          coords_type =1 ! default project coords
      
          if (coords_type.eq.1) then
             lat_coord_standard_name = "projection_x_coordinate"
             lon_coord_standard_name = "projection_y_coordinate"
             lat_str_len = 23
             lon_str_len = 23
          else
             lat_coord_standard_name = "latitude"
             lon_coord_standard_name = "longitude"
             lat_str_len = 8
             lon_str_len = 9
          endif

      
    !     define variables
          time_dims(1) = ntime_dim
          iret=nf_def_var(ncid,'time',NF_DOUBLE,1,time_dims,itime_id)
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

     
    !     define mesh
          iret=nf_def_var(ncid,'SCHISM_hgrid',NF_INT,0,dummy_dim,ihgrid_id)
          iret=nf_put_att_text(ncid,ihgrid_id,'long_name',37,'Topology data of 2d unstructured mesh')
          tempint_array(1)=2
          iret=nf_put_att_int(ncid,ihgrid_id,'topology_dimension',NF_INT,1,tempint_array)
          iret=nf_put_att_text(ncid,ihgrid_id,'cf_role',13,'mesh_topology')
          iret=nf_put_att_text(ncid,ihgrid_id,'node_coordinates',39,'SCHISM_hgrid_node_x SCHISM_hgrid_node_y')
          iret=nf_put_att_text(ncid,ihgrid_id,'face_node_connectivity',23,'SCHISM_hgrid_face_nodes')
          iret=nf_put_att_text(ncid,ihgrid_id,'edge_coordinates',39,'SCHISM_hgrid_edge_x SCHISM_hgrid_edge_y')
          iret=nf_put_att_text(ncid,ihgrid_id,'face_coordinates',39,'SCHISM_hgrid_face_x SCHISM_hgrid_face_y')
          iret=nf_put_att_text(ncid,ihgrid_id,'edge_node_connectivity',23,'SCHISM_hgrid_edge_nodes')

          iret=nf_def_var(ncid,'Mesh3D',NF_INT,0,dummy_dim,i3Dgrid_id)
          iret=nf_put_att_text(ncid,i3Dhgrid_id,'long_name',37,'Topology data of 3d unstructured mesh')
          tempint_array(1)=3
          iret=nf_put_att_int(ncid,i3Dhgrid_id,'topology_dimension',NF_INT,1,tempint_array)
          iret=nf_put_att_text(ncid,i3Dhgrid_id,'cf_role',13,'mesh_topology')
     
          ele_dims(2)=nele_dim; ele_dims(1)=nfour_dim
          iret=nf_def_var(ncid,'SCHISM_hgrid_face_nodes',NF_INT,2,ele_dims,iele_id)
          iret=nf_put_att_text(ncid,iele_id,'long_name',35,'Horizontal Triangular Element Table')
          iret=nf_put_att_text(ncid,iele_id,'units',15,'non-dimensional')
          iret=nf_put_att_text(ncid,iele_id,'cf_role',22,'face_node_connectivity')
          int_buffer(1) = 1  ! fortran indexing starts at 1
          iret=nf_put_att_int(ncid,iele_id,'start_index',NF_INT,1,int_buffer)
          int_buffer(1) = invalid_index
          iret=nf_put_att_int(ncid,iele_id,'_FillValue',NF_INT,1,int_buffer)
      
          ele_dims(2)=nedge_dim; ele_dims(1)=ntwo_dim
          iret=nf_def_var(ncid,'SCHISM_hgrid_edge_nodes',NF_INT,2,ele_dims,iedge_id)
          iret=nf_put_att_text(ncid,iedge_id,'long_name',48,'Map every edge to the two nodes that it connects')
          iret=nf_put_att_text(ncid,iedge_id,'units',15,'non-dimensional')
          iret=nf_put_att_text(ncid,iedge_id,'cf_role',22,'edge_node_connectivity')
          int_buffer(1) = 1  ! fortran indexing starts at 1
          iret=nf_put_att_int(ncid,iedge_id,'start_index',NF_INT,1,int_buffer)
      
          x_dims(1)=node_dim
          iret=nf_def_var(ncid,'SCHISM_hgrid_node_x',NF_FLOAT,1,x_dims,ix_id)
          iret=nf_put_att_text(ncid,ix_id,'long_name',17,'node x-coordinate')
          iret=nf_put_att_text(ncid,ix_id,'standard_name',lat_str_len,lat_coord_standard_name)
          iret=nf_put_att_text(ncid,ix_id,'units',1,'m')
          iret=nf_put_att_text(ncid,ix_id,'mesh',12,'SCHISM_hgrid')
      
          iret=nf_def_var(ncid,'SCHISM_hgrid_node_y',NF_FLOAT,1,x_dims,iy_id)
          iret=nf_put_att_text(ncid,iy_id,'long_name',17,'node y-coordinate')
          iret=nf_put_att_text(ncid,iy_id,'standard_name',lon_str_len,lon_coord_standard_name)
          iret=nf_put_att_text(ncid,iy_id,'units',1,'m')
          iret=nf_put_att_text(ncid,iy_id,'mesh',12,'SCHISM_hgrid')
      
      
      
          iret=nf_def_var(ncid,'node_bottom_index',NF_INT,1,x_dims,inode_bottom_id)
          iret=nf_put_att_text(ncid,inode_bottom_id,'long_name',31,'bottom level index at each node')
          iret=nf_put_att_text(ncid,inode_bottom_id,'units',15,'non-dimensional')
          iret=nf_put_att_text(ncid,inode_bottom_id,'mesh',12,'SCHISM_hgrid')
          iret=nf_put_att_text(ncid,inode_bottom_id,'location',4,'node')
          int_buffer(1) = 1  
          iret=nf_put_att_int(ncid,inode_bottom_id,'start_index',NF_INT,1,int_buffer)
      
      
          x_dims(1)=nele_dim
          iret=nf_def_var(ncid,'SCHISM_hgrid_face_x',NF_FLOAT,1,x_dims,ix_face_id)
          iret=nf_put_att_text(ncid,ix_face_id,'long_name',28,'x_coordinate of 2D mesh face')
           iret=nf_put_att_text(ncid,ix_face_id,'standard_name',lat_str_len,lat_standard_name)
          iret=nf_put_att_text(ncid,ix_face_id,'units',1,'m')
          iret=nf_put_att_text(ncid,ix_face_id,'mesh',12,'SCHISM_hgrid')
      

          iret=nf_def_var(ncid,'SCHISM_hgrid_face_y',NF_FLOAT,1,x_dims,iy_face_id)
          iret=nf_put_att_text(ncid,iy_face_id,'long_name',28,'y_coordinate of 2D mesh face')
          iret=nf_put_att_text(ncid,iy_face_id,'standard_name',lon_str_len,lon_coord_standard_name)
          iret=nf_put_att_text(ncid,iy_face_id,'units',1,'m')
          iret=nf_put_att_text(ncid,iy_face_id,'mesh',12,'SCHISM_hgrid')
      
      
          iret=nf_def_var(ncid,'ele_bottom_index',NF_INT,1,x_dims,iele_bottom_id)
          iret=nf_put_att_text(ncid,iele_bottom_id,'long_name',34,'bottom level index at each element')
          iret=nf_put_att_text(ncid,iele_bottom_id,'units',15,'non-dimensional')
          iret=nf_put_att_text(ncid,iele_bottom_id,'mesh',12,'SCHISM_hgrid')
          iret=nf_put_att_text(ncid,iele_bottom_id,'location',4,'elem')
          int_buffer(1) = 1  
          iret=nf_put_att_int(ncid,iele_bottom_id,'start_index',NF_INT,1,int_buffer)
      
      
          x_dims(1)=nedge_dim
          iret=nf_def_var(ncid,'SCHISM_hgrid_edge_x',NF_FLOAT,1,x_dims,ix_edge_id)
          iret=nf_put_att_text(ncid,ix_edge_id,'long_name',28,'x_coordinate of 2D mesh edge')
          iret=nf_put_att_text(ncid,ix_edge_id,'standard_name',lat_str_len, lat_coord_standard_name)
          iret=nf_put_att_text(ncid,ix_edge_id,'units',1,'m')
          iret=nf_put_att_text(ncid,ix_edge_id,'mesh',12,'SCHISM_hgrid')
      

          iret=nf_def_var(ncid,'SCHISM_hgrid_edge_y',NF_FLOAT,1,x_dims,iy_edge_id)
          iret=nf_put_att_text(ncid,iy_edge_id,'long_name',28,'y_coordinate of 2D mesh edge')
          iret=nf_put_att_text(ncid,iy_edge_id,'standard_name',lon_str_len,lon_coord_standard_name)
          iret=nf_put_att_text(ncid,iy_edge_id,'units',1,'m')
          iret=nf_put_att_text(ncid,iy_edge_id,'mesh',12,'SCHISM_hgrid')
    
      
          iret=nf_def_var(ncid,'edge_bottom_index',NF_INT,1,x_dims,iedge_bottom_id)
          iret=nf_put_att_text(ncid,iedge_bottom_id,'long_name',31,'bottom level index at each edge')
          iret=nf_put_att_text(ncid,iedge_bottom_id,'units',15,'non-dimensional')
          iret=nf_put_att_text(ncid,iedge_bottom_id,'mesh',12,'SCHISM_hgrid')
          iret=nf_put_att_text(ncid,iedge_bottom_id,'location',4,'edge')
          int_buffer(1) = 1  
          iret=nf_put_att_int(ncid,iedge_bottom_id,'start_index',NF_INT,1,int_buffer)
      
      
          x_dims(1)=node_dim
          iret=nf_def_var(ncid,'depth',NF_FLOAT,1,x_dims,idepth_id)
          iret=nf_put_att_text(ncid,idepth_id,'long_name',10,'Bathymetry')
          iret=nf_put_att_text(ncid,idepth_id,'units',6,'meters')
          iret=nf_put_att_text(ncid,idepth_id,'positive',4,'down')
          iret=nf_put_att_text(ncid,idepth_id,'mesh',12,'SCHISM_hgrid')
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
          iret=nf_put_att_text(ncid,ihc_id,'units',6,'meters')
          iret=nf_put_att_text(ncid,ihc_id,'positive',4,'down')

          var1d_dims(1)=none_dim
          iret=nf_def_var(ncid,'sigma_theta_b',NF_FLOAT,1,var1d_dims,itheta_b_id)
          cbuffer='ocean_s_coordinate theta_b constant'
          iret=nf_put_att_text(ncid,itheta_b_id,'long_name',len(trim(cbuffer)),cbuffer)

          var1d_dims(1)=none_dim
          iret=nf_def_var(ncid,'sigma_theta_f',NF_FLOAT,1,var1d_dims,itheta_f_id)
          cbuffer='ocean_s_coordinate theta_f constant'
          iret=nf_put_att_text(ncid,itheta_f_id,'long_name',len(trim(cbuffer)),cbuffer)

          var1d_dims(1)=none_dim
          iret=nf_def_var(ncid,'sigma_maxdepth',NF_FLOAT,1,var1d_dims,ihs_id)
          cbuffer='ocean_s_coordinate maximum depth cutoff (mixed s over z boundary'
          iret=nf_put_att_text(ncid,ihs_id,'long_name',len(trim(cbuffer)),cbuffer)
          iret=nf_put_att_text(ncid,ihs_id,'units',6,'meters')
          iret=nf_put_att_text(ncid,ihs_id,'positive',4,'down')

       
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
      
          call get_netcdf_var_names(file63(1:lfile63-3),&
                                    &file63((lfile63-1):lfile63),&
                                    &netcdf_var_name,&
                                    &netcdf_var_long_name,&
                                    &netcdf_var_standard_name,&
                                    &netcdf_out_location,&
                                    &netcdf_level_location,&
                                    &netcdf_var_dim,&
                                    &found_netcdf_var)
          if (found_netcdf_var.eq.(.false.)) then
              print *, "unsupported data type :", file63(1:lfile63-3)
              exit
          endif
      
          iret = nf_put_att_text(ncid, NF_GLOBAL, 'data_horizontal_center', 4, netcdf_out_location)
          iret = nf_put_att_text(ncid, NF_GLOBAL, 'data_vertical_center', 4, netcdf_level_location)
       
          if (netcdf_out_location.eq.'node') then
               var2d_dims(1)=node_dim; var2d_dims(2)=ntime_dim
          else if (netcdf_out_location.eq.'elem') then
                 var2d_dims(1)=nele_dim; var2d_dims(2)=ntime_dim
          else
                 var2d_dims(1)=nedge_dim; var2d_dims(2)=ntime_dim
          endif
      
          iret=nf_def_var(ncid,"eta",NF_FLOAT,2,var2d_dims,isurface_id)
          iret=nf_put_att_text(ncid,isurface_id,'long_name',&
               & len_trim("water surface height above reference datum"),&
               & "water surface height above reference datum")
          iret=nf_put_att_text(ncid,isurface_id,'standard_name',&
               & len_trim("sea_surface_height_above_reference_ellipsoid"),&
               & "sea_surface_height_above_reference_ellipsoid")
          iret=nf_put_att_text(ncid,isurface_id,'missing_value',5,'-999.')
          iret=nf_put_att_text(ncid,isurface_id,'mesh',12,'SCHISM_hgrid')
          iret=nf_put_att_text(ncid,isurface_id,'location',4,netcdf_out_location)
      
      
          if((i23d==2).or.(i23d==10).or.(i23d==4).or.(i23d==5)) then
          
            if (netcdf_out_location.eq.'node') then
               var2d_dims(1)=node_dim; var2d_dims(2)=ntime_dim
            else if (netcdf_out_location.eq.'elem') then
                 var2d_dims(1)=nele_dim; var2d_dims(2)=ntime_dim
            else
                 var2d_dims(1)=nedge_dim; var2d_dims(2)=ntime_dim
            endif
        
            iret=nf_def_var(ncid,netcdf_var_name(1),NF_FLOAT,2,var2d_dims,iu_id)
            iret=nf_put_att_text(ncid,iu_id,'long_name',len_trim(netcdf_var_long_name(1)),netcdf_var_long_name(1))
            iret=nf_put_att_text(ncid,iu_id,'standard_name',len_trim(netcdf_var_standard_name(1)),netcdf_var_standard_name(1))
            iret=nf_put_att_text(ncid,iu_id,'missing_value',5,'-999.')
            iret=nf_put_att_text(ncid,iu_id,'mesh',12,'SCHISM_hgrid')
            iret=nf_put_att_text(ncid,iu_id,'location',4,netcdf_out_location)
            iret=nf_put_att_text(ncid,iu_id,'level',5,netcdf_level_location)
        
            if (netcdf_var_dim.eq.2) then
                 iret=nf_def_var(ncid,netcdf_var_name(2),NF_FLOAT,2,var2d_dims,iv_id)
                 iret=nf_put_att_text(ncid,iv_id,'long_name',len_trim(netcdf_var_long_name(2)),netcdf_var_long_name(2))
                 iret=nf_put_att_text(ncid,iv_id,'standard_name',len_trim(netcdf_var_standard_name(2)),netcdf_var_standard_name(2))
                 iret=nf_put_att_text(ncid,iv_id,'missing_value',5,'-999.')
                 iret=nf_put_att_text(ncid,iv_id,'mesh',12,'SCHISM_hgrid')
                 iret=nf_put_att_text(ncid,iv_id,'location',4,netcdf_out_location)  
                  iret=nf_put_att_text(ncid,iv_id,'level',5,netcdf_level_location)
            endif
        
          else !3D
            var3d_dims(2)=nv_dim
            var3d_dims(3)=ntime_dim
            chunks(2)=nvrt
            chunks(3)=1
            if (netcdf_out_location.eq.'node') then
               var3d_dims(1)=node_dim
               chunks(1)=np_global
            else if (netcdf_out_location.eq.'elem') then
               var3d_dims(1)=nele_dim
               chunks(1)=ne_global  
            else        
               var3d_dims(1)=nedge_dim
               chunks(1)=ns_global
            endif
        
            iret=nf_def_var(ncid,netcdf_var_name(1),NF_FLOAT,3,var3d_dims,iu_id)
#ifdef   NETCDF_4
            iret=nf_def_var_chunking(ncid, iu_id, 0, chunks)
            iret=nf_def_var_deflate(ncid,iu_id,0,1,4)
#endif
            iret=nf_put_att_text(ncid,iu_id,'long_name',len_trim(netcdf_var_long_name(1)),netcdf_var_long_name(1))
            iret=nf_put_att_text(ncid,iu_id,'standard_name',len_trim(netcdf_var_standard_name(1)),netcdf_var_standard_name(1))
            iret=nf_put_att_text(ncid,iu_id,'missing_value',5,'-999.')
            iret=nf_put_att_text(ncid,iu_id,'mesh',6,'Mesh3D')
            iret=nf_put_att_text(ncid,iu_id,'location',4,netcdf_out_location)
            iret=nf_put_att_text(ncid,iu_id,'level',5,netcdf_level_location)
        
            if (netcdf_var_dim.eq.2) then
                 iret=nf_def_var(ncid,netcdf_var_name(2),NF_FLOAT,3,var3d_dims,iv_id)
#ifdef   NETCDF_4
                 iret=nf_def_var_chunking(ncid, iv_id, 0, chunks)
                 iret=nf_def_var_deflate(ncid,iv_id,0,1,4)
#endif
                 iret=nf_put_att_text(ncid,iv_id,'long_name',len_trim(netcdf_var_long_name(2)),netcdf_var_long_name(2))
                 iret=nf_put_att_text(ncid,iv_id,'standard_name',len_trim(netcdf_var_standard_name(2)),netcdf_var_standard_name(2))
                 iret=nf_put_att_text(ncid,iv_id,'missing_value',5,'-999.')
                 iret=nf_put_att_text(ncid,iv_id,'mesh',6,'Mesh3D')
                 iret=nf_put_att_text(ncid,iv_id,'location',4,netcdf_out_location)    
                 iret=nf_put_att_text(ncid,iv_id,'level',5,netcdf_level_location)
            endif
        
    
          endif !i23d

          iret = nf_put_att_text(ncid, NF_GLOBAL, 'Conventions', 6, 'CF-1.0')
    !     leave define mode
          iret = nf_enddef(ncid)

    !     Write mode (header part only)
          data_start_2d(1:2)=1
          data_count_2d(1)=4; data_count_2d(2)=ne_global
          iret=nf_put_vara_int(ncid,iele_id,data_start_2d,data_count_2d,elnode)
      
      
          if(iret.ne.NF_NOERR) then
             print*, nf_strerror(iret); stop
          endif
      
          data_start_2d(1:2)=1
          data_count_2d(1)=2; data_count_2d(2)=ns_global
          iret=nf_put_vara_int(ncid,iedge_id,data_start_2d,data_count_2d,isidenode)
      
          if(iret.ne.NF_NOERR) then
             print*, nf_strerror(iret); stop
          endif
 
          ! fill bottom index
          iret=nf_put_vara_int(ncid,inode_bottom_id,1,np_global,kbp00)
          iret=nf_put_vara_int(ncid,iedge_bottom_id,1,ns_global,kbs)
          iret=nf_put_vara_int(ncid,iele_bottom_id,1,ne_global,kbe)
      
      
          ! fill node,side center coords
          iret=nf_put_vara_real(ncid,ix_id,1,np_global,x)
          iret=nf_put_vara_real(ncid,iy_id,1,np_global,y)
          iret=nf_put_vara_real(ncid,ix_edge_id,1,ns_global,xcj)
          iret=nf_put_vara_real(ncid,iy_edge_id,1,ns_global,ycj)
          
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
         &(((outb(iplg(irank,i),k,m),m=1,ivs),k=max0(1,kbp01(irank,i)),nvrt),i=1,np(irank))
            else if(i23d==4) then !2D side 
              read(63,rec=ispool)time,it,itmp,(eta2(iplg(irank,i)),i=1,np(irank)), &
         &(outsb(islg(irank,i),1,1:ivs),i=1,ns(irank))
            else if(i23d==5) then !2D elem
              read(63,rec=ispool)time,it,itmp,(eta2(iplg(irank,i)),i=1,np(irank)), &
         &(outeb(ielg(irank,i),1,1:ivs),i=1,ne(irank))
            else if(i23d==6.or.i23d==7) then !3D side and whole/half level
              read(63,rec=ispool)time,it,itmp,(eta2(iplg(irank,i)),i=1,np(irank)), &
         &(((outsb(islg(irank,i),k,m),m=1,ivs),k=1,nvrt),i=1,ns(irank))
            else if(i23d==8.or.i23d==9) then !3D element at whole/half levels
              read(63,rec=ispool)time,it,itmp,(eta2(iplg(irank,i)),i=1,np(irank)), &
         &(((outeb(ielg(irank,i),k,m),m=1,ivs),k=1,nvrt),i=1,ne(irank))
            else if(i23d==10) then !2D node (n.s.)
              read(63,rec=ispool)time,it,itmp,(eta2(iplg(irank,i)),i=1,np(irank)), &
         &((outb(iplg(irank,i),1,m),m=1,ivs),i=1,np(irank))
            else !3D node
              read(63,rec=ispool)time,it,itmp,(eta2(iplg(irank,i)),i=1,np(irank)), &
         &(((outb(iplg(irank,i),k,m),m=1,ivs),k=1,nvrt),i=1,np(irank))
            endif
            ! Close input file
            close(63)
          enddo !irank

          !Output
          !print*, 'it=',it,' time=',time/86400.0    
           !Compute eta2 at sides/elem.
            do i=1,ns_global
              n1=isidenode(1,i); n2=isidenode(2,i)
              eta2s(i)=(eta2(n1)+eta2(n2))/2.d0
            enddo !i
            do i=1,ne_global
              eta2e(i)=sum(eta2(elnode(1:i34(i),i)))/i34(i)
            enddo !i

          if(inc==0) then !binary
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
            data_start_2d(1)=1; data_start_2d(2)=ispool
            data_count_2d(2)=1
            if (netcdf_out_location.eq.'node') then
                  data_count_2d(1)=np_global
                  iret=nf_put_vara_real(ncid,isurface_id,data_start_2d,data_count_2d,eta2) 
            else if (netcdf_out_location.eq.'elem') then
                  data_count_2d(1)=ne_global;
                  iret=nf_put_vara_real(ncid,isurface_id,data_start_2d,data_count_2d,eta2e) 
            else
                  data_count_2d(1)=ns_global;
                  iret=nf_put_vara_real(ncid,isurface_id,data_start_2d,data_count_2d,eta2s) 
            endif
       
          
            if((i23d==2).or.(i23d==10).or.(i23d==4).or.(i23d==5)) then
         
             iret=nf_put_vara_real(ncid,iu_id,data_start_2d,data_count_2d,outb(:,1,1))
             if (netcdf_var_dim.eq.2) then
                 iret=nf_put_vara_real(ncid,iv_id,data_start_2d,data_count_2d,outb(:,1,2))
             endif
      
            else !3D
              data_start_3d(1:2)=1; data_start_3d(3)=ispool
              data_count_3d(2)=nvrt; data_count_3d(3)=1
              if (netcdf_out_location.eq.'node') then
                  data_count_3d(1)=np_global
                  iret=nf_put_vara_real(ncid,iu_id,data_start_3d,data_count_3d,outb(:,:,1))
                  if (netcdf_var_dim.eq.2) then
                    iret=nf_put_vara_real(ncid,iv_id,data_start_3d,data_count_3d,outb(:,:,2))
                  endif
              else if (netcdf_out_location.eq.'elem') then
                  data_count_3d(1)=ne_global;
                  iret=nf_put_vara_real(ncid,iu_id,data_start_3d,data_count_3d,outeb(:,:,1))
                  if (netcdf_var_dim.eq.2) then
                    iret=nf_put_vara_real(ncid,iv_id,data_start_3d,data_count_3d,outeb(:,:,2))
                  endif
              else
                  data_count_3d(1)=ns_global;
                  iret=nf_put_vara_real(ncid,iu_id,data_start_3d,data_count_3d,outsb(:,:,1))
                  if (netcdf_var_dim.eq.2) then
                    iret=nf_put_vara_real(ncid,iv_id,data_start_3d,data_count_3d,outsb(:,:,2))
                  endif
              endif
        
            endif !i23d
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

  end do ! loop over file names 
  stop
end subroutine combine_output9



subroutine combine_output9_input(files,nfile,ibgn,iend,inetcdf)
use argparse
implicit none
integer, parameter                   :: nfilemax = 20
integer                              :: nfile   !< number of file base names
character(len=30),dimension(nfilemax) :: files   !< base names (e.g. elev.61) of files to combine
integer                              :: ibgn    !< first output file index to process
integer                              :: iend    !< last output file index to process
integer                              :: inetcdf !< netcdf flag
character(len=80) :: infile
character(len=30) :: cfile = ""

! local
integer :: comcount
integer :: ifile

infile = ""
cfile = ""
files=""

cmd_name = "combine_output9"
call cla_init(cmd_name,"Combine time blocked per-processor binary outputs (e.g. '1_000*_elev.61') into time blocked global outputs ('1_elev.61')")

call cla_register('-i','--in', 'input file (e.g. combine_input.in) containing options (overridden by command line specs)', cla_char,'')
call cla_register('-b','--begin', 'start day', cla_int,'-1')
call cla_register('-e','--end','end day', cla_int  ,'-1')
call cla_register('-n','--nc','combine to NetCDF format (1) or ordinary binary (0)', cla_int  ,'0')
call cla_register('-f','--file','base file name like elev.61',  cla_char,'') 
call cla_validate


    
comcount = command_argument_count()
if (comcount == 0) then
  infile = "combine_output.in"
else
  call cla_get("--in",infile)
end if



! possibly read files from input file if there is one
if (len_trim(infile)>1) then
  ! read inputs from combine_output.in or other provided file
  open(10,file=infile,status='old')
  read(10,'(a12)') cfile !e.g. 'hvel.64'
  read(10,*) ibgn,iend  
  read(10,*) inetcdf       ! netcdf option
  close(10)
! command line takes precedence
  if (nfile == 0) then
    nfile = 1
    files(1) = trim(cfile)
  end if 
  if (cla_key_present("--begin"))call cla_get("--begin",ibgn)
  if (cla_key_present("--end"))  call cla_get("--end",iend)
  if (cla_key_present("--nc"))   call cla_get("--nc",inetcdf)
  if (cla_key_present("--file")) call cla_get("--file",cfile)

else
  call cla_get("--begin",ibgn)
  call cla_get("--end",iend)
  call cla_get("--nc",inetcdf)
  call cla_get("--file", cfile)
end if

if (len_trim(cfile) > 1) then
  nfile = 1
  files(1) = cfile
else
  print*, "No base files specified."
end if

! validate and prioritize
if (ibgn > iend) then
  stop("Beginning index after end. Check input file or -b and -e options")
end if
print '("Begin: ",i4," End: ",i4," NetCDF enabled: ", i4)',ibgn,iend,inetcdf
print*,"File:"
do ifile = 1,nfile
  print*, files(ifile)
end do

return
end subroutine


!===============================================================================
program combine_output
implicit none
integer, parameter                   :: nfilemax = 20
integer                              :: nfile  !< number of file base names
character(len=30),dimension(nfilemax) :: files  !< base names (e.g. elev.61) of files to combine
integer                              :: ibgn   !< first output file index to process
integer                              :: iend   !< last output file index to process
integer                              :: inetcdf !< netcdf flag
call combine_output9_input(files,nfile,ibgn,iend,inetcdf)
call combine_output9(files,nfile,ibgn,iend,inetcdf)
end program



