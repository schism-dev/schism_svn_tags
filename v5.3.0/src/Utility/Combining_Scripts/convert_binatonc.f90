!===============================================================================
! Read in binary outputs (*_elev.61) and convert them into
! netcdf format(*_elev.nc). 
! Run this program inside the directory outputs/, where some of the input files below
! can be found.

! Inputs:
!        *_elev.61 files (from SCHISM outputs/); 
!        combine_output.in (1st line: elev.61 etc; 2nd line: start and end file indices;
!                          3rd line: inc=1 (netcdf option)); 
! Output: netcdf format(*_elev.nc)
!
!  ifort -Bstatic -O2 -assume byterecl -o combine_output8 combine_output8.f90 ...../schism_geometry.f90 -Vaxlib -I/share/apps/netcdf/include/ -L/share/apps/netcdf/lib/ -lnetcdf

!  On TYPHOON:
!  pgf90 -O2 -mcmodel=medium  -Bstatic -o convert_binatonc convert_binatonc.f90 -Mbounds -I/usr/local/amd64a/opteron/pgi/netcdf-3.6.2/include -L/usr/local/amd64a/opteron/pgi/netcdf-3.6.2/lib -lnetcdf

!  History: 
!           2016-01-19 Hai Huang: from combined binary to nc
!           2016-01-19 drf: 
!===============================================================================

program convert_binatonc
!-------------------------------------------------------------------------------
!  use typeSizes
!  use netcdf

!-------------------------------------------------------------read_output7_allnodes.f90
      parameter(nbyte=4)
      character*30 file63
      character*12 it_char
      character*48 start_time,version,variable_nm,variable_dim
      character*48 data_format
  
      integer,allocatable :: i34(:),elnode(:,:)
      allocatable :: sigma(:),cs(:),ztot(:),x(:),y(:),dp(:),kbp00(:),kfp(:)
      allocatable :: out(:,:,:),icum(:,:,:),eta2(:),ztmp(:,:)
      allocatable :: xctr(:),yctr(:),dpe(:),kbe(:)
      allocatable :: zs(:,:),ze(:,:),idry(:),outs(:,:,:),oute(:,:,:)
      allocatable :: kbp(:),sigma_lcl(:,:)
!-------------------------------------------------------------read_output7_allnodes.f90

!  implicit real(4)(a-h,o-z),integer(i-n)
  include 'netcdf.inc'
  allocatable outb(:,:,:)

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

!-------------------------------------------------------------read_output7_allnodes.f90
  open(10,file='combine_output.in',status='old') 
  read(10,'(a30)') file63 !e.g. 'hvel.64'
  read(10,*) ibgn,iend
  read(10,*) inc !inc=1: netcdf option
  close(10)
  file63=adjustl(file63)
  lfile63=len_trim(file63)

!...  Header
! Loop over input files
  do iinput=ibgn,iend 
      write(it_char,'(i12)')iinput 
    it_char=adjustl(it_char)  !place blanks at end
    it_len=len_trim(it_char)  !length without trailing blanks
!      open(63,file=it_char(1:it_len)//'_'//'elev.61',status='old',access='direct',recl=nbyte)
      open(63,file=it_char(1:it_len)//'_'//file63(1:lfile63),status='old',access='direct',recl=nbyte)
!      open(63,file='31_elev.61',status='old',access='direct',recl=nbyte)
      irec=0
      do m=1,48/nbyte 
        read(63,rec=irec+m) data_format(nbyte*(m-1)+1:nbyte*m) 
      enddo
      if(data_format.ne.'DataFormat v5.0') then
        print*, 'This code reads only v5.0:  ',data_format
        stop
      endif
      irec=irec+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) version(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) start_time(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) variable_nm(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) variable_dim(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte

      write(*,'(a48)')data_format
      write(*,'(a48)')version
      write(*,'(a48)')start_time
      write(*,'(a48)')variable_nm
      write(*,'(a48)')variable_dim

      read(63,rec=irec+1) nrec 
      read(63,rec=irec+2) dtout
      read(63,rec=irec+3) nspool
      read(63,rec=irec+4) ivs
      read(63,rec=irec+5) i23d 
      irec=irec+5

      print*, 'i23d=',i23d,' nrec= ',nrec

!     Vertical grid (obsolete)
      read(63,rec=irec+1) nvrt 
      read(63,rec=irec+2) kz
      read(63,rec=irec+3) h0   
      read(63,rec=irec+4) h_s
      read(63,rec=irec+5) h_c
      read(63,rec=irec+6) theta_b
      read(63,rec=irec+7) theta_f
      irec=irec+7
      allocate(ztot(nvrt),sigma(nvrt),cs(nvrt),stat=istat)
      if(istat/=0) stop 'falied to allocate (1)'

      do k=1,kz-1
        read(63,rec=irec+k) ztot(k)
      enddo
      do k=kz,nvrt
        kin=k-kz+1
        read(63,rec=irec+k) sigma(kin)
        cs(kin)=(1-theta_b)*sinh(theta_f*sigma(kin))/sinh(theta_f)+ &
     &theta_b*(tanh(theta_f*(sigma(kin)+0.5))-tanh(theta_f*0.5))/2/tanh(theta_f*0.5)
      enddo

!     Horizontal grid
      irec=irec+nvrt 
      read(63,rec=irec+1) np 
      read(63,rec=irec+2) ne 
      irec=irec+2
      allocate(x(np),y(np),dp(np),kbp00(np),kfp(np),i34(ne),elnode(4,ne),out(np,nvrt,2),idry(np), &
     !&icum(np,nvrt,2),eta2(np),stat=istat)
     &eta2(np),xctr(ne),yctr(ne),dpe(ne),kbe(ne),ztmp(nvrt,np),ze(nvrt,ne), &
     &oute(2,nvrt,ne),stat=istat)
      if(istat/=0) stop 'Falied to allocate (2)'

      do m=1,np 
        read(63,rec=irec+1)x(m)
        read(63,rec=irec+2)y(m)
        read(63,rec=irec+3)dp(m)
        read(63,rec=irec+4)kbp00(m) 
        irec=irec+4
      enddo !m=1,np

      do m=1,ne
        read(63,rec=irec+1)i34(m) 
        irec=irec+1
        do mm=1,i34(m)
          read(63,rec=irec+1)elnode(mm,m) 
          irec=irec+1
        enddo !mm
      enddo !m
      irec0=irec


      print*, 'last element',elnode(1:i34(ne),ne)

!--------------------------------------------------------------Write header to output files 

!!     enter define mode
      fname=it_char(1:it_len)//'_'//file63(1:lfile63-3)//'.nc' 
!      fname=it_char(1:it_len)//'_'//'elev'//'.nc' 
!      fname='31_elev.nc' 
      iret = nf_create(trim(fname), NF_CLOBBER, ncid)
!     define dimensions
      iret = nf_def_dim(ncid, 'node',np, node_dim) 
      iret = nf_def_dim(ncid, 'nele',ne, nele_dim)
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

      if(i23d==2) then
        var2d_dims(1)=node_dim; var2d_dims(2)=ntime_dim
        if(file63((lfile63-1):lfile63).eq.'62') then
          iret=nf_def_var(ncid,'dahv_u',NF_FLOAT,2,var2d_dims,iu_id)
          iret=nf_put_att_text(ncid,iu_id,'long_name',27,'Eastward Velocity or stress')
          iret=nf_put_att_text(ncid,iu_id,'missing_value',4,'-99.')
          iret=nf_put_att_text(ncid,iu_id,'mesh',4,'mesh')
          iret=nf_put_att_text(ncid,iu_id,'location',4,'node')
          iret=nf_def_var(ncid,'dahv_v',NF_FLOAT,2,var2d_dims,iv_id)
          iret=nf_put_att_text(ncid,iv_id,'long_name',28,'Northward Velocity or stress')
          iret=nf_put_att_text(ncid,iv_id,'missing_value',4,'-99.')
          iret=nf_put_att_text(ncid,iv_id,'mesh',4,'mesh')
          iret=nf_put_att_text(ncid,iv_id,'location',4,'node')
        else
!          iret=nf_def_var(ncid,file63(1:lfile63-3),NF_FLOAT,2,var2d_dims,ivar2d_id)
          iret=nf_def_var(ncid,'elev',NF_FLOAT,2,var2d_dims,ivar2d_id)
          iret=nf_put_att_text(ncid,ivar2d_id,'long_name',lfile63,file63(1:lfile63)) 
          iret=nf_put_att_text(ncid,ivar2d_id,'mesh',4,'mesh')
          iret=nf_put_att_text(ncid,ivar2d_id,'location',4,'node')
        endif !file63
      else !3D
        var3d_dims(1)=node_dim; var3d_dims(2)=nv_dim; var3d_dims(3)=ntime_dim
        if(file63((lfile63-1):lfile63).eq.'64') then
          iret=nf_def_var(ncid,'u',NF_FLOAT,3,var3d_dims,iu_id)
          iret=nf_put_att_text(ncid,iu_id,'long_name',23,'Eastward Water Velocity')
          iret=nf_put_att_text(ncid,iu_id,'missing_value',4,'-99.')
          iret=nf_put_att_text(ncid,iu_id,'mesh',4,'mesh')
          iret=nf_put_att_text(ncid,iu_id,'location',4,'node')
          iret=nf_def_var(ncid,'v',NF_FLOAT,3,var3d_dims,iv_id)
          iret=nf_put_att_text(ncid,iv_id,'long_name',24,'Northward Water Velocity')
          iret=nf_put_att_text(ncid,iv_id,'missing_value',4,'-99.')
          iret=nf_put_att_text(ncid,iv_id,'mesh',4,'mesh')
          iret=nf_put_att_text(ncid,iv_id,'location',4,'node')
        else
          iret=nf_def_var(ncid,file63(1:lfile63-3),NF_FLOAT,3,var3d_dims,ivar3d_id)
          iret=nf_put_att_text(ncid,ivar3d_id,'long_name',lfile63,file63(1:lfile63))
          iret=nf_put_att_text(ncid,ivar3d_id,'missing_value',4,'-99.')
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

      iret = nf_put_att_text(ncid, imesh_id, 'Conventions', 6, 'CF-1.0')
!     leave define mode
      iret = nf_enddef(ncid) 

!     Write mode (header part only)
      data_start_2d(1:2)=1
      data_count_2d(1)=3; data_count_2d(2)=ne
      iret=nf_put_vara_int(ncid,iele_id,data_start_2d,data_count_2d,elnode(1:3,1:ne))
      if(iret.ne.NF_NOERR) then
         print*, nf_strerror(iret); stop
      endif
      iret=nf_put_vara_real(ncid,ix_id,1,np,x)
      iret=nf_put_vara_real(ncid,iy_id,1,np,y)
      iret=nf_put_vara_real(ncid,idepth_id,1,np,dp)
      iret=nf_put_vara_real(ncid,isigma_id,1,nvrt-kz+1,sigma)
      iret=nf_put_vara_real(ncid,ics_id,1,nvrt-kz+1,cs)
      iret=nf_put_vara_real(ncid,ihc_id,1,1,h_c)
      iret=nf_put_vara_real(ncid,itheta_b_id,1,1,theta_b)
      iret=nf_put_vara_real(ncid,itheta_f_id,1,1,theta_f)
      iret=nf_put_vara_real(ncid,ihs_id,1,1,h_s)
      if(kz/=1) iret=nf_put_vara_real(ncid,iz_id,1,kz-1,ztot)
!--------------------------------------------------------------Write header to output files 

!...  Time iteration
!...
      out=-99 !init.
      ztmp=-99

      irec=irec0 
      do ispool=1,nrec 
        read(63,rec=irec+1) time
        read(63,rec=irec+2) it
        irec=irec+2

        print*, 'time in days =',time/86400

        do i=1,np
          read(63,rec=irec+i) eta2(i) 
        enddo !i
        irec=irec+np


          do i=1,np 
            do m=1,ivs 
              read(63,rec=irec+1) out(i,1,m)
              irec=irec+1
            enddo !m
          enddo !i


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!-------------------------------------------------------------read_output7_allnodes.f90

      
!-------------------------------------------------------------------------------combine_output8.f90




    ! Loop over output spools in file


      !Output
      !print*, 'it=',it,' time=',time/86400.0
        iret=nf_put_vara_real(ncid,itime_id,ispool,1,time)
        if(i23d==2) then
          data_start_2d(1)=1; data_start_2d(2)=ispool
          data_count_2d(1)=np; data_count_2d(2)=1
          if(file63((lfile63-1):lfile63).eq.'62') then
            iret=nf_put_vara_real(ncid,iu_id,data_start_2d,data_count_2d,outb(:,1,1))
            iret=nf_put_vara_real(ncid,iv_id,data_start_2d,data_count_2d,outb(:,1,2))
          else
            iret=nf_put_vara_real(ncid,ivar2d_id,data_start_2d,data_count_2d,out(:,1,1))
          endif
        else !3D
          data_start_3d(1:2)=1; data_start_3d(3)=ispool
          data_count_3d(1)=np; data_count_3d(2)=nvrt; data_count_3d(3)=1
          if(file63((lfile63-1):lfile63).eq.'64') then
            iret=nf_put_vara_real(ncid,iu_id,data_start_3d,data_count_3d,outb(:,:,1))
            iret=nf_put_vara_real(ncid,iv_id,data_start_3d,data_count_3d,outb(:,:,2))
          else
            iret=nf_put_vara_real(ncid,ivar3d_id,data_start_3d,data_count_3d,outb(:,:,1))
          endif
        endif !i23d
        print*, 'done record # ',ispool

    enddo !ispool=1,nrec
    ! Close output file
       iret = nf_close(ncid)

  deallocate(ztot,sigma,cs,x,y,dp,kbp00,kfp,i34,elnode,out,idry, &
     !&icum(np,nvrt,2),eta2(np),stat=istat)
     &eta2,xctr,yctr,dpe,kbe,ztmp,ze, &
     &oute)

  enddo ! do iinput=ibgn,iend
  stop
end program convert_binatonc
