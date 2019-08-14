! Read NRL/NCOM netcdf files and output in SELFE binary output format (v5.0).
! Sub-sample original NCOM grid.
! Currently hvel.64 has only surface values for all levels, and elevations
! and vertical vel. are 0.
!
! pturner 12-2004
!
! X = longitude dimension is sampled in the salt and temp data files.
! Y = latitude dimension is sampled in the salt and temp data files.
! The level dimensions goes from 1 = surface, 40 = bottom.
!
! Need to have the following modules symlinked into the current directory
!
! NETCDF.mod -> /usr/local/netcdf/include/NETCDF.mod
! TYPESIZES.mod -> /usr/local/netcdf/include/TYPESIZES.mod
!

!zyl
!   Input: 
!     (1) vgrid.in (SELFE or ELCIRC test02k4);
!     (2) date.in: year, month, day, # of days; 2nd line:  slam0, sfea0 (center of CPP projection used
!                  for output)
!     (3) screen input: itype (1: elev.61; 2: hvel.64; 3: vert.63)
!   Output: 
!     (1) binary.out (1_hvel.64, or 1_elev.61, or 1_vert.63 depending on itype input);
!     (2) hgrid.gr3 (sub-domain of original NCOM domain);
!     (3) hgrid.ll;

      program readNCOM
!     netcdf modules from symlinks
      use typeSizes
      use netcdf

!      implicit none

!      integer, parameter :: debug=1
      integer, parameter :: mnp=100000
      integer, parameter :: mne=200000
      integer, parameter :: mns=130000
      integer, parameter :: mnv=10
      integer, parameter :: mnei=20 !neighbor
      integer, parameter :: nbyte=4
      real, parameter :: small1=1.e-3 !used to check area ratios
  
!     netcdf related variables
      integer :: hid, latid, lonid, zmid, sid, tid ! Netcdf file IDs
      integer :: latvid, lonvid, zmvid, hvid ! positional variables
      integer :: xdid, xvid ! longitude index
      integer :: ydid, yvid ! latitude index
      integer :: ldid, lvid ! vertical level, 1 is top
      integer :: svid, tvid ! uu2 & vv2 variable IDs
      integer, dimension(nf90_max_var_dims) :: dids
             
!     Local variables for data
      integer, dimension(:), allocatable :: xind, yind, lind ! sampled indexes
!     Lat, lon, bathymetry
      real (kind = FourByteReal), dimension(:,:), allocatable :: lat, lon, h
!     Vertical postion
      real (kind = FourByteReal), dimension(:,:,:), allocatable :: zm, uu2, vv2
      integer, dimension(:,:), allocatable :: kbp,ihope
!     File names for netcdf files
      character(len=1024) :: hnc, latnc, lonnc, zmnc, sfile, tfile
!     Command line arguments
      character(len=1024) :: s, yr, md
      character(len=4) :: iyear_char
      character(len=1) :: char1,char2
      character(len=2) :: char3,char4
!     external function for number of command line arguments
      integer :: iargc

      integer :: status ! netcdf local status variable
      integer :: ier ! allocate error return.
      integer :: ixlen, iylen, ilen ! sampled lengths in each coordinate direction
      integer :: ixlen1, iylen1,ixlen2, iylen2 !reduced indices for CORIE grid to speed up interpolation
!     integer :: i,j,k,i1,i2,i3,j1,j2,j3
      character(len=1) :: imodel

      character(len=48) :: start_time,version,data_format
      character(len=48) :: variable_nm,variable_dim
      integer :: elnode(4,mne)
      dimension xout(mnp),yout(mnp),dp(mnp),i34(mne),kbp00(mnp)
      dimension xl(mnp),yl(mnp)
      dimension ztot(0:mnv),sigma(mnv),cs(mnv),z(mnp,mnv)
      dimension vel(mnp,mnv,2)
      dimension wild(100),wild2(100,2),month_day(12)

!     First statement
      pi=3.1415926
      open(10,file='date.in',status='old')
      read(10,*) iyear,month0,iday0,ndays
      read(10,*) slam0,sfea0
      slam0=slam0*pi/180
      sfea0=sfea0*pi/180
      close(10)

      print*, 'Which file do you want: (1: elev.61; 2: hvel.64; 3: vert.63)'
!'
      read(*,*) itype

      open(63,file='binary.out',access='direct',recl=nbyte)

!     Read in vgrid
      open(19,file='vgrid.in',status='old')
      read(19,*) nvrt,kz,h_s !kz>=1
      if(nvrt>mnv.or.nvrt<3) then
        write(11,*)'nvrt > mnv or nvrt<4'
        stop
      endif
      if(kz<1.or.kz>nvrt-2) then
        write(11,*)'Wrong kz:',kz
        stop
      endif
      if(h_s<10) then
        write(11,*)'h_s needs to be larger:',h_s
        stop
      endif

!     # of z-levels excluding "bottom" at h_s
      read(19,*) !for adding comment "Z levels"
      do k=1,kz-1
        read(19,*)j,ztot(k)
        if(k>1.and.ztot(k)<=ztot(k-1).or.ztot(k)>=-h_s) then
          write(11,*)'z-level inverted:',k
          stop
        endif
      enddo !k
      read(19,*) !level kz       
!     In case kz=1, there is only 1 ztot(1)=-h_s
      ztot(kz)=-h_s

      nsig=nvrt-kz+1 !# of S levels (including "bottom" & f.s.)
      read(19,*) !for adding comment "S levels"
      read(19,*)h_c,theta_b,theta_f
      if(h_c<5) then !large h_c to avoid 2nd type abnormaty
        write(11,*)'h_c needs to be larger:',h_c
        stop
      endif
      if(theta_b<0.or.theta_b>1) then
        write(11,*)'Wrong theta_b:',theta_b
        stop
      endif
      if(theta_f<=0) then 
        write(11,*)'Wrong theta_f:',theta_f 
        stop
      endif
!     Pre-compute constants
      s_con1=sinh(theta_f)

      sigma(1)=-1 !bottom
      sigma(nsig)=0 !surface
      read(19,*) !level kz
      do k=kz+1,nvrt-1
        kin=k-kz+1
        read(19,*) j,sigma(kin)
        if(sigma(kin)<=sigma(kin-1).or.sigma(kin)>=0) then
          write(11,*)'Check sigma levels at:',k,sigma(kin),sigma(kin-1)
          stop
        endif
      enddo
      read(19,*) !level nvrt
      close(19)

!     Compute C(s) and C'(s)
      do k=1,nsig
        cs(k)=(1-theta_b)*sinh(theta_f*sigma(k))/sinh(theta_f)+ &
     &theta_b*(tanh(theta_f*(sigma(k)+0.5))-tanh(theta_f*0.5))/2/tanh(theta_f*0.5)
      enddo !k=1,nvrt

!     Read NCOM results
!     A check to netcdf library to make sure all is OK. Not needed 
!     but left here anyway.
      if(.not. byteSizesOK()) then
        print *, "Compiler does not appear to support required kinds of variables."
!"
        stop
      end if

!
! geometry files - should never change.
!
      hnc = '/home/workspace/ccalmr6/nrldata/model_h.nc'
      latnc = '/home/workspace/ccalmr6/nrldata/model_lat.nc'
      lonnc = '/home/workspace/ccalmr6/nrldata/model_lon.nc'
      zmnc = '/home/workspace/ccalmr6/nrldata/model_zm.nc'
  
!      if (iargc() /= 4) then
!        call getarg(0,s)
!        write(*,*) 'Usage: ', trim(s), ' [YYYY] [Month Day] [# of days]' 
!        write(*,*) 'Example: ', trim(s), ' 2004 2 4 7' 
!        stop
!      end if

!!    Get year and date from the command line
!     call getarg(1, yr)
!     call getarg(2, month_char)
!     call getarg(3, day_char)
!     call getarg(4, ndays_char)

      write(iyear_char,'(i4)') iyear
      yr=iyear_char
      month_day(1:7:2)=31
      month_day(4:6:2)=30
      month_day(8:12:2)=31
      month_day(9:11:2)=30
      if(mod(iyear,4)==0) then
        month_day(2)=29
      else
        month_day(2)=28
      endif

!
!     Loop over all days
!
      month=month0
      iday=iday0-1
      do ifile=0,ndays
!-----------------------------------------------------------------------------------------------------------

      iday=iday+1
      if(iday>month_day(month)) then
        iday=iday-month_day(month)
        month=month+1
        if(month>12) then
          write(11,*)'Beyond one year'
          stop
        endif
      endif

      if(month<10) then
        write(char1,'(i1)') month
        char3='0'//char1
      else
        write(char3,'(i2)') month
      endif
      if(iday<10) then
        write(char2,'(i1)') iday
        char4='0'//char2
      else
        write(char4,'(i2)') iday 
      endif
      md=char3//char4

      write(*,*)'doing '//iyear_char//' '//trim(md)

!     Create u3d and v3d file names for this date.
!     Compile spits out a warning here, dunno how to fix.
!     Read in from symlinks in ccalmr (to fill gaps manually) 
      sfile = '/home/workspace/ccalmr6/nrldata/'//trim(yr)//'/u3d/u3d.glb8_2f_' // trim(yr) // trim(md) // '00.nc'
      tfile = '/home/workspace/ccalmr6/nrldata/'//trim(yr)//'/v3d/v3d.glb8_2f_' // trim(yr) // trim(md) // '00.nc'

!     Open files 
      status = nf90_open(trim(sfile), nf90_nowrite, sid)
      call check(status)
      status = nf90_open(trim(tfile), nf90_nowrite, tid)
      call check(status)

!     Get index information for sampled grid, assumed the same for vv2.
!     These will index into the static lat, lon, bathymetry, and zm files.
      status = nf90_inq_dimid(sid, "X_Index", xdid)
      call check(status)
      status = nf90_inq_dimid(sid, "Y_Index", ydid)
      call check(status)
      status = nf90_inq_dimid(sid, "level", ldid)
      call check(status)

      status = nf90_inq_varid(sid, "X_Index", xvid)
      call check(status)
      status = nf90_inq_varid(sid, "Y_Index", yvid)
      call check(status)
      status = nf90_inq_varid(sid, "level", lvid)
      call check(status)

!     get the variable ids 
      status = nf90_inq_varid(sid, "U_Velocity", svid)
      call check(status)
      status = nf90_inq_varid(tid, "V_Velocity", tvid)
      call check(status)

!     Get the lengths of the dimensions 
      status = nf90_Inquire_Variable(sid, svid, dimids = dids)
      call check(status)
      status = nf90_Inquire_Dimension(sid, dids(1), len = ixlen)
      call check(status)
      status = nf90_Inquire_Dimension(sid, dids(2), len = iylen)
      call check(status)
      status = nf90_Inquire_Dimension(sid, dids(3), len = ilen)
      call check(status)

!     allocate memory.
      if(ifile==0) then
        allocate( xind(ixlen),stat=ier)
        allocate( yind(iylen),stat=ier)
        allocate( lind(ilen),stat=ier)
        allocate( uu2(ixlen, iylen, ilen),stat=ier)
        if(ier /= 0) then
          write(*,*) ' Could not allocate uu2'
          stop 'Allocate uu2'
        endif
        allocate( vv2(ixlen, iylen, ilen),stat=ier)
        if(ier /= 0) then
          write(*,*) ' Could not allocate vv2'
          stop 'Allocate vv2'
        endif
      endif
  
!     get the index values in all directions
      status = nf90_get_var(sid, xvid, xind)
      call check(status)
      status = nf90_get_var(sid, yvid, yind)
      call check(status)
      status = nf90_get_var(sid, lvid, lind)
      call check(status)
!     Read uu2
      status = nf90_get_var(sid,svid,uu2,start=(/1,1,1/),count=(/ixlen, iylen, ilen/),stride=(/1,1,1/))
      call check(status)
!     Read vv2
      status = nf90_get_var(tid,tvid,vv2,start=(/1,1,1/),count=(/ixlen, iylen, ilen/),stride=(/1,1,1/))
      call check(status)

      status = nf90_close(sid)
      call check(status)
      status = nf90_close(tid)
      call check(status)

!     Extracting a subset of all lat,lon,h,vertical pos data
      if(ifile==0) then
        allocate( lat(ixlen,iylen))
        allocate( lon(ixlen,iylen))
        allocate( zm(ixlen, iylen, ilen))
        allocate( h(ixlen,iylen))
        allocate( kbp(ixlen,iylen))
        allocate( ihope(ixlen,iylen))
      endif

      status = nf90_open(trim(hnc), nf90_nowrite, hid)
      call check(status)
      status = nf90_inq_varid(hid, "bathymetry", hvid)
      call check(status)
      status = nf90_get_var(hid,hvid,h,start=(/xind(1),yind(1)/),count=(/ixlen, iylen/),stride=(/1,1/))
      call check(status)
      status = nf90_close(hid)
!     TODO  call read2dvar(hnc, 'bathymetry', xind(1), ixlen, yind(1), iylen, h)

      status = nf90_open(trim(latnc), nf90_nowrite, latid)
      call check(status)
      status = nf90_inq_varid(latid, "Lat", latvid)
      call check(status)
      status = nf90_get_var(latid,latvid,lat,start=(/xind(1),yind(1)/),count=(/ixlen, iylen/),stride=(/1,1/))
      call check(status)
      status = nf90_close(latid)
!     TODO  call read2dvar(latnc, 'Lat', xind(1), ixlen, yind(1), iylen, lat)

      status = nf90_open(trim(lonnc), nf90_nowrite, lonid)
      call check(status)
      status = nf90_inq_varid(lonid, "Long", lonvid)
      call check(status)
      status = nf90_get_var(lonid,lonvid,lon,start=(/xind(1),yind(1)/),count=(/ixlen, iylen/),stride=(/1,1/))
      call check(status)
      status = nf90_close(lonid)
!     TODO  call read2dvar(lonnc, 'Long', xind(1), ixlen, yind(1), iylen, lon)

!     Read level depths information
      status = nf90_open(trim(zmnc), nf90_nowrite, zmid)
      call check(status)
      status = nf90_inq_varid(zmid, "zm", zmvid)
      call check(status)
      status = nf90_get_var(zmid,zmvid,zm,start=(/xind(1),yind(1),lind(1)/),count=(/ixlen, iylen, ilen/),stride=(/1,1,1/))
      call check(status)
      status = nf90_close(zmid)
!     TODO  call read3dvar(zmnc, 'zm', xind(1), ixlen, yind(1), iylen, lind(1), ilen, zm)

!
!     At this point all variables have been read, you may proceed with processing.
!     top layer is 1, ilen is bottom layer.
!
      lon=lon-360 !convert to our long.
!      print*, ixlen,iylen,ilen
!     Compute bottom indices and extend S,T data below bottom
      do i=1,ixlen
        do j=1,iylen
          kbp(i,j)=-99
          do k=ilen,1,-1
            if(zm(i,j,k)>-1.e20) then
              kbp(i,j)=k
              exit
            endif
          enddo !k

          do k=1,ilen
            if(uu2(i,j,k)<-99) uu2(i,j,k)=0
            if(vv2(i,j,k)<-99) vv2(i,j,k)=0
          enddo !k
        enddo !j
      enddo !i  

!     Write a test file of build points
!      write(15, *)'Test NCOM output.'
!      write(15, *)200*iylen
!      do i=1,200 !ixlen 
!        do j=1,iylen 
!!         top layer is 1, ilen is bottom layer.
!          write(15, *) (i-1)*iylen+j,lon(i,j),lat(i,j),vv2(i,j,2) !-h(i,j) 
!        enddo
!      enddo
!      stop

!     Sub-sample data to cut down domain
      ix1=92; ix2=200
      iy1=206; iy2=iylen
 
!     Output hvel.64 etc
      if(ifile==0) then !header

!       Compute hgrid
        np=(ix2-ix1+1)*(iy2-iy1+1)
        ne=2*(ix2-ix1)*(iy2-iy1)
        if(np>mnp.or.ne>mne) stop 'Increase mnp/mne'
        do i=ix1,ix2
          do j=iy1,iy2
            nd=(i-ix1)*(iy2-iy1+1)+j-iy1+1
            xl(nd)=lon(i,j)*pi/180
            yl(nd)=lat(i,j)*pi/180
            call cpp(xout(nd),yout(nd),xl(nd),yl(nd),slam0,sfea0)
            dp(nd)=-h(i,j) !reverse sign
            kbp00(nd)=1 !not kbp(i,j) as NCOM vertical order is reversed (1 is top)
          enddo !j
        enddo !i
      
        do i=ix1,ix2-1
          do j=iy1,iy2-1
            ibox=(i-ix1)*(iy2-iy1)+j-iy1+1
            n1=(i-ix1)*(iy2-iy1+1)+j-iy1+1
            n2=n1+iy2-iy1+1
            n3=n2+1
            n4=n1+1
            if(n3>np) stop 'Overflow 1'
            if(2*ibox>ne) stop 'Overflow 2'
            elnode(1,2*ibox-1)=n1
            elnode(2,2*ibox-1)=n2
            elnode(3,2*ibox-1)=n3
            elnode(1,2*ibox)=n1
            elnode(2,2*ibox)=n3
            elnode(3,2*ibox)=n4
          enddo !j
        enddo !i
!       Output hgrid.gr3 and hgrid.ll
        open(13,file='hgrid.ll')
        open(14,file='hgrid.gr3')
        write(13,*)
        write(13,*)ne,np
        write(14,*)
        write(14,*)ne,np
        do i=1,np
          write(13,*)i,xl(i),yl(i),dp(i)
          write(14,*)i,xout(i),yout(i),dp(i)
        enddo !i
        do i=1,ne
          write(13,*)i,3,(elnode(l,i),l=1,3)
          write(14,*)i,3,(elnode(l,i),l=1,3)
        enddo !i
        close(13)
        close(14)

!       Compute z coordinates
        do i=1,np
          do k=kz,nvrt
            kin=k-kz+1
            hmod2=max(0.1,min(dp(i),h_s))
            if(hmod2<=h_c) then
              z(i,k)=sigma(kin)*hmod2
            else
              z(i,k)=h_c*sigma(kin)+(hmod2-h_c)*cs(kin)
            endif
          enddo !k

!         Z-levels; shallow pts have junk values
          do k=1,kz-1
            z(i,k)=ztot(k)
          enddo !k
        enddo !i

        data_format='DataFormat v5.0'
        version='v1.5h'
        start_time='01/28/2007  00:00:00 PST'
        variable_nm='velocity'
        variable_dim='3D vector'
        irec=0
        do m=1,48/nbyte
          write(63,rec=irec+m) data_format(nbyte*(m-1)+1:nbyte*m)
        enddo
        irec=irec+48/nbyte
        do m=1,48/nbyte
          write(63,rec=irec+m) version(nbyte*(m-1)+1:nbyte*m)
        enddo
        irec=irec+48/nbyte
        do m=1,48/nbyte
          write(63,rec=irec+m) start_time(nbyte*(m-1)+1:nbyte*m)
        enddo
        irec=irec+48/nbyte
        do m=1,48/nbyte
          write(63,rec=irec+m) variable_nm(nbyte*(m-1)+1:nbyte*m)
        enddo
        irec=irec+48/nbyte
        do m=1,48/nbyte
          write(63,rec=irec+m) variable_dim(nbyte*(m-1)+1:nbyte*m)
        enddo
        irec=irec+48/nbyte

        if(itype==1) then
          ivs=1
          i23d=2
        else if(itype==2) then
          ivs=2
          i23d=3
        else
          ivs=1
          i23d=3
        endif
        write(63,rec=irec+1) ndays+1
        write(63,rec=irec+2) 86400. !dtout
        write(63,rec=irec+3) 1 !nspool
        write(63,rec=irec+4) ivs
        write(63,rec=irec+5) i23d
        irec=irec+5


!       Vertical grid
        h_0=0.1
        write(63,rec=irec+1) nvrt
        write(63,rec=irec+2) kz
        write(63,rec=irec+3) h_0
        write(63,rec=irec+4) h_s
        write(63,rec=irec+5) h_c
        write(63,rec=irec+6) theta_b
        write(63,rec=irec+7) theta_f
        irec=irec+7
        do k=1,kz-1
          write(63,rec=irec+k) ztot(k)
        enddo
        do k=kz,nvrt
          kin=k-kz+1
          write(63,rec=irec+k) sigma(kin)
        enddo
        irec=irec+nvrt

!       Horizontal grid
        write(63,rec=irec+1) np
        write(63,rec=irec+2) ne
        irec=irec+2
        do m=1,np
          write(63,rec=irec+1)xout(m)
          write(63,rec=irec+2)yout(m)
          write(63,rec=irec+3)dp(m)
          write(63,rec=irec+4)kbp00(m)

!          write(88,*)m,kbp00(m)

          irec=irec+4
        enddo !m=1,np
        do m=1,ne
          write(63,rec=irec+1)3
          irec=irec+1
          do mm=1,3
            write(63,rec=irec+1)elnode(mm,m)
            irec=irec+1
          enddo !mm
        enddo !m
      endif !ifile=0; header

!     Non-header part
      do i=ix1,ix2
        do j=iy1,iy2
          nd=(i-ix1)*(iy2-iy1+1)+j-iy1+1
          do k=1,nvrt
            vel(nd,k,1)=uu2(i,j,1) !surface vel.
            vel(nd,k,2)=vv2(i,j,1)
          enddo !k
        enddo !j
      enddo !i

      write(63,rec=irec+1) ifile*86400.
      write(63,rec=irec+2) ifile
      irec=irec+2
      do i=1,np
        write(63,rec=irec+i) 0.
      enddo !i
      irec=irec+np

      do i=1,np
        if(itype==1) then
          write(63,rec=irec+1) 0.
          irec=irec+1
        else if(itype==2) then
          do k=max0(1,kbp00(i)),nvrt
            write(63,rec=irec+1) vel(i,k,1)
            write(63,rec=irec+2) vel(i,k,2)
            irec=irec+2
          enddo !k
        else !vert.63
          do k=max0(1,kbp00(i)),nvrt
            write(63,rec=irec+1) 0.
            irec=irec+1
          enddo !k
        endif
      enddo !i
      
!-----------------------------------------------------------------------------------------------------------
      enddo !ifile=0,ndays

      print*, 'Finished'

!     End of main
!     Subroutines
      contains
!     Internal subroutine - checks error status after each netcdf, 
!     prints out text message each time an error code is returned. 
      subroutine check(status)
      integer, intent ( in) :: status
    
      if(status /= nf90_noerr) then 
        print *, trim(nf90_strerror(status))
      end if
      end subroutine check  

!     A simpler interface for reading 1, 2, 3d variables
!     Not activated right now.
!
!     fname = Netcdf filename
!     vname = variable name
!     nix = starting index to read >= 1
!     n = number to read
!     var = return value
      subroutine read1dvar(fname, vname, nix, n, var)
      character(len=*), intent (in) :: fname
      character(len=*), intent (in) :: vname
      integer, intent ( in) :: nix, n
      real (kind = FourByteReal), dimension(n), intent(out) :: var
      integer ncid, nvid, status
      status = nf90_open(trim(fname), nf90_nowrite, ncid)
      call check(status)
      status = nf90_inq_varid(ncid, vname, nvid)
      call check(status)
      status = nf90_get_var(ncid,nvid,var,start=(/nix/),count=(/n/),stride=(/1/))
      call check(status)
      end subroutine read1dvar

!     fname = Netcdf filename
!     vname = variable name
!     nix = starting index to read >= 1
!     n = number to read
!     mix = starting index to read >= 1
!     m = number to read
!     var = return value
      subroutine read2dvar(fname, vname, nix, n, mix, m, var)
      character(len=*), intent (in) :: fname
      character(len=*), intent (in) :: vname
      integer, intent ( in) :: nix, n, mix, m
      real (kind = FourByteReal), dimension(n, m), intent(out) :: var
      integer ncid, nvid, status
      status = nf90_open(trim(fname), nf90_nowrite, ncid)
      call check(status)
      status = nf90_inq_varid(ncid, vname, nvid)
      call check(status)
      status = nf90_get_var(ncid,nvid,var,start=(/nix,mix/),count=(/n,m/),stride=(/1,1/))
      call check(status)
      end subroutine read2dvar

!     fname = Netcdf filename
!     vname = variable name
!     nix = starting index to read >= 1
!     n = number to read
!     mix = starting index to read >= 1
!     m = number to read
!     lix = starting index to read >= 1
!     l = number to read
!     var = return value
      subroutine read3dvar(fname, vname, nix, n, mix, m, lix, l, var)
      character(len=*), intent (in) :: fname
      character(len=*), intent (in) :: vname
      integer, intent ( in) :: n, m, l, nix, mix, lix
      real (kind = FourByteReal), dimension(n, m, l), intent(out) :: var
      integer ncid, nvid, status
      status = nf90_open(trim(fname), nf90_nowrite, ncid)
      call check(status)
      status = nf90_inq_varid(ncid, vname, nvid)
      call check(status)
      status = nf90_get_var(ncid,nvid,var,start=(/nix,mix,lix/),count=(/n,m,l/),stride=(/1,1,1/))
      call check(status)
      status = nf90_close(ncid)
      end subroutine read3dvar

      end program readNCOM

!      function signa(x1,x2,x3,y1,y2,y3)
!
!      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2
!
!      return
!      end

      subroutine cpp(x,y,rlambda,phi,rlambda0,phi0)
!      implicit real(kind=dbl_kind1)(a-h,o-z), integer(i-n)

      r=6378206.4
      x=r*(rlambda-rlambda0)*cos(phi0)
      y=phi*r

      return
      end
