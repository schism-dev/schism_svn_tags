! Read NRL/NCOM_CCS netcdf files and do analyses (time series).
!
! X = longitude dimension is sampled in the salt and temp data files.
! Y = latitude dimension is sampled in the salt and temp data files.
! The level dimensions goes from 1 = surface, 40 = bottom.
!
! Compilation on amb6402:
!
! ifort -O3 -Bstatic -o analyze_ncomccs1 analyze_ncomccs1.f90 -Vaxlib  -I/usr/local/netcdf/include/ -L/usr/local/netcdf/lib -lnetcdf
!
!   Input: 
!     (1) station.bp: depth is z (>=0 from F.S.); x,y are lon., lat in degrees;
!     (2) date.in: 1st line: year, month, day, starting CORIE day, # of days; 2nd line: ifiletype (1: T,S; 2: UV).
!   Output: (1) fort.18 (T or u), fort.19 (S or v).
!           (2) fort.11: fatal and non-fatal errors;

      program readNCOM
!     netcdf modules from symlinks
      use typeSizes
      use netcdf

!      implicit real*8(a-h,o-z)

!      integer, parameter :: debug=1
      integer, parameter :: mnp=50000
      integer, parameter :: mne=100000
      integer, parameter :: mns=150000
      integer, parameter :: mnv=60
      integer, parameter :: mnei=20 !neighbor
      integer, parameter :: nbyte=4
      real, parameter :: small1=1.e-3 !used to check area ratios
  
!     netcdf related variables
      integer :: hid, latid, lonid, zmid, sid, tid ! Netcdf file IDs
      integer :: latvid, lonvid, zmvid, hvid ! positional variables
      integer :: xdid, xvid ! longitude index
      integer :: ydid, yvid ! latitude index
      integer :: ldid, lvid ! vertical level, 1 is top
      integer :: svid, tvid ! salt & temp variable IDs
      integer, dimension(nf90_max_var_dims) :: dids
             
!     Lat, lon, bathymetry
      real (kind = FourByteReal), dimension(:), allocatable :: xind, yind
      real (kind = FourByteReal), dimension(:,:), allocatable :: lat, lon , depths0, depths
!     Vertical postion, salinity, and temperature
      real (kind = FourByteReal), dimension(:,:,:), allocatable :: zm0,zm,salt, temp
      integer, dimension(:,:), allocatable :: kbp
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
!     integer :: i,j,k,i1,i2,i3,j1,j2,j3

      dimension xl(mnp),yl(mnp),nm(mne,4),dp(mnp),i34(mne)
      dimension ixy(mnp,2),arco(4)
      dimension wild(100),wild2(100,2)
      dimension tempout(mnp), saltout(mnp),month_day(12)

!     First statement
      open(10,file='date.in',status='old')
      read(10,*) iyear,month0,iday0,corieday,ndays
      read(10,*) ifiletype
      if(ifiletype/=1.and.ifiletype/=2) stop ' Wrong ifiletype'
      close(10)

!     Read in station.bp
      open(14,file='station.bp',status='old')
      read(14,*)
      read(14,*)np
      if(np.gt.mnp) then
        write(*,*)'Increase mnp'
        stop
      endif
      do i=1,np
        read(14,*)j,xl(i),yl(i),dp(i)
        if(dp(i)<0) stop 'Depth must be non-negative' 
      enddo !i
      close(14)

!     Read NCOM results
! geometry files - should never change.
!
      hnc = '/home/workspace/ccalmr/nrldata_links/ncom_ccs/model_h.h002.nc'
      zmnc = '/home/workspace/ccalmr/nrldata_links/ncom_ccs/model_zm.h002.nc'
!'
  
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

      do it=0,1 !2 records per day
!==============================================================================
      timeout=corieday+ifile+it/2.0-1./3 !UTC to PST

!     Create salt and temp (or UV) file names for this date.
      if(it==0) then
        char3='00'
      else
        char3='12'
      endif
      if(ifiletype==1) then !S,T
        sfile='/home/workspace/ccalmr8/frolovs/ncomccs/output/s3d/s3d.ncom.h002.401.nest1.'//trim(yr)//trim(md)//char3//'_00000000.nc'
        tfile='/home/workspace/ccalmr8/frolovs/ncomccs/output/t3d/t3d.ncom.h002.401.nest1.'//trim(yr)//trim(md)//char3//'_00000000.nc'
      else !UV
        sfile='/home/workspace/ccalmr8/frolovs/ncomccs/output/u3d/u3d.ncom.h002.401.nest1.'//trim(yr)//trim(md)//char3//'_00000000.nc'
        tfile='/home/workspace/ccalmr8/frolovs/ncomccs/output/v3d/v3d.ncom.h002.401.nest1.'//trim(yr)//trim(md)//char3//'_00000000.nc'
      endif

!     Open files for salinity and temperature.
      status = nf90_open(trim(sfile), nf90_nowrite, sid)
      call check(status)
      if(status /= nf90_noerr) cycle
      status = nf90_open(trim(tfile), nf90_nowrite, tid)
      call check(status)
      if(status /= nf90_noerr) cycle

!     Get index information for sampled grid, assumed the same for temperature.
      status = nf90_inq_dimid(sid, "Longitude", xdid)
      call check(status)
      status = nf90_inq_dimid(sid, "Latitude", ydid)
      call check(status)
      status = nf90_inq_dimid(sid, "Depth", ldid)
      call check(status)

      status = nf90_inq_varid(sid, "Longitude", xvid)
      call check(status)
      status = nf90_inq_varid(sid, "Latitude", yvid)
      call check(status)
      status = nf90_inq_varid(sid, "Depth", lvid)
      call check(status)

!     get the variable ids from salt and temperature files
      if(ifiletype==1) then
        status = nf90_inq_varid(sid, "Salinity", svid)
        call check(status)
        status = nf90_inq_varid(tid, "Temperature", tvid)
        call check(status)
      else
        status = nf90_inq_varid(sid, "U3D", svid)
        call check(status)
        status = nf90_inq_varid(tid, "V3D", tvid)
        call check(status)
      endif

!     Get the lengths of the dimensions from the salt file.
!     Assumed same as in temperature
      status = nf90_Inquire_Variable(sid, svid, dimids = dids)
      call check(status)
      status = nf90_Inquire_Dimension(sid, dids(1), len = ixlen)
      call check(status)
      status = nf90_Inquire_Dimension(sid, dids(2), len = iylen)
      call check(status)
      status = nf90_Inquire_Dimension(sid, dids(3), len = ilen)
      call check(status)

!     allocate memory.
      if(ifile==0.and.it==0) then
        allocate( xind(ixlen),stat=ier)
        allocate( yind(iylen),stat=ier)
        allocate( salt(ixlen, iylen, ilen),stat=ier) !also used for v
        allocate( lat(ixlen,iylen))
        allocate( lon(ixlen,iylen))
        allocate( depths0(ixlen,iylen))
        allocate( depths(ixlen,iylen))
        allocate( zm0(ixlen, iylen, ilen))
        allocate( zm(ixlen, iylen, ilen))
        allocate( kbp(ixlen,iylen))
        if(ier /= 0) then
          write(*,*) ' Could not allocate salt'
          stop 'Allocate salt'
        endif
        allocate( temp(ixlen, iylen, ilen),stat=ier) !also used for u
        if(ier /= 0) then
          write(*,*) ' Could not allocate temp'
          stop 'Allocate temp'
        endif
      endif

!     get lon, lat 
      status = nf90_get_var(sid, xvid, xind) !,start=(/1/),count=(/ixlen/),stride=(/1/))
      call check(status)
      status = nf90_get_var(sid, yvid, yind)
      call check(status)
!      status = nf90_get_var(sid, lvid, sigmaz)
!      call check(status)

!     Bathymetry (not changing with time)
!     Note: order of y reversed!
      status = nf90_open(trim(hnc), nf90_nowrite, hid)
      call check(status)
      status = nf90_inq_varid(hid, "bathymetry", hvid)
      call check(status)
      status = nf90_get_var(hid,hvid,depths0)
      call check(status)
      status = nf90_close(hid)
      depths(:,1:iylen)=-depths0(:,iylen:1:-1)

!     Read level depths information (not changing with time)
!     Note: order of y reversed!
      status = nf90_open(trim(zmnc), nf90_nowrite, zmid)
      call check(status)
      status = nf90_inq_varid(zmid, "zm", zmvid)
      call check(status)
      status = nf90_get_var(zmid,zmvid,zm0)
      call check(status)
      status = nf90_close(zmid)
      zm(:,1:iylen,:)=zm0(:,iylen:1:-1,:)

!      print*, 'Step5:',xind(ixlen),yind(iylen)
!      stop

!     Read salinity/U
      status = nf90_get_var(sid,svid,salt) !,start=(/1,1,1/),count=(/ixlen, iylen, ilen/),stride=(/1,1,1/))
      call check(status)
!     Read temperature/V
      status = nf90_get_var(tid,tvid,temp) !,start=(/1,1,1/),count=(/ixlen, iylen, ilen/),stride=(/1,1,1/))
      call check(status)

      status = nf90_close(sid)
      call check(status)
      status = nf90_close(tid)
      call check(status)

!     At this point all variables have been read, you may proceed with processing.
!     Variables: [salt,temp,zm](1:ixlen,1:iylen,1:ilen) (level 1 is surface).
!     Compute lon, lat
!      h_s=150 !demarcation depth between sigmna and z
      do i=1,ixlen
        do j=1,iylen
          lon(i,j)=xind(i)
          lat(i,j)=yind(j)
        enddo !j
      enddo !i

!      write(16,*)zm(1,1,1:ilen)

!     Compute bottom indices and extend S,T data below bottom
      do i=1,ixlen
        do j=1,iylen
          kbp(i,j)=-99
          do k=ilen,1,-1
            if(salt(i,j,k)>-1.e20) then
              kbp(i,j)=k
              exit
            endif
          enddo !k

          if(kbp(i,j)/=-99) then !valid pts (wet)
            do k=1,ilen
              if(k<=kbp(i,j)) then
                if(temp(i,j,k)<-99) then
                  write(11,*)'Fatal: no valid S,T:',i,j,k,salt(i,j,k),temp(i,j,k)
                  stop
                endif
              else !extend
                salt(i,j,k)=salt(i,j,kbp(i,j))
                temp(i,j,k)=temp(i,j,kbp(i,j))
              endif
            enddo !k
          endif
        enddo !j
      enddo !i  

!     Write a test file of build points
!     top layer is 1, ilen is bottom layer.
!      write(15, *)'Test NCOM output.'
!      write(15, *)ixlen * iylen
!      do i=1,ixlen
!        do j=1,iylen
!!         write(15, *) (i-1)*iylen+j,lon(i,j),lat(i,j),depths(i,j) 
!          if(kbp(i,j)==-99) then
!            ibot=ilen
!          else
!            ibot=kbp(i,j)
!          endif
!!          write(15, *) (i-1)*iylen+j,lon(i,j),lat(i,j),salt(i,j,1)
!          if(salt(i,j,1)>-99) write(15, *) lon(i,j),lat(i,j),salt(i,j,1),temp(i,j,1)
!        enddo !j
!      enddo !i
!      stop

!     Find parent elements and levels for station.bp, and do interpolation
      tempout=-99; saltout=-99
      ixy=0
      loop4: do ix=1,ixlen-1 
        do iy=1,iylen-1 
          x1=lon(ix,iy); x2=lon(ix+1,iy); x3=lon(ix+1,iy+1); x4=lon(ix,iy+1)
          y1=lat(ix,iy); y2=lat(ix+1,iy); y3=lat(ix+1,iy+1); y4=lat(ix,iy+1)
          b1=abs(signa(x1,x2,x3,y1,y2,y3))
          b2=abs(signa(x1,x3,x4,y1,y3,y4))
          if(b1+b2==0) cycle
          ialldone=1 !to see if all pts have been interpolated
          do i=1,np
            if(ixy(i,1)/=0) cycle
            ialldone=0
            a1=abs(signa(xl(i),x1,x2,yl(i),y1,y2))
            a2=abs(signa(xl(i),x2,x3,yl(i),y2,y3))
            a3=abs(signa(xl(i),x3,x4,yl(i),y3,y4))
            a4=abs(signa(xl(i),x4,x1,yl(i),y4,y1))
            rat=abs(a1+a2+a3+a4-b1-b2)/(b1+b2)

            !write(9,*)ix,iy,i,xl(i),yl(i),x1,y1,x2,y2,x3,y3,x4,y4
!            write(9,*)ix,iy,i,rat,a1,a2,a3,a4,b1+b2

            if(rat<small1) then
              ixy(i,1)=ix; ixy(i,2)=iy
!             Find a triangle
              in=0 !flag
              do l=1,2
                ap=abs(signa(xl(i),x1,x3,yl(i),y1,y3))
                if(l==1) then !nodes 1,2,3
                  bb=abs(signa(x1,x2,x3,y1,y2,y3))
                  wild(l)=abs(a1+a2+ap-bb)/bb
                  if(wild(l)<small1*5) then
                    in=1
                    arco(1)=min(1.,a2/bb)
                    arco(2)=min(1.,ap/bb)
                    arco(3)=max(0.,min(1.,1-arco(1)-arco(2)))
                    arco(4)=0
                    exit
                  endif
                else !nodes 1,3,4
                  bb=abs(signa(x1,x3,x4,y1,y3,y4))
                  wild(l)=abs(a3+a4+ap-bb)/bb
                  if(wild(l)<small1*5) then
                    in=2
                    arco(1)=min(1.,a3/bb)
                    arco(3)=min(1.,a4/bb)
                    arco(4)=max(0.,min(1.,1-arco(1)-arco(3)))
                    arco(2)=0
                    exit
                  endif
                endif
              enddo !l=1,2
              if(in==0) then
                write(11,*)'Cannot find a triangle:',(wild(l),l=1,2)
                stop
              endif
              if(sum(arco)/=1) then
                write(11,*)'Sum of coordinates/=1:',arco(1:4)
                stop
              endif

!             Find vertical level
              if(kbp(ix,iy)==-99.or.kbp(ix+1,iy)==-99.or.kbp(ix,iy+1)==-99.or.kbp(ix+1,iy+1)==-99) then
                if(ifile==0) then
                  write(11,*)'One of parent point is dry:',i,ix,iy,kbp(ix,iy),kbp(ix+1,iy),kbp(ix,iy+1),kbp(ix+1,iy+1)
                  write(11,*)(ix-1)*iylen+iy,ix*iylen+iy,(ix-1)*iylen+iy+1,ix*iylen+iy+1
                endif
                lev=1; vrat=0
!               Find first wet pt
                ifound=0
                arco=0
                do l=1,4
                  if(l==1) then
                    in1=ix; in2=iy
                  else if(l==2) then
                    in1=ix+1; in2=iy
                  else if(l==3) then
                    in1=ix+1; in2=iy+1
                  else
                    in1=ix; in2=iy+1
                  endif   
                  if(kbp(in1,in2)/=-99) then
                    ifound=l
                    arco(l)=1
                    exit
                  endif
                enddo !l
                if(ifound==0) then
                  write(11,*)'All parents dry:',i,ix,iy
                  stop
                endif
              else !all 4 pts wet
                in1=ix; in2=iy
              endif

              if(-dp(i)<zm(in1,in2,kbp(in1,in2))) then
                lev=kbp(in1,in2)-1; vrat=1
              else if(-dp(i)>zm(in1,in2,1)) then !above f.s.
                lev=1; vrat=0
              else
                lev=-99 !flag
                do kk=1,kbp(in1,in2)-1
                  if(-dp(i)<=zm(in1,in2,kk).and.-dp(i)>=zm(in1,in2,kk+1)) then
                    lev=kk
                    vrat=(zm(in1,in2,kk)+dp(i))/(zm(in1,in2,kk)-zm(in1,in2,kk+1))
                    if(vrat<0.or.vrat>1) then
                      write(11,*)'vrat out of bound:',vrat
                      stop
                    endif
                    exit
                  endif
                enddo !kk
                if(lev==-99) then
                  write(11,*)'Cannot find a level:',i,k,-dp(i),(zm(in1,in2,l),l=1,kbp(in1,in2))
                  stop
                endif
              endif
           
              wild2(1,1)=temp(ix,iy,lev)*(1-vrat)+temp(ix,iy,lev+1)*vrat
              wild2(1,2)=salt(ix,iy,lev)*(1-vrat)+salt(ix,iy,lev+1)*vrat
              wild2(2,1)=temp(ix+1,iy,lev)*(1-vrat)+temp(ix+1,iy,lev+1)*vrat
              wild2(2,2)=salt(ix+1,iy,lev)*(1-vrat)+salt(ix+1,iy,lev+1)*vrat
              wild2(3,1)=temp(ix+1,iy+1,lev)*(1-vrat)+temp(ix+1,iy+1,lev+1)*vrat
              wild2(3,2)=salt(ix+1,iy+1,lev)*(1-vrat)+salt(ix+1,iy+1,lev+1)*vrat
              wild2(4,1)=temp(ix,iy+1,lev)*(1-vrat)+temp(ix,iy+1,lev+1)*vrat
              wild2(4,2)=salt(ix,iy+1,lev)*(1-vrat)+salt(ix,iy+1,lev+1)*vrat
              tempout(i)=wild2(1,1)*arco(1)+wild2(2,1)*arco(2)+wild2(3,1)*arco(3)+wild2(4,1)*arco(4)
              saltout(i)=wild2(1,2)*arco(1)+wild2(2,2)*arco(2)+wild2(3,2)*arco(3)+wild2(4,2)*arco(4)

!             Debug
              if(ifile==0.and.it==0) then
                write(16,*)i,kbp(ix,iy),kbp(ix+1,iy),kbp(ix+1,iy+1),kbp(ix,iy+1)
                write(16,*)i,depths(ix,iy)*arco(1)+depths(ix+1,iy)*arco(3)+ &
     &depths(ix+1,iy+1)*arco(3)+depths(ix,iy+1)*arco(4)
              endif

              if(abs(tempout(i))>98.or.abs(saltout(i))>98) then
                write(11,*)'Failed to interpolate:',i,tempout(i),saltout(i)
                stop
              endif 
            endif !rat<small1
          enddo !i=1,np
          if(ialldone==1) exit loop4 
        enddo !iy=1,iylen-1
      enddo loop4 !ix=1,ixlen-1

      do i=1,np
        if(ixy(i,1)==0.or.ixy(i,2)==0) write(*,*)'Cannot find a parent element for pt:',i,ixy(i,1:2)
!'
      enddo !i
      write(18,'(e13.6,6000(1x,f12.3))')timeout,(saltout(i),i=1,np) !S or u; time in PST
      write(19,'(e13.6,6000(1x,f12.3))')timeout,(tempout(i),i=1,np) !T or v
!==============================================================================
      enddo !it=0,1

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

      function signa(x1,x2,x3,y1,y2,y3)

      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2

      return
      end

