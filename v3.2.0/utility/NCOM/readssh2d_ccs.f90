! Read NRL/NCOM netcdf files
!
! pturner 12-2004
!
! X = longitude dimension is sampled in the salt and temp data files.
! Y = latitude dimension is sampled in the salt and temp data files.
! The level dimensions goes from 1 = surface, 40 = bottom.
!
! Compilation on canopus:
!
! ifort -g -O3 -Bstatic -o readssh2d_ccs_canopus readssh2d_ccs.f90 -Vaxlib -I/usr/local/include -L/usr/local/lib -lnetcdf
!
! Need to have the following modules symlinked into the current directory
!
! NETCDF.mod -> /usr/local/netcdf/include/NETCDF.mod
! TYPESIZES.mod -> /usr/local/netcdf/include/TYPESIZES.mod
!

!zyl
!   Input: 
!     (1) hgrid.gr3 (only for depth and b.c.);
!     (2) hgrid.ll;
!     (3) date.in: year, month, day, # of days (must=0); 2nd line: ireduce; 
!                  3rd line: # of open bnds; list of open bnds
!   Output: 
!     (1) fort.99: original SSH on the lat/lon from NRL; build pts, at the beginning and end 
!                  of the period specified. This can be used to "extrapolate" ssh
!                  at invalid pts
!     (2) Z0.out: Z0 at the specified open bnds. "Invalid" if parent nodes are all dry 

      program readNCOM
!     netcdf modules from symlinks
      use typeSizes
      use netcdf

!      implicit none

!      integer, parameter :: debug=1
      integer, parameter :: mnp=60000
      integer, parameter :: mne=100000
      integer, parameter :: mns=130000
      integer, parameter :: mnv=90
      integer, parameter :: mnope=9 !# of open bnd segements
      integer, parameter :: mnond=1000 !max. # of open-bnd nodes on each segment
      integer, parameter :: mnei=20 !neighbor
      integer, parameter :: nbyte=4
      real, parameter :: small1=1.e-3 !used to check area ratios
  
!     netcdf related variables
      integer :: hid, latid, lonid, sid ! Netcdf file IDs
      integer :: latvid, lonvid,  hvid ! positional variables
      integer :: xdid, xvid ! longitude index
      integer :: ydid, yvid ! latitude index
!      integer :: ldid, lvid ! vertical level, 1 is top
      integer :: svid ! elevation variable IDs
      integer, dimension(nf90_max_var_dims) :: dids
             
!     Local variables for data
      real (kind = FourByteReal), dimension(:), allocatable :: xind, yind 
!     Lat, lon, bathymetry
      real (kind = FourByteReal), dimension(:,:), allocatable :: lat, lon, h, ssh
!     Vertical postion, salinity, and temperature
!      real (kind = FourByteReal), dimension(:,:,:), allocatable :: zm, salt, temp
!      integer, dimension(:,:), allocatable :: kbp,ihope
!     File names for netcdf files
      character(len=1024) :: hnc, latnc, lonnc, zmnc, sfile !, tfile
!     Command line arguments
      character(len=1024) :: s, yr, md
      character(len=4) :: iyear_char
      character(len=1) :: char1,char2
      character(len=2) :: char3,char4
!     external function for number of command line arguments
      integer :: iargc

      integer :: status ! netcdf local status variable
      integer :: ier ! allocate error return.
      integer :: ixlen, iylen !, ilen ! sampled lengths in each coordinate direction
      integer :: ixlen1, iylen1,ixlen2, iylen2 !reduced indices for CORIE grid to speed up interpolation
!     integer :: i,j,k,i1,i2,i3,j1,j2,j3

      dimension xl(mnp),yl(mnp),nm(mne,3),dp(mnp),xcor(mnp),ycor(mnp)
!      dimension ztot(mnv),sigma(mnv),cs(mnv),z(mnp,mnv),iest(mnp)
      dimension wild(100),wild2(100,2),ixy(mnp,2),arco(3)
      dimension etaout(mnp),month_day(12)
!      dimension tsd(mns,mnv),ssd(mns,mnv)
      dimension nne(mnp),ine(mnp,mnei),ic3(mne,3),nx(3,2),js(mne,3),is(mns,2),isidenode(mns,2)
      dimension xcj(mns),ycj(mns)
      dimension nond(mnope),iond(mnope,mnond),isbnd(mnp),iob(mnope),z0_out(2,mnope,mnond)
      dimension ival(mnp),ival_out(2,mnope,mnond)

!     First statement
!     Read in hgrid and vgrid
!      open(17,file='estuary.gr3',status='old')
      open(16,file='hgrid.ll',status='old')
      open(14,file='hgrid.gr3',status='old') !only need depth and bnd info 
!      open(20,file='ssh.ll')
      read(14,*)
      read(14,*)ne,np
      if(np.gt.mnp.or.ne.gt.mne) then
        write(*,*)'Increase mnp/mne'
        stop
      endif
      read(16,*)
      read(16,*)
      do i=1,np
        read(14,*)j,xcor(i),ycor(i),dp(i) !xcor, ycor for debugging only 
        read(16,*)j,xl(i),yl(i) !,dp(i)
      enddo !i
      do i=1,ne
        read(14,*)j,k,(nm(i,l),l=1,3)
      enddo !i

!     Open bnds
      read(14,*) nope
      if(nope>mnope) then
        write(11,*) 'nope > mnope'
        stop
      endif

      read(14,*) neta
      ntot=0
      do k=1,nope
        read(14,*) nond(k)
        if(nond(k)>mnond) then
          write(11,*) 'nond(k) > mnond'
          stop
        endif
        do i=1,nond(k)
          read(14,*) iond(k,i)
          isbnd(iond(k,i))=k
        enddo
        if(iond(k,1)==iond(k,nond(k))) then
          write(11,*)'Looped open bnd:',k
          stop
        endif
        ntot=ntot+nond(k)
      enddo

      close(14)
      close(16)

!     Compute geometry
      do i=1,3
        do j=1,2
          nx(i,j)=i+j
          if(nx(i,j)>3) nx(i,j)=nx(i,j)-3
          if(nx(i,j)<1.or.nx(i,j)>3) then
            write(*,*)'nx wrong',i,j,nx(i,j)
            stop
          endif
        enddo !j
      enddo !i

      do i=1,np
        nne(i)=0
      enddo

      do i=1,ne
        do j=1,3
          nd=nm(i,j)
          nne(nd)=nne(nd)+1
          if(nne(nd)>mnei) then
            write(11,*)'Too many neighbors',nd
            stop
          endif
          ine(nd,nne(nd))=i
        enddo
      enddo

!     Compute ball info; this won't be affected by re-arrangement below
      do i=1,ne
        do j=1,3
          ic3(i,j)=0 !index for bnd sides
          nd1=nm(i,nx(j,1))
          nd2=nm(i,nx(j,2))
          do k=1,nne(nd1)
            ie=ine(nd1,k)
            if(ie/=i.and.(nm(ie,1)==nd2.or.nm(ie,2)==nd2.or.nm(ie,3)==nd2)) ic3(i,j)=ie
          enddo !k
        enddo !j
      enddo !i

      ns=0 !# of sides
      do i=1,ne
        do j=1,3
          nd1=nm(i,nx(j,1))
          nd2=nm(i,nx(j,2))
          if(ic3(i,j)==0.or.i<ic3(i,j)) then !new sides
            ns=ns+1
            if(ns>mns) then
              write(11,*)'Too many sides'
              stop
            endif
            js(i,j)=ns
            is(ns,1)=i
            isidenode(ns,1)=nd1
            isidenode(ns,2)=nd2
            xcj(ns)=(xl(nd1)+xl(nd2))/2
            ycj(ns)=(yl(nd1)+yl(nd2))/2
!            dps(ns)=(dp(nd1)+dp(nd2))/2
!            distj(ns)=dsqrt((x(nd2)-x(nd1))**2+(y(nd2)-y(nd1))**2)
!            if(distj(ns)==0) then
!              write(11,*)'Zero side',ns
!              stop
!            endif
!            thetan=datan2(x(nd1)-x(nd2),y(nd2)-y(nd1))
!            snx(ns)=dcos(thetan)
!            sny(ns)=dsin(thetan)

            is(ns,2)=ic3(i,j) !bnd element => bnd side
!           Corresponding side in element ic3(i,j)
            if(ic3(i,j)/=0) then !old internal side
              iel=ic3(i,j)
              index=0
              do k=1,3
                if(ic3(iel,k)==i) then
                  index=k
                  exit
                endif
              enddo !k
              if(index==0) then
                write(11,*)'Wrong ball info',i,j
                stop
              endif
              js(iel,index)=ns
            endif !ic3(i,j).ne.0
          endif !ic3(i,j)==0.or.i<ic3(i,j)
        enddo !j=1,3
      enddo !i=1,ne

      if(ns<ne.or.ns<np) then
        write(11,*)'Weird grid with ns < ne or ns < np',np,ne,ns
        stop
      endif

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
      hnc = '/home/workspace/ccalmr/nrldata_links/ncom_ccs/model_h.h002.nc'
!      latnc = '/home/workspace/ccalmr/nrldata_links/ncom_ccs/model_lat.h002.nc'
!      lonnc = '/home/workspace/ccalmr/nrldata_links/ncom_ccs/model_lon.h002.nc'
!      zmnc = '/home/workspace/ccalmr/nrldata_links/ncom_ccs/model_zm.h002.nc'
  
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

      open(10,file='date.in',status='old')
      read(10,*) iyear,month0,iday0,ndays
      if(ndays/=0) then
        write(11,*)'ndays /=0'
        stop
      endif
!     ireduce=0: original set; =1: reduced bg set for quicker interpolation (edit ixlen[1,2],iylen[1,2])
      read(10,*) ireduce 
      read(10,*) nob,(iob(l),l=1,nob)
      if(nob>mnope) then
        write(11,*)'nob>mnope'
        stop
      endif
      close(10)

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
          month = 1
          iday = 1
          iyear = iyear + 1
          if(mod(iyear,4)==0) then
              month_day(2)=29
          else
              month_day(2)=28
          endif
          write(iyear_char,'(i4)') iyear
          yr=iyear_char
          write(*,*)'Beyond one year - going to next year', iyear,month,iday
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

!     Create salt and temp file names for this date.
!     Compile spits out a warning here, dunno how to fix.
!     Read in from symlinks in ccalmr (to fill gaps manually) 
!     Use sfile to denote ssh
!      sfile = '/home/workspace/ccalmr/nrldata_links/'//trim(yr)//'/ssh/ssh.glb8_2f_'//trim(yr)//trim(md)//'00.nc'
      sfile='/home/workspace/ccalmr8/frolovs/ncomccs/output/ssh/ssh.ncom.h002.401.nest1.'//trim(yr)//trim(md)//'00_00000000.nc'
!'
      print*, 'opening ',trim(sfile)
!     Open files for ssh
      status = nf90_open(trim(sfile), nf90_nowrite, sid)
      call check(status)

!     Get index information for sampled grid, 
!     These will index into the static lat, lon, bathymetry, and zm files.
      status = nf90_inq_dimid(sid, "Longitude", xdid)
      call check(status)
      status = nf90_inq_dimid(sid, "Latitude", ydid)
      call check(status)
!      status = nf90_inq_dimid(sid, "level", ldid)
!      call check(status)

      status = nf90_inq_varid(sid, "Longitude", xvid)
      call check(status)
      status = nf90_inq_varid(sid, "Latitude", yvid)
      call check(status)
!      status = nf90_inq_varid(sid, "level", lvid)
!      call check(status)

!     get the variable id
      status = nf90_inq_varid(sid, "SSH", svid)
      call check(status)

!     Get the lengths of the dimensions from the ssh file.
      status = nf90_Inquire_Variable(sid, svid, dimids = dids)
      call check(status)
      status = nf90_Inquire_Dimension(sid, dids(1), len = ixlen)
      call check(status)
      status = nf90_Inquire_Dimension(sid, dids(2), len = iylen)
      call check(status)
!      status = nf90_Inquire_Dimension(sid, dids(3), len = ilen)
!      call check(status)

!     allocate memory.
      if(ifile==0) then
        allocate( xind(ixlen),stat=ier)
        allocate( yind(iylen),stat=ier)
        allocate( ssh(ixlen, iylen),stat=ier)
        allocate( lat(ixlen,iylen))
        allocate( lon(ixlen,iylen))
        allocate( h(ixlen,iylen))

        if(ier /= 0) then
          write(*,*) ' Could not allocate ssh'
          stop 'Allocate salt'
        endif
      endif
  
!     get lon, lat
      status = nf90_get_var(sid, xvid, xind)
      call check(status)
      status = nf90_get_var(sid, yvid, yind)
      call check(status)
!      status = nf90_get_var(sid, lvid, lind)
!      call check(status)
!     Read ssh
      status = nf90_get_var(sid,svid,ssh,start=(/1,1/),count=(/ixlen, iylen/),stride=(/1,1/))
      call check(status)

      status = nf90_close(sid)
      call check(status)

!     Extracting depths
!     Note: order of y reversed!
      status = nf90_open(trim(hnc), nf90_nowrite, hid)
      call check(status)
      status = nf90_inq_varid(hid, "bathymetry", hvid)
      call check(status)
      status = nf90_get_var(hid,hvid,h,start=(/1,1/),count=(/ixlen, iylen/),stride=(/1,1/))
      call check(status)
      status = nf90_close(hid)

      h(:,1:iylen)=-h(:,iylen:1:-1) !reverse sign

      do i=1,ixlen
        do j=1,iylen
          lon(i,j)=xind(i)
          lat(i,j)=yind(j)
        enddo !j
      enddo !i

!     Reduce indices ixlen, iylen for CORIE grid to speed up inetrpolation
      if(ireduce==0) then
        ixlen1=1
        ixlen2=ixlen
        iylen1=1
        iylen2=iylen
      else
        ixlen1=49
        ixlen2=80
        iylen1=279
        iylen2=370
        if(ixlen2<ixlen1.or.ixlen2>ixlen.or.iylen2<iylen1.or.iylen2>iylen) then
          write(11,*)'Wrong reduced indices:',ixlen1,ixlen2,iylen1,iylen2,ixlen,iylen
          stop
        endif
      endif

!
!     At this point all variables have been read, you may proceed with processing.
!
!     zyl: interpolate to CORIE grid
!      lon=lon-360 !convert to our long.

!     Write a test file of build points
      if(ifile==0.or.ifile==ndays) then
        write(99, *)'SSh from NRL',ifile
        write(99, *)(ixlen2-ixlen1+1)*(iylen2-iylen1+1)
        do i=ixlen1,ixlen2
          do j=iylen1,iylen2
            write(99, *) (i-1)*iylen1+j,lon(i,j),lat(i,j),h(i,j),ssh(i,j)
          enddo
        enddo
        write(99, *)'+++++++++++++++++++++++++++'
      endif

!     Find parent elements and levels for hgrid.ll, and do interpolation
      loop4: do i=1,np
        etaout(i)=-9999 !flag
        ifl=0
        do l=1,nob
          if(isbnd(i)==iob(l)) ifl=1
        enddo !l
        if(ifl==0) cycle loop4

        ixy(i,1)=0; ixy(i,2)=0
        do ix=ixlen1,ixlen2-1 
          do iy=iylen1,iylen2-1 
            x1=lon(ix,iy); x2=lon(ix+1,iy); x3=lon(ix+1,iy+1); x4=lon(ix,iy+1)
            y1=lat(ix,iy); y2=lat(ix+1,iy); y3=lat(ix+1,iy+1); y4=lat(ix,iy+1)
            a1=abs(signa(xl(i),x1,x2,yl(i),y1,y2))
            a2=abs(signa(xl(i),x2,x3,yl(i),y2,y3))
            a3=abs(signa(xl(i),x3,x4,yl(i),y3,y4))
            a4=abs(signa(xl(i),x4,x1,yl(i),y4,y1))
            b1=abs(signa(x1,x2,x3,y1,y2,y3))
            b2=abs(signa(x1,x3,x4,y1,y3,y4))
            rat=abs(a1+a2+a3+a4-b1-b2)/(b1+b2)
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
                    exit
                  endif
                else !nodes 1,3,4
                  bb=abs(signa(x1,x3,x4,y1,y3,y4))
                  wild(l)=abs(a3+a4+ap-bb)/bb
                  if(wild(l)<small1*5) then
                    in=2
                    arco(1)=min(1.,a3/bb)
                    arco(2)=min(1.,a4/bb)
                    exit
                  endif
                endif
              enddo !l=1,2
              if(in==0) then
                write(11,*)'Cannot find a triangle:',(wild(l),l=1,2)
                stop
              endif
              arco(3)=max(0.,min(1.,1-arco(1)-arco(2)))

!             Debug
!              if(i==9351) then
!                if(in==1) then
!                  print*, x1,x2,x3
!                  print*, y1,y2,y3
!                else
!                  print*, x1,x3,x4
!                  print*, y1,y3,y4
!                endif
!                print*, arco(1)+arco(2)+arco(3)
!              endif

              mcount=0
              ssh_av=0 !average
              ssh_int=0 !interpolated
              do jj=1,3
                if(jj==1) then
                  id1=ix; id2=iy
                else if(jj==2) then
                  if(in==1) then
                    id1=ix+1; id2=iy
                  else
                    id1=ix+1; id2=iy+1
                  endif
                else !jj=3
                  if(in==1) then
                    id1=ix+1; id2=iy+1
                  else
                    id1=ix; id2=iy+1
                  endif
                endif

                if(h(id1,id2)>0) then
                  mcount=mcount+1
                  ssh_av=ssh_av+ssh(id1,id2)
                  ssh_int=ssh_int+ssh(id1,id2)*arco(jj)
                endif
              enddo !jj
              if(mcount==0) then
                ival(i)=0 !flag
                etaout(i)=-99 
              else if(mcount==3) then
                ival(i)=1
                etaout(i)=ssh_int
              else !mcount/=0
                ival(i)=0
                etaout(i)=ssh_av/mcount
              endif

              cycle loop4
            endif !rat<small1
          enddo !iy=iylen1,iylen2-1
        enddo !ix=ixlen1,ixlen2-1
        if(ixy(i,1)==0.or.ixy(i,2)==0) then
          write(11,*)'Cannot find a parent element:',i
          stop
        endif
      end do loop4 !i=1,np

!     Output 
      nond0=0
      do l=1,nob
        k=iob(l)
        nond0=nond0+nond(k)
        do i=1,nond(k)
          nd=iond(k,i)
          if(etaout(nd)<-9998) then
            write(11,*)'Flag wrong:',nd,k,i
            stop
          endif
        
          z0_out(1,l,i)=etaout(nd)
          ival_out(1,l,i)=ival(nd)

!          if(etaout(nd)<-98) then
!            write(11,*)'Invalid'
!            stop
!          else
!            if(ival(nd)==0) then
!              write(22,*)k,i,etaout(nd),' Invalid'
!            else
!              write(22,*)k,i,etaout(nd)
!            endif
!          endif
        enddo !i
      enddo !l
!      write(22,*)ifile*86400.,nond0,((etaout(iond(iob(l),i)),i=1,nond(iob(l))),l=1,nob)

!     Debug
!      write(98,*)ifile,etaout(2),etaout(5),etaout(15086),etaout(25969)

!      write(20,*)'SSH on hgrid.ll (-99 for invalid pts)',ifile
!      write(20,*)ne,np
!      do i=1,np
!        write(20,*)i,xl(i),yl(i),etaout(i)
!      enddo !i
!      do i=1,ne
!        write(20,*)i,3,(nm(i,l),l=1,3)
!      enddo !i
!      write(20,*)'+++++++++++++++++++++++++'
!-----------------------------------------------------------------------------------------------------------
      enddo !ifile=0,ndays

!     Output
      open(22,file='Z0.out')
      do l=1,nob
        write(22,*)'Z0'
        k=iob(l)
        do i=1,nond(k)
          if(z0_out(1,l,i)<-98) then !.or.z0_out(2,l,i)<-98) then
            write(22,*)'Invalid'
          else if(ival_out(1,l,i)==1) then
            write(22,*)z0_out(1,l,i),0.
          else
            write(22,*)z0_out(1,l,i),0.,' Invalid'
          endif
        enddo !i
      enddo !l

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

