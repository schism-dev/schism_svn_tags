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

! Read sflux/ netcdf files and do analyses (time series). All vars are 2D.
! Run in sflux/ dir.
! Tested so far on NARR and NAM.
!
! Compilation on amb6402:
!
! ifort -O3 -Bstatic -o analyze_sflux1 analyze_sflux1.f90 -Vaxlib  -I/usr/local/netcdf/include/ -L/usr/local/netcdf/lib -lnetcdf
!
!   Input: 
!     (1) screen inputs;
!     (2) .nc files;
!     (3) station.bp: x,y are lon., lat in degrees; z is not used;
!   Output: (1) fort.18 (time series)
!           (2) fort.11: fatal and non-fatal errors;

      program readNCOM
!     netcdf modules from .../include/
      use typeSizes
      use netcdf

      integer, parameter :: nbyte=4
      real, parameter :: small1=1.e-3 !used to check area ratios
  
!     netcdf related variables
      integer :: sid ! Netcdf file IDs
      integer :: svid,uvid ! variable IDs
      integer, dimension(nf90_max_var_dims) :: dids
      real (kind = FourByteReal), dimension(:,:), allocatable :: lat,lon
      allocatable :: varin(:,:,:)

!     File names for netcdf files
      character(len=3) :: day_char
      character(len=5) :: fname
      character(len=18) :: sfile,varname

      integer :: status ! netcdf local status variable
      integer :: ier ! allocate error return.

      allocatable :: xl(:),yl(:),ixy(:,:),varout(:,:),timeout(:)
      dimension arco(4),wild(100),wild2(100,2),month_day(12)

!     First statement
      print*, 'Input nc file name (air_[12], rad_[12] or prc_[12]):'
      read*, fname
!     Possible choices of varname: 
!          In sflux_air*: uwind, vwind, prmsl (pressure), stmp (airt in Kelvin),
!                         spfh (humidity)
!          In sflux_rad*: dlwrf (Downward Long Wave Radiation Flux)
!                         dswrf (Downward Short Wave Radiation Flux i.e. solar)
!          In sflux_prc*: prate (precip), 
      print*, 'Input var. name:'
      read*, varname
      varname=adjustl(varname)
      ivarname=len_trim(varname)
      print*, 'Input # of stacks:'
      read*, ndays
      print*, 'Input offset in hours to GMT (used for output only):'
      read*, offset

!     Read in station.bp
      open(14,file='station.bp',status='old') !only need depth info and connectivity
      read(14,*)
      read(14,*)np
      allocate(xl(np),yl(np),ixy(np,2),stat=istat)
      if(istat/=0) stop 'Failed to allocate (1)'
      do i=1,np
        read(14,*)j,xl(i),yl(i)
      enddo !i
      close(14)

!     Loop over all days
      do iday=1,ndays
!------------------------------------------------------------------------------------------------
      day_char='000'
      write(day_char(1:3),'(i3.3)') iday

!     Create nc file names for this date.
      sfile='sflux_'//fname//'.'//day_char//'.nc'
      write(*,*)'doing ',sfile

      status = nf90_open(trim(sfile), nf90_nowrite, sid)
      call check(status)

!     Get index information 
      status = nf90_inq_varid(sid, varname(1:ivarname), svid) !var. name is after, say, 'float'
      call check(status)
      status = nf90_Inquire_Variable(sid, svid, dimids = dids)
      call check(status)
      status = nf90_Inquire_Dimension(sid, dids(1), len = nx) !last index in nc file
      call check(status)
      status = nf90_Inquire_Dimension(sid, dids(2), len = ny)
      call check(status)
      status = nf90_Inquire_Dimension(sid, dids(3), len = nstep)
      call check(status)

!     allocate memory.
      allocate(lat(nx,ny),lon(nx,ny),varin(nx,ny,nstep),varout(np,nstep), &
     &timeout(nstep),stat=istat)
      if(istat/=0) stop 'Failed to allocate (2)'

!     get the index values in all directions
      status = nf90_get_var(sid,svid,varin,start=(/1,1,1/),count=(/nx, ny, nstep/),stride=(/1,1,1/))
      call check(status)
      status = nf90_inq_varid(sid, "lon", uvid)
      call check(status)
      status = nf90_get_var(sid,uvid,lon,start=(/1,1/),count=(/nx, ny/),stride=(/1,1/))
      call check(status)
      status = nf90_inq_varid(sid, "lat", uvid)
      call check(status)
      status = nf90_get_var(sid,uvid,lat,start=(/1,1/),count=(/nx, ny/),stride=(/1,1/))
      call check(status)
      !time in fraction of days
      status = nf90_inq_varid(sid, "time", uvid)
      call check(status)
      status = nf90_get_var(sid,uvid,timeout,start=(/1/),count=(/nstep/),stride=(/1/))
      call check(status)
      timeout=timeout+iday-1-offset/24

!     At this point all variables have been read, you may proceed with processing.
!      lon=lon-360 !convert to our long.

!     Write a test file of build points
!      write(15, *)'Test NCOM output.'
!      write(15, *)nx*ny
!      do i=1,nx
!        do j=1,ny
!         write(15, *) (i-1)*ny+j,lon(i,j),lat(i,j),varin(i,j,1)
!        enddo
!      enddo
!      stop

!     Find parent elements and levels for hgrid.ll, and do interpolation
      varout=-99
      ixy=0
      loop4: do ix=1,nx-1 
        do iy=1,ny-1 
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
                write(11,*)'Sum of coordinates/=1:',arco(1:4),sum(arco)
                stop
              endif

              varout(i,:)=varin(ix,iy,:)*arco(1)+varin(ix+1,iy,:)*arco(2)+ &
     &varin(ix+1,iy+1,:)*arco(3)+varin(ix,iy+1,:)*arco(4)
            endif !rat<small1
          enddo !i=1,np
          if(ialldone==1) exit loop4 
        enddo !iy
      enddo loop4 !ix

      do i=1,np
        if(ixy(i,1)==0.or.ixy(i,2)==0)  &
     &write(*,*)'Cannot find a parent element for pt:',i,ixy(i,1:2)
      enddo !i

      do j=1,nstep
        write(18,'(e13.6,6000(1x,f12.3))')timeout(j),(varout(i,j),i=1,np)
      enddo !i

      deallocate(lat,lon,varin,varout,timeout)
!-----------------------------------------------------------------------------------------------------------
      enddo !iday

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

      end program readNCOM

!      function signa(x1,x2,x3,y1,y2,y3)
!
!      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2
!
!      return
!      end

