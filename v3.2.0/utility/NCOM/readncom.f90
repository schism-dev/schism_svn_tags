! Read NRL/NCOM netcdf files
!
! pturner 12-2004
!
! X = longitude dimension is sampled in the salt and temp data files.
! Y = latitude dimension is sampled in the salt and temp data files.
! The level dimensions goes from 1 = surface, 40 = bottom.
!
! Compilation on amb1001:
!
! ifc -g -O3 -Bstatic -o readncom readncom.f90 -Vaxlib -L/usr/local/netcdf/lib -lnetcdf
!
! Need to have the following modules symlinked into the current directory
!
! NETCDF.mod -> /usr/local/netcdf/include/NETCDF.mod
! TYPESIZES.mod -> /usr/local/netcdf/include/TYPESIZES.mod
!

!zyl
!   Input: hgrid.gr3 (only for depth), hgrid.ll, vgrid.in (SELFE); estuary.gr3 (flags)
program readNCOM
! netcdf modules from symlinks
  use typeSizes
  use netcdf

!  implicit none

  integer, parameter :: debug=1
  integer, parameter :: mnp=60000
  integer, parameter :: mne=100000
  integer, parameter :: mnv=60
  real, parameter :: small1=1.e-3 !used to check area ratios
  
! netcdf related variables
  integer :: hid, latid, lonid, zmid, sid, tid ! Netcdf file IDs
  integer :: latvid, lonvid, zmvid, hvid ! positional variables
  integer :: xdid, xvid ! longitude index
  integer :: ydid, yvid ! latitude index
  integer :: ldid, lvid ! vertical level, 1 is top
  integer :: svid, tvid ! salt & temp variable IDs
  integer, dimension(nf90_max_var_dims) :: dids
             
! Local variables for data
  integer, dimension(:), allocatable :: xind, yind, lind ! sampled indexes
! Lat, lon, bathymetry
  real (kind = FourByteReal), dimension(:,:), allocatable :: lat, lon, h
! Vertical postion, salinity, and temperature
  real (kind = FourByteReal), dimension(:,:,:), allocatable :: zm, salt, temp
  integer, dimension(:,:), allocatable :: kbp,ihope
! File names for netcdf files
  character(len=1024) :: hnc, latnc, lonnc, zmnc, sfile, tfile
! Command line arguments
  character(len=1024) :: s, yr, mdh
! external function for number of command line arguments
  integer :: iargc

  integer :: status ! netcdf local status variable
  integer :: ier ! allocate error return.
  integer :: ixlen, iylen, ilen ! sampled lengths in each coordinate direction
  integer :: ixlen1, iylen1,ixlen2, iylen2 !reduced indices for CORIE grid to speed up interpolation
!  integer :: i,j,k,i1,i2,i3,j1,j2,j3

  dimension xl(mnp),yl(mnp),nm(mne,3),dp(mnp)
  dimension sigma(mnv),cs(mnv),z(mnp,mnv),iest(mnp),ixy(mnp,2),arco(3)
  dimension wild(100),wild2(100,2)
  dimension tempout(mnp,mnv), saltout(mnp,mnv)

! First statement
! A check to netcdf library to make sure all is OK. Not needed 
! but left here anyway.
  if(.not. byteSizesOK()) then
    print *, "Compiler does not appear to support required kinds of variables."
    stop
  end if

!
! geometry files - should never change.
!
  hnc = '/home/workspace/ccalmr6/nrldata/model_h.nc'
  latnc = '/home/workspace/ccalmr6/nrldata/model_lat.nc'
  lonnc = '/home/workspace/ccalmr6/nrldata/model_lon.nc'
  zmnc = '/home/workspace/ccalmr6/nrldata/model_zm.nc'
  
  if (iargc() /= 2) then
    call getarg(0,s)
    write(*,*) 'Usage: ', trim(s), ' [YYYY] [MMDDHH]' 
    write(*,*) 'Example: ', trim(s), ' 2004 122400' 
    stop
  end if

! Get year and date from the command line
  call getarg(1, yr)
  call getarg(2, mdh)

! Create salt and temp file names for this date.
! Compile spits out a warning here, dunno how to fix.
  sfile = '/home/workspace/ccalmr6/nrldata/' // trim(yr) // '/s3d/s3d.glb8_2f_' // trim(yr) // trim(mdh) // '.nc'
  tfile = '/home/workspace/ccalmr6/nrldata/' // trim(yr) // '/t3d/t3d.glb8_2f_' // trim(yr) // trim(mdh) // '.nc'

!
! Open files for salinity and temperature.
  status = nf90_open(trim(sfile), nf90_nowrite, sid)
  call check(status)
  status = nf90_open(trim(tfile), nf90_nowrite, tid)
  call check(status)

! Get index information for sampled grid, assumed the same for temperature.
! These will index into the static lat, lon, bathymetry, and zm files.
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

! get the variable ids from salt and temperature files
  status = nf90_inq_varid(sid, "Salinity", svid)
  call check(status)
  status = nf90_inq_varid(tid, "Temperature", tvid)
  call check(status)

! Get the lengths of the dimensions from the salt file.
! Assumed same as in temperature
  status = nf90_Inquire_Variable(sid, svid, dimids = dids)
  call check(status)
  status = nf90_Inquire_Dimension(sid, dids(1), len = ixlen)
  call check(status)
  status = nf90_Inquire_Dimension(sid, dids(2), len = iylen)
  call check(status)
  status = nf90_Inquire_Dimension(sid, dids(3), len = ilen)
  call check(status)

! allocate memory.
  allocate( xind(ixlen),stat=ier)
  allocate( yind(iylen),stat=ier)
  allocate( lind(ilen),stat=ier)
  allocate( salt(ixlen, iylen, ilen),stat=ier)
  if (ier /= 0) then
    write(*,*) ' Could not allocate salt'
    stop 'Allocate salt'
  endif
  allocate( temp(ixlen, iylen, ilen),stat=ier)
  if (ier /= 0) then
    write(*,*) ' Could not allocate temp'
    stop 'Allocate temp'
  endif
  
! get the index values in all directions
  status = nf90_get_var(sid, xvid, xind)
  call check(status)
  status = nf90_get_var(sid, yvid, yind)
  call check(status)
  status = nf90_get_var(sid, lvid, lind)
  call check(status)
! Read salinity
  status = nf90_get_var(sid,svid,salt,start=(/1,1,1/),count=(/ixlen, iylen, ilen/),stride=(/1,1,1/))
  call check(status)
! Read temperature
  status = nf90_get_var(tid,tvid,temp,start=(/1,1,1/),count=(/ixlen, iylen, ilen/),stride=(/1,1,1/))
  call check(status)

! TODO  call read3dvar(sfile, 'Salinity', 1, ixlen, 1, iylen, 1, ilen, salt)
! TODO  call read3dvar(tfile, 'Temperature', 1, ixlen, 1, iylen, 1, ilen, salt)

  status = nf90_close(sid)
  call check(status)
  status = nf90_close(tid)
  call check(status)

! Extracting a subset of all lat,lon,h,vertical pos data
  allocate( lat(ixlen,iylen))
  allocate( lon(ixlen,iylen))
  allocate( zm(ixlen, iylen, ilen))
  allocate( h(ixlen,iylen))
  allocate( kbp(ixlen,iylen))
  allocate( ihope(ixlen,iylen))

  status = nf90_open(trim(hnc), nf90_nowrite, hid)
  call check(status)
  status = nf90_inq_varid(hid, "bathymetry", hvid)
  call check(status)
  status = nf90_get_var(hid,hvid,h,start=(/xind(1),yind(1)/),count=(/ixlen, iylen/),stride=(/1,1/))
  call check(status)
  status = nf90_close(hid)
! TODO  call read2dvar(hnc, 'bathymetry', xind(1), ixlen, yind(1), iylen, h)

  status = nf90_open(trim(latnc), nf90_nowrite, latid)
  call check(status)
  status = nf90_inq_varid(latid, "Lat", latvid)
  call check(status)
  status = nf90_get_var(latid,latvid,lat,start=(/xind(1),yind(1)/),count=(/ixlen, iylen/),stride=(/1,1/))
  call check(status)
  status = nf90_close(latid)
! TODO  call read2dvar(latnc, 'Lat', xind(1), ixlen, yind(1), iylen, lat)

  status = nf90_open(trim(lonnc), nf90_nowrite, lonid)
  call check(status)
  status = nf90_inq_varid(lonid, "Long", lonvid)
  call check(status)
  status = nf90_get_var(lonid,lonvid,lon,start=(/xind(1),yind(1)/),count=(/ixlen, iylen/),stride=(/1,1/))
  call check(status)
  status = nf90_close(lonid)
! TODO  call read2dvar(lonnc, 'Long', xind(1), ixlen, yind(1), iylen, lon)

! Read level depths information
  status = nf90_open(trim(zmnc), nf90_nowrite, zmid)
  call check(status)
  status = nf90_inq_varid(zmid, "zm", zmvid)
  call check(status)
  status = nf90_get_var(zmid,zmvid,zm,start=(/xind(1),yind(1),lind(1)/),count=(/ixlen, iylen, ilen/),stride=(/1,1,1/))
  call check(status)
  status = nf90_close(zmid)
! TODO  call read3dvar(zmnc, 'zm', xind(1), ixlen, yind(1), iylen, lind(1), ilen, zm)

! Reduce indices ixlen, iylen for CORIE grid to speed up inetrpolation
  ixlen1=49
  ixlen2=80
  iylen1=279
  iylen2=370
  if(ixlen2<ixlen1.or.ixlen2>ixlen.or.iylen2<iylen1.or.iylen2>iylen) then
    write(11,*)'Wrong reduced indices:',ixlen1,ixlen2,iylen1,iylen2,ixlen,iylen
    stop
  endif

!
! At this point all variables have been read, you may proceed with processing.
!
!zyl: interpolate to CORIE grid
  lon=lon-360 !convert to our long.
!  print*, ixlen,iylen,ilen
!Compute bottom indices and extend S,T data below bottom
  do i=ixlen1,ixlen2
    do j=iylen1,iylen2
      kbp(i,j)=-99
      do k=ilen,1,-1
        if(zm(i,j,k)>-1.e20) then
          kbp(i,j)=k
          exit
        endif
      enddo !k

      if(kbp(i,j)/=-99) then !valid pts
        do k=1,ilen
          if(k<=kbp(i,j)) then
            if(salt(i,j,k)<-99.or.temp(i,j,k)<-99) then
              write(11,*)'Fatal: no valid S,T:',i,j,k,salt(i,j,k),temp(i,j,k)
              stop
            endif
          else !extend
            salt(i,j,k)=salt(i,j,kbp(i,j))
            temp(i,j,k)=temp(i,j,kbp(i,j))
          endif
        enddo !k
      endif

!      write(12,*)(i-1)*iylen + j, lon(i,j), lat(i,j), kbp(i,j)
    enddo !j
  enddo !i  

! Write a test file of build points
  if (debug == 1) then
      write(15, *)'Test NCOM output.'
      write(15, *)(ixlen2-ixlen1+1)*(iylen2-iylen1+1) !ixlen * iylen
      do i=ixlen1,ixlen2 !1,ixlen
        do j=iylen1,iylen2 !1,iylen
! top layer is 1, ilen is bottom layer.
!           write(15, *) (i-1)*iylen+j,lon(i,j),lat(i,j),salt(i,j,max(kbp(i,j),1)) !-h(i,j) 
           write(15, *) (i-1)*iylen+j,lon(i,j),lat(i,j),temp(i,j,2) !-h(i,j) 
        enddo
      enddo
   endif

!Coompute S,T @ invalid pts based on nearest neighbor
  ihope=1
  do i=ixlen1,ixlen2
    do j=iylen1,iylen2
      if(kbp(i,j)==-99) then !invalid pts  
!       Search along x-direction first
        i1=-99999
        do k=i-1,ixlen1,-1
          if(kbp(k,j)/=-99) then
            i1=k
            exit
          endif
        enddo !k
        i2=-99999
        do k=i+1,ixlen2
          if(kbp(k,j)/=-99) then
            i2=k
            exit
          endif
        enddo !k
        if(iabs(i-i1)>=iabs(i-i2)) then
          i3=i2 !could be -99999
        else
          i3=i1
        endif

!       Search along y-direction 
        j1=-99999
        do k=j-1,iylen1,-1
          if(kbp(i,k)/=-99) then
            j1=k
            exit
          endif
        enddo !k
        j2=-99999
        do k=j+1,iylen2
          if(kbp(i,k)/=-99) then
            j2=k
            exit
          endif
        enddo !k
        if(iabs(j-j1)>=iabs(j-j2)) then
          j3=j2 !could be -99999
        else
          j3=j1
        endif

        if(iabs(i-i3)>=iabs(j-j3)) then
          i1=i; j1=j3
        else
          i1=i3; j1=j
        endif

        if(i1==-99999.or.j1==-99999) then !hopeless pts
          ihope(i,j)=0
        else
!          ihope(i,j)=1
          if(i1<ixlen1.or.i1>ixlen2.or.j1<iylen1.or.j1>iylen2) then
            write(11,*)'Indices out of bound (1):',i,j,i1,j1
            stop
          endif
          salt(i,j,1:ilen)=salt(i1,j1,1:ilen)
          temp(i,j,1:ilen)=temp(i1,j1,1:ilen)
        endif
      endif !kbp(i,j)==-99

!      write(13, *)(i-1)*iylen+j,lon(i,j),lat(i,j),salt(i,j,ilen) !-h(i,j) 
!      write(15, *) (i-1)*iylen + j, lon(i,j), lat(i,j), ihope(i,j) !-h(i,j) 
    enddo !j=iylen1,iylen2
  enddo !i=ixlen1,ixlen2

! Read in hgrid and vgrid
  open(17,file='estuary.gr3',status='old')
  open(16,file='hgrid.ll',status='old')
  open(14,file='hgrid.gr3',status='old') !only need depth info
  open(19,file='vgrid.in',status='old')
  read(14,*)
  read(14,*)ne,np
  if(np.gt.mnp.or.ne.gt.mne) then
    write(*,*)'Increase mnp/mne'
    stop
  endif
  read(16,*)
  read(16,*)
  read(17,*)
  read(17,*)
  icount=0 !# of anchor pts for transition (should be 1 or 2)
  icount2=0
  ianchor2=-99
  do i=1,np
    read(14,*)j,xtmp,ytmp,dp(i)
    read(16,*)j,xl(i),yl(i) !,dp(i)
    read(17,*)j,xtmp,ytmp,iest(i)
    if(abs(iest(i))>1.and.iest(i)/=-2) then
      write(11,*)'Estuary flag wrong:',i,iest(i)
      stop
    endif
    if(iest(i)==-1) then
      icount=icount+1
      ianchor1=i
    endif
    if(iest(i)==-2) then
      icount2=icount2+1
      ianchor2=i
    endif
  enddo !i
  if(icount/=1) then
    write(11,*)'# of anchor pts/=1:',icount
    stop
  endif
  if(icount2>1) then
    write(11,*)'# of anchor pts>1:',icount2
    stop
  endif
  nanchor=icount+icount2
  if(nanchor==2.and.ianchor2==-99) then
    write(11,*)'Impossible 1'
    stop
  endif
!  do i=1,ne
!    read(14,*)j,k,(nm(i,l),l=1,3)
!  enddo !i
  close(14)
  close(16)

  read(19,*) ivcor
  if(ivcor==1) then !traditional sigma
    h_c=0
  else if(ivcor==2) then !s1
    read(19,*)h_c,theta_b,theta_f
    if(h_c<0) then
      write(11,*)'Negative h_c:',h_c
      stop
    endif
    if(theta_b<0.or.theta_b>1) then
      write(11,*)'Wrong theta_b:',theta_b
      stop
    endif
    if(theta_f<=0.or.theta_f>20) then
      write(11,*)'Wrong theta_f:',theta_f 
      stop
    endif
  else
    write(11,*)'Unknown sigma:',ivcor
    stop
  endif

  read(19,*) nvrt
  if(nvrt>mnv.or.nvrt<3) then
    write(11,*)'nvrt > mnv or nvrt<3'
    stop
  endif

  sigma(1)=-1 !bottom
  sigma(nvrt)=0 !surface
  do k=2,nvrt-1
    read(19,*) j,sigma(k)
    if(sigma(k)<=sigma(k-1).or.sigma(k)>=0) then
      write(11,*)'Check sigma levels at:',k
      stop
    endif
  enddo
  close(19)

! Compute C(s)
  do k=1,nvrt
    if(ivcor==1) then
      cs(k)=sigma(k)
    else if(ivcor==2) then
      cs(k)=(1-theta_b)*sinh(theta_f*sigma(k))/sinh(theta_f)+ &
     &theta_b*(tanh(theta_f*(sigma(k)+0.5))-tanh(theta_f*0.5))/2/tanh(theta_f*0.5)
    endif
  enddo !k=1,nvrt

  do i=1,np
    do k=1,nvrt
      if(dp(i)<=h_c) then
        z(i,k)=max(dp(i),0.1)*sigma(k)
      else
        z(i,k)=h_c*sigma(k)+(dp(i)-h_c)*cs(k)
      endif
    enddo !k
  enddo !i

! Find parent elements and levels for hgrid.ll, and do interpolation
  tempout=-99; saltout=-99
  loop4: do i=1,np
    if(iest(i)==1.or.iest(i)==-2) cycle loop4

!   Non-estuary nodes
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
!        if((ix-1)*iylen+iy==27837) then
!          print*, rat,a1+a2+a3+a4,b1+b2
!          print*, x1,x2,x3,x4
!          print*, y1,y2,y3,y4
!          print*, xl(i),yl(i)
!        endif
        if(rat<small1) then
          ixy(i,1)=ix; ixy(i,2)=iy
!         Find a triangle
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

!         Debug
!          if(i==9351) then
!            if(in==1) then
!              print*, x1,x2,x3
!              print*, y1,y2,y3
!            else
!              print*, x1,x3,x4
!              print*, y1,y3,y4
!            endif
!            print*, arco(1)+arco(2)+arco(3)
!          endif

          if(ihope(ix,iy)==0.or.ihope(ix+1,iy)==0.or.ihope(ix+1,iy+1)==0.or.ihope(ix,iy+1)==0) then
            write(11,*)'Hopeless parent:',i,ix,iy
            stop
          else !do interpolation
!         Find vertical level
            do k=1,nvrt
              if(kbp(ix,iy)==-99) then
                lev=ilen-1; vrat=1
              else if(z(i,k)<zm(ix,iy,kbp(ix,iy))) then
                lev=kbp(ix,iy)-1; vrat=1
              else if(z(i,k)>zm(ix,iy,1)) then
                lev=1; vrat=0
              else
                lev=-99 !flag
                do kk=1,kbp(ix,iy)-1
                  if(z(i,k)<=zm(ix,iy,kk).and.z(i,k)>=zm(ix,iy,kk+1)) then
                    lev=kk
                    vrat=(zm(ix,iy,kk)-z(i,k))/(zm(ix,iy,kk)-zm(ix,iy,kk+1))
                    exit
                  endif
                enddo !kk
                if(lev==-99) then
                  write(11,*)'Cannot find a level:',i,k,z(i,k),(zm(ix,iy,l),l=1,kbp(ix,iy))
                  stop
                endif
              endif
           
!             write(18,*)i,k,ix,iy,lev,vrat,kbp(ix,iy)
              wild2(1,1)=temp(ix,iy,lev)*(1-vrat)+temp(ix,iy,lev+1)*vrat
              wild2(1,2)=salt(ix,iy,lev)*(1-vrat)+salt(ix,iy,lev+1)*vrat
              wild2(2,1)=temp(ix+1,iy,lev)*(1-vrat)+temp(ix+1,iy,lev+1)*vrat
              wild2(2,2)=salt(ix+1,iy,lev)*(1-vrat)+salt(ix+1,iy,lev+1)*vrat
              wild2(3,1)=temp(ix+1,iy+1,lev)*(1-vrat)+temp(ix+1,iy+1,lev+1)*vrat
              wild2(3,2)=salt(ix+1,iy+1,lev)*(1-vrat)+salt(ix+1,iy+1,lev+1)*vrat
              wild2(4,1)=temp(ix,iy+1,lev)*(1-vrat)+temp(ix,iy+1,lev+1)*vrat
              wild2(4,2)=salt(ix,iy+1,lev)*(1-vrat)+salt(ix,iy+1,lev+1)*vrat
              if(in==1) then
                tempout(i,k)=wild2(1,1)*arco(1)+wild2(2,1)*arco(2)+wild2(3,1)*arco(3)
                saltout(i,k)=wild2(1,2)*arco(1)+wild2(2,2)*arco(2)+wild2(3,2)*arco(3)
              else
                tempout(i,k)=wild2(1,1)*arco(1)+wild2(3,1)*arco(2)+wild2(4,1)*arco(3)
                saltout(i,k)=wild2(1,2)*arco(1)+wild2(3,2)*arco(2)+wild2(4,2)*arco(3)
              endif

!             enforce lower bound for salt
              if(z(i,k)<-4) saltout(i,k)=max(saltout(i,k),33.)
            enddo !k=1,nvrt
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

! Estuary nodes
  if(tempout(ianchor1,1)==-99) then
    write(11,*)'Anchor node not assigned'
    stop
  endif
  do i=1,np
    if(iest(i)==1.or.iest(i)==-2) then
    
      tempout(i,1:nvrt)=tempout(ianchor1,1:nvrt)
      if(nanchor==1) then
        saltout(i,1:nvrt)=saltout(ianchor1,1:nvrt)
      else !=2
        if(xl(ianchor2)-xl(ianchor1)==0) then
          write(11,*)'Wrong anchor pts:',ianchor1,ianchor2
          stop
        endif
        xrat=(xl(i)-xl(ianchor1))/(xl(ianchor2)-xl(ianchor1))
        xrat=max(0.,min(1.,xrat))
        do k=1,nvrt
          saltout(i,k)=saltout(ianchor1,k)*(1-xrat) !+0
        enddo !k 
      endif
    endif !inside estuary
  enddo !i
    
  do i=1,np
    write(20,*)i,xl(i),yl(i),tempout(i,1)
    write(21,*)i,xl(i),yl(i),tempout(i,nvrt)
  enddo !i

  print*, 'Finished'

! End of main
! Subroutines
contains
  ! Internal subroutine - checks error status after each netcdf, 
  ! prints out text message each time an error code is returned. 
  subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
    end if
  end subroutine check  

  ! A simpler interface for reading 1, 2, 3d variables
  ! Not activated right now.
  !
  ! fname = Netcdf filename
  ! vname = variable name
  ! nix = starting index to read >= 1
  ! n = number to read
  ! var = return value
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

  ! fname = Netcdf filename
  ! vname = variable name
  ! nix = starting index to read >= 1
  ! n = number to read
  ! mix = starting index to read >= 1
  ! m = number to read
  ! var = return value
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

  ! fname = Netcdf filename
  ! vname = variable name
  ! nix = starting index to read >= 1
  ! n = number to read
  ! mix = starting index to read >= 1
  ! m = number to read
  ! lix = starting index to read >= 1
  ! l = number to read
  ! var = return value
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

