! Read NRL/NCOM netcdf files
!
! pturner 12-2004
!
! X = longitude dimension is sampled in the salt and temp data files.
! Y = latitude dimension is sampled in the salt and temp data files.
! The level dimensions goes from 1 = surface, 40 = bottom.
!
! Compilation on amb1004 or amb6402:
! ifort -g -O3 -Bstatic -o readncom8c_intel readncom8c.f90 -Vaxlib -I/usr/local/netcdf/include -L/usr/local/netcdf/lib -lnetcdf
!
! Compilation on canopus01:
! ifort -g -O3 -Bstatic -o readncom8c_canopus readncom8c.f90 -Vaxlib -I/usr/local/include -L/usr/local/lib -lnetcdf

!zyl
!   Input: 
!     (1) hgrid.gr3;
!     (2) hgrid.ll;
!     (3) vgrid.in (SELFE or ELCIRC test02k4);
!     (4) estuary.gr3 (flags): depth=0: outside; =1: inside; =-1: 1st anchor; =-2: 2nd anchor (optional)
!     (5) date.in: 1st line: ivartype (1: T,S; 2; u,v), nob, (iob(i),i=1,nob) (# of open bnd segment 
!                            that has ifltype=4 or -4; segment # follow bctides.in; not used for T,S);
!                  2nd line: year, month, day, ndays (# of days); 
!                  (note: input "0" in month or day wouldn't cause any difference); 
!                  3rd line: ireduce; 
!                  4th line: "SELFE" or "ELCIRC"; 
!                  5th line: depth (ht), Smin to be imposed when z <= -ht (for S,T only);
!                  6th line: ntracers (note: tracer conc. =0 uniformly);
!                  7th line: h_mean, xl_mean,yl_mean. Depth (h_mean) below which a profile 
!                            is used to generate T,S for hotstart.in only; xl_mean and yl_mean 
!                            are used to find a nearest node (in NCOM grid) for this profile.
!     (6) temp.th: needed if there are 2 anchors. Same format as CORIE runs. Only first river temp.
!                  on the first bnd is used
!   Output: hotstart.in (itur=3 for ELCIRC; unformatted for MPI SELFE) when ndays=0 together
!           with ts.ic.0 (mean profile of T,S);
!           or salt_nu.in and temp_nu.in (unformatted) when ndays>0 and ivartype=1;
!           or uv_bcc.th (ASCII;  nodes on bnds nob; use gen_uv3D.f90 to generate uv3D.th)
!           when ndays>0 and ivartype=2.
!           Note that the unformatted outputs are different between Intel and AMD!!!

      program readNCOM
!     netcdf modules from symlinks
      use typeSizes
      use netcdf

!      implicit none

!      integer, parameter :: debug=1
      integer, parameter :: mnp=70000
      integer, parameter :: mne=140000
      integer, parameter :: mns=210000
      integer, parameter :: mnv=60
      integer, parameter :: mnope=10 !max # of open bnd segments
      integer, parameter :: mnond=1000 !max # of open bnd nodes in each segment
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
             
!     Local variables for data
      integer, dimension(:), allocatable :: xind, yind, lind ! sampled indexes
!     Lat, lon, bathymetry
      real (kind = FourByteReal), dimension(:,:), allocatable :: lat, lon, h
!     Vertical postion, salinity, and temperature
      real (kind = FourByteReal), dimension(:,:,:), allocatable :: zm, salt, temp
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
      integer :: elnode(4,mne),ic3(4,mne),elside(4,mne)
      integer :: isdel(2,mns)
      dimension xl(mnp),yl(mnp),dp(mnp),i34(mne),kbp2(mnp)
      dimension ztot(0:mnv),sigma(mnv),cs(mnv),z(mnp,mnv),iest(mnp),ixy(mnp,2),arco(3)
      dimension wild(100),wild2(100,2),swild(mnv,2),swild2(mnv),zs(mnv),ze(mnv)
      dimension tempout(mnp,mnv), saltout(mnp,mnv),month_day(12)
      dimension tsd(mns,mnv),ssd(mns,mnv),tsel(mnv,mne,2)
      dimension nne(mnp),indel(mnei,mnp),nx(4,4,3),isidenode(2,mns)
      dimension xcj(mns),ycj(mns),nond(mnope),iond(mnope,mnond),iob(mnope),iond2(mnope*mnond)
      allocatable :: z_r(:),temp_bar(:),salt_bar(:),indx_var(:),cspline_ypp(:,:)

!     First statement
      itur=3 !for ELCIRC only
      open(10,file='date.in',status='old')
!     ivartype=1: T&S; 2: u,v; if 2, then nob,(iob(i),i=1,nob) are the open bnd segment # (not used for T,S)
      read(10,*) ivartype,nob,(iob(i),i=1,nob) 
      if(ivartype/=1.and.ivartype/=2) stop 'Check ivartype in date.in'
      if(nob>mnope) stop 'Increase mnope'
      read(10,*) iyear,month0,iday0,ndays
      if(ndays==0.and.ivartype==2) stop 'ivartype=2 and ndays==0'
!     ireduce=0: original set; =1: reduced bg set for quicker interpolation (edit ixlen[1,2],iylen[1,2])
      read(10,*) ireduce 
      read(10,*) imodel !SELFE or ELCIRC
      if(imodel.ne."S".and.imodel.ne."E") then
        print*, 'Unknown model:',imodel
        stop
      endif
      read(10,*) ht,smin
      read(10,*) ntracers
      read(10,*) h_mean,xl_mean,yl_mean
      close(10)

!     Read in hgrid and vgrid
      open(17,file='estuary.gr3',status='old')
      open(16,file='hgrid.ll',status='old')
      open(14,file='hgrid.gr3',status='old') !only need depth info and connectivity
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
        if(iest(i)<-2.or.iest(i)>1) then
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
      do i=1,ne
        read(14,*)j,i34(i),(elnode(l,i),l=1,i34(i))
      enddo !i
!     Open bnds
      read(14,*) nope
      read(14,*) neta
      ntot=0
      if(nope>mnope) stop 'Increase mnope (2)'
      do k=1,nope
        read(14,*) nond(k)
        if(nond(k)>mnond) stop 'Increase mnond'
        do i=1,nond(k)
          read(14,*) iond(k,i)
        enddo
      enddo

      nond0=0
      do i=1,nob
        ibnd=iob(i)
        do j=1,nond(ibnd)
          nond0=nond0+1
          iond2(nond0)=iond(ibnd,j)
        enddo !j
      enddo !i

      close(14)
      close(16)

!     Calculate the extent of hgrid.ll for averaging NCOM for hotstart.in
      xl_max=maxval(xl); xl_min=minval(xl)
      yl_max=maxval(yl); yl_min=minval(yl)

!     Open temp.th for 2 anchors
!     Since the nudging is not done in estuary, only the river temp. at the beginning is read in
      if(nanchor==2) then
        open(9,file='temp.th',status='old')
        read(9,*)dt_th,tempth
        close(9)
      endif

!     V-grid
      if(imodel.eq."S") then !SELFE
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

!       # of z-levels excluding "bottom" at h_s
        read(19,*) !for adding comment "Z levels"
        do k=1,kz-1
          read(19,*)j,ztot(k)
          if(k>1.and.ztot(k)<=ztot(k-1).or.ztot(k)>=-h_s) then
            write(11,*)'z-level inverted:',k
            stop
          endif
        enddo !k
        read(19,*) !level kz       
!       In case kz=1, there is only 1 ztot(1)=-h_s
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
!       Pre-compute constants
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

!       Compute C(s) and C'(s)
        do k=1,nsig
          cs(k)=(1-theta_b)*sinh(theta_f*sigma(k))/sinh(theta_f)+ &
     &theta_b*(tanh(theta_f*(sigma(k)+0.5))-tanh(theta_f*0.5))/2/tanh(theta_f*0.5)
        enddo !k=1,nvrt

!       Z-coordinates
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

!         Z-levels
          if(dp(i)<=h_s) then
            kbp2(i)=kz
          else !has z layers
            kbp2(i)=0 !flag
            do k=1,kz-1
              if(-dp(i)>=ztot(k).and.-dp(i)<ztot(k+1)) then
                kbp2(i)=k
                exit
              endif
            enddo !k
            if(kbp2(i)==0) then
              write(11,*)'Cannot find a bottom level for node (3):',i
              stop
            endif

          do k=kbp2(i)+1,kz-1
            z(i,k)=ztot(k)
          enddo !k
          z(i,kbp2(i))=-dp(i)
          endif !dp(i)
        enddo !i
      else !ELCIRC
        read(19,*) nvrt,zmsl
        ztot(0)=-zmsl
        do k=1,nvrt
          read(19,*) i,tmp,ztot(k)
          ztot(k)=ztot(k)-zmsl
        enddo !k
        close(19)
       
        do i=1,np
          kbot=0; kf=0
          do k=0,nvrt-1
            if(-dp(i)>=ztot(k).and.-dp(i)<ztot(k+1)) then
              kbot=k+1
              exit
            endif
          enddo !k
          do k=0,nvrt-1
            if(0>ztot(k).and.0<=ztot(k+1)) then
              kf=k+1
              exit
            endif
          enddo !k
          if(kbot==0.or.kf==0) then
            write(11,*)'Cannot find surface/bottom:',i,kbot,kf
            stop
          endif

          do k=1,nvrt
            if(k<kbot.or.k>kf) then
              z(i,k)=ztot(k)
            else if(k==kbot) then
              z(i,k)=(ztot(k)-dp(i))/2
            else if(k==kf) then
              z(i,k)=ztot(k-1)/2
            else !normal
              z(i,k)=(ztot(k-1)+ztot(k))/2
            endif
          enddo !k
        enddo !i
      endif !SELFE or ELCIRC

!     Compute geometry
      do k=3,4
        do i=1,k
          do j=1,k-1
            nx(k,i,j)=i+j
            if(nx(k,i,j)>k) nx(k,i,j)=nx(k,i,j)-k
            if(nx(k,i,j)<1.or.nx(k,i,j)>k) then
              write(*,*)'nx wrong',i,j,k,nx(k,i,j)
              stop
            endif
          enddo !j
        enddo !i
      enddo !k

      do i=1,np
        nne(i)=0
      enddo

      do i=1,ne
        do j=1,i34(i)
          nd=elnode(j,i)
          nne(nd)=nne(nd)+1
          if(nne(nd)>mnei) then
            write(11,*)'Too many neighbors',nd
            stop
          endif
          indel(nne(nd),nd)=i
        enddo
      enddo

!     Compute ball info; this won't be affected by re-arrangement below
      do i=1,ne
        do j=1,i34(i)
          ic3(j,i)=0 !index for bnd sides
          nd1=elnode(nx(i34(i),j,1),i)
          nd2=elnode(nx(i34(i),j,2),i)
          do k=1,nne(nd1)
            ie=indel(k,nd1)
            if(ie/=i.and.(elnode(1,ie)==nd2.or.elnode(2,ie)==nd2.or.elnode(3,ie)==nd2.or.(i34(ie)==4.and.elnode(4,ie)==nd2))) ic3(j,i)=ie
          enddo !k
        enddo !j
      enddo !i

      ns=0 !# of sides
      do i=1,ne
        do j=1,i34(i)
          nd1=elnode(nx(i34(i),j,1),i)
          nd2=elnode(nx(i34(i),j,2),i)
          if(ic3(j,i)==0.or.i<ic3(j,i)) then !new sides
            ns=ns+1
            if(ns>mns) then
              write(11,*)'Too many sides'
              stop
            endif
            elside(j,i)=ns
            isdel(1,ns)=i
            isidenode(1,ns)=nd1
            isidenode(2,ns)=nd2
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

            isdel(2,ns)=ic3(j,i) !bnd element => bnd side
!           Corresponding side in element ic3(j,i)
            if(ic3(j,i)/=0) then !old internal side
              iel=ic3(j,i)
              index=0
              do k=1,i34(iel)
                if(ic3(iel,k)==i) then
                  index=k
                  exit
                endif
              enddo !k
              if(index==0) then
                write(11,*)'Wrong ball info',i,j
                stop
              endif
              elside(index,iel)=ns
            endif !ic3(j,i).ne.0
          endif !ic3(j,i)==0.or.i<ic3(j,i)
        enddo !j=1,i34
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

!     Open nudging files
      if(ndays/=0) then
        if(ivartype==1) then
          open(37,file='temp_nu.in',form='unformatted')
          open(35,file='salt_nu.in',form='unformatted')
        else
          open(35,file='uv_bcc.th')
        endif
      endif

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
!      sfile = '/home/workspace/ccalmr6/nrldata/'//trim(yr)//'/s3d/s3d.glb8_2f_' // trim(yr) // trim(md) // '00.nc'
!      tfile = '/home/workspace/ccalmr6/nrldata/'//trim(yr)//'/t3d/t3d.glb8_2f_' // trim(yr) // trim(md) // '00.nc'

      if(ivartype==1) then !T,S
        sfile='/home/workspace/ccalmr/nrldata_links/'//trim(yr)//'/s3d/s3d.glb8_2f_'//trim(yr)//trim(md)//'00.nc'
        tfile='/home/workspace/ccalmr/nrldata_links/'//trim(yr)//'/t3d/t3d.glb8_2f_'//trim(yr)//trim(md)//'00.nc'
      else !U,V
        sfile='/home/workspace/ccalmr/nrldata_links/'//trim(yr)//'/u3d/u3d.glb8_2f_'//trim(yr)//trim(md)//'00.nc'
        tfile='/home/workspace/ccalmr/nrldata_links/'//trim(yr)//'/v3d/v3d.glb8_2f_'//trim(yr)//trim(md)//'00.nc'
      endif

      print*, 'Trying to open nc files:',trim(sfile),trim(tfile)      

!     Open files for salinity and temperature (or u,v).
      status = nf90_open(trim(sfile), nf90_nowrite, sid)
      call check(status)
      status = nf90_open(trim(tfile), nf90_nowrite, tid)
      call check(status)

!     Get index information for sampled grid, assumed the same for temperature.
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

!     get the variable ids from salt and temperature files
      if(ivartype==1) then !S&T
        status = nf90_inq_varid(sid, "Salinity", svid)
        call check(status)
        status = nf90_inq_varid(tid, "Temperature", tvid)
        call check(status)
      else !U,V
        status = nf90_inq_varid(sid, "U_Velocity", svid)
        call check(status)
        status = nf90_inq_varid(tid, "V_Velocity", tvid)
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
!      if(ifile==0) then
      allocate( xind(ixlen),stat=ier)
      allocate( yind(iylen),stat=ier)
      allocate( lind(ilen),stat=ier)
      allocate( salt(ixlen, iylen, ilen),stat=ier)
      if(ier /= 0) then
        write(*,*) ' Could not allocate salt'
        stop 'Allocate salt'
      endif
      allocate( temp(ixlen, iylen, ilen),stat=ier)
      if(ier /= 0) then
        write(*,*) ' Could not allocate temp'
        stop 'Allocate temp'
      endif
!      endif
  
!     get the index values in all directions
      status = nf90_get_var(sid, xvid, xind)
      call check(status)
      status = nf90_get_var(sid, yvid, yind)
      call check(status)
      status = nf90_get_var(sid, lvid, lind)
      call check(status)
!     Read salinity/U
      status = nf90_get_var(sid,svid,salt,start=(/1,1,1/),count=(/ixlen, iylen, ilen/),stride=(/1,1,1/))
      call check(status)
!     Read temperature/V
      status = nf90_get_var(tid,tvid,temp,start=(/1,1,1/),count=(/ixlen, iylen, ilen/),stride=(/1,1,1/))
      call check(status)

!     TODO  call read3dvar(sfile, 'Salinity', 1, ixlen, 1, iylen, 1, ilen, salt)
!     TODO  call read3dvar(tfile, 'Temperature', 1, ixlen, 1, iylen, 1, ilen, salt)

      status = nf90_close(sid)
      call check(status)
      status = nf90_close(tid)
      call check(status)

!     Extracting a subset of all lat,lon,h,vertical pos data
!      if(ifile==0) then
      allocate( lat(ixlen,iylen))
      allocate( lon(ixlen,iylen))
      allocate( zm(ixlen, iylen, ilen))
      allocate( h(ixlen,iylen))
      allocate( kbp(ixlen,iylen))
      allocate( ihope(ixlen,iylen))
!      endif

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
      lon=lon-360 !convert to our long.
!      print*, ixlen,iylen,ilen
!     Compute bottom indices and extend S,T data below bottom
!     Level 1 is surface, ilen is towards bottom
      hncom_max=-1 !max. depth in NCOM
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
            if(lon(i,j)>=xl_min.and.lon(i,j)<=xl_max.and.lat(i,j)>=yl_min.and. &
     &lat(i,j)<=yl_max.and.-zm(i,j,kbp(i,j))>hncom_max) hncom_max=-zm(i,j,kbp(i,j))
            do k=1,ilen
              if(k<=kbp(i,j)) then
                if(salt(i,j,k)<-99.or.temp(i,j,k)<-99) then
                  write(11,*)'Fatal: no valid S,T:',i,j,k,salt(i,j,k),temp(i,j,k)
                  stop
                endif
              else !extend
                if(ivartype==1) then
                  salt(i,j,k)=salt(i,j,kbp(i,j))
                  temp(i,j,k)=temp(i,j,kbp(i,j))
                else
                  salt(i,j,k)=0
                  temp(i,j,k)=0
                endif
              endif
            enddo !k
          endif

!          write(12,*)(i-1)*iylen + j, lon(i,j), lat(i,j), kbp(i,j)
        enddo !j
      enddo !i  

!     Calculate a profile at a node closest to (xl_mean,yl_mean) for hotstart.in
      if(ndays==0.and.ivartype==1) then
!       Find a node in NCOM grid for profile
        distm=1.e25 !mion .distance
        i_mean=0; j_mean=0 !init. indices
        do i=ixlen1,ixlen2
          do j=iylen1,iylen2
            dist=(lon(i,j)-xl_mean)**2+(lat(i,j)-yl_mean)**2
            if(dist<distm) then
              distm=dist
              i_mean=i; j_mean=j
            endif
          enddo !j
        enddo !i
        if(i_mean==0) stop 'Failed to find a node in NCOM'
        if(kbp(i_mean,j_mean)==-99) stop 'Found a dry node in NCOM'
        print*, 'Min. distance found = ',distm,'; node indices=',i_mean,j_mean

        nz_r=kbp(i_mean,j_mean)
        allocate(z_r(nz_r),temp_bar(nz_r),salt_bar(nz_r),indx_var(nz_r), &
     &cspline_ypp(nz_r,2),stat=ier)
        if(ier /= 0) stop 'Failed to allocate (6)'
        temp_bar(1:nz_r)=temp(i_mean,j_mean,nz_r:1:-1)
        salt_bar(1:nz_r)=salt(i_mean,j_mean,nz_r:1:-1)
        z_r(1:nz_r)=zm(i_mean,j_mean,nz_r:1:-1)

!       Sort T,S to avoid inversion
        call bubble_sort(1,nz_r,temp_bar,indx_var) 
        call bubble_sort(-1,nz_r,salt_bar,indx_var) 

        open(21,file='ts.ic.0',status='replace')
        do k=1,nz_r
          write(21,*)k,z_r(k),temp_bar(k),salt_bar(k),indx_var(k)
        enddo !k
        close(21)
      endif !ndays==0.and.ivartype==1

!     Write a test file of build points
!      write(15, *)'Test NCOM output.'
!      write(15, *)(ixlen2-ixlen1+1)*(iylen2-iylen1+1) !ixlen * iylen
!      do j=iylen1,iylen2 !1,iylen
!        do i=ixlen1,ixlen2 !1,ixlen
!!         top layer is 1, ilen is bottom layer.
!          if(lat(i,j)>43.73.and.lat(i,j)<44.40.and.lon(i,j)>-126.13.and.lon(i,j)<-124.52) then
!            write(15,*)'Point:',i,j
!            do k=kbp(i,j),1,-1
!              write(15,'(i2,5(1x,f12.4))')k,lon(i,j),lat(i,j),zm(i,j,k),salt(i,j,k),temp(i,j,k)
!            enddo !k
!          endif
!!          write(15, *) (i-1)*iylen+j,lon(i,j),lat(i,j),-h(i,j) 
!!          write(15,*) (i-1)*iylen+j,lon(i,j),lat(i,j),salt(i,j,max(kbp(i,j),1)) !-h(i,j) 
!        enddo
!      enddo
!      stop

!     Compute S,T @ invalid pts based on nearest neighbor
      ihope=1
      do i=ixlen1,ixlen2
        do j=iylen1,iylen2
          if(kbp(i,j)==-99) then !invalid pts  
!         Search along x-direction first
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

!           Search along y-direction 
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
!              ihope(i,j)=1
              if(i1<ixlen1.or.i1>ixlen2.or.j1<iylen1.or.j1>iylen2) then
                write(11,*)'Indices out of bound (1):',i,j,i1,j1
                stop
              endif
              salt(i,j,1:ilen)=salt(i1,j1,1:ilen)
              temp(i,j,1:ilen)=temp(i1,j1,1:ilen)
            endif
          endif !kbp(i,j)==-99

!          write(13, *)(i-1)*iylen+j,lon(i,j),lat(i,j),salt(i,j,ilen) !-h(i,j) 
!          write(15, *) (i-1)*iylen + j, lon(i,j), lat(i,j), ihope(i,j) !-h(i,j) 
        enddo !j=iylen1,iylen2
      enddo !i=ixlen1,ixlen2

!     Find parent elements and levels for hgrid.ll, and do interpolation
      tempout=-99; saltout=-99
      if(ivartype==1) then
        limit=np
      else
        limit=nond0 
      endif

      loop4: do ii=1,limit
        if(ivartype==1) then
          i=ii
        else
          i=iond2(ii)
        endif

        if(iest(i)==1.or.iest(i)==-2) cycle loop4

!       Non-estuary nodes
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
!            if((ix-1)*iylen+iy==27837) then
!              print*, rat,a1+a2+a3+a4,b1+b2
!              print*, x1,x2,x3,x4
!              print*, y1,y2,y3,y4
!              print*, xl(i),yl(i)
!            endif
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

              if(ihope(ix,iy)==0.or.ihope(ix+1,iy)==0.or.ihope(ix+1,iy+1)==0.or.ihope(ix,iy+1)==0) then
                write(11,*)'Hopeless parent:',i,ix,iy
                stop
              else !do interpolation
!               Find vertical level
                do k=kbp2(i),nvrt
                  if(kbp(ix,iy)==-99) then
                    lev=ilen-1; vrat=1
                  else if(z(i,k)<zm(ix,iy,kbp(ix,iy))) then
                    lev=kbp(ix,iy)-1; vrat=1
                  else if(z(i,k)>zm(ix,iy,1)) then !above f.s.
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
           
!                  write(18,*)i,k,ix,iy,lev,vrat,kbp(ix,iy)
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

!                 Enforce lower bound for salt (this is the only occurence in the code)
                  if(ivartype==1.and.z(i,k)<=-ht) saltout(i,k)=max(saltout(i,k),smin)

                enddo !k=kbp2(i),nvrt
 
!               Extend below bottom
                do k=1,kbp2(i)-1
                  tempout(i,k)=tempout(i,kbp2(i))
                  saltout(i,k)=saltout(i,kbp2(i))
                enddo !k
              endif !ihope(ix,iy)
              cycle loop4
            endif !rat<small1
          enddo !iy=iylen1,iylen2-1
        enddo !ix=ixlen1,ixlen2-1
        if(ixy(i,1)==0.or.ixy(i,2)==0) then
          write(11,*)'Cannot find a parent element:',i
          stop
        endif
      end do loop4 !i=1,limit

!     Estuary nodes (for S,T only)
      if(ivartype==1) then
        if(tempout(ianchor1,1)==-99) then
          write(11,*)'Anchor node not assigned'
          stop
        endif
        do i=1,np
          if(iest(i)==1.or.iest(i)==-2) then
            if(nanchor==1) then
              tempout(i,1:nvrt)=tempout(ianchor1,1:nvrt)
              saltout(i,1:nvrt)=saltout(ianchor1,1:nvrt)
            else !=2
              if(xl(ianchor2)-xl(ianchor1)==0) then
                write(11,*)'Wrong anchor pts:',ianchor1,ianchor2
                stop
              endif
              xrat=(xl(i)-xl(ianchor1))/(xl(ianchor2)-xl(ianchor1))
              xrat=max(0.,min(1.,xrat))
              do k=1,nvrt
!               Since the nudging is not done in estuary, only the river temp. at the beginning is used
                tempout(i,k)=tempout(ianchor1,k)*(1-xrat)+tempth*xrat
                saltout(i,k)=saltout(ianchor1,k)*(1-xrat) !+0
              enddo !k 
            endif
          endif !inside estuary
        enddo !i
      endif !ivartype==1
    
!     Output hotstart or nudging files
      if(ndays==0) then !ivartype must be 1
!       hotstart.in
!       Note: although I swapped order of indices in SELFE for S,T,
!             I didn't do the same for them in this code; as long as 
!             the order of writing is correct in hotstart.in, it does not matter
        do i=1,ns
          n1=isidenode(1,i)
          n2=isidenode(2,i)
          do k=1,nvrt
            tsd(i,k)=(tempout(n1,k)+tempout(n2,k))/2
            ssd(i,k)=(saltout(n1,k)+saltout(n2,k))/2
          enddo !k
!          write(88,*)i,xcj(i),ycj(i),ssd(i,1),ssd(i,nvrt)
        enddo !i

        do i=1,ne
          n1=elnode(1,i)
          n2=elnode(2,i)
          n3=elnode(3,i)
          do k=2,nvrt
            tsel(k,i,1)=(tempout(n1,k)+tempout(n2,k)+tempout(n3,k)+tempout(n1,k-1)+tempout(n2,k-1)+tempout(n3,k-1))/6
            tsel(k,i,2)=(saltout(n1,k)+saltout(n2,k)+saltout(n3,k)+saltout(n1,k-1)+saltout(n2,k-1)+saltout(n3,k-1))/6
          enddo !k
          tsel(1,i,1)=tsel(2,i,1) !mainly for hotstart format
          tsel(1,i,2)=tsel(2,i,2)
        enddo !i

        if(imodel.eq."S") then
!         MPI SELFE
!         Overwrite T,S with a profile below a depth for hotstart.in only
!         Cubic spline coefficients (save for interpolation later)
          call cubic_spline(nz_r,z_r,temp_bar,0.,0.,cspline_ypp(:,1))
          call cubic_spline(nz_r,z_r,salt_bar,0.,0.,cspline_ypp(:,2))

!         T,S @ nodes
          do i=1,np
            if(dp(i)>=h_mean) then

!             Wet nodes with depth > h_mean
              if(z(i,kbp2(i))<z_r(1)) then !.or.z(nvrt,i)>z_r(nz_r)) then
                stop 'MAIN: node depth too big for ts.ic'
              endif 
              call eval_cubic_spline(nz_r,z_r,temp_bar,cspline_ypp(:,1),nvrt-kbp2(i)+1,z(i,kbp2(i):nvrt), &
     &0,z_r(1),z_r(nz_r),swild(kbp2(i):nvrt,1))
              call eval_cubic_spline(nz_r,z_r,salt_bar,cspline_ypp(:,2),nvrt-kbp2(i)+1,z(i,kbp2(i):nvrt), &
     &0,z_r(1),z_r(nz_r),swild(kbp2(i):nvrt,2))
           
              do k=kbp2(i),nvrt
                if(z(i,k)<-h_mean) then
                  tempout(i,k)=swild(k,1)
                  saltout(i,k)=swild(k,2)
                endif !z(i,k)<-h_mean
              enddo !k
!           Impose no slip b.c. to be consistent with ELM transport
!          if(Cdp(i)/=0) then
!            tem0(kbp(i),i)=tem0(kbp(i)+1,i)
!            sal0(kbp(i),i)=sal0(kbp(i)+1,i)
!          endif

!             Extend
              do k=1,kbp2(i)-1
                tempout(i,k)=tempout(i,kbp2(i))
                saltout(i,k)=saltout(i,kbp2(i))
              enddo !k
            endif !dp(i)>=h_mean
          enddo !i=1,npa

!         T,S @ sides 
          do i=1,ns
            n1=isidenode(1,i); n2=isidenode(2,i)
            dps=(dp(n1)+dp(n2))/2
            if(dps>=h_mean) then
              kbs=max(kbp2(n1),kbp2(n2))
              zs(kbs:nvrt)=(z(n1,kbs:nvrt)+z(n2,kbs:nvrt))/2
              if(zs(kbs)<z_r(1)) then !.or.zs(nvrt,i)>z_r(nz_r)) then
                stop 'MAIN: side depth too big for ts.ic'
              endif 
              call eval_cubic_spline(nz_r,z_r,temp_bar,cspline_ypp(:,1),nvrt-kbs+1,zs(kbs:nvrt), &
     &0,z_r(1),z_r(nz_r),swild(kbs:nvrt,1))
              call eval_cubic_spline(nz_r,z_r,salt_bar,cspline_ypp(:,2),nvrt-kbs+1,zs(kbs:nvrt), &
     &0,z_r(1),z_r(nz_r),swild(kbs:nvrt,2))
             
              do k=kbs,nvrt
                if(zs(k)<-h_mean) then
                  tsd(i,k)=swild(k,1)
                  ssd(i,k)=swild(k,2)
                endif                
              enddo !k

!             Extend
              do k=1,kbs-1
                tsd(i,k)=tsd(i,kbs)
                ssd(i,k)=ssd(i,kbs)
              enddo !k
            endif !dps>=h_mean
          enddo !i=1,nsa

!         T,S @ elements
          do i=1,ne
            n1=elnode(1,i)
            n2=elnode(2,i)
            n3=elnode(3,i)
            dpe=(dp(n1)+dp(n2)+dp(n3))/3
            if(dpe>=h_mean) then
              kbe=max(kbp2(n1),kbp2(n2),kbp2(n3))
              ze(kbe:nvrt)=(z(n1,kbe:nvrt)+z(n2,kbe:nvrt)+z(n3,kbe:nvrt))/3
              if(ze(kbe)<z_r(1)) then !.or.ze(nvrt,i)>z_r(nz_r)) then
                stop 'MAIN: ele. depth too big for ts.ic'
              endif 

              do k=kbe+1,nvrt
                swild2(k)=(ze(k)+ze(k-1))/2
              enddo !k
              call eval_cubic_spline(nz_r,z_r,temp_bar,cspline_ypp(:,1),nvrt-kbe,swild2(kbe+1:nvrt), &
     &0,z_r(1),z_r(nz_r),swild(kbe+1:nvrt,1))
              call eval_cubic_spline(nz_r,z_r,salt_bar,cspline_ypp(:,2),nvrt-kbe,swild2(kbe+1:nvrt), &
     &0,z_r(1),z_r(nz_r),swild(kbe+1:nvrt,2))

              do k=kbe+1,nvrt
                if(swild2(k)<-h_mean) then
                  tsel(k,i,1:2)=swild(k,1:2)
                endif
              enddo !k

!             Extend
              do k=1,kbe
                tsel(k,i,1:2)=tsel(kbe+1,i,1:2)
              enddo !k
            endif !dpe>=h_mean
          enddo !i=1,nea

       
          open(36,file='hotstart.in',form='unformatted',status='replace')
          write(36) 0.d0,0,1
          do i=1,ne
            write(36) i,0,(0.d0,dble(tsel(j,i,1:2)),(0.d0,0.d0,l=1,ntracers),j=1,nvrt)
          enddo !i
          do i=1,ns
            write(36) i,0,(0.d0,0.d0,dble(tsd(i,j)),dble(ssd(i,j)),j=1,nvrt)
          enddo !i
          do i=1,np
            write(36) i,0.d0,0,(dble(tempout(i,j)),dble(saltout(i,j)), &
                      dble(tempout(i,j)),dble(saltout(i,j)),0.d0,0.d0, &
                      0.d0,0.d0,0.d0,0.d0,0.d0,j=1,nvrt)
          enddo !i
          close(36)

        else !ELCIRC
          ihot_len=nbyte*(3+4*ne+2*ne*(nvrt+1)+4*ne*nvrt+4*ns*(nvrt+1)+4*ns*nvrt+ &
     &         3*np+7*np*(nvrt+1)+8*np*nvrt+1)+12
          if(itur.eq.3) ihot_len=ihot_len+nbyte*4*ns*(nvrt+1)
          open(36,file='hotstart.in',access='direct',recl=ihot_len)
          if(itur.eq.3) then
            write(36,rec=1)0.d0,0,(0.d0,0.d0,(0.d0,j=0,nvrt),(dble(tsel(j,i,1)),dble(tsel(j,i,2)),j=1,nvrt),i=1,ne), &
     &((0.d0,0.d0,j=0,nvrt),(dble(tsd(i,j)),dble(ssd(i,j)),j=1,nvrt),i=1,ns) &
     &,(0.d0,0,(0.d0,0.d0,0.d0,0,j=0,nvrt), &
     &(dble(tempout(i,j)),dble(saltout(i,j)),dble(tempout(i,j)),dble(saltout(i,j)),j=1,nvrt),i=1,np), &
     &((0.d0,0.d0,j=0,nvrt),i=1,ns),1,'1           '
          else
            write(36,rec=1)0.d0,0,(0.d0,0.d0,(0.d0,j=0,nvrt),(dble(tsel(j,i,1)),dble(tsel(j,i,2)),j=1,nvrt),i=1,ne), &
     &((0.d0,0.d0,j=0,nvrt),(dble(tsd(i,j)),dble(ssd(i,j)),j=1,nvrt),i=1,ns) &
     &,(0.d0,0,(0.d0,0.d0,0.d0,0,j=0,nvrt), &
     &(dble(tempout(i,j)),dble(saltout(i,j)),dble(tempout(i,j)),dble(saltout(i,j)),j=1,nvrt),i=1,np), &
     &1,'1           '
          endif
        endif

      else !ndays/=0
        if(ivartype==1) then !ST
!         nudging files
          write(37)ifile*86400.
          write(35)ifile*86400.
          do i=1,np
            write(37)(tempout(i,j),j=1,nvrt)
            write(35)(saltout(i,j),j=1,nvrt)
          enddo !i
        else !UV
          write(35,*)ifile*86400.,nond0,nvrt,(saltout(iond2(i),1:nvrt), &
     &tempout(iond2(i),1:nvrt),i=1,nond0)
        endif

!       Debug
!        do i=1,nond0
!          nd=iond2(i)
!          if(i==1.or.i==2.or.i==15) write(40+i,*)ifile,tempout(nd,nvrt)
!        enddo !i

      endif !output
 
!     Debug
!      if(ifile==0) then
!!          write(20,*)i,xl(i),yl(i),saltout(i,1)
!        do i=1,ne
!          write(20,*)i,(tsel(k,i,1),k=1,nvrt)
!        enddo !i
!      endif

      deallocate(lat,lon,zm,h,kbp,ihope,xind,yind,lind,salt,temp)
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
        print *, 'failed to open nc files'
        stop
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

!     Inefficient bubble sort routine for sorting into ascending or descending order
!     If there are equal entries in the list the sorted indices (indx_var)
!     is the same as tje original indices for those entries.
!     Input
!     sort_type: 1: ascending order; -1: descending order;
!     ndim:  dimension parameter;
!     In/Out
!     var(ndim): list to be sorted
!     Output
!     indx_var(ndim):  list of original indices after sorting.      
      subroutine bubble_sort(sort_type,ndim,var,indx_var)
      implicit none
      integer, parameter :: rkind=4
      integer, intent(in) :: sort_type,ndim
      real(rkind), intent(inout) :: var(ndim)
      integer, intent(out) :: indx_var(ndim)

      integer :: i,j,itmp
      real(rkind) :: tmp
      
      if(iabs(sort_type)/=1) stop 'bubble_sort: Wrong sort_type'
      indx_var(:)=(/(i,i=1,ndim)/)
      do i=1,ndim-1
        do j=i+1,ndim
          if(sort_type*var(i)>sort_type*var(j)) then !swap
            tmp=var(i)
            var(i)=var(j)
            var(j)=tmp
            itmp=indx_var(i)
            indx_var(i)=indx_var(j)
            indx_var(j)=itmp
          endif
        enddo !j
      enddo !i

      end subroutine bubble_sort

!===============================================================================
!     Do interpolation with cubic spline
!     Needs coefficients from routine cubic_spline()
!     Inputs: 
!            npts: dimesion for xcor etc; # of data pts;
!            xcor(npts),yy(npts): x and y coordinates of the original function 
!                                 (same as in cubic_spline()); xcor in ascending order;
!            ypp(npts): 2nd deriavtives (output from cubic_spline);
!            npts2: # of output pts;
!            xout(npts2): x coordinates of the output pts (no ordering required);
!            xmax: if xout>xmax, it is reset to xmax;
!            ixmin (0 or 1): bottom option;
!            xmin: if xout<xcor(1), it is either reset to xmin (ixmin=0), or 
!                  to xcor(1) (ixmin=1), i.e. yyout takes the value of yy(1), and
!                  xmin is not used except for debugging messages.
!     Output: 
!            yyout(npts2): output y values; if xmin>xmax, yyout=yy(1).
!===============================================================================
      subroutine eval_cubic_spline(npts,xcor,yy,ypp,npts2,xout,ixmin,xmin,xmax,yyout)
      implicit real(4)(a-h,o-z), integer(i-n)
      integer, intent(in) :: npts,npts2,ixmin
      real(4), intent(in) :: xcor(npts),yy(npts),ypp(npts),xout(npts2),xmin,xmax
      real(4), intent(out) :: yyout(npts2)

      if(xmin>xmax) then
        yyout=yy(1); return
      endif

      do i=1,npts2
        ifl=0 !flag
        xtmp=min(xout(i),xmax)
        if(ixmin==0) then
          xtmp=max(xtmp,xmin)
        else
          if(xout(i)<xcor(1)) then
            yyout(i)=yy(1); cycle
          endif
        endif

        do j=1,npts-1
          if(xtmp>=xcor(j).and.xtmp<=xcor(j+1)) then
            ifl=1
            aa=(xcor(j+1)-xtmp)/(xcor(j+1)-xcor(j))
            bb=1-aa
            cc=(aa*aa*aa-aa)*(xcor(j+1)-xcor(j))/6
            dd=(bb*bb*bb-bb)*(xcor(j+1)-xcor(j))/6
            yyout(i)=aa*yy(j)+bb*yy(j+1)+cc*ypp(j)+dd*ypp(j+1)
            exit
          endif
        enddo !j
        if(ifl==0) then
          write(*,*)'EVAL_CUBIC: Falied to find:',i,xtmp,xmin,xmax
          stop
        endif
      enddo !i=1,npts2

      end subroutine eval_cubic_spline

!===============================================================================
!     Generate coefficients (2nd derivatives) for cubic spline for interpolation later
!     Inputs: 
!            npts: dimesion for xcor etc; # of data pts;
!            xcor(npts): x coordinates; must be in ascending order (and distinctive);
!            yy(npts): functional values; 
!            yp1 and yp2: 1st derivatives at xcor(1) and xcor(npts);
!     Output: 
!            ypp(npts): 2nd deriavtives used in interpolation.
!===============================================================================
      subroutine cubic_spline(npts,xcor,yy,yp1,yp2,ypp)
      implicit real(4)(a-h,o-z), integer(i-n)
      integer, intent(in) :: npts
      real(4), intent(in) :: xcor(npts),yy(npts),yp1,yp2
      real(4), intent(out) :: ypp(npts)
  
      real(4) :: alow(npts),bdia(npts),cupp(npts),rrhs(npts,1),gam(npts)

      do k=1,npts
        if(k==1) then
          bdia(k)=(xcor(k+1)-xcor(k))/3
          if(bdia(k)==0) then
            write(*,*)'CUBIC_SP: bottom problem:',xcor(k+1),xcor(k)
            stop
          endif
          cupp(k)=bdia(k)/2
          rrhs(k,1)=(yy(k+1)-yy(k))/(xcor(k+1)-xcor(k))-yp1
        else if(k==npts) then
          bdia(k)=(xcor(k)-xcor(k-1))/3
          if(bdia(k)==0) then
            write(*,*)'CUBIC_SP: surface problem:',xcor(k),xcor(k-1)
            stop
          endif
          alow(k)=bdia(k)/2
          rrhs(k,1)=-(yy(k)-yy(k-1))/(xcor(k)-xcor(k-1))+yp2
        else
          bdia(k)=(xcor(k+1)-xcor(k-1))/3
          alow(k)=(xcor(k)-xcor(k-1))/6
          cupp(k)=(xcor(k+1)-xcor(k))/6
          if(alow(k)==0.or.cupp(k)==0) then
            write(*,*)'CUBIC_SP: middle problem:',xcor(k),xcor(k-1),xcor(k+1)
            stop
          endif
          rrhs(k,1)=(yy(k+1)-yy(k))/(xcor(k+1)-xcor(k))-(yy(k)-yy(k-1))/(xcor(k)-xcor(k-1))
        endif
      enddo !k

      call tridag(npts,1,npts,1,alow,bdia,cupp,rrhs,ypp,gam)
!      ypp(:)=soln(:,1)

      end subroutine cubic_spline

!===============================================================================
subroutine tridag(nmax,nvec,n,nc,a,b,c,r,u,gam)
!-------------------------------------------------------------------------------
! This program solves a tridiagonal system. It was adapted from "Numerical 
! Recipes in FORTRAN (pp.43 ).
!
! a,b,c,r: input vectors and are not modified by this program.
! b is the main diagonal, a below and c above. a(1) and c(n) are not used.
! r is the r.h.s.
! (nmax,nvec) is the dimension of r() _and_ u() in the driving routine.
! n: actual rank of the system.
! nc: input; actual # of columns of rhs (<= nvec).
! u: output with nc columns (depending on input nc).
! gam: a working array.
!-------------------------------------------------------------------------------
  implicit real(4)(a-h,o-z), integer(i-n)

  integer, intent(in) :: nmax,nvec,n,nc
  real(4), dimension(nmax), intent(in) :: a,b,c
  real(4), dimension(nmax,nvec), intent(in) :: r
  real(4), dimension(nmax), intent(out) :: gam
  real(4), dimension(nmax,nvec), intent(out) :: u

  if(n<1) stop 'TRIDAG: n must be >= 1'
  if(nc>nvec) stop 'TRIDAG: Increase # of columns'
  if(b(1)==0) stop 'TRIDAG:  b(1)=0'

  bet=b(1)
  do i=1,nc
    u(1,i)=r(1,i)/bet
  enddo

  do j=2,n
    gam(j)=c(j-1)/bet
    bet=b(j)-a(j)*gam(j)
    if(bet.eq.0) stop 'TRIDAG: failed'
    do i=1,nc
      u(j,i)=(r(j,i)-a(j)*u(j-1,i))/bet
    enddo !i
  enddo !j

! Backsubstitution
  do j=n-1,1,-1
    do i=1,nc
      u(j,i)=u(j,i)-gam(j+1)*u(j+1,i)
    enddo
  enddo
  
end subroutine tridag

