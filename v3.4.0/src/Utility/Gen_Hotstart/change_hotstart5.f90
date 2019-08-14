! Generate hotstart.in (MPI SELFE) from gridded data (nc file)
! Assume elev=0 in 3D variables
! History: added ivcor=1 option
!   Input: 
!     (1) hgrid.gr3;
!     (2) hgrid.ll;
!     (3) vgrid.in (SELFE or ELCIRC test02k4);
!     (4) estuary.gr3 (flags for extrapolating S,T, vel.): depth=0: outside; =1: inside
!     (5) obs_T.gr3: 1: obs. region (near Fehm & Darss); 0: otherwise. T is corrected in the obs region
!                    (search for 'T_obs');
!     (6) obs_S.gr3: 1: obs. region (near Fehm); 2: region around Oder; 3: rest of Baltic; 
!                    0: otherwise. S is corrected in the obs region (search for 'S_obs');
!     (7) gen_hot.in: 1st line: name of nc file;
!                    2nd line: 1: include vel and elev. in hotstart.in (not working at the moment); 0: only T,S
!                    3rd line: T,S values for estuary points defined in estuary.gr3;
!                    4th line: T,S values for pts outside bg grid in nc
!                    5th line: month # (1-12)
!     (8) nc file (prescribed in (7));
!   Output: hotstart.in ; fort.11 (fatal errors); fort.2[013], fort.9[89] (debug outputs)

! pgf90 -O2 -mcmodel=medium -Mbounds -Bstatic -o gen_hot7 gen_hot7.f90 -I/usr/local/amd64a/opteron/pgi/netcdf-3.6.2/include -L/usr/local/amd64a/opteron/pgi/netcdf-3.6.2/lib -lnetcdf

      program gen_hot
!     netcdf modules from symlinks
      use typeSizes
      use netcdf

!      implicit none

!      integer, parameter :: debug=1
      integer, parameter :: mnp=350000
      integer, parameter :: mne=700000
      integer, parameter :: mns=1050000
      integer, parameter :: mnv=70
      integer, parameter :: mnope=40 !max # of open bnd segments
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
      integer :: uvid, vvid ! vel variable IDs
      integer, dimension(nf90_max_var_dims) :: dids
             
!     Local variables for data
      real (kind = FourByteReal), dimension(:), allocatable :: xind, yind, lind 
!     Lat, lon, bathymetry
      real (kind = FourByteReal), dimension(:,:), allocatable :: lat, lon, hnc
!     Vertical postion, salinity, and temperature
      real (kind = FourByteReal), dimension(:,:,:), allocatable :: zm
      real (kind = FourByteReal), dimension(:,:,:,:), allocatable :: salt,temp
      integer, dimension(:,:), allocatable :: kbp,ihope
!     File names for netcdf files
      character(len=1024) :: sfile
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

      dimension xl(mnp),yl(mnp),nm(mne,4),i34(mne)
      dimension ztot(0:mnv),sigma(mnv),cs(mnv),z(mnp,mnv),iest(mnp),ixy(mnp,2),arco(3)
      dimension wild(100),wild2(100,2),iobs_T(mnp),iobs_S(mnp),T_obs(mnp,mnv),S_obs1(mnp,mnv)
      dimension tempout(mnp,mnv), saltout(mnp,mnv),month_day(12),S_obs2(mnp,mnv)
      dimension tsd(mns,mnv),ssd(mns,mnv),tsel(mnv,mne,2)
      dimension nne(mnp),ine(mnp,mnei),ic3(mne,4),nx(4,4,3),js(mne,4),is(mns,2),isidenode(mns,2)
      dimension xcj(mns),ycj(mns),nond(mnope),iond(mnope,mnond),iob(mnope),iond2(mnope*mnond)
      dimension kbp2(mnp),sigma_lcl(mnv,mnp)
      real*8 :: hsm(mnv),dp(mnp),h_tran1,h_tran2,z0,z_1,sp

!     First statement
!     Constant
      tempmin=-10 !used to check non-junk values from nc files

      open(10,file='gen_hot.in',status='old')
      read(10,*) sfile !nc file name
      read(10,*) iuv !1: include vel and elev. in hotstart.in; 0: only T,S
      read(10,*) tem_es,sal_es !T,S values for estuary points defined in estuary.gr3
      read(10,*) tem_outside,sal_outside !T,S values for pts outside bg grid in nc
      read(10,*) mon0 !month #
      close(10)

!     Read in hgrid and vgrid
      open(13,file='obs_T.gr3',status='old')
      open(15,file='obs_S.gr3',status='old')
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
      read(13,*)
      read(13,*)
      read(15,*)
      read(15,*)
      do i=1,np
        read(14,*)j,xtmp,ytmp,dp(i)
        read(16,*)j,xl(i),yl(i) !,dp(i)
        read(17,*)j,xtmp,ytmp,tmp
        iest(i)=tmp
        if(iest(i)/=0.and.iest(i)/=1) then
          write(11,*)'Estuary flag wrong:',i,iest(i)
          stop
        endif
        read(13,*)j,xtmp,ytmp,tmp1
        read(15,*)j,xtmp,ytmp,tmp2
        iobs_T(i)=tmp1
        iobs_S(i)=tmp2
        if(iobs_T(i)/=0.and.iobs_T(i)/=1.or.iobs_S(i)<0.or.iobs_S(i)>3) then
          write(11,*)'Obs. flag wrong:',i,iobs_T(i),iobs_S(i)
          stop
        endif
      enddo !i
      close(13); close(15); close(17); close(16);

      do i=1,ne
        read(14,*)j,i34(i),(nm(i,l),l=1,i34(i))
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

!     V-grid
      read(19,*)ivcor

      if(ivcor==2) then !SZ
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

!       Z-coord.
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
          do k=1,kz-1
            z(i,k)=ztot(k)
          enddo !k
        enddo !i=1,np

      else if(ivcor==1) then !localized
        read(19,*)nvrt
        if(nvrt>mnv.or.nvrt<3) then
          write(11,*)'nvrt > mnv or nvrt<4'
          stop
        endif

        do i=1,np
          read(19,*)j,kbp2(i),sigma_lcl(kbp2(i):nvrt,i)

          z(i,kbp2(i))=-dp(i) !to avoid underflow
          z(i,nvrt)=0
          do k=kbp2(i)+1,nvrt-1
            z(i,k)=sigma_lcl(k,i)*dp(i)
          enddo !k

          !Check
          do k=kbp2(i)+1,nvrt
            if(z(i,k)-z(i,k-1)<=0) then
              write(11,*)'Inverted Z:',z(i,k),z(i,k-1),k,i
              stop
            endif
          enddo !k

          !Extend bottom (as these are used)
          z(i,1:kbp2(i)-1)=z(i,kbp2(i))
        enddo !i=1,np

      else
        write(11,*)'Unknown ivcor'
        stop
      endif !ivcor

!     T, S profiles according to obs
!     shallow pts have junk values
      do i=1,np
        do k=1,nvrt
          S_obs2(i,k)=7.06 !Oder
          if(dp(i)<=0.1) then
            T_obs(i,k)=17
            S_obs1(i,k)=10 !Fehm region
          else !z() is valid
            if(z(i,k)<=-10) then
              S_obs1(i,k)=max(28.,min(31.,28-3*(z(i,k)+10)/10))
            else
              S_obs1(i,k)=max(10.,min(28.,10-18*(z(i,k)+6)/4))
            endif

            if(z(i,k)<=-23) then
              T_obs(i,k)=7.74
            else if(z(i,k)<=-20) then
              T_obs(i,k)=7.74+(8.28-7.74)*(z(i,k)+23)/3
            else if(z(i,k)<=-15) then
              T_obs(i,k)=8.28+(11.27-8.28)*(z(i,k)+20)/5
            else if(z(i,k)<=-10) then
              T_obs(i,k)=11.27+(15.43-11.27)*(z(i,k)+15)/5
            else if(z(i,k)<=-6) then
              T_obs(i,k)=15.43+(17-15.43)*(z(i,k)+10)/4
            else
              T_obs(i,k)=17
            endif
          endif
        enddo !k
      enddo !i=1,np

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
        do j=1,i34(i)
          ic3(i,j)=0 !index for bnd sides
          nd1=nm(i,nx(i34(i),j,1))
          nd2=nm(i,nx(i34(i),j,2))
          do k=1,nne(nd1)
            ie=ine(nd1,k)
            if(ie/=i.and.(nm(ie,1)==nd2.or.nm(ie,2)==nd2.or.nm(ie,3)==nd2.or.(i34(ie)==4.and.nm(ie,4)==nd2))) ic3(i,j)=ie
          enddo !k
        enddo !j
      enddo !i

      ns=0 !# of sides
      do i=1,ne
        do j=1,i34(i)
          nd1=nm(i,nx(i34(i),j,1))
          nd2=nm(i,nx(i34(i),j,2))
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
              js(iel,index)=ns
            endif !ic3(i,j).ne.0
          endif !ic3(i,j)==0.or.i<ic3(i,j)
        enddo !j=1,i34
      enddo !i=1,ne

      if(ns<ne.or.ns<np) then
        write(11,*)'Weird grid with ns < ne or ns < np',np,ne,ns
        stop
      endif
     print*, 'Done computing geometry...'

!     The following section needs to be updated for different nc format
!     Open nc file 
      status = nf90_open(trim(sfile), nf90_nowrite, sid)
      call check(status)

!      status = nf90_inq_dimid(sid, "lonc", xdid)
!      call check(status)
!      status = nf90_inq_dimid(sid, "latc", ydid)
!      call check(status)
!      status = nf90_inq_dimid(sid, "sigma", ldid)
!      call check(status)
!      print*, 'Done reading dimensions...'

      status = nf90_inq_varid(sid, "lon", xvid)
      call check(status)
      status = nf90_inq_varid(sid, "lat", yvid)
      call check(status)
      status = nf90_inq_varid(sid, "zax", lvid)
      call check(status)
      print*, 'Done reading grid IDs...'

      status = nf90_inq_varid(sid, "salt", svid)
      call check(status)
      status = nf90_inq_varid(sid, "temp", tvid)
      call check(status)
      print*, 'Done reading variable IDs'

!     Get the lengths of the dimensions from the salt file.
!     Assumed same as for T, u,v
!     WARNING: indices reversed from ncdump!
      status = nf90_Inquire_Variable(sid, svid, dimids = dids)
      call check(status)
      status = nf90_Inquire_Dimension(sid, dids(1), len = ixlen)
      call check(status)
      status = nf90_Inquire_Dimension(sid, dids(2), len = iylen)
      call check(status)
      status = nf90_Inquire_Dimension(sid, dids(3), len = ilen)
      call check(status)
      status = nf90_Inquire_Dimension(sid, dids(4), len = ntime)
      call check(status)

      print*, 'ixlen,iylen,ilen,ntime= ',ixlen,iylen,ilen,ntime

!     allocate memory.
      allocate(xind(ixlen),stat=ier)
      allocate(yind(iylen),stat=ier)
      allocate(lind(ilen),stat=ier)
      allocate(lat(ixlen,iylen))
      allocate(lon(ixlen,iylen))
      allocate(zm(ixlen, iylen, ilen))
!      allocate(hnc(ixlen,iylen))
      allocate(kbp(ixlen,iylen))
      allocate(ihope(ixlen,iylen))
!      allocate(uvel(ixlen,iylen,ilen),stat=ier)
!      allocate(vvel(ixlen,iylen,ilen),stat=ier)
      allocate(salt(ixlen,iylen,ilen,ntime),stat=ier)
      allocate(temp(ixlen,iylen,ilen,ntime),stat=ier)
  
!     get the index values in all directions
      status = nf90_get_var(sid, xvid, xind)
      call check(status)
      status = nf90_get_var(sid, yvid, yind)
      call check(status)
      !lind may be sigma coord.
      status = nf90_get_var(sid, lvid, lind)
      call check(status)

      do i=1,ixlen
        lon(i,:)=xind(i)
      enddo !i
      do j=1,iylen
        lat(:,j)=yind(j)
      enddo !j
!      lon=lon-360 !convert to our long.

!     Compute z-coord. (assuming eta=0)
!     In zm(), 1 is bottom; ilen is surface (i.e. follow SELFE convention)
      do i=1,ixlen
        do j=1,iylen
          do k=1,ilen
            zm(i,j,k)=-lind(ilen+1-k)
          enddo !k
        enddo !j
      enddo !i

!     Read T,S,u,v
      status = nf90_get_var(sid,svid,salt(:,:,ilen:1:-1,:)) 
      call check(status)
      status = nf90_get_var(sid,tvid,temp(:,:,ilen:1:-1,:)) 
      call check(status)

      status = nf90_close(sid)
      call check(status)
!     End reading nc file

!     Test
!      icount=0
!      do i=1,ixlen
!        do j=1,iylen
!          icount=icount+1
!          write(98,*)icount,lon(i,j),lat(i,j),salt(i,j,1,mon0)
!          write(99,*)icount,lon(i,j),lat(i,j),temp(i,j,1,mon0)
!        enddo !j
!      enddo !i

!     At this point all variables have been read, you may proceed with processing.
!
!     Compute bottom indices for nc files
      kbp=1 !no Z layers
      do i=1,ixlen
        do j=1,iylen
          if(salt(i,j,ilen,mon0)<0) kbp(i,j)=-1 !dry

          if(kbp(i,j)==1) then !wet
          !Extend near bottom
            klev0=-1 !flag
            do k=1,ilen
              if(salt(i,j,k,mon0)>=0) then
                klev0=k; exit
              endif
            enddo !k
            if(klev0<=0) then
              write(11,*)'Impossible (1):',i,j,salt(i,j,ilen,mon0)
              stop
            endif !klev0
            salt(i,j,1:klev0-1,mon0)=salt(i,j,klev0,mon0)
            temp(i,j,1:klev0-1,mon0)=temp(i,j,klev0,mon0)

            !Check
            do k=1,ilen
              if(salt(i,j,k,mon0)<0.or.temp(i,j,k,mon0)<tempmin) then
               write(11,*)'Fatal: no valid S,T:',i,j,k,salt(i,j,k,mon0), &
     &temp(i,j,k,mon0)
                stop
              endif
            enddo !k
          endif !kbp
        enddo !j
      enddo !i  

!     Extend S,T (from nc) @ invalid pts based on nearest neighbor
!     Search around neighborhood of a pt
      do i=1,ixlen
        do j=1,iylen
          if(kbp(i,j)==-1) then !invalid pts  
            !Compute max possible tier #
            mmax=max(i-1,ixlen-i,j-1,iylen-j)

            m=0 !tier #
            loop6: do
              m=m+1
              do ii=max(-m,1-i),min(m,ixlen-i)
                i3=max(1,min(ixlen,i+ii))
                do jj=max(-m,1-j),min(m,iylen-j)
                  j3=max(1,min(iylen,j+jj))
                  if(kbp(i3,j3)==1) then !found
                    i1=i3; j1=j3
                    exit loop6   
                  endif
                enddo !jj
              enddo !ii

              if(m==mmax) then
                write(11,*)'Max. exhausted:',i,j,mmax
                write(11,*)'kbp'
                do ii=1,ixlen
                  do jj=1,iylen
                    write(11,*)ii,jj,kbp(ii,jj)
                  enddo !jj
                enddo !ii
                stop
              endif
            end do loop6

            salt(i,j,1:ilen,mon0)=salt(i1,j1,1:ilen,mon0)
            temp(i,j,1:ilen,mon0)=temp(i1,j1,1:ilen,mon0)
          endif !kbp(i,j)==-1
        enddo !j=iylen1,iylen2
      enddo !i=ixlen1,ixlen2

!     Test
      icount=0
      do i=1,ixlen
        do j=1,iylen
          icount=icount+1
          write(98,*)icount,lon(i,j),lat(i,j),salt(i,j,1,mon0)
          write(99,*)icount,lon(i,j),lat(i,j),temp(i,j,1,mon0)
        enddo !j
      enddo !i

!     Find parent elements and levels for hgrid.ll, and do interpolation
!     Error: have not added elev., u,v
      tempout=-99; saltout=-99
      loop4: do i=1,np
        if(iest(i)==1) then
          tempout(i,:)=tem_es
          saltout(i,:)=sal_es
          cycle loop4
        endif

!       Non-estuary nodes
        ixy(i,1)=0; ixy(i,2)=0
        do ix=1,ixlen-1 
          do iy=1,iylen-1 
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

!             Find vertical level
              do k=1,nvrt
                if(kbp(ix,iy)==-1) then
                  lev=ilen-1; vrat=1
                else if(z(i,k)<=zm(ix,iy,kbp(ix,iy))) then
                  lev=kbp(ix,iy); vrat=0
                else if(z(i,k)>=zm(ix,iy,ilen)) then !above f.s.
                  lev=ilen-1; vrat=1
                else
                  lev=-99 !flag
                  do kk=1,ilen-1
                    if(z(i,k)>=zm(ix,iy,kk).and.z(i,k)<=zm(ix,iy,kk+1)) then
                      lev=kk
                      vrat=(zm(ix,iy,kk)-z(i,k))/(zm(ix,iy,kk)-zm(ix,iy,kk+1))
                      exit
                    endif
                  enddo !kk
                  if(lev==-99) then
                    write(11,*)'Cannot find a level:',i,k,z(i,k),zm(ix,iy,:)
                    stop
                  endif
                endif
          
                wild2(1,1)=temp(ix,iy,lev,mon0)*(1-vrat)+temp(ix,iy,lev+1,mon0)*vrat
                wild2(1,2)=salt(ix,iy,lev,mon0)*(1-vrat)+salt(ix,iy,lev+1,mon0)*vrat
                wild2(2,1)=temp(ix+1,iy,lev,mon0)*(1-vrat)+temp(ix+1,iy,lev+1,mon0)*vrat
                wild2(2,2)=salt(ix+1,iy,lev,mon0)*(1-vrat)+salt(ix+1,iy,lev+1,mon0)*vrat
                wild2(3,1)=temp(ix+1,iy+1,lev,mon0)*(1-vrat)+temp(ix+1,iy+1,lev+1,mon0)*vrat
                wild2(3,2)=salt(ix+1,iy+1,lev,mon0)*(1-vrat)+salt(ix+1,iy+1,lev+1,mon0)*vrat
                wild2(4,1)=temp(ix,iy+1,lev,mon0)*(1-vrat)+temp(ix,iy+1,lev+1,mon0)*vrat
                wild2(4,2)=salt(ix,iy+1,lev,mon0)*(1-vrat)+salt(ix,iy+1,lev+1,mon0)*vrat
                if(in==1) then !nodes 1-3
                  tempout(i,k)=wild2(1,1)*arco(1)+wild2(2,1)*arco(2)+wild2(3,1)*arco(3)
                  saltout(i,k)=wild2(1,2)*arco(1)+wild2(2,2)*arco(2)+wild2(3,2)*arco(3)
                else
                  tempout(i,k)=wild2(1,1)*arco(1)+wild2(3,1)*arco(2)+wild2(4,1)*arco(3)
                  saltout(i,k)=wild2(1,2)*arco(1)+wild2(3,2)*arco(2)+wild2(4,2)*arco(3)
                endif

                !Check
                if(i==224547) then
                  write(87,*)'Node 224547:',k,lev,ix,iy,(ix-1)*iylen+iy,salt(ix,iy,lev,mon0), &
     &salt(ix+1,iy,lev,mon0),salt(ix+1,iy+1,lev,mon0),salt(ix,iy+1,lev,mon0), &
     &salt(ix,iy,lev+1,mon0)
                endif      
 
                if(tempout(i,k)<tempmin.or.saltout(i,k)<0) then
                  write(11,*)'T,S<0:',i,k,tempout(i,k),saltout(i,k),temp(ix,iy,lev,mon0)
                  stop
                endif

                !Enforce lower bound for temp. for eqstate
                tempout(i,k)=max(0.,tempout(i,k))

                !Correct with obs. T,S near Danish Str.
                if(iobs_T(i)==1.and.z(i,k)<-15) tempout(i,k)=min(tempout(i,k),T_obs(i,k))
                if(iobs_S(i)==1.and.z(i,k)<=-10) saltout(i,k)=max(saltout(i,k),S_obs1(i,k)) !Fehm region
                if(iobs_S(i)==2) saltout(i,k)=S_obs2(i,k) !Oder region
                if(iobs_S(i)==3.and.z(i,k)>=-10) saltout(i,k)=max(0.,saltout(i,k)-1) !Baltic

              enddo !k=1,nvrt
              cycle loop4
            endif !rat<small1
          enddo !iy=iylen1,iylen2-1
        enddo !ix=ixlen1,ixlen2-1
        if(ixy(i,1)==0.or.ixy(i,2)==0) then
          write(*,*)'Cannot find a parent element:',i
          tempout(i,:)=tem_outside
          saltout(i,:)=sal_outside
        endif
      end do loop4 !i=1,np
    
!     Note: although I swapped order of indices in SELFE for S,T,
!           I didn't do the same for them in this code; as long as 
!           the order of writing is correct in hotstart.in, it does not matter
      do i=1,ns
        n1=isidenode(i,1)
        n2=isidenode(i,2)
        do k=1,nvrt
          tsd(i,k)=(tempout(n1,k)+tempout(n2,k))/2
          ssd(i,k)=(saltout(n1,k)+saltout(n2,k))/2
        enddo !k
!       write(88,*)i,xcj(i),ycj(i),ssd(i,1),ssd(i,nvrt)
      enddo !i

      do i=1,ne
        n1=nm(i,1)
        n2=nm(i,2)
        n3=nm(i,3)
        do k=2,nvrt
          tsel(k,i,1)=(tempout(n1,k)+tempout(n2,k)+tempout(n3,k)+tempout(n1,k-1)+tempout(n2,k-1)+tempout(n3,k-1))/6
          tsel(k,i,2)=(saltout(n1,k)+saltout(n2,k)+saltout(n3,k)+saltout(n1,k-1)+saltout(n2,k-1)+saltout(n3,k-1))/6
        enddo !k
        tsel(1,i,1)=tsel(2,i,1) !mainly for hotstart format
        tsel(1,i,2)=tsel(2,i,2)
      enddo !i

!     Debug
      do i=1,np
        write(20,*)i,xl(i),yl(i),saltout(i,1)
        write(23,*)i,xl(i),yl(i),tempout(i,1)
      enddo !i
      do i=1,ne
        write(21,*)i,tsel(1,i,1)
      enddo !i

!     Output hotstart 
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

!      deallocate(lat,lon,zm,h,kbp,ihope,xind,yind,lind,salt,temp)

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

      end program gen_hot

      function signa(x1,x2,x3,y1,y2,y3)

      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2

      return
      end

