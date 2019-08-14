! Generate hotstart.in with 0 elev. and vel (i.e. S,T only; no tracers)
! Need to update T,S below.
!
! ifort -O2 -Bstatic -o change_hotstart5 change_hotstart5.f90

!   Input: 
!     (1) hgrid.gr3;
!     (2) vgrid.in; 
!     (3) estuary.gr3 (flags): depth=0: outside; =1: inside
!   Output: hotstart.in 

!      implicit none

!      integer, parameter :: debug=1
      integer, parameter :: mnp=50000
      integer, parameter :: mne=100000
      integer, parameter :: mns=150000
      integer, parameter :: mnv=60
      integer, parameter :: mnope=10 !max # of open bnd segments
      integer, parameter :: mnond=1000 !max # of open bnd nodes in each segment
      integer, parameter :: mnei=20 !neighbor
      integer, parameter :: nbyte=4
      real, parameter :: small1=1.e-3 !used to check area ratios
  
!     Vertical postion, salinity, and temperature
      integer, dimension(:,:), allocatable :: kbp

      integer :: ier ! allocate error return.

      dimension xnd(mnp),ynd(mnp),nm(mne,4),dp(mnp),i34(mne)
      dimension ztot(0:mnv),sigma(mnv),cs(mnv),z(mnp,mnv),iest(mnp),ixy(mnp,2),arco(3)
      dimension wild(100),wild2(100,2)
      dimension tempout(mnp,mnv), saltout(mnp,mnv),month_day(12)
      dimension tsd(mns,mnv),ssd(mns,mnv),tsel(mnv,mne,2)
      dimension nne(mnp),ine(mnp,mnei),ic3(mne,4),nx(4,4,3),js(mne,4),is(mns,2),isidenode(mns,2)
      dimension xcj(mns),ycj(mns),nond(mnope),iond(mnope,mnond),iob(mnope),iond2(mnope*mnond)

!     Read in hgrid and vgrid
      open(17,file='estuary.gr3',status='old')
      open(14,file='hgrid.gr3',status='old') !only need depth info and connectivity
      open(19,file='vgrid.in',status='old')
      read(14,*)
      read(14,*)ne,np
      if(np.gt.mnp.or.ne.gt.mne) then
        write(*,*)'Increase mnp/mne'
        stop
      endif
      read(17,*)
      read(17,*)
      do i=1,np
        read(14,*)j,xnd(i),ynd(i),dp(i)
        read(17,*)j,xtmp,ytmp,iest(i)
        if(iest(i)/=0.and.iest(i)/=1) then
          write(11,*)'Estuary flag wrong:',i,iest(i)
          stop
        endif
      enddo !i
      do i=1,ne
        read(14,*)j,i34(i),(nm(i,l),l=1,3)
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

!       Z-levels; shallow pts have junk values
        do k=1,kz-1
          z(i,k)=ztot(k)
        enddo !k
      enddo !i

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
            xcj(ns)=(xnd(nd1)+xnd(nd2))/2
            ycj(ns)=(ynd(nd1)+ynd(nd2))/2

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

!     Compute node T,S
      tempout=10
      y1=4109420 !mouth
      y2=4200499 !Potomac
      y3=4324141 !CB3.3
      y4=4363309 !upper
      do i=1,np
        if(dp(i)<=0) then
          saltout(i,:)=0; cycle
        endif

        if(iest(i)==0) then !outside Bay
          S_bot=34.7
          S_sur=34.7       
        else if(ynd(i)<=y2) then
          S_bot=min(34.7,34.7-4.7*(ynd(i)-y1)/(y2-y1))
          S_sur=min(34.7,34.7-6.7*(ynd(i)-y1)/(y2-y1))
        else if(ynd(i)<=y3) then
          S_bot=30-3*(ynd(i)-y2)/(y3-y2)
          S_sur=28-8*(ynd(i)-y2)/(y3-y2)
        else
          S_bot=max(17.,27-10*(ynd(i)-y3)/(y4-y3))
          S_sur=max(10.,20-10*(ynd(i)-y3)/(y4-y3))
        endif
 
        do k=1,nvrt
          saltout(i,k)=S_sur-(S_bot-S_sur)*z(i,k)/dp(i)
        enddo !k
        write(27,*)i,saltout(i,1:nvrt)
      enddo !i=1,np
      write(*,*)'S max/min=',maxval(saltout(1:np,1:nvrt)),minval(saltout(1:np,1:nvrt))

!     Output hotstart 
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

!     MPI SELFE
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

      stop
      end
