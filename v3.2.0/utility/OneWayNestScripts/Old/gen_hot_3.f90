!     Generate hotstart.in for estuary only runs (MPI version 3.0b and up).
!     New from previous version: added non-hydrostatic pressure

!   Input: 
!     (1) nvrt,nond,ntracers,start_slope(1:2): # of vertical levels; # of nodes at mouth 
!                                              or used in salt3D.th); # of tracers; S-T slope
!         info for starting CORIE day from st_slope.txt.
!     (2) hgrid.gr3
!     (3) estuary.gr3 (must have 2 anchors): depth=0: outside; =1: inside; =-1: 1st anchor; =-2: 2nd anchor ;
!                                            1st anchor is near the mouth; 2nd in river.
!     (4) salt3D.th (binary): mouth S;
!     (5) temp.th (ASCII): river temp;
!   Output: hotstart.in (unformatted)

!   ifort -Bstatic -O3 -assume byterecl -o gen_hot_3_canopus gen_hot_3.f90

      program gen_hot
      implicit real*8(a-h,o-z)

!      integer, parameter :: mnp=60000
!      integer, parameter :: mne=100000
!      integer, parameter :: mns=130000
!      integer, parameter :: mnv=90
!      integer, parameter :: mnei=20 !neighbor
      integer, parameter :: nbyte=4
      real*8, parameter :: small1=1.e-3 !used to check area ratios
  
!      dimension x(mnp),y(mnp),nm(mne,4),dp(mnp),i34(mne),iest(mnp)
!      dimension tempout(mnp,mnv), saltout(mnp,mnv)
!      dimension tsd(mns,mnv),ssd(mns,mnv),tsel(mnv,mne,2)
!      dimension nne(mnp),ine(mnp,mnei),ic3(mne,4),nx(4,4,3),js(mne,4),is(mns,2),isidenode(mns,2)
!      dimension xcj(mns),ycj(mns),start_slope(2)
!      real*4 dt_th,s3D(mnp,mnv)

      dimension nx(4,4,3),start_slope(2)
      allocatable :: x(:),y(:),nm(:,:),dp(:),i34(:),iest(:)
      allocatable :: tempout(:,:), saltout(:,:),tsel(:,:,:),nne(:),ic3(:,:),js(:,:)
      allocatable :: tsd(:,:),ssd(:,:),ine(:,:),is(:,:),isidenode(:,:),xcj(:),ycj(:)
      real*4, allocatable :: s3D(:,:)
      real*4 dt_th

      print*, 'Input nvrt,nond, and starting S-T slope and intersection:'
!'
      read*, nvrt,nond,start_slope(1:2)
      ntracers=0
!      if(nvrt>mnv.or.nond>mnp) then
!        write(*,*)'Increase mnv'
!        stop
!      endif

!     Open temp.th for river temp.
      open(9,file='temp.th',status='old')
      read(9,*)dt_th,temp_ri
      close(9)

!     Read in hgrid
      open(14,file='hgrid.gr3',status='old')
      open(17,file='estuary.gr3',status='old')
      read(14,*)
      read(14,*)ne,np

      if(nond>np) stop 'Impossible'

!     Allocate arrays
      allocate(x(np),y(np),nm(ne,4),dp(np),i34(ne),iest(np),tempout(np,nvrt), &
     &saltout(np,nvrt),tsel(nvrt,ne,2),nne(np),ic3(ne,4),js(ne,4),s3D(np,nvrt),stat=istat)
      if(istat/=0) stop 'Failed to allocate (1)'

      read(17,*)
      read(17,*)
      icount=0 !# of 1st anchor pt
      icount2=0
      do i=1,np
        read(14,*)j,x(i),y(i),dp(i)
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

      if(icount/=1.or.icount2/=1) then
        write(11,*)'# of 1st or 2nd anchor pt/=1:',icount,icount2
        stop
      endif
      if(ianchor1==ianchor2) then
        write(11,*)'Duplicate anchors'
        stop
      endif

      do i=1,ne
        read(14,*)j,i34(i),(nm(i,l),l=1,i34(i))
        if(i34(i)/=3) then
          write(*,*)'Triangles only please:',i34(i)
          stop
        endif
      enddo !i
      close(14)
      close(16)

!     Open temp3D.th & salt3D.th for mouth T,S
      nrecl_te=nbyte*(1+nond*nvrt)
!      open(9,file='temp3D.th',access='direct',recl=nrecl_te,status='old')
      open(8,file='salt3D.th',access='direct',recl=nrecl_te,status='old')
!      read(9,rec=1)dt_th,((t3D(i,k),k=1,nvrt),i=1,nond)
      read(8,rec=1)dt_th,((s3D(i,k),k=1,nvrt),i=1,nond)
      ianchor0=max(1,nond/2)
!      do i=1,nond
!        read(9,*)node,(t3D(node,k),k=1,nvrt)
!        read(8,*)node,(s3D(node,k),k=1,nvrt)
!        if(node>mnp) then
!          write(11,*)'Out of bound:',node
!          stop
!        endif
!       Choose mid node as anchor pt for interpolating S,T
!        if(i==nond/2.or.nond==1) ianchor0=node
!      enddo !i
      close(8)
      close(9)

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

      nne=0
      do i=1,ne
        do j=1,i34(i)
          nd=nm(i,j)
          nne(nd)=nne(nd)+1
        enddo
      enddo
      mnei=maxval(nne)
      print*, 'mnei = ',mnei
 
      allocate(ine(np,mnei),stat=istat)
      if(istat/=0) stop 'Failed to allocate (2)'

      nne=0
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
          endif !ic3(i,j)==0.or.i<ic3(i,j)
        enddo !j=1,i34
      enddo !i=1,ne

      if(ns<ne.or.ns<np) then
        write(11,*)'Weird grid with ns < ne or ns < np',np,ne,ns
        stop
      endif

!     Allocate
      ns0=ns
      allocate(tsd(ns,nvrt),ssd(ns,nvrt),is(ns,2),isidenode(ns,2),xcj(ns),ycj(ns),stat=istat)
      if(istat/=0) stop 'Failed to allocate (3)'

      ns=0 !# of sides
      do i=1,ne
        do j=1,i34(i)
          nd1=nm(i,nx(i34(i),j,1))
          nd2=nm(i,nx(i34(i),j,2))
          if(ic3(i,j)==0.or.i<ic3(i,j)) then !new sides
            ns=ns+1
            if(ns>ns0) then
              write(11,*)'Too many sides'
              stop
            endif
            js(i,j)=ns
            is(ns,1)=i
            isidenode(ns,1)=nd1
            isidenode(ns,2)=nd2
            xcj(ns)=(x(nd1)+x(nd2))/2
            ycj(ns)=(y(nd1)+y(nd2))/2
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

!     Compute initial S,T
      do i=1,np
        if(x(ianchor2)-x(ianchor1)==0) then
          write(11,*)'Wrong anchor pts:',ianchor1,ianchor2
          stop
        endif
        xrat=(x(i)-x(ianchor1))/(x(ianchor2)-x(ianchor1))
        xrat=max(0.,min(1.,xrat))
        do k=1,nvrt
          saltout(i,k)=s3D(ianchor0,k)*(1-xrat) !+0
          temp_ocean=start_slope(1)*s3D(ianchor0,k)+start_slope(2)
          tempout(i,k)=temp_ocean*(1-xrat)+temp_ri*xrat
        enddo !k 
      enddo !i
    
!       Note: although I swapped order of indices in SELFE for S,T,
!             I didn't do the same for them in this code; as long as 
!             the order of writing is correct in hotstart.in, it does not matter
        do i=1,ns
          n1=isidenode(i,1)
          n2=isidenode(i,2)
          do k=1,nvrt
            tsd(i,k)=(tempout(n1,k)+tempout(n2,k))/2
            ssd(i,k)=(saltout(n1,k)+saltout(n2,k))/2
          enddo !k
!          write(88,*)i,xcj(i),ycj(i),ssd(i,1),ssd(i,nvrt)
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

      print*, 'Finished'

      stop
      end program gen_hot

