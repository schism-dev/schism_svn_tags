!     Compute geometric arrays used in SELFE (side info etc)
!     In the driver routine, the following 2 routines should be called
!     in the following fashion (all float arguments are real*4):
!
!     First get # of sides 'ns' with inputs: np (# of nodes),ne (# of elements), and elnode (connectivity table)
!     call compute_nside(np,ne,elnode,ns) 
!     Allocate side-related arrays
!     allocate(ic3(3,ne),elside(3,ne),isdel(2,ns),isidenode(2,ns),xcj(ns),ycj(ns))
!     Then compute the rest of side related arrays with additional inputs (xnd,ynd) (x,y coordinates of each node)
!     call selfe_geometry(np,ne,ns,xnd,ynd,elnode,ic3,elside,isdel,isidenode,xcj,ycj)
!     The outputs:
!                ns: # of sides
!                ic3(3,ne): 3 neighboring elements of an element
!                elside(3,ne): 3 sides of an element
!                isdel(2,ns): 2 adjacent elements of a side
!                isidenode(2,ns): 2 end nodes of a side
!                xcj(ns),ycj(ns): x,y of each side center

!     Common data used in 2 routines
      module selfe_geometry_mod
      implicit none
      public
      integer,save :: nx(4,4,3)
      integer,save,allocatable :: nne(:),indel(:,:),ic3(:,:)
      end module selfe_geometry_mod

      subroutine compute_nside(np,ne,elnode,ns)
      use selfe_geometry_mod
      implicit real(4)(a-h,o-z),integer(i-n)
      integer, intent(in) :: np,ne,elnode(3,ne)
      integer, intent(out) :: ns !,ic3(3,ne)
!      integer, allocatable :: indel(:,:)

      allocate(nne(np),ic3(3,ne),stat=istat)
      if(istat/=0) stop 'Failed to alloc. nne'

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
        do j=1,3
          nd=elnode(j,i)
          nne(nd)=nne(nd)+1
!          indel(nne(nd),nd)=i
        enddo
      enddo
      mnei=maxval(nne)

      allocate(indel(mnei,np),stat=istat)
      if(istat/=0) stop 'Failed to alloc. indel'
      nne=0
      do i=1,ne
        do j=1,3
          nd=elnode(j,i)
          nne(nd)=nne(nd)+1
          if(nne(nd)>mnei) then
            write(*,*)'Too many neighbors',nd
            stop
          endif
          indel(nne(nd),nd)=i
        enddo
      enddo !i

!     Compute ball info; this won't be affected by re-arrangement below
      do i=1,ne
        do j=1,3 
          ic3(j,i)=0 !index for bnd sides
          nd1=elnode(nx(3,j,1),i)
          nd2=elnode(nx(3,j,2),i)
          do k=1,nne(nd1)
            ie=indel(k,nd1)
            if(ie/=i.and.(elnode(1,ie)==nd2.or.elnode(2,ie)==nd2.or.elnode(3,ie)==nd2)) ic3(j,i)=ie
          enddo !k
        enddo !j
      enddo !i

      ns=0
      do ie=1,ne
        do j=1,3 !visit each side associated with element ie
          if(ic3(j,ie)==0.or.ie<ic3(j,ie)) then !new side
            ns=ns+1
          endif
        enddo !j
      enddo !ie

      end subroutine compute_nside

      subroutine selfe_geometry(np,ne,ns0,xnd,ynd,elnode,ic3_out,&
     &elside,isdel,isidenode,xcj,ycj)
      use selfe_geometry_mod
      implicit real(4)(a-h,o-z),integer(i-n)
      integer, intent(in) :: np,ne,ns0,elnode(3,ne)
      real, intent(in) :: xnd(np),ynd(np)
      integer, intent(out) :: ic3_out(3,ne),elside(3,ne),isdel(2,ns0),isidenode(2,ns0)
      real, intent(out) :: xcj(ns0),ycj(ns0)
      
      ic3_out=ic3

      ns=0 !# of sides
      do i=1,ne
        do j=1,3
          nd1=elnode(nx(3,j,1),i)
          nd2=elnode(nx(3,j,2),i)
          if(ic3(j,i)==0.or.i<ic3(j,i)) then !new sides
            ns=ns+1
            if(ns>ns0) then
              write(*,*)'Too many sides'
              stop
            endif
            elside(j,i)=ns
            isdel(1,ns)=i
            isidenode(1,ns)=nd1
            isidenode(2,ns)=nd2
            xcj(ns)=(xnd(nd1)+xnd(nd2))/2
            ycj(ns)=(ynd(nd1)+ynd(nd2))/2

            isdel(2,ns)=ic3(j,i) !bnd element => bnd side
!           Corresponding side in element ic3(j,i)
            if(ic3(j,i)/=0) then !old internal side
              iel=ic3(j,i)
              index=0
              do k=1,3
                if(ic3(k,iel)==i) then
                  index=k
                  exit
                endif
              enddo !k
              if(index==0) then
                write(*,*)'Wrong ball info',i,j
                stop
              endif
              elside(index,iel)=ns
            endif !ic3(j,i).ne.0
          endif !ic3(j,i)==0.or.i<ic3(j,i)
        enddo !j
      enddo !i=1,ne
      if(ns/=ns0) stop 'Side count mismatch'

      end subroutine selfe_geometry
