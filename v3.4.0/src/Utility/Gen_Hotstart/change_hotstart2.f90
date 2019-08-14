!     Manipulate hotstart.in
!     Input: hotstart.old; hgrid.gr3
!     Output: hotstart.new
      program change_hotstart
      implicit real*8(a-h,o-z)
      integer, parameter :: mnp=40000
      integer, parameter :: mne=80000
      integer, parameter :: mns=90000
      integer, parameter :: mnv=70
      integer, parameter :: nbyte=4
      character(len=12) :: ifile_char
      dimension tem1(mnv),sal1(mnv),idry_e(mne),we(mnv,mne),idry_s(mns),eta2(mnp),idry(mnp)
      dimension tsel(mnv,mne,2),x(mnp),y(mnp),elnode(3,mne)
      real*8, dimension(mnv,mnp) :: tnd,snd,tem0,sal0
      real*8, dimension(mnv,mns) :: su2,sv2,tsd,ssd
      real*8, dimension(mnp,mnv) :: q2,xl,dfv,dfh,dfq1,dfq2

!     Input np etc.
      ns=59884
      nvrt=54
      open(14,file='hgrid.gr3',status='old')
      read(14,*)
      read(14,*) ne,np    
      if(np>mnp.or.ne>mne.or.ns>mns.or.nvrt>mnv) then
        print*, 'Increase dimensions'
        stop
      endif
      do i=1,np
        read(14,*) j,x(i),y(i) !,dp(i)
      enddo !i
      do i=1,ne
        read(14,*)j,l,(elnode(k,i),k=1,3)
      enddo !i
      close(14)

!     Old hotstart.in
      ihot_len=nbyte*(3+(2*nvrt+1)*ne+(8*nvrt+1)*ns+3*np+20*np*nvrt+1)+12
      open(36,file='hotstart.old',access='direct',recl=ihot_len)

      read(36,rec=1)iths,time,(idry_e(i),(we(j,i),j=1,nvrt),i=1,ne), &
     &(idry_s(i),(su2(j,i),sv2(j,i),tsd(j,i),ssd(j,i),j=1,nvrt),i=1,ns), &
     &(eta2(i),idry(i),(tnd(j,i),snd(j,i),tem0(j,i),sal0(j,i),q2(i,j),xl(i,j), &
     &dfv(i,j),dfh(i,j),dfq1(i,j),dfq2(i,j),j=1,nvrt),i=1,np),ifile,ifile_char

!     Compute tsel
      do i=1,ne
        n1=elnode(1,i)
        n2=elnode(2,i)
        n3=elnode(3,i)
        do k=2,nvrt
          tsel(k,i,1)=(tem0(k,n1)+tem0(k,n2)+tem0(k,n3)+tem0(k-1,n1)+tem0(k-1,n2)+tem0(k-1,n3))/6
          tsel(k,i,2)=(sal0(k,n1)+sal0(k,n2)+sal0(k,n3)+sal0(k-1,n1)+sal0(k-1,n2)+sal0(k-1,n3))/6
        enddo !k
        tsel(1,i,1)=tsel(2,i,1) !mainly for hotstart format
        tsel(1,i,2)=tsel(2,i,2)
      enddo !i

!     Output
      ihot_len=ihot_len+4*nvrt*ne*nbyte
      open(37,file='hotstart.new',access='direct',recl=ihot_len)
      write(37,rec=1)iths,time,(idry_e(i),(we(j,i),tsel(j,i,1),tsel(j,i,2),j=1,nvrt),i=1,ne), &
     &(idry_s(i),(su2(j,i),sv2(j,i),tsd(j,i),ssd(j,i),j=1,nvrt),i=1,ns), &
     &(eta2(i),idry(i),(tnd(j,i),snd(j,i),tem0(j,i),sal0(j,i),q2(i,j),xl(i,j), &
     &dfv(i,j),dfh(i,j),dfq1(i,j),dfq2(i,j),j=1,nvrt),i=1,np),ifile,ifile_char
      stop
      end

