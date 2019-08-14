!     Manipulate hotstart.in
!     Input: hotstart.old
!     Output: hotstart.new
!     ifort -Bstatic -O3 -o change_hotstart change_hotstart.f90
      program change_hotstart
      implicit real*8(a-h,o-z)
!      integer, parameter :: mnp=20000
!      integer, parameter :: mne=36000
!      integer, parameter :: mns=55000
!      integer, parameter :: mnv=50
      integer, parameter :: nbyte=4
      character(len=12) :: ifile_char
      allocatable :: z_r(:),tem1(:),sal1(:),idry_e(:),we(:,:),idry_s(:),eta2(:),idry(:)
      allocatable :: tnd(:,:),snd(:,:),tem0(:,:),sal0(:,:),q2(:,:),xl(:,:),su2(:,:),sv2(:,:),tsd(:,:),ssd(:,:)
      !allocatable :: z_r(mnv),tem1(mnv),sal1(mnv),idry_e(mne),we(mne,mnv),idry_s(mns),eta2(mnp),idry(mnp)
      !real*8, dimension(mnp,mnv) :: tnd,snd,tem0,sal0,q2,xl
      !real*8, dimension(mns,mnv) :: su2,sv2,tsd,ssd

!     Input np etc.
      np=18757
      ne=35234
      ns=54005
      nvrt=26

      allocate(z_r(nvrt),tem1(nvrt),sal1(nvrt),idry_e(ne),we(ne,nvrt),idry_s(ns),eta2(np), &
     &idry(np),tnd(np,nvrt),snd(np,nvrt),tem0(np,nvrt),sal0(np,nvrt),q2(np,nvrt),xl(np,nvrt), &
     &su2(ns,nvrt),sv2(ns,nvrt),tsd(ns,nvrt),ssd(ns,nvrt),stat=istat)
     if(istat/=0) stop 'Failed to allocate'
!      if(np>mnp.or.ne>mne.or.ns>mns.or.nvrt>mnv) then
!        print*, 'Increase dimensions'
!        stop
!      endif

!     Record length for hot start files (double precision for all reals)
      ihot_len=nbyte*(5+6*mnv+(2*nvrt+1)*ne+ns+(8*nvrt+1)*ns+3*np+12*np*nvrt+1)+12

      open(36,file='hotstart.old',access='direct',recl=ihot_len)
      read(36,rec=1)iths,time,icst,nz_r,(z_r(k),tem1(k),sal1(k),k=1,mnv),(idry_e(i),(we(i,j),j=1,nvrt),i=1,ne), &
     &(idry_s(i),(su2(i,j),sv2(i,j),tsd(i,j),ssd(i,j),j=1,nvrt),i=1,ns), &
     &(eta2(i),idry(i),(tnd(i,j),snd(i,j),tem0(i,j),sal0(i,j),q2(i,j),xl(i,j),j=1,nvrt),i=1,np), &
     &ifile,ifile_char

!     New record length for output
!      ihot_len=nbyte*(3+(2*nvrt+1)*ne+ns+(8*nvrt+1)*ns+3*np+12*np*nvrt+1)+12
      tnd=10; tsd=10; tem0=10
      open(37,file='hotstart.new',access='direct',recl=ihot_len)
      write(37,rec=1)iths,time,icst,nz_r,(z_r(k),tem1(k),sal1(k),k=1,mnv),(idry_e(i),(we(i,j),j=1,nvrt),i=1,ne), &
     &(idry_s(i),(su2(i,j),sv2(i,j),tsd(i,j),ssd(i,j),j=1,nvrt),i=1,ns), &
     &(eta2(i),idry(i),(tnd(i,j),snd(i,j),tem0(i,j),sal0(i,j),q2(i,j),xl(i,j),j=1,nvrt),i=1,np), &
     &ifile,ifile_char

      stop
      end

