!===============================================================================
! Manipulate hotstart.in for MPI versions (v3.0a and up): only vgrid (pure S) is
! allowed to change.

! Inputs:
!        hotstart.old (unformatted binary); screen inputs.
! Output: hotstart.in (unformatted binary). This format is different
!         between Intel and AMD!
!
!  ifort -Bstatic -O3 -assume byterecl -o change_hotstart4 change_hotstart4.f90

!===============================================================================

program combine_hotstart1
!-------------------------------------------------------------------------------
  implicit real(8)(a-h,o-z),integer(i-n)
  parameter(nbyte=4)
!  character(12) :: it_char
!  character(72) :: fgb,fgb2,fdb  ! Processor specific global output file name
!  integer :: lfgb,lfdb       ! Length of processor specific global output file name
!  allocatable ner(:),npr(:),nsr(:)
!  allocatable ielg(:,:),iplg(:,:),islg(:,:)
  allocatable idry_e(:),we(:,:),tsel(:,:,:),idry_s(:),su2(:,:),sv2(:,:)
  allocatable tsd(:,:),ssd(:,:),idry(:),eta2(:),tnd(:,:),snd(:,:)
  allocatable tem0(:,:),sal0(:,:),q2(:,:),xl(:,:),dfv(:,:),dfh(:,:)
  allocatable dfq1(:,:),dfq2(:,:),qnon(:,:),trel0(:,:,:),trel(:,:,:)
  allocatable intv(:),zrat(:),swild(:,:)
!-------------------------------------------------------------------------------
      
  print*, 'Input ne, np, ns, nvrt, and ntracers:'
  read*, ne_global,np_global,ns_global,nvrt,ntracers
  print*, 'Input new # of levels:'
  read*, nvrt1

  allocate(idry_e(ne_global),we(nvrt,ne_global),tsel(2,nvrt,ne_global), &
           idry_s(ns_global),su2(nvrt,ns_global),sv2(nvrt,ns_global), &
           tsd(nvrt,ns_global),ssd(nvrt,ns_global), &
           idry(np_global),eta2(np_global),tnd(nvrt,np_global),snd(nvrt,np_global), &
           tem0(nvrt,np_global),sal0(nvrt,np_global),q2(np_global,nvrt), &
           xl(np_global,nvrt),dfv(np_global,nvrt),dfh(np_global,nvrt), &
           dfq1(np_global,nvrt),dfq2(np_global,nvrt),qnon(nvrt,np_global), &
           swild(nvrt1,max(11,3+ntracers)),intv(nvrt1),zrat(nvrt1),stat=istat)
  if(istat/=0) stop 'Allocation error (2)'
  if(ntracers==0) then
    allocate(trel0(1,1,1),trel(1,1,1),stat=istat)
  else
    allocate(trel0(ntracers,nvrt,ne_global),trel(ntracers,nvrt,ne_global),stat=istat)
  endif
  if(istat/=0) stop 'Allocation error (3)'

! Calculate interplation coefficients, assuming equal distance in vertical
  do k=1,nvrt1
   sig=-1+(k-1.)/(nvrt1-1)
   tmp=(sig+1)*(nvrt-1)+1
   intv(k)=int(tmp)
   if(intv(k)==nvrt) intv(k)=nvrt-1
   zrat(k)=tmp-intv(k) !from level intv(k)
  enddo !k

!-------------------------------------------------------------------------------
! Read hotstart files
!-------------------------------------------------------------------------------
  open(36,file='hotstart.old',form='unformatted',status='old')
  open(37,file='hotstart.in',form='unformatted',status='replace')
  read(36) time,it,ifile
  write(37) time,it,ifile
  do i=1,ne_global
    read(36) j,idry_e(i),(we(j,i),tsel(1:2,j,i),(trel0(l,j,i),trel(l,j,i),l=1,ntracers),j=1,nvrt)
    do j=1,nvrt1
      k=intv(j)
      swild(j,1)=we(k,i)*(1-zrat(k))+we(k+1,i)*zrat(k)
      swild(j,2:3)=tsel(1:2,k,i)*(1-zrat(k))+tsel(1:2,k+1,i)*zrat(k)
!      swild(j,4:3+ntracers)=trel(1:ntracers,k,i)*(1-zrat(k))+trel(1:ntracers,k+1,i)*zrat(k)
    enddo !j
    write(37) i,idry_e(i),(swild(j,1:3),j=1,nvrt1)
    !write(37) i,idry_e(i),(swild(j,1),swild(j,2:3),swild(j,4:3+ntracers),swild(j,4:3+ntracers),j=1,nvrt)
  enddo !i
  do i=1,ns_global
    read(36) j,idry_s(i),(su2(j,i),sv2(j,i),tsd(j,i),ssd(j,i),j=1,nvrt)
    do j=1,nvrt1
      k=intv(j)
      swild(j,1)=su2(k,i)*(1-zrat(k))+su2(k+1,i)*zrat(k)
      swild(j,2)=sv2(k,i)*(1-zrat(k))+sv2(k+1,i)*zrat(k)
      swild(j,3)=tsd(k,i)*(1-zrat(k))+tsd(k+1,i)*zrat(k)
      swild(j,4)=ssd(k,i)*(1-zrat(k))+ssd(k+1,i)*zrat(k)
    enddo !j
    write(37) i,idry_s(i),(swild(j,1:4),j=1,nvrt1)
  enddo !i
  do i=1,np_global
    read(36) j,eta2(i),idry(i),(tnd(j,i),snd(j,i),tem0(j,i),sal0(j,i),q2(i,j),xl(i,j), &
             dfv(i,j),dfh(i,j),dfq1(i,j),dfq2(i,j),qnon(j,i),j=1,nvrt)
    do j=1,nvrt1
      k=intv(j)
      swild(j,1)=tnd(k,i)*(1-zrat(k))+tnd(k+1,i)*zrat(k)
      swild(j,2)=snd(k,i)*(1-zrat(k))+snd(k+1,i)*zrat(k)
      swild(j,3)=tem0(k,i)*(1-zrat(k))+tem0(k+1,i)*zrat(k)
      swild(j,4)=sal0(k,i)*(1-zrat(k))+sal0(k+1,i)*zrat(k)
      swild(j,5)=q2(i,k)*(1-zrat(k))+q2(i,k+1)*zrat(k)
      swild(j,6)=xl(i,k)*(1-zrat(k))+xl(i,k+1)*zrat(k)
      swild(j,7)=dfv(i,k)*(1-zrat(k))+dfv(i,k+1)*zrat(k)
      swild(j,8)=dfh(i,k)*(1-zrat(k))+dfh(i,k+1)*zrat(k)
      swild(j,9)=dfq1(i,k)*(1-zrat(k))+dfq1(i,k+1)*zrat(k)
      swild(j,10)=dfq2(i,k)*(1-zrat(k))+dfq2(i,k+1)*zrat(k)
      swild(j,11)=qnon(k,i)*(1-zrat(k))+qnon(k+1,i)*zrat(k)
    enddo !j
    write(37) i,eta2(i),idry(i),(swild(j,1:11),j=1,nvrt1)
  enddo !i
  close(36)
  close(37)

! Output

end program combine_hotstart1
