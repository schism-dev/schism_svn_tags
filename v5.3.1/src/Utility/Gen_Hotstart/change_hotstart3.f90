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

!===============================================================================
! Manipulate hotstart.in

! Inputs:
!        hotstart.old (unformatted binary); screen inputs.
! Output: hotstart.in (unformatted binary). This format is different
!         between Intel and AMD!
!
!  ifort -Bstatic -O3 -assume byterecl -o change_hotstart3 change_hotstart3.f90

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
!-------------------------------------------------------------------------------
      
  print*, 'Input ne, np, ns, nvrt, and ntracers:'
  read*, ne_global,np_global,ns_global,nvrt,ntracers
!  ne_global=118707
!  np_global=60250
!  ns_global=178963
!  nvrt=18
!  ntracers=0

  allocate(idry_e(ne_global),we(nvrt,ne_global),tsel(2,nvrt,ne_global), &
           idry_s(ns_global),su2(nvrt,ns_global),sv2(nvrt,ns_global), &
           tsd(nvrt,ns_global),ssd(nvrt,ns_global), &
           idry(np_global),eta2(np_global),tnd(nvrt,np_global),snd(nvrt,np_global), &
           tem0(nvrt,np_global),sal0(nvrt,np_global),q2(np_global,nvrt), &
           xl(np_global,nvrt),dfv(np_global,nvrt),dfh(np_global,nvrt), &
           dfq1(np_global,nvrt),dfq2(np_global,nvrt),qnon(nvrt,np_global),stat=istat)
  if(istat/=0) stop 'Allocation error (2)'
  if(ntracers==0) then
    allocate(trel0(1,1,1),trel(1,1,1),stat=istat)
  else
    allocate(trel0(ntracers,nvrt,ne_global),trel(ntracers,nvrt,ne_global),stat=istat)
  endif
  if(istat/=0) stop 'Allocation error (3)'
!-------------------------------------------------------------------------------
! Read hotstart files
!-------------------------------------------------------------------------------
  open(36,file='hotstart.old',form='unformatted',status='old')
  open(37,file='hotstart.in',form='unformatted',status='replace')
  read(36) time,it,ifile
  write(37) time,it,ifile
  do i=1,ne_global
    read(36) j,idry_e(i),(we(j,i),tsel(1:2,j,i),(trel0(l,j,i),trel(l,j,i),l=1,ntracers),j=1,nvrt)
    idry_e(i)=0
    tsel(1,:,i)=max(tsel(1,:,i),10.d0)
    tsel(2,:,i)=max(tsel(2,:,i),0.d0)
    write(37) i,idry_e(i),(we(j,i),tsel(1:2,j,i),(trel0(l,j,i),trel(l,j,i),l=1,ntracers),j=1,nvrt)
  enddo !i
  do i=1,ns_global
    read(36) j,idry_s(i),(su2(j,i),sv2(j,i),tsd(j,i),ssd(j,i),j=1,nvrt)
    idry_s(i)=0
    tsd(:,i)=max(tsd(:,i),10.d0)
    ssd(:,i)=max(ssd(:,i),0.d0)
    write(37) i,idry_s(i),(su2(j,i),sv2(j,i),tsd(j,i),ssd(j,i),j=1,nvrt)
  enddo !i
  do i=1,np_global
    read(36) j,eta2(i),idry(i),(tnd(j,i),snd(j,i),tem0(j,i),sal0(j,i),q2(i,j),xl(i,j), &
             dfv(i,j),dfh(i,j),dfq1(i,j),dfq2(i,j),qnon(j,i),j=1,nvrt)
    idry(i)=0
    tnd(:,i)=max(tnd(:,i),10.d0)
    tem0(:,i)=max(tem0(:,i),10.d0)
    snd(:,i)=max(snd(:,i),0.d0)
    sal0(:,i)=max(sal0(:,i),0.d0)
    eta2(i)=max(eta2(i),-2.d0)
    write(37) i,eta2(i),idry(i),(tnd(j,i),snd(j,i),tem0(j,i),sal0(j,i),q2(i,j),xl(i,j), &
             dfv(i,j),dfh(i,j),dfq1(i,j),dfq2(i,j),qnon(j,i),j=1,nvrt)
  enddo !i
  close(36)
  close(37)

! Output

end program combine_hotstart1
