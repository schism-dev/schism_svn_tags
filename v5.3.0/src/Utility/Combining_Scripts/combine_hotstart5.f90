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
! Read in rank-specific hotstart outputs and combine them into hotstart.in.
! Gobal-local mappings are read in from separate files.

! Inputs:
!        rank-specific hotstart files; 
!        local_to_global_*; 
!        screen: nproc, ntracers; it_char (iteration #)
! Output: hotstart.in (unformatted binary). This format is different
!         between Intel and AMD!
!
!  PGI compiler:
!  pgf90 -O2 -mcmodel=medium  -Bstatic -o combine_hotstart5 ../UtilLib/argparse.F90 combine_hotstart5.f90

! Revisions: v4 with SED2D & SED parts added
!================================================================================
program combine_hotstart
integer :: istep
integer :: nproc
integer :: ntracer
call combine_hotstart5_input(istep,nproc,ntracer)
call combine_hotstart5(istep,nproc,ntracer)
end program combine_hotstart

subroutine combine_hotstart5_input(istep, nproc, ntracer)
use argparse
implicit none
integer       :: istep   !< global iteration/step 
integer       :: nproc   !< number of processors in original simulation
integer       :: ntracer !< number of tracers 

cmd_name = "combine_hotstart5"
call cla_init(cmd_name)
call cla_register('-i','--iteration', 'global iteration to combine (before _00*_hotstart', cla_int,'')
call cla_register('-p','--nproc', 'number of procs to combine', cla_int,'1')
call cla_register('-t','--ntracer','number of tracers', cla_int  ,'0')
call cla_validate
        
call cla_get("--iteration",istep)
call cla_get("--nproc",nproc)
call cla_get("--ntracer",ntracer)
end subroutine


subroutine combine_hotstart5(istep, nproc, ntracers)
!-------------------------------------------------------------------------------
  implicit real(8)(a-h,o-z),integer(i-n)
  parameter(nbyte=4)
  character(12) :: it_char
  character(72) :: fgb,fgb2,fdb  ! Processor specific global output file name
  integer :: lfgb,lfdb       ! Length of processor specific global output file name
  allocatable ner(:),npr(:),nsr(:) !resident (no ghosts)
  allocatable ielg(:,:),iplg(:,:),islg(:,:)
  allocatable idry_e(:),we(:,:),tsel(:,:,:),idry_s(:),su2(:,:),sv2(:,:)
  allocatable tsd(:,:),ssd(:,:),idry(:),eta2(:),tnd(:,:),snd(:,:)
  allocatable tem0(:,:),sal0(:,:),q2(:,:),xl(:,:),dfv(:,:),dfh(:,:)
  allocatable dfq1(:,:),dfq2(:,:),qnon(:,:),trel0(:,:,:),trel(:,:,:)
  allocatable iwild(:),swild(:,:,:),wild2(:),wild3(:),wild4(:)
  allocatable ihotstp(:),bed(:,:,:),bed_frac(:,:,:),rough_p(:)
!-------------------------------------------------------------------------------
      
!-------------------------------------------------------------------------------
! Aquire user inputs
!-------------------------------------------------------------------------------

  allocate(ner(0:nproc-1),npr(0:nproc-1),nsr(0:nproc-1),ihotstp(0:nproc-1),stat=istat) 
  if(istat/=0) stop 'Allocation error'

! Read mapping info
  fdb='local_to_global_0000'; fdb=adjustl(fdb)
  lfdb=len_trim(fdb)
  mxner=0; mxnpr=0; mxnsr=0
  do irank=0,nproc-1
    write(fdb(lfdb-3:lfdb),'(i4.4)') irank
    open(10,file=fdb,status='old')
    read(10,*)ne_global,np_global,nvrt,ntmp
    if(ntmp/=nproc) then
      print*, '# of proc mistmatch!',ntmp,nproc,ne_global,np_global,nvrt
      stop
    endif
    read(10,*) 
    read(10,*)ner(irank); do i=1,ner(irank); read(10,*); enddo;
    read(10,*)npr(irank); do i=1,npr(irank); read(10,*); enddo;
    read(10,*)nsr(irank); do i=1,nsr(irank); read(10,*); enddo;
    mxner=max0(mxner,ner(irank))
    mxnpr=max0(mxnpr,npr(irank))
    mxnsr=max0(mxnsr,nsr(irank))
    close(10)
  enddo !irank

  allocate(ielg(0:nproc-1,mxner),iplg(0:nproc-1,mxnpr),islg(0:nproc-1,mxnsr), &
          stat=istat)
  if(istat/=0) stop 'Allocation error (3)'

  do irank=0,nproc-1
    write(fdb(lfdb-3:lfdb),'(i4.4)') irank
    open(10,file=fdb,status='old')
    read(10,*); read(10,*)
    read(10,*)ner(irank)
    do i=1,ner(irank)
      read(10,*)j,ielg(irank,i)
    enddo
    read(10,*)npr(irank)
    do i=1,npr(irank)
      read(10,*)j,iplg(irank,i)
    enddo
    read(10,*)nsr(irank)
    do i=1,nsr(irank)
      read(10,*)j,islg(irank,i)
    enddo
    close(10)
  enddo !irank

! Compute global side #
  ns_global=0
  do irank=0,nproc-1
    do i=1,nsr(irank)
      ns_global=max0(ns_global,islg(irank,i))
    enddo !i 
  enddo !irank
  print*, 'Global quantities:',ne_global,np_global,ns_global

!-------------------------------------------------------------------------------
! Read hotstart files
!-------------------------------------------------------------------------------

  ! Open file
  write(it_char,'(i)') istep
  it_char=adjustl(it_char)  !place blanks at end
  it_len=len_trim(it_char)  !length without trailing blanks
  fgb=it_char(1:it_len)//'_0000'; lfgb=len_trim(fgb);
  print*, fgb
  !Gather all ranks
  ! Output
  open(37,file='hotstart.in',form='unformatted',status='replace')

  !Read/write elements first
  allocate(iwild(ne_global),swild(3+2*ntracers,nvrt,ne_global),stat=istat)
  if(istat/=0) stop 'Allocation error (7)'
  do irank=0,nproc-1
    fgb2=fgb
    write(fgb2(lfgb-3:lfgb),'(i4.4)') irank
    !Reserve 8 bytes for all integers as well
    ihot_len=8*(7+((3+2*ntracers)*nvrt+1)*ner(irank))
    open(36,file=fgb2(1:lfgb)//'_hotstart',access='direct',recl=ihot_len,status='old')
    read(36,rec=1)icm,ised2d,ised3d,iha,time,it,ifile,(iwild(ielg(irank,i)),(swild(:,j,ielg(irank,i)),j=1,nvrt),i=1,ner(irank))
    close(36)

    !Debug
    !write(98,*)icm,ised2d,iha,time,it,ifile
  enddo !irank

  write(37) time,it,ifile
  do i=1,ne_global
    write(37) i,iwild(i),(swild(:,j,i),j=1,nvrt) !(we(j,i),tsel(1:2,j,i),(trel0(l,j,i),trel(l,j,i),l=1,ntracers),j=1,nvrt)
  enddo !i

  !Debug
  write(98,*)'Tracer:'
  do i=1,ne_global
    do l=1,ntracers
      write(98,*)i,2*l-1,swild(3+2*l-1,:,i)
      write(98,*)i,2*l,swild(3+2*l,:,i)
    enddo !l
  enddo !i

  deallocate(iwild,swild)

  !Read/write sides
  allocate(iwild(ns_global),swild(4,nvrt,ns_global),wild2((3+2*ntracers)*nvrt),stat=istat)
  if(istat/=0) stop 'Allocation error (7)'
  do irank=0,nproc-1
    fgb2=fgb
    write(fgb2(lfgb-3:lfgb),'(i4.4)') irank
    ihot_len=8*(7+((3+2*ntracers)*nvrt+1)*ner(irank)+(4*nvrt+1)*nsr(irank))
    open(36,file=fgb2(1:lfgb)//'_hotstart',access='direct',recl=ihot_len,status='old')
    read(36,rec=1)icm,ised2d,ised3d,iha,time,it,ifile,(itmp,wild2(:),i=1,ner(irank)), &
     &(iwild(islg(irank,i)),(swild(:,j,islg(irank,i)),j=1,nvrt),i=1,nsr(irank))
    close(36)
  enddo !irank

  do i=1,ns_global
    write(37) i,iwild(i),(swild(:,j,i),j=1,nvrt) !idry_s(i),(su2(j,i),sv2(j,i),tsd(j,i),ssd(j,i),j=1,nvrt)
  enddo !i
  deallocate(iwild,swild)

  !Read/write nodes
  allocate(iwild(np_global),swild(11,nvrt,np_global),wild3(4*nvrt),wild4(np_global),stat=istat)
  if(istat/=0) stop 'Allocation error (7)'
  do irank=0,nproc-1
    fgb2=fgb
    write(fgb2(lfgb-3:lfgb),'(i4.4)') irank
    ihot_len=8*(7+((3+2*ntracers)*nvrt+1)*ner(irank)+ &
              &(4*nvrt+1)*nsr(irank)+(2+11*nvrt)*npr(irank))
    !Save the record # at the end of Hydro part of hotstart.in (assuming 8-byte)
    ihotstp(irank)=ihot_len/8
 
    open(36,file=fgb2(1:lfgb)//'_hotstart',access='direct',recl=ihot_len,status='old')
    read(36,rec=1)icm,ised2d,ised3d,iha,time,it,ifile,(itmp,wild2(:),i=1,ner(irank)), &
     &(itmp,wild3(:),i=1,nsr(irank)), &
     &(wild4(iplg(irank,i)),iwild(iplg(irank,i)), &
     &(swild(:,j,iplg(irank,i)),j=1,nvrt),i=1,npr(irank))
    close(36)
  enddo !irank

  do i=1,np_global
    write(37) i,wild4(i),iwild(i),(swild(:,j,i),j=1,nvrt) 
  enddo !i
  !Debug
  !do i=1,np_global
  !  write(96,*)i,real(wild4(i))
  !enddo !i
  !write(98,*)'Last record before modules:',swild(11,:,np_global)

  deallocate(iwild,swild,wild2,wild3,wild4)

  !Modules data
  if(icm==1) then
    allocate(swild(1,40,ne_global),stat=istat)
    if(istat/=0) stop 'Allocation error (8)'
    do irank=0,nproc-1
      fgb2=fgb
      write(fgb2(lfgb-3:lfgb),'(i4.4)') irank
      !Change to 8-byte recl!!
      open(36,file=fgb2(1:lfgb)//'_hotstart',access='direct',recl=8,status='old')

      !Debug
!      read(36,rec=ihotstp(irank))tmp
!      print*, 'B4 ICM:',tmp

      do i=1,ner(irank)
        do j=1,40
          read(36,rec=ihotstp(irank)+j)swild(1,j,ielg(irank,i))
!          print*, 'i,j:',i,j,swild(1,j,ielg(irank,i))
        enddo !j
        ihotstp(irank)=ihotstp(irank)+40
      enddo !i
      close(36)
    enddo !irank

    do i=1,ne_global
      write(37) i,swild(1,1:40,i)
      !Debug
      write(97,*)'ICM_SED:',i,swild(1,1:40,i)
    enddo !i
    deallocate(swild)
  endif !icm==1

  if(ised2d==1) then !SED2D
    allocate(wild2(np_global),stat=istat)
    if(istat/=0) stop 'Allocation error (9)'
    do irank=0,nproc-1
      fgb2=fgb
      write(fgb2(lfgb-3:lfgb),'(i4.4)') irank
      !Change to 8-byte recl!!
      open(36,file=fgb2(1:lfgb)//'_hotstart',access='direct',recl=8,status='old')

      !Debug
!      read(36,rec=ihotstp(irank))tmp
!      print*, 'B4 ICM:',tmp

      do i=1,npr(irank)
        read(36,rec=ihotstp(irank)+1)wild2(iplg(irank,i))
!       print*, 'i,j:',i,j,wild2(iplg(irank,i))
        ihotstp(irank)=ihotstp(irank)+1
      enddo !i
      close(36)
    enddo !irank

    do i=1,np_global
      write(37) i,wild2(i)
      !Debug
      write(97,*)'SED:',i,wild2(i)
    enddo !i
    deallocate(wild2)
  endif !ised2d==1

  if(ised3d==1) then !SED3D
    allocate(wild2(np_global),stat=istat)
    if(istat/=0) stop 'Allocation error (9)'
    do irank=0,nproc-1
      fgb2=fgb
      write(fgb2(lfgb-3:lfgb),'(i4.4)') irank
      !Change to 8-byte recl!!
      open(36,file=fgb2(1:lfgb)//'_hotstart',access='direct',recl=8,status='old')

      !Debug
!      read(36,rec=ihotstp(irank))tmp
!      print*, 'B4 ICM:',tmp

      read(36,rec=ihotstp(irank)+1)MBEDP
      read(36,rec=ihotstp(irank)+2)Nbed
      ihotstp(irank)=ihotstp(irank)+2
!      print*, 'SED3D dim=',MBEDP,Nbed,ntracers,ne_global
      if(irank==0) then !allocate
        allocate(bed(Nbed,ne_global,MBEDP),rough_p(np_global), &
     &bed_frac(Nbed,ne_global,ntracers),stat=istat)
        if(istat/=0) stop 'Allocation error (12)'
        !Test memory
        bed(Nbed,ne_global,MBEDP)=0 
        bed_frac(Nbed,ne_global,ntracers)=0
      endif !irank  

      do i=1,npr(irank)
        read(36,rec=ihotstp(irank)+1)wild2(iplg(irank,i))
        read(36,rec=ihotstp(irank)+2)rough_p(iplg(irank,i))
!       print*, 'i,j:',i,j,wild2(iplg(irank,i))
        ihotstp(irank)=ihotstp(irank)+2
      enddo !i

      do i=1,MBEDP
        do j=1,ner(irank)
          do k=1,Nbed
            read(36,rec=ihotstp(irank)+1)bed(k,ielg(irank,j),i)
            ihotstp(irank)=ihotstp(irank)+1
          enddo !k
        enddo !j
      enddo !i

      do i=1,ntracers
        do k=1,ner(irank)
          do m=1,Nbed
            read(36,rec=ihotstp(irank)+1)bed_frac(m,ielg(irank,k),i)
            ihotstp(irank)=ihotstp(irank)+1
          enddo !m
        enddo !k
      enddo !i

      close(36)
    enddo !irank

    !Write
    do i=1,np_global
      write(37) i,wild2(i) !depth
      write(37) i,rough_p(i)
      !Debug
      write(97,*)'SED3D:',i,wild2(i)
    enddo !i
    deallocate(wild2)

    do i=1,ne_global    
      !bed(Nbed,nea,3), bed_frac(Nbed,nea,ntracers)
      write(37) i,bed_frac(:,i,:),bed(:,i,:)
      write(97,*)'SED3D bed:',i,bed_frac(:,i,:)
    enddo !i
    deallocate(bed,bed_frac,rough_p)
  endif !ised3d==1

  if(iha==1) then
    !HA part
  endif !iha==1

  close(37)

  print*, 'time,it,ifile:',time,it,ifile

  stop
end subroutine combine_hotstart5

