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
!  ifort -Bstatic -O3 -assume byterecl -o combine_hotstart3 combine_hotstart3.f90
!  PGI compiler:
!  pgf90 -O2 -mcmodel=medium  -Bstatic -o combine_hotstart3 combine_hotstart3.f90

! Revisions: non-hydrostatic pressure added (v3.0a up).
!===============================================================================

program combine_hotstart1
!-------------------------------------------------------------------------------
  implicit real(8)(a-h,o-z),integer(i-n)
  parameter(nbyte=4)
  character(12) :: it_char
  character(72) :: fgb,fgb2,fdb  ! Processor specific global output file name
  integer :: lfgb,lfdb       ! Length of processor specific global output file name
  allocatable ner(:),npr(:),nsr(:)
  allocatable ielg(:,:),iplg(:,:),islg(:,:)
  allocatable idry_e(:),we(:,:),tsel(:,:,:),idry_s(:),su2(:,:),sv2(:,:)
  allocatable tsd(:,:),ssd(:,:),idry(:),eta2(:),tnd(:,:),snd(:,:)
  allocatable tem0(:,:),sal0(:,:),q2(:,:),xl(:,:),dfv(:,:),dfh(:,:)
  allocatable dfq1(:,:),dfq2(:,:),qnon(:,:),trel0(:,:,:),trel(:,:,:)
  allocatable iwild(:),swild(:,:,:),wild2(:),wild3(:),wild4(:)
!-------------------------------------------------------------------------------
      
!-------------------------------------------------------------------------------
! Aquire user inputs
!-------------------------------------------------------------------------------

!  open(10,file='combine_hotstart1.in',status='old')
  write(*,*) 'Input nproc, ntracers:'
  read(*,*) nproc,ntracers
  write(*,*) 'Input iteration # of the hotstart file (before _00*_hotstart) :'
  read(*,'(a)')it_char 
  write(*,*) 'Do you want to use memorys-saving mode (1 - yes):'
  read(*,*)imem
  if(imem/=0.and.imem/=1) stop 'Unknown imem'

  allocate(ner(0:nproc-1),npr(0:nproc-1),nsr(0:nproc-1),stat=istat) 
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

  if(imem==0) then
  endif !imem==0

!-------------------------------------------------------------------------------
! Read hotstart files
!-------------------------------------------------------------------------------

  ! Open file
  it_char=adjustl(it_char)  !place blanks at end
  it_len=len_trim(it_char)  !length without trailing blanks
  fgb=it_char(1:it_len)//'_0000'; lfgb=len_trim(fgb);

  !Gather all ranks
  ! Output
  open(37,file='hotstart.in',form='unformatted',status='replace')

  if(imem==0) then !non memory saving mode
!============================================================================
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

    do irank=0,nproc-1
      !Open input file
      fgb2=fgb
      write(fgb2(lfgb-3:lfgb),'(i4.4)') irank
      ihot_len=nbyte*(4+((6+4*ntracers)*nvrt+1)*ner(irank)+ &
                &(8*nvrt+1)*nsr(irank)+(3+22*nvrt)*npr(irank))
      open(36,file=fgb2(1:lfgb)//'_hotstart',access='direct',recl=ihot_len,status='old')

      read(36,rec=1)time,it,ifile,(idry_e(ielg(irank,i)),(we(j,ielg(irank,i)),tsel(1:2,j,ielg(irank,i)), &
     &(trel0(l,j,ielg(irank,i)),trel(l,j,ielg(irank,i)),l=1,ntracers),j=1,nvrt),i=1,ner(irank)), &
     &(idry_s(islg(irank,i)),(su2(j,islg(irank,i)),sv2(j,islg(irank,i)), &
     &tsd(j,islg(irank,i)),ssd(j,islg(irank,i)),j=1,nvrt),i=1,nsr(irank)), &
     &(eta2(iplg(irank,i)),idry(iplg(irank,i)),(tnd(j,iplg(irank,i)),snd(j,iplg(irank,i)), &
     &tem0(j,iplg(irank,i)),sal0(j,iplg(irank,i)),q2(iplg(irank,i),j),xl(iplg(irank,i),j), &
     &dfv(iplg(irank,i),j),dfh(iplg(irank,i),j),dfq1(iplg(irank,i),j), &
     &dfq2(iplg(irank,i),j),qnon(j,iplg(irank,i)),j=1,nvrt),i=1,npr(irank)) 
      close(36)
    enddo !irank

    write(37) time,it,ifile
    do i=1,ne_global
      write(37) i,idry_e(i),(we(j,i),tsel(1:2,j,i),(trel0(l,j,i),trel(l,j,i),l=1,ntracers),j=1,nvrt)
    enddo !i
    do i=1,ns_global
      write(37) i,idry_s(i),(su2(j,i),sv2(j,i),tsd(j,i),ssd(j,i),j=1,nvrt)
    enddo !i
    do i=1,np_global
      write(37) i,eta2(i),idry(i),(tnd(j,i),snd(j,i),tem0(j,i),sal0(j,i),q2(i,j),xl(i,j), &
               dfv(i,j),dfh(i,j),dfq1(i,j),dfq2(i,j),qnon(j,i),j=1,nvrt)
    enddo !i
!============================================================================
  else !memory saving mode
!    allocate(idry_e(mxner),we(nvrt,mxner),tsel(2,nvrt,mxner), &
!             idry_s(mxnsr),su2(nvrt,mxnsr),sv2(nvrt,mxnsr), &
!             tsd(nvrt,mxnsr),ssd(nvrt,mxnsr), &
!             idry(mxnpr),eta2(mxnpr),tnd(nvrt,mxnpr),snd(nvrt,mxnpr), &
!             tem0(nvrt,mxnpr),sal0(nvrt,mxnpr),q2(mxnpr,nvrt), &
!             xl(mxnpr,nvrt),dfv(mxnpr,nvrt),dfh(mxnpr,nvrt), &
!             dfq1(mxnpr,nvrt),dfq2(mxnpr,nvrt),qnon(nvrt,mxnpr),stat=istat)
!    if(istat/=0) stop 'Allocation error (5)'
!    if(ntracers==0) then
!      allocate(trel0(1,1,1),trel(1,1,1),stat=istat)
!    else
!      allocate(trel0(ntracers,nvrt,mxner),trel(ntracers,nvrt,mxner),stat=istat)
!    endif
!    if(istat/=0) stop 'Allocation error (6)'

    !Read/write elements first
    allocate(iwild(ne_global),swild(3+2*ntracers,nvrt,ne_global),stat=istat)
    if(istat/=0) stop 'Allocation error (7)'
    do irank=0,nproc-1
      fgb2=fgb
      write(fgb2(lfgb-3:lfgb),'(i4.4)') irank
      ihot_len=nbyte*(4+((6+4*ntracers)*nvrt+1)*ner(irank)) !+ &
!                &(8*nvrt+1)*nsr(irank)+(3+22*nvrt)*npr(irank))
      open(36,file=fgb2(1:lfgb)//'_hotstart',access='direct',recl=ihot_len,status='old')
      read(36,rec=1)time,it,ifile,(iwild(ielg(irank,i)),(swild(:,j,ielg(irank,i)),j=1,nvrt),i=1,ner(irank))
      close(36)
    enddo !irank

    write(37) time,it,ifile
    do i=1,ne_global
      write(37) i,iwild(i),(swild(:,j,i),j=1,nvrt) !(we(j,i),tsel(1:2,j,i),(trel0(l,j,i),trel(l,j,i),l=1,ntracers),j=1,nvrt)
    enddo !i
    deallocate(iwild,swild)

    !Read/write sides
    allocate(iwild(ns_global),swild(4,nvrt,ns_global),wild2((3+2*ntracers)*nvrt),stat=istat)
    if(istat/=0) stop 'Allocation error (7)'
    do irank=0,nproc-1
      fgb2=fgb
      write(fgb2(lfgb-3:lfgb),'(i4.4)') irank
      ihot_len=nbyte*(4+((6+4*ntracers)*nvrt+1)*ner(irank)+ &
               &(8*nvrt+1)*nsr(irank)) !+(3+22*nvrt)*npr(irank))
      open(36,file=fgb2(1:lfgb)//'_hotstart',access='direct',recl=ihot_len,status='old')
      read(36,rec=1)time,it,ifile,(itmp,wild2(:),i=1,ner(irank)), &
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
      ihot_len=nbyte*(4+((6+4*ntracers)*nvrt+1)*ner(irank)+ &
                &(8*nvrt+1)*nsr(irank)+(3+22*nvrt)*npr(irank))
      open(36,file=fgb2(1:lfgb)//'_hotstart',access='direct',recl=ihot_len,status='old')
      read(36,rec=1)time,it,ifile,(itmp,wild2(:),i=1,ner(irank)), &
     &(itmp,wild3(:),i=1,nsr(irank)), &
     &(wild4(iplg(irank,i)),iwild(iplg(irank,i)), &
     &(swild(:,j,iplg(irank,i)),j=1,nvrt),i=1,npr(irank))
      close(36)
    enddo !irank

    do i=1,np_global
      write(37) i,wild4(i),iwild(i),(swild(:,j,i),j=1,nvrt) 
!eta2(i),idry(i),(tnd(j,i),snd(j,i),tem0(j,i),sal0(j,i),q2(i,j),xl(i,j), &
!               dfv(i,j),dfh(i,j),dfq1(i,j),dfq2(i,j),qnon(j,i),j=1,nvrt)
    enddo !i
!============================================================================
  endif !imem

  close(37)

  print*, 'time,it,ifile:',time,it,ifile

  stop
end program combine_hotstart1
