!===============================================================================
! Read in *.gr3-like (rank-specific) outputs from SELFE and combine them into one global output
! e.g. maxelev.gr3; may have multiple scalar fields

! Inputs:
!        (0) screen: filenm - name of file (e.g. maxelev);
!        (1) hgrid.gr3;
!        (2) outputs/<filenm>_0*
! Output: <filenm>.gr3
!
!  Compile on amb64xx:
!  ifort -Bstatic -O3 -assume byterecl -o combine_gr3 combine_gr3.f90
!  PGI compiler:
!  pgf90 -O2 -mcmodel=medium  -Bstatic -o combine_gr3 combine_gr3.f90
!===============================================================================

program combine_gr3
!-------------------------------------------------------------------------------

  implicit real(8)(a-h,o-z),integer(i-n)
  character(36) :: fdb,filenm 
  integer :: lfdb,lfilenm
  allocatable x(:),y(:),elevmax(:,:)
      
!-------------------------------------------------------------------------------
! inputs
!-------------------------------------------------------------------------------

  print*, 'Input file name (e.g.: maxelev):'
  read*, filenm
  filenm=adjustl(filenm); lfilenm=len_trim(filenm)
  print*, 'Input # of scalar fields:'
  read*, nscal
  if(nscal<=0) stop 'Wrong nscal'

  open(14,file='hgrid.gr3',status='old')
  read(14,*); read(14,*)ne,np
  allocate(x(np),y(np),elevmax(nscal,np),stat=istat)
  if(istat/=0) stop 'Allocation error: x,y'

  open(10,file='outputs/'//filenm(1:lfilenm)//'_0000',status='old')
  read(10,*)icount,nproc
  close(10)

!-------------------------------------------------------------------------------
! Combine
!-------------------------------------------------------------------------------
  fdb=filenm(1:lfilenm)//'_0000'
  lfdb=len_trim(fdb)

  do irank=0,nproc-1
    write(fdb(lfdb-3:lfdb),'(i4.4)') irank
    open(10,file='outputs/'//fdb,status='old')
    read(10,*)icount !,nproc
    do i=1,icount
      read(10,*)nd,xtmp,ytmp,elevmax(:,nd)
    enddo !i
  enddo !irank

  open(13,file=filenm(1:lfilenm)//'.gr3',status='replace')
  write(13,*); write(13,*)ne,np
  do i=1,np
    read(14,*)j,xtmp,ytmp
    write(13,'(i10,100(1x,e22.11))')i,xtmp,ytmp,elevmax(:,i)
  enddo !i
  do i=1,ne
    read(14,*)j,k,n1,n2,n3
    write(13,*)j,k,n1,n2,n3
  enddo !i

end program combine_gr3
