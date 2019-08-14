!===============================================================================
! Read in binary outputs (rank-specific) from parallel code and combine them into
! one global output in v5.0 format. Works for partial outputs.
! Gobal-local mappings are read in from separate files.

! Inputs:
!        rank-specific binary files; 
!        local_to_global_*; 
!        combine_output.in (1st line: elev.61 etc; 2nd line: start and end file indices); 
! Output: combined binary file.
!
!  ifort -Bstatic -O3 -assume byterecl -o combine_output4 combine_output4.f90
!
!  PGI compiler:
!  pgf90 -O2 -mcmodel=medium  -Bstatic -o combine_output4 combine_output4.f90

!  History: (1) rank-specific binary files stripped of headers (v2.0d and up) for more efficient combining.
!           (2) Reduced the dimensions of iplg, ielg, kbp01
!===============================================================================

program read_iwrite1
!-------------------------------------------------------------------------------
  implicit real(4)(a-h,o-z),integer(i-n)
  parameter(nbyte=4)
  character*30 file63
  character*12 it_char
  character*48 start_time,version,variable_nm,variable_dim
  character*48 data_format
  character(72) :: fgb,fgb2,fdb  ! Processor specific global output file name
  integer :: lfgb,lfdb       ! Length of processor specific global output file name
  character(len= 4) :: a_4
  allocatable ne(:),np(:),nsproc(:),ihot_len(:)
  allocatable ztot(:),sigma(:),outb(:,:,:),eta2(:)
  allocatable i34(:),nm(:,:),nm2(:,:),js(:,:),xctr(:),yctr(:),dpe(:)
  allocatable x(:),y(:),dp(:),kbp00(:),iplg(:,:),ielg(:,:),kbp01(:,:)
  allocatable xcj(:),ycj(:),dps(:)
!-------------------------------------------------------------------------------
      
!-------------------------------------------------------------------------------
! Aquire user inputs
!-------------------------------------------------------------------------------

  open(10,file='combine_output.in',status='old')
  read(10,'(a30)') file63
  read(10,*) ibgn,iend
  close(10)

! Read local_to_global_0000 for global info
  open(10,file='local_to_global_0000',status='old')
  read(10,*)ne_global,np_global,nvrt,nproc
  close(10)

  allocate(x(np_global),y(np_global),dp(np_global),kbp00(np_global), &
           np(0:nproc-1),ne(0:nproc-1),nm(ne_global,3),nm2(ne_global,3),eta2(np_global), &
           ztot(nvrt),sigma(nvrt),outb(np_global,nvrt,2),ihot_len(0:nproc-1),stat=istat)
  if(istat/=0) stop 'Allocation error: x,y'

!-------------------------------------------------------------------------------
! Read rank-specific local_to_global*
!-------------------------------------------------------------------------------

  ! Compute ivs and i23d
  lfile63=len_trim(file63)
  if(file63((lfile63-1):lfile63).eq.'61'.or.file63((lfile63-1):lfile63).eq.'63') then
    ivs=1
  else if(file63((lfile63-1):lfile63).eq.'62'.or.file63((lfile63-1):lfile63).eq.'64') then
    ivs=2
  else
    print*, 'Unknown file type:',file63
  endif
  if(file63((lfile63-1):lfile63).eq.'61'.or.file63((lfile63-1):lfile63).eq.'62') then
    i23d=2
  else
    i23d=3
  endif
   
  ! Read in local-global mappings from all ranks
  fdb='local_to_global_0000'
  lfdb=len_trim(fdb)

  !Find max. for dimensioning
  np_max=0; ne_max=0
  do irank=0,nproc-1
    write(fdb(lfdb-3:lfdb),'(i4.4)') irank
    open(10,file=fdb,status='old')
    read(10,*) !global info
    read(10,*) !info
    read(10,*)ne(irank)
    do i=1,ne(irank)
      read(10,*)!j,ielg(irank,i)
    enddo !i
    read(10,*)np(irank)
    close(10)
    np_max=max(np_max,np(irank))
    ne_max=max(ne_max,ne(irank))
  enddo !irank

  allocate(iplg(0:nproc-1,np_max),kbp01(0:nproc-1,np_max), &
     &ielg(0:nproc-1,ne_max),stat=istat)
  if(istat/=0) stop 'Allocation error (2)'

  !Re-read
  do irank=0,nproc-1
    write(fdb(lfdb-3:lfdb),'(i4.4)') irank
    open(10,file=fdb,status='old')

    read(10,*) !global info

    read(10,*) !info
    read(10,*)ne(irank)
    do i=1,ne(irank)
      read(10,*)j,ielg(irank,i)
    enddo !i
    read(10,*)np(irank)
    do i=1,np(irank)
      read(10,*)j,iplg(irank,i)
    enddo
    read(10,*)itmp !sides
    do i=1,itmp
      read(10,*)
    enddo
!    print*, 'Mapping:',irank,ielg(irank,ne(irank))

    read(10,*) !'Header:'
    read(10,'(a)')data_format,version,start_time
    read(10,*)nrec,dtout,nspool,nvrt,kz,h0,h_s,h_c,theta_b,theta_f
    read(10,*)(ztot(k),k=1,kz-1),(sigma(k),k=1,nvrt-kz+1)
    read(10,*)np(irank),ne(irank),(x(iplg(irank,m)),y(iplg(irank,m)),dp(iplg(irank,m)),kbp01(irank,m),m=1,np(irank)), &
    &         (ntmp,(nm2(m,mm),mm=1,3),m=1,ne(irank))

    close(10)

!   Compute kbp00 (to avoid mismatch of indices)
    do m=1,np(irank)
      ipgb=iplg(irank,m)
      kbp00(ipgb)=kbp01(irank,m)
    enddo !m
 
!   Reconstruct connectivity table
    do m=1,ne(irank)
      iegb=ielg(irank,m)
      if(iegb>ne_global) stop 'Overflow!'
      do mm=1,3
        itmp=nm2(m,mm)
        if(itmp>np(irank).or.itmp<=0) then
          write(*,*)'Overflow:',m,mm,itmp
          stop
        endif
        nm(iegb,mm)=iplg(irank,itmp)
      enddo !mm
    enddo !m
  enddo !irank=0,nproc-1

! Compute record length for each rank-specific binary output per time step
  do irank=0,nproc-1
    ihot_len(irank)=nbyte*(2+np(irank))
    if(i23d==2) then
      ihot_len(irank)=ihot_len(irank)+nbyte*ivs*np(irank)
    else
      do i=1,np(irank)
        do k=max0(1,kbp01(irank,i)),nvrt
          do m=1,ivs
            ihot_len(irank)=ihot_len(irank)+nbyte
          enddo !m
        enddo !k
      enddo !i
    endif
  enddo !irank

!-------------------------------------------------------------------------------
! Time iteration -- select "node" data
!-------------------------------------------------------------------------------

! Loop over input files
  do iinput=ibgn,iend
    write(it_char,'(i12)')iinput
    it_char=adjustl(it_char)  !place blanks at end
    it_len=len_trim(it_char)  !length without trailing blanks
    fgb=it_char(1:it_len)//'_0000'; lfgb=len_trim(fgb);

    ! Read actual number of spools in this file
!    open(63,file=fgb(1:lfgb)//'_'//file63,access='direct',recl=nbyte,status='old')
!    read(63,rec=irec_nrec) nrec
!    close(63)

    !Write header to output files
    open(65,file=it_char(1:it_len)//'_'//file63,status='replace')
    data_format='DataFormat v5.0'
    variable_nm=file63 !not important
    variable_dim=file63

    write(65,'(a48)',advance="no") data_format
    write(65,'(a48)',advance="no") version
    write(65,'(a48)',advance="no") start_time
    write(65,'(a48)',advance="no") variable_nm
    write(65,'(a48)',advance="no") variable_dim

    a_4 = transfer(source=nrec,mold=a_4)
    write(65,"(a4)",advance="no") a_4
    a_4 = transfer(source=dtout,mold=a_4)
    write(65,"(a4)",advance="no") a_4
    a_4 = transfer(source=nspool,mold=a_4)
    write(65,"(a4)",advance="no") a_4
    a_4 = transfer(source=ivs,mold=a_4)
    write(65,"(a4)",advance="no") a_4
    a_4 = transfer(source=i23d,mold=a_4)
    write(65,"(a4)",advance="no") a_4

    !Vertical grid
    a_4 = transfer(source=nvrt,mold=a_4)
    write(65,'(a4)',advance="no") a_4
    a_4 = transfer(source=kz,mold=a_4)
    write(65,'(a4)',advance="no") a_4
    a_4 = transfer(source=h0,mold=a_4)
    write(65,'(a4)',advance="no") a_4
    a_4 = transfer(source=h_s,mold=a_4)
    write(65,'(a4)',advance="no") a_4
    a_4 = transfer(source=h_c,mold=a_4)
    write(65,'(a4)',advance="no") a_4
    a_4 = transfer(source=theta_b,mold=a_4)
    write(65,'(a4)',advance="no") a_4
    a_4 = transfer(source=theta_f,mold=a_4)
    write(65,'(a4)',advance="no") a_4

    do k=1,kz-1
      a_4 = transfer(source=ztot(k),mold=a_4)
      write(65,'(a4)',advance="no") a_4
    enddo
    do k=kz,nvrt
      kin=k-kz+1
      a_4 = transfer(source=sigma(kin),mold=a_4)
      write(65,'(a4)',advance="no") a_4
    enddo !k

    !Horizontal grid
    a_4 = transfer(source=np_global,mold=a_4)
    write(65,'(a4)',advance="no") a_4
    a_4 = transfer(source=ne_global,mold=a_4)
    write(65,'(a4)',advance="no") a_4

    do m=1,np_global
      a_4 = transfer(source=x(m),mold=a_4)
      write(65,'(a4)',advance="no") a_4
      a_4 = transfer(source=y(m),mold=a_4)
      write(65,'(a4)',advance="no") a_4
      a_4 = transfer(source=dp(m),mold=a_4)
      write(65,'(a4)',advance="no") a_4
      a_4 = transfer(source=kbp00(m),mold=a_4)
      write(65,'(a4)',advance="no") a_4
    enddo !m=1,np
    do m=1,ne_global
      a_4 = transfer(source=3,mold=a_4)
      write(65,'(a4)',advance="no") a_4
      do mm=1,3
        a_4 = transfer(source=nm(m,mm),mold=a_4)
        write(65,'(a4)',advance="no") a_4
      enddo !mm
    enddo !m

    !print*, 'Last element:',nm(ne_global,1:3)
    !end output header

    ! Loop over output spools in file
    do ispool=1,nrec

      !Gather all ranks
      do irank=0,nproc-1
        !Open input file
        fgb2=fgb
        write(fgb2(lfgb-3:lfgb),'(i4.4)') irank

        open(63,file=fgb2(1:lfgb)//'_'//file63,access='direct',recl=ihot_len(irank),status='old')
 
        if(i23d==2) then
          read(63,rec=ispool)time,it,(eta2(iplg(irank,i)),i=1,np(irank)),((outb(iplg(irank,i),1,m),m=1,ivs),i=1,np(irank))
        else !3D
          read(63,rec=ispool)time,it,(eta2(iplg(irank,i)),i=1,np(irank)), &
     &                       (((outb(iplg(irank,i),k,m),m=1,ivs),k=max0(1,kbp01(irank,i)),nvrt),i=1,np(irank))
        endif
    
        ! Close input file
        close(63)
      enddo !irank

      !Output
      !print*, 'it=',it,' time=',time/86400.0

      a_4 = transfer(source=time,mold=a_4)
      write(65,"(a4)",advance="no") a_4
      a_4 = transfer(source=it,mold=a_4)
      write(65,"(a4)",advance="no") a_4

      do i=1,np_global
        a_4 = transfer(source=eta2(i),mold=a_4)
        write(65,"(a4)",advance="no") a_4
      enddo !i
      
      do i=1,np_global
        if(i23d==2) then
          do m=1,ivs
            a_4 = transfer(source=outb(i,1,m),mold=a_4)
            write(65,"(a4)",advance="no") a_4
          enddo !m
        else !i23d=3 
          do k=max0(1,kbp00(i)),nvrt
            do m=1,ivs
              a_4 = transfer(source=outb(i,k,m),mold=a_4)
              write(65,"(a4)",advance="no") a_4
            enddo !m
          enddo !k
        endif !i23d
      enddo !i

    enddo !ispool=1,nrec
  enddo !iinput=1,ninput_files

! Close output file
  close(65)

  stop
end program read_iwrite1
