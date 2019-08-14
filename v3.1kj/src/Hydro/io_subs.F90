!===============================================================================
!===============================================================================
! ELCIRC FILE I/O SUBROUTINES
!
! subroutine write_header1 (not used due to combine_output*.F90)
! subroutine write_header0 (not used)
! subroutine write_obe
! subroutine report_timers
! subroutine get_param

!===============================================================================
!===============================================================================

      subroutine write_header1
!-------------------------------------------------------------------------------
! Write header data to global output files (iwrite=1). Separate file from each processor.
!-------------------------------------------------------------------------------
#ifdef USE_MPIMODULE
      use mpi
#endif
      use elfe_glbl
      use elfe_msgp
      implicit none
#ifndef USE_MPIMODULE
      include 'mpif.h'
#endif
!      integer,intent(in) :: iths,it !nt
      integer :: i,k,m,mm,ivs,i23d,ipgb,iegb,irecm,kin
      character(72) :: fgb  ! Processor specific global output file name
      integer :: lfgb       ! Length of processor specific global output file name
      character(len=48) :: variable_out
!      integer :: nrec
!      real(rkind) :: cwtmp
!      logical :: cwtime
!-------------------------------------------------------------------------------

! Actual number of spools to occur
!      nrec=min(nt-it,ihfskip)/nspool

! Global output files
      do i=1,noutput

        irec(i)=0 !reset record # for binary
        ichan(i)=100+i !output channel #
        if(i>=13.and.i<=15.or.i==26) then
          ivs=2
        else
          ivs=1
        endif
        if(i<=15) then
          i23d=2 !2 or 3D
        else
          i23d=3
        endif

        if(iof(i)==1) then

          !Open new processor specific file
          fgb=ifile_char(1:ifile_len)//'_0000'; lfgb=len_trim(fgb);
          write(fgb(lfgb-3:lfgb),'(i4.4)') myrank
          open(ichan(i),file='outputs/'//(fgb(1:lfgb)//'_'//outfile(i)),status='replace')

          write(ichan(i),'(a48)',advance="no") data_format
          write(ichan(i),'(a48)',advance="no") version
          write(ichan(i),'(a48)',advance="no") start_time
          write(ichan(i),'(a48)',advance="no") variable_nm(i)
          write(ichan(i),'(a48)',advance="no") variable_dim(i)

          irec(i)=irec(i)+48/nbyte

          a_4 = transfer(source=nrec,mold=a_4)
          write(ichan(i),"(a4)",advance="no") a_4
          a_4 = transfer(source=real(dt*nspool),mold=a_4)
          write(ichan(i),"(a4)",advance="no") a_4
          a_4 = transfer(source=nspool,mold=a_4)
          write(ichan(i),"(a4)",advance="no") a_4
          a_4 = transfer(source=ivs,mold=a_4)
          write(ichan(i),"(a4)",advance="no") a_4
          a_4 = transfer(source=i23d,mold=a_4)
          write(ichan(i),"(a4)",advance="no") a_4

          irec(i)=irec(i)+5

          !Vertical grid size and data
          a_4 = transfer(source=nvrt,mold=a_4)
          write(ichan(i),'(a4)',advance="no") a_4
          a_4 = transfer(source=kz,mold=a_4)
          write(ichan(i),'(a4)',advance="no") a_4
          a_4 = transfer(source=real(h0),mold=a_4)
          write(ichan(i),'(a4)',advance="no") a_4
          a_4 = transfer(source=real(h_s),mold=a_4)
          write(ichan(i),'(a4)',advance="no") a_4
          a_4 = transfer(source=real(h_c),mold=a_4)
          write(ichan(i),'(a4)',advance="no") a_4
          a_4 = transfer(source=real(theta_b),mold=a_4)
          write(ichan(i),'(a4)',advance="no") a_4
          a_4 = transfer(source=real(theta_f),mold=a_4)
          write(ichan(i),'(a4)',advance="no") a_4

          irec(i)=irec(i)+7
          do k=1,kz-1
            a_4 = transfer(source=real(ztot(k)),mold=a_4)
            write(ichan(i),'(a4)',advance="no") a_4
          enddo
          do k=kz,nvrt
            kin=k-kz+1
            a_4 = transfer(source=real(sigma(kin)),mold=a_4)
            write(ichan(i),'(a4)',advance="no") a_4
          enddo !k
          irec(i)=irec(i)+nvrt

          !Horizontal grid size (resident)
          a_4 = transfer(source=np,mold=a_4)
          write(ichan(i),'(a4)',advance="no") a_4
          a_4 = transfer(source=ne,mold=a_4)
          write(ichan(i),'(a4)',advance="no") a_4

          irec(i)=irec(i)+2

          !Write resident node data
          do m=1,np
!            ipgb=iplg(m) !global node index
!            write(ichan(i),rec=irec(i)+1)ipgb
            a_4 = transfer(source=real(xnd(m)),mold=a_4)
            write(ichan(i),'(a4)',advance="no") a_4
            a_4 = transfer(source=real(ynd(m)),mold=a_4)
            write(ichan(i),'(a4)',advance="no") a_4
            a_4 = transfer(source=real(dp00(m)),mold=a_4)
            write(ichan(i),'(a4)',advance="no") a_4
            a_4 = transfer(source=kbp00(m),mold=a_4)
            write(ichan(i),'(a4)',advance="no") a_4
            irec(i)=irec(i)+4
          enddo !m=1,np

          !Write resident element connectivity data 
          do m=1,ne
!            iegb=ielg(m) !global element index
            a_4 = transfer(source=3,mold=a_4)
            write(ichan(i),'(a4)',advance="no") a_4
            irec(i)=irec(i)+1
            do mm=1,3
              a_4 = transfer(source=nm(m,mm),mold=a_4)
              write(ichan(i),'(a4)',advance="no") a_4
            enddo !mm
            irec(i)=irec(i)+3
          enddo !m=1,ne

          !Close file
!          close(ichan(i))

!         Estimate total # of records
          irecm=48/nbyte*5+5+7+nvrt+4*np+4*ne
          if(i23d==2) then !2D
            irecm=irecm+(2+np+np*ivs)*nrec
          else !3D
            irecm=irecm+(2+np+np*nvrt*ivs)*nrec
          endif

          if(irecm>mirec) then
            write(errmsg,*)'Output file too large',i,irecm
            call parallel_abort(errmsg)
          endif      
        endif !iof(i)=1

      enddo !i=1,noutput

! Test output
      igmp=0
      if(noutgm==1) then

        !All ranks initialize record counter
        igmp=(4+3+3)*8/nbyte

        fgb=ifile_char(1:ifile_len)//'_0000'; lfgb=len_trim(fgb);
        write(fgb(lfgb-3:lfgb),'(i4.4)') myrank
        open(100,file='outputs/'//(fgb(1:lfgb)//'_test.60'), &
             access='direct',recl=nbyte,status='replace')

! TODO: mapping from local to global
        write(100,rec=igmp+1) nrec
        write(100,rec=igmp+2) ns
        write(100,rec=igmp+3) real(dt*nspool)
        write(100,rec=igmp+4) nspool
        write(100,rec=igmp+5) 1
        igmp=igmp+5
        close(100)
      endif !noutgm.eq.1

      end subroutine write_header1

!===============================================================================
!===============================================================================

      subroutine write_header0(iths,it)
!-------------------------------------------------------------------------------
! Write header data to global output files (iwrite=0).
!-------------------------------------------------------------------------------
#ifdef USE_MPIMODULE
      use mpi
#endif
      use elfe_glbl
      use elfe_msgp
      implicit none
#ifndef USE_MPIMODULE
      include 'mpif.h'
#endif
      integer,intent(in) :: iths,it !nt
      integer :: i,k,m,mm,ivs,i23d,ipgb,iegb,ipgb_rec,iegb_rec,irecm,kin
!      integer :: nrec
      real(rkind) :: cwtmp
      character(len=48) :: variable_out
      logical :: cwtime
!-------------------------------------------------------------------------------

! Flag for comm timing
      cwtime=it/=iths

! Actual number of spools to occur
!      nrec=min(nt-it,ihfskip)/nspool

! Global output files
      do i=1,noutput

        irec(i)=0 !reset record # for binary
        ichan(i)=100+i !output channel #
        if(i>=13.and.i<=15.or.i==26) then
          ivs=2
        else
          ivs=1
        endif
        if(i<=15) then
          i23d=2 !2 or 3D
        else
          i23d=3
        endif

        if(iof(i)==1) then

          !Synchronization
#ifdef INCLUDE_TIMING
!          if(cwtime) cwtmp=mpi_wtime()
#endif
          call parallel_rrsync(1)
#ifdef INCLUDE_TIMING
!          if(cwtime) wtimer(14,2)=wtimer(14,2)+mpi_wtime()-cwtmp
#endif

          !Rank 0 open new file all others open existing file
          if(myrank==0) then
            open(ichan(i),file=ifile_char(1:ifile_len)//'_'//outfile(i), &
                 access='direct',recl=nbyte,status='replace')
          else
            open(ichan(i),file=ifile_char(1:ifile_len)//'_'//outfile(i), &
                 access='direct',recl=nbyte,status='old')
          endif

          if(myrank==0) then
            do m=1,48/nbyte
              write(ichan(i),rec=irec(i)+m) data_format(nbyte*(m-1)+1:nbyte*m)
            enddo
            do m=1,48/nbyte
              write(ichan(i),rec=irec(i)+48/nbyte+m) version(nbyte*(m-1)+1:nbyte*m)
            enddo
            do m=1,48/nbyte
              write(ichan(i),rec=irec(i)+2*48/nbyte+m) start_time(nbyte*(m-1)+1:nbyte*m)
            enddo
            variable_out=variable_nm(i)
            do m=1,48/nbyte
              write(ichan(i),rec=irec(i)+3*48/nbyte+m) variable_out(nbyte*(m-1)+1:nbyte*m)
            enddo
            variable_out=variable_dim(i)
            do m=1,48/nbyte
              write(ichan(i),rec=irec(i)+4*48/nbyte+m) variable_out(nbyte*(m-1)+1:nbyte*m)
            enddo
          endif !myrank==0
          irec(i)=irec(i)+5*48/nbyte
          if(myrank==0) then
            write(ichan(i),rec=irec(i)+1) nrec
            write(ichan(i),rec=irec(i)+2) real(dt*nspool)
            write(ichan(i),rec=irec(i)+3) nspool
            write(ichan(i),rec=irec(i)+4) ivs
            write(ichan(i),rec=irec(i)+5) i23d
          endif !myrank==0
          irec(i)=irec(i)+5

          !Vertical grid size and data
          if(myrank==0) then
            write(ichan(i),rec=irec(i)+1) nvrt
            write(ichan(i),rec=irec(i)+2) kz
            write(ichan(i),rec=irec(i)+3) real(h0)
            write(ichan(i),rec=irec(i)+4) real(h_s)
            write(ichan(i),rec=irec(i)+5) real(h_c)
            write(ichan(i),rec=irec(i)+6) real(theta_b)
            write(ichan(i),rec=irec(i)+7) real(theta_f)
          endif
          irec(i)=irec(i)+7
          if(myrank==0) then
            do k=1,kz-1
              write(ichan(i),rec=irec(i)+k) real(ztot(k))
            enddo
            do k=kz,nvrt
              kin=k-kz+1
              write(ichan(i),rec=irec(i)+k) real(sigma(kin))
            enddo !k
          endif !myrank==0
          irec(i)=irec(i)+nvrt

          !Horizontal grid size
          if(myrank==0) then
            write(ichan(i),rec=irec(i)+1) np_global
            write(ichan(i),rec=irec(i)+2) ne_global
          endif
          irec(i)=irec(i)+2

          !Write resident node data according to global node index
          do m=1,np
            ipgb=iplg(m) !global node index
            if(associated(ipgl(ipgb)%next)) then !interface node
              if(ipgl(ipgb)%next%rank<myrank) cycle !already written so skip
            endif
            ipgb_rec=irec(i)+4*(ipgb-1) !record before ipgb
            !Node coordinates and depth
            write(ichan(i),rec=ipgb_rec+1)real(xnd(m))
            write(ichan(i),rec=ipgb_rec+2)real(ynd(m))
            write(ichan(i),rec=ipgb_rec+3)real(dp(m))
            write(ichan(i),rec=ipgb_rec+4)kbp00(m)
          enddo !m=1,np
          irec(i)=irec(i)+4*np_global !last node record written

          !Write resident element connectivity data using global indices
          do m=1,ne
            iegb=ielg(m) !global element index
            iegb_rec=irec(i)+4*(iegb-1) !record before iegb
            write(ichan(i),rec=iegb_rec+1) 3
            do mm=1,3
              write(ichan(i),rec=iegb_rec+1+mm)iplg(nm(m,mm))
            enddo !mm
          enddo !m=1,ne
          irec(i)=irec(i)+4*ne_global !last element record written

          !Close file
          close(ichan(i))

          !Synchronization
#ifdef INCLUDE_TIMING
!          if(cwtime) cwtmp=mpi_wtime()
#endif
          call parallel_rrsync(-1)
#ifdef INCLUDE_TIMING
!          if(cwtime) wtimer(14,2)=wtimer(14,2)+mpi_wtime()-cwtmp
#endif

!         Estimate total # of records
          irecm=48/nbyte*5+5+7+nvrt+4*np_global+4*ne_global
          if(i23d==2) then !2D
            irecm=irecm+(2+np_global+np_global*ivs)*nrec
          else !3D
            irecm=irecm+(2+np_global+np_global*nvrt*ivs)*nrec
          endif

          if(irecm>mirec) then
            write(errmsg,*)'Output file too large',i,irecm
            call parallel_abort(errmsg)
          endif      
        endif !iof(i)=1

      enddo !i=1,noutput

! Test output
      igmp=0
      if(noutgm==1) then

        !All ranks initialize record counter
        igmp=(4+3+3)*8/nbyte

        !Rank 0 open new file
        if(myrank==0) then
          open(100,file=ifile_char(1:ifile_len)//'_test.60', &
               access='direct',recl=nbyte,status='replace')
          write(100,rec=igmp+1) nrec
          write(100,rec=igmp+2) ns_global
          write(100,rec=igmp+3) real(dt*nspool)
          write(100,rec=igmp+4) nspool
          write(100,rec=igmp+5) 1
          close(100)
        endif !myrank==0

        !All ranks increment record counter
        igmp=igmp+5
      endif !noutgm.eq.1

      !All ranks wait here
#ifdef INCLUDE_TIMING
!      if(cwtime) cwtmp=mpi_wtime()
#endif
      call parallel_barrier
#ifdef INCLUDE_TIMING
!      if(cwtime) wtimer(14,2)=wtimer(14,2)+mpi_wtime()-cwtmp
#endif

      end subroutine write_header0

!===============================================================================
!===============================================================================

      subroutine write_obe
!-------------------------------------------------------------------------------
! Output centers.bp and sidecenters.bp
! NOTE: Valid for single processor only!
!-------------------------------------------------------------------------------
      use elfe_glbl
      use elfe_msgp
      implicit none
      integer :: i,j

      !Output sidecenters.bp
      open(32,file='sidecenters.bp')
      if(iwrite.eq.0) then
        write(32,*) 'Sidegrid'
        write(32,*) ns
      else !evm
        write(32,"(a)",advance="no") 'Sidegrid\n'
        write(32,"(i12,a)",advance="no") ns,'\n'
      endif
      do i=1,ns
        if(iwrite.eq.0) then
          write(32,*) i,xcj(i),ycj(i),real(dps(i))
        else !evm
          write(32,"(i12,a,f19.9,a,f19.9,a,f12.6,a)",advance="no") &
             i," ",xcj(i),"      ",ycj(i),"    ",real(dps(i)),'\n'
        endif
      enddo !i
      close(32)

      !Output centers.bp
      open(32,file='centers.bp')
      if(iwrite.eq.0) then
        write(32,*) 'centers pts'
        write(32,*) ne
      else !evm
        write(32,"(a)",advance="no") 'centers pts\n'
        write(32,"(i12,a)",advance="no") ne,'\n'
      endif
      do i=1,ne
        if(iwrite.eq.0) then
          write(32,*) i,xctr(i),yctr(i),real(dpe(i))
        else !evm
          write(32,"(i12,a,f19.9,a,f19.9,a,f12.6,a)",advance="no") &
             i," ",xctr(i),"      ",yctr(i),"    ",real(dpe(i)),'\n'
        endif
      enddo !i
      close(32)

!      write(*,*)'Pre-processing completed successfully!'
!      call parallel_finalize

      end subroutine write_obe


!===============================================================================
!===============================================================================

      subroutine report_timers
!-------------------------------------------------------------------------------
! Write timing data for all tasks to timer.out file
!-------------------------------------------------------------------------------
#ifdef USE_MPIMODULE
      use mpi
#endif
      use elfe_glbl, only : rkind,mxtimer,wtimer
      use elfe_msgp
      implicit none
#ifndef USE_MPIMODULE
      include 'mpif.h'
#endif
      integer :: i
      real(rkind) :: wavg1(0:mxtimer),wavg2(0:mxtimer)
      real(rkind) :: wbuf(2,0:mxtimer)
      real(rkind) :: wmin1(2,0:mxtimer),wmin2(2,0:mxtimer)
      real(rkind) :: wmax1(2,0:mxtimer),wmax2(2,0:mxtimer)
!-------------------------------------------------------------------------------

      ! Sum communication time for timestepping section
      do i=3,13; wtimer(2,2)=wtimer(2,2)+wtimer(i,2); enddo;

      ! Total communication time
      wtimer(0,2)=wtimer(1,2)+wtimer(2,2)

      ! Make computation time excluding communication time
      wtimer(:,1)=wtimer(:,1)-wtimer(:,2)

      ! Compute average time over all tasks
      call mpi_allreduce(wtimer(0,1),wavg1,mxtimer+1,rtype,MPI_SUM,comm,ierr)
      wavg1=wavg1/dble(nproc)
      call mpi_allreduce(wtimer(0,2),wavg2,mxtimer+1,rtype,MPI_SUM,comm,ierr)
      wavg2=wavg2/dble(nproc)

      ! Compute min & max computation time over all tasks
      wbuf(1,:)=wtimer(:,1); wbuf(2,:)=myrank;
      call mpi_allreduce(wbuf,wmin1,mxtimer+1,MPI_2DOUBLE_PRECISION,MPI_MINLOC,comm,ierr)
      call mpi_allreduce(wbuf,wmax1,mxtimer+1,MPI_2DOUBLE_PRECISION,MPI_MAXLOC,comm,ierr)

      ! Compute min & max communication time over all tasks
      wbuf(1,:)=wtimer(:,2); wbuf(2,:)=myrank;
      call mpi_allreduce(wbuf,wmin2,mxtimer+1,MPI_2DOUBLE_PRECISION,MPI_MINLOC,comm,ierr)
      call mpi_allreduce(wbuf,wmax2,mxtimer+1,MPI_2DOUBLE_PRECISION,MPI_MAXLOC,comm,ierr)

      ! Rank 0 create new file and write header and avg time
      if(myrank==0) then
        ! Open new file
        open(36,file='timer.out',form='formatted',status='replace')

        ! Write ledger
        write(36,'(2a)') '# ','********** Timer Index Mapping **********'
        write(36,'(2a)') '# ','00 -- Total'
        write(36,'(2a)') '# ','01 -- Init Section'
        write(36,'(2a)') '# ','02 -- Timestepping Section'
        write(36,'(2a)') '# ','03 -- Forcings & Prep Section'
        write(36,'(2a)') '# ','04 -- Backtracking Section'
        write(36,'(2a)') '# ','05 -- Turbulence Closure Section'
        write(36,'(2a)') '# ','06 -- Matrix Preparation Section'
        write(36,'(2a)') '# ','07 -- Wave-Cont. Solver Section'
        write(36,'(2a)') '# ','08 -- Momentum Eqs. Solve Section'
        write(36,'(2a)') '# ','09 -- Transport Section'
        write(36,'(2a)') '# ','10 -- Recomputing Levels Section'
        write(36,'(2a)') '# ','11 -- Conservation Check Section'
        write(36,'(2a)') '# ','12 -- Global Output Section'
        write(36,'(2a)') '# ','13 -- Hotstart Section'
!'

        ! Write average, min & max times
        write(36,'(/)')
        write(36,'(2a)') '# ','********** Average, Min & Max Times in secs **********'
!'
        write(36,'(11a)') 'ID', &
     '        CompT','     MinCompT',' RankMinCompT','     MaxCompT',' RankMaxCompT', &
     '        CommT','     MinCommT',' RankMinCommT','     MaxCommT',' RankMaxCommT'
        do i=0,13
          write(36,'(i2.2,2(e13.6,2(e13.6,i13)))') i, &
          wavg1(i),wmin1(1,i),int(wmin1(2,i)),wmax1(1,i),int(wmax1(2,i)), &
          wavg2(i),wmin2(1,i),int(wmin2(2,i)),wmax2(1,i),int(wmax2(2,i))
        enddo

        ! Close file
        if(nproc>1) close(36)

      endif !myrank=0

      ! Initiate round-robin synchronization
      call parallel_rrsync(1)

      ! Open file to append
      if(nproc>1) &
        open(36,file='timer.out',form='formatted',status='old',position='append')

      ! Task 0 write next header
      if(myrank==0) then
        write(36,'(/)')
        write(36,'(a)') '# ********** Computation Times (sec) For Each MPI Task **********'
        write(36,'(a)') '# ********** Rows = Ranks; Columns = Timers      **********'
        write(36,'(a,20i13)') '# Rank',(i,i=0,13)
      endif

      ! Each task write its own timer data
      write(36,'(a,i4.4,20e13.6)') '# ',myrank,(wtimer(i,1),i=0,13)

      ! Close file
      if(nproc>1) close(36)

      ! Round-robin synchronization step
      call parallel_rrsync(2)

      ! Open file to append
      if(nproc>1) &
        open(36,file='timer.out',form='formatted',status='old',position='append')

      ! Task 0 write next header
      if(myrank==0) then
        write(36,'(/)')
        write(36,'(a)') '# ********** Communication Times For Each MPI Task **********'
        write(36,'(a)') '# ********** Rows = Ranks; Columns = Timers        **********'
        write(36,'(a,20i13)') '# Rank',(i,i=0,13)
      endif

      ! Each task write its own timer data
      write(36,'(a,i4.4,20e13.6)') '# ',myrank,(wtimer(i,2),i=0,13)

      ! Close file
      if(nproc>1) close(36)

      ! Round-robin synchronization final
      call parallel_rrsync(-2)

      end subroutine report_timers


!===============================================================================
!===============================================================================
! Routine to read in param.in
!===============================================================================
!===============================================================================
      subroutine get_param(varname,vartype,ivarvalue,varvalue1,varvalue2)
! Get a parameter from param.in
! Inputs:
!        varname: parameter name (string no longer than 90)
!        vartype: parameter value type (0: 2-char string; 1: integer; 2: float)
! Outputs:
!        ivarvalue: integer output;
!        varvalue1: float output;
!        varvalue2: 2-char string output.
! Format rules for param.in:
! (1) Lines beginning with "!" are comments; blank lines are ignored;
! (2) one line for each parameter in the format: keywords= value;
!     keywords are case sensitive; spaces allowed between keywords and "=" and value;
!     comments starting with "!"  allowed after value;
! (3) value is an integer, double, or 2-char string; for double, any of the format is acceptable:
!     40 40. 4.e1
!     Use of decimal point in integers is OK but discouraged.
      use elfe_glbl, only : rkind,errmsg
      use elfe_msgp, only : parallel_abort,myrank
      implicit real(rkind)(a-h,o-z), integer(i-n)

      character(*),intent(in) :: varname
      integer,intent(in) :: vartype
      integer,intent(out) :: ivarvalue
      real(rkind),intent(out) :: varvalue1
      character(len=2),intent(out) :: varvalue2

      character(len=90) :: line_str,str_tmp,str_tmp2

      str_tmp2=adjustl(varname)
      lstr_tmp2=len_trim(str_tmp2)
!      print*, varname !,str_tmp2(1:lstr_tmp2)

      ! Scan param.in
      open(15,file='param.in',status='old')
      rewind(15)
      line=0
      do
        line=line+1
        read(15,'(a)',end=99)line_str
        line_str=adjustl(line_str) !place blanks at end
        len_str=len_trim(line_str)
        if(len_str==0.or.line_str(1:1)=='!') cycle

        loc=index(line_str,'=')
        loc2=index(line_str,'!')
        if(loc2/=0.and.loc2-1<loc+1) call parallel_abort('READ_PARAM: comments before =')
!'

        str_tmp=''
        str_tmp(1:loc-1)=line_str(1:loc-1) !keyword
        str_tmp=trim(str_tmp)
        lstr_tmp=len_trim(str_tmp)
    
        if(str_tmp(1:lstr_tmp)==str_tmp2(1:lstr_tmp2)) then
          if(loc2/=0) then
            str_tmp2=line_str(loc+1:loc2-1)
          else
            str_tmp2=line_str(loc+1:len_str)
          endif
          str_tmp2=adjustl(str_tmp2)
          str_tmp2=trim(str_tmp2)
          if(vartype==0) then !string
            varvalue2=str_tmp2(1:2)
#ifdef DEBUG
            if(myrank==0) write(86,*)varname,' = ',varvalue2
#endif
          else if(vartype==1) then !integer
            read(str_tmp2,*)ivarvalue
#ifdef DEBUG
            if(myrank==0) write(86,*)varname,' = ',ivarvalue
#endif
          else if(vartype==2) then !float
            read(str_tmp2,*)varvalue1
#ifdef DEBUG
            if(myrank==0) write(86,*)varname,' = ',real(varvalue1)
#endif
          else
            write(errmsg,*)'read_param: unknown type:',vartype
            call parallel_abort(errmsg)
          endif
          exit
        endif
      enddo !scan param.in
  
!     print*, 'Found it on line: ',line
      close(15)
      return

99    close(15)
      write(errmsg,*)'Failed to find parameter:',varname
      call parallel_abort(errmsg)

      end subroutine get_param

!===============================================================================
!===============================================================================
! END ELFE FILE I/O SUBROUTINES
!===============================================================================
!===============================================================================
