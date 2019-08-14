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
!===============================================================================
! SCHISM FILE I/O SUBROUTINES
!
! subroutine write_obe
! subroutine report_timers

!===============================================================================
!===============================================================================

      subroutine write_obe
!-------------------------------------------------------------------------------
! Output centers.bp and sidecenters.bp
! NOTE: Valid for single processor only!
!-------------------------------------------------------------------------------
      use schism_glbl
      use schism_msgp
      implicit none
      integer :: i,j
      real(rkind) ::  tmp1,tmp2

      !Output sidecenters.bp
      open(32,file='sidecenters.bp')
      write(32,*) 'Sidegrid'
      write(32,*) ns
      do i=1,ns
        if(ics==1) then
          write(32,*) i,real(xcj(i)),real(ycj(i)),real(dps(i))
        else
          tmp1=sum(xlon(isidenode(1:2,i)))/2
          tmp2=sum(ylat(isidenode(1:2,i)))/2
          write(32,*) i,real(tmp1),real(tmp2),real(dps(i))
        endif
      enddo !ics
      close(32)

      !Output centers.bp
      open(32,file='centers.bp')
      write(32,*) 'centers pts'
      write(32,*) ne
      do i=1,ne
        if(ics==1) then
          write(32,'(i12,2(1x,e22.14),e12.4)') i,xctr(i),yctr(i),dpe(i)
        else
          tmp1=sum(xlon(elnode(1:i34(i),i)))/i34(i)
          tmp2=sum(ylat(elnode(1:i34(i),i)))/i34(i)
          write(32,*) i,real(tmp1),real(tmp2),real(dpe(i))
        endif !ics
      enddo !i
      close(32)

      end subroutine write_obe


!===============================================================================
!===============================================================================

      subroutine report_timers
!-------------------------------------------------------------------------------
! Write timing data for all tasks to timer.out file
!-------------------------------------------------------------------------------
!#ifdef USE_MPIMODULE
!      use mpi
!#endif
      use schism_glbl, only : rkind,mxtimer,wtimer
      use schism_msgp
      implicit none
!#ifndef USE_MPIMODULE
      include 'mpif.h'
!#endif
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
! END FILE I/O SUBROUTINES
!===============================================================================
!===============================================================================
