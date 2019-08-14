      subroutine selfe_finalize
      use elfe_glbl
      use elfe_msgp
      USE hydraulic_structures
#ifdef USE_HA
      USE harm
#endif
      implicit none
      include 'mpif.h'

      integer :: i,iwbl_itmax_gb
      character(len=72) :: it_char
      character(len=8) :: date
      character(len=10) :: timestamp
      real(rkind) :: wtmp2

!...  Output max. diffusion number reached
      if(difnum_max_l2>0.5) write(12,*)'MAIN: max. diffusion # exceeds 0.5:',difnum_max_l2
!'

!...  Output max. # of iterations for all ranks for WBL (Grant-Madsen formulation)
#ifdef USE_WWM
        if(iwbl==1) then
           call mpi_reduce(iwbl_itmax,iwbl_itmax_gb,1,itype,MPI_MAX,0,comm,ierr)
           if(myrank==0) write(16,*)'Max. iteration for Grant-Madsen = ',iwbl_itmax_gb
        endif !iwbl
#endif /*USE_WWM*/

!...  Output max. elevations & dahv
      fdb='maxelev_0000'
      lfdb=len_trim(fdb)
      write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
      open(10,file='outputs/'//fdb,status='replace')
      write(10,*)np,nproc
      do i=1,np
        write(10,*)iplg(i),real(xnd(i)),real(ynd(i)),real(elevmax(i))
      enddo !i
      close(10)

      fdb='maxdahv_0000'
      lfdb=len_trim(fdb)
      write(fdb(lfdb-3:lfdb),'(i4.4)') myrank
      open(10,file='outputs/'//fdb,status='replace')
      write(10,*)np,nproc
      do i=1,np
        write(10,*)iplg(i),real(xnd(i)),real(ynd(i)),real(dav_maxmag(i)),real(dav_max(1:2,i))
      enddo !i
      close(10)

#ifdef USE_HA
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! SOLVE THE HARMONIC ANALYSIS PROBLEM (Adapted from ADCIRC)
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!     Use 12 to write from each process
      write(12,*)'HA',iharind,it_main,ITHAS,CHARMV
      IF ((iharind.EQ.1).AND.(it_main.GT.ITHAS)) THEN
!...Compute means and variances for checking the harmonic analysis results
!...Accumulate mean and variance at each node.
        if (CHARMV) then
          IF (FMV.NE.0.) THEN
            DO I=1,np
              ELAV(I)   = ELAV(I)/NTSTEPS
              XVELAV(I) = XVELAV(I)/NTSTEPS
              YVELAV(I) = YVELAV(I)/NTSTEPS
              ELVA(I)   = ELVA(I)/NTSTEPS   - ELAV(I)*ELAV(I)
              XVELVA(I) = XVELVA(I)/NTSTEPS - XVELAV(I)*XVELAV(I)
              YVELVA(I) = YVELVA(I)/NTSTEPS - YVELAV(I)*YVELAV(I)
            END DO
            TIMEBEG=ITMV*dt
            write(it_char(1:4),'(i4.4)') myrank
            open(55,file='outputs/harme.55'//it_char(1:4))
            WRITE(55,*) np
          ENDIF
        endif

!......Fill out and decompose the LHS harmonic analaysis matrix

        CALL FULSOL(0)

!......Solve the harmonic analysis problem and write the output

        write(12,*)'myrank=',myrank,TIMEBEG,ITMV,dt
        IF(NHAGE.EQ.1) CALL LSQSOLEG(np,ELAV,ELVA,CHARMV,myrank,TIMEBEG,DT,FMV,NTSTEPS,ITMV)
        IF(NHAGV.EQ.1) CALL LSQSOLVG(np,XVELAV,YVELAV,XVELVA,YVELVA,CHARMV,myrank,TIMEBEG,DT,FMV,NTSTEPS,ITMV)
      ENDIF
#endif /*USE_HA*/

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
! Finalize
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      if(ihydraulics/=0) call finalize_hydraulic_structures

#ifdef INCLUDE_TIMING
      wtmp2=mpi_wtime()
      wtimer(2,1)=wtmp2-wtimer(2,1) !time-stepping section
      wtimer(0,1)=wtmp2-wtimer(0,1) !total
!     Report timing
      call report_timers

!     Report wall-clock time of some targeted routines
      do i=1,size(timer_ns)
        write(12,*)'Custom timers in hours:',i,myrank,real(timer_ns(i)/3600)
      enddo !i
#endif

      call date_and_time(date,timestamp)
      if(myrank==0) write(16,'(/4a)') 'Run completed successfully at ',date,', ',timestamp

!     Close file handles - may have problem on some systems
      do i=1,2000
        close(i)
      enddo !i

#ifdef USE_WWM
      call TERMINATE_WWM
#endif /*USE_WWM*/

      end subroutine selfe_finalize
