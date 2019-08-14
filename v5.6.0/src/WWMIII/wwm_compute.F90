#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_EXPLICIT
        USE DATAPOOL
        IMPLICIT NONE

        REAL(rkind)       :: VEC2RAD
        REAL(rkind)       :: SSBRL(NUMSIG,NUMDIR), SSLIM(NUMSIG,NUMDIR)
#ifdef TIMINGS
        REAL(rkind)       :: TIME1, TIME2, TIME3, TIME4, TIME5
        REAL(rkind)       :: TIME6, TIME7, TIME8
#endif

         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'START COMPUTE COMPUTE_SIMPLE_EXPLICIT'
         FLUSH(STAT%FHNDL)
#ifdef TIMINGS
         CALL WAV_MY_WTIME(TIME1)
#endif
         IF (LNANINFCHK) THEN
           WRITE(DBG%FHNDL,*) ' AFTER ENTERING COMPUTE ',  SUM(AC2)
           IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN COMPUTE 1') 
           IF (MINVAL(AC2) .LT. ZERO) CALL WWM_ABORT(' NEGATIVE IN COMPUTE 1')
         ENDIF

         IF (.NOT. LSTEA .AND. .NOT. LQSTEA) THEN
           DT4A = MAIN%DELT
           DT4S = DT4A
           DT4D = 0.5_rkind*DT4A
           DT4F = 0.5_rkind*DT4A 
         ELSE IF (LQSTEA) THEN
           DT4A = DT_ITER
           DT4S = DT4A
           DT4D = 0.5_rkind*DT4A
           DT4F = 0.5_rkind*DT4A
         END IF

#ifdef TIMINGS
         CALL WAV_MY_WTIME(TIME2)
#endif
         CALL COMPUTE_DIFFRACTION
#ifdef DEBUG
         CALL Print_SumAC2("After COMPUTE_DIFFRACTION")
#endif
#ifdef TIMINGS
         CALL WAV_MY_WTIME(TIME3)
#endif
         IF (FMETHOD .GT. 0) CALL COMPUTE_FREQUENCY
#ifdef DEBUG
         CALL Print_SumAC2("After COMPUTE_FREQUENCY 1")
#endif
         IF (DMETHOD .GT. 0) CALL COMPUTE_DIRECTION
#ifdef DEBUG
         CALL Print_SumAC2("After COMPUTE_DIRECTION 1")
#endif
         IF (LNANINFCHK) THEN
           WRITE(DBG%FHNDL,*) ' AFTER DIRECTION AND FREQUENCY -1- ',  SUM(AC2)
           IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN COMPUTE 2')
           IF (MINVAL(AC2) .LT. ZERO) CALL WWM_ABORT(' NEGATIVE IN COMPUTE 2')
         ENDIF
  
#ifdef TIMINGS
         CALL WAV_MY_WTIME(TIME4)
#endif
         IF (AMETHOD .GT. 0) CALL COMPUTE_SPATIAL
#ifdef DEBUG
         CALL Print_SumAC2("After COMPUTE_SPATIAL")
#endif

         IF (LNANINFCHK) THEN
           WRITE(DBG%FHNDL,*) ' AFTER SPATIAL ',  SUM(AC2)
           IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN COMPUTE 3') 
           IF (MINVAL(AC2) .LT. ZERO) CALL WWM_ABORT(' NEGATIVE IN COMPUTE 3')
         ENDIF

#ifdef TIMINGS
         CALL WAV_MY_WTIME(TIME5)
#endif
         IF (FMETHOD .GT. 0) CALL COMPUTE_FREQUENCY
#ifdef DEBUG
         CALL Print_SumAC2("After COMPUTE_FREQUENCY 2")
#endif
         IF (DMETHOD .GT. 0) CALL COMPUTE_DIRECTION
#ifdef DEBUG
         CALL Print_SumAC2("After COMPUTE_DIRECTION 2")
#endif

         IF (LNANINFCHK) THEN
           WRITE(DBG%FHNDL,*) ' AFTER DIRECTION AND FREQUENCY -2-  ',  SUM(AC2)
           IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN COMPUTE 4')
           IF (MINVAL(AC2) .LT. ZERO) CALL WWM_ABORT(' NEGATIVE IN COMPUTE 4')
         ENDIF

#ifdef TIMINGS
         CALL WAV_MY_WTIME(TIME6)
#endif

#ifdef DEBUG
         CALL Print_SumAC2("Before COMPUTE_SOURCES")
#endif
         AC1 = AC2
         IF (SMETHOD .GT. 0) CALL COMPUTE_SOURCES
#ifdef DEBUG
         CALL Print_SumAC2(" After COMPUTE_SOURCES")
#endif
         IF (LNANINFCHK) THEN
           WRITE(DBG%FHNDL,*) ' AFTER SOURCES ',  SUM(AC2)
           IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN COMPUTE 5')
           IF (MINVAL(AC2) .LT. ZERO) CALL WWM_ABORT(' NEGATIVE IN COMPUTE 5')
         ENDIF

#ifdef TIMINGS
         CALL WAV_MY_WTIME(TIME7)
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-----SIMPLE SPLITTING SCHEME-----'
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'DIFFRACTION                      ', TIME3-TIME2
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'SOURCES                          ', TIME7-TIME6
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS ADVEKTION            ', TIME5-TIME4
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS SPECTRAL SPACE       ', TIME6-TIME5 + TIME4-TIME3
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU MICHE LIMITER                ', TIME6-TIME5 ! ???????
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS TOTAL TIME           ', TIME7-TIME1
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-------------TIMINGS-------------'
#endif
         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'FINISHED COMPUTE COMPUTE_SIMPLE_EXPLICIT'
         FLUSH(STAT%FHNDL)
!        CALL PLOT_SHADED_CONTOUR_POLAR(SPSIG/PI2,SPDIR*RADDEG,NUMSIG,NUMDIR,AC2(137,:,:),10,NUMSIG,NUMDIR,'BEFORE ANY CALL')
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_SEMI_IMPLICIT
        USE DATAPOOL
        IMPLICIT NONE

        REAL(rkind), SAVE  :: TIME1, TIME2, TIME3, TIME4, TIME5, TIME6, TIME7, TIME8, GTEMP1, GTEMP2
        REAL(rkind)        :: SSBRL(NUMSIG,NUMDIR), SSLIM(NUMSIG,NUMDIR)


        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING COMPUTE COMPUTE_SEMI_IMPLICIT'
        FLUSH(STAT%FHNDL)

        IF (.NOT. LSTEA .AND. .NOT. LQSTEA) THEN
          DT4A = MAIN%DELT
          DT4S = DT4A
          DT4D = 0.5_rkind*DT4A
          DT4F = 0.5_rkind*DT4A
        ELSE IF (LQSTEA) THEN
          DT4A = DT_ITER
          DT4S = DT4A
          DT4D = 0.5_rkind*DT4A
          DT4F = 0.5_rkind*DT4A
        END IF

        AC1 = AC2

#ifdef TIMINGS
        CALL WAV_MY_WTIME(TIME1)
#endif

        CALL COMPUTE_DIFFRACTION

#ifdef TIMINGS
        CALL WAV_MY_WTIME(TIME2)
#endif

        IF (DMETHOD .GT. 0) CALL COMPUTE_DIRECTION
        IF (FMETHOD .GT. 0) CALL COMPUTE_FREQUENCY
        IF (DMETHOD .GT. 0) CALL COMPUTE_DIRECTION
        IF (FMETHOD .GT. 0) CALL COMPUTE_FREQUENCY

#ifdef TIMINGS
        CALL WAV_MY_WTIME(TIME3)
#endif

        IF (SMETHOD .GT. 0) CALL COMPUTE_SOURCES 

#ifdef TIMINGS
        CALL WAV_MY_WTIME(TIME5)
#endif
        IF (AMETHOD .GT. 0) CALL COMPUTE_SPATIAL
!
#ifdef TIMINGS
        CALL WAV_MY_WTIME(TIME6)
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-----IMPLICIT SPLITTING SCHEME-----'
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'DIFFRACTION                      ', TIME2-TIME1
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS ADVEKTION            ', TIME6-TIME5
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS SPECTRAL SPACE       ', TIME3-TIME2
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS SOURCES              ', TIME5-TIME4
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS TOTAL TIME           ', TIME8-TIME1
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-------------TIMINGS-------------'
        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'FINISHED COMPUTE COMPUTE_SEMI_IMPLICIT'
        FLUSH(STAT%FHNDL)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_IMPLICIT
        USE DATAPOOL
#ifdef PETSC
        USE PETSC_BLOCK, ONLY : EIMPS_PETSC_BLOCK
#endif
        IMPLICIT NONE
#ifdef TIMINGS
        REAL(rkind)       :: TIME1, TIME2, TIME3, TIME4, TIME5
        REAL(rkind)       :: TIME6, TIME7, TIME8
#endif
        REAL(rkind)       :: SSBRL(NUMSIG,NUMDIR), SSLIM(NUMSIG,NUMDIR)

       IF (.NOT. LSTEA .AND. .NOT. LQSTEA) THEN
         DT4A = MAIN%DELT
         DT4S = DT4A
         DT4D = DT4A
         DT4F = DT4A
       ELSE IF (LQSTEA) THEN
         DT4A = DT_ITER
         DT4S = DT4A
         DT4D = DT4A
         DT4F = DT4A
       END IF

       AC1 = AC2

       IF (LNANINFCHK) THEN
         WRITE(DBG%FHNDL,*) ' AFTER ENTERING COMPUTE ',  SUM(AC2)
         IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN COMPUTE 1')
       ENDIF

#ifdef TIMINGS
        CALL WAV_MY_WTIME(TIME1)
#endif
        CALL COMPUTE_DIFFRACTION

        IF (LNANINFCHK) THEN
          WRITE(DBG%FHNDL,*) ' AFTER DIFFRACTION',  SUM(AC2)
          IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN COMPUTE 2')
        ENDIF
#ifdef TIMINGS
        CALL WAV_MY_WTIME(TIME2)
#endif
#ifdef TIMINGS
        CALL WAV_MY_WTIME(TIME3)
#endif
        CALL SOURCES_IMPLICIT

        IF (LNANINFCHK) THEN
          WRITE(DBG%FHNDL,*) ' AFTER SOURCES',  SUM(AC2)
          IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN COMPUTE 3a')
          IF (SUM(PHIA) .NE. SUM(PHIA)) CALL WWM_ABORT('NAN IN COMPUTE 3b')
          IF (SUM(DPHIDNA) .NE. SUM(DPHIDNA)) CALL WWM_ABORT('NAN IN COMPUTE 3c')
        ENDIF
#ifdef TIMINGS
        CALL WAV_MY_WTIME(TIME4)
#endif
        IF (AMETHOD .eq.5) THEN
#ifdef PETSC
          CALL EIMPS_PETSC_BLOCK
#endif
        ELSE IF (AMETHOD .eq. 7) THEN
#ifdef WWM_SOLVER
          CALL EIMPS_TOTAL_JACOBI_ITERATION
#endif
        END IF

        IF (LMAXETOT) CALL BREAKING_LIMITER_GLOBAL 

        IF (LNANINFCHK) THEN
          WRITE(DBG%FHNDL,*) ' AFTER ADVECTION',  SUM(AC2)
          IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN COMPUTE 4')
        ENDIF

        IF (LMAXETOT) CALL BREAKING_LIMITER_GLOBAL

#ifdef TIMINGS
        CALL WAV_MY_WTIME(TIME5)
#endif
#ifdef TIMINGS
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-----IMPLICIT -----'
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'DIFFRACTION                      ', TIME2-TIME1
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS SOLVER               ', TIME5-TIME4
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS SOURCES              ', TIME4-TIME3
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS TOTAL TIME           ', TIME5-TIME1
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-------------TIMINGS-------------'
        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'FINISHED COMPUTE COMPUTE_IMPLICIT'
#endif

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_SPATIAL
        USE DATAPOOL
        IMPLICIT NONE

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING COMPUTE COMPUTE_SPATIAL'
        FLUSH(STAT%FHNDL)

        IF (DIMMODE == 1) THEN
          CALL COMPUTE_ADVECTION1D_QUICKEST_A
        ELSE IF (DIMMODE == 2) THEN
          IF(ICOMP == 0) THEN
            CALL FLUCT_EXPLICIT
          ELSE IF(ICOMP == 1) THEN
            CALL FLUCT_IMP_EXP_SOURCES
          ELSE IF(ICOMP == 2) THEN
            CALL FLUCT_IMP_SOURCES
          ELSE IF(ICOMP == 3) THEN 
            CALL FLUCT_IMP_ALL
          ENDIF
          IF ( ICOMP .GE. 1 .AND. (AMETHOD .EQ. 2 .OR. AMETHOD .EQ. 3 )) CALL RESCALE_SPECTRUM
        END IF

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'FINISHED COMPUTE_SPATIAL'
        FLUSH(STAT%FHNDL)

      END SUBROUTINE COMPUTE_SPATIAL
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_SOURCES
        USE DATAPOOL
        IMPLICIT NONE

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING COMPUTE_SOURCES_EXP'
        FLUSH(STAT%FHNDL)

        IF (ICOMP < 2) THEN
          CALL SOURCES_EXPLICIT
        ELSEIF (ICOMP  .GE. 2) THEN
          CALL SOURCES_IMPLICIT
        ENDIF 

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'FINISHED COMPUTE_SOURCES_EXP'
        FLUSH(STAT%FHNDL)

      END SUBROUTINE COMPUTE_SOURCES
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CFLSPEC
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER              :: IP

         REAL(rkind)                 :: TMPCFLCAD(MNP), TMPCAD(MNP)
         REAL(rkind)                 :: TMPCFLCAS(MNP), TMPCAS(MNP)
         REAL(rkind)                 :: CAS(NUMSIG,NUMDIR), CAD(NUMSIG,NUMDIR)

         OPEN(310, FILE='cflcad.bin', FORM = 'UNFORMATTED', STATUS = 'UNKNOWN')
         OPEN(311, FILE='cflcas.bin', FORM = 'UNFORMATTED', STATUS = 'UNKNOWN')

         TMPCFLCAS = 0.
         TMPCFLCAD = 0.
         TMPCAS    = 0. 
         TMPCAD    = 0.

         DO IP = 1, MNP
           CALL PROPTHETA(IP,CAD)
           CALL PROPSIGMA(IP,CAS)
           TMPCAD(IP)    = MAXVAL(ABS(CAD))
! 0.5 since the directional and frequency intergration is split in two parts .... no of course not. Why should it be so, u solve 2 times a 1d equation and this has just the normal cfl 
           TMPCFLCAD(IP) = TMPCAD(IP)*MAIN%DELT/DDIR
           TMPCAS(IP)    = MAXVAL(ABS(CAS))
! absolute max. value ... lies on the secure side ... to do ...
           TMPCFLCAS(IP) = TMPCAS(IP)*MAIN%DELT/MINVAL(DS_INCR)
           CFL_CASD(1,IP)=TMPCAS(IP)
           CFL_CASD(2,IP)=TMPCAD(IP)
           CFL_CASD(3,IP)=TMPCFLCAS(IP)
           CFL_CASD(4,IP)=TMPCFLCAD(IP)
         END DO

         MAXCFLCAD = MAXVAL(TMPCAD)
         MAXCFLCAS = MAXVAL(TMPCAS)

         WRITE (310) SNGL(RTIME)
         WRITE (310) (SNGL(TMPCAD(IP)), SNGL(TMPCAD(IP)), SNGL(TMPCFLCAD(IP)), IP = 1, MNP)
         WRITE (311) SNGL(RTIME)
         WRITE (311) (SNGL(TMPCAS(IP)), SNGL(TMPCAS(IP)), SNGL(TMPCFLCAS(IP)), IP = 1, MNP)
         
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
