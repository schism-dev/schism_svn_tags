#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_SIMPLE_EXPLICIT(SSBR_DUMON_ALL,QB_DUMON_ALL)
        USE DATAPOOL
#ifdef MPI_PARALL_GRID
        use elfe_msgp
#endif

        IMPLICIT NONE

        REAL              :: TIME1, TIME2, TIME3, TIME4, TIME5
        REAL              :: TIME6, TIME7, TIME8, TIME9, TIME10, TIME11, TIME12, TIME13


!!!!!! modif AD
        REAL(rkind),INTENT(OUT)  :: SSBR_DUMON_ALL(MNP,MSC,MDC),QB_DUMON_ALL(MNP)
!!!!!! end modif AD


         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'START COMPUTE COMPUTE_SIMPLE_EXPLICIT'
         CALL FLUSH(STAT%FHNDL)

         AC1 = AC2 
         IF (LNANINFCHK) THEN
           WRITE(DBG%FHNDL,*) ' AFTER ENTERING COMPUTE ',  SUM(AC2)
           IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN COMPUTE 1') 
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

         CALL CPU_TIME(TIME1)

         CALL COMPUTE_DIFFRACTION

         CALL CPU_TIME(TIME2)

         IF (FMETHOD .GT. 0 .AND. (LSECU .OR. LSTCU .OR. LSEWL) ) CALL COMPUTE_FREQUENCY()
         IF (DMETHOD .GT. 0) CALL COMPUTE_DIRECTION()

         IF (LNANINFCHK) THEN
           WRITE(DBG%FHNDL,*) ' AFTER DIRECTION AND FREQUENCY -1- ',  SUM(AC2)
           IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN COMPUTE 2')
         ENDIF
  
         CALL CPU_TIME(TIME3)

         IF (FMETHOD .GT. 0 .AND. (LSECU .OR. LSTCU .OR. LSEWL) ) CALL COMPUTE_FREQUENCY()
         IF (DMETHOD .GT. 0) CALL COMPUTE_DIRECTION()

         IF (LNANINFCHK) THEN
           WRITE(DBG%FHNDL,*) ' AFTER DIRECTION AND FREQUENCY -2- ',  SUM(AC2)
           IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN COMPUTE 3') 
         ENDIF

         CALL CPU_TIME(TIME4)

         IF (AMETHOD .GT. 0) CALL COMPUTE_SPATIAL()

         IF (LNANINFCHK) THEN
           WRITE(DBG%FHNDL,*) ' AFTER SPATIAL ',  SUM(AC2)
           IF (SUM(AC2) .NE. SUM(AC2)) STOP 'NAN IN COMPUTE 4'
         ENDIF

         CALL CPU_TIME(TIME5)
!!!!!! modif AD
         IF (SMETHOD .GT. 0) CALL COMPUTE_SOURCES_EXP_DUMON(SSBR_DUMON_ALL,QB_DUMON_ALL)
!!!!!! end modif AD
         IF (LNANINFCHK) THEN
           WRITE(DBG%FHNDL,*) ' AFTER SOURCES ',  SUM(AC2)
           IF (SUM(AC2) .NE. SUM(AC2)) STOP 'NAN IN COMPUTE 5'
         ENDIF

         CALL CPU_TIME(TIME6)

         IF (LMAXETOT .AND. SMETHOD .EQ. 0) CALL BREAK_LIMIT_ALL ! Miche for no source terms ... may cause oscilations ...

         IF (LNANINFCHK) THEN
           WRITE(DBG%FHNDL,*) ' AFTER BREAK LIMIT ',  SUM(AC2)
           IF (SUM(AC2) .NE. SUM(AC2)) STOP 'NAN IN COMPUTE 6'
         ENDIF

         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-----SIMPLE SPLITTING SCHEME-----'
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'DIFFRACTION                      ', TIME2-TIME1
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'SOURCES                          ', TIME6-TIME5
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS ADVEKTION            ', TIME5-TIME4
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS THETA SPACE          ', TIME4-TIME3
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS SIGMA SPACE          ', TIME3-TIME2
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU MICHE LIMITER                ', TIME6-TIME5
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS TOTAL TIME           ', TIME6-TIME1
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-------------TIMINGS-------------'

         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'FINISHED COMPUTE COMPUTE_SIMPLE_EXPLICIT'
         CALL FLUSH(STAT%FHNDL)

        IF (.NOT. LDIFR) LCALC = .FALSE.

!        CALL PLOT_SHADED_CONTOUR_POLAR(SPSIG/PI2,SPDIR*RADDEG,MSC,MDC,AC2(137,:,:),10,MSC,MDC,'BEFORE ANY CALL')    
        RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_DOUBLE_STRANG_EXPLICIT()
        USE DATAPOOL

#ifdef MPI_PARALL_GRID
        use elfe_msgp
#endif
        IMPLICIT NONE

        REAL       :: TIME1, TIME2, TIME3, TIME4, TIME5
        REAL       :: TIME6, TIME7, TIME8, TIME9, TIME10
        REAL       :: TIME11, TIME12, TIME13, TIME14, TIME15, TIME16, TIME17

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'START COMPUTE'

        IF (.NOT. LSTEA .AND. .NOT. LQSTEA) THEN
          DT4A = 0.5*MAIN%DELT
          DT4S = 0.5*DT4A
          DT4D = ONETHIRD*MAIN%DELT
          DT4F = DT4D 
        ELSE IF (LQSTEA) THEN
          DT4A = 0.5*DT_ITER
          DT4S = DT4A * 0.25
          DT4D = ONETHIRD*DT_ITER
          DT4F = DT4D 
        END IF

        CALL CPU_TIME(TIME1)

        CALL COMPUTE_DIFFRACTION

        IF (FMETHOD .GT. 0 .AND. (LSECU .OR. LSTCU .OR. LSEWL) ) CALL COMPUTE_FREQUENCY
        IF (DMETHOD .GT. 0) CALL COMPUTE_DIRECTION
! ---- 1st spectra 
        IF (AMETHOD .GT. 0) CALL COMPUTE_SPATIAL
        IF (SMETHOD .GT. 0) CALL COMPUTE_SOURCES_EXP
        IF (DMETHOD .GT. 0) CALL COMPUTE_DIRECTION
        IF (FMETHOD .GT. 0 .AND. (LSECU .OR. LSTCU .OR. LSEWL) ) CALL COMPUTE_FREQUENCY
! ---- 2nd spectra 
        IF (AMETHOD .GT. 0) CALL COMPUTE_SPATIAL()
        IF (FMETHOD .GT. 0 .AND. (LSECU .OR. LSTCU .OR. LSEWL) ) CALL COMPUTE_FREQUENCY
        IF (DMETHOD .GT. 0) CALL COMPUTE_DIRECTION()
! ---- 3rd spectra 
        IF (AMETHOD .GT. 0) CALL COMPUTE_SPATIAL()
        IF (SMETHOD .GT. 0) CALL COMPUTE_SOURCES_EXP()

        CALL CPU_TIME(TIME17)

        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '------DOUBLE STRANG SPLITTING----'
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'DIFFRACTION                      ', TIME2-TIME1 + TIME5-TIME4 + TIME8-TIME7 + TIME14+TIME13
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'SOURCES                          ', TIME7-TIME6 + TIME13-TIME12
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS ADVEKTION            ', TIME6-TIME5 + TIME12-TIME11
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS THETA SPACE          ', TIME4-TIME3 + TIME9-TIME8 + TIME16-TIME15
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS SIGMA SPACE          ', TIME3-TIME2 + TIME10-TIME9 + TIME15-TIME14 
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS TOTAL TIME           ', TIME17-TIME1
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-------------TIMINGS-------------'

        IF (.NOT. LDIFR) LCALC = .FALSE.

        RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_IMPLICIT
        USE DATAPOOL
        IMPLICIT NONE

        REAL, SAVE       :: TIME1, TIME2, TIME3, TIME4, TIME5, TIME6, TIME7, TIME8, TIME9, TIME10, TIME11
        INTEGER          :: IP, IT


        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING COMPUTE COMPUTE_IMPLICIT'
        CALL FLUSH(STAT%FHNDL)

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

        AC1  = AC2

        CALL CPU_TIME(TIME1)
        CALL COMPUTE_DIFFRACTION
        CALL CPU_TIME(TIME2)
        IF (DMETHOD .GT. 0) CALL COMPUTE_DIRECTION
        IF (FMETHOD .GT. 0 .AND. (LSECU .OR. LSTCU .OR. LSEWL) ) CALL COMPUTE_FREQUENCY
        IF (DMETHOD .GT. 0) CALL COMPUTE_DIRECTION
        IF (FMETHOD .GT. 0 .AND. (LSECU .OR. LSTCU .OR. LSEWL) ) CALL COMPUTE_FREQUENCY
        CALL CPU_TIME(TIME3)
        IF (LMAXETOT) CALL BREAK_LIMIT_ALL ! Enforce Miche
        CALL CPU_TIME(TIME4)
        IF (SMETHOD .GT. 0) CALL COMPUTE_SOURCES_IMP
        CALL CPU_TIME(TIME5)
        IF (AMETHOD .GT. 0) CALL COMPUTE_SPATIAL
        CALL CPU_TIME(TIME6)
        IF (LLIMT .AND. SMETHOD .GT. 0) CALL ACTION_LIMITER
        CALL CPU_TIME(TIME7)
        IF (LMAXETOT) CALL BREAK_LIMIT_ALL ! Enforce Miche  
        CALL CPU_TIME(TIME8)

        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-----IMPLICIT SPLITTING SCHEME-----'
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'DIFFRACTION                      ', TIME2-TIME1
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS ADVEKTION            ', TIME6-TIME5
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS SPECTRAL SPACE       ', TIME3-TIME2
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS SOURCES              ', TIME5-TIME4
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'ACTION LIMITER                   ', TIME7-TIME6
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'MICHE LIMITER                    ', TIME8-TIME7+TIME4-TIME3
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS TOTAL TIME           ', TIME8-TIME1
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-------------TIMINGS-------------'

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'FINISHED COMPUTE COMPUTE_IMPLICIT'
        CALL FLUSH(STAT%FHNDL)

        RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE COMPUTE_ITERATIVE_SPLITTING()
        USE DATAPOOL
        IMPLICIT NONE

        INTEGER          :: ITER, IS, ID, IP

        REAL, SAVE       :: TIME1, TIME2

         CALL CPU_TIME(TIME1)

         DT4A = MAIN%DELT 
         DT4S = DT4A
         DT4D = DT4A
         DT4F = DT4A

! Set DAC's to Zero ... 

         CALL CPU_TIME(TIME1)

         DAC_THE = 0.
         DAC_SIG = 0.
         DAC_SOU = 0.
         DAC_ADV = 0.

         AC1 = AC2

! 1st step ...

         IITERSPLIT = 0

         CALL COMPUTE_SPATIAL()
         CALL COMPUTE_FREQUENCY
         CALL COMPUTE_DIRECTION()
         CALL COMPUTE_SOURCES_EXP()

         WRITE(DBG%FHNDL,*) SUM(DAC_ADV), SUM(DAC_SOU), SUM(DAC_THE)

         AC2 = AC1

         IITERSPLIT = 1

! iteration ...

         CALL COMPUTE_SPATIAL()
         CALL COMPUTE_FREQUENCY
         CALL COMPUTE_DIRECTION()
         CALL COMPUTE_SOURCES_EXP()

         WRITE(DBG%FHNDL,*) SUM(DAC_ADV), SUM(DAC_SOU), SUM(DAC_THE)

         AC2 = AC1

! final step ... 

         CALL COMPUTE_SPATIAL()
         CALL COMPUTE_FREQUENCY
         CALL COMPUTE_DIRECTION()
         CALL COMPUTE_SOURCES_EXP()

         CALL CPU_TIME(TIME2)

        RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_DIFFRACTION
        USE DATAPOOL
        IF (LDIFR) THEN
          IF (IDIFFR == 1 ) THEN
            CALL DIFFRA_SIMPLE
          ELSE IF (IDIFFR == 2) THEN
            CALL DIFFRA_EXTENDED
          END IF
        END IF
      END SUBROUTINE COMPUTE_DIFFRACTION
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_SPATIAL()
        USE DATAPOOL
        IMPLICIT NONE

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING COMPUTE COMPUTE_SPATIAL'
        CALL FLUSH(STAT%FHNDL)

        IF (DIMMODE == 1) THEN
          CALL COMPUTE_ADVECTION1D_QUICKEST_A()
        ELSE IF (DIMMODE == 2) THEN
          IF (LVECTOR) THEN
            CALL FLUCT_3
          ELSE
            CALL FLUCT_1
! don't forget to uncomment FLUCT* in wwm_fluctsplit
!             IF(ICOMP == 0) THEN
!               CALL FLUCT_EXPLICIT()
!             ELSE IF(ICOMP == 1) THEN
!               CALL FLUCT_SEMIIMPLICIT()
!             ELSE IF(ICOMP == 2) THEN
!               CALL FLUCT_IMPLICIT()
!             ENDIF
          END IF
          IF ( ICOMP .GE. 1 .AND. (AMETHOD .EQ. 2 .OR. AMETHOD .EQ. 3 )) CALL RESCALE_SPECTRUM
        END IF

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'FINISHED COMPUTE_SPATIAL'
        CALL FLUSH(STAT%FHNDL)

      END SUBROUTINE COMPUTE_SPATIAL
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_DIRECTION()
        USE DATAPOOL
        IMPLICIT NONE

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING COMPUTE_DIRECTION'
        CALL FLUSH(STAT%FHNDL)
 
        IF (DMETHOD > 0) THEN
          IF (DMETHOD == 1) THEN
            CALL COMPUTE_DIRECTION_CNTG_A()
          ELSE IF (DMETHOD == 2) THEN
            CALL COMPUTE_DIRECTION_QUICKEST_A()
          ELSE IF (DMETHOD == 3) THEN
            CALL COMPUTE_DIRECTION_WENO_A()
          ELSE IF (DMETHOD == 4) THEN
            CALL COMPUTE_DIRECTION_UPWIND_A()
          END IF
        END IF

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'FINISHED COMPUTE_DIRECTION'
        CALL FLUSH(STAT%FHNDL)

        IF ( DMETHOD == 1) CALL RESCALE_SPECTRUM

      END SUBROUTINE COMPUTE_DIRECTION
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_FREQUENCY()
        USE DATAPOOL
        IMPLICIT NONE

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING COMPUTE_FREQUENCY'
        CALL FLUSH(STAT%FHNDL)

        IF (FMETHOD == 1) THEN
          CALL COMPUTE_FREQUENCY_QUICKEST_A()
        END IF

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'FINISHED COMPUTE_FREQUENCY'
        CALL FLUSH(STAT%FHNDL)

      END SUBROUTINE COMPUTE_FREQUENCY
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_SOURCES_EXP_DUMON(SSBR_DUMON_ALL,QB_DUMON_ALL)
        USE DATAPOOL
        IMPLICIT NONE

!!!!! modif AD
        REAL(rkind), INTENT(OUT) :: SSBR_DUMON_ALL(MNP,MSC,MDC),QB_DUMON_ALL(MNP) 
!!!!! end modif AD

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING COMPUTE_SOURCES_EXP'
        CALL FLUSH(STAT%FHNDL)

        IF (ICOMP < 2 .AND. SMETHOD > 0) THEN
!!!!! modif AD
          CALL SOURCE_INT_EXP(SSBR_DUMON_ALL,QB_DUMON_ALL)
!!!!! modif AD
        END IF

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'FINISHED COMPUTE_SOURCES_EXP'
        CALL FLUSH(STAT%FHNDL)

      END SUBROUTINE COMPUTE_SOURCES_EXP_DUMON
!!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_SOURCES_EXP()
        USE DATAPOOL
        IMPLICIT NONE

!!!!! modif AD
        REAL(rkind) :: SSBR_DUMON_ALL(MNP,MSC,MDC),QB_DUMON_ALL(MNP) 
!!!!! end modif AD

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING COMPUTE_SOURCES_EXP'
        CALL FLUSH(STAT%FHNDL)

        IF (ICOMP < 2 .AND. SMETHOD > 0) THEN
!!!!! modif AD
          CALL SOURCE_INT_EXP(SSBR_DUMON_ALL,QB_DUMON_ALL)
!!!!! modif AD
        END IF

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'FINISHED COMPUTE_SOURCES_EXP'
        CALL FLUSH(STAT%FHNDL)

      END SUBROUTINE COMPUTE_SOURCES_EXP
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_SOURCES_IMP()
        USE DATAPOOL
        IMPLICIT NONE

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING COMPUTE_SOURCES_IMP'
        CALL FLUSH(STAT%FHNDL)

        IF (ICOMP >= 2  .AND. SMETHOD > 0) THEN
          CALL SOURCE_INT_IMP()
        END IF

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'FINISHED COMPUTE_SOURCES_IMP'
        CALL FLUSH(STAT%FHNDL)

      END SUBROUTINE COMPUTE_SOURCES_IMP
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CFLSPEC()
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER              :: IP

         REAL(rkind)                 :: TMPCFLCAD(MNP), TMPCAD(MNP)
         REAL(rkind)                 :: TMPCFLCAS(MNP), TMPCAS(MNP)
         REAL(rkind)                 :: CAS(MSC,MDC), CAD(MSC,MDC)

         OPEN(310, FILE='cflcad.bin', FORM = 'UNFORMATTED', STATUS = 'UNKNOWN')
         OPEN(311, FILE='cflcas.bin', FORM = 'UNFORMATTED', STATUS = 'UNKNOWN')

         TMPCFLCAS = 0.
         TMPCFLCAD = 0.
         TMPCAS    = 0. 
         TMPCAD    = 0.

         DO IP = 1, MNP
           IF (DEP(IP) .GT. DMIN) THEN
             CALL PROPTHETA(IP,CAD)
             CALL PROPSIGMA(IP,CAS)
             TMPCAD(IP)    = MAXVAL(ABS(CAD))
! 0.5 since the directional and frequency intergration is split in two parts ....
             TMPCFLCAD(IP) = 0.5 * TMPCAD(IP)*MAIN%DELT/DDIR
             TMPCAS(IP)    = MAXVAL(ABS(CAS))
! absolute max. value ... lies on the secure side ... to do ...
             TMPCFLCAS(IP) = 0.5 * TMPCAS(IP)*MAIN%DELT/MINVAL(DS_INCR)
           ELSE
             CALL PROPTHETA(IP,CAD)
             CALL PROPSIGMA(IP,CAS)
             TMPCFLCAD(IP) = 0.
             TMPCAD(IP)    = 0.
             TMPCFLCAS(IP) = 0.
             TMPCAS(IP)    = 0.
           END IF
         END DO

         MAXCFLCAD = MAXVAL(TMPCAD)
         MAXCFLCAS = MAXVAL(TMPCAS)

         WRITE (310) RTIME
         WRITE (310) (TMPCAD(IP), TMPCAD(IP), TMPCFLCAD(IP), IP = 1, MNP)
         WRITE (311) RTIME
         WRITE (311) (TMPCAS(IP), TMPCAS(IP), TMPCFLCAS(IP), IP = 1, MNP)

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
