!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_SIMPLE_EXPLICIT()
        USE DATAPOOL

#ifdef SELFE
        use elfe_msgp
#endif 
        IMPLICIT NONE

        INTEGER    :: IS, ID, IP
        REAL       :: ACLOC(MSC,MDC)
        REAL       :: SSBRL(MSC,MDC)
        REAL       :: TIME1, TIME2, TIME3, TIME4, TIME5
        REAL       :: TIME6, TIME7, TIME8, TIME9, TIME10, TIME11, TIME12, TIME13

#ifdef SELFE
!        REAL*8     :: WILD(MNP)
#endif 

         AC1 = AC2

         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'START COMPUTE'

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
 
         CALL CPU_TIME(TIME1)

         IF (LDIFR) THEN
           IF (IDIFFR == 1 ) THEN
             CALL DIFFRA_SIMPLE
           ELSE IF (IDIFFR == 2) THEN
             CALL DIFFRA_EXTENDED
           END IF 
         END IF 

         CALL CPU_TIME(TIME2)

         IF (FMETHOD .GT. 0 .AND. (LSECU .OR. LSTCU .OR. LSEWL) ) CALL COMPUTE_FREQUENCY()

         IF (LDIFR) THEN
           IF (IDIFFR == 1 ) THEN
             CALL DIFFRA_SIMPLE
           ELSE IF (IDIFFR == 2) THEN
             CALL DIFFRA_EXTENDED
           END IF
         END IF

         CALL CPU_TIME(TIME3)

         IF (DMETHOD .GT. 0) CALL COMPUTE_DIRECTION()

         CALL CPU_TIME(TIME4)

         IF (LDIFR) THEN
           IF (IDIFFR == 1 ) THEN
             CALL DIFFRA_SIMPLE
           ELSE IF (IDIFFR == 2) THEN
             CALL DIFFRA_EXTENDED
           END IF
         END IF

         IF (AMETHOD .GT. 0) CALL COMPUTE_SPATIAL()

         IF (LDIFR) THEN
           IF (IDIFFR == 1 ) THEN
             CALL DIFFRA_SIMPLE
           ELSE IF (IDIFFR == 2) THEN
             CALL DIFFRA_EXTENDED
           END IF
         END IF

         CALL CPU_TIME(TIME5)

         IF (SMETHOD .GT. 0) CALL COMPUTE_SOURCES_EXP()

         CALL CPU_TIME(TIME6)

         IF (LMAXETOT .AND. SMETHOD .EQ. 0) CALL BREAK_LIMIT_ALL ! Miche for no source terms

         CALL CPU_TIME(TIME7)

         !WRITE(*,*) SUM(AC2)

#ifdef SELFE
      IF (myrank == 0) THEN
#endif
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-----SIMPLE SPLITTING SCHEME-----'
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'DIFFRACTION                      ', TIME2-TIME1
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'SOURCES                          ', TIME6-TIME5
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS ADVEKTION            ', TIME5-TIME4
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS THETA SPACE          ', TIME4-TIME3
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS SIGMA SPACE          ', TIME3-TIME2
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU MICHE LIMITER                ', TIME7-TIME6
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS TOTAL TIME           ', TIME7-TIME1
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-------------TIMINGS-------------'
        CALL FLUSH(STAT%FHNDL)
#ifdef SELFE
      ENDIF
#endif 

        IF (.NOT. LDIFR) LCALC = .FALSE.

!        CALL PLOT_SHADED_CONTOUR_POLAR(SPSIG/PI2,SPDIR*RADDEG,MSC,MDC,AC2(137,:,:),10,MSC,MDC,'BEFORE ANY CALL')    
        RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_DOUBLE_STRANG_EXPLICIT()
        USE DATAPOOL

#ifdef SELFE
        use elfe_msgp
#endif 
        IMPLICIT NONE

        INTEGER    :: IS, ID, IP
        REAL       :: ACLOC(MSC,MDC)
        REAL       :: TIME1, TIME2, TIME3, TIME4, TIME5
        REAL       :: TIME6, TIME7, TIME8, TIME9, TIME10
        REAL       :: TIME11, TIME12, TIME13, TIME14, TIME15, TIME16, TIME17

#ifdef SELFE
!        REAL*8     :: WILD(MNP)
#endif 

         !WRITE(*,'("+TRACE...",A)') 'START COMPUTER'

         IF (.NOT. LSTEA .AND. .NOT. LQSTEA) THEN
           DT4A = 0.5*MAIN%DELT
           DT4S = DT4A * 0.5
           DT4D = ONETHIRD*MAIN%DELT
           DT4F = DT4D 
         ELSE IF (LQSTEA) THEN
           DT4A = 0.5*DT_ITER
           DT4S = DT4A * 0.5
           DT4D = ONETHIRD*DT_ITER
           DT4F = DT4D 
         END IF

         CALL CPU_TIME(TIME1)

         IF (LDIFR) THEN
           IF (IDIFFR == 1 ) THEN
             CALL DIFFRA_SIMPLE
           ELSE IF (IDIFFR == 3) THEN
             CALL DIFFRA_EXTENDED
           END IF 
         END IF 

         CALL CPU_TIME(TIME2)

         IF (FMETHOD .GT. 0 .AND. (LSECU .OR. LSTCU .OR. LSEWL) ) CALL COMPUTE_FREQUENCY

         CALL CPU_TIME(TIME3)

         IF (DMETHOD .GT. 0) CALL COMPUTE_DIRECTION()

         CALL CPU_TIME(TIME4)

         !WRITE(*,'("+TRACE...",A)') 'FINISHED SPECTRAL PART -1-'

         IF (LDIFR) THEN
           IF (IDIFFR == 1 ) THEN
             CALL DIFFRA_SIMPLE
           ELSE IF (IDIFFR == 2) THEN
             CALL DIFFRA_EXTENDED
           END IF
         END IF

         CALL CPU_TIME(TIME5)

         IF (SMETHOD .GT. 0) CALL COMPUTE_SOURCES_EXP()


         IF (AMETHOD .GT. 0) CALL COMPUTE_SPATIAL()

         CALL CPU_TIME(TIME6)

         IF (SMETHOD .GT. 0) CALL COMPUTE_SOURCES_EXP()

         CALL CPU_TIME(TIME7)

         !WRITE(*,'("+TRACE...",A)') 'FINISHED SPATIAL PART AND SOURCE -1-'

         IF (LDIFR) THEN
           IF (IDIFFR == 1 ) THEN
             CALL DIFFRA_SIMPLE
           ELSE IF (IDIFFR == 2) THEN
             CALL DIFFRA_EXTENDED
           END IF
         END IF

         CALL CPU_TIME(TIME8)

         IF (DMETHOD .GT. 0) CALL COMPUTE_DIRECTION()

         CALL CPU_TIME(TIME9)

         IF (FMETHOD .GT. 0 .AND. (LSECU .OR. LSTCU .OR. LSEWL) ) CALL COMPUTE_FREQUENCY

         !WRITE(*,'("+TRACE...",A)') 'FINISHED SPECTRAL PART -2-'

         CALL CPU_TIME(TIME10)

         IF (LDIFR) THEN
           IF (IDIFFR == 1 ) THEN
             CALL DIFFRA_SIMPLE
           ELSE IF (IDIFFR == 2) THEN
             CALL DIFFRA_EXTENDED
           END IF
         END IF

         CALL CPU_TIME(TIME11)

         IF (SMETHOD .GT. 0) CALL COMPUTE_SOURCES_EXP()

         IF (AMETHOD .GT. 0) CALL COMPUTE_SPATIAL()

         CALL CPU_TIME(TIME12)

         IF (SMETHOD .GT. 0) CALL COMPUTE_SOURCES_EXP()

         CALL CPU_TIME(TIME13)

         !WRITE(*,'("+TRACE...",A)') 'FINISHED SPATIAL PART AND SOURCE -2-'

         IF (LDIFR) THEN
           IF (IDIFFR == 1 ) THEN
             CALL DIFFRA_SIMPLE
           ELSE IF (IDIFFR == 2) THEN
             CALL DIFFRA_EXTENDED
           END IF
         END IF

         CALL CPU_TIME(TIME14)

         IF (FMETHOD .GT. 0 .AND. (LSECU .OR. LSTCU .OR. LSEWL) ) CALL COMPUTE_FREQUENCY

         CALL CPU_TIME(TIME15)

         IF (DMETHOD .GT. 0) CALL COMPUTE_DIRECTION()

         CALL CPU_TIME(TIME16)

         !WRITE(*,'("+TRACE...",A)') 'FINISHED SPECTRAL PART -3-'

#ifdef SELFE
!       DO IS = 1, MSC
!         DO ID = 1, MDC
!           WILD = DBLE(AC2(:,IS,ID))
!           CALL EXCHANGE_P2D(WILD)
!           AC2(:,IS,ID) = REAL(WILD)
!         END DO ! ID
!       END DO ! IS
#endif 

         CALL CPU_TIME(TIME17)

#ifdef SELFE
      IF (myrank == 0) THEN
#endif
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '------DOUBLE STRANG SPLITTING----'
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'DIFFRACTION                      ', &
     &                            TIME2-TIME1 + TIME5-TIME4 + TIME8-TIME7 + TIME14+TIME13
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'SOURCES                          ', &
     &                            TIME7-TIME6 + TIME13-TIME12
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS ADVEKTION            ', &
     &                            TIME6-TIME5 + TIME12-TIME11
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS THETA SPACE          ', &
     &                            TIME4-TIME3 + TIME9-TIME8 + TIME16-TIME15
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS SIGMA SPACE          ', &
     &                            TIME3-TIME2 + TIME10-TIME9 + TIME15-TIME14
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS TOTAL TIME           ', &
     &                            TIME17-TIME1
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-------------TIMINGS-------------'
#ifdef SELFE
      ENDIF
#endif 

!        WRITE(*,'("+TRACE...",A,F15.6)') 'CPU TIMINGS TOTAL TIME           ', TIME2-TIME1

        IF (.NOT. LDIFR) LCALC = .FALSE.

!        DO ID = 1, MDC
!          WRITE(*,*) SPDIR(ID)*RADDEG, SUM(AC2(137,:,ID)), IOBPD(ID,137), IOBPDS(ID,137)
!        END DO
!        CALL PLOT_SHADED_CONTOUR_POLAR(SPSIG/PI2,SPDIR*RADDEG,MSC,MDC,AC2(137,:,:),10,MSC,MDC,'BEFORE ANY CALL')    
        RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_IMP()
        USE DATAPOOL
#ifdef SELFE
        use elfe_msgp
#endif 
        IMPLICIT NONE

        REAL, SAVE       :: TIME1, TIME2, TIME3, TIME4, TIME5, TIME6, TIME7
        INTEGER          :: IP, IT

!       Three level splitting to reduce splitting errors when using implicit-explicit schemes

         DT4A = MAIN%DELT
         DT4S = DT4A
         DT4D = DT4A*0.5 
         DT4F = DT4A*0.5

         AC1  = AC2
         CALL CPU_TIME(TIME1)
         IF (LDIFR) THEN
           IF (IDIFFR == 1 ) THEN
             CALL DIFFRA_SIMPLE
           ELSE IF (IDIFFR == 2) THEN
             CALL DIFFRA_EXTENDED
           END IF
         END IF
         CALL CPU_TIME(TIME2)
         IF (FMETHOD .GT. 0 .AND. (LSECU .OR. LSTCU .OR. LSEWL) ) CALL COMPUTE_FREQUENCY
         IF (DMETHOD .GT. 0) CALL COMPUTE_DIRECTION
         IF (FMETHOD .GT. 0 .AND. (LSECU .OR. LSTCU .OR. LSEWL) ) CALL COMPUTE_FREQUENCY
         IF (DMETHOD .GT. 0) CALL COMPUTE_DIRECTION
         CALL CPU_TIME(TIME4)
         IF (SMETHOD .GT. 0) CALL COMPUTE_SOURCES_IMP
         IF (AMETHOD .GT. 0) CALL COMPUTE_SPATIAL
         CALL CPU_TIME(TIME5)
         IF (LLIMT .AND. SMETHOD .GT. 0) CALL ACTION_LIMITER
         IF (.NOT. LDIFR) LCALC = .FALSE.
         CALL CPU_TIME(TIME7)
         IF (LMAXETOT) CALL BREAK_LIMIT_ALL ! Enforce Miche  


#ifdef SELFE
      IF (myrank == 0) THEN
#endif
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-----IMPLICIT SPLITTING SCHEME-----'
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'DIFFRACTION                      ', TIME2-TIME1
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'SOURCES                          ', TIME5-TIME4
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS ADVEKTION            ', TIME6-TIME5
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS THETA SPACE          ', TIME4-TIME3
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS SIGMA SPACE          ', TIME3-TIME2
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CPU TIMINGS TOTAL TIME           ', TIME7-TIME1
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-------------TIMINGS-------------'
#ifdef SELFE
      ENDIF
#endif 

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

         WRITE(*,*) SUM(DAC_ADV), SUM(DAC_SOU), SUM(DAC_THE)

         AC2 = AC1

         IITERSPLIT = 1

! iteration ...

         CALL COMPUTE_SPATIAL()
         CALL COMPUTE_FREQUENCY
         CALL COMPUTE_DIRECTION()
         CALL COMPUTE_SOURCES_EXP()

         WRITE(*,*) SUM(DAC_ADV), SUM(DAC_SOU), SUM(DAC_THE)

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
      SUBROUTINE COMPUTE_SPATIAL()
        USE DATAPOOL
        IMPLICIT NONE

        IF (DIMMODE == 1) THEN
          CALL COMPUTE_ADVECTION1D_QUICKEST_A()
        ELSE IF (DIMMODE == 2) THEN
          IF (LVECTOR) THEN
            CALL FLUCT_3
          ELSE
            CALL FLUCT_1
          END IF
          IF ( ICOMP .GE. 1 .AND. (AMETHOD .EQ. 2 .OR. AMETHOD .EQ. 3 )) CALL RESCALE_SPECTRUM
        END IF

      END SUBROUTINE COMPUTE_SPATIAL
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_DIRECTION()
        USE DATAPOOL
        IMPLICIT NONE
 
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

      END SUBROUTINE COMPUTE_DIRECTION
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_FREQUENCY()
        USE DATAPOOL
        IMPLICIT NONE

        IF (FMETHOD == 1) THEN
          CALL COMPUTE_FREQUENCY_QUICKEST_A()
        END IF

      END SUBROUTINE COMPUTE_FREQUENCY
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_SOURCES_EXP()
        USE DATAPOOL
        IMPLICIT NONE

          IF (ICOMP < 2 .AND. SMETHOD > 0) THEN
            CALL SOURCE_INT_EXP()
          END IF

      END SUBROUTINE COMPUTE_SOURCES_EXP
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_SOURCES_IMP()
        USE DATAPOOL
        IMPLICIT NONE

          IF (ICOMP >= 2  .AND. SMETHOD > 0) THEN
            CALL SOURCE_INT_IMP()
          END IF

      END SUBROUTINE COMPUTE_SOURCES_IMP
!**********************************************************************
!*                                                                    *
!**********************************************************************
