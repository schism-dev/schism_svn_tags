!:
! __      __  __      __  _____  .___.___.___ 
!/  \    /  \/  \    /  \/     \ |   |   |   |
!\   \/\/   /\   \/\/   /  \ /  \|   |   |   |
! \        /  \        /    Y    \   |   |   |
!  \__/\  /    \__/\  /\____|__  /___|___|___|
!       \/          \/         \/             
!
!
! WWM-III (Wind Wave Model) source code 
! 
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef WWMONLY
      PROGRAM WWMII

         USE DATAPOOL

         IMPLICIT NONE

         INTEGER :: K

         CALL INITIALIZE_WWM

         DO K = 1, MAIN%ISTP
            IF (LQSTEA) THEN
              CALL QUASI_STEADY(K)
            ELSE
              CALL UN_STEADY(K)
            END IF
         END DO

         STOP 'END OF SIMULATION'

      END PROGRAM
#endif 
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef SELFE
      SUBROUTINE WWM_II(IT_SELFE,icou_elfe_wwm,DT_SELFE0,NSTEP_WWM0)

         USE DATAPOOL
         use elfe_msgp!, only : myrank,parallel_abort,itype,comm,ierr
         use elfe_glbl, only : iplg,ielg

         IMPLICIT NONE

         INTEGER, INTENT(IN)   :: NSTEP_WWM0, icou_elfe_wwm
         REAL*8, INTENT(IN)    :: DT_SELFE0

         REAL, SAVE  :: SIMUTIME
         REAL*8      :: T1, T2
         REAL        :: TIME1, TIME2, TIME3, TIME4, TIME5, TIME6, TIME7

         INTEGER     :: I, IP, IT_SELFE, K, IFILE, IT, ITMP, IPGL
         REAL        :: OUTPAR(OUTVARS), OUTWINDPAR(WINDVARS), ACLOC(MSC,MDC)
         REAL*8      :: DTMP

         LOGICAL, PARAMETER :: LNANCHECK = .TRUE.

         CALL CPU_TIME(TIME1) 

!zyl: check dimension
!ar: obsolete since already checked in init_wwm, after model has been tested this should be removed
         if(WINDVARS/=size(WIND_INTPAR,2)) call parallel_abort('Dimension mismatch: OUTWINDPAR and out_wwm_windpar')
         if(OUTVARS/=size(OUTT_INTPAR,2)) call parallel_abort('Dimension mismatch: OUTPAR and out_wwm')

         NSTEPWWM = NSTEP_WWM0

         DT_SELFE  = SNGL(DT_SELFE0)

         T1 = DBLE(IT_SELFE-NSTEPWWM)*DT_SELFE0 ! Beginn time step ...
         T2 = DBLE(IT_SELFE)*DT_SELFE0          ! End of time time step ...

         MAIN%DELT = NSTEPWWM*DT_SELFE 

         SIMUTIME = SIMUTIME + MAIN%DELT

         IF (icou_elfe_wwm == 1) THEN ! Full coupling 
           WLDEP       = REAL(DEP8)
           WATLEV      = REAL(ETA2)
           WATLEVOLD   = REAL(ETA1)
           DEP         = MAX(0.,WLDEP + WATLEV)
           CURTXY(:,1) = REAL(UU2(NVRT,:))
           CURTXY(:,2) = REAL(VV2(NVRT,:))
           IF (.NOT. LWINDFROMWWM) THEN
             WINDXY(:,1) = REAL(WINDX0)
             WINDXY(:,2) = REAL(WINDY0)
           END IF
           LSECU       = .TRUE.
           LSEWL       = .TRUE.
           LCALC       = .TRUE.
         ELSE IF (icou_elfe_wwm == 0) THEN ! No interaction at all 
           WLDEP       = SNGL(DEP8)
           WATLEV      = 0. 
           WATLEVOLD   = 0.
           DEP         = WLDEP
           CURTXY(:,1) = 0.!REAL(UU2(NVRT,:))
           CURTXY(:,2) = 0.!REAL(VV2(NVRT,:))
           IF (.NOT. LWINDFROMWWM) THEN
             WINDXY(:,1) = REAL(WINDX0)
             WINDXY(:,2) = REAL(WINDY0)
           END IF
           LSECU       = .FALSE.
           LSEWL       = .FALSE.
           LCALC       = .TRUE. 
         ELSE IF (icou_elfe_wwm == 2) THEN ! current effect in wwm but no radiation stress in SELFE 
           WLDEP       = REAL(DEP8)
           WATLEV      = REAL(ETA2)
           WATLEVOLD   = REAL(ETA1)
           DEP         = MAX(0.,WLDEP + WATLEV)
           CURTXY(:,1) = REAL(UU2(NVRT,:))
           CURTXY(:,2) = REAL(VV2(NVRT,:))
           IF (.NOT. LWINDFROMWWM) THEN
             WINDXY(:,1) = REAL(WINDX0)
             WINDXY(:,2) = REAL(WINDY0)
           END IF
           LSECU       = .TRUE.
           LSEWL       = .TRUE.
           LCALC       = .TRUE.
         ELSE IF (icou_elfe_wwm == 3) THEN ! No current effect in wwm but radiation stress in SELFE 
           WLDEP       = SNGL(DEP8)
           WATLEV      = 0.
           WATLEVOLD   = 0.
           DEP         = WLDEP
           CURTXY(:,1) = 0.!REAL(UU2(NVRT,:))
           CURTXY(:,2) = 0.!REAL(VV2(NVRT,:))
           IF (.NOT. LWINDFROMWWM) THEN
             WINDXY(:,1) = REAL(WINDX0)
             WINDXY(:,2) = REAL(WINDY0)
           END IF
           LSECU       = .FALSE.
           LSEWL       = .FALSE.
           LCALC       = .TRUE.
         ELSE IF (icou_elfe_wwm == 4) THEN ! No current effect in wwm and no radiation stress in SELFE 
           WLDEP       = SNGL(DEP8)
           WATLEV      = 0.
           WATLEVOLD   = 0.
           DEP         = WLDEP
           CURTXY(:,1) = 0.!REAL(UU2(NVRT,:))
           CURTXY(:,2) = 0.!REAL(VV2(NVRT,:))
           IF (.NOT. LWINDFROMWWM) THEN
             WINDXY(:,1) = REAL(WINDX0)
             WINDXY(:,2) = REAL(WINDY0)
           END IF
           LSECU       = .FALSE.
           LSEWL       = .FALSE.
           LCALC       = .TRUE.
         END IF

!         CALL CGSPEC 

         IFILE = 1
         IT    = 1

         IF (LBCSE .AND. (LBCWA .OR. LBCSP)) THEN
           IF ( MAIN%TMJD > SEBO%TMJD-1.E-8 .AND. MAIN%TMJD < SEBO%EMJD ) THEN
           !WRITE(*,*) MAIN%TMJD, SEBO%TMJD-1.E-8, SEBO%EMJD
             IF (IBOUNDFORMAT == 3) THEN ! Find the right position in the file ...
               DTMP = (MAIN%TMJD-BND_TIME_ALL_FILES(1,1)) * DAY2SEC
               !WRITE(*,*) DTMP, MAIN%BMJD, BND_TIME_ALL_FILES(1,1)
               ITMP  = 0
               DO IFILE = 1, NUM_NETCDF_FILES_BND
                 ITMP = ITMP + NDT_BND_FILE(IFILE)
                 IF (ITMP .GT. INT(DTMP/SEBO%DELT)) EXIT
               END DO
               ITMP = SUM(NDT_BND_FILE(1:IFILE-1))
               IT   = NINT(DTMP/SEBO%DELT) - ITMP + 1
               IF (LBINTER) IT = IT + 1
               IF (IT .GT. NDT_BND_FILE(IFILE)) THEN
                 IFILE = IFILE + 1
                 IT    = 1
               ENDIF
             END IF ! BOUNDFORMAT
             !WRITE(*,*) IFILE, IT, SUM(NDT_BND_FILE(1:IFILE-1)), NINT(DTMP/SEBO%DELT), SEBO%DELT
             IF (LBINTER) THEN 
               IF (IBOUNDFORMAT == 3) THEN
                 WRITE(STAT%FHNDL,*) 'CALL TO WAVE_BOUNDARY_CONDITION FROM SELFE', IFILE, IT, LBINTER
                 CALL WAVE_BOUNDARY_CONDITION(IFILE,IT,WBACNEW)
               ELSE
                 WRITE(STAT%FHNDL,*) 'CALL TO WAVE_BOUNDARY_CONDITION FROM SELFE', IFILE, IT, LBINTER
                 CALL WAVE_BOUNDARY_CONDITION(1,1,WBACNEW)
               END IF
               DSPEC   = (WBACNEW-WBACOLD)/SEBO%DELT*MAIN%DELT
               WBAC    =  WBACOLD
               WBACOLD =  WBACNEW
               !WRITE(*,*) SUM(WBAC), SUM(WBACNEW), SUM(WBACOLD)
             ELSE ! .NOT. LBINTER
               IF (IBOUNDFORMAT == 3) THEN
                 WRITE(STAT%FHNDL,*) 'CALL TO WAVE_BOUNDARY_CONDITION FROM SELFE', IFILE, IT, LBINTER
                 CALL WAVE_BOUNDARY_CONDITION(IFILE,IT,WBAC)
               ELSE
                 WRITE(STAT%FHNDL,*) 'CALL TO WAVE_BOUNDARY_CONDITION FROM SELFE', IFILE, IT, LBINTER
                 CALL WAVE_BOUNDARY_CONDITION(1,1,WBAC)
               END IF
             END IF

             SEBO%TMJD = SEBO%TMJD + DBLE(SEBO%DELT*SEC2DAY) ! Increment boundary time line ...

           ELSE ! Interpolate in time ...

             IF (LBINTER) THEN
               WBAC = WBAC + DSPEC
             END IF
 
           END IF

           IF (LINHOM) THEN
             DO IP = 1, IWBMNP
               IPGL = IWBNDLC(IP)
               AC2(IPGL,:,:) = DBLE(WBAC(:,:,IP))
             END DO
           ELSE
             DO IP = 1, IWBMNP
               IPGL = IWBNDLC(IP)
               AC2(IPGL,:,:) = DBLE(WBAC(:,:,1))
             END DO
           ENDIF

         END IF ! LBCWA .OR. LBCSP

         IF (LFIRSTSTEP) THEN
           CALL INITIAL_CONDITION(IFILE,IT-1)
           LFIRSTSTEP = .FALSE.
           LCALC      = .TRUE.
         END IF

         CALL CPU_TIME(TIME2) 

         IF (LQSTEA) THEN
            CALL QUASI_STEADY(KKK)
         ELSE
            CALL UN_STEADY(KKK)
         END IF

         CALL CPU_TIME(TIME3)
!
! Write output when requested ...
!
         IF (myrank == 0) WRITE(STAT%FHNDL,'("+TRACE...",A,F15.4)') 'nth call to WWM', SIMUTIME
         IF ( (MAIN%TMJD .GE. OUTF%TMJD-1.E-8) .AND. (MAIN%TMJD .LE. OUTF%EMJD)) THEN
           IF (IOUTP .GT. 0) THEN
             WRITE(STAT%FHNDL,'("+TRACE...",A,F15.4)') 'WRITING OUTPUT', MOD(T2,OUTF%DELT) 
           ENDIF
           CALL OUTPUT(REAL(RTIME*DAY2SEC),.FALSE.)
           OUTF%TMJD = OUTF%TMJD + DBLE(OUTF%DELT*SEC2DAY)
         ENDIF 
!
! Compute WWM output parameters for SELFE, this must be done after each call to WWMII
!
         CALL CPU_TIME(TIME4)
  
         OUTT_INTPAR = 0.
         WIND_INTPAR = 0.
         DO IP = 1, MNP
           ACLOC = AC2(IP,:,:)
           IF (DEP(IP) .GT. DMIN) THEN
             CALL INTPAR(IP, MSC, ACLOC, OUTPAR)
             OUTT_INTPAR(IP,:) = OUTPAR
             CALL WINDPAR(IP,OUTWINDPAR)
             WIND_INTPAR(IP,:) = OUTWINDPAR
             IF (LMONO_OUT) THEN
               OUTT_INTPAR(IP,1) = OUTT_INTPAR(IP,1) / SQRT(2.)
             ELSE
               OUTT_INTPAR(IP,1) = OUTT_INTPAR(IP,1) 
             END IF
           ELSE
             OUTT_INTPAR(IP,:) = 0.
             WIND_INTPAR(IP,:) = 0.
           END IF
         END DO

         CALL CPU_TIME(TIME5)
!
! Compute radiation stress ...
!
         IF (icou_elfe_wwm == 0 .OR. icou_elfe_wwm == 2 .OR. icou_elfe_wwm == 4) THEN
           WWAVE_FORCE = 0.
         ELSE 
           CALL RADIATION_STRESS
         END IF 

         CALL CPU_TIME(TIME6)
 
         IF (LNANCHECK) THEN
           DO IP = 1, MNP
             IF (SUM(OUTT_INTPAR(IP,:)) .NE. SUM(OUTT_INTPAR(IP,:))) THEN
               DO I = 1, SIZE(OUTT_INTPAR(IP,:))
                 WRITE(DBG%FHNDL,*) 'NaN in OUTT_INTPAR', IP, I, OUTT_INTPAR(IP,I)
                 CALL FLUSH(DBG%FHNDL)
               END DO
             END IF
             IF (SUM(WIND_INTPAR(IP,:)) .NE. SUM(WIND_INTPAR(IP,:))) THEN
               DO I = 1, SIZE(WIND_INTPAR(IP,:))
                 WRITE(DBG%FHNDL,*) 'NaN in WIND_INTPAR', IP, I, WIND_INTPAR(IP,I)
                 CALL FLUSH(DBG%FHNDL)
               END DO
             END IF
             IF (SUM(WWAVE_FORCE(:,IP,:)) .NE. SUM(WWAVE_FORCE(:,IP,:))) THEN
               DO I = 1, SIZE(WWAVE_FORCE(:,IP,1)) ! loop over layers ...
                 WRITE(DBG%FHNDL,*) 'NaN in WWAVE_FORCE', IP, I, WWAVE_FORCE(I,IP,1), WWAVE_FORCE(I,IP,2)
                 CALL FLUSH(DBG%FHNDL)
               END DO
             END IF 
           END DO
         END IF!LNANCHECK

         MAIN%TMJD = MAIN%TMJD + DBLE(MAIN%DELT*SEC2DAY)
         RTIME = REAL( MAIN%TMJD - MAIN%BMJD  )
         KKK = KKK + 1

         CALL CPU_TIME(TIME7)

#ifdef SELFE
      IF (myrank == 0) THEN
#endif
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-----TOTAL TIMINGS-----'
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'PREPARATION        ', TIME2-TIME1
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'INTEGRATION        ', TIME3-TIME2
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'OUTPUT WWM         ', TIME4-TIME3
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'OUTPUT TO SELFE    ', TIME5-TIME4
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'RADIATION STRESSES ', TIME6-TIME5
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'NAN CHECK          ', TIME7-TIME6
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'TOTAL TIME         ', TIME7-TIME1
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '------END-TIMINGS-  ---'
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.4)') 'FINISHED WITH WWM', SIMUTIME
        CALL FLUSH(STAT%FHNDL)
#ifdef SELFE
      ENDIF
#endif
 
      END SUBROUTINE WWM_II
#endif 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE UN_STEADY(K)

         USE DATAPOOL
#ifdef SELFE
         use elfe_msgp
#endif

         IMPLICIT NONE

         INTEGER, INTENT(IN) :: K

         REAL    :: CONV1, CONV2, CONV3, CONV4, CONV5
         REAL    :: TIME1, TIME2, TIME3, TIME4, TIME5, TIME6

         INTEGER :: IP

         CALL CPU_TIME(TIME1)

         CALL IO_1(K)

         CALL CPU_TIME(TIME2)

#ifdef WWMONLY
         IF (LCFL) CALL CFLSPEC
#endif

         CALL CPU_TIME(TIME3)

         IF (ICOMP .EQ. 0) THEN
           IF (LITERSPLIT) THEN
            CALL COMPUTE_DOUBLE_STRANG_EXPLICIT
             !CALL COMPUTE_ITERATIVE_SPLITTING
           ELSE
             CALL COMPUTE_SIMPLE_EXPLICIT
           END IF
         ELSE IF (ICOMP .EQ. 1) THEN 
           IF (LITERSPLIT) THEN
            CALL COMPUTE_DOUBLE_STRANG_EXPLICIT
           ELSE
             CALL COMPUTE_SIMPLE_EXPLICIT
           END IF
         ELSE IF (ICOMP .EQ. 2) THEN 
           CALL COMPUTE_IMP
         END IF

         CALL CPU_TIME(TIME4)

#ifdef WWMONLY
         MAIN%TMJD = MAIN%BMJD + DBLE(K*MAIN%DELT*SEC2DAY)
         RTIME = REAL( MAIN%TMJD - MAIN%BMJD  )
         WRITE(*,101)  K, MAIN%ISTP, RTIME
#endif 

         CALL IO_2(K)

         CALL CPU_TIME(TIME5)

         IF (LCONV) THEN
            CALL CHECK_STEADY(RTIME,CONV1,CONV2,CONV3,CONV4,CONV5)
         END IF

         CALL CPU_TIME(TIME6)

         IF (.NOT. LDIFR) LCALC = .FALSE.

#ifdef SELFE
      IF (myrank == 0) THEN
#endif
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-----TOTAL RUN TIMES-----'
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'PREPROCESSING                    ', TIME2-TIME1
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CHECK CFL                        ', TIME3-TIME2
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'INTEGRATION                      ', TIME4-TIME3
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'POSTPROCESSING                   ', TIME5-TIME4
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CHECK STEADY                     ', TIME6-TIME5
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-------------TIMINGS-------------'
        CALL FLUSH(STAT%FHNDL)
#ifdef SELFE
      ENDIF
#endif 

101      FORMAT ('+STEP = ',I10,'/',I10,' ( TIME = ',F15.4,' DAYS)')

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE QUASI_STEADY(K)

#ifdef SELFE
         USE elfe_msgp
#endif
         USE DATAPOOL

         IMPLICIT NONE
         INTEGER, INTENT(IN) :: K

         INTEGER :: IT, ITER
         REAL    :: ITERTIME
         REAL    :: CONV1, CONV2, CONV3, CONV4, CONV5

         CALL IO_1(K)

         IF (LCFL) CALL CFLSPEC()

#ifdef SELFE
         NQSITER = NSTEPWWM
#endif 

         DO IT = 1, NQSITER

           DT_ITER = MAIN%DELT/REAL(NQSITER)

           IF (ICOMP .LT. 2) THEN
             IF (LITERSPLIT) THEN
               CALL COMPUTE_DOUBLE_STRANG_EXPLICIT
             ELSE 
               CALL COMPUTE_SIMPLE_EXPLICIT
             END IF
           ELSE
             CALL COMPUTE_IMP
           END IF

           ITERTIME = RTIME*DAY2SEC+IT*DT_ITER

           IF (LCHKCONV) THEN
             CALL CHECK_STEADY(ITERTIME, CONV1, CONV2, CONV3, CONV4, CONV5)
             IF ( (CONV1 .GT. 100.*QSCONV1 .AND. &
     &           CONV2 .GT. 100.*QSCONV2 .AND. & 
     &           CONV3 .GT. 100.*QSCONV3 .AND. &
     &           CONV4 .GT. 100.*QSCONV4 .AND. &
     &           CONV5 .GT. 100.*QSCONV5 .AND. &
     &           K .NE. 1) .OR.  & 
     &           IT .EQ. NQSITER ) THEN
#ifdef WWMONLY
               WRITE(QSTEA%FHNDL,'(3I10,5F15.8)') K, IT, NQSITER, CONV1, CONV2, CONV3, CONV4, CONV5
#elif SELFE
               if (myrank == 0) WRITE(QSTEA%FHNDL,'(3I10,5F15.8)') K, IT, NQSITER, CONV1, CONV2, CONV3, CONV4, CONV5
#endif 
               EXIT 
             END IF
           END IF
           IF (LOUTITER) CALL OUTPUT(ITERTIME,.FALSE.)
         END DO

#ifdef SELFE 
         IF (myrank == 0) WRITE(STAT%FHNDL,101)  K, MAIN%ISTP, RTIME*DAY2SEC
#elif WWMONLY
         WRITE(STAT%FHNDL,101)  K, MAIN%ISTP, RTIME*DAY2SEC
         MAIN%TMJD = MAIN%BMJD + DBLE(K)*DBLE(MAIN%DELT)*SEC2DAY
         RTIME = REAL( MAIN%TMJD - MAIN%BMJD )
#endif

         CALL IO_2(K)

         IF (.NOT. LDIFR) LCALC = .FALSE.

101      FORMAT ('+STEP = ',I5,'/',I5,' ( TIME = ',F15.4,'HR )')
         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE IO_1(K)
#ifdef NCDF
         USE NETCDF 
#endif
         USE DATAPOOL
#ifdef SELFE
         USE elfe_msgp
#endif 
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: K
         REAL                :: TEST, TMP(MNP,2)
         REAL*8              :: DTMP, WI(3)
         CHARACTER(LEN=20)   :: STIME

         INTEGER             :: IP, ISTAT, IT, IFILE, ITMP, IPGL, FORECASTHOURS 

         FORECASTHOURS = 0 

         IF (LWINDFROMWWM) THEN
           WRITE(STAT%FHNDL,*) MAIN%TMJD, SEWI%TMJD-1.E-8, MAIN%TMJD > SEWI%TMJD-1.E-8, MAIN%TMJD < SEWI%EMJD, SEWI%EMJD
           IF ( LSEWD .AND. (MAIN%TMJD > SEWI%TMJD-1.E-8) .AND. (MAIN%TMJD < SEWI%EMJD) ) THEN
             IF (IWINDFORMAT == 1) THEN
               CALL CSEVAL( WIN%FHNDL, WIN%FNAME, LINTERWD, LWINFILE, MNP, 2, WINDXY, SEWI%DELT, MAIN%DELT, DVWIND )
#ifdef NCDF
             ELSE IF (IWINDFORMAT == 2) THEN ! DWD_NETCDF
               DTMP = (MAIN%TMJD-WIND_TIME_ALL_FILES(1)) * DAY2SEC
! This are the indices that are used the open the right NETCDF file at the right time ...
               IFILE = INT(DTMP/SEWI%DELT/DBLE(NDT_WIND_FILE)) + 1
               IT    = NINT(DTMP/SEWI%DELT)-(IFILE-1)*NDT_WIND_FILE + 1
! Since we have already read the first time step and we are just reading now the next one for interpolation ... we add + 1
               IT = IT + 1 
               IF (IT .GT. NDT_WIND_FILE) THEN
                 IFILE = IFILE + 1
                 IT    = 1
               ENDIF
               IF (IFILE .GT. NUM_NETCDF_FILES) STOP 'SOMETHING IS WRONG WE RUN OUT OF WIND TIME' 
               ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(IFILE), NF90_NOWRITE, WIND_NCID)
               CALL READ_NETCDF_DWD(IT)
               ISTAT = NF90_CLOSE(WIND_NCID)
               CALL INTER_STRUCT_DATA(NDX_WIND,NDY_WIND,DX_WIND,DY_WIND,OFFSET_X_WIND,OFFSET_Y_WIND,WIND_X,TMP(:,1))
               CALL INTER_STRUCT_DATA(NDX_WIND,NDY_WIND,DX_WIND,DY_WIND,OFFSET_X_WIND,OFFSET_Y_WIND,WIND_Y,TMP(:,2))
!               CALL INTER_STRUCT_DATA(NDX_WIND,NDY_WIND,DX_WIND,DY_WIND,OFFSET_X_WIND,OFFSET_Y_WIND,ATMO_PRESS,TMP(:,1))
               DVWIND = (TMP-WINDXY)/SEWI%DELT*MAIN%DELT
             ELSE IF (IWINDFORMAT == 3) THEN ! NOAA CFRS ...
               DTMP = (MAIN%TMJD-WIND_TIME_ALL_FILES(1)) * DAY2SEC
               IFILE = INT(DTMP/SEWI%DELT/DBLE(NDT_WIND_FILE)) + 1
               IT    = NINT(DTMP/SEWI%DELT)-(IFILE-1)*NDT_WIND_FILE + NINT((MAIN%TMJD-WIND_TIME_ALL_FILES(1)) * 4) 
               !WRITE(*,*) (MAIN%TMJD-WIND_TIME_ALL_FILES(1))*4, NINT(DTMP/SEWI%DELT)-(IFILE-1)*NDT_WIND_FILE
! Since we have already read the first time step and we are just reading now the next one for interpolation ... we add + 1
!AR: Same error CFRS is hourly ...
               !WRITE(*,*) IT 
               IF (IT .GT. NDT_WIND_FILE) THEN
                 IFILE = IFILE + 1
                 IT    = 1
               ENDIF
               IF (WIND_TIME_ALL_FILES(IT+(IFILE-1)*NDT_WIND_FILE) .EQ. WIND_TIME_ALL_FILES(IT+1+(IFILE-1)*NDT_WIND_FILE)) THEN
                 !WRITE(*,*) IT, IFILE, WIND_TIME_ALL_FILES(IT+(IFILE-1)*NDT_WIND_FILE), WIND_TIME_ALL_FILES(IT+1+(IFILE-1)*NDT_WIND_FILE)
                 IT = IT + 1 ! Skip the initital condition 
               END IF
               WRITE(STAT%FHNDL,*) 'STEPPINGS IN NARR_NETCDF'
               WRITE(STAT%FHNDL,*) IT, SEWI%DELT, MAIN%TMJD, NDT_WIND_FILE, FORECASTHOURS
               IF (IFILE .GT. NUM_NETCDF_FILES) STOP 'SOMETHING IS WRONG WE RUN OUT OF WIND TIME'
               ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(IFILE), NF90_NOWRITE, WIND_NCID)
               WRITE(STAT%FHNDL,*) 'CALLING READ_NETCDF_NARR FROM MAIN', IT
               CALL READ_NETCDF_NARR(IT)
               ISTAT = NF90_CLOSE(WIND_NCID)
               CALL INTER_STRUCT_DATA(NDX_WIND,NDY_WIND,DX_WIND,DY_WIND,OFFSET_X_WIND,OFFSET_Y_WIND,WIND_X,TMP(:,1))
               CALL INTER_STRUCT_DATA(NDX_WIND,NDY_WIND,DX_WIND,DY_WIND,OFFSET_X_WIND,OFFSET_Y_WIND,WIND_Y,TMP(:,2))
!               CALL INTER_STRUCT_DATA(NDX_WIND,NDY_WIND,DX_WIND,DY_WIND,OFFSET_X_WIND,OFFSET_Y_WIND,ATMO_PRESS,TMP(:,1))
               DVWIND = (TMP-WINDXY)/SEWI%DELT*MAIN%DELT
             ELSE IF (IWINDFORMAT == 4) THEN ! NOAA NARR ...
               DTMP = (MAIN%TMJD-WIND_TIME_ALL_FILES(1)) * DAY2SEC
               IFILE = INT(DTMP/SEWI%DELT/DBLE(NDT_WIND_FILE)) + 1
               IT    = NINT(DTMP/SEWI%DELT)-(IFILE-1)*NDT_WIND_FILE + 1
! Since we have already read the first time step and we are just reading now the next one for interpolation ... we add + 1
               IT = IT + 1
               IF (IT .GT. NDT_WIND_FILE) THEN
                 IFILE = IFILE + 1
                 IT    = 1
               ENDIF
               WRITE(STAT%FHNDL,*) 'STEPPINGS IN NARR_NETCDF'
               WRITE(STAT%FHNDL,*) IT, SEWI%DELT, MAIN%TMJD, NDT_WIND_FILE, FORECASTHOURS
               IF (IFILE .GT. NUM_NETCDF_FILES) STOP 'SOMETHING IS WRONG WE RUN OUT OF WIND TIME'
               ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(IFILE), NF90_NOWRITE, WIND_NCID)
               WRITE(STAT%FHNDL,*) 'CALLING READ_NETCDF_NARR FROM MAIN', IT
               CALL READ_NETCDF_NARR2(IT)
               ISTAT = NF90_CLOSE(WIND_NCID)
               DO IP = 1, MNP
                 CALL INTELEMENT(XYPWIND(1,INE_WIND(:,WIND_ELE(IP))),XYPWIND(2,INE_WIND(:,WIND_ELE(IP))),UWND_NARR(INE_WIND(:,WIND_ELE(IP))),XP(IP),YP(IP),Wi,TMP(IP,1),.FALSE.)
                 CALL INTELEMENT(XYPWIND(1,INE_WIND(:,WIND_ELE(IP))),XYPWIND(2,INE_WIND(:,WIND_ELE(IP))),VWND_NARR(INE_WIND(:,WIND_ELE(IP))),XP(IP),YP(IP),Wi,TMP(IP,2),.FALSE.)
                 !WRITE(*,*) IP, XYPWIND(1,INE_WIND(:,WIND_ELE(IP)))
                 !WRITE(*,*) IP, XYPWIND(2,INE_WIND(:,WIND_ELE(IP)))
                 !WRITE(*,*) IP, UWND_NARR(INE_WIND(:,WIND_ELE(IP))) 
                 !WRITE(*,*) IP, VWND_NARR(INE_WIND(:,WIND_ELE(IP)))
               END DO
               DVWIND = (TMP-WINDXY)/SEWI%DELT*MAIN%DELT
#endif
             END IF
             SEWI%TMJD = SEWI%TMJD + DBLE(SEWI%DELT)*SEC2DAY
           END IF
           IF (LWINDSWAN) THEN
             WRITE(3333,*) SEWI%TMJD
             WRITE(3333,*) WINDXY(:,1)
             WRITE(3333,*) WINDXY(:,2)
           END IF
         END IF

#ifdef WWMONLY
         IF ( LSECU .AND. (MAIN%TMJD > SECU%TMJD-1.E-8) .AND. (MAIN%TMJD < SECU%EMJD)) THEN
           CALL CSEVAL( CUR%FHNDL, CUR%FNAME, LINTERCU, LCURFILE, MNP, 2, CURTXY,SECU%DELT, MAIN%DELT, DVCURT)
           SECU%TMJD = SECU%TMJD + DBLE(SECU%DELT)*SEC2DAY
           LCALC = .TRUE.
         END IF
         IF ( LSEWL .AND. (MAIN%TMJD > SEWL%TMJD-1.E-8) .AND. (MAIN%TMJD < SEWL%EMJD)) THEN
           CALL CSEVAL( WAT%FHNDL, WAT%FNAME, LINTERWL, LWATLFILE, MNP, 1, WATLEV, SEWL%DELT, MAIN%DELT, DVWALV )
           SEWL%TMJD = SEWL%TMJD + DBLE(SEWL%DELT)*SEC2DAY
           LCALC = .TRUE.
         END IF

         IF (LBCSE .AND. (LBCWA .OR. LBCSP)) THEN 
           IF ( MAIN%TMJD > SEBO%TMJD-1.E-8 .AND. MAIN%TMJD < SEBO%EMJD ) THEN
             IF (IBOUNDFORMAT == 3) THEN ! Find the right position in the file ...
               DTMP = (MAIN%TMJD-BND_TIME_ALL_FILES(1,1)) * DAY2SEC
               !WRITE(*,*) DTMP, MAIN%BMJD, BND_TIME_ALL_FILES(1,1)
               ITMP  = 0
               DO IFILE = 1, NUM_NETCDF_FILES_BND
                 ITMP = ITMP + NDT_BND_FILE(IFILE)
                 IF (ITMP .GT. INT(DTMP/SEBO%DELT)) EXIT
               END DO
               ITMP = SUM(NDT_BND_FILE(1:IFILE-1))
               IT   = NINT(DTMP/SEBO%DELT) - ITMP + 1
               IF (LBINTER) IT = IT + 1
               IF (IT .GT. NDT_BND_FILE(IFILE)) THEN
                 IFILE = IFILE + 1
                 IT    = 1
               ENDIF
             END IF ! BOUNDFORMAT
             !WRITE(*,*) IFILE, IT, SUM(NDT_BND_FILE(1:IFILE-1)), NINT(DTMP/SEBO%DELT), SEBO%DELT
             IF (LBINTER) THEN
               IF (IBOUNDFORMAT == 3) THEN
                 CALL WAVE_BOUNDARY_CONDITION(IFILE,IT,WBACNEW)
               ELSE
                 CALL WAVE_BOUNDARY_CONDITION(1,1,WBACNEW)
               END IF
               DSPEC   = (WBACNEW-WBACOLD)/SEBO%DELT*MAIN%DELT
               WBAC    =  WBACOLD
               WBACOLD =  WBACNEW
               !WRITE(*,*) SUM(WBAC), SUM(WBACNEW), SUM(WBACOLD)
             ELSE ! .NOT. LBINTER
               IF (IBOUNDFORMAT == 3) THEN
                 CALL WAVE_BOUNDARY_CONDITION(IFILE,IT,WBAC)
               ELSE
                 CALL WAVE_BOUNDARY_CONDITION(1,1,WBAC)
               END IF
             END IF
             SEBO%TMJD = SEBO%TMJD + DBLE(SEBO%DELT*SEC2DAY) ! Increment boundary time line ...

           ELSE ! Interpolate in time ...

             IF (LBINTER) THEN
               WBAC = WBAC + DSPEC
             END IF

           END IF

           IF (LINHOM) THEN
             DO IP = 1, IWBMNP
               IPGL = IWBNDLC(IP)
               AC2(IPGL,:,:) = DBLE(WBAC(:,:,IP))
             END DO
           ELSE
             DO IP = 1, IWBMNP
               IPGL = IWBNDLC(IP)
               AC2(IPGL,:,:) = DBLE(WBAC(:,:,1))
             END DO
           ENDIF

         END IF ! LBCWA .OR. LBCSP 

         IF (LSEWL .OR. LSECU .AND. LCALC) THEN
           WATLEVOLD = WATLEV
         END IF
#endif 

         IF (LSEWD .AND. LWINDFROMWWM) THEN
           WINDXY = WINDXY + DVWIND
           WINDXY = WINDXY * WINDFAC
         END IF
!
!      *** coupling via pipe *** read pipe
!
#ifdef WWMONLY
         IF (LCPL .AND. LTIMOR) THEN
           CALL PIPE_TIMOR_IN(K)
           LCALC = .TRUE.
         ELSE IF (LCPL .AND. LSHYFEM) THEN
           CALL PIPE_SHYFEM_IN(K)
           LCALC = .TRUE.
         ELSE IF (LCPL .AND. LROMS) THEN
           CALL PIPE_ROMS_IN(K,IFILE,IT)
           LCALC = .TRUE.
         END IF
#endif 
!
!      *** recalculate water level and current related values 
!
         IF (LSEWL .OR. LSECU .OR. LCPL) THEN ! LCPL makes sure that when the model is coupled it gets into this part for 100%
           IF (.NOT. LCPL) THEN
             WATLEV = WATLEVOLD + DVWALV
             DEP    = MAX(0.,WLDEP + WATLEV) ! d = a + h  if -h .gt. a set d to zero
             DEPDT  = DVWALV / MAIN%DELT
             CURTXY = CURTXY + DVCURT
           ELSE
             DEPDT(:) = ( WATLEV(:) - WATLEVOLD(:) ) / MAIN%DELT
           END IF
           DEP  = MAX(0.,WLDEP + WATLEV) ! d = a + h  if -h .gt. a set d to zero
           CALL GRADDEP
           CALL WAVE_K_C_CG
           CALL GRADCURT
           CALL SET_IOBPD
           CALL SET_DEP
           IF (LCFL) CFLCXY = 0.
           CALL SETSHALLOW
           IF (LMAXETOT .AND. MESBR == 0) CALL SET_HMAX
         END IF

!         IF (MAIN%TMJD .GT. 51545.) WINDXY = 0.

        END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE IO_2(K)

         USE DATAPOOL
         IMPLICIT NONE
!
         INTEGER, INTENT(IN) :: K
!
#ifdef WWMONLY
!
         IF ( (MAIN%TMJD .GE. OUTF%TMJD-1.E-8) .AND. (MAIN%TMJD .LE. OUTF%EMJD)) THEN
           WRITE(STAT%FHNDL,'("+TRACE...",A,4F15.4)') 'WRITING OUTPUT INTERNAL TIME', RTIME, MAIN%TMJD, OUTF%TMJD-1.E-8, OUTF%EMJD
           CALL OUTPUT(REAL(RTIME*DAY2SEC),.FALSE.)
           OUTF%TMJD = OUTF%TMJD + OUTF%DELT*DBLE(SEC2DAY)
         END IF
!
         IF (LHOTF) THEN
           IF ( (MAIN%TMJD .GE. HOTF%TMJD-1.E-8) .AND. (MAIN%TMJD .LE. HOTF%EMJD)) THEN
             WRITE(STAT%FHNDL,'("+TRACE...",A,F15.4)') 'WRITING HOTFILE INTERNAL TIME', RTIME
#ifdef NCDF
             CALL OUTPUT_HOTFILE_NETCDF()
#endif
             HOTF%TMJD = HOTF%TMJD + HOTF%DELT*DBLE(SEC2DAY)
           END IF
         END IF
!
!      *** coupling via pipe *** write pipe
!
         IF (LCPL .AND. LTIMOR) THEN
           CALL PIPE_TIMOR_OUT(K)           
         ELSE IF (LCPL .AND. LSHYFEM) THEN
           CALL PIPE_SHYFEM_OUT(K)
         ELSE IF (LCPL .AND. LROMS) THEN
           CALL PIPE_ROMS_OUT(K)
         END IF
!
#endif

        END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PIPE_TIMOR_IN(K)
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER, INTENT(IN)  :: K

         IF ( K-INT(K/MAIN%ICPLT)*MAIN%ICPLT .EQ. 0 ) THEN
           WRITE(*,'("+TRACE...",A)') 'READING PIPE'
           WRITE(*,'("+TRACE...",A)') 'END READING PIPE'
         END IF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PIPE_SHYFEM_IN(K)
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER, INTENT(IN)  :: K
         INTEGER              :: IP, IL

         IF ( K-INT(K/MAIN%ICPLT)*MAIN%ICPLT .EQ. 0 ) THEN

           WRITE(STAT%FHNDL,'("+TRACE...",A)') 'READING PIPE'
           IF (.NOT. LWINDWWM) THEN
!      *** WIND FROM SHYFEM *** ccf
             DO IP = 1, MNP
               READ(1000) CURTXY(IP,1)
               READ(1001) CURTXY(IP,2)
               READ(1002) WATLEV(IP)
               READ(1003) WLDEP(IP)
               READ(1003) NLEV(IP)
               DO IL = 1,NLVT
                 READ(1004) SHYFZETA(IL,IP)
               END DO
               READ(112) WINDXY(IP,1)
               READ(113) WINDXY(IP,2)
             END DO
            ELSE
!      *** WIND FROM WWM *** ccf
             DO IP = 1, MNP
               READ(1000) CURTXY(IP,1)
               READ(1001) CURTXY(IP,2)
               READ(1002) WATLEV(IP)
               READ(1003) WLDEP(IP)
               READ(1003) NLEV(IP)
               DO IL = 1,NLVT
                 READ(1004) SHYFZETA(IL,IP)
               END DO
             END DO
             WRITE(STAT%FHNDL,*) 'CHECK MAX UX,UY,H'
             WRITE(STAT%FHNDL,*) MAXVAL(CURTXY(:,1)), MAXVAL(CURTXY(:,2)), MAXVAL(WATLEV)
             WRITE(STAT%FHNDL,*) 'CHECK MIN UX,UY,H'
             WRITE(STAT%FHNDL,*) MINVAL(CURTXY(:,1)), MINVAL(CURTXY(:,2)), MINVAL(WATLEV)
             WRITE(2001,*) 'CHECK MAX UX,UY,H'
             WRITE(2001,*) MAXVAL(CURTXY(:,1)), MAXVAL(CURTXY(:,2)), MAXVAL(WATLEV)
             WRITE(2001,*) 'CHECK MIN UX,UY,H'
             WRITE(2001,*) MINVAL(CURTXY(:,1)), MINVAL(CURTXY(:,2)), MINVAL(WATLEV)
           endif
           WRITE(STAT%FHNDL,'("+TRACE...",A)') 'END READ PIPE WWM'
           LCALC = .TRUE.

         END IF

      END SUBROUTINE

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PIPE_ROMS_IN(K,IFILE,IT)
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER, INTENT(IN)  :: K,IFILE,IT
         INTEGER              :: IP

         IF ( K-INT(K/MAIN%ICPLT)*MAIN%ICPLT .EQ. 0 ) THEN
           WRITE(*,'("+TRACE...",A)') 'READING PIPE'
           DO IP = 1, MNP
             READ(1000) WINDXY(IP,1), WINDXY(IP,2), CURTXY(IP,1), CURTXY(IP,2), WATLEV(IP)                   
           END DO
           WRITE(*,'("+TRACE...",A)') 'END READING PIPE'
         END IF
!
!2do make a initialization section for ROMS and WWM
!
         IF (K == 1) CALL INITIAL_CONDITION(IFILE,IT)

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PIPE_TIMOR_OUT(K)
         IMPLICIT NONE

         INTEGER, INTENT(IN)  :: K

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PIPE_ROMS_OUT(K)
         USE DATAPOOL
         IMPLICIT NONE
!
         INTEGER, INTENT(IN)  :: K 

         INTEGER              :: IP
         REAL                 :: TheDISSIP
         REAL                 :: ACLOC(MSC,MDC)
         REAL                 :: HS,WLM,KMWAM2,LPP,FPP,CPP,BOTEXPER
         REAL                 :: SME01,SME10,KME01,KMWAM,URSELL,UBOT,TM,TM01
         REAL                 :: TMBOT,ETOT,FP,CP,KPP,DM,DSPR,ORBITAL,ETOTS,ETOTC,WNPP,TPP,CGPP,TPPD,KPPD,CGPD,CPPD
         REAL                 :: PEAKDSPR, PEAKDM, KME, SME, UR, HSWE, HSLIM, TM02, KLM, DPEAK

         TheDISSIP=10
         Print *, 'export WWM: begining of writing data'
         IF ( K-INT(K/MAIN%ICPLT)*MAIN%ICPLT .EQ. 0 ) THEN
           DO IP = 1, MNP
             ACLOC = AC2(IP,:,:)
             CALL MEAN_PARAMETER(IP,ACLOC,MSC,HS,TM01,TM02,KLM,WLM)
             CALL WAVE_CURRENT_PARAMETER(IP,ACLOC,UBOT,ORBITAL,BOTEXPER,TMBOT)
             CALL MEAN_DIRECTION_AND_SPREAD(IP,ACLOC,MSC,ETOTS,ETOTC,DM,DSPR)
             CALL PEAK_PARAMETER(IP,ACLOC,MSC,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD)
! HS, HSWE, HSLIM  ! - Significant wave height (m) -- HS
! DM  ! - Wave direction (degrees)
! TPP ! - Surface wave relative peak period (s) -- TP
! WLM, KME, SME ! - Average Wave Length [m] - LME
! ORBITAL(IP) ! - Wave bottom orbital velocity (m/s)
! TMBOT ! - Bottom wave period (s)
! TheDissip ! - Wave energy dissipation (W/m2)
! QBLOCAL(IP) ! - Percent of breakig waves (nondimensional)
! DSPR ! - directional spreading
! PEAKDSPR ! - peak directional spreading
! PEAKDM ! - Peak direction
!n
!AR: what for you need this HSWE and HSLIM and so on ... what is exactly SME in your definition ...
!AR: I have deleted them ...
!
             WRITE(101)  HS, DM, TPP, WLM, KLM, TM01, ORBITAL, TMBOT, &
      &                  TheDissip, QBLOCAL(IP), DSPR, PEAKDSPR, PEAKDM, TM02
             CALL FLUSH(101)
           END DO
         END IF
         WRITE(*,*) 'export WWM: ending of writing data'
      END SUBROUTINE

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PIPE_SHYFEM_OUT(K)
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER, INTENT(IN)  :: K

         INTEGER              :: IL, IP
         REAL                 :: ACLOC(MSC,MDC)
         REAL                 :: HS,WLM,LPP
         REAL                 :: SME01,SME10,KME01,KMWAM,KMWAM2,URSELL,UBOT
         REAL                 :: TMBOT,ETOT,FP,CP,KPP,DM,DSPR,PEAKDSPR,PEAKDM
         REAL                 :: ORBITAL, USTOKES, USTOKES_X, USTOKES_Y
         REAL                 :: BOTEXPER, ETOTS, ETOTC, DPEAK
         REAL                 :: FPP, TPP, CPP, WNPP, CGPP, TPPD,KPPD,CGPD,CPPD 

           IF ( K-INT(K/MAIN%ICPLT)*MAIN%ICPLT .EQ. 0 ) THEN

             WRITE(STAT%FHNDL,'("+TRACE...",A)') 'WRITING PIPE'
!
!            *** COMPUTE RADIATION STRESSES 2D OR 3D *** ccf
!
             IF (LCPL) THEN
               CALL RADIATION_STRESS
             END IF

             IF (.NOT. LWINDWWM) THEN
               DO IP = 1, MNP
                 DO IL = 1, NLVT
                   WRITE(101)  SXX3D(IL,IP)             !ccf
                   WRITE(102)  SYY3D(IL,IP)             !ccf
                   WRITE(142)  SXY3D(IL,IP)             !ccf
                 END DO
                 CALL FLUSH(101)
                 CALL FLUSH(102)
                 CALL FLUSH(142)                        !ccf
                 ACLOC = AC2(IP,:,:)
!2do check ustokes ...
                 CALL MEAN_WAVE_PARAMETER(IP,ACLOC,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2)
                 CALL WAVE_CURRENT_PARAMETER(IP,ACLOC,UBOT,ORBITAL,BOTEXPER,TMBOT)
                 CALL MEAN_DIRECTION_AND_SPREAD(IP,ACLOC,MSC,ETOTS,ETOTC,DM,DSPR)
                 CALL PEAK_PARAMETER(IP,ACLOC,MSC,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD)

                 WRITE(103) HS
                 CALL FLUSH(103)
                 WRITE(104)  SME01
                 CALL FLUSH(104)
                 WRITE(105)  DM
                 CALL FLUSH(105)
                 WRITE(106)  KMWAM
                 CALL FLUSH(106)
                 WRITE(107)  TPP
                 CALL FLUSH(107)
                 WRITE(108)  KPP
                 CALL FLUSH(108)
                 WRITE(109)  ORBITAL
                 CALL FLUSH(109)
!AR: Delete this stuff 
                 !WRITE(110)  USTOKES_X
                 CALL FLUSH(110)
                 !WRITE(111)  USTOKES_Y
                 CALL FLUSH(111)
               END DO
             ELSE
               DO IP = 1, MNP
                 DO IL = 1, NLVT
                   WRITE(101)  SXX3D(IL,IP)             !ccf
                   WRITE(102)  SYY3D(IL,IP)             !ccf
                   WRITE(142)  SXY3D(IL,IP)             !ccf
                 END DO
                 CALL FLUSH(101)
                 CALL FLUSH(102)
                 CALL FLUSH(142)                        !ccf
                 CALL MEAN_WAVE_PARAMETER(IP,ACLOC,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2)
                 CALL WAVE_CURRENT_PARAMETER(IP,ACLOC,UBOT,ORBITAL,BOTEXPER,TMBOT)
                 CALL MEAN_DIRECTION_AND_SPREAD(IP,ACLOC,MSC,ETOTS,ETOTC,DM,DSPR)
                 CALL PEAK_PARAMETER(IP,ACLOC,MSC,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD)
                 WRITE(103) HS
                 CALL FLUSH(103)
                 WRITE(104)  SME01
                 CALL FLUSH(104)
                 WRITE(105)  DM
                 CALL FLUSH(105)
                 WRITE(106)  KMWAM
                 CALL FLUSH(106)
                 WRITE(107)  TPP
                 CALL FLUSH(107)
                 WRITE(108)  KPP
                 CALL FLUSH(108)
                 WRITE(109)  ORBITAL
                 CALL FLUSH(109)
!AR: Delete this stuff 
                 !WRITE(110)  USTOKES_X
                 CALL FLUSH(110)
                 !WRITE(111)  USTOKES_Y
                 CALL FLUSH(111)
                 WRITE(112)  WINDXY(IP,1)
                 CALL FLUSH(112)
                 WRITE(113)  WINDXY(IP,2)
                 CALL FLUSH(113)
               END DO
             END IF
             WRITE(STAT%FHNDL,'("+TRACE...",A)') 'END WRITING PIPE'
          END IF

      END SUBROUTINE

!**********************************************************************
!*                                                                    *
!**********************************************************************
