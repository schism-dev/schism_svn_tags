#include "wwm_functions.h"
!:
! __      __  __      __  _____  .___.___. 
!/  \    /  \/  \    /  \/     \ |   |   |  
!\   \/\/   /\   \/\/   /  \ /  \|   |   |   
! \        /  \        /    Y    \   |   |   
!  \__/\  /    \__/\  /\____|__  /___|___|
!       \/          \/         \/             
!
! WWM-II (Wind Wave Model) source code 
! 
! The source code is entirely rewritten with respect to WWM (Hsu et al., 2005). The numerics have been completely revised (Roland, 2008)
! The code included various source term packages (see Manual) and can be coupled to various ocean models on structured and unstructured grids
! Parallelization is done using OpenMP or MPI. For coupling to certain models we are using either Pipes (Roland et al. 2009), 
! coupling libraries (PGMCL, Dutour-Sikiric et al. 2013) or tightly coupled with SELFE (Roland et al. 2012).  
! 
! Developers:                                                   
! Lead: Aron Roland (IT&E, Frankfurt, Z&P, Hannover), Yinglong Joseph Zhang (VIMS), Mathieu-Dutour Sikiric (IRB, Zagreb), Ulrich Zanke (Z&P, Hannover) 
!
! Contributors: Christian Ferrarin (ISMAR-CNR), Fabrice Ardhuin (IFREMER), Yaron Toledo (Tel-Aviv University), Thomas Huxhorn (TUD), Ivica Janekovic (IRB, Zagreb), 
! Will Perrie (Fisheries, Canada), Bash Toulany (Fisheries, Canada), Harry Wang (VIMS), Andrea Fortunato (LNEC), Guillaume Dodet (LNEC), Kai Li (XXXX)
! Andreas Wurpts (Forschungsstelle Küste, Norderney), Michael Glozman (Cameri, Technion), Boris Lehmann (TUD, Darmstadt) 
!				
! Copyright: 2008 - XXXX Z&P (Aron Roland, IT&E, Frankfurt, Zanke&Partner, Hannover, Germany)
! All Rights Reserved                                     
!
! License: Redistribution of any files contained in this package is strictly prohibited
! Any kind usage only allowed only with permission of Zanke & Partner (aaronroland@gmx.de)
! This includes commerical as well as academic usage. Developers and Contributers 
! are not subject to this licese condition and can use this code as they wish. 
! For any kind of questions or licence inquries please contact: aaronroland@gmx.de
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef SELFE
 !!!     SUBROUTINE WWM_II(IT_SELFE,icou_elfe_wwm,DT_SELFE0,NSTEP_WWM0)
      SUBROUTINE WWM_II(IT_SELFE,icou_elfe_wwm,DT_SELFE0,NSTEP_WWM0,HS_DUMON,SBR_X,SBR_Y,STOKES_X,STOKES_Y,JPRESS,RADFLAG2)   !! modif AD



         USE DATAPOOL
         use elfe_msgp!, only : myrank,parallel_abort,itype,comm,ierr
         use elfe_glbl, only : iplg,ielg

         IMPLICIT NONE

         INTEGER, INTENT(IN)   :: NSTEP_WWM0, icou_elfe_wwm
         REAL(rkind), INTENT(IN)    :: DT_SELFE0
         CHARACTER(LEN=3), INTENT(OUT) :: RADFLAG2

         REAL(rkind), SAVE  :: SIMUTIME
         REAL(rkind)        :: T1, T2, mpi_wtime
         REAL(rkind)        :: TIME1, TIME2, TIME3, TIME4, TIME5, TIME6, TIME7

         INTEGER     :: I, IP, IT_SELFE, K, IFILE, IT, IPGL
         REAL(rkind) :: DT_PROVIDED
         REAL(rkind) :: OUTPAR(OUTVARS), OUTWINDPAR(WINDVARS), ACLOC(MSC,MDC)
         character(LEN=15) :: CALLFROM


!!!!!!!!!!!!!!!!!! modif AD
          REAL(rkind) :: SSBR_DUMON_ALL(MNP,MSC,MDC),QB_DUMON_ALL(MNP)
          REAL(rkind) :: SBR_X(MNP),SBR_Y(MNP)
          REAL(rkind) :: STOKES_X(NVRT,MNP),STOKES_Y(NVRT,MNP),JPRESS(MNP)
          REAL(rkind) :: HS_DUMON(MNP),ETOT_DUMON(MNP)
!!!!!!!!!!!!!!!!!!  end modif AD
      


         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING WWM_II'
         CALL FLUSH(STAT%FHNDL)

         TIME1 = mpi_wtime()

         IF (LNANINFCHK) THEN
           WRITE(DBG%FHNDL,*) ' STARTING WWM FROM SELFE ',  SUM(AC2)
           IF (SUM(AC2) .NE. SUM(AC2)) call wwm_abort('NAN IN MAIN 1')
         ENDIF

         if(WINDVARS/=size(WIND_INTPAR,2)) call wwm_abort('Dimension mismatch: OUTWINDPAR and out_wwm_windpar')
         if(OUTVARS/=size(OUTT_INTPAR,2))  call wwm_abort('Dimension mismatch: OUTPAR and out_wwm')

         NSTEPWWM = NSTEP_WWM0

         DT_SELFE      = DT_SELFE0
         DELTAT_WATLEV = DT_SELFE0

         T1 = MyREAL(IT_SELFE-NSTEPWWM)*DT_SELFE0 ! Beginn time step ...
         T2 = MyREAL(IT_SELFE)*DT_SELFE0          ! End of time time step ...
         DT_PROVIDED=NSTEPWWM*DT_SELFE

         IF (abs(MAIN%DELT - DT_PROVIDED).gt.THR) THEN
           WRITE(DBG%FHNDL,*) 'MAIN%DELT=', MAIN%DELT, ' in wwminput.nml'
           WRITE(DBG%FHNDL,*) 'But nstep_wwm*dt=', DT_PROVIDED
           WRITE(DBG%FHNDL,*) 'nstep_wwm=', NSTEPWWM
           WRITE(DBG%FHNDL,*) '       dt=', DT_SELFE
           CALL WWM_ABORT('Correct coupled model time-steppings')
         ENDIF

         IF (LNANINFCHK) THEN
           WRITE(DBG%FHNDL,*) ' FIRST SUM IN MAIN ',  SUM(AC2)
           IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN MAIN 2')
         ENDIF

         WRITE(STAT%FHNDL,'("+TRACE...",A)') ' ---- ALL CHECKS DONE'
         CALL FLUSH(STAT%FHNDL)

         SIMUTIME = SIMUTIME + MAIN%DELT

         IF (icou_elfe_wwm == 1) THEN ! Full coupling 
           WLDEP       = DEP8
           WATLEV      = ETA2
           WATLEVOLD   = ETA1
           DEP         = MAX(ZERO,WLDEP + WATLEV)
           CURTXY(:,1) = UU2(NVRT,:)
           CURTXY(:,2) = VV2(NVRT,:)
           IF (.NOT. LWINDFROMWWM) THEN
             WINDXY(:,1) = WINDX0
             WINDXY(:,2) = WINDY0
           END IF
           LSECU       = .TRUE.
           LSEWL       = .TRUE.
           LCALC       = .TRUE.
         ELSE IF (icou_elfe_wwm == 0) THEN ! No interaction at all 
           WLDEP       = DEP8
           WATLEV      = ZERO 
           WATLEVOLD   = ZERO
           DEP         = WLDEP
           CURTXY(:,1) = ZERO !REAL(rkind)(UU2(NVRT,:))
           CURTXY(:,2) = ZERO !REAL(rkind)(VV2(NVRT,:))
           IF (.NOT. LWINDFROMWWM) THEN
             WINDXY(:,1) = WINDX0
             WINDXY(:,2) = WINDY0
           END IF
           LSECU       = .FALSE.
           LSEWL       = .FALSE.
           LCALC       = .TRUE. 
         ELSE IF (icou_elfe_wwm == 2) THEN ! Currents and water levels in wwm but no radiation stress in SELFE 
           WLDEP       = DEP8
           WATLEV      = ETA2
           WATLEVOLD   = ETA1
           DEP         = MAX(ZERO, WLDEP + WATLEV)
           CURTXY(:,1) = UU2(NVRT,:)
           CURTXY(:,2) = VV2(NVRT,:)
           IF (.NOT. LWINDFROMWWM) THEN
             WINDXY(:,1) = WINDX0
             WINDXY(:,2) = WINDY0
           END IF
           LSECU       = .TRUE.
           LSEWL       = .TRUE.
           LCALC       = .TRUE.
         ELSE IF (icou_elfe_wwm == 3) THEN ! No current and no water levels in wwm but radiation stress in SELFE 
           WLDEP       = DEP8
           WATLEV      = ZERO
           WATLEVOLD   = ZERO
           DEP         = WLDEP
           CURTXY(:,1) = ZERO !REAL(rkind)(UU2(NVRT,:))
           CURTXY(:,2) = ZERO !REAL(rkind)(VV2(NVRT,:))
           IF (.NOT. LWINDFROMWWM) THEN
             WINDXY(:,1) = WINDX0
             WINDXY(:,2) = WINDY0
           END IF
           LSECU       = .FALSE.
           LSEWL       = .FALSE.
           LCALC       = .TRUE.
         ELSE IF (icou_elfe_wwm == 4) THEN ! No current but water levels in wwm and radiation stresss in selfe
           WLDEP       = DEP8
           WATLEV      = ETA2
           WATLEVOLD   = ETA1
           DEP         = WLDEP
           CURTXY(:,1) = 0.!UU2(NVRT,:) 
           CURTXY(:,2) = 0.!UU2(NVRT,:) 
           IF (.NOT. LWINDFROMWWM) THEN
             WINDXY(:,1) = WINDX0
             WINDXY(:,2) = WINDY0
           END IF
           LSECU       = .FALSE.
           LSEWL       = .TRUE.
           LCALC       = .TRUE.
         ELSE IF (icou_elfe_wwm == 5) THEN ! No current but water levels in wwm and no radiation stress in selfe  
           WLDEP       = DEP
           WATLEV      = ETA2
           WATLEVOLD   = ETA1
           DEP         = WLDEP
           CURTXY(:,1) = 0.!UU2(NVRT,:) 
           CURTXY(:,2) = 0.!UU2(NVRT,:) 
           IF (.NOT. LWINDFROMWWM) THEN
             WINDXY(:,1) = WINDX0
             WINDXY(:,2) = WINDY0
           END IF
           LSECU       = .FALSE.
           LSEWL       = .TRUE.
           LCALC       = .TRUE.
         ELSE IF (icou_elfe_wwm == 6) THEN ! Currents but no water levels in wwm and radiation stress in selfe  
           WLDEP       = DEP
           WATLEV      = ZERO 
           WATLEVOLD   = ZERO 
           DEP         = WLDEP
           CURTXY(:,1) = UU2(NVRT,:) 
           CURTXY(:,2) = VV2(NVRT,:) 
           IF (.NOT. LWINDFROMWWM) THEN
             WINDXY(:,1) = WINDX0
             WINDXY(:,2) = WINDY0
           END IF
           LSECU       = .TRUE.
           LSEWL       = .FALSE.
           LCALC       = .TRUE.
         ELSE IF (icou_elfe_wwm == 7) THEN ! Currents but no water levels in wwm and no radiation stress in selfe  
           WLDEP       = DEP
           WATLEV      = ZERO
           WATLEVOLD   = ZERO
           DEP         = WLDEP
           CURTXY(:,1) = UU2(NVRT,:)
           CURTXY(:,2) = UU2(NVRT,:)
           IF (.NOT. LWINDFROMWWM) THEN
             WINDXY(:,1) = WINDX0
             WINDXY(:,2) = WINDY0
           END IF
           LSECU       = .TRUE.
           LSEWL       = .FALSE.
           LCALC       = .TRUE.
         END IF

         IF (LNANINFCHK) THEN
           CALL SELFE_NANCHECK_INPUT_A
         END IF

         IFILE = 1
         IT    = 1

         CALL SET_WAVE_BOUNDARY_CONDITION

         IF (LNANINFCHK) THEN
           WRITE(DBG%FHNDL,*) ' AFTER SETTING BOUNDARY CONDITION IN MAIN ',  SUM(AC2)
           IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN MAIN 3')
         ENDIF

         IF (LFIRSTSTEP) THEN
           IF (INITSTYLE == 1) CALL INITIAL_CONDITION(IFILE,IT)!We need to call for the case of wind dependent intiial guess this call since before we have no wind from SELFE
           LFIRSTSTEP = .FALSE.
           LCALC      = .TRUE.
         END IF

         TIME2 = mpi_wtime() 

         IF (LNANINFCHK) THEN
           WRITE(DBG%FHNDL,*) ' BEFORE COMPUTE ',  SUM(AC2)
           IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN MAIN 4')
         ENDIF

         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING COMPUTE'
         CALL FLUSH(STAT%FHNDL)

         CALLFROM='SELFE'
         IF (LQSTEA) THEN
            CALL QUASI_STEADY(KKK)
         ELSE
!!!!! modif AD
            CALL UN_STEADY(KKK,CALLFROM,SSBR_DUMON_ALL,QB_DUMON_ALL)
!!!!! modif AD
         END IF

         IF (LNANINFCHK) THEN
           WRITE(DBG%FHNDL,*) ' AFTER COMPUTE ',  SUM(AC2)
           IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT('NAN IN MAIN 5') 
         ENDIF

         TIME3 = mpi_wtime()

         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.4)') 'FINISHED COMPUTE nth call to WWM', SIMUTIME
         CALL FLUSH(STAT%FHNDL)

         DO IP = 1, MNP
           ACLOC = AC2(IP,:,:)
           write(44444,*) IP, DEP8(IP)
           IF (DEP(IP) .GT. DMIN) THEN
             CALL INTPAR(IP, MSC, ACLOC, OUTPAR)
             OUTT_INTPAR(IP,:) = OUTPAR
             CALL WINDPAR(IP,OUTWINDPAR)
             WIND_INTPAR(IP,:) = OUTWINDPAR
             IF (LMONO_OUT) THEN
               OUTT_INTPAR(IP,1) = OUTT_INTPAR(IP,1) / SQRT(2.)
             END IF
           ELSE
             OUTT_INTPAR(IP,:) = ZERO
             WIND_INTPAR(IP,:) = ZERO
           END IF
         END DO

         TIME4 = mpi_wtime()
!
! Compute radiation stress ...
!

!!!!!! modif AD : if RADFLAG=VOR , then coupling with selfe will gives stokes_velocity (Eq. 17 from Bennis 2011), Wave-induced pressure (Eq. 20) and source momentums (Eq.21) 
         RADFLAG2=RADFLAG !for output into SELFE
         IF (icou_elfe_wwm == 0 .OR. icou_elfe_wwm == 2 .OR. icou_elfe_wwm == 5 .OR. icou_elfe_wwm == 7) THEN
           WWAVE_FORCE = ZERO
           STOKES_X=ZERO
           STOKES_Y=ZERO
           JPRESS=ZERO
           SBR_X=ZERO
           SBR_Y=ZERO
         ELSE 
           IF (RADFLAG == 'VOR') THEN
             CALL STOKES_STRESS_INTEGRAL_SELFE(SSBR_DUMON_ALL,STOKES_X,STOKES_Y,JPRESS,SBR_X,SBR_Y,HS_DUMON,ETOT_DUMON)
           ELSE
             CALL RADIATION_STRESS_SELFE
           ENDIF
         END IF 
!!!!!! end modif AD


         TIME5 = mpi_wtime()
 
         IF (LNANINFCHK) THEN
           CALL SELFE_NANCHECK_INPUT_B
         END IF

         KKK = KKK + 1

         TIME6 = mpi_wtime()

         IF (LNANINFCHK) THEN
           WRITE(DBG%FHNDL,*) ' END OF MAIN ',  SUM(AC2)
           IF (SUM(AC2) .NE. SUM(AC2)) CALL WWM_ABORT ('NAN IN MAIN 5')
         ENDIF

         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-----TOTAL TIMINGS-----'
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'PREPARATION        ', TIME2-TIME1
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'INTEGRATION        ', TIME3-TIME2
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'OUTPUT TO SELFE    ', TIME4-TIME3
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'RADIATION STRESSES ', TIME5-TIME4
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'NAN CHECK          ', TIME6-TIME5
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'TOTAL TIME         ', TIME6-TIME1
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '------END-TIMINGS-  ---'

         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.4)') 'FINISHED WITH WWM', SIMUTIME
         CALL FLUSH(STAT%FHNDL)
 
      END SUBROUTINE WWM_II
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SELFE_NANCHECK_INPUT_A
      USE DATAPOOL
      implicit none
      integer IP
      DO IP = 1, MNP
        IF (WINDXY(IP,1) .NE. WINDXY(IP,1)) THEN
          WRITE(DBG%FHNDL,*) 'NaN in WINDX', IP, WINDXY(IP,1) 
          CALL FLUSH(DBG%FHNDL)
        END IF
        IF (WINDXY(IP,2) .NE. WINDXY(IP,2)) THEN
          WRITE(DBG%FHNDL,*) 'NaN in WINDY', IP, WINDXY(IP,2) 
          CALL FLUSH(DBG%FHNDL)
        END IF
        IF (WATLEV(IP) .NE. WATLEV(IP)) THEN
          WRITE(DBG%FHNDL,*) 'NaN in WATLEV', IP, WATLEV(IP) 
          CALL FLUSH(DBG%FHNDL)
        END IF
        IF (WATLEVOLD(IP) .NE. WATLEVOLD(IP)) THEN
          WRITE(DBG%FHNDL,*) 'NaN in WATLEV', IP, WATLEV(IP)
          CALL FLUSH(DBG%FHNDL)
        END IF
        IF (CURTXY(IP,1) .NE. CURTXY(IP,1)) THEN
          WRITE(DBG%FHNDL,*) 'NaN in CURTX', IP, CURTXY(IP,1)
          CALL FLUSH(DBG%FHNDL)
        END IF
        IF (CURTXY(IP,2) .NE. CURTXY(IP,2)) THEN
          WRITE(DBG%FHNDL,*) 'NaN in CURTY', IP, CURTXY(IP,2)
          CALL FLUSH(DBG%FHNDL)
        END IF
      END DO
      END SUBROUTINE SELFE_NANCHECK_INPUT_A
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SELFE_NANCHECK_INPUT_B
      USE DATAPOOL
      implicit none
      integer IP, I
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
!Error: force defined at side centers
        IF (SUM(WWAVE_FORCE(:,IP,:)) .NE. SUM(WWAVE_FORCE(:,IP,:))) THEN
          DO I = 1, SIZE(WWAVE_FORCE(1,:,IP)) ! loop over layers ...
            WRITE(DBG%FHNDL,*) 'NaN in WWAVE_FORCE', IP, I, WWAVE_FORCE(1,I,IP), WWAVE_FORCE(2,I,IP)
            CALL FLUSH(DBG%FHNDL)
          END DO
        END IF 
      END DO
      END SUBROUTINE SELFE_NANCHECK_INPUT_B
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE UN_STEADY(K,CALLFROM,SSBR_DUMON_ALL,QB_DUMON_ALL)

         USE DATAPOOL
#ifdef MPI_PARALL_GRID
         use elfe_msgp, only : myrank
#endif

         IMPLICIT NONE

!!!!! modif AD
         REAL(rkind), INTENT(OUT)  :: SSBR_DUMON_ALL(MNP,MSC,MDC),QB_DUMON_ALL(MNP)
!!!!! end modif AD

         INTEGER, INTENT(IN) :: K

         REAL(rkind)    :: CONV1, CONV2, CONV3, CONV4, CONV5
         REAL           :: TIME1, TIME2, TIME3, TIME4, TIME5, TIME6

         CHARACTER(LEN=15)   :: CTIME,CALLFROM

         CALL CPU_TIME(TIME1)
       
         CALL IO_1(K)

         CALL CPU_TIME(TIME2)

#if !defined MPI_PARALL_GRID
         IF (LCFL) CALL CFLSPEC
#endif

         CALL CPU_TIME(TIME3)

         IF (ICOMP .EQ. 0) THEN
           IF (LITERSPLIT) THEN
            CALL COMPUTE_DOUBLE_STRANG_EXPLICIT
             !CALL COMPUTE_ITERATIVE_SPLITTING
           ELSE
!!!!! modif AD
             CALL COMPUTE_SIMPLE_EXPLICIT(SSBR_DUMON_ALL,QB_DUMON_ALL)
!!!!! end modif AD
           END IF
         ELSE IF (ICOMP .EQ. 1) THEN 
           IF (LITERSPLIT) THEN
            CALL COMPUTE_DOUBLE_STRANG_EXPLICIT
           ELSE
!!!!! modif AD
             CALL COMPUTE_SIMPLE_EXPLICIT(SSBR_DUMON_ALL,QB_DUMON_ALL)
!!!!! end modif AD
            END IF
         ELSE IF (ICOMP .EQ. 2) THEN 
           CALL COMPUTE_IMPLICIT
         END IF

         CALL CPU_TIME(TIME4)

         MAIN%TMJD = MAIN%BMJD + MyREAL(K)*MAIN%DELT*SEC2DAY
         RTIME = MAIN%TMJD - MAIN%BMJD
#ifndef SELFE
# if defined WWM_MPI
         IF (myrank.eq.0) THEN
# endif
         WRITE(*,101)  K, MAIN%ISTP, RTIME
# if defined WWM_MPI
         ENDIF
# endif
#endif
         CALL IO_2(K)

         CALL CPU_TIME(TIME5)

         IF (LCONV) THEN
            CALL CHECK_STEADY(RTIME,CONV1,CONV2,CONV3,CONV4,CONV5)
         END IF

         CALL CPU_TIME(TIME6)

         IF (.NOT. LDIFR) LCALC = .FALSE.

         CALL MJD2CT(MAIN%TMJD, CTIME)

#ifdef MPI_PARALL_GRID
      IF (myrank == 0) THEN
#endif
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6,A20)') '-----SIMULATION TIME-----        ', MAIN%TMJD, CTIME
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-----TOTAL RUN TIMES-----        '
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'PREPROCESSING                    ', TIME2-TIME1
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CHECK CFL                        ', TIME3-TIME2
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'INTEGRATION                      ', TIME4-TIME3
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'POSTPROCESSING                   ', TIME5-TIME4
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'CHECK STEADY                     ', TIME6-TIME5
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') '-------------TIMINGS-------------'
        CALL FLUSH(STAT%FHNDL)
#ifdef MPI_PARALL_GRID
      ENDIF
#endif

101      FORMAT ('+STEP = ',I10,'/',I10,' ( TIME = ',F15.4,' DAYS)')

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE QUASI_STEADY(K)

#ifdef MPI_PARALL_GRID
         USE elfe_msgp
#endif
         USE DATAPOOL

         IMPLICIT NONE
         INTEGER, INTENT(IN) :: K

         INTEGER :: IT
         REAL(rkind)    :: ITERTIME
         REAL(rkind)    :: CONV1, CONV2, CONV3, CONV4, CONV5

!!!!! modif AD
         REAL(rkind)  :: SSBR_DUMON_ALL(MNP,MSC,MDC),QB_DUMON_ALL(MNP)
!!!!! end modif AD



         CALL IO_1(K)

         IF (LCFL) CALL CFLSPEC()

         IF (LCHKCONV) IP_IS_STEADY = 0 ! Reset local convergence indicators ...
         IF (LCHKCONV) IE_IS_STEADY = 0

#ifdef MPI_PARALL_GRID
!         NQSITER = NSTEPWWM ! this is not very flexible!
#endif

         DO IT = 1, NQSITER

           DT_ITER = MAIN%DELT/MyREAL(NQSITER)

           IF (ICOMP .LT. 2) THEN
 !!!!! modif AD
             CALL COMPUTE_SIMPLE_EXPLICIT(SSBR_DUMON_ALL,QB_DUMON_ALL)
!!!!! end modif AD
           ELSE
             CALL COMPUTE_IMPLICIT
           END IF

           ITERTIME = RTIME*DAY2SEC+IT*DT_ITER

           IF (LCHKCONV) THEN
             CALL CHECK_STEADY(ITERTIME,CONV1,CONV2,CONV3,CONV4,CONV5)
!             DO IP = 1, MNP
!               IF (IP_IS_STEADY(IP) .GE. 1) AC2(IP,:,:) = AC1(IP,:,:)
!             ENDDO
             IF ( (CONV1 .GT. 100._rkind*QSCONV1 .AND.                  &
     &             CONV2 .GT. 100._rkind*QSCONV2 .AND.                    &
     &             CONV3 .GT. 100._rkind*QSCONV3 .AND.                    &
     &             CONV4 .GT. 100._rkind*QSCONV4 .AND.                    &
     &             CONV5 .GT. 100._rkind*QSCONV5 .AND.                    &
     &             K .NE. 1) .OR.                                         &
     &             IT .EQ. NQSITER ) THEN
#ifndef SELFE
               WRITE(QSTEA%FHNDL,'(3I10,5F15.8)') K, IT, NQSITER, CONV1, CONV2, CONV3, CONV4, CONV5
#elif SELFE
               if (myrank == 0) WRITE(QSTEA%FHNDL,'(3I10,5F15.8)') K, IT, NQSITER, CONV1, CONV2, CONV3, CONV4, CONV5
#endif
               CALL FLUSH(QSTEA%FHNDL)
               EXIT 
             END IF
           END IF
           IF (LOUTITER) CALL WWM_OUTPUT(ITERTIME,.FALSE.)
         END DO

#ifdef MPI_PARALL_GRID
         MAIN%TMJD = MAIN%BMJD + MyREAL(K)*MAIN%DELT*SEC2DAY
         RTIME = MAIN%TMJD - MAIN%BMJD
         IF (myrank == 0) WRITE(STAT%FHNDL,101)  K, MAIN%ISTP, RTIME*DAY2SEC
#else
         MAIN%TMJD = MAIN%BMJD + MyREAL(K)*MAIN%DELT*SEC2DAY
         RTIME = MAIN%TMJD - MAIN%BMJD
         WRITE(STAT%FHNDL,101)  K, MAIN%ISTP, RTIME*DAY2SEC
#endif
         CALL FLUSH(STAT%FHNDL)

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
#ifdef MPI_PARALL_GRID
         USE elfe_msgp
#endif
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: K
         REAL(rkind)  :: TMP_CUR(MNP,2), TMP_WAT(MNP)
         INTEGER             :: IT, IFILE

! update wind ...
         IF (LWINDFROMWWM) THEN
           CALL UPDATE_WIND(K)
         END IF

#ifndef SELFE
         IF ( LSECU .AND. (MAIN%TMJD > SECU%TMJD-1.E-8) .AND. (MAIN%TMJD < SECU%EMJD)) THEN
           CALL CSEVAL( CUR%FHNDL, CUR%FNAME, .TRUE., 2, TMP_CUR)
           DVCURT=(TMP_CUR - CURTXY)/SECU%DELT*MAIN%DELT
           SECU%TMJD = SECU%TMJD + SECU%DELT*SEC2DAY
           LCALC = .TRUE.
         END IF
         IF ( LSEWL .AND. (MAIN%TMJD > SEWL%TMJD-1.E-8) .AND. (MAIN%TMJD < SEWL%EMJD)) THEN
           CALL CSEVAL( WAT%FHNDL, WAT%FNAME, .TRUE., 1, TMP_WAT)
           DVWALV=(TMP_WAT - WATLEV)/SEWL%DELT*MAIN%DELT
           SEWL%TMJD = SEWL%TMJD + SEWL%DELT*SEC2DAY
           LCALC = .TRUE.
         END IF

         CALL SET_WAVE_BOUNDARY_CONDITION

         IF (LSEWL .OR. LSECU .AND. LCALC) THEN
           DELTAT_WATLEV = MAIN%DELT
           WATLEVOLD = WATLEV
         END IF
#endif
!
!      *** coupling via pipe *** read pipe
!
#if !defined SELFE && !defined PGMCL_COUPLING
         IF (LCPL .AND. LTIMOR) THEN
           CALL PIPE_TIMOR_IN(K)
           LCALC = .TRUE.
# ifdef SHYFEM_COUPLING
         ELSE IF (LCPL .AND. LSHYFEM) THEN
           CALL PIPE_SHYFEM_IN(K)
           LCALC = .TRUE.
# endif
         ELSE IF (LCPL .AND. LROMS) THEN
           CALL PIPE_ROMS_IN(K,IFILE,IT)
           LCALC = .TRUE.
         END IF
#endif
#ifdef PGMCL_COUPLING
         IF ( K-INT(K/MAIN%ICPLT)*MAIN%ICPLT .EQ. 0 ) THEN
           CALL PGMCL_ROMS_IN(K,IFILE,IT)
# ifdef ANALYTICAL_WIND_CURR
           DO IP=1,MNP
             WINDXY(IP,1)=ZERO
             WINDXY(IP,2)=XP(IP)/5.0_rkind
             CURTXY(IP,1)=XP(IP)/30.0_rkind
             CURTXY(IP,2)=ZERO
             WATLEV(IP)=ZERO
             WATLEVOLD(IP)=ZERO
           END DO
# endif
         END IF
         IF (K == 1) CALL INITIAL_CONDITION(IFILE,IT)
#endif
!
!      *** recalculate water level and current related values 
!
         IF (LSEWL .OR. LSECU .OR. LCPL) THEN ! LCPL makes sure that when the model is coupled it gets into this part for 100%
           IF (.NOT. LCPL) THEN
             WATLEV = WATLEVOLD + DVWALV
             DEP    = MAX(ZERO,WLDEP + WATLEV) ! d = a + h  if -h .gt. a set d to zero
             DEPDT  = DVWALV / MAIN%DELT
             CURTXY = CURTXY + DVCURT
           ELSE
             DEPDT(:) = ( WATLEV(:) - WATLEVOLD(:) ) / DELTAT_WATLEV
           END IF
           DEP  = MAX(ZERO,WLDEP + WATLEV) ! d = a + h  if -h .gt. a set d to zero
           CALL SETSHALLOW
           CALL GRADDEP
           CALL GRAD_CG_K
           CALL WAVE_K_C_CG
           CALL GRADCURT
           CALL SET_IOBPD
           CALL SET_IOBPD_BY_DEP
           IF (LCFL) CFLCXY = ZERO
           IF (LMAXETOT .AND. MESBR == 0) CALL SET_HMAX
         END IF

!         IF (MAIN%TMJD .GT. 51545.) WINDXY = ZERO

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE IO_2(K)
         USE WWM_HOTFILE_MOD
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: K

         CALL OUTPUT_HISTORY_AND_STATION

         IF (LHOTF) THEN
           IF ( (MAIN%TMJD .GE. HOTF%TMJD-1.E-8) .AND. (MAIN%TMJD .LE. HOTF%EMJD)) THEN
             WRITE(STAT%FHNDL,'("+TRACE...",A,F15.4)') 'WRITING HOTFILE INTERNAL TIME', RTIME
             CALL FLUSH(STAT%FHNDL)
             CALL OUTPUT_HOTFILE
             HOTF%TMJD = HOTF%TMJD + HOTF%DELT*SEC2DAY
           END IF
         END IF
#ifndef SELFE
# if !defined PGMCL_COUPLING
         IF (LCPL .AND. LTIMOR) THEN
           CALL PIPE_TIMOR_OUT(K)
#  ifdef SHYFEM_COUPING
         ELSE IF (LCPL .AND. LSHYFEM) THEN
           CALL PIPE_SHYFEM_OUT(K)
#  endif
         ELSE IF (LCPL .AND. LROMS) THEN
           CALL PIPE_ROMS_OUT(K)
         END IF
# endif
# ifdef PGMCL_COUPLING
         IF ( K-INT(K/MAIN%ICPLT)*MAIN%ICPLT .EQ. 0 ) THEN
           CALL PGMCL_ROMS_OUT(K)
         END IF
# endif
#endif
        END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#if !defined SELFE
      SUBROUTINE SET_WWMINPULNML
        USE DATAPOOL, only : INP
!This is all to special stuff if PGMCL then the input file has other name and all this special shit that makes you the life easier and fucks up the code and nobody can reproduce it. 
# ifndef PGMCL_COUPLING
        IMPLICIT NONE
        INTEGER nbArg
        nbArg=command_argument_count()
        IF (nbArg > 1) THEN
          CALL WWM_ABORT('Number of argument is 0 or 1')
        ENDIF
        IF (nbArg.eq.0) THEN
          INP%FNAME  = 'wwminput.nml'
        ELSE
          CALL GET_COMMAND_ARGUMENT(1, INP%FNAME)
        ENDIF
# else
        USE mod_coupler, only : Iwaves, INPname
        IMPLICIT NONE
        INP%FNAME=INPname(Iwaves)
# endif
      END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
#if !defined SELFE
# if defined MPI_PARALL_GRID
      SUBROUTINE SIMPLE_PRE_READ
      USE DATAPOOL
      USE ELFE_GLBL, only : msc2, mdc2, ics
      IMPLICIT NONE
      CHARACTER(LEN=20) :: BEGTC, UNITC, ENDTC
      REAL(rkind) DELTC
         NAMELIST /PROC/ PROCNAME, DIMMODE, LSTEA, LQSTEA, LSPHE,       &
     &      LNAUTIN, LNAUTOUT, LMONO_OUT, LMONO_IN,                     &
     &      BEGTC, DELTC, UNITC, ENDTC, DMIN 
         NAMELIST /GRID/ LCIRD, LSTAG, MINDIR, MAXDIR, MDC, FRLOW,      &
     &      FRHIGH, MSC, FILEGRID, IGRIDTYPE, LSLOP, SLMAX, LVAR1D,     &
     &      LOPTSIG
      INTEGER FHNDL
      !
      FHNDL=12
      CALL TEST_FILE_EXIST_DIE("Missing input file : ", TRIM(INP%FNAME))
      OPEN(FHNDL, FILE = TRIM(INP%FNAME))
      READ(FHNDL, NML = PROC)
      IF (LSPHE) THEN
        ics=2
      ELSE
        ics=1
      ENDIF
      READ(FHNDL, NML = GRID)
      msc2=MSC
      mdc2=MDC
      CLOSE(FHNDL)
      END SUBROUTINE
# endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
# ifdef PGMCL_COUPLING
      SUBROUTINE WWMIII_MPI(MyCOMM)
# else
      PROGRAM WWMIII_MPI
# endif
# ifdef PGMCL_COUPLING
      USE mod_coupler, only : WAV_COMM_WORLD
# endif
      USE DATAPOOL, only: MAIN, SEBO,                                  &
     &      NDT_BND_FILE, IWBNDLC, AC2, WBAC, STAT, RTIME,             &
     &      bnd_time_all_files, LSPHE, WLDEP, DEP, SMALL, KKK,         &
     &      WATLEV, LBCSE, LBCWA, LBCSP, IWBMNP, IWBNDLC, AC2, WBAC,   &
     &      WBACOLD, WBACNEW, DSPEC, LBINTER, LFIRSTSTEP, LQSTEA,      &
     &      LINHOM, IBOUNDFORMAT, LCALC, DAY2SEC, SEC2DAY,             &
     &      NUM_NETCDF_FILES_BND, LSECU
# ifdef WWM_MPI
      use elfe_glbl
      use elfe_msgp
# else
      use datapool, only : MyRank, nproc, comm
# endif
      implicit none
# if defined MPI_PARALL_GRID || defined PGMCL_COUPLING
      include 'mpif.h'
# endif
# ifdef PGMCL_COUPLING
      integer, intent(in) :: MyCOMM
# endif
# ifdef MPI_PARALL_GRID
      integer i, j
# endif
      integer :: k
      character(len=15) CALLFROM
# if !defined PGMCL_COUPLING && defined WWM_MPI
      call mpi_init(ierr)
      if(ierr/=MPI_SUCCESS) call wwm_abort('Error at mpi_init')
# endif
      CALL SET_WWMINPULNML
# ifdef PGMCL_COUPLING
      comm=MyCOMM
      WAV_COMM_WORLD=MyCOMM
# else
#  ifdef WWM_MPI
      call mpi_comm_dup(MPI_COMM_WORLD,comm,ierr)
      if(ierr/=MPI_SUCCESS) call wwm_abort('Error at mpi_comm_dup')
#  endif
# endif
# ifdef MPI_PARALL_GRID
      call mpi_comm_size(comm,nproc,ierr)
      if(ierr/=MPI_SUCCESS) call wwm_abort('Error at mpi_comm_size')
      call mpi_comm_rank(comm,myrank,ierr)
      if(ierr/=MPI_SUCCESS) call wwm_abort('Error at mpi_comm_rank')
      CALL SIMPLE_PRE_READ
      do i=1,3
        do j=1,2
          nx(i,j)=i+j
          if(nx(i,j)>3) nx(i,j)=nx(i,j)-3
          if(nx(i,j)<1.or.nx(i,j)>3) then
            write(errmsg,*)'MAIN: nx wrong',i,j,nx(i,j)
            call wwm_abort(errmsg)
          endif
        enddo
      enddo
      call partition_hgrid
      call aquire_hgrid(.true.)

      call msgp_tables
      call msgp_init
      call parallel_barrier
      CALL INITIALIZE_WWM
      CALLFROM='WWM_MPI'

      DO K = 1, MAIN%ISTP
        IF (LQSTEA) THEN
          CALL QUASI_STEADY(K)
        ELSE
          CALL UN_STEADY(K,CALLFROM)
        END IF
      END DO
# else
      CALL INITIALIZE_WWM
      CALLFROM='WWM'
      DO K = 1, MAIN%ISTP
        IF (LQSTEA) THEN
          CALL QUASI_STEADY(K)
        ELSE
          CALL UN_STEADY(K,CALLFROM)
        END IF
      END DO
# endif
# if defined MPI_PARALL_GRID && !defined PGMCL_COUPLING
      call parallel_finalize
# endif
# ifdef PGMCL_COUPLING
      END SUBROUTINE
# else
      END PROGRAM
# endif
#endif
