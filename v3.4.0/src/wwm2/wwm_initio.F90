#include "wwm_functions.h"
#undef DEBUG
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE INIT_ARRAYS

         USE DATAPOOL
         IMPLICIT NONE

#ifdef MPI_PARALL_GRID
         INTEGER :: IE, IP
#endif
         integer istat

         ALLOCATE( DX1(0:MNP+1), DX2(0:MNP+1), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 1')
         DX1 = zero
         DX2 = zero

         ALLOCATE( XP(MNP), YP(MNP), DEP(MNP), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 2')
         DEP = zero

         ALLOCATE( INVSPHTRANS(MNP,2), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 3')
         INVSPHTRANS = zero

         ALLOCATE( INE(3,MNE), IEN(6,MNE), TRIA(MNE), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 4')
         INE = 0
         IEN = zero
         TRIA = zero
#ifdef MPI_PARALL_GRID
         DO IE = 1, MNE
           INE(:,IE) = INETMP(:,IE)   !todo: this used to be an index swap in wwm, but now they are the same. Still copy or use selfe as-is?
         END DO
#else
         XP  = zero
         YP  = zero
#endif
!
! spectral grid - shared
!
         ALLOCATE(MSC_HF(MNP), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 5')
         MSC_HF = MSC
!
! action densities and source terms - shared
!
         ALLOCATE (AC2(MNP,MSC,MDC), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 8')
         AC2 = zero

         ALLOCATE (AC1(MNP,MSC,MDC), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 9')
         AC1 = zero

         IF (ICOMP .GE. 2) THEN
           ALLOCATE (IMATRAA(MNP,MSC,MDC), IMATDAA(MNP,MSC,MDC), stat=istat)
           IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 10')
           IMATRAA = zero
           IMATDAA = zero
         END IF

         IF (LITERSPLIT) THEN
!           ALLOCATE (AC1(MNP,MSC,MDC)); AC1 = zero
           ALLOCATE (DAC_ADV(2,MNP,MSC,MDC), DAC_THE(2,MNP,MSC,MDC), DAC_SIG(2,MNP,MSC,MDC), DAC_SOU(2,MNP,MSC,MDC), stat=istat)
           IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 11')
           DAC_ADV = zero
           DAC_THE = zero
           DAC_SIG = zero
           DAC_SOU = zero
         END IF

#ifdef SHYFEM_COUPLING
         IF (LSHYFEM) THEN
           ALLOCATE(SHYFZETA(NLVT,MNP),NLEV(MNP), stat=istat)
           IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 12')
           SHYFZETA = ZERO
           NLEV = 0
           ALLOCATE(STOKES_X(NLVT,MNP), STOKES_Y(NLVT,MNP), JPRESS(MNP), stat=istat)
           IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 13')
           STOKES_X = ZERO; STOKES_Y = ZERO; JPRESS = ZERO
         END IF
#endif
!
! WAM Cycle 4.5 - shared
!
         IF (MESIN .EQ. 2) THEN
           ALLOCATE ( TAUHFT(0:IUSTAR,0:IALPHA,1), TAUT(0:ITAUMAX,0:JUMAX,1),stat=istat)
           IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 14')
           TAUHFT = zero; TAUT = zero
         ELSE IF (MESIN .EQ. 6) THEN
           ALLOCATE(TAUHFT(0:IUSTAR,0:IALPHA,MSC), stat=istat)
           IF (istat/=0) CALL WWM_ABORT('wwm_ecmwf, allocate error 14a')
           ALLOCATE(TAUT(0:ITAUMAX,0:JUMAX,JPLEVT), stat=istat)
           IF (istat/=0) CALL WWM_ABORT('wwm_ecmwf, allocate error 14b')
           TAUHFT = zero; TAUT = zero
           INQUIRE(FILE='fort.5010',EXIST=LPRECOMP_EXIST)
           INQUIRE(FILE='fort.5011',EXIST=LPRECOMP_EXIST)
           IF (LPRECOMP_EXIST) THEN
             READ(5010) DELU, DELTAUW
             READ(5010) TAUT
             READ(5011) DELUST, DELALP 
             READ(5011) TAUHFT
           ENDIF
         ENDIF

! Boundary conditions - shared
!
         ALLOCATE( IOBDP(MNP), IOBPD(MDC,MNP), IOBP(MNP), IOBWB(MNP), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 15')
         IOBDP  = 1
         IOBPD  = 0
         IOBP   = 0
         IOBWB  = 1
!
! phase velocity, wave number, group velocity, dwdh, kh
!
         ALLOCATE( WK(MNP,MSC), CG(MNP,MSC), DCGDX(MNP,MSC), WC(MNP,MSC), DWKDX(MNP,MSC), DCGDY(MNP,MSC), DWKDY(MNP,MSC), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 16')
         WK = ZERO 
         CG = ZERO 
         DWKDX = ZERO
         DCGDX = ZERO
#ifdef SELFE
!         ALLOCATE( CGX(MNP,MSC,MDC) ); CGX = zero
!         ALLOCATE( CGY(MNP,MSC,MDC) ); CGY = zero
#endif
!
         ALLOCATE( TABK (0:IDISPTAB), TABCG(0:IDISPTAB), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 17')
         TABK  = zero
         TABCG = zero
!
! diffraction parameter - shared
!
         IF (LDIFR) THEN
           ALLOCATE ( DIFRM(MNP), DIFRX(MNP), DIFRY(MNP), stat=istat)
           IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 18')
           DIFRM = zero
           DIFRX = zero
           DIFRY = zero
         END IF
!
! water level, currents and depths ...
!
         ALLOCATE( WINDXY(MNP,2), PRESSURE(MNP), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 19')
         WINDXY = zero
         PRESSURE = zero
         ALLOCATE( DVWIND(MNP,2), CURTXY(MNP,2), DVCURT(MNP,2), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 20')
         DVWIND = zero
         CURTXY = zero
         DVCURT = zero
         ALLOCATE( DDEP(MNP,2), DCUX(MNP,2), DCUY(MNP,2), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 21')
         DDEP = zero
         DCUX = zero
         DCUY = zero
         ALLOCATE( WATLEV(MNP), WATLEVOLD(MNP), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 22')
         WATLEV = zero
         WATLEVOLD = zero
         ALLOCATE( DVWALV(MNP), WLDEP(MNP), DEPDT(MNP), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 23')
         DVWALV = zero
         WLDEP = zero
         DEPDT = zero
!
!  convergence analysis - shared
!
         IF (LCONV .OR. (LQSTEA .AND. LCHKCONV)) THEN
           ALLOCATE ( SUMACOLD(MNP), HSOLD(MNP), KHSOLD(MNP), TM02OLD(MNP), IP_IS_STEADY(MNP), IE_IS_STEADY(MNE), STAT2D(MSC,MDC), stat=istat)
           IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 24')
           IP_IS_STEADY = 0
           IE_IS_STEADY = 0
           SUMACOLD = zero
           HSOLD  = zero
           TM02OLD = zero
           KHSOLD  = zero
         END IF
!
!  output - shared
!
         ALLOCATE( QBLOCAL(MNP), DISSIPATION(MNP), AIRMOMENTUM(MNP), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 25')
         QBLOCAL = zero
         DISSIPATION = zero
         AIRMOMENTUM = zero

         ALLOCATE( UFRIC(MNP), ALPHA_CH(MNP), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 26')
         UFRIC = zero
         ALPHA_CH = zero

         ALLOCATE( TAUW(MNP), TAUTOT(MNP), TAUWX(MNP), TAUWY(MNP), TAUHF(MNP), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 27')
         TAUW = zero
         TAUTOT = zero
         TAUWX = zero
         TAUWY = zero
         TAUHF = zero

         ALLOCATE( Z0(MNP), CD(MNP), USTDIR(MNP), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 28')
         Z0 = zero
         CD = zero
         USTDIR = zero

         ALLOCATE( RSXX(MNP), RSXY(MNP), RSYY(MNP), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 29')
         RSXX = zero
         RSXY = zero
         RSYY = zero

#ifdef SELFE
         ALLOCATE( SXX3D(NVRT,MNP), SXY3D(NVRT,MNP), SYY3D(NVRT,MNP), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 30')
         SXX3D = zero
         SXY3D = zero
         SYY3D = zero
#endif
#ifndef SELFE
         IF (LCPL) THEN
           IF (LTIMOR.or.LSHYFEM) THEN
             ALLOCATE( SXX3D(NLVT,MNP), SXY3D(NLVT,MNP), SYY3D(NLVT,MNP), stat=istat)
             IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 31')
             SXX3D = zero
             SXY3D = zero
             SYY3D = zero
           END IF
         END IF
#endif
         ALLOCATE(HMAX(MNP), ISHALLOW(MNP), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 32')
         HMAX = zero
         ISHALLOW = zero
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE DEALLOC_ARRAYS

         USE DATAPOOL
         IMPLICIT NONE

#ifdef MPI_PARALL_GRID
         INTEGER :: IE, IP
#endif

         DEALLOCATE( DX1, DX2, XP, YP, INVSPHTRANS, DEP, INE, IEN, TRIA)
!
! spectral grid - shared
!
         DEALLOCATE( SPSIG, SPDIR, FR, COSTH, SINTH, COS2TH, SIN2TH)
         DEALLOCATE( SINCOSTH, SIGPOW, DS_BAND, DS_INCR)
!
! action densities and source terms - shared
!
         DEALLOCATE (AC2)

         DEALLOCATE (AC1)
         IF (ICOMP .GE. 2) THEN
           DEALLOCATE (IMATRAA, IMATDAA)
         END IF

         IF (LITERSPLIT) THEN
           DEALLOCATE (DAC_ADV, DAC_THE, DAC_SIG, DAC_SOU)
         END IF

#ifdef SHYFEM_COUPLING
         IF (LSHYFEM) THEN
           DEALLOCATE(SHYFZETA, NLEV)
         END IF
#endif
!
! WAM Cycle 4.5 - shared
!
         IF (MESIN .EQ. 2) THEN
           DEALLOCATE ( TAUHFT, TAUT)
         ENDIF
!
! Boundary conditions - shared
!
         DEALLOCATE( IOBPD, IOBP, IOBWB)
!
! phase velocity, wave number, group velocity, dwdh, kh
!
         DEALLOCATE( WK, CG, TABK, TABCG)
!
! diffraction parameter - shared
!
         IF (LDIFR) THEN
           DEALLOCATE ( DIFRM, DIFRX, DIFRY )
         END IF
!
! water level, currents and depths ...
!
         DEALLOCATE( WINDXY, PRESSURE, DVWIND, CURTXY, DVCURT)
         DEALLOCATE( DDEP, DCUX, DCUY, WATLEV, WATLEVOLD)
         DEALLOCATE( DVWALV, WLDEP, DEPDT)
!
!  convergence analysis - shared
!
         IF (LCONV .OR. (LQSTEA .AND. LCHKCONV)) THEN
           DEALLOCATE ( SUMACOLD, HSOLD, KHSOLD, TM02OLD)
         END IF
!
!  output - shared
!
         DEALLOCATE( QBLOCAL, DISSIPATION, AIRMOMENTUM, UFRIC, ALPHA_CH)
         DEALLOCATE( TAUW, TAUTOT, TAUWX, TAUWY, TAUHF)
         DEALLOCATE( Z0, CD, USTDIR)
         DEALLOCATE( RSXX, RSXY, RSYY)
#ifdef SELFE
         DEALLOCATE( SXX3D, SXY3D, SYY3D)
#endif
#ifndef SELFE
         IF (LCPL) THEN
           IF (LTIMOR.or.LSHYFEM) THEN
             DEALLOCATE( SXX3D, SXY3D, SYY3D)
           END IF
         END IF
#endif
         DEALLOCATE(HMAX, ISHALLOW)
         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INITIALIZE_WWM
         USE DATAPOOL
#ifdef WWM_SOLVER
# ifdef MPI_PARALL_GRID
         USE WWM_PARALL_SOLVER, only : WWM_SOLVER_INIT
# endif
#endif
#ifdef PETSC
         USE PETSC_CONTROLLER, ONLY : PETSC_INIT
#endif
#ifdef ST41
         USE W3SRC4MD_OLD
#endif
#ifdef ST42
         USE W3SRC4MD_NEW
#endif
#ifdef MPI_PARALL_GRID
         use elfe_glbl, only : iplg, np_global, ne_global
         use elfe_msgp, only : myrank, nproc, parallel_abort, comm, ierr
         use elfe_msgp, only : rtype, itype
#endif
         IMPLICIT NONE
         integer istat
#ifdef MPI_PARALL_GRID
         integer FHNDLspec
         CHARACTER(LEN=40) :: FILEspec
#endif
         REAL           :: TIME1, TIME2
         INTEGER        :: IT, IFILE
         CALL CPU_TIME(TIME1)

         CALL INIT_FILE_HANDLES
         CALL READ_WWMINPUT
         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'DONE READING NAMELIST'
         CALL FLUSH(STAT%FHNDL)

#ifdef VDISLIN
         CALL INIT_DISLIN()
         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'INIT DISLIN                '
         CALL FLUSH(STAT%FHNDL)
#endif

         CALL INIT_ARRAYS
         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ARRAY INITIALIZATION'
         CALL FLUSH(STAT%FHNDL)

#ifdef MPI_PARALL_GRID 
         DEP  = DEP8
         WLDEP  = DEP
         IF (LSPHE) THEN
           XP = XLON*RADDEG
           YP = YLAT*RADDEG
         ELSE
           XP = XPTMP
           YP = YPTMP
         END IF
#endif
         CALL CHECK_LOGICS
         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'CHECK LOGICS                '
         CALL FLUSH(STAT%FHNDL)

#ifndef MPI_PARALL_GRID
         CALL READ_SPATIAL_GRID
         WLDEP = DEP
         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'READ SPATIAL GRID'
         CALL FLUSH(STAT%FHNDL)
#endif
         CALL SPATIAL_GRID
         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'INIT SPATIAL GRID'
         CALL FLUSH(STAT%FHNDL)

         IF (LADVTEST) THEN
           ALLOCATE(UTEST(MNP), stat=istat)
           IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 33')
           UTEST = 0.
           CALL ADVTEST(UTEST)
           AC2(:,1,1) = UTEST
           CALL CHECKCONS(UTEST,SUMACt0)
!AR: check the advections test if it still works ...
           DEALLOCATE(UTEST)
         END IF
#ifdef MPI_PARALL_GRID
         CALL BUILD_WILD_ARRAY
#endif
#ifdef MPI_PARALL_GRID
         NP_TOTAL=np_global
         NE_TOTAL=ne_global
#else
         NP_TOTAL=MNP
         NE_TOTAL=MNE
#endif
         CALL SET_IOBPD_BY_DEP
         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'SET DEPTH POINTER'
         CALL FLUSH(STAT%FHNDL)

         IF (DIMMODE .EQ. 2) THEN
           WRITE(STAT%FHNDL,'("+TRACE...",A)') 'THE FLUCTUATION SPLITTING PREPROCESSOR HAS STARTED'
           CALL FLUSH(STAT%FHNDL)
           CALL INIT_FLUCT_ARRAYS
           CALL INIT_FLUCT

           IF (AMETHOD .EQ. 4 .OR. AMETHOD .EQ. 5) THEN
#ifdef PETSC
             CALL PETSC_INIT
#endif
           END IF
           IF (AMETHOD .EQ. 6) THEN
#if defined WWM_SOLVER && defined MPI_PARALL_GRID
             CALL WWM_SOLVER_INIT
#endif
           END IF
           WRITE(STAT%FHNDL,'("+TRACE...",A)') 'THE FLUCTUATION SPLITTING PREPROCESSOR HAS ENDED'
           CALL FLUSH(STAT%FHNDL)
         END IF

         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'INITIALIZE SPECTRAL GRID'
         CALL FLUSH(STAT%FHNDL)
         CALL INIT_SPECTRAL_GRID
 
         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'INITIALIZE BOUNDARY POINTER 1/2'
         CALL FLUSH(STAT%FHNDL)
#if defined SELFE 
         DMIN = DMIN_SELFE
         CALL SET_IOBP_SELFE
#else
         CALL SET_IOBP_NEXTGENERATION
#endif
         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'INITIALIZE BOUNDARY POINTER 2/2'
         CALL FLUSH(STAT%FHNDL)
         CALL SET_IOBPD

#ifndef PGMCL_COUPLING
         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'INITIALIZE WIND CURRENT WATERLEVEL'
         CALL FLUSH(STAT%FHNDL)
         IF (LWINDFROMWWM) THEN
           CALL INIT_WIND_INPUT
         ENd IF
         IF (.NOT. LCPL) THEN
           CALL INIT_CURRENT_INPUT
           CALL INIT_WATLEV_INPUT
         END IF
#endif
         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'COMPUTE THE WAVE PARAMETER'
         CALL FLUSH(STAT%FHNDL)
         CALL INITIATE_WAVE_PARAMETER
         CALL SETSHALLOW
         CALL SET_HMAX

         !CALL NLWEIGT(MSC,MDC)

         IF ( (MESIN .EQ. 1 .OR. MESDS .EQ. 1) .AND. SMETHOD .GT. 0) THEN
           WRITE(STAT%FHNDL,'("+TRACE...",A)') 'INIT ARDHUIN et al.'
           CALL FLUSH(STAT%FHNDL)
#ifdef ST41
           CALL PREPARE_ARDHUIN_OLD
#elif ST42
           CALL PREPARE_ARDHUIN_NEW
#endif
         ELSE IF (MESIN .EQ. 2 .AND. .NOT. LPRECOMP_EXIST) THEN
           WRITE(STAT%FHNDL,'("+TRACE...",A)') 'INIT CYCLE4 WIND INPUT'
           CALL FLUSH(STAT%FHNDL)
           CALL PREPARE_SOURCE
           CALL FLUSH(STAT%FHNDL)
         ELSE IF (MESIN .EQ. 6 .AND. .NOT. LPRECOMP_EXIST) THEN
           WRITE(STAT%FHNDL,'("+TRACE...",A)') 'INIT CYCLE4 WIND INPUT'
           CALL PREPARE_SOURCE_ECMWF 
           CALL FLUSH(STAT%FHNDL)
         END IF

         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'SET THE INITIAL WAVE BOUNDARY CONDITION'
         CALL FLUSH(STAT%FHNDL)
         CALL INIT_WAVE_BOUNDARY_CONDITION(IFILE,IT)
         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'SET THE INITIAL CONDITION'
         CALL FLUSH(STAT%FHNDL)
         CALL INITIAL_CONDITION(IFILE,IT)
         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'SET BOUNDARY CONDITIONS'
         CALL FLUSH(STAT%FHNDL)
         CALL SET_WAVE_BOUNDARY

         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'INIT STATION OUTPUT'
         CALL FLUSH(STAT%FHNDL)
         CALL INIT_STATION_OUTPUT

         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'WRITING INITIAL TIME STEP'
         CALL FLUSH(STAT%FHNDL)
         CALL WWM_OUTPUT(ZERO,.TRUE.)

         IF (LWXFN) THEN
           CALL WRINPGRD_XFN
         ELSE IF(LWSHP) THEN
           CALL WRINPGRD_SHP
         END IF

#if !defined SELFE && !defined PGMCL_COUPLING
         IF (LCPL) THEN
           WRITE(STAT%FHNDL,'("+TRACE...",A)') 'OPEN PIPES FOR COUPLING'
           CALL FLUSH(STAT%FHNDL)
           IF (LTIMOR) THEN
             CALL INIT_PIPES_TIMOR()
#ifdef SHYFEM_COUPLING
           ELSE IF (LSHYFEM) THEN
             CALL INIT_PIPES_SHYFEM()
#endif
           ELSE IF (LROMS) THEN
             CALL INIT_PIPES_ROMS()
           END IF
         END IF
#endif
#ifdef PGMCL_COUPLING
         CALL ROMS_COUPL_INITIALIZE
#endif
!#ifdef MPI_PARALL_GRID 
!         IF (size(OUTT_INTPAR(1,:)) .NE. OUTVARS) THEN
!           call parallel_abort('OUTPUTSIZE IS NOT RIGHT CHECK OUTVARS')
!         END IF
!         IF (size(WIND_INTPAR(1,:)) .NE. WINDVARS) THEN
!           call parallel_abort('OUTPUTSIZE IS NOT RIGHT CHECK VARS')
!         END IF
!#endif
         CALL CPU_TIME(TIME2)

#if defined SELFE
         IF (MSC_SELFE .NE. MSC .OR. MDC_SELFE .NE. MDC) THEN
           WRITE(DBG%FHNDL,*) 'MSC_SELFE', MSC_SELFE
           WRITE(DBG%FHNDL,*) 'MSC', MSC
           WRITE(DBG%FHNDL,*) 'MDC_SELFE', MDC_SELFE
           WRITE(DBG%FHNDL,*) 'MDC', MDC
           CALL FLUSH(DBG%FHNDL)
           CALL PARALLEL_ABORT('THERE IS AND ERROR IN MSC2 OR MDC2 IN PARAM.IN')
         END IF
#endif

         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'INITIALIZE_WWM'
         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.4)') 'CPU Time for the preprocessing', TIME2-TIME1
         CALL FLUSH(STAT%FHNDL)

         AC1 = AC2
       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE TERMINATE_WWM
         USE DATAPOOL
#ifdef ST41
         USE W3SRC4MD_OLD
#endif
#ifdef ST42
         USE W3SRC4MD_NEW
#endif

#ifdef PETSC
         USE PETSC_CONTROLLER, ONLY : PETSC_FINALIZE
#endif

         CALL CLOSE_FILE_HANDLES
         CALL DEALLOC_ARRAYS
#ifdef MPI_PARALL_GRID
         CALL DEALLOC_WILD_ARRAY
#endif
         IF (DIMMODE .EQ. 2) THEN
           CALL DEALLOC_FLUCT_ARRAYS
           CALL DEALLOC_FLUCT

           IF (AMETHOD .EQ. 4 .OR. AMETHOD .EQ. 5) THEN
#ifdef PETSC
             CALL PETSC_FINALIZE
#endif
           END IF
         END IF
         CALL DEALLOC_SPECTRAL_GRID
         CALL CLOSE_IOBP

         IF ( (MESIN .EQ. 1 .OR. MESDS .EQ. 1) .AND. SMETHOD .GT. 0) THEN
#ifdef ST41
!           CALL PREPARE_ARDHUIN_OLD
#elif ST42
!           CALL PREPARE_ARDHUIN_NEW
#endif
         ELSE IF (MESIN .EQ. 2) THEN
!           CALL PREPARE_SOURCE
         ELSE IF (MESIN .EQ. 6) THEN
!           CALL PREPARE_SOURCE_ECMWF 
         END IF
!         CALL SET_WAVE_BOUNDARY

         CALL TERMINATE_STATION_OUTPUT
#if !defined SELFE && !defined PGMCL_COUPLING
         IF (LCPL) THEN
           IF (LTIMOR) THEN
             CALL TERMINATE_PIPES_TIMOR()
# ifdef SHYFEM_COUPLING
           ELSE IF (LSHYFEM) THEN
             CALL TERMINATE_PIPES_SHYFEM()
# endif
           ELSE IF (LROMS) THEN
             CALL TERMINATE_PIPES_ROMS()
           END IF
         END IF
#endif
#ifdef PGMCL_COUPLING
         CALL ROMS_COUPL_DEALLOCATE
#endif

       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef MPI_PARALL_GRID 
       SUBROUTINE BUILD_WILD_ARRAY
       USE DATAPOOL, only : MNP, NP_RES, ONE, rkind
       USE DATAPOOL, only : nwild_gb, nwild_loc, nwild_loc_res
       USE elfe_glbl, only : iplg, np_global
       USE elfe_msgp, only : istatus, itype, comm, ierr, myrank, rtype, nproc
       IMPLICIT NONE
       include 'mpif.h'
       integer, allocatable :: nwild_i(:), nwild_gbi(:)
       integer :: iProc, istat
       INTEGER :: IP
       real(rkind), allocatable :: nwild_gb_res(:)
       ALLOCATE(nwild_i(np_global), nwild_gbi(np_global), stat=istat)
       IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 33')
       !
       nwild_i=0
       nwild_gbi=0
       DO IP=1,MNP
         nwild_i(IPLG(IP))=1
       END DO
       call mpi_reduce(nwild_i,nwild_gbi,NP_GLOBAL,itype,MPI_SUM,0,comm,ierr)
       ALLOCATE(nwild_gb(NP_GLOBAL), stat=istat)
       IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 34')
       IF (myrank.eq.0) THEN
         DO IP=1,NP_GLOBAL
           nwild_gb(IP)=ONE/MyREAL(nwild_gbi(IP))
         END DO
         DO iProc=2,nproc
           CALL MPI_SEND(nwild_gb,np_global,rtype, iProc-1, 2037, comm, ierr)
         END DO
       ELSE
         CALL MPI_RECV(nwild_gb,np_global,rtype, 0, 2037, comm, istatus, ierr)
       END IF
       ALLOCATE(nwild_loc(MNP), stat=istat)
       IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 35')
       DO IP=1,MNP
         nwild_loc(IP)=nwild_gb(iplg(IP))
       END DO
       IF (myrank.ne.0) THEN
         deallocate(nwild_gb)
       END IF
       !
       nwild_i=0
       nwild_gbi=0
       DO IP=1,NP_RES
         nwild_i(IPLG(IP))=1
       END DO
       call mpi_reduce(nwild_i,nwild_gbi,NP_GLOBAL,itype,MPI_SUM,0,comm,ierr)
       ALLOCATE(nwild_gb_res(NP_GLOBAL), stat=istat)
       IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 36')
       IF (myrank.eq.0) THEN
         DO IP=1,NP_GLOBAL
           nwild_gb_res(IP)=ONE/MyREAL(nwild_gbi(IP))
         END DO
         DO iProc=2,nproc
           CALL MPI_SEND(nwild_gb_res,np_global,rtype, iProc-1, 277, comm, ierr)
         END DO
       ELSE
         CALL MPI_RECV(nwild_gb_res,np_global,rtype, 0, 277, comm, istatus, ierr)
       END IF
       ALLOCATE(nwild_loc_res(NP_RES), stat=istat)
       IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 37')
       DO IP=1,NP_RES
         nwild_loc_res(IP)=nwild_gb_res(iplg(IP))
       END DO
       deallocate(nwild_gb_res)
       !
       DEALLOCATE(nwild_i)
       DEALLOCATE(nwild_gbi)
       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE DEALLOC_WILD_ARRAY
!AR: what is this kind of shit ?
       USE DATAPOOL
       implicit none
       DEALLOCATE(nwild_loc)
       DEALLOCATE(nwild_loc_res)
       END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE BASIC_PARAMETER
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER :: IP
         REAL(rkind)    :: PPTAIL
!
!2do: Check with WAM parameterization and WW3
!
         PQUAD(1) = 0.25
         IF (MESIN .EQ. 1) THEN
           PQUAD(2) = 2.5E7
         ELSE
           PQUAD(2) = 2.78E7
         END IF
         PQUAD(3) = 5.5_rkind
         PQUAD(4) = 6._rkind/7._rkind
         PQUAD(5) = -1.25_rkind
         PQUAD(6) = zero
!
! High frequency tail integration factors as defined in the SWAN model.
!
         PTAIL(1) = 4._rkind
         PTAIL(2) = 2.5_rkind
         PTAIL(3) = PTAIL(1) + 1._rkind
         IF (MESIN .EQ. 2) THEN
           PTAIL(1) = 5._rkind
           PTAIL(3) = PTAIL(1) + 1._rkind
         ENDIF
         PTAIL(4) = 3._rkind
         DO IP = 0, 3
           PPTAIL = PTAIL(1) - MyREAL(IP)
           PTAIL(5+IP) = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
         ENDDO
!
! Output Variables ... LVARS, NVARS
!
         OUTT_VARNAMES(1)  = 'HS'
         OUTT_VARNAMES(2)  = 'TM01'
         OUTT_VARNAMES(3)  = 'TM02'
         OUTT_VARNAMES(4)  = 'TM10'
         OUTT_VARNAMES(5)  = 'KLM'
         OUTT_VARNAMES(6)  = 'WLM'
         OUTT_VARNAMES(7)  = 'ETOTS'
         OUTT_VARNAMES(8)  = 'ETOTC'
         OUTT_VARNAMES(9)  = 'DM'
         OUTT_VARNAMES(10)  = 'DSPR'
         OUTT_VARNAMES(11) = 'TPPD'
         OUTT_VARNAMES(12) = 'TPP'
         OUTT_VARNAMES(13) = 'CPP'
         OUTT_VARNAMES(14) = 'WNPP'
         OUTT_VARNAMES(15) = 'CGPP'
         OUTT_VARNAMES(16) = 'KPP'
         OUTT_VARNAMES(17) = 'LPP'
         OUTT_VARNAMES(18) = 'PEAKD'
         OUTT_VARNAMES(19) = 'PEAKDSPR'
         OUTT_VARNAMES(20) = 'DPEAK'
         OUTT_VARNAMES(21) = 'UBOT'
         OUTT_VARNAMES(22) = 'ORBITAL'
         OUTT_VARNAMES(23) = 'BOTEXPER'
         OUTT_VARNAMES(24) = 'TMBOT'
         OUTT_VARNAMES(25) = 'URSELL'
         OUTT_VARNAMES(26) = 'USTAR'
         OUTT_VARNAMES(27) = 'ALPHA'
         OUTT_VARNAMES(28) = 'Z0'
         OUTT_VARNAMES(29) = 'WIND-X'
         OUTT_VARNAMES(30) = 'WIND-Y'
         OUTT_VARNAMES(31) = 'CD'
         OUTT_VARNAMES(32) = 'CURR-X'
         OUTT_VARNAMES(33) = 'CURR-Y'
         OUTT_VARNAMES(34) = 'DEPTH'
         OUTT_VARNAMES(35) = 'ELEVATION'

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INITIATE_WAVE_PARAMETER()
         USE DATAPOOL
         USE M_CONSTANTS
         USE M_XNLDATA
         USE M_FILEIO

         IMPLICIT NONE
         INTEGER IQGRID, INODE, IERR

         WRITE(STAT%FHNDL,*) 'START WAVE PARAMETER'
         CALL FLUSH(STAT%FHNDL)
         CALL GRADDEP
         WRITE(STAT%FHNDL,*) 'GRADDEP'
         CALL FLUSH(STAT%FHNDL)
         IF (LSTCU .OR. LSECU) CALL GRADCURT
         WRITE(STAT%FHNDL,*) 'GRADCURT'
         CALL FLUSH(STAT%FHNDL)
         CALL BASIC_PARAMETER
         WRITE(STAT%FHNDL,*) 'BASIC'
         CALL FLUSH(STAT%FHNDL)
         CALL MAKE_WAVE_TABLE
         CALL WAVE_K_C_CG
         WRITE(STAT%FHNDL,*) 'WAVEKCG'
         CALL FLUSH(STAT%FHNDL)


         IF (MESNL .LT. 4) THEN
           CALL PARAMETER4SNL
         ELSE IF (MESNL .EQ. 5) THEN
           CALL XNL_INIT(REAL(SPSIG),REAL(SPDIR),MSC,MDC,-4.0,REAL(G9),REAL(DEP),MNP,1,IQGRID,INODE,IERR)
           CALL INIT_CONSTANTS()
           IQGRID = 3
           INODE  = 1
           CALL XNL_INIT(REAL(SPSIG),REAL(SPDIR),MSC,MDC,-4.0,REAL(G9),REAL(DEP),MNP,1,IQGRID,INODE,IERR)
           IF (IERR .GT. 0) CALL WWM_ABORT('IERR XNL_INIT')
         ENDIF

         IF (MESTR .GT. 5) CALL GRAD_CG_K 

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INITIAL_CONDITION(IFILE,IT)
         USE WWM_HOTFILE_MOD
         USE DATAPOOL
#ifdef NCDF
         USE NETCDF
#endif

#ifdef MPI_PARALL_GRID
         use elfe_glbl, only : iplg, np_global
         use elfe_msgp, only : myrank
#endif
         IMPLICIT NONE

         INTEGER         :: IP, IFILE, IT
         REAL(rkind)     :: SPPAR(8)
         REAL(rkind)     :: MS
         REAL(rkind)     :: HS, TP, HSLESS, TPLESS, FDLESS
         REAL(rkind)     :: WIND10, WINDTH, VEC2DEG
         REAL(rkind)     :: WINDX, WINDY
         REAL(rkind)     :: ACLOC(MSC,MDC)
         REAL(rkind)     :: DEG
         REAL(rkind)     :: TMPPAR(8,MNP), SSBRL(MSC,MDC)

         CHARACTER(len=25) :: CALLFROM

         TMPPAR = 0.

         IF (.NOT. LHOTR .AND. LINID) THEN
            IF (INITSTYLE == 2 .AND. IBOUNDFORMAT == 3) THEN
#ifdef NCDF
              CALLFROM = 'CALL FROM INIT. CONDITION'
              CALL READ_NETCDF_WW3(IFILE,IT,CALLFROM)
#else
              CALL WWM_ABORT('compile with DNCDF PPFLAG')
#endif
              CALL INTER_STRUCT_DOMAIN(NDX_BND,NDY_BND,DX_BND,DY_BND,OFFSET_X_BND,OFFSET_Y_BND,TMPPAR)
!test the initial conditions ...
              WRITE(2001) 1.
              WRITE(2001) (TMPPAR(3,IP), TMPPAR(2,IP), TMPPAR(1,IP), IP = 1, MNP)
            END IF
            DO IP = 1, MNP
               ACLOC = 0.
               WINDX  = WINDXY(IP,1)
               WINDY  = WINDXY(IP,2)
               WIND10 = SQRT(WINDX**2+WINDY**2)
!AR: SELFEWWM: WE HAVE TO CHECK WHY THE INITIAL WIND FIELD IS ZERO WHEN COUPLING TO SELFE
               IF(DEP(IP) .GT. DMIN) THEN
                 IF (DIMMODE .EQ. 1 .AND. IP .EQ. 1) CYCLE
                 IF (INITSTYLE == 1) THEN
                   IF (WIND10 .GT. 1.) THEN
                     WINDTH = VEC2DEG(WINDX,WINDY)
                     CALL DEG2NAUT(WINDTH, DEG, LNAUTIN)
                     FDLESS = G9*AVETL/WIND10**2
                     HSLESS = MIN(0.21_rkind, 0.00288_rkind*FDLESS**0.45_rkind)
                     TPLESS = MIN(ONE/0.13_rkind, 0.46_rkind*FDLESS**0.27_rkind)
                     HS = HSLESS*WIND10**2/G9
                     TP = TPLESS*WIND10/G9
                     MS = 2.0_rkind
                     SPPAR(1) = HS
                     SPPAR(2) = TP ! Check whether TP is inside of spectal representation !!!
                     SPPAR(3) = DEG
                     SPPAR(4) = MS
                     SPPAR(5) = 2.
                     SPPAR(6) = 2.
                     SPPAR(7) = 0.1
                     SPPAR(8) = 3.3
                     CALL SPECTRAL_SHAPE(SPPAR,ACLOC,.FALSE.,'INITIAL CONDITION PARA', .FALSE.)
                     AC2(IP,:,:) = ACLOC
                   ELSE
                     ACLOC = 1.E-8
                   END IF
                 ELSE IF (INITSTYLE == 2 .AND. IBOUNDFORMAT == 3) THEN
                   TMPPAR(5,IP) = 2.
                   TMPPAR(6,IP) = 1.
                   TMPPAR(7,IP) = 0.1
                   TMPPAR(8,IP) = 3.3
                   CALL SPECTRAL_SHAPE(TMPPAR(:,IP),ACLOC,.FALSE.,'INITIAL CONDITION WW3', .FALSE.)
                   AC2(IP,:,:) = ACLOC
                 END IF ! INITSTYLE
                 IF (LMAXETOT .AND. ISHALLOW(IP) .EQ. 1) CALL BREAK_LIMIT(IP,ACLOC,SSBRL) ! Miche for initial cond.
               ELSE
                 FDLESS      = 0.
                 HS          = 0.
                 TP          = 0.
                 WINDTH      = 0.
                 AC2(IP,:,:) = 0.
               END IF ! DEP(IP) .GT. DMIN .AND. WIND10 .GT. SMALL
            END DO ! IP
         END IF
         IF (LHOTR) THEN
           CALL INPUT_HOTFILE
         END IF
         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_FILE_HANDLES()
#ifdef MPI_PARALL_GRID
         use elfe_msgp
#endif
         USE DATAPOOL
         IMPLICIT NONE
         CHARACTER (LEN = 30) :: FDB
         INTEGER              :: LFDB

#ifdef SELFE
            INP%FNAME  = 'wwminput.nml'
#endif
            CHK%FNAME  = 'wwmcheck.nml'
           QSTEA%FNAME = 'qstea.out'
         WINDBG%FNAME  = 'winddbg.out'
         IOBPOUT%FNAME = 'iobp.out'
        IOBPDOUT%FNAME = 'iobpd.out'
!
!2do ... dinstinguish between binary and ascii stuff ...
!
             BND%FHNDL  = STARTHNDL + 1
             WIN%FHNDL  = STARTHNDL + 2
             CUR%FHNDL  = STARTHNDL + 3
             WAT%FHNDL  = STARTHNDL + 4
             WAV%FHNDL  = STARTHNDL + 5
             CHK%FHNDL  = STARTHNDL + 6
           HOTIN%FHNDL  = STARTHNDL + 7
          HOTOUT%FHNDL  = STARTHNDL + 8
             INP%FHNDL  = STARTHNDL + 9
             GRD%FHNDL  = STARTHNDL + 10
          GRDCOR%FHNDL  = STARTHNDL + 11

           QSTEA%FHNDL  = STARTHNDL + 12

         IOBPOUT%FHNDL  = STARTHNDL + 13
         IOBPDOUT%FHNDL = STARTHNDL + 14

         DBG%FHNDL      = STARTHNDL + 15 
         STAT%FHNDL     = STARTHNDL + 16 
         WINDBG%FHNDL   = STARTHNDL + 17 

#ifndef MPI_PARALL_GRID
         open(DBG%FHNDL,file='wwmdbg.out',status='unknown') !non-fatal errors
         open(STAT%FHNDL,file='wwmstat.out',status='unknown') !non-fatal errors
         open(WINDBG%FHNDL,file='windbg.out',status='unknown') !non-fatal errors
#else
# ifdef SELFE
         FDB  ='wwmdbg_0000'
         LFDB =len_trim(FDB)
         write(FDB(LFDB-3:LFDB),'(i4.4)') MYRANK
         open(DBG%FHNDL,file='outputs/'//fdb,status='replace') 
         FDB  ='wwmstat_0000'
         LFDB =len_trim(FDB)
         write(FDB(LFDB-3:LFDB),'(i4.4)') MYRANK
         open(STAT%FHNDL,file='outputs/'//fdb,status='replace') 
         FDB  ='windbg_0000'
         LFDB =len_trim(FDB)
         write(FDB(LFDB-3:LFDB),'(i4.4)') MYRANK
         open(WINDBG%FHNDL,file='outputs/'//fdb,status='replace') 
# else
         FDB  ='wwmdbg_0000'
         LFDB =len_trim(FDB)
         write(FDB(LFDB-3:LFDB),'(i4.4)') MYRANK
         open(DBG%FHNDL,file=fdb,status='replace') 
         FDB  ='wwmstat_0000'
         LFDB =len_trim(FDB)
         write(FDB(LFDB-3:LFDB),'(i4.4)') MYRANK
         open(STAT%FHNDL,file=fdb,status='replace') 
         FDB  ='windbg_0000'
         LFDB =len_trim(FDB)
         write(FDB(LFDB-3:LFDB),'(i4.4)') MYRANK
         open(WINDBG%FHNDL,file=fdb,status='replace')
# endif
#endif
         WRITE(DBG%FHNDL, *) 'THR=', THR
         WRITE(DBG%FHNDL, *) 'THR8=', THR8
         CALL TEST_FILE_EXIST_DIE("Missing input file : ", TRIM(INP%FNAME))
#ifdef MPIP_PARALL_GRID
         IF (myrank == 0) THEN
#endif
         WRITE(STAT%FHNDL,*) 'Input Filename   =', TRIM(INP%FNAME)
         WRITE(STAT%FHNDL,*) 'Check Filename   =', TRIM(CHK%FNAME)
         WRITE(STAT%FHNDL,*) 'Qstea Filename   =', TRIM(QSTEA%FNAME)
         WRITE(STAT%FHNDL,*) 'Iobp Filename    =', TRIM(IOBPOUT%FNAME)
         WRITE(STAT%FHNDL,*) 'Iobpd Filename   =', TRIM(IOBPDOUT%FNAME)
         WRITE(STAT%FHNDL,*) 'WindDbg Filename =', TRIM(WINDBG%FNAME)
#ifdef MPIP_PARALL_GRID
         ENDIF
#endif
         CALL TEST_FILE_EXIST_DIE("Missing input file : ", INP%FNAME)
         OPEN( INP%FHNDL,      FILE = TRIM(INP%FNAME))
         OPEN( CHK%FHNDL,      FILE = TRIM(CHK%FNAME))
         OPEN( QSTEA%FHNDL,    FILE = TRIM(QSTEA%FNAME))
         OPEN( IOBPOUT%FHNDL,  FILE = TRIM(IOBPOUT%FNAME))
         OPEN( IOBPDOUT%FHNDL, FILE = TRIM(IOBPDOUT%FNAME))

           OUT1D%FHNDL = STARTHNDL + 18 
            MISC%FHNDL = STARTHNDL + 19 
         OUTSP1D%FHNDL = STARTHNDL + 20 
         OUTPARM%FHNDL = STARTHNDL + 21 
         OUTSP2D%FHNDL = STARTHNDL + 22 

         OUT%FHNDL     = STARTHNDL + 23 

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CLOSE_FILE_HANDLES()
#ifdef MPI_PARALL_GRID
         use elfe_msgp
#endif
         USE DATAPOOL
         IMPLICIT NONE
         CHARACTER (LEN = 30) :: FDB
         INTEGER              :: LFDB
         close(DBG%FHNDL)
         close(STAT%FHNDL)
         close( QSTEA%FHNDL)
         close( IOBPOUT%FHNDL)
         close( IOBPDOUT%FHNDL)
         close( WINDBG%FHNDL)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_STATION_OUTPUT()
      USE DATAPOOL
#ifdef MPI_PARALL_GRID
      USE elfe_msgp
      USE elfe_glbl, only : iplg, ielg
#endif
      IMPLICIT NONE
#ifdef MPI_PARALL_GRID
      include 'mpif.h'
#endif
      INTEGER           :: I, NI(3), IP, IS
      integer istat
      REAL(rkind)              :: XYTMP(2,MNP)
#ifdef MPI_PARALL_GRID
      integer :: iProc
      integer, allocatable :: rbuf_int(:)
#endif
!
!    set the site output
!
      IF (LOUTS) THEN
        ALLOCATE (ACLOC_STATIONS(IOUTS,MSC,MDC), CDLOC_STATIONS(IOUTS), Z0LOC_STATIONS(IOUTS), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 38')

        ALLOCATE (ALPHALOC_STATIONS(IOUTS), WINDXLOC_STATIONS(IOUTS), WINDYLOC_STATIONS(IOUTS), USTARLOC_STATIONS(IOUTS), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 39')

        ALLOCATE (DEPLOC_STATIONS(IOUTS), WKLOC_STATIONS(IOUTS,MSC), CURTXYLOC_STATIONS(IOUTS,2), WATLEVLOC_STATIONS(IOUTS), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 40')
        ACLOC_STATIONS      = 0.
        DEPLOC_STATIONS     = 0.
        WKLOC_STATIONS      = 0.
        CURTXYLOC_STATIONS  = 0.
        USTARLOC_STATIONS   = 0.
        CDLOC_STATIONS      = 0.
        Z0LOC_STATIONS      = 0.
        WINDXLOC_STATIONS   = 0.
        WINDYLOC_STATIONS   = 0.
        WATLEVLOC_STATIONS  = 0.
      END IF
      IF (LOUTS .and. (DIMMODE .EQ. 2)) THEN
#ifdef MPI_PARALL_GRID
        XYTMP(1,:) = XP
        XYTMP(2,:) = YP
        WRITE(DBG%FHNDL,*) 'SEARCHING FOR STATION ACROSS RANKS', myrank
        DO I = 1, IOUTS
          CALL FIND_ELE ( MNE,MNP,INE,XYTMP,STATION(I)%XCOORD, STATION(I)%YCOORD,STATION(I)%ELEMENT )
          IF (STATION(I)%ELEMENT .GT. 0) THEN
            STATION(I)%IFOUND  = 1
            NI                 = INE(:,STATION(I)%ELEMENT)
            STATION(I)%XELE(:) = XP(NI)
            STATION(I)%YELE(:) = YP(NI)
            CALL INTELEMENT_COEF(XP(NI),YP(NI), STATION(I)%XCOORD,STATION(I)%YCOORD, STATION(I)%WI)
            WRITE(DBG%FHNDL,'(A10,I10,A20,I10,A15,2I10)') 'MYRANK', MYRANK, 'STATION =',I, 'IN ELEMENT =', IELG(STATION(I)%ELEMENT), STATION(I)%IFOUND
            CALL FLUSH(DBG%FHNDL)
          ELSE
            STATION(I)%IFOUND  = 0
            NI                 = 0
            STATION(I)%XELE(:) = 0.
            STATION(I)%YELE(:) = 0.
            WRITE(DBG%FHNDL,'(A10,I10,A20,I10,A15,2I10)') 'MYRANK', MYRANK, 'STATION =',I, 'IN ELEMENT =', STATION(I)%ELEMENT, STATION(I)%IFOUND
            CALL FLUSH(DBG%FHNDL)
          END IF
        END DO
        allocate(rbuf_int(1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 41')
        DO I = 1, IOUTS
          CALL MPI_REDUCE(STATION(I)%IFOUND,STATION(I)%ISUM,1, itype,MPI_SUM,0,COMM,IERR)
          IF (myrank == 0) THEN
            WRITE(DBG%FHNDL,'(A30,3I10)') 'SUM OF THE FOUND STATIONS MYRANK', MYRANK, I, STATION(I)%ISUM
            CALL FLUSH(DBG%FHNDL)
            rbuf_int(1)=STATION(I)%ISUM
            DO iProc=2,nproc
              CALL MPI_SEND(rbuf_int,1,itype, iProc-1, 144, COMM, ierr)
            ENDDO
          ELSE
            CALL MPI_RECV(rbuf_int,1,itype, 0, 144, COMM, istatus, ierr)
            STATION(I)%ISUM=rbuf_int(1)
          END IF
        END DO
        deallocate(rbuf_int)
        IF (myrank == 0) THEN
          DO I = 1, IOUTS
            IF (STATION(I)%ISUM .EQ. 0) THEN
              WRITE(DBG%FHNDL,'(A20,I10,A10,2F15.8)') 'STATION NOT FOUND', I, STATION(I)%NAME, STATION(I)%XCOORD, STATION(I)%YCOORD
            ELSE
              WRITE(DBG%FHNDL,'(A25,I10,A10,2F15.8)') 'STATION FOUND    ', I, STATION(I)%NAME, STATION(I)%XCOORD, STATION(I)%YCOORD
            END IF
          END DO
        END IF

        ALLOCATE (DEPLOC_SUM(IOUTS), WKLOC_SUM(IOUTS,MSC), CURTXYLOC_SUM(IOUTS,2), ACLOC_SUM(IOUTS,MSC,MDC), USTAR_SUM(IOUTS), ALPHA_SUM(IOUTS), WINDY_SUM(IOUTS), WINDX_SUM(IOUTS), Z0_SUM(IOUTS), CD_SUM(IOUTS), WATLEVLOC_SUM(IOUTS), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 42')
        ACLOC_SUM           = 0.
        WKLOC_SUM           = 0.
        DEPLOC_SUM          = 0.
        CURTXYLOC_SUM       = 0.
        USTAR_SUM           = 0.
        CD_SUM              = 0.
        Z0_SUM              = 0.
        WINDX_SUM           = 0.
        WINDY_SUM           = 0.
        WATLEVLOC_SUM = 0.
#else
        XYTMP(1,:) = XP
        XYTMP(2,:) = YP
        IF (DIMMODE .EQ. 2) THEN
          WRITE(STAT%FHNDL,*) 'FINDING ELEMENT CONNECTED TO STATION'
          DO I = 1, IOUTS
            CALL FIND_ELE ( MNE,MNP,INE,XYTMP,STATION(I)%XCOORD, STATION(I)%YCOORD,STATION(I)%ELEMENT )
            IF (STATION(I)%ELEMENT == 0) THEN
              STATION(I)%IFOUND = 0
              WRITE(STAT%FHNDL,*) STATION(I)%NAME, STATION(I)%XCOORD, STATION(I)%YCOORD, ' is out of mesh !'
            ELSE
              STATION(I)%IFOUND = 1
              NI = INE(:,STATION(I)%ELEMENT)
              STATION(I)%XELE(:) = XP(NI)
              STATION(I)%YELE(:) = YP(NI)
              CALL INTELEMENT_COEF(XP(NI),YP(NI), STATION(I)%XCOORD,STATION(I)%YCOORD, STATION(I)%WI)
              WRITE(STAT%FHNDL,*)'Site    ',STATION(I)%NAME, STATION(I)%XCOORD,STATION(I)%YCOORD,STATION(I)%IFOUND
            END IF
          END DO
        END IF
#endif
        STATION(:)%ISMAX = MSC
        IF (LSIGMAX) THEN
          DO IP = 1, IOUTS
            DO IS = 1, MSC
              IF (SPSIG(IS)/PI2 .GT. STATION(IP)%CUTOFF) THEN
                STATION(IP)%ISMAX = IS - 1
                EXIT
              END IF
            END DO
            WRITE(DBG%FHNDL,*) 'CUT-OFF FREQ. OF STATION =', IP, STATION(IP)%CUTOFF, 'RAD - IS =', STATION(IP)%ISMAX
          END DO
        END IF
      END IF
      WRITE(STAT%FHNDL,'("+TRACE...",A)')'FINISHED WITH INIT_STATION_OUTPUT'
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE TERMINATE_STATION_OUTPUT
      USE DATAPOOL
      IMPLICIT NONE
         IF (LOUTS) THEN
           DEALLOCATE (ACLOC_STATIONS)
           DEALLOCATE (CDLOC_STATIONS)
           DEALLOCATE (Z0LOC_STATIONS)
           DEALLOCATE (ALPHALOC_STATIONS)
           DEALLOCATE (WINDXLOC_STATIONS)
           DEALLOCATE (WINDYLOC_STATIONS)
           DEALLOCATE (USTARLOC_STATIONS)
           DEALLOCATE (DEPLOC_STATIONS)
           DEALLOCATE (WKLOC_STATIONS)
           DEALLOCATE (CURTXYLOC_STATIONS)
           DEALLOCATE (WATLEVLOC_STATIONS)
#ifdef MPI_PARALL_GRID
           DEALLOCATE (DEPLOC_SUM)
           DEALLOCATE (WKLOC_SUM)
           DEALLOCATE (CURTXYLOC_SUM)
           DEALLOCATE (ACLOC_SUM)
           DEALLOCATE (USTAR_SUM, ALPHA_SUM)
           DEALLOCATE (WINDY_SUM, WINDX_SUM)
           DEALLOCATE (Z0_SUM, CD_SUM, WATLEVLOC_SUM)
#endif
         END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READWAVEPARWWM()
#ifdef MPI_PARALL_GRID
      USE elfe_glbl, only: ipgl
      USE elfe_msgp
#endif
      USE DATAPOOL

      IMPLICIT NONE

      INTEGER         :: IP
#ifdef MPI_PARALL_GRID
      INTEGER         :: IPP
      REAL(rkind)     :: RTMP
#endif

!     SPPARM(1): Hs, sign. wave height
!     SPPARM(2): Wave period given by user (either peak or mean)
!     SPPARM(3): average direction
!     SPPARM(4): directional spread
!     SPPARM(5): spectral shape (1-4),
!                (1 - Pierson-Moskowitz,
!                 2 - JONSWAP,
!                 3 - BIN,
!                 4 - Gauss)
!                     negative peak (+) 
!                     or mean frequency (-)
!     SPPARM(6): directional spreading in degree (1) or exponent (2)
!     SPPARM(7): gaussian width for the gauss spectrum 0.1
!     SPPARM(8): peak enhancement factor for the JONSWAP spectra 3.3

      IF (LINHOM) THEN
        READ(WAV%FHNDL,*)
      END IF

#ifdef MPI_PARALL_GRID
      IPP = 0
      IF (LINHOM) THEN
        DO IP = 1, IWBMNPGL
          IF(ipgl(IWBNDGL(IP))%rank == myrank) THEN ! if boundary nodes belong to local domain ...
            IPP = IPP + 1
            READ (WAV%FHNDL, *) SPPARM(:,IPP) ! ... read values into boundary array
          ELSE
            READ (WAV%FHNDL, *) RTMP, RTMP, RTMP, RTMP, RTMP, RTMP, RTMP, RTMP ! ... else ... throw them away
          ENDIF
        END DO
      ELSE
        READ (WAV%FHNDL, *) SPPARM(:,1)
        DO IP = 2, IWBMNPGL
          SPPARM(:,IP) = SPPARM(:,1)
        END DO
      END IF
#else 
      IF (LINHOM) THEN
        DO IP = 1, IWBMNP
          READ (WAV%FHNDL, *) SPPARM(:,IP)
        END DO
      ELSE
        READ (WAV%FHNDL, *) SPPARM(:,1)
      END IF
#endif

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
!AR: check this for wave boundary nodes ....
      SUBROUTINE READWAVEPARFVCOM()
#ifdef MPI_PARALL_GRID
         USE elfe_glbl, only: ipgl
         USE elfe_msgp
#endif
      USE DATAPOOL

      IMPLICIT NONE
 
      INTEGER         :: IP
#ifdef MPI_PARALL_GRID
      REAL(rkind)     :: RTMP
      INTEGER         :: IPP
#endif
!     SPPARM(1): Hs, sign. wave height
!     SPPARM(2): Wave period given by user (either peak or mean)
!     SPPARM(3): average direction
!     SPPARM(4): directional spread
!     SPPARM(5): spectral shape (1-4),
!                (1 - Pierson-Moskowitz,
!                 2 - JONSWAP,
!                 3 - BIN,
!                 4 - Gauss)
!                     negative peak (+) 
!                     or mean frequency (-)
!     SPPARM(6): directional spreading in degree (1) or exponent (2)
!     SPPARM(7): gaussian width for the gauss spectrum 0.1
!     SPPARM(8): peak enhancement factor for the JONSWAP spectra 3.3

      SPPARM(4,:) = 20.
      SPPARM(5,:) = 2.
      SPPARM(6,:) = 2.
      SPPARM(7,:) = 0.1
      SPPARM(8,:) = 3.3

      IF (LINHOM) THEN
        READ(WAV%FHNDL,*)
      END IF

#ifdef MPI_PARALL_GRID
      IPP = 0
      IF (LINHOM) THEN
        DO IP = 1, IWBMNPGL
          IF(ipgl(IWBNDGL(IP))%rank == myrank) THEN ! IF boundary nodes belong to local domain read values into boundary array
            IPP = IPP + 1
            READ (WAV%FHNDL, *) SPPARM(1,IPP), SPPARM(2,IPP), SPPARM(3,IPP)
          ELSE
            READ (WAV%FHNDL, *) RTMP, RTMP, RTMP ! ELSE ... throw them away ...
          ENDIF
        END DO
      ELSE
        READ (WAV%FHNDL, *) SPPARM(1,1), SPPARM(2,1), SPPARM(3,1)
        DO IP = 1, IWBMNPGL
          SPPARM(1:3,IP) = SPPARM(1:3,1)
        END DO
      END IF
#else
      IF (LINHOM) THEN
        DO IP = 1, IWBMNP
          READ (WAV%FHNDL, *) SPPARM(1,IP), SPPARM(2,IP), SPPARM(3,IP)
        END DO
      ELSE
        READ (WAV%FHNDL, *) SPPARM(1,1), SPPARM(2,1), SPPARM(3,1)
        DO IP = 1, IWBMNP
          SPPARM(1:3,IP) = SPPARM(1:3,1)
        END DO
      END IF
#endif

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READSPEC1D(LFIRST)
      USE DATAPOOL
      IMPLICIT NONE

      INTEGER              :: IP, IS
      LOGICAL, INTENT(IN)  :: LFIRST
!
!     Read Spectrum 1-D File ...
!     Second Line ... number of frequencies
!     Third  Line ... number of directions
!
      IF (LFIRST) THEN
        OPEN(WAV%FHNDL, FILE = TRIM(WAV%FNAME), STATUS = 'OLD')   
        READ (WAV%FHNDL,*) WBMSC
        READ (WAV%FHNDL,*) WBMDC
        !WRITE(*,*) WBMSC, WBMDC
        CALL ALLOC_SPEC_BND
        DO IS = 1, WBMSC
          READ (WAV%FHNDL,*) SFRQ(IS,1)
        END DO
      END IF

      IF (LINHOM) THEN
        DO IP = 1, IWBMNP
          DO IS = 1, WBMSC
            READ (WAV%FHNDL, *) SPEG(IS,1,IP), SDIR(IS,1), SPRD(IS,1)
          END DO
        END DO
      ELSE
        DO IS = 1, WBMSC
          READ (WAV%FHNDL, *) SPEG(IS,1,1), SDIR(IS,1), SPRD(IS,1)
          !WRITE(*,*) SPEG(IS,1,1), SDIR(IS,1), SPRD(IS,1)
        END DO
      END IF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READSPEC2D
      USE DATAPOOL
      IMPLICIT NONE

      END SUBROUTINE
!**********************************************************************
!* READSPEC2D_WW3INIT1
!* Read the header of a WAVEWATCHIII binary spectral file to get 
!* dimensions of the spectral grid and number of output locations 
!* required for dynamic allocation.
!*
!* Called by GETWW3SPECTRA
!*
!* Authors: Aron Roland 
!*          Kai Li
!*          Guillaume Dodet (18/12/2012)  
!**********************************************************************
      SUBROUTINE READSPEC2D_WW3_INIT_SPEC
      USE DATAPOOL
      IMPLICIT NONE
      CHARACTER(LEN=30) :: GNAME
      CHARACTER(LEN=21) :: LABEL
        OPEN(WAV%FHNDL,FILE=WAV%FNAME, STATUS='OLD',CONVERT='BIG_ENDIAN',FORM='UNFORMATTED')
        READ(WAV%FHNDL)LABEL, MSC_WW3,MDC_WW3, NP_WW3, GNAME
        WRITE(STAT%FHNDL,*) 'START READSPEC2D_WW3_INIT_SPEC'
        WRITE(STAT%FHNDL,*) 'LABEL, MSC_WW3,MDC_WW3, NP_WW3, GNAME'
        WRITE(STAT%FHNDL,*) LABEL, MSC_WW3,MDC_WW3, NP_WW3, GNAME
        CLOSE(WAV%FHNDL)
        WRITE(STAT%FHNDL,*)'DIRECTION NUMBER IN WW3 SPECTRUM:',MDC_WW3
        WRITE(STAT%FHNDL,*)'FREQUENCY NUMBER IN WW3 SPECTRUM:',MSC_WW3
        WRITE(STAT%FHNDL,'("+TRACE...",A)')'DONE READSPEC2D_WW3_INIT_SPEC'
      END SUBROUTINE
!**********************************************************************
!* READSPEC2D_WW3INIT2
!* Read the header of a WAVEWATCHIII binary spectral file to get
!* frequencies and directions of the spectral grid and starting time
!*
!* Called by GETWW3SPECTRA
!*
!* Authors: Aron Roland
!*          Kai Li
!*          Guillaume Dodet (18/12/2012)
!**********************************************************************
      SUBROUTINE READSPEC2D_WW3_INIT_TIME 
      USE DATAPOOL
        IMPLICIT NONE
        REAL                 :: SPECDMP(MSC_WW3,MDC_WW3)
        INTEGER              :: TMP, IFLAG, IP, TIME(2), IT
        INTEGER, ALLOCATABLE :: ITIME(:,:)
        CHARACTER(LEN=30)    :: GNAME
        CHARACTER(LEN=21)    :: LABEL
        CHARACTER(LEN=10)    :: PID
        CHARACTER(LEN=20)    :: CTIME1, CTIME2
        CHARACTER(LEN=15)    :: TIMESTRING

        REAL :: TMP1(MSC_WW3),TMP2(MDC_WW3) !GD: in ww3 binary file, reals 
        REAL :: TMPR1, TMPR2, TMPR3, TMPR4, TMPR5, TMPR6, TMPR7
        integer istat

        WRITE(STAT%FHNDL,*)'START READSPEC2D_WW3_INIT_TIME'
   
        ALLOCATE(FQ_WW3(MSC_WW3), DR_WW3(MDC_WW3), XP_WW3(NP_WW3), YP_WW3(NP_WW3), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 45')
        FQ_WW3=ZERO
        DR_WW3=ZERO
        XP_WW3=ZERO
        YP_WW3=ZERO
 
        OPEN(WAV%FHNDL, FILE = WAV%FNAME, STATUS = 'OLD',CONVERT='BIG_ENDIAN',FORM='UNFORMATTED')
        READ(WAV%FHNDL)LABEL,TMP,TMP,TMP,GNAME
        READ(WAV%FHNDL)TMP1
        FQ_WW3 = TMP1
        READ(WAV%FHNDL)TMP2
        DR_WW3 = TMP2
        MAXSTEP_WW3 = 0

        DO 
          READ(WAV%FHNDL,IOSTAT=IFLAG) TIME
          IF (IFLAG .GT. 0) THEN
            CALL WWM_ABORT('IFLAG incorrectly set 1')
          ELSE IF (IFLAG .LT. 0) THEN
            WRITE(STAT%FHNDL,*) 'END OF FILE REACHED AT 1, WHICH IS NICE'
            EXIT
          END IF
          DO IP = 1, NP_WW3 
            READ(WAV%FHNDL,IOSTAT=IFLAG) PID,TMPR1,TMPR2,TMPR3,TMPR4,TMPR5,TMPR6,TMPR7
!            WRITE(STAT%FHNDL,'(A10,7F15.4)') PID,TMPR1,TMPR2,TMPR3,TMPR4,TMPR5,TMPR6,TMPR7
            IF (IFLAG .GT. 0) THEN
              CALL WWM_ABORT('IFLAG incorrectly set 2')
            ELSE IF (IFLAG .LT. 0) THEN
              WRITE(STAT%FHNDL,*) 'END OF FILE REACHED AT 2, WHICH IS NOT GOOD'
              EXIT
            END IF
            READ(WAV%FHNDL,IOSTAT=IFLAG) SPECDMP(:,:)
            IF (IFLAG .GT. 0) THEN
              CALL WWM_ABORT('IFLAG incorrectly set 3')
            ELSE IF (IFLAG .LT. 0) THEN
              WRITE(STAT%FHNDL,*) 'END OF FILE REACHED AT 3, WHICH IS NOT GOOD'  
              EXIT
            END IF
          ENDDO
          MAXSTEP_WW3 = MAXSTEP_WW3 + 1
          IF (MAXSTEP_WW3 == 1) TSTART_WW3 = TIME
        END DO

        WRITE(STAT%FHNDL,*) 'NUMBER OF BUOYS', NP_WW3
        WRITE(STAT%FHNDL,*) 'NUMBER OF TIME STEPS IN FILE', MAXSTEP_WW3 
 
        REWIND(WAV%FHNDL)

        ALLOCATE(ITIME(MAXSTEP_WW3,2), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 46')
        ITIME = 0

        READ(WAV%FHNDL)LABEL,TMP,TMP,TMP,GNAME
        READ(WAV%FHNDL)TMP1
        READ(WAV%FHNDL)TMP2
        DO IT = 1, MAXSTEP_WW3 
          READ(WAV%FHNDL,IOSTAT=IFLAG) ITIME(IT,:)
          DO IP = 1, NP_WW3
            READ(WAV%FHNDL,IOSTAT=IFLAG) PID,TMPR1,TMPR2,TMPR3,TMPR4,TMPR5,TMPR6,TMPR7
            READ(WAV%FHNDL,IOSTAT=IFLAG) SPECDMP(:,:)
          ENDDO
        END DO

        ALLOCATE (BND_TIME_ALL_FILES(1,MAXSTEP_WW3), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 47')
        BND_TIME_ALL_FILES = ZERO

        DO IT = 1, MAXSTEP_WW3 
          WRITE(CTIME1,*) ITIME(IT,1) 
          CALL LEADINGZERO(ITIME(IT,2),CTIME2)
          TIMESTRING = ADJUSTL(TRIM(CTIME1)//'.'//ADJUSTL(CTIME2))
          CALL CT2MJD(TIMESTRING, BND_TIME_ALL_FILES(1,IT))
          WRITE(STAT%FHNDL,*) IT, TIMESTRING, BND_TIME_ALL_FILES(1,IT)
        END DO

        DTBOUND_WW3 = (BND_TIME_ALL_FILES(1,2)-BND_TIME_ALL_FILES(1,1))*DAY2SEC
        SEBO%DELT = DTBOUND_WW3
        SEBO%BMJD = BND_TIME_ALL_FILES(1,1) 
        SEBO%EMJD = BND_TIME_ALL_FILES(1,MAXSTEP_WW3) 

        CLOSE(WAV%FHNDL)
        WRITE(STAT%FHNDL,*)'MIN. FREQ. IN WW3 SPECTRUM:',FQ_WW3(1)
        WRITE(STAT%FHNDL,*)'MAX. FREQ. IN WW3 SPECTRUM:',FQ_WW3(MSC_WW3)
        WRITE(STAT%FHNDL,*)'NUMBER OF TIME STEPS',MAXSTEP_WW3
        WRITE(STAT%FHNDL,*)'TIME INCREMENT IN SPECTRAL FILE', DTBOUND_WW3
        WRITE(STAT%FHNDL,*)'FIRST TIME STEP IN WW3 SPECTRUM FILE:',BND_TIME_ALL_FILES(1,1)
        WRITE(STAT%FHNDL,*)'BEGING TIME, END TIME and DELT of wave boundary', SEBO%BMJD, SEBO%EMJD, SEBO%DELT
        WRITE(STAT%FHNDL,*)'BEGING TIME, END TIME and DELT of simulation', MAIN%BMJD, MAIN%EMJD, MAIN%DELT
        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'DONE READSPEC2D_WW3INIT2'

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE LEADINGZERO(INTIN,STRING)
! AR: some fortran magic with leading zero's ...
      INTEGER, INTENT(IN) :: INTIN
      CHARACTER(LEN=6), INTENT(OUT) :: STRING
        write( string, '(I6)' ) INTIN 
        string = repeat('0',6-len_trim(adjustl(string)))//adjustl(string)
      END SUBROUTINE
!**********************************************************************
!* READSPEC2D_WW3
!* Reads spectra in WAVEWATCHIII binary spectral file 
!*
!* Called by GETWW3SPECTRA
!*
!* Authors: Aron Roland
!*          Kai Li
!*          Guillaume Dodet (18/12/2012)
!*
!* Remarks: GD: Need to be modified for time interpolation. Input
!*              arguments should include date for interpolation. If
!*              constant forcing is required, this time can be set to 
!*              zero and only the first step is read.  
!* Remakrs: AR: At this time the whole file is read but this should be 
!*              replaced by a direct access call to the binary file 
!*              Gulliaume can you do this? It is not that urgent.
!**********************************************************************
      SUBROUTINE READ_SPEC_WW3(ISTEP,SPECOUT)
      USE DATAPOOL
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: ISTEP

      INTEGER :: TMP
      INTEGER :: IP, IT, TIME(2)

      REAL(rkind), INTENT(OUT) :: SPECOUT(MSC_WW3,MDC_WW3,NP_WW3)

      REAL :: SPECOUT_SGLE(MSC_WW3,MDC_WW3)
      REAL :: FQ_WW3_SNGL(MSC_WW3), DR_WW3_SNGL(MDC_WW3)
      REAL :: XP_WW3_SGLE(NP_WW3), YP_WW3_SGLE(NP_WW3)
      REAL :: D(NP_WW3),UA(NP_WW3),UD(NP_WW3),CA(NP_WW3),CD2(NP_WW3)
      REAL :: m0,m1,m2,df
      INTEGER :: is,id

      CHARACTER(LEN=30) :: GNAME
      CHARACTER(LEN=21) :: LABEL
      CHARACTER(LEN=10) :: PID(NP_WW3)

      OPEN(WAV%FHNDL, FILE = WAV%FNAME, STATUS = 'OLD',CONVERT='BIG_ENDIAN',FORM='UNFORMATTED')
      READ(WAV%FHNDL) LABEL,TMP,TMP,TMP,GNAME
      READ(WAV%FHNDL) FQ_WW3_SNGL
      READ(WAV%FHNDL) DR_WW3_SNGL
      IF(LBCSE) THEN ! non-stationary ...
        DO IT=1,MAXSTEP_WW3
          READ(WAV%FHNDL) TIME
          DO IP = 1, NP_WW3
            READ(WAV%FHNDL)PID(IP),YP_WW3_SGLE(IP),XP_WW3_SGLE(IP),D(IP),UA(IP),UD(IP),CA(IP),CD2(IP) ! As if XP and YP would change in time ... well i leave it as it is ... 
            YP_WW3(IP) = YP_WW3_SGLE(IP)*DEGRAD
            XP_WW3(IP) = XP_WW3_SGLE(IP)*DEGRAD
            READ(WAV%FHNDL)SPECOUT_SGLE(:,:)
!            write(*,*) sum(SPECOUT_SGLE)
            SPECOUT(:,:,IP) = SPECOUT_SGLE
            m0 = 0.; m1 = 0.; m2 = 0.
            DO ID = 1,MDC_WW3-1
              DO IS = 1,MSC_WW3-1
                DF = FQ_WW3(IS+1)-FQ_WW3(IS)
                M0 = M0+((SPECOUT_SGLE(IS+1,ID)+SPECOUT_SGLE(IS,ID))/TWO)*DF*DDIR_WW3
                M1 = M1+((FQ_WW3(IS+1)-FQ_WW3(IS))/TWO)*((SPECOUT_SGLE(IS+1,ID)+SPECOUT_SGLE(IS,ID))/TWO)*DF*DDIR_WW3
                M2 = M2+(((FQ_WW3(IS+1)-FQ_WW3(IS))/TWO)**TWO)*((SPECOUT_SGLE(IS+1,ID)+SPECOUT_SGLE(IS,ID))/TWO)*DF*DDIR_WW3
!                write(*,*) 4*sqrt(m0), m1, m2, df
              ENDDO
            ENDDO
          ENDDO ! IP
          IF (IT == ISTEP) EXIT ! Read the certain timestep indicated by ISTEP ...
        ENDDO ! IT 
      ELSE ! stationary ... 
        READ(WAV%FHNDL) TIME
        DO IP = 1, NP_WW3
          READ(WAV%FHNDL)PID(IP),YP_WW3_SGLE(IP),XP_WW3_SGLE(IP),D(IP), UA(IP),UD(IP),CA(IP),CD2(IP)
          YP_WW3(IP) = YP_WW3_SGLE(IP)*DEGRAD
          XP_WW3(IP) = XP_WW3_SGLE(IP)*DEGRAD
          READ(WAV%FHNDL)SPECOUT_SGLE(:,:)
          SPECOUT(:,:,IP) = SPECOUT_SGLE
        END DO
      ENDIF
      CLOSE(WAV%FHNDL)
      WRITE(STAT%FHNDL,'("+TRACE...",A)') 'DONE READSPEC2D_WW3'
      END SUBROUTINE
!**********************************************************************
!* GETWW3SPECTRA
!* Read a WAVEWATCHIII binary spectral file and do time and space 
!* interpolation if required.
!*
!* Called by WAVE_BOUNDARY_CONDITIONS
!*
!* Authors: Aron Roland
!*          Kai Li
!*          Guillaume Dodet (18/12/2012)
!**********************************************************************
        SUBROUTINE GET_BINARY_WW3_SPECTRA (ISTEP,WBACOUT)

        USE DATAPOOL
        IMPLICIT NONE

        INTEGER, INTENT(IN)      :: ISTEP
        REAL(rkind), INTENT(OUT) :: WBACOUT(MSC,MDC,IWBMNP)


        INTEGER     :: I,J,IB,IPGL,IBWW3,STEP,TIME(2),IS
        REAL(rkind) :: SPEC_WW3(MSC_WW3,MDC_WW3,NP_WW3),SPEC_WWM(MSC,MDC,NP_WW3)
        REAL(rkind) :: DIST(NP_WW3),TMP(NP_WW3), INDBWW3(NP_WW3)
        REAL(rkind) :: SPEC_WW3_TMP(MSC_WW3,MDC_WW3,NP_WW3),SPEC_WW3_UNSORT(MSC_WW3,MDC_WW3,NP_WW3)
        REAL(rkind) :: JUNK(MDC_WW3),DR_WW3_UNSORT(MDC_WW3),DR_WW3_TMP(MDC_WW3)
        REAL(rkind) :: XP_WWM,YP_WWM
!
! Read spectra in file
!
        CALL READ_SPEC_WW3(ISTEP,SPEC_WW3_UNSORT)

        !DO IBWW3=1,NP_WW3
        !  WRITE(STAT%FHNDL,*) 'ORIG WW3 SUM SPEC', IBWW3, SUM(SPEC_WW3_UNSORT(:,:,IBWW3))
        !END DO

!
! Sort directions and carries spectra along (ww3 directions are not
! montonic)
!
        DO IBWW3 = 1, NP_WW3
          DO IS = 1,MSC_WW3
            DR_WW3_TMP=DR_WW3
            SPEC_WW3_TMP(IS,:,IBWW3) = SPEC_WW3_UNSORT(IS,:,IBWW3)
            CALL SSORT2 (DR_WW3_TMP,SPEC_WW3_TMP(IS,:,IBWW3),JUNK,MDC_WW3,2)
            SPEC_WW3(IS,:,IBWW3) = SPEC_WW3_TMP(IS,:,IBWW3)
            DR_WW3 = DR_WW3_TMP
          ENDDO
          DDIR_WW3 = DR_WW3(2) - DR_WW3(1)
!          WRITE(STAT%FHNDL,*) 'AFTER SORTING', IBWW3, SUM(SPEC_WW3(:,:,IBWW3))
        ENDDO ! IBWW3
!
! Interpolate ww3 spectra on wwm frequency grid
! GD: at the moment 360º spanning grids are mandatory
!
        IF((FQ_WW3(1).GT.FRLOW).OR.(FQ_WW3(MSC_WW3).LT.FRHIGH)) THEN
!          WRITE(STAT%FHNDL,*)'WW3 FMIN = ',FQ_WW3(1),'WWM FMIN = ',FRLOW
!          WRITE(STAT%FHNDL,*)'WW3 FMAX = ',FQ_WW3(MSC_WW3),'WWM FMAX = ', FRHIGH
!          WRITE(STAT%FHNDL,*)'WW3 spectra does not encompass the whole WWM spectra, please carefully check if this makes sense for your simulations'
          CALL SPECTRALINT(SPEC_WW3,SPEC_WWM)
        ELSE
!          WRITE(STAT%FHNDL,*)'WW3 FMIN = ',FQ_WW3(1),'WWM FMIN = ',FRLOW
!          WRITE(STAT%FHNDL,*)'WW3 FMAX = ',FQ_WW3(MSC_WW3),'WWM FMAX = ', FRHIGH
          CALL SPECTRALINT(SPEC_WW3,SPEC_WWM)
        ENDIF

        IF (LNANINFCHK) THEN
          WRITE(DBG%FHNDL,*) 'SUMS AFTER INTERPOLATION', SUM(SPEC_WW3),SUM(SPEC_WWM)
        ENDIF
!
! Interpolate ww3 spectra on wwm boundary nodes
! GD: ww3 forcing works until here. Some more debugging is needed
!
        TMP = 0
        IF(LINHOM) THEN !nearest-neighbour interpolation
          DO IB=1,IWBMNP
            IPGL = IWBNDLC(IB)
#ifdef SELFE
            XP_WWM=XLON(IPGL)! * RADDEG 
            YP_WWM=YLAT(IPGL)! * RADDEG 
#else
            XP_WWM = XP(IPGL)
            YP_WWM = YP(IPGL)
#endif
            IF (NP_WW3 .GT. 1) THEN
              DO IBWW3=1,NP_WW3
                !WRITE(STAT%FHNDL,*)'XP_WWM =',XP_WWM,'XP_WW3 =',XP_WW3(IBWW3)
                !WRITE(STAT%FHNDL,*)'YP_WWM =',YP_WWM,'YP_WW3 =',YP_WW3(IBWW3)
                DIST(IBWW3)=SQRT((XP_WWM-XP_WW3(IBWW3))**2+(YP_WWM-YP_WW3(IBWW3))**2)
                INDBWW3(IBWW3)=IBWW3
                !WRITE(STAT%FHNDL,*) 'orig', IBWW3, INDBWW3(IBWW3), DIST(IBWW3)
              ENDDO
              CALL SSORT2 (DIST, INDBWW3, TMP, NP_WW3, 2)
              !DO IBWW3=1,NP_WW3
              !  WRITE(STAT%FHNDL,*) 'sorted', IBWW3, INDBWW3(IBWW3), DIST(IBWW3)
              !END DO
              CALL SHEPARDINT2D(2, 1./DIST(1:2),MSC,MDC,SPEC_WWM(:,:,INT(INDBWW3(1:2))), WBACOUT(:,:,IB), 1)
              !WRITE(STAT%FHNDL,'(A20, 2F20.5,3F30.10)') ' AFTER INTERPOLATION ', INDBWW3(1), INDBWW3(2), sum(SPEC_WWM(:,:,INT(INDBWW3(1)))), sum(SPEC_WWM(:,:,INT(INDBWW3(2)))), SUM(WBACOUT(:,:,IB))
            ELSE
              WBACOUT(:,:,IB) = SPEC_WWM(:,:,1)
            ENDIF
          ENDDO ! IB 
        ELSE
          DO IB=1,IWBMNP
            WBACOUT(:,:,IB) = SPEC_WWM(:,:,1) 
          ENDDO
        ENDIF

        !DO IB = 1, IWBMNP
        !  WRITE(STAT%FHNDL,*) 'SUM OF WBAC', IB, SUM(WBACOUT(:,:,IB)) 
        !ENDDO 

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'DONE GETWW3SPECTRA'
        END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
        SUBROUTINE INIT_BINARY_WW3_SPECTRA 

        USE DATAPOOL
        IMPLICIT NONE

        INTEGER :: I,J
!
! Read header to get grid dimension, frequencies and directions, 
! output locations and first time step
!
        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'START WITH INIT_BINARY_WW3_SPECTRA'

        CALL READSPEC2D_WW3_INIT_SPEC
        CALL READSPEC2D_WW3_INIT_TIME

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'DONE WITH INIT_BINARY_WW3_SPECTRA'
        END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
        SUBROUTINE SHEPARDINT2D(NP,WEIGHT,D1,D2,Z,ZINT,P)
        USE DATAPOOL
        IMPLICIT NONE
!AR: Kai Li, please carefully describe the method and comment on the input parameters ...

        INTEGER, INTENT(IN) :: NP,P,D1,D2
        REAL(rkind), INTENT(IN)  :: WEIGHT(NP), Z(D1,D2,NP)
        REAL(rkind), INTENT(OUT) :: ZINT(D1,D2)
        INTEGER :: IP
        REAL :: SW

        SW=0
        ZINT=0
        DO IP=1,NP
          SW=SW+WEIGHT(IP)**P
          ZINT(:,:)=ZINT(:,:)+Z(:,:,IP)*(WEIGHT(IP)**P)
        ENDDO
        ZINT=ZINT/SW

        END SUBROUTINE
!**********************************************************************
!* SPECTRALINT
!* Interpolate spectrum on WWM spectral grid.
!*
!* Called by GETWW3SPECTRA
!*
!* Authors: Aron Roland
!*          Kai Li
!*          Guillaume Dodet (18/12/2012)
!*
!* Remarks: GD)The spectral grid needs to be define over 360º and
!*             the min (resp. max) frequencies need to be smaller
!*             (resp. higher) than in WWM frequency grid.
!**********************************************************************
        SUBROUTINE SPECTRALINT(SPEC_WW3,SPEC_WWM)
        USE DATAPOOL
        IMPLICIT NONE

        REAL(rkind), INTENT(IN)  :: SPEC_WW3(MSC_WW3,MDC_WW3,NP_WW3)
        REAL(rkind), INTENT(OUT) :: SPEC_WWM(MSC,MDC,NP_WW3)

        REAL(rkind) :: SPEC_WW3_TMP(MSC_WW3,MDC,NP_WW3),tmp(msc)
        REAL(rkind) :: DF, M0_WW3, M1_WW3, M2_WW3, HSWW3, M0_WWM, M1_WWM, M2_WWM

        INTEGER     :: IP,IS,ID

        REAL(rkind) :: JACOBIAN(MSC), AM, SM, OUTPAR(OUTVARS)

        WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ENTERING SPECTRALINT'

        JACOBIAN = ONE/(SPSIG*PI2)! ENERGY / HZ -> ACTION / RAD
 
        SPEC_WW3_TMP = ZERO
        SPEC_WWM     = ZERO

        IF (LNANINFCHK) THEN
          WRITE(DBG%FHNDL,'(A20,I10,3F30.2)') 'BEFORE INTERPOLATION', IP, SUM(SPEC_WW3), SUM(SPEC_WW3_TMP), SUM(SPEC_WWM) 
        ENDIF

        DO IP=1,NP_WW3
          DO IS=1,MSC_WW3
            CALL INTERLIND(MDC_WW3,MDC,DR_WW3,SPDIR,SPEC_WW3(IS,:,IP),SPEC_WW3_TMP(IS,:,IP))
          ENDDO
          DO ID=1,MDC 
            CALL INTERLIN (MSC_WW3,MSC,FQ_WW3,FR,SPEC_WW3_TMP(:,ID,IP),SPEC_WWM(:,ID,IP))
          ENDDO
          M0_WW3 = ZERO; M1_WW3 = ZERO; M2_WW3 = ZERO
          DO ID = 1,MDC_WW3
            DO IS = 1,MSC_WW3-1
              DF = FQ_WW3(IS+1)-FQ_WW3(IS)
              AM = (SPEC_WW3(IS+1,ID,IP)+SPEC_WW3(IS,ID,IP))/TWO
              SM = (FQ_WW3(IS+1)+FQ_WW3(IS))/TWO
              M0_WW3 =M0_WW3+AM*DF*DDIR_WW3
              M1_WW3 =M1_WW3+AM*SM*DF*DDIR_WW3
              M2_WW3 =M2_WW3+AM*SM**2*DF*DDIR_WW3
            ENDDO
          ENDDO
          M0_WWM = ZERO; M1_WWM = ZERO; M2_WWM = ZERO 
          DO ID = 1,MDC
            DO IS = 1,MSC-1
              DF = FR(IS+1)-FR(IS)
              AM = (SPEC_WWM(IS+1,ID,IP)+SPEC_WWM(IS,ID,IP))/TWO
              SM = (FR(IS+1)+FR(IS))/TWO
              M0_WWM =M0_WWM+AM*DF*DDIR
              M1_WWM =M1_WWM+AM*SM*DF*DDIR
              M2_WWM =M2_WWM+AM*SM**2*DF*DDIR
            ENDDO
          ENDDO
          WRITE(STAT%FHNDL,*) 'POINT NUMBER', IP
          WRITE(STAT%FHNDL,'(A10,2F20.10,A10,2F20.10)') 'M1 = ',M1_WW3, M1_WWM, 'M2 = ',M2_WW3, M2_WWM
        END DO

        IF (LNANINFCHK) THEN
          WRITE(DBG%FHNDL,'(A20,I10,3F30.2)') 'AFTER INTERPOLATION', IP, SUM(SPEC_WW3), SUM(SPEC_WW3_TMP), SUM(SPEC_WWM)
        ENDIF

! Do jacobian

        DO IP = 1, NP_WW3
          DO ID = 1, MDC
            SPEC_WWM(:,ID,IP) = SPEC_WWM(:,ID,IP) * JACOBIAN(:) ! convert to wave action in sigma,theta space
          END DO              
        END DO

        WRITE(STAT%FHNDL,*)'CHECKING INTEGRATED PARAMETERS AFTER JACOBIAN'
        DO IP = 1, NP_WW3
          M0_WW3 = ZERO; M1_WW3 = ZERO; M2_WW3 = ZERO
          DO ID = 1,MDC_WW3
            DO IS = 1,MSC_WW3-1
              DF = FQ_WW3(IS+1)-FQ_WW3(IS)
              AM = (SPEC_WW3(IS+1,ID,IP)+SPEC_WW3(IS,ID,IP))/TWO
              SM = (FQ_WW3(IS+1)+FQ_WW3(IS))/TWO
              M0_WW3 =M0_WW3+AM*DF*DDIR_WW3
              M1_WW3 =M1_WW3+AM*SM*DF*DDIR_WW3
              M2_WW3 =M2_WW3+AM*SM**2*DF*DDIR_WW3
            ENDDO
          ENDDO
          M0_WWM = ZERO; M1_WWM = ZERO; M2_WWM = ZERO
          DO ID = 1,MDC
            DO IS = 1,MSC-1
              DF = SPSIG(IS+1)-SPSIG(IS)
              SM = (SPSIG(IS+1)+SPSIG(IS))/TWO
              AM = (SPEC_WWM(IS+1,ID,IP)+SPEC_WWM(IS,ID,IP))/TWO * SM
              M0_WWM =M0_WWM+AM*DF*DDIR
              M1_WWM =M1_WWM+AM*SM*DF*DDIR
              M2_WWM =M2_WWM+AM*SM**2*DF*DDIR
            ENDDO
          ENDDO
          WRITE(STAT%FHNDL,*) 'POINT NUMBER', IP
          WRITE(STAT%FHNDL,'(A10,2F20.10,A10,2F20.10)') 'M1 = ',M1_WW3, M1_WWM, 'M2 = ',M2_WW3, M2_WWM
        END DO

        IF (LNANINFCHK) THEN
          WRITE(DBG%FHNDL,'(A20,I10,3F30.2)') 'AFTER JACOBIAN', IP, SUM(SPEC_WW3), SUM(SPEC_WW3_TMP), SUM(SPEC_WWM)
        ENDIF

        END SUBROUTINE
!**********************************************************************
!* INTERLIND
!* Interpolate vector on a 2-pi periodic axis (directions in wwm grid).
!*
!* Called by SPECTRALINT
!*
!* Authors: Aron Roland
!*          Kai Li
!*          Guillaume Dodet (18/12/2012)
!*
!* Remarks: GD)This routine should be togeher with INTERLIN either here 
!*             or in wwm_aux.F90.
!**********************************************************************
      SUBROUTINE INTERLIND (NX1, NX2, X1, X2, Y1, Y2)
      USE DATAPOOL, ONLY : THR, rkind,PI2
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NX1, NX2

      REAL(rkind), INTENT(IN)    :: X1(NX1), Y1(NX1)
      REAL(rkind), INTENT(IN)    :: X2(NX2)
      REAL(rkind), INTENT(OUT)   :: Y2(NX2)

      INTEGER             :: I, J
      REAL(rkind)         :: DX1(NX1)

      DO I = 1, NX1 - 1
        DX1(I) = X1(I+1)-X1(I)
      END DO
      DX1(NX1) = X1(1)+PI2-X1(NX1)

      DO I = 1, NX2
        DO J = 1,NX1 - 1
          IF (ABS(X2(I)-X1(J)) .LT. THR) THEN
            Y2(I) = Y1(J)
            EXIT
          ELSE IF (X2(I) .GT. X1(J) .AND. X2(I) .LT. X1(J+1)) THEN
            Y2(I) = Y1(J)+(Y1(J+1)-Y1(J))/DX1(J)*(X2(I)-X1(J))
            EXIT
          ELSE IF (X2(I) .GT. X1(NX1)) THEN
            Y2(I) = Y1(NX1)+(Y1(1)-Y1(NX1))/DX1(NX1)*(X2(I)-X1(NX1))
            EXIT
          ELSE IF (X2(I) .LT. X1(1)) THEN
            Y2(I) = Y1(NX1)+(Y1(1)-Y1(NX1))/DX1(NX1)*(X2(I)+PI2-X1(NX1))
            EXIT            
          END IF
        END DO
      END DO

      END SUBROUTINE INTERLIND
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE  ALLOC_SPEC_BND()
      USE DATAPOOL
      IMPLICIT NONE
      integer istat
      IF (LINHOM) THEN
        ALLOCATE (SFRQ(WBMSC,IWBMNP), SPEG(WBMSC,WBMDC,IWBMNP), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 48')
        ALLOCATE (SDIR(WBMSC,IWBMNP), SPRD(WBMSC,IWBMNP), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 49')
      ELSE
        ALLOCATE (SFRQ(WBMSC,1), SPEG(WBMSC,WBMDC,1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 50')
        ALLOCATE (SDIR(WBMSC,1), SPRD(WBMSC,1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_initio, allocate error 51')
      END IF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READWAVEPARWW3()
#ifdef MPI_PARALL_GRID
         USE elfe_glbl, only: ipgl
         USE elfe_msgp
#endif
      USE DATAPOOL

      IMPLICIT NONE

      INTEGER         :: IP
#ifdef MPI_PARALL_GRID
      INTEGER         :: IPP
      REAL(rkind)     :: RTMP
#endif

!     SPPARM(1): Hs, sign. wave height
!     SPPARM(2): Wave period given by user (either peak or mean)
!     SPPARM(3): average direction
!     SPPARM(4): directional spread
!     SPPARM(5): spectral shape (1-4),
!                (1 - Pierson-Moskowitz,
!                 2 - JONSWAP,
!                 3 - BIN,
!                 4 - Gauss)
!                     negative peak (+) 
!                     or mean frequency (-)
!     SPPARM(6): directional spreading in degree (1) or exponent (2)
!     SPPARM(7): gaussian width for the gauss spectrum 0.1
!     SPPARM(8): peak enhancement factor for the JONSWAP spectra 3.3

      SPPARM(4,:) = 20.
      SPPARM(5,:) = 2.
      SPPARM(6,:) = 2.
      SPPARM(7,:) = 0.1
      SPPARM(8,:) = 3.3

      IF (LINHOM) THEN
        READ(WAV%FHNDL,*)
      END IF

#ifdef MPI_PARALL_GRID
      IPP = 0
      IF (LINHOM) THEN
        DO IP = 1, IWBMNPGL
          IF(ipgl(IWBNDGL(IP))%rank == myrank) THEN ! IF boundary nodes belong to local domain read values into boundary array
            IPP = IPP + 1
            READ (WAV%FHNDL, *) SPPARM(1,IPP), SPPARM(2,IPP), SPPARM(3,IPP)
          ELSE
            READ (WAV%FHNDL, *) RTMP, RTMP, RTMP ! ELSE ... throw them away ...
          ENDIF
        END DO
      ELSE
        READ (WAV%FHNDL, *) SPPARM(1,1), SPPARM(2,1), SPPARM(3,1)
        DO IP = 1, IWBMNPGL
          SPPARM(1:3,IP) = SPPARM(1:3,1)
        END DO
      END IF
#else
      IF (LINHOM) THEN
        DO IP = 1, IWBMNP
          READ (WAV%FHNDL, *) SPPARM(1,IP), SPPARM(2,IP), SPPARM(3,IP)
        END DO
      ELSE
        READ (WAV%FHNDL, *) SPPARM(1,1), SPPARM(2,1), SPPARM(3,1)
        DO IP = 1, IWBMNP
          SPPARM(1:3,IP) = SPPARM(1:3,1)
        END DO
      END IF
#endif

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_WAVE_BOUNDARY_CONDITION(IFILE, IT)
        USE DATAPOOL
        IMPLICIT NONE

        INTEGER, INTENT(OUT) :: IFILE, IT

        REAL(rkind)            :: DTMP
        INTEGER                :: ITMP
        CHARACTER(len=25)     :: CHR

!TODO: Makes sure initial condition work also when no wave boundary is set ...

        IF (LBCWA .OR. LBCSP) THEN
          IF (IBOUNDFORMAT == 3) THEN
            IF (LBCSP) THEN
              CALL INIT_BINARY_WW3_SPECTRA
              DTMP = (MAIN%BMJD-BND_TIME_ALL_FILES(1,1)) * DAY2SEC
              IT   = NINT(DTMP/SEBO%DELT) + 1
              CHR = 'FROM INIT_WAVE_BOUNDARY 1'
              CALL WAVE_BOUNDARY_CONDITION(1,IT,WBAC,CHR)
              IF (LBINTER) WBACOLD = WBAC
              WRITE(STAT%FHNDL,*) 'INITIALIZING WAVE BOUNDARY IT =', IT
              WRITE(STAT%FHNDL,*) 'SUM OF WAVE ACTION', SUM(WBACOLD), SUM(WBAC) 
            ELSE IF (LBCWA) THEN 
#ifdef NCDF
              CALL INIT_NETCDF_WW3_WAVEPARAMETER
#else
              CALL WWM_ABORT('Compile with NCDF For WW3 bdcons')
#endif
              DTMP = (MAIN%BMJD-BND_TIME_ALL_FILES(1,1)) * DAY2SEC
              ITMP  = 0
              DO IFILE = 1, NUM_NETCDF_FILES_BND
                ITMP = ITMP + NDT_BND_FILE(IFILE)
                IF (ITMP .GT. INT(DTMP/SEBO%DELT)) EXIT
              END DO
              ITMP = SUM(NDT_BND_FILE(1:IFILE-1))
              IT   = NINT(DTMP/SEBO%DELT) - ITMP + 1
              IF (IT .GT. NDT_BND_FILE(IFILE)) THEN
                IFILE = IFILE + 1
                IT    = 1
              ENDIF
              WRITE(STAT%FHNDL,*) IFILE, IT, SUM(NDT_BND_FILE(1:IFILE-1)), NINT(DTMP/SEBO%DELT), SEBO%DELT
              CHR = 'FROM INIT_WAVE_BOUNDARY 2'
              CALL WAVE_BOUNDARY_CONDITION(IFILE,IT,WBAC,CHR)
              IF (LBINTER) WBACOLD = WBAC
            END IF
          ELSE ! BOUNDFORMAT
            CHR = 'FROM INIT_WAVE_BOUNDARY 4'
            CALL WAVE_BOUNDARY_CONDITION(1,1,WBAC,CHR)
            IF (LBINTER) WBACOLD = WBAC
          END IF
        END IF

      END SUBROUTINE INIT_WAVE_BOUNDARY_CONDITION
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_HMAX
        USE DATAPOOL

        HMAX = BRHD * DEP

        IF (LMONO_IN) HMAX = HMAX * SQRT(2.)

      END SUBROUTINE SET_HMAX
!**********************************************************************
!*                                                                    *
!**********************************************************************
