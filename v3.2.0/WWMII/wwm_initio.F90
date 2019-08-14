!     Last change:  1    10 Mar 2004   11:54 pm
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INITIALIZE_WWM()
         USE DATAPOOL
#ifdef PETSC
#ifdef SELFE
         USE PETSCPOOL, ONLY : PETSC_INIT_PARALLEL
#else
         USE PETSCPOOL, ONLY : PETSC_INIT
#endif
#endif
#ifdef SELFE
         use elfe_msgp!, only : myrank,parallel_abort,itype,comm,ierr
#endif
         IMPLICIT NONE
#ifdef SELFE
         include 'mpif.h'
#endif

         REAL    :: TIME1, TIME2
         INTEGER :: IP, ITMP, IT, IFILE, IPGL
         REAL*8  :: DTMP

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

#ifdef SELFE
         DEP  = REAL(DEP8)
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

#ifdef WWMONLY
         CALL READ_SPATIAL_GRID
         WLDEP = DEP ! Set reference depth ...
         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'READ SPATIAL GRID'
         CALL FLUSH(STAT%FHNDL)
#endif
         CALL SPATIAL_GRID
         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'INIT SPATIAL GRID'
         CALL FLUSH(STAT%FHNDL)

         IF (LADVTEST) THEN
           ALLOCATE(UTEST(MNP)); UTEST = 0.
           CALL ADVTEST(UTEST)
           AC2(:,1,1) = UTEST
           CALL CHECKCONS(UTEST,SUMACt0)
         END IF

         CALL SET_DEP
         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'SET DEPTH POINTER'
         CALL FLUSH(STAT%FHNDL)

         IF (DIMMODE .EQ. 2) THEN
           WRITE(STAT%FHNDL,'("+TRACE...",A)') 'THE FLUCTUATION SPLITTING PREPROCESSOR HAS STARTED'
           CALL INIT_FLUCT_ARRAYS
           CALL INIT_FLUCT

          IF (AMETHOD .EQ. 4) THEN
#ifdef SELFE
#ifdef PETSC
            CALL PETSC_INIT_PARALLEL
#endif PETSC
#else
#ifdef PETSC
            CALL PETSC_INIT
#endif PETSC
#endif
          END IF

          WRITE(STAT%FHNDL,'("+TRACE...",A)') 'THE FLUCTUATION SPLITTING PREPROCESSOR HAS ENDED'
          CALL FLUSH(STAT%FHNDL)
         END IF

         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'INITIALIZE SPECTRAL GRID'
         CALL FLUSH(STAT%FHNDL)
         CALL SPECTRAL_GRID

#ifdef WWMONLY
         CALL SET_IOBP
#elif SELFE
         CALL SET_IOBP_SELFE
         DMIN = DMIN_SELFE
#endif 
         CALL SET_IOBPD

         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'INITIALIZE WIND CURRENT WATERLEVEL'
         CALL FLUSH(STAT%FHNDL)
         CALL INIT_WIND_CURRENT_WATLEV()

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
         ELSE IF (MESIN .EQ. 2) THEN
           WRITE(STAT%FHNDL,'("+TRACE...",A)') 'INIT CYCLE4 WIND INPUT'
           CALL PREPARE_SOURCE
           CALL FLUSH(STAT%FHNDL)
         END IF

         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'SET THE INITIAL WAVE BOUNDARY CONDITION'
         CALL FLUSH(STAT%FHNDL)
         CALL INIT_WAVE_BOUNDARY(IFILE,IT)

         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'SET THE INITIAL CONDITION'
         CALL FLUSH(STAT%FHNDL)
         CALL INITIAL_CONDITION(IFILE,IT)

         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'SET BOUNDARY CONDITIONS'
         CALL FLUSH(STAT%FHNDL)
         CALL SET_WAVE_BOUNDARY

         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'Writing Initial Condition to Output and initialize Station Output'
         CALL FLUSH(STAT%FHNDL)
         CALL INIT_STATION_OUTPUT

         CALL OUTPUT(0.,.TRUE.)
         OUTF%TMJD = OUTF%TMJD + DBLE(OUTF%DELT*SEC2DAY)

         IF (LWXFN) THEN
           CALL WRINPGRD
         END IF

#ifdef WWMONLY
         IF (LCPL) THEN
           WRITE(STAT%FHNDL,'("+TRACE...",A)') 'OPEN PIPES FOR COUPLING'
           CALL FLUSH(STAT%FHNDL)
           IF (LTIMOR) THEN
             CALL INIT_PIPES_TIMOR()
           ELSE IF (LSHYFEM) THEN
             CALL INIT_PIPES_SHYFEM()
           ELSE IF (LROMS) THEN
             CALL INIT_PIPES_ROMS()
           END IF
         END IF
#endif

#ifdef SELFE
         IF (size(OUTT_INTPAR(1,:)) .NE. OUTVARS) THEN
           call parallel_abort('OUTPUTSIZE IS NOT RIGHT CHECK OUTVARS')
         END IF
         IF (size(WIND_INTPAR(1,:)) .NE. WINDVARS) THEN
           call parallel_abort('OUTPUTSIZE IS NOT RIGHT CHECK VARS')
         END IF
#endif
         CALL CPU_TIME(TIME2)

#ifdef SELFE
         IF (MSC_SELFE .NE. MSC .OR. MDC_SELFE .NE. MDC) THEN
           WRITE(*,*) 'MSC_SELFE', MSC_SELFE
           WRITE(*,*) 'MSC', MSC
           WRITE(*,*) 'MDC_SELFE', MDC_SELFE
           WRITE(*,*) 'MDC', MDC
           CALL PARALLEL_ABORT('THERE IS AND ERROR IN MSC2 OR MDC2 IN PARAM.IN')
         END IF
#endif

         WRITE(STAT%FHNDL,'("+TRACE...",A,F15.4)') 'CPU Time for the preprocessing', TIME2-TIME1
         CALL FLUSH(STAT%FHNDL)

       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_PIPES_TIMOR()
         USE DATAPOOL
         IMPLICIT NONE
!
! open pipe data files for coupling
!
           LSEWL = .TRUE.
           LSECU = .TRUE.
           WRITE(*,'("+TRACE...",A)') 'OPEN PIPE'
!          Pipes that are read by the wave model
           OPEN(1000,file='p_velx.dat'  ,form='unformatted', action='read')
           OPEN(1001,file='p_vely.dat'  ,form='unformatted', action='read')
           OPEN(1002,file='p_lev.dat'   ,form='unformatted', action='read')
           OPEN(1003,file='p_bot.dat'   ,form='unformatted', action='read')
           OPEN(1004,file='p_time.dat'  ,form='unformatted', action='write')
           OPEN(1005,file='p_dthyd.dat',form='unformatted', action='read')
!          Pipes that are written by the wave modell
           OPEN(101 ,file='p_stressx.dat' ,form='unformatted', action='write')
           OPEN(102 ,file='p_stressy.dat' ,form='unformatted', action='write')
           OPEN(103 ,file='p_waveh.dat'   ,form='unformatted', action='write')
           OPEN(104 ,file='p_wavet.dat'   ,form='unformatted', action='write')
           OPEN(105 ,file='p_waved.dat'   ,form='unformatted', action='write')
           OPEN(106 ,file='p_wavekm.dat'  ,form='unformatted', action='write')
           OPEN(107 ,file='p_wavetp.dat'  ,form='unformatted', action='write')
           OPEN(108 ,file='p_wavekp.dat'  ,form='unformatted', action='write')
           OPEN(109 ,file='p_orbit.dat'   ,form='unformatted', action='write')
           OPEN(110 ,file='p_stokesx.dat' ,form='unformatted', action='write')
           OPEN(111 ,file='p_stokesy.dat' ,form='unformatted', action='write')
           OPEN(112 ,file='p_windx.dat'   ,form='unformatted', action='write')
           OPEN(113 ,file='p_windy.dat'   ,form='unformatted', action='write')
!          Pipes writen as ergzus.bin for XFN
           OPEN(2003, file='pipe_wave.bin' ,form='unformatted' ,status = 'unknown')
           OPEN(2004, file='pipe_orbi.bin' ,form='unformatted' ,status = 'unknown')
           OPEN(2005, file='pipe_stok.bin' ,form='unformatted' ,status = 'unknown')
           WRITE(*,'("+TRACE...",A)') 'END OPEN PIPE'
!
! write coupling tyme step
!
           WRITE(1004) MAIN%DTCOUP
           CALL FLUSH(1004)
!
! read coupling tyme step
!
           READ(1005)  MAIN%DTCUR

           IF(ABS(MAIN%DTCOUP-INT(MAIN%DTCOUP/MAIN%DTCUR)*MAIN%DTCUR).GT. .001) THEN
             write(*,*) 'TIME STEP OF THE hydraulic flow MODEL CANNOT BE DIVIDIED WITHOUT A REST'
             write(*,*)'dt Stroemung (s) =',MAIN%DTCUR, ',  dt Kopplung (s) = ',MAIN%DTCOUP
!             STOP 'dt Stroemungsmodell muss gerades Vielfaches des Kopplungs-dt sein'
           END IF

           WRITE(*,'("+TRACE... DTCUR and DTCOUP",A)') MAIN%DTCUR, MAIN%DTCOUP

       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_PIPES_SHYFEM()

         USE DATAPOOL
         IMPLICIT NONE
!
! open pipe data files for coupling
!
           LSEWL = .TRUE.
           LSECU = .TRUE.

           WRITE(*,'("+TRACE...",A)') 'OPEN PIPE'
!          Pipes that are read by the wave model
           OPEN(1000,file='p_velx.dat'    ,form='unformatted', action='read')
           OPEN(1001,file='p_vely.dat'    ,form='unformatted', action='read')
           OPEN(1002,file='p_lev.dat'     ,form='unformatted', action='read')
           OPEN(1003,file='p_bot.dat'     ,form='unformatted', action='read')
           OPEN(1004,file='p_zeta3d.dat'  ,form='unformatted', action='read')
!          Pipes that are written by the wave model
           OPEN(101 ,file='p_stressx.dat' ,form='unformatted', action='write')
           OPEN(102 ,file='p_stressy.dat' ,form='unformatted', action='write')
           OPEN(142 ,file='p_stresxy.dat' ,form='unformatted', action='write')  !ccf
           OPEN(103 ,file='p_waveh.dat'   ,form='unformatted', action='write')
           OPEN(104 ,file='p_wavet.dat'   ,form='unformatted', action='write')
           OPEN(105 ,file='p_waved.dat'   ,form='unformatted', action='write')
           OPEN(106 ,file='p_wavekm.dat'  ,form='unformatted', action='write')
           OPEN(107 ,file='p_wavetp.dat'  ,form='unformatted', action='write')
           OPEN(108 ,file='p_wavekp.dat'  ,form='unformatted', action='write')
           OPEN(109 ,file='p_orbit.dat'   ,form='unformatted', action='write')
           OPEN(110 ,file='p_stokesx.dat' ,form='unformatted', action='write')
           OPEN(111 ,file='p_stokesy.dat' ,form='unformatted', action='write')

           IF (.NOT. LWINDWWM) THEN
!            *** WIND FROM SHYFEM *** ccf
             OPEN(112,file='p_windx.dat'  ,form='unformatted', action='read')
             OPEN(113,file='p_windy.dat'  ,form='unformatted', action='read')
            ELSE
!            *** WIND FROM WWM *** ccf
             OPEN(112 ,file='p_windx.dat',form='unformatted', action='write')
             OPEN(113 ,file='p_windy.dat',form='unformatted', action='write')
           END IF

           WRITE(*,'("+TRACE...",A)') 'END OPEN PIPE'

       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_PIPES_ROMS()
         USE DATAPOOL
         IMPLICIT NONE
!
! open pipe data files for coupling
!
! Mathieu can you please edit this part ... or I can do it during this weekend just tell me ... :)
           LSEWL = .TRUE.
           LSECU = .TRUE.
           WRITE(*,'("+TRACE...",A)') 'OPEN PIPE ROMS'
!          Pipes that are read by the wave model
           OPEN(1000,file='pipe/ExchRW'  ,form='unformatted', action='read')
           Print *, 'WWM: open pipe ExchImport'
!          Pipes that are written by the wave modell
           OPEN(101 ,file='pipe/ExchWR' ,form='unformatted', action='write')
           Print *, 'WWM: open pipe ExchExport'
           WRITE(*,'("+TRACE...",A)') 'END OPEN PIPE ROMS'

       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_ARRAYS()

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER :: IP

#ifdef SELFE
         INTEGER :: IE
#endif

         ALLOCATE( DX1(0:MNP+1) ); DX1 = 0.
         ALLOCATE( DX2(0:MNP+1) ); DX2 = 0.
         ALLOCATE( XP(MNP) )
         ALLOCATE( YP(MNP) )
         ALLOCATE( INVSPHTRANS(MNP,2)); INVSPHTRANS = 0.
         ALLOCATE( DEP(MNP) ); DEP = 0.
         ALLOCATE( INE(3,MNE) ); INE = 0
         ALLOCATE( IEN(6,MNE) );  IEN = 0.d0 
         ALLOCATE( TRIA(MNE) ); TRIA = 0.

#ifdef SELFE
         DO IE = 1, MNE
           INE(:,IE) = INETMP(IE,:)
         END DO
#elif WWMONLY
         XP  = 0.
         YP  = 0.
#endif
!
! spectral grid - shared
!
         ALLOCATE( SPSIG(MSC) ); SPSIG = 0.
         ALLOCATE( SPDIR(MDC) ); SPDIR = 0.
         ALLOCATE( FR(MSC) ); FR    = 0.
         ALLOCATE( COSTH(MDC) ); COSTH = 0.
         ALLOCATE( SINTH(MDC) ); SINTH = 0.
         ALLOCATE( COS2TH(MDC) ); COS2TH = 0.
         ALLOCATE( SIN2TH(MDC) ); SIN2TH = 0.
         ALLOCATE( SINCOSTH(MDC) ); SINCOSTH = 0.
         ALLOCATE( SIGPOW(MSC,6) ); SIGPOW = 0.
         ALLOCATE( DS_BAND(0:MSC+1) ); DS_BAND = 0.
         ALLOCATE( DS_INCR(0:MSC+1) ); DS_INCR = 0.
!
! action densities and source terms - shared
!
         ALLOCATE (AC2(MNP,MSC,MDC)); AC2 = 0.

         !IF (ICOMP .GE. 2) THEN
           ALLOCATE (AC1(MNP,MSC,MDC)); AC1 = 0.
         IF (ICOMP .GE. 2) THEN
           ALLOCATE (IMATRAA(MNP,MSC,MDC)); IMATRAA = 0.
           ALLOCATE (IMATDAA(MNP,MSC,MDC)); IMATDAA = 0.
         END IF

         IF (LITERSPLIT) THEN
           ALLOCATE (AC1(MNP,MSC,MDC)); AC1 = 0.
           ALLOCATE (DAC_ADV(2,MNP,MSC,MDC)); DAC_ADV = 0.d0
           ALLOCATE (DAC_THE(2,MNP,MSC,MDC)); DAC_THE = 0.d0
           ALLOCATE (DAC_SIG(2,MNP,MSC,MDC)); DAC_SIG = 0.d0
           ALLOCATE (DAC_SOU(2,MNP,MSC,MDC)); DAC_SOU = 0.d0
         END IF

         IF (LSHYFEM) THEN
           ALLOCATE(SHYFZETA(NLVT,MNP),NLEV(MNP))
         END IF
!
! WAM Cycle 4.5 - shared
!
         IF (MESIN .EQ. 2) THEN
           ALLOCATE ( TAUHFT(0:IUSTAR,0:IALPHA) ); TAUHFT   = 0.
           ALLOCATE ( TAUT  (0:ITAUMAX,0:JUMAX) ); TAUT = 0.
         ENDIF
!
! Boundary conditions - shared
!
         ALLOCATE( IOBPD(MDC,MNP) ); IOBPD  = 1
         ALLOCATE( IOBP(MNP) ); IOBP   = 0
         ALLOCATE( IOBWB(MNP) ); IOBWB  = 1
!
! phase velocity, wave number, group velocity, dwdh, kh
!
         ALLOCATE( WK(MNP,MSC) ); WK = 0.
         ALLOCATE( CG(MNP,MSC) ); CG = 0.

#ifdef SELFE
!         ALLOCATE( CGX(MNP,MSC,MDC) ); CGX = 0.
!         ALLOCATE( CGY(MNP,MSC,MDC) ); CGY = 0.
#endif
!
         ALLOCATE( TABK (0:IDISPTAB) ); TABK = 0.
         ALLOCATE( TABCG(0:IDISPTAB) ); TABCG = 0.
!
! diffraction parameter - shared
!
         IF (LDIFR) THEN
           ALLOCATE ( DIFRM(MNP), DIFRX(MNP), DIFRY(MNP) )
           DIFRM = 0.
           DIFRX = 0.
           DIFRY = 0.
         END IF
!
! water level, currents and depths ...
!
         ALLOCATE( WINDXY(MNP,2) ); WINDXY = 0.
         ALLOCATE( PRESSURE(MNP) ); PRESSURE = 0.
         ALLOCATE( DVWIND(MNP,2) ); DVWIND = 0.
         ALLOCATE( CURTXY(MNP,2) ); CURTXY = 0.
         ALLOCATE( DVCURT(MNP,2) ); DVCURT = 0.
         ALLOCATE( DDEP(MNP,2) ); DDEP = 0.
         ALLOCATE( DCUX(MNP,2) ); DCUX = 0.
         ALLOCATE( DCUY(MNP,2) ); DCUY = 0. 
         ALLOCATE( WATLEV(MNP) ); WATLEV = 0.
         ALLOCATE( WATLEVOLD(MNP) ); WATLEVOLD = 0.
         ALLOCATE( DVWALV(MNP) ); DVWALV = 0.
         ALLOCATE( WLDEP(MNP) ); WLDEP = 0.
         ALLOCATE( DEPDT(MNP) ); DEPDT = 0.
!
!  convergence analysis - shared
!
         IF (LCONV .OR. (LQSTEA .AND. LCHKCONV)) THEN
           ALLOCATE ( SUMACOLD(MNP), HSOLD(MNP), KHSOLD(MNP), TM02OLD(MNP) )
           SUMACOLD = 0.
           HSOLD  = 0.
           TM02OLD = 0.
           KHSOLD  = 0.
         END IF
!
!  output - shared
!
         ALLOCATE( QBLOCAL(MNP) ); QBLOCAL = 0.
         ALLOCATE( DISSIPATION(MNP) ); DISSIPATION = 0.
         ALLOCATE( AIRMOMENTUM(MNP) ); AIRMOMENTUM = 0.
         ALLOCATE( UFRIC(MNP) ); UFRIC = 0.
         ALLOCATE( ALPHA_CH(MNP) ); ALPHA_CH = 0.
         ALLOCATE( TAUW(MNP) ); TAUW = 0.
         ALLOCATE( TAUTOT(MNP) ); TAUTOT = 0.
         ALLOCATE( TAUWX(MNP) ); TAUWX = 0.
         ALLOCATE( TAUWY(MNP) ); TAUWY = 0.
         ALLOCATE( TAUHF(MNP) ); TAUHF = 0.
         ALLOCATE( Z0(MNP) ); Z0 = 0.
         ALLOCATE( CD(MNP) ); CD = 0.
         ALLOCATE( USTDIR(MNP) ); USTDIR = 0.
         ALLOCATE( RSXX(MNP)); RSXX = 0.
         ALLOCATE( RSXY(MNP)); RSXY = 0.
         ALLOCATE( RSYY(MNP)); RSYY = 0.
!
#ifdef SELFE
!
         ALLOCATE( SXX3D(NVRT,MNP), SXY3D(NVRT,MNP), SYY3D(NVRT,MNP) )
         SXX3D = 0.
         SXY3D = 0.
         SYY3D = 0.
!
#elif WWMONLY
!
         IF (LCPL) THEN
           ALLOCATE( SXX3D(NLVT,MNP), SXY3D(NLVT,MNP), SYY3D(NLVT,MNP) )
           SXX3D = 0.
           SXY3D = 0.
           SYY3D = 0.
         END IF
!
#endif
!
         ALLOCATE(HMAX(MNP)); HMAX = 0.
         ALLOCATE(ISHALLOW(MNP)); ISHALLOW = 0.

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_FLUCT_ARRAYS()

         USE DATAPOOL
         IMPLICIT NONE

         ALLOCATE( CCON(MNP) ); CCON = 0
         ALLOCATE( SI(MNP) ); SI = 0.d0
         ALLOCATE( ITER_EXP(MSC,MDC) ); ITER_EXP = 0
         ALLOCATE( ITER_EXPD(MSC) ); ITER_EXPD = 0
         IF (ICOMP .GE. 1) THEN
           ALLOCATE( I_DIAG(MNP) ); I_DIAG = 0
         END IF
         IF (LCFL) THEN
           ALLOCATE (CFLCXY(3,MNP))
           CFLCXY = 0.
         END IF

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE BASIC_PARAMETER()
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER :: IP
         REAL    :: PPTAIL
!
!2do: Check with WAM parameterization and WW3
!
         PQUAD(1) = 0.25
         IF (MESIN .EQ. 1) THEN
           PQUAD(2) = 2.5E7
         ELSE
           PQUAD(2) = 2.78E7
         END IF
         PQUAD(3) = 5.5
         PQUAD(4) = 6.0/7.0
         PQUAD(5) = -1.25
         PQUAD(6) = 0.0
!
! High frequency tail integration factors as defined in the SWAN model.
!
         PTAIL(1) = 4.
         PTAIL(2) = 2.5
         PTAIL(3) = PTAIL(1) + 1.
         IF (MESIN .EQ. 2) THEN
           PTAIL(1) = 5.
           PTAIL(3) = PTAIL(1) + 1.
         ENDIF 
         PTAIL(4) = 3.
         DO IP = 0, 3
           PPTAIL = PTAIL(1) - REAL(IP)
           PTAIL(5+IP) = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
         ENDDO
!
! Output Variables ... LVARS, NVARS
!
         OUTT_VARNAMES(1)  = 'HS'
         OUTT_VARNAMES(2)  = 'TM01'
         OUTT_VARNAMES(3)  = 'TM02'
         OUTT_VARNAMES(4)  = 'KLM'
         OUTT_VARNAMES(5)  = 'WLM'
         OUTT_VARNAMES(6)  = 'ETOTS'
         OUTT_VARNAMES(7)  = 'ETOTC'
         OUTT_VARNAMES(8)  = 'DM'
         OUTT_VARNAMES(9)  = 'DSPR'
         OUTT_VARNAMES(10) = 'TPPD'
         OUTT_VARNAMES(11) = 'TPP'
         OUTT_VARNAMES(12) = 'CPP'
         OUTT_VARNAMES(13) = 'WNPP'
         OUTT_VARNAMES(14) = 'CGPP'
         OUTT_VARNAMES(15) = 'KPP'
         OUTT_VARNAMES(16) = 'LPP'
         OUTT_VARNAMES(17) = 'PEAKD'
         OUTT_VARNAMES(18) = 'PEAKDSPR'
         OUTT_VARNAMES(19) = 'DPEAK'
         OUTT_VARNAMES(20) = 'UBOT'
         OUTT_VARNAMES(21) = 'ORBITAL'
         OUTT_VARNAMES(22) = 'BOTEXPER'
         OUTT_VARNAMES(23) = 'TMBOT'
         OUTT_VARNAMES(24) = 'URSELL'
         OUTT_VARNAMES(25) = 'USTAR'
         OUTT_VARNAMES(26) = 'ALPHA'
         OUTT_VARNAMES(27) = 'Z0'
         OUTT_VARNAMES(28) = 'WIND-X'
         OUTT_VARNAMES(29) = 'WIND-Y'
         OUTT_VARNAMES(30) = 'CD'
         OUTT_VARNAMES(31) = 'CURR-X'
         OUTT_VARNAMES(32) = 'CURR-Y'
         OUTT_VARNAMES(33) = 'DEPTH'
         OUTT_VARNAMES(34) = 'ELEVATION'

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INITIATE_WAVE_PARAMETER()
         USE DATAPOOL
         IMPLICIT NONE

         WRITE(STAT%FHNDL,*) 'START WAVE PARAMETER'
         CALL GRADDEP
         WRITE(STAT%FHNDL,*) 'GRADDEP'
         IF (LSTCU .OR. LSECU) CALL GRADCURT
         WRITE(STAT%FHNDL,*) 'GRADCURT'
         CALL BASIC_PARAMETER
         WRITE(STAT%FHNDL,*) 'BASIC'
         CALL MAKE_WAVE_TABLE
         CALL WAVE_K_C_CG
         WRITE(STAT%FHNDL,*) 'WAVEKCG'
         IF (MESNL .GT. 0) CALL PARAMETER4SNL
         WRITE(STAT%FHNDL,*) 'SNL4'
         !IF (MESTR .GT. 0) CALL INIT_TRIADSWAN
         !WRITE(STAT%FHNDL,*) 'SNL3'

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INITIAL_CONDITION(IFILE,IT)

         USE DATAPOOL
#ifdef NCDF
         USE NETCDF
#endif
#ifdef SELFE
         use elfe_msgp, only : myrank
#endif
         IMPLICIT NONE

         INTEGER :: IP, IS, IFILE, IT
         REAL    :: SPPAR(8)
         REAL    :: MS
         REAL    :: HS, TP, HSLESS, TPLESS, FDLESS
         REAL    :: WIND10, WINDTH, VEC2DEG
         REAL    :: WINDX, WINDY
         REAL    :: ACLOC(MSC,MDC)
         REAL    :: DEG, HOT_FRL, HOT_FRH
         REAL    :: ETOT, TMPPAR(8,MNP), SSBRL(MSC,MDC)
         INTEGER :: iret, ncid, ac_id

         IF (.NOT. LHOTR .AND. LINID) THEN
            IF (INITSTYLE == 2 .AND. IBOUNDFORMAT == 3) THEN
#ifdef NCDF
              CALL READ_NETCDF_WW3(IFILE,IT)
#else
                 STOP 'compile with DNCDF PPFLAG'
#endif
              CALL INTER_STRUCT_DOMAIN(NDX_BND,NDY_BND,DX_BND,DY_BND,OFFSET_X_BND,OFFSET_Y_BND,TMPPAR)
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
                     FDLESS = G9*REAL(AVETL)/WIND10**2
                     HSLESS = MIN(0.21, 0.00288*FDLESS**0.45)
                     TPLESS = MIN(1.0/0.13, 0.46*FDLESS**0.27)
                     HS = HSLESS*WIND10**2/G9
                     TP = TPLESS*WIND10/G9
                     MS = 2.0
                     SPPAR(1) = HS
                     SPPAR(2) = TP ! Check whether TP is inside of spectal representation !!!
                     SPPAR(3) = DEG
                     SPPAR(4) = MS
                     SPPAR(5) = 2.
                     SPPAR(6) = 2.
                     SPPAR(7) = 0.1
                     SPPAR(8) = 3.3
                     CALL SPECTRAL_SHAPE(SPPAR,ACLOC)
                     AC2(IP,:,:) = ACLOC
                   ELSE
                     ACLOC = 1.E-8
                   END IF
                 ELSE IF (INITSTYLE == 2 .AND. IBOUNDFORMAT == 3) THEN
                   CALL SPECTRAL_SHAPE(TMPPAR(:,IP),ACLOC)
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

#ifdef SELFE
         IF (MYRANK == 0) THEN
#endif
           IF (LHOTR) THEN
             IF (HOTSTYLE == 1) THEN

               OPEN(HOTIN%FHNDL, FILE = TRIM(HOTIN%FNAME), STATUS = 'OLD', FORM = 'UNFORMATTED')
               READ(HOTIN%FHNDL,*) HMNP, HMNE
               READ(HOTIN%FHNDL,*) HMSC, HMDC, HFRLOW, HFRHIGH
               IF ( HMNP .NE. MNP .OR. HMNE .NE. MNE .OR. HMSC .NE. MSC .OR. HFRLOW .NE. FRLOW .OR. HFRHIGH .NE. FRHIGH ) THEN
                 CALL ERROR_MSG('THE HOTFILE GEOMETRY DOES NOT FIT THE INPUT FILE','INITIAL_CONDITION')
               ELSE
                 READ(HOTIN%FHNDL,*) AC2
                 CLOSE(HOTIN%FHNDL)
               END IF

             ELSE IF (HOTSTYLE == 2) THEN
#ifdef NCDF
               iret=nf90_open(HOTIN%FNAME, nf90_nowrite, ncid)
               iret=nf90_inq_varid(ncid, "ac", ac_id)
               iret=nf90_get_var(ncid,ac_id,AC2,start = (/1, 1, 1, IHOTPOS/), count = (/MNP, MSC, MDC, 1 /))
               iret=nf90_close(ncid)
#endif
             END IF
           END IF

#ifdef SELFE
         END IF
#endif

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_FILE_HANDLES()
#ifdef SELFE
         use elfe_msgp
#endif
         USE DATAPOOL

         CHARACTER (LEN = 30) :: FDB

            INP%FNAME  = 'wwminput.nml'
            CHK%FNAME  = 'wwmcheck.nml'
           QSTEA%FNAME = 'qstea.out'
         IOBPOUT%FNAME = 'iobp.out'
        IOBPDOUT%FNAME = 'iobpd.out'
!
!2do ... dinstinguish between binary and ascii stuff ...
!
             BND%FHNDL = STARTHNDL + 1
             WIN%FHNDL = STARTHNDL + 2
             CUR%FHNDL = STARTHNDL + 3
             WAT%FHNDL = STARTHNDL + 4
             WAV%FHNDL = STARTHNDL + 5
             CHK%FHNDL = STARTHNDL + 6
           HOTIN%FHNDL = STARTHNDL + 7
          HOTOUT%FHNDL = STARTHNDL + 8
             INP%FHNDL = STARTHNDL + 9
             GRD%FHNDL = STARTHNDL + 10
          GRDCOR%FHNDL = STARTHNDL + 11

           QSTEA%FHNDL = STARTHNDL + 12
        IOBPOUT%FHNDL  = STARTHNDL + 13
        IOBPDOUT%FHNDL = STARTHNDL + 14

           OUT1D%FHNDL = STARTHNDL + 11
            MISC%FHNDL = STARTHNDL + 12
         OUTSP1D%FHNDL = STARTHNDL + 13
         OUTSP2D%FHNDL = STARTHNDL + 14
         OUTPARM%FHNDL = STARTHNDL + 15

#ifdef SELFE
         DBG%FHNDL  = STARTHNDL + 20
         FDB  ='wwmdbg_0000'
         LFDB =len_trim(FDB)
         write(FDB(LFDB-3:LFDB),'(i4.4)') MYRANK
         open(DBG%FHNDL,file='outputs/'//fdb,status='replace') ! errors
         STAT%FHNDL = STARTHNDL + 21 
         FDB  ='wwmstat_0000'
         LFDB =len_trim(FDB)
         write(FDB(LFDB-3:LFDB),'(i4.4)') MYRANK
         open(STAT%FHNDL,file='outputs/'//fdb,status='replace') ! errors
#elif WWMONLY
         DBG%FHNDL  = STARTHNDL + 20
         STAT%FHNDL = STARTHNDL + 21
         open(DBG%FHNDL,file='wwmdbg.out',status='unknown') !non-fatal errors
         open(STAT%FHNDL,file='wwmstat.out',status='unknown') !non-fatal errors
#endif
         OPEN( INP%FHNDL,      FILE = TRIM(INP%FNAME))
         OPEN( CHK%FHNDL,      FILE = TRIM(CHK%FNAME))
         OPEN( QSTEA%FHNDL,    FILE = TRIM(QSTEA%FNAME))
         OPEN( IOBPOUT%FHNDL,  FILE = TRIM(IOBPOUT%FNAME))
         OPEN( IOBPDOUT%FHNDL, FILE = TRIM(IOBPDOUT%FNAME))

!ascii output files ...
           OUT1D%FHNDL = STARTHNDL + 30 
            MISC%FHNDL = STARTHNDL + 31 
         OUTSP1D%FHNDL = STARTHNDL + 32 
         OUTPARM%FHNDL = STARTHNDL + 33 

         OUTSP2D%FHNDL = STARTHNDL + 40 
         OUT%FHNDL     = STARTHNDL + 41 

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_STATION_OUTPUT()
         USE DATAPOOL
#ifdef SELFE
         USE elfe_msgp
         USE elfe_glbl, only : iplg, ielg
#endif 
         IMPLICIT NONE
#ifdef SELFE
         include 'mpif.h'
#endif
         INTEGER           :: I, NI(3), IP, IS
         REAL              :: XYTMP(2,MNP)
         REAL              :: DELTC, DEG
!
!    set the site output
!
#ifdef WWMONLY

         !WRITE(*,*) SIZE(INE(1,:)), SIZE(INE(:,1))
         !WRITE(*,*) MAXVAL(INE), MINVAL(INE)

         IF (LOUTS) THEN
           XYTMP(1,:) = XP
           XYTMP(2,:) = YP
           IF (DIMMODE .EQ. 2) THEN
              WRITE(STAT%FHNDL,*) 'FINDING ELEMENT CONNECTED TO STATION'
              DO I = 1, IOUTS ! Loop over stations ...
                CALL FIND_ELE ( MNE,MNP,INE,XYTMP,STATION(I)%XCOORD,STATION(I)%YCOORD,STATION(I)%ELEMENT ) ! Fine element of output locations ...
                IF (STATION(I)%ELEMENT == 0) THEN
                  STATION(I)%IFOUND = 0
                  WRITE(STAT%FHNDL,*) STATION(I)%NAME, STATION(I)%XCOORD, STATION(I)%YCOORD, ' is out of mesh !'
                ELSE
                  STATION(I)%IFOUND = 1
                  NI = INE(:,STATION(I)%ELEMENT)
                  STATION(I)%XELE(:) = XP(NI)
                  STATION(I)%YELE(:) = YP(NI)
                  WRITE(STAT%FHNDL,*)'Site    ',STATION(I)%NAME,STATION(I)%XCOORD,STATION(I)%YCOORD,STATION(I)%IFOUND
                END IF
             END DO
           END IF

          IF (LOUTS) THEN
            STATION(:)%ISMAX = MSC
            IF (LOUTS) THEN
              IF (LSIGMAX) THEN
                DO IP = 1, IOUTS
                 DO IS = 1, MSC
                   IF (SPSIG(IS)/PI2 .GT. STATION(IP)%CUTOFF) THEN
                     STATION(IP)%ISMAX = IS - 1
                     EXIT
                     END IF
                   END DO
                 END DO
               END IF
             END IF
           END IF

         END IF
!
#elif SELFE
!
         IF (LOUTS) THEN

           XYTMP(1,:) = XP
           XYTMP(2,:) = YP

           IF (DIMMODE .EQ. 2) THEN

              WRITE(DBG%FHNDL,*) 'SEARCHING FOR STATION ACROSS RANKS', myrank 
              DO I = 1, IOUTS ! Loop over stations ...
                CALL FIND_ELE ( MNE,MNP,INE,XYTMP,STATION(I)%XCOORD,STATION(I)%YCOORD,STATION(I)%ELEMENT ) ! Find element of output locations ...
                IF (STATION(I)%ELEMENT .GT. 0) THEN
                  STATION(I)%IFOUND  = 1
                  NI                 = INE(:,STATION(I)%ELEMENT)
                  STATION(I)%XELE(:) = XP(NI)
                  STATION(I)%YELE(:) = YP(NI)
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

              DO I = 1, IOUTS
                CALL MPI_REDUCE(STATION(I)%IFOUND,STATION(I)%ISUM,1,MPI_INTEGER,MPI_SUM,0,COMM,IERR)
                IF (myrank == 0) THEN
                  WRITE(DBG%FHNDL,'(A30,3I10)') 'SUM OF THE FOUND STATIONS MYRANK', MYRANK, I, STATION(I)%ISUM
                  CALL FLUSH(DBG%FHNDL)
                END IF
              END DO

              IF (myrank == 0) THEN
                DO I = 1, IOUTS ! Loop over stations ...
                  IF (STATION(I)%ISUM .EQ. 0) THEN
                    WRITE(DBG%FHNDL,'(A20,I10,A10,2F15.8)') 'STATION NOT FOUND', I, STATION(I)%NAME, STATION(I)%XCOORD, STATION(I)%YCOORD
                  ELSE
                    WRITE(DBG%FHNDL,'(A25,I10,A10,2F15.8)') 'STATION FOUND    ', I, STATION(I)%NAME, STATION(I)%XCOORD, STATION(I)%YCOORD
                  END IF
                END DO
              END IF

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

           ALLOCATE (ACLOC_STATIONS(IOUTS,MSC,MDC), ACLOC_SUM(IOUTS,MSC,MDC))
           ALLOCATE (DEPLOC_STATIONS(IOUTS), WKLOC_STATIONS(IOUTS,MSC), CURTXYLOC_STATIONS(IOUTS,2))
           ALLOCATE (DEPLOC_SUM(IOUTS), WATLEVLOC_STATIONS(IOUTS), WKLOC_SUM(IOUTS,MSC), CURTXYLOC_SUM(IOUTS,2))
           ALLOCATE (USTARLOC_STATIONS(IOUTS), CDLOC_STATIONS(IOUTS), Z0LOC_STATIONS(IOUTS)) 
           ALLOCATE (ALPHALOC_STATIONS(IOUTS), WINDXLOC_STATIONS(IOUTS), WINDYLOC_STATIONS(IOUTS))
           ALLOCATE (USTAR_SUM(IOUTS), WINDY_SUM(IOUTS), WINDX_SUM(IOUTS), ALPHA_SUM(IOUTS))
           ALLOCATE (Z0_SUM(IOUTS), CD_SUM(IOUTS), WATLEVLOC_SUM(IOUTS))

           ACLOC_STATIONS      = 0.
           ACLOC_SUM           = 0.
           DEPLOC_STATIONS     = 0.
           WKLOC_SUM           = 0.
           WKLOC_STATIONS      = 0.
           DEPLOC_SUM          = 0.
           CURTXYLOC_SUM       = 0.
           CURTXYLOC_STATIONS  = 0.
           USTARLOC_STATIONS   = 0.
           USTAR_SUM           = 0.
           CDLOC_STATIONS      = 0.
           CD_SUM              = 0.
           Z0LOC_STATIONS      = 0.
           Z0_SUM              = 0.
           WINDXLOC_STATIONS      = 0.
           WINDX_SUM           = 0.
           WINDYLOC_STATIONS      = 0.
           WINDY_SUM           = 0.
           WATLEVLOC_STATIONS = 0.
           WATLEVLOC_SUM = 0.
         END IF
#endif
!
!     *** output the testing message
!
         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READWAVEPARWWM()
#ifdef SELFE
         USE elfe_glbl, only: ipgl
         USE elfe_msgp
#endif
      USE DATAPOOL

      IMPLICIT NONE

      INTEGER :: IP, IPP, IPPP
      REAL    :: RTMP

!     SPPARM(1): Hs, sign. wave height
!     SPPARM(2): Wave period given by user (either peak or mean)
!     SPPARM(3): average direction
!     SPPARM(4): directional spread
!     SPPARM(5): spectral shape (1-4), (1 - Pierson-Moskowitz, 2 - JONSWAP, 3 - BIN, 4 - Gauss) negative peak (+) or mean frequency (-)
!     SPPARM(6): directional spreading in degree (1) or exponent (2)
!     SPPARM(7): gaussian width for the gauss spectrum 0.1
!     SPPARM(8): peak enhancement factor for the JONSWAP spectra 3.3

      IF (LINHOM) THEN
        READ(WAV%FHNDL,*)
      END IF
#ifdef SELFE
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
#elif WWMONLY
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
#ifdef SELFE
         USE elfe_glbl, only: ipgl
         USE elfe_msgp
#endif
      USE DATAPOOL

      IMPLICIT NONE

      INTEGER :: IP, IPP, IPPP
      REAL    :: RTMP

!     SPPARM(1): Hs, sign. wave height
!     SPPARM(2): Wave period given by user (either peak or mean)
!     SPPARM(3): average direction
!     SPPARM(4): directional spread
!     SPPARM(5): spectral shape (1-4), (1 - Pierson-Moskowitz, 2 - JONSWAP, 3 - BIN, 4 - Gauss) negative peak (+) or mean frequency (-)
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
#ifdef SELFE
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
#elif WWMONLY
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

      INTEGER              :: I, IP, IS, TMP
      LOGICAL, INTENT(IN)  :: LFIRST
!
!     Read Spectrum 1-D File ...
!     Second Line ... number of frequencies
!     Third  Line ... number of directions
!
      IF (LFIRST) THEN
        READ (WAV%FHNDL,*) WBMSC
        READ (WAV%FHNDL,*) WBMDC
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
        END DO
      END IF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READSPEC2D (LFIRST)
      USE DATAPOOL
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: LFIRST

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READSPEC2D_WW3 (LFIRST)
      USE DATAPOOL
      IMPLICIT NONE
      LOGICAL, INTENT(IN) :: LFIRST

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE  ALLOC_SPEC_BND()
      USE DATAPOOL
      IMPLICIT NONE

      IF (LINHOM) THEN
        ALLOCATE (SFRQ(WBMSC,IWBMNP))
        ALLOCATE (SPEG(WBMSC,WBMDC,IWBMNP))
        ALLOCATE (SDIR(WBMSC,IWBMNP))
        ALLOCATE (SPRD(WBMSC,IWBMNP))
      ELSE
        ALLOCATE (SFRQ(WBMSC,1))
        ALLOCATE (SPEG(WBMSC,WBMDC,1))
        ALLOCATE (SDIR(WBMSC,1))
        ALLOCATE (SPRD(WBMSC,1))
      END IF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READWAVEPARWW3()
#ifdef SELFE
         USE elfe_glbl, only: ipgl
         USE elfe_msgp
#endif
      USE DATAPOOL

      IMPLICIT NONE

      INTEGER :: IP, IPP, IPPP
      REAL    :: RTMP

!     SPPARM(1): Hs, sign. wave height
!     SPPARM(2): Wave period given by user (either peak or mean)
!     SPPARM(3): average direction
!     SPPARM(4): directional spread
!     SPPARM(5): spectral shape (1-4), (1 - Pierson-Moskowitz, 2 - JONSWAP, 3 - BIN, 4 - Gauss) negative peak (+) or mean frequency (-)
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
#ifdef SELFE
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
#elif WWMONLY
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
      SUBROUTINE PREPARE_ARDHUIN_OLD()

        USE DATAPOOL
        USE W3SRC4MD_OLD

        IMPLICIT  NONE

        INTEGER :: IK, ISP, ITH, ITH0
        REAL    :: SIGMA, FR1, RTH0

        NK    = MSC
        MK    = NK  ! ?????????????????????????????
        NTH   = MDC
        MTH   = NTH ! ?????????????????????????????
        NSPEC = MSC * MDC
        MSPEC = NSPEC

        ALLOCATE(SIG(0:MSC+1), SIG2(NSPEC), DSIP(0:MSC+1), TH(MDC))
        ALLOCATE(ESIN(MSPEC+MTH), ECOS(MSPEC+MTH), EC2(MSPEC+MTH), ES2(MSPEC+MTH),ESC(MSPEC+MTH) )
        ALLOCATE(DSII(MSC), DDEN(MSC), DDEN2(NSPEC))

        DTH   = DDIR
        FR1   = SPSIG(1)/PI2
        TH    = SPDIR

        RTH0 = 0.
        DO ITH=1, NTH
          TH  (ITH) = DTH * ( RTH0 + REAL(ITH-1) )
          ESIN(ITH) = SIN ( TH(ITH) )
          ECOS(ITH) = COS ( TH(ITH) )
          IF ( ABS(ESIN(ITH)) .LT. 1.E-5 ) THEN
            ESIN(ITH) = 0.
            IF ( ECOS(ITH) .GT. 0.5 ) THEN
              ECOS(ITH) =  1.
            ELSE
              ECOS(ITH) = -1.
              END IF
          END IF
          IF ( ABS(ECOS(ITH)) .LT. 1.E-5 ) THEN
            ECOS(ITH) = 0.
            IF ( ESIN(ITH) .GT. 0.5 ) THEN
              ESIN(ITH) =  1.
            ELSE
              ESIN(ITH) = -1.
            END IF
          END IF
          ES2 (ITH) = ESIN(ITH)**2
          EC2 (ITH) = ECOS(ITH)**2
          ESC (ITH) = ESIN(ITH)*ECOS(ITH)
        END DO

        DO IK=2, NK+1
          ITH0 = (IK-1)*NTH
          DO ITH=1, NTH
            ESIN(ITH0+ITH) = ESIN(ITH)
            ECOS(ITH0+ITH) = ECOS(ITH)
            ES2 (ITH0+ITH) = ES2 (ITH)
            EC2 (ITH0+ITH) = EC2 (ITH)
            ESC (ITH0+ITH) = ESC (ITH)
          END DO
        END DO

        WRITE(5001,*) 'DTH'
        WRITE(5001,*) DTH
        WRITE(5001,*) 'FR1'
        WRITE(5001,*) FR1
        WRITE(5001,*) 'TH'
        WRITE(5001,*) TH
        WRITE(5001,*) 'ESIN, ECOS, EC2'
        WRITE(5001,*) ESIN, ECOS, EC2

!        TAUWSHELTER = 1. ! This maybe even too big ... wave supportesed are ...

        WNMEANP = 0.5
        WNMEANPTAIL = -0.5
        WRITE(5001,*) 'WNMEANP, WNMEANPTAIL'
        WRITE(5001,*) WNMEANP, WNMEANPTAIL

        XFR = EXP(FRINTF) ! Check with Fabrice ... should be 1.1

        SIGMA   = FR1 * TPI / XFR**2 ! What is going on here ?
        SXFR    = 0.5 * (XFR-1./XFR)

        WRITE(5001,*) 'XFR, SIGMA, SXFR'
        WRITE(5001,*) XFR, SIGMA, SXFR

        DO IK=0, NK+1
         SIGMA    = SIGMA * XFR ! What is going on here ...
         SIG (IK) = SIGMA
         DSIP(IK) = SIGMA * SXFR
        END DO

        WRITE(5001,*) 'SIGMA'
        WRITE(5001,*)  SIGMA
        WRITE(5001,*) 'SIG'
        WRITE(5001,*)  SIG
        WRITE(5001,*) 'DSIP'
        WRITE(5001,*)  DSIP

        DSII(1) = 0.5 * SIG( 1) * (XFR-1.)
        DO IK = 2, NK - 1
          DSII(IK) = DSIP(IK)
        END DO
        DSII(NK) = 0.5 * SIG(NK) * (XFR-1.) / XFR

        DDEN = DTH * DSII(:) * SIG(:)

        DO ISP=1, NSPEC
          IK         = 1 + (ISP-1)/NTH
          SIG2 (ISP) = SIG (IK)
          DDEN2(ISP) = DDEN(IK)
        END DO

        WRITE(5001,*) 'SIG2'
        WRITE(5001,*) SIG2
        WRITE(5001,*) 'DSII'
        WRITE(5001,*) DSII
        WRITE(5001,*) 'DDEN'
        WRITE(5001,*) DDEN
        WRITE(5001,*) 'DDEN2'
        WRITE(5001,*) DDEN2

        FTE = 0.25 * SIG(NK) * DTH * SIG(NK)
        FTF = 0.20           * DTH * SIG(NK)

        FACHF  = 5.
        FACHFE = XFR**(-FACHF)

        STXFTFTAIL  = 1./(FACHF-1.-WNMEANPTAIL*2)
        STXFTF      = 1./(FACHF-1.-WNMEANP*2)
        STXFTWN     = 1./(FACHF-1.-WNMEANP*2) * SIG(NK)**(2)

        WRITE(5001,*) 'FTE, FTF, FACHF, FACHFE'
        WRITE(5001,*) FTE, FTF, FACHF, FACHFE

        SSWELLF(1) = 0.8
        SSWELLF(2) = -0.018
        SSWELLF(3) = 0.015
        SSWELLF(4) = 1.E5
        SSWELLF(5) = 1.2
        SSWELLF(6) = 0.
        SSWELLF(7) = 0.

        WRITE(5001,*) 'SSWELLF'
        WRITE(5001,*) SSWELLF

        AALPHA = 0.0095
        BBETA  = 1.31 ! NOMAD for ECMWF 1.54
        ZZALP   = 0.006
        ZZWND   = 10.

        SWELLFPAR = 1
        SSDSTH     = 80.
        SSDSCOS    = 2.
        SSDSBRF1   = 0.5
        SSDSHCK    = 1.
        SSDSC3     = -0.80
        SSDSBCK    = 0.
        SSDSBINT   = 0.3
        SSDSPBK    = 4.
        SSDSABK    = 1.5

        SSDSBR     = 0.90E-3
        SSDSBRFDF  = 0
        SSDSBRF1   = 0.5
        SSDSBRF2   = 0.
        SSDSBR2    = 0.8

        SSDSP      = 2.
        SSDSPBK    = 4.

        SSDSISO     = 2

        SSDSBM(0)  = 1.
        SSDSBM(1)  = 0.
        SSDSBM(2)  = 0.
        SSDSBM(3)  = 0.
        SSDSBM(4)  = 0.

        ZZ0MAX = 0.0

!        CICE0 = 0.25
!        CICEN = 0.75
!        FLAGTR =4
!        P2SF = 1
!        BSSUBGRID = 0.1
        ZZ0MAX = 1.0020
        SSINTHP = 2.0
        SSWELLFPAR = 3
        SSWELLF(1) = 0.80
        TTAUWSHELTER = 1.0
        SSWELLF(2) = -0.018
        SSWELLF(3) = 0.015
        ZZ0RAT = 0.04
        SSWELLF(4) = 100000
        SSWELLF(5) = 1.2
        SSDSC1 = -4.2
!        SDSLF = 0.
!        SDSHF = 0.
        SSDSBR = 0.00090
        SSDSC2 = -2.2E-5
        SSDSC3 = -0.80
        SSDSC4 = 1.
        SSDSC5 = 0.
        SSDSC6 = 0.30
        SSDSDTH = 80.
!        FXFM3 = 9.9
        SSDSCOS = 2.
        SSDSISO = 2
!        NLPROP = 2.5E7


        SDSNTH  = MIN(NINT(SSDSDTH/(DTH*RADDEG)),NTH/2-1)
        ALLOCATE(IKTAB(MK,2000)); IKTAB = 0
        ALLOCATE(SATINDICES(MTH,2*SDSNTH+1)); SATINDICES = 0
        ALLOCATE(SATWEIGHTS(MTH,2*SDSNTH+1)); SATWEIGHTS = 0.
        ALLOCATE(CUMULW(MK,MTH,MK,MTH)); CUMULW = 0.
        ALLOCATE(DCKI(NKHS,NKD)); DCKI = 0.
        ALLOCATE(LLWS(NSPEC)); LLWS = .FALSE.
        !ALLOCATE(CD(MNP),Z0(MNP),USTDIR(MNP)) ! Already allocated before
        TAUWX = 0.; TAUWY = 0.; CD = 0.; Z0 = 0.; USTDIR = 0.

        IF (LPRECOMPST4) THEN
          CALL INSIN4_OLD(.TRUE.)
        ELSE
          CALL READ_INSIN4_OLD
        END IF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_INSIN4_OLD()
        USE DATAPOOL, ONLY : LPRECOMPST4
        USE W3SRC4MD_OLD
        IF (.NOT. LPRECOMPST4) THEN
          READ (5002)                &
          ZZWND, AALPHA, ZZ0MAX, BBETA, SSINTHP, ZZALP,    &
          TTAUWSHELTER, SSWELLFPAR, SSWELLF,               &
          ZZ0RAT, SSDSC1, SSDSC2, SSDSC3, SSDSC4, SSDSC5,  &
          SSDSC6, SSDSISO, SSDSBR, SSDSBR2, SSDSBM, SSDSP, &
          SSDSCOS, SSDSDTH, WWNMEANP, WWNMEANPTAIL, SSTXFTF, &
          SSTXFTFTAIL, SSTXFTWN, SSTXFTF, SSTXFTWN,        &
          SSDSBRF1, SSDSBRF2, SSDSBRFDF,SSDSBCK, SSDSABK,  &
          SSDSPBK, SSDSBINT, &
          SSDSHCK, DELUST, DELTAIL, DELTAUW, &
          DELU, DELALP, DELAB, TAUT, TAUHFT, TAUHFT2,      &
          SWELLFT, IKTAB, DCKI, SATINDICES, SATWEIGHTS, &
          DIKCUMUL, CUMULW
        END IF
     END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PREPARE_ARDHUIN_NEW()

        USE DATAPOOL
        USE W3SRC4MD_NEW

        IMPLICIT  NONE

        INTEGER :: IK, ISP, ITH, ITH0
        REAL    :: SIGMA, FR1, RTH0

        NK    = MSC
        MK    = NK  ! ? 
        NTH   = MDC
        MTH   = NTH ! ?
        NSPEC = MSC * MDC
        MSPEC = NSPEC

        ALLOCATE(SIG(0:MSC+1), SIG2(NSPEC), DSIP(0:MSC+1), TH(MDC))
        ALLOCATE(ESIN(MSPEC+MTH), ECOS(MSPEC+MTH), EC2(MSPEC+MTH), ES2(MSPEC+MTH),ESC(MSPEC+MTH) )
        ALLOCATE(DSII(MSC), DDEN(MSC), DDEN2(NSPEC))

        DTH   = DDIR
        FR1   = SPSIG(1)/PI2
        TH    = SPDIR

        RTH0 = 0.
        DO ITH=1, NTH
          TH  (ITH) = DTH * ( RTH0 + REAL(ITH-1) )
          ESIN(ITH) = SIN ( TH(ITH) )
          ECOS(ITH) = COS ( TH(ITH) )
          IF ( ABS(ESIN(ITH)) .LT. 1.E-5 ) THEN
            ESIN(ITH) = 0.
            IF ( ECOS(ITH) .GT. 0.5 ) THEN
              ECOS(ITH) =  1.
            ELSE
              ECOS(ITH) = -1.
              END IF
          END IF
          IF ( ABS(ECOS(ITH)) .LT. 1.E-5 ) THEN
            ECOS(ITH) = 0.
            IF ( ESIN(ITH) .GT. 0.5 ) THEN
              ESIN(ITH) =  1.
            ELSE
              ESIN(ITH) = -1.
            END IF
          END IF
          ES2 (ITH) = ESIN(ITH)**2
          EC2 (ITH) = ECOS(ITH)**2
          ESC (ITH) = ESIN(ITH)*ECOS(ITH)
        END DO
!
        DO IK=2, NK+1
          ITH0 = (IK-1)*NTH
          DO ITH=1, NTH
            ESIN(ITH0+ITH) = ESIN(ITH)
            ECOS(ITH0+ITH) = ECOS(ITH)
            ES2 (ITH0+ITH) = ES2 (ITH)
            EC2 (ITH0+ITH) = EC2 (ITH)
            ESC (ITH0+ITH) = ESC (ITH)
          END DO
        END DO

        WRITE(5001,*) 'DTH'
        WRITE(5001,*) DTH
        WRITE(5001,*) 'FR1'
        WRITE(5001,*) FR1
        WRITE(5001,*) 'TH'
        WRITE(5001,*) TH
        WRITE(5001,*) 'ESIN, ECOS, EC2'
        WRITE(5001,*) ESIN, ECOS, EC2

        WNMEANP = 0.5
        WNMEANPTAIL = -0.5
        WRITE(5001,*) 'WNMEANP, WNMEANPTAIL'
        WRITE(5001,*) WNMEANP, WNMEANPTAIL

        XFR = EXP(FRINTF) ! Check with Fabrice ... should be 1.1

        SIGMA   = FR1 * TPI / XFR**2 ! What is going on here ?
        SXFR    = 0.5 * (XFR-1./XFR)

        WRITE(5001,*) 'XFR, SIGMA, SXFR'
        WRITE(5001,*) XFR, SIGMA, SXFR

        DO IK=0, NK+1
         SIGMA    = SIGMA * XFR ! What is going on here ...
         SIG (IK) = SIGMA
         DSIP(IK) = SIGMA * SXFR
        END DO

        WRITE(5001,*) 'SIGMA'
        WRITE(5001,*)  SIGMA
        WRITE(5001,*) 'SIG'
        WRITE(5001,*)  SIG
        WRITE(5001,*) 'DSIP'
        WRITE(5001,*)  DSIP

        DSII(1) = 0.5 * SIG( 1) * (XFR-1.)
        DO IK = 2, NK - 1
          DSII(IK) = DSIP(IK)
        END DO
        DSII(NK) = 0.5 * SIG(NK) * (XFR-1.) / XFR

        DDEN = DTH * DSII(:) * SIG(:)

        DO ISP=1, NSPEC
          IK         = 1 + (ISP-1)/NTH
          SIG2 (ISP) = SIG (IK)
          DDEN2(ISP) = DDEN(IK)
        END DO

        WRITE(5001,*) 'SIG2'
        WRITE(5001,*) SIG2
        WRITE(5001,*) 'DSII'
        WRITE(5001,*) DSII
        WRITE(5001,*) 'DDEN'
        WRITE(5001,*) DDEN
        WRITE(5001,*) 'DDEN2'
        WRITE(5001,*) DDEN2

        FTE = 0.25 * SIG(NK) * DTH * SIG(NK)
        FTF = 0.20           * DTH * SIG(NK)

        FACHF  = 5.
        FACHFE = XFR**(-FACHF)

        STXFTFTAIL  = 1./(FACHF-1.-WNMEANPTAIL*2)
        STXFTF      = 1./(FACHF-1.-WNMEANP*2)
        STXFTWN     = 1./(FACHF-1.-WNMEANP*2) * SIG(NK)**(2)

        WRITE(5001,*) 'FTE, FTF, FACHF, FACHFE'
        WRITE(5001,*) FTE, FTF, FACHF, FACHFE

        SSWELLF(1) = 0.8
        SSWELLF(2) = -0.018
        SSWELLF(3) = 0.015
        SSWELLF(4) = 1.E5
        SSWELLF(5) = 1.2
        SSWELLF(6) = 0.
        SSWELLF(7) = 0.

        WRITE(5001,*) 'SSWELLF'
        WRITE(5001,*) SSWELLF

        AALPHA = 0.0095
        BBETA  = 1.31 ! 1.52 for ECMWF
        ZZALP   = 0.006
        ZZWND   = 10.

        SWELLFPAR = 1

        SSDSBRF1   = 0.5
        SSDSHCK    = 1.
        SSDSBCK    = 0.
        SSDSBINT   = 0.3
        SSDSPBK    = 4.
        SSDSABK    = 1.5

        SSDSBR     = 0.90E-3
        SSDSBRFDF  = 0
        SSDSBRF1   = 0.5
        SSDSBRF2   = 0.
        SSDSBR2    = 0.8

        SSDSP      = 2.
        SSDSPBK    = 4.

        SSDSISO     = 2

        SSDSBM(0)  = 1.
        SSDSBM(1)  = 0.
        SSDSBM(2)  = 0.
        SSDSBM(3)  = 0.
        SSDSBM(4)  = 0.

        ZZ0MAX     = 0. 
        SSINTHP    = 2.0
        SSWELLFPAR = 3

        TTAUWSHELTER = 1.0
        ZZ0RAT       = 0.04

        SSDSC1 = 0. 
        SSDSC2 = -2.2E-5 ! AR: was originally 2.2
        !SSDSC3 = -0.80  ! overwritten by SSDSCUM 
        SSDSC4 = 1.
        SSDSC5 = 0.
        SSDSC6 = 0.30

        WHITECAPWIDTH = 0.8

        SSDSCUM = -0.40344
        SSDSDTH = 80. ! not used ...
        SSDSCOS = 2.
        SSDSISO = 2 !AR: changes to isotropic breaking ... 

        SSDSC(1)   = SSDSC1
        SSDSC(2)   = SSDSC2
        SSDSC(3)   = SSDSCUM
        SSDSC(4)   = SSDSC4
        SSDSC(5)   = SSDSC5
        SSDSC(6)   = SSDSC6
        SSDSC(7)   = WHITECAPWIDTH

        SDSNTH  = MIN(NINT(SSDSDTH/(DTH*RADDEG)),NTH/2-1)
        DELAB   = (ABMAX-ABMIN)/REAL(SIZEFWTABLE)

        ALLOCATE(IKTAB(MK,2000)); IKTAB = 0
        ALLOCATE(SATINDICES(2*SDSNTH+1,MTH)); SATINDICES = 0
        ALLOCATE(SATWEIGHTS(2*SDSNTH+1,MTH)); SATWEIGHTS = 0.
        ALLOCATE(CUMULW(MK*MTH,MK*MTH)); CUMULW = 0.
        ALLOCATE(DCKI(NKHS,NKD)); DCKI = 0.
        ALLOCATE(LLWS(NSPEC)); LLWS = .FALSE.
        ALLOCATE(QBI(NKHS,NKD)); QBI = 0.

        TAUWX = 0.; TAUWY = 0.; CD = 0.; Z0 = 0.; USTDIR = 0.

        IF (LPRECOMPST4) THEN
          CALL INSIN4_NEW(.TRUE.)
        ELSE
          CALL READ_INSIN4_NEW
        END IF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_INSIN4_NEW()
        USE DATAPOOL, ONLY : LPRECOMPST4
        USE W3SRC4MD_NEW
        IF (.NOT. LPRECOMPST4) THEN
          READ (5002)                &
        & ZZWND, AALPHA, ZZ0MAX, BBETA, SSINTHP, ZZALP,    &
        & TTAUWSHELTER, SSWELLFPAR, SSWELLF,               &
        & ZZ0RAT, SSDSC1, SSDSC2, SSDSC3, SSDSC4, SSDSC5,  &
        & SSDSC6, SSDSISO, SSDSBR, SSDSBR2, SSDSBM, SSDSP, &
        & SSDSCOS, SSDSDTH, WNMEANP, WNMEANPTAIL, SSTXFTF, &
        & SSTXFTFTAIL, SSTXFTWN, SSTXFTF, SSTXFTWN,        &
        & SSDSBRF1, SSDSBRF2, SSDSBRFDF,SSDSBCK, SSDSABK,  &
        & SSDSPBK, SSDSBINT, &
        & SSDSHCK, DELUST, DELTAIL, DELTAUW, &
        & DELU, DELALP, DELAB, TAUT, TAUHFT, TAUHFT2,      &
        & SWELLFT, IKTAB, DCKI, SATINDICES, SATWEIGHTS, &
        & DIKCUMUL, CUMULW, QBI
        END IF
     END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_WAVE_BOUNDARY(IFILE, IT)
        USE DATAPOOL
        IMPLICIT NONE
      
        INTEGER, INTENT(INOUT) :: IFILE, IT
 
        REAL*8                 :: DTMP
        INTEGER                :: ITMP

!TODO: Makes sure initial condition work also when no wave boundary is set ...
        IF (LBCWA .OR. LBCSP) THEN
          IF (IBOUNDFORMAT == 3) THEN
            IF (LBCSP) THEN
              WRITE(DBG%FHNDL,*) 'SPECTRAL WW3 NOT READY YET'
              STOP 'SPECTRAL WW3 NOT READY YET'
            END IF
            CALL INIT_NETCDF_WW3
            DTMP = (MAIN%BMJD-BND_TIME_ALL_FILES(1,1)) * DAY2SEC
            !WRITE(*,*) DTMP, MAIN%BMJD, BND_TIME_ALL_FILES(1,1)
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
            !WRITE(*,*) IFILE, IT, SUM(NDT_BND_FILE(1:IFILE-1)), NINT(DTMP/SEBO%DELT), SEBO%DELT
            CALL WAVE_BOUNDARY_CONDITION(IFILE,IT,WBAC)
            IF (LBINTER) WBACOLD = WBAC
          ELSE ! BOUNDFORMAT
            CALL WAVE_BOUNDARY_CONDITION(1,1,WBAC)
            IF (LBINTER) WBACOLD = WBAC
          END IF
        END IF

      END SUBROUTINE INIT_WAVE_BOUNDARY
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SET_WAVE_BOUNDARY
        USE DATAPOOL

         IF (LBCWA .OR. LBCSP) THEN
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
         END IF
      END SUBROUTINE SET_WAVE_BOUNDARY
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
