!    25 Mar 2004    4:34 pm
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_WWMINPUT()
         USE DATAPOOL
#ifdef SELFE
         use elfe_msgp, only : myrank,parallel_abort
         use elfe_glbl, only : errmsg
#endif 
         IMPLICIT NONE

         LOGICAL           :: LFLIVE

         INTEGER           :: I, NI(3)
         INTEGER, PARAMETER :: INUMOUTS = 50

         CHARACTER(LEN=20) :: BEGTC, UNITC, ENDTC, NOUTS(INUMOUTS), NLOUTS(INUMOUTS)
         INTEGER           :: NLPOINTS(INUMOUTS)
         REAL*8            :: XOUTS(INUMOUTS), YOUTS(INUMOUTS), CUTOFF(INUMOUTS)
         REAL*8            :: XLOUTS(2*INUMOUTS), YLOUTS(2*INUMOUTS), DELTC
         REAL              :: DEG

         NAMELIST /PROC/ PROCNAME, DIMMODE, LSTEA, LQSTEA, LSPHE, LNAUTIN, &
     &                   BEGTC, DELTC, UNITC, ENDTC, DMIN, DMINTRIAD

         NAMELIST /COUPL/ LCPL, LROMS, LTIMOR, LSHYFEM, RADFLAG, LETOT, NLVT, LWINDWWM, DTCOUP

#ifdef WWMONLY
         NAMELIST /GRID/ MNP, MNE, LCIRD, LSTAG, MINDIR, MAXDIR, MDC, FRLOW, FRHIGH, &
     &                   MSC, FILEGRID, IGRIDTYPE, LSLOP, SLMAX, LVAR1D, LOPTSIG
#elif SELFE
         NAMELIST /GRID/ LCIRD, LSTAG, MINDIR, MAXDIR, MDC, FRLOW, FRHIGH, MSC,  &
     &                   FILEGRID, IGRIDTYPE, LSLOP, SLMAX, LVAR1D, LOPTSIG
#endif
         NAMELIST /INIT/ LHOTR, IHOTPOS, FILEHOT, LINID, INITSTYLE

         NAMELIST /BOUC/ LBCSE, LBCWA, LBCSP, LINHOM, LBSP1D, LBSP2D, LBINTER, BEGTC, DELTC, UNITC, ENDTC, &
     &                   FILEBOUND, IBOUNDFORMAT, FILEWAVE, LMONO_IN, LINDSPRDEG, LPARMDIR, &
     &                   WBHS, WBTP, WBDM, WBDS, WBSS, WBDSMS, WBGAUSS, WBPKEN, &
     &                   NCDF_HS_NAME, NCDF_DIR_NAME, NCDF_SPR_NAME, NCDF_FP_NAME, NCDF_F02_NAME 

         NAMELIST /WIND/ LSEWD, LSTWD, LCWIN, LWDIR, BEGTC, DELTC, UNITC, ENDTC, LINTERWD, &
     &                   WDIR, WVEL, CWINDX, CWINDY, FILEWIND, WINDFAC, IWINDFORMAT, LWINDFROMWWM

         NAMELIST /CURR/ LSECU, BEGTC, DELTC, UNITC, ENDTC, LINTERCU, &
     &                   LSTCU, LCCUR, CCURTX, CCURTY, FILECUR, LERGINP, CURFAC, ICURRFORMAT

         NAMELIST /WALV/ LSEWL, BEGTC, DELTC, UNITC, ENDTC, LINTERWL, &
     &                   LSTWL, LCWLV, CWATLV, FILEWATL, LERGINP, WALVFAC, IWATLVFORMAT

         NAMELIST /ENGS/ MESNL, MESIN, IFRIC, MESBF, FRICC, MESBR, ICRIT, ALPBJ, BRHD, &
     &                   LMAXETOT, MESDS, MESTR, TRICO, TRIRA, TRIURS, LPRECOMPST4

         NAMELIST /NUMS/ ICOMP, AMETHOD, SMETHOD, DMETHOD, LITERSPLIT, LFILTERTH, & 
     &                   MAXCFLTH, FMETHOD, LFILTERCXY, MAXCFLCXY, LFILTERSIG, MAXCFLSIG, &
     &                   LLIMT, LIMFAK, MELIM, LDIFR, IDIFFR, LADVTEST, &
     &                   LCONV, LCFL, LCHKCONV, NQSITER, QSCONV1, QSCONV2, QSCONV3, QSCONV4, QSCONV5, & 
     &                   EPSH1, EPSH2, EPSH3, EPSH4, EPSH5, RTHETA, LEXPIMP, FREQEXP, LVECTOR,IVECTOR, &
     &                   DTMIN_DYN, NDYNITER, DTMIN_SIN, DTMIN_SNL4, DTMIN_SDS, DTMIN_SNL3, DTMIN_SBR, DTMIN_SBF, &
     &                   NDYNITER_SIN, NDYNITER_SNL4, NDYNITER_SDS, NDYNITER_SBR, NDYNITER_SNL3, NDYNITER_SBF

         NAMELIST /OUTP/ BEGTC, DELTC, UNITC, ENDTC, OUTSTYLE, FILEOUT, LOUTITER, LOUTS, &
     &                   IOUTS, NOUTS, XOUTS, YOUTS, CUTOFF, LSIGMAX, LSP1D, LSP2D, LENERGY, &
     &                   LLOUTS, ILOUTS, NLOUTS, XLOUTS, YLOUTS, LNAUTOUT, LWXFN, LMONO_OUT

         NAMELIST /HOTFILE/ BEGTC, DELTC, UNITC, ENDTC, LHOTF, LCYCLEHOT, FILEHOT

         XOUTS = 0.
         YOUTS = 0.
         XLOUTS = 0.
         YLOUTS = 0.
         NOUTS = ''
         NLOUTS = ''
         CUTOFF = 0.

         READ( INP%FHNDL,  NML = PROC)
#ifdef SELFE
         IF (myrank == 0) THEN
#endif 
         WRITE(CHK%FHNDL,  NML = PROC)
         CALL FLUSH(CHK%FHNDL)
#ifdef SELFE
         ENDIF
#endif

         READ( INP%FHNDL,  NML = COUPL)
#ifdef SELFE
         IF (myrank == 0) THEN
#endif 
         WRITE(CHK%FHNDL,  NML = COUPL)
         CALL FLUSH(CHK%FHNDL)
#ifdef SELFE
         ENDIF
#endif

!
!    *** Estimate various timings ...
!
         MAIN%BEGT   = BEGTC
         MAIN%DELT   = DELTC
         MAIN%DTCOUP = DTCOUP
         MAIN%UNIT   = UNITC
         MAIN%ENDT   = ENDTC
         CALL CT2MJD(MAIN%BEGT, MAIN%BMJD)
         CALL CT2MJD(MAIN%ENDT, MAIN%EMJD)
         CALL CU2SEC(MAIN%UNIT, MAIN%DELT)
         MAIN%TOTL = REAL( (MAIN%EMJD - MAIN%BMJD) * DAY2SEC )
         MAIN%ISTP = NINT( MAIN%TOTL / DBLE(MAIN%DELT) )
         MAIN%TMJD = MAIN%BMJD

#ifdef SELFE
         IF (DIMMODE .NE. 2 .OR. .NOT. LCPL) THEN
           WRITE(errmsg,*)'You are running in less than 1d or LCPL = .F.',DIMMODE,LCPL
           call parallel_abort(errmsg)
         endif
         MAIN%DELT   = DT_WWM   !DELTC
         MAIN%DTCOUP = DT_SELFE !coupling time step
#endif
!
!     *** GRID section
!
         READ (INP%FHNDL,   NML = GRID)
#ifdef SELFE
         IF (myrank == 0) THEN
#endif 
         WRITE(CHK%FHNDL,   NML = GRID)
         CALL FLUSH(CHK%FHNDL)
#ifdef SELFE
         ENDIF
#endif

         GRD%FNAME = FILEGRID

         IF (LCIRD) THEN
            MINDIR = 0.0
            MAXDIR = PI2
         ELSE
           IF (LNAUTIN) THEN
             CALL DEG2NAUT (MINDIR, DEG, LNAUTIN)
             MINDIR = DEG
             MINDIR = MINDIR * DEGRAD
             CALL DEG2NAUT (MAXDIR, DEG, LNAUTIN)
             MAXDIR = DEG
             MAXDIR = MAXDIR * DEGRAD
           ELSE
             MINDIR = MINDIR*DEGRAD
             MAXDIR = MAXDIR*DEGRAD
           END IF
         END IF

         IF (FRLOW > FRHIGH) THEN
#ifdef WWMONLY
            WRITE(*,*) 'error, the FRHIG must be greater than FRLOW'
            STOP 'READ WWM INPUT'
#elif SELFE
            WRITE(STAT%FHNDL,*) 'error, the FRHIG must be greater than FRLOW'
            call parallel_abort('SELFEWWM: error, the FRHIG must be greater than FRLOW ... :(....')
#endif

         END IF
!
!     *** INIT section
!
         READ(INP%FHNDL,  NML = INIT)
#ifdef SELFE
         IF (myrank == 0) THEN
#endif 
         WRITE(CHK%FHNDL, NML = INIT)
         CALL FLUSH(CHK%FHNDL)
#ifdef SELFE
         ENDIF
#endif
         HOTIN%FNAME = FILEHOT 

         IF (LHOTR) THEN
#ifdef SELFE
           call parallel_abort('SELFEWWM: no hotstart implemented yet ... :(....')
#endif 
            INQUIRE( FILE = TRIM(HOTIN%FNAME), EXIST = LFLIVE )
            IF (.NOT. LFLIVE) THEN
#ifdef SELFE
               call parallel_abort('cant open the hotfile')
#elif WWMONLY
               WRITE(*,*) "can't open the hotfile"
               STOP 'READ WWM INPUT' 
#endif 
            ELSE
               WRITE(STAT%FHNDL,'("+TRACE...",A)') 'HOTFILE is used as Initital Condition'
            END IF
         END IF
!
!     *** BOUNDARY CONDITIONS section
!
         READ(INP%FHNDL,  NML = BOUC )
#ifdef SELFE
         IF (myrank == 0) THEN
#endif 
         WRITE(CHK%FHNDL, NML = BOUC )
         CALL FLUSH(CHK%FHNDL)
#ifdef SELFE
         ENDIF
#endif

         BND%FNAME = FILEBOUND
         WAV%FNAME = FILEWAVE

         WRITE(STAT%FHNDL,'("+TRACE...",A10,I5)') 'BOUNDARY FILE FORMAT IS', IBOUNDFORMAT

         SEBO%BEGT = BEGTC
         SEBO%DELT = DELTC
         SEBO%UNIT = UNITC
         SEBO%ENDT = ENDTC

         CALL CT2MJD(SEBO%BEGT, SEBO%BMJD)
         CALL CT2MJD(SEBO%ENDT, SEBO%EMJD)

         CALL CU2SEC(SEBO%UNIT, SEBO%DELT)

         SEBO%TMJD = SEBO%BMJD
!
!     *** WIND section
!
         READ(INP%FHNDL, NML = WIND)
#ifdef WWMONLY
         WRITE(CHK%FHNDL, NML = WIND)
#elif SELFE
         IF (myrank == 0) THEN
           WRITE(CHK%FHNDL, NML = WIND)
           CALL FLUSH(CHK%FHNDL)
         ENDIF
#endif 

         WIN%FNAME = TRIM(FILEWIND)

         SEWI%BEGT = BEGTC
         SEWI%DELT = DELTC
         SEWI%UNIT = UNITC
         SEWI%ENDT = ENDTC

         CALL CT2MJD(SEWI%BEGT, SEWI%BMJD)
         CALL CT2MJD(SEWI%ENDT, SEWI%EMJD)

         CALL CU2SEC(SEWI%UNIT, SEWI%DELT)

         SEWI%TMJD = 0.0
!
!     *** CURR section
!
         READ(INP%FHNDL, NML = CURR)
#ifdef SELFE
         IF (myrank == 0) THEN
#endif 
         WRITE(CHK%FHNDL, NML = CURR)
         CALL FLUSH(CHK%FHNDL)
#ifdef SELFE
         ENDIF
#endif
         CUR%FNAME = TRIM(FILECUR)

         IF (LSECU .AND. LSTCU) THEN
            WRITE(DBG%FHNDL,*) 'Error: LSECU and LSTCU are .TRUE. value'
            STOP 'READ WWM INPUT'
         END IF

         SECU%BEGT = BEGTC
         SECU%DELT = DELTC
         SECU%UNIT = UNITC
         SECU%ENDT = ENDTC

         CALL CT2MJD(SECU%BEGT, SECU%BMJD) ! convert string to internal time ... in double ...
         CALL CT2MJD(SECU%ENDT, SECU%EMJD)

         CALL CU2SEC(SECU%UNIT, SECU%DELT)

         SECU%TMJD = 0.0
!
!     *** water level section
!
         READ(INP%FHNDL, NML = WALV)
#ifdef SELFE
         IF (myrank == 0) THEN
#endif 
         WRITE(CHK%FHNDL, NML = WALV)
         CALL FLUSH(CHK%FHNDL)
#ifdef SELFE
         ENDIF
#endif

         WAT%FNAME = FILEWATL

         IF (LSEWL .AND. LSTWL) THEN
            WRITE(DBG%FHNDL,*) 'Error: LSEWL and LSTWL have .TRUE. value'
            STOP 'READ WWM INPUT'
         END IF

         SEWL%BEGT = BEGTC
         SEWL%DELT = DELTC
         SEWL%UNIT = UNITC
         SEWL%ENDT = ENDTC

         CALL CT2MJD(SEWL%BEGT, SEWL%BMJD)
         CALL CT2MJD(SEWL%ENDT, SEWL%EMJD)

         CALL CU2SEC(SEWL%UNIT, SEWL%DELT)
         SEWL%TMJD = 0.0
!
!     *** ENGS section
!
         READ(INP%FHNDL, NML = ENGS)
#ifdef SELFE
         IF (myrank == 0) THEN
#endif 
         WRITE(CHK%FHNDL, NML = ENGS)
         CALL FLUSH(CHK%FHNDL)
#ifdef SELFE
         ENDIF
#endif
!
!     *** NUMS section
!
         READ(INP%FHNDL, NML = NUMS)
#ifdef SELFE
         IF (myrank == 0) THEN
#endif 
         WRITE(CHK%FHNDL, NML = NUMS)
         CALL FLUSH(CHK%FHNDL)
#ifdef SELFE
         ENDIF
#endif
!
!     **** OUTP section
!
         READ(INP%FHNDL, NML = OUTP)
#ifdef SELFE
         IF (myrank == 0) THEN
#endif 
         WRITE(CHK%FHNDL, NML = OUTP)
         CALL FLUSH(CHK%FHNDL)
#ifdef SELFE
         ENDIF
#endif

         OUT%FNAME = FILEOUT

         OUTF%BEGT = BEGTC
         OUTF%DELT = DELTC
         OUTF%UNIT = UNITC
         OUTF%ENDT = ENDTC

         CALL CT2MJD(OUTF%BEGT, OUTF%BMJD)
         CALL CT2MJD(OUTF%ENDT, OUTF%EMJD)
         CALL CU2SEC(OUTF%UNIT, OUTF%DELT)

         OUTF%TOTL = REAL( (OUTF%EMJD - OUTF%BMJD) * DAY2SEC )
         OUTF%ISTP = NINT( OUTF%TOTL / DBLE(OUTF%DELT) ) + 1
         OUTF%TMJD = OUTF%BMJD

         IF (LOUTS) THEN

           ALLOCATE ( STATION(IOUTS) )
#ifdef SELFE
           STATION%IFOUND = 0
           STATION%ISUM   = 0 
#endif
           DO I = 1, IOUTS
             STATION(I)%XELE = 0.
             STATION(I)%YELE = 0.
             STATION(I)%ZELE = 0.
             STATION(I)%WI   = 0.d0
             STATION(I)%OUTPAR_NODE = 0.
           END DO
           STATION(1:IOUTS)%NAME   = NOUTS(1:IOUTS)
           STATION(1:IOUTS)%XCOORD = XOUTS(1:IOUTS)
           STATION(1:IOUTS)%YCOORD = YOUTS(1:IOUTS)

           IF (LSIGMAX) THEN
             STATION(1:IOUTS)%CUTOFF = CUTOFF(1:IOUTS)
           ELSE
             STATION(1:IOUTS)%CUTOFF = FRHIGH
           END IF

           WRITE(DBG%FHNDL,*) 'STATION X and Y Coordinates'  
           WRITE(DBG%FHNDL,*) STATION%XCOORD
           WRITE(DBG%FHNDL,*) STATION%YCOORD
           WRITE(DBG%FHNDL,*) 'STATION Names'
           WRITE(DBG%FHNDL,*) STATION%NAME
           CALL FLUSH(DBG%FHNDL)

        END IF

         IF (LLOUTS) THEN

          ALLOCATE ( LINES(IOUTS) )
#ifdef SELFE
           STATION%IFOUND = 0
           STATION%ISUM   = 0
#endif
           DO I = 1, IOUTS
             STATION(I)%XELE = 0.
             STATION(I)%YELE = 0.
             STATION(I)%ZELE = 0.
             STATION(I)%WI   = 0.d0
             STATION(I)%OUTPAR_NODE = 0.
           END DO
           STATION(1:IOUTS)%NAME   = NOUTS(1:IOUTS)
           STATION(1:IOUTS)%XCOORD = XOUTS(1:IOUTS)
           STATION(1:IOUTS)%YCOORD = YOUTS(1:IOUTS)

           IF (LSIGMAX) THEN
             STATION(1:IOUTS)%CUTOFF = CUTOFF(1:IOUTS)
           ELSE
             STATION(1:IOUTS)%CUTOFF = FRHIGH
           END IF

           WRITE(DBG%FHNDL,*) 'STATION X and Y Coordinates'
           WRITE(DBG%FHNDL,*) STATION%XCOORD
           WRITE(DBG%FHNDL,*) STATION%YCOORD
           WRITE(DBG%FHNDL,*) 'STATION Names'
           WRITE(DBG%FHNDL,*) STATION%NAME
           CALL FLUSH(DBG%FHNDL)

        END IF
!
! set the output flag
!
         IF (     TRIM(OUTSTYLE) == 'NO') THEN
            IOUTP = 0
         ELSE IF (TRIM(OUTSTYLE) == 'XFN') THEN
            IOUTP = 1
         ELSE IF (TRIM(OUTSTYLE) == 'NC') THEN
            IOUTP = 2
         ELSE IF (TRIM(OUTSTYLE) == 'GMT') THEN
            IOUTP = 3 
         ELSE IF (TRIM(OUTSTYLE) == 'SHP') THEN
            IOUTP = 4 
         END IF
!
!     **** HOTFILE section
!
         READ(INP%FHNDL, NML = HOTFILE)
#ifdef SELFE
         IF (myrank == 0) THEN
#endif 
         WRITE(CHK%FHNDL, NML = HOTFILE)
         CALL FLUSH(CHK%FHNDL)
#ifdef SELFE
         ENDIF
#endif
         HOTOUT%FNAME = FILEHOT

         HOTF%BEGT = BEGTC
         HOTF%DELT = DELTC
         HOTF%UNIT = UNITC
         HOTF%ENDT = ENDTC

         CALL CT2MJD(HOTF%BEGT, HOTF%BMJD)
         CALL CT2MJD(HOTF%ENDT, HOTF%EMJD)

         CALL CU2SEC(HOTF%UNIT, HOTF%DELT)

         HOTF%TOTL = REAL( (HOTF%EMJD - HOTF%BMJD) * DAY2SEC )
         HOTF%ISTP = NINT( HOTF%TOTL / DBLE(HOTF%DELT) ) + 1
         HOTF%TMJD = HOTF%BMJD

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_WIND_CURRENT_WATLEV()
#ifdef NCDF
         USE NETCDF
#endif
         USE DATAPOOL
#ifdef SELFE
         USE ELFE_MSGP, ONLY : myrank, comm, ierr
         USE ELFE_GLBL, ONLY : ipgl
#endif
         IMPLICIT NONE

         INTEGER :: IP, ISTAT, IT, IFILE, FORECASTHOURS
         LOGICAL :: LFLIVE
         REAL    :: WDIRT
         REAL*8  :: DTMP, WI(3)
         CHARACTER(LEN=20) :: TIMESTRING

         FORECASTHOURS = 0
  
!
!     *** wind field ...
!
         IF (LWINDFROMWWM) THEN

           WINDXY(:,:) = 0.0

           IF (LSTWD) THEN
             IF (LCWIN) THEN
               WRITE(STAT%FHNDL,'("+TRACE...",A)') 'HOMOGENOUS STEADY WIND FIELD IS USED'
               WRITE(STAT%FHNDL,'("+TRACE...",A,I10)') 'WIND IS COMING FROM WWM - WINDFORMAT', IWINDFORMAT
               IF (LWDIR) THEN
                 CALL DEG2NAUT(WDIR, WDIRT, LNAUTIN)
                 DO IP = 1, MNP
                   WINDXY(IP,1) =  WVEL * COS(WDIRT * DEGRAD)
                   WINDXY(IP,2) =  WVEL * SIN(WDIRT * DEGRAD)
                 END DO
               ELSE
                 DO IP = 1, MNP
                   WINDXY(IP,1) = CWINDX
                   WINDXY(IP,2) = CWINDY
                 END DO
               END IF
             ELSE ! LCWIN
               IF (IWINDFORMAT == 1) THEN
                 WRITE(STAT%FHNDL,'("+TRACE...",A)') 'SPATIAL VARIABLE WIND FIELD IS USED        '
                 WRITE(STAT%FHNDL,'("+TRACE...",A,I10)') 'WIND IS COMING FROM WWM - WINDFORMAT', IWINDFORMAT
                 INQUIRE( FILE = TRIM(WIN%FNAME), EXIST = LFLIVE )
                 IF ( .NOT. LFLIVE ) THEN
                   CALL ERROR_MSG("can't open the wind velocity file", 'PREPARE')
                 ELSE
                   OPEN(WIN%FHNDL, FILE = TRIM(WIN%FNAME), STATUS = 'OLD')
                   READ(WIN%FHNDL,*)
                   READ(WIN%FHNDL, *, IOSTAT = ISTAT) WINDXY(:,1)
                   IF ( ISTAT > 0 ) CALL ERROR_MSG('error in the wind velocity file', 'PREPARE')
                   READ(WIN%FHNDL, *, IOSTAT = ISTAT) WINDXY(:,2)
                   IF ( ISTAT > 0 ) CALL ERROR_MSG('error in the wind velocity file', 'PREPARE')
                   CLOSE(WIN%FHNDL)
                 END IF
#ifdef NCDF
               ELSE IF (IWINDFORMAT == 2) THEN ! NETCDF created using ncl_convert2nc using DWD grib 
                 WRITE(*,'("+TRACE...",A)') 'SPATIAL VARIABLE WIND FIELD IS USED DWD_NETCDF'
                 WRITE(STAT%FHNDL,'("+TRACE...",A,I10)') 'WIND IS COMING FROM WWM - WINDFORMAT', IWINDFORMAT
                 INQUIRE( FILE = TRIM(WIN%FNAME), EXIST = LFLIVE )
                 IF ( .NOT. LFLIVE ) THEN
                 CALL ERROR_MSG("can't open the wind velocity file", 'PREPARE')
                 ELSE
                   CALL INIT_NETCDF_DWD
                   DTMP = SEWI%BMJD-WIND_TIME_ALL_FILES(1)
                   IF (DTMP .GT. -THR8) THEN
                     IFILE = NINT(DTMP/SEWI%DELT/DBLE(NDT_WIND_FILE)) + 1
                     IT    = INT(DTMP/SEWI%DELT)-(IFILE-1)*NDT_WIND_FILE + 1
                     ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(IFILE), NF90_NOWRITE, WIND_NCID)
                     CALL READ_NETCDF_DWD(IT)
                     ISTAT = NF90_CLOSE(WIND_NCID)
                     CALL INTER_STRUCT_DATA(NDX_WIND,NDY_WIND,DX_WIND,DY_WIND,OFFSET_X_WIND,OFFSET_Y_WIND,WIND_X,WINDXY(:,1))
                     CALL INTER_STRUCT_DATA(NDX_WIND,NDY_WIND,DX_WIND,DY_WIND,OFFSET_X_WIND,OFFSET_Y_WIND,WIND_Y,WINDXY(:,2))
                   ELSE
                     WRITE(STAT%FHNDL,*) 'DATA OF STEADY WIND FIELD OUT OF TIME HISTORY OF AVAILABLE WIND DATA' 
                   END IF
                 END IF
               ELSE IF (IWINDFORMAT == 3) THEN ! NETCDF created using cdo -f nc copy file.grb file.nc 
                 WRITE(STAT%FHNDL,'("+TRACE...",A)') 'SPATIAL VARIABLE WIND FIELD IS USED NARR_NETCDF'
                 WRITE(STAT%FHNDL,'("+TRACE...",A,I10)') 'WIND IS COMING FROM WWM - WINDFORMAT', IWINDFORMAT
                 INQUIRE( FILE = TRIM(WIN%FNAME), EXIST = LFLIVE )
                 IF ( .NOT. LFLIVE ) THEN
                   CALL ERROR_MSG("can't open the wind velocity file", 'PREPARE')
                 ELSE
                   CALL INIT_NETCDF_NARR
                   DTMP = SEWI%BMJD-WIND_TIME_ALL_FILES(1)
                   IF (DTMP .GT. -THR8) THEN
                     IFILE = NINT(DTMP/SEWI%DELT/DBLE(NDT_WIND_FILE)) + 1
                     IT    = INT(DTMP/SEWI%DELT)-(IFILE-1)*NDT_WIND_FILE + 1
                     ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(IFILE), NF90_NOWRITE, WIND_NCID)
                     WRITE(STAT%FHNDL,*) 'CALLING READ_NETCDF_NARR FROM INIT_WIND_CURRENT_WATLEV', IT
                     CALL READ_NETCDF_NARR(IT)
                     ISTAT = NF90_CLOSE(WIND_NCID)
                     CALL INTER_STRUCT_DATA(NDX_WIND,NDY_WIND,DX_WIND,DY_WIND,OFFSET_X_WIND,OFFSET_Y_WIND,WIND_X,WINDXY(:,1))
                     CALL INTER_STRUCT_DATA(NDX_WIND,NDY_WIND,DX_WIND,DY_WIND,OFFSET_X_WIND,OFFSET_Y_WIND,WIND_Y,WINDXY(:,2))
                   ELSE
                     WRITE(STAT%FHNDL,*) 'DATA OF STEADY WIND FIELD OUT OF TIME HISTORY OF AVAILABLE WIND DATA'
                   END IF
                 END IF
               ELSE IF (IWINDFORMAT == 4) THEN ! NETCDF NARR downloaded from NOMAD 
                 WRITE(STAT%FHNDL,'("+TRACE...",A)') 'SPATIAL VARIABLE WIND FIELD IS USED NARR_NETCDF'
                 WRITE(STAT%FHNDL,'("+TRACE...",A,I10)') 'WIND IS COMING FROM WWM - WINDFORMAT', IWINDFORMAT
                 INQUIRE( FILE = TRIM(WIN%FNAME), EXIST = LFLIVE )
                 IF ( .NOT. LFLIVE ) THEN
                   CALL ERROR_MSG("can't open the wind velocity file", 'PREPARE')
                 ELSE
                   CALL INIT_NETCDF_NARR2
                   DTMP = SEWI%BMJD-WIND_TIME_ALL_FILES(1)
                   IF (DTMP .GT. -THR8) THEN
                     IFILE = NINT(DTMP/SEWI%DELT/DBLE(NDT_WIND_FILE)) + 1
                     IT    = INT(DTMP/SEWI%DELT)-(IFILE-1)*NDT_WIND_FILE + 1
                     ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(IFILE), NF90_NOWRITE, WIND_NCID)
                     WRITE(STAT%FHNDL,*) 'CALLING READ_NETCDF_NARR FROM INIT_WIND_CURRENT_WATLEV', IT
                     CALL READ_NETCDF_NARR2(IT)
                     ISTAT = NF90_CLOSE(WIND_NCID)
                     CALL INTER_STRUCT_DATA(NDX_WIND,NDY_WIND,DX_WIND,DY_WIND,OFFSET_X_WIND,OFFSET_Y_WIND,WIND_X,WINDXY(:,1))
                     CALL INTER_STRUCT_DATA(NDX_WIND,NDY_WIND,DX_WIND,DY_WIND,OFFSET_X_WIND,OFFSET_Y_WIND,WIND_Y,WINDXY(:,2))
                   ELSE
                     WRITE(STAT%FHNDL,*) 'DATA OF STEADY WIND FIELD OUT OF TIME HISTORY OF AVAILABLE WIND DATA'
                   END IF
                 END IF
#endif
               END IF
             ENDIF
           ELSE IF (LSEWD) THEN

             IF (LCWIN) THEN
               STOP 'LSEWD + LCWIN NOT READY'
             ELSE
               IF (IWINDFORMAT == 1) THEN
                 WRITE(*,'("+TRACE...",A)') 'NONSTATIONARY WIND FIELD IS USED        '
                 WRITE(STAT%FHNDL,'("+TRACE...",A)') 'WIND IS COMING FROM WWM - WINDFORMAT', IWINDFORMAT
                 SEWI%TOTL = REAL( (SEWI%EMJD - SEWI%BMJD) * DAY2SEC )
                 SEWI%ISTP = NINT( SEWI%TOTL / DBLE(SEWI%DELT) ) + 1
                 SEWI%TMJD = SEWI%BMJD
                 WRITE(STAT%FHNDL,*) 'Serial Wind Condition -----------'
                 WRITE(STAT%FHNDL,*) SEWI%BEGT, SEWI%ENDT, SEWI%ISTP, SEWI%TOTL/3600.0, SEWI%DELT
                 INQUIRE( FILE = TRIM(WIN%FNAME), EXIST = LFLIVE )
                 IF ( .NOT. LFLIVE ) THEN
                   CALL ERROR_MSG("can't open the wind velocity file", 'PREPARE')
                 ELSE
                   OPEN(WIN%FHNDL, FILE = TRIM(WIN%FNAME), STATUS = 'OLD')
                   CALL CSEVAL( WIN%FHNDL, TRIM(WIN%FNAME), .FALSE., LWINFILE, MNP, 2, WINDXY, SEWI%DELT, MAIN%DELT, DVWIND )
                 ENDIF
#ifdef NCDF
               ELSE IF (IWINDFORMAT == 2) THEN ! NETCDF created using ncl_convert2nc using DWD grib 
                 WRITE(STAT%FHNDL,'("+TRACE...",A)') 'NONSTATIONARY WIND FIELD IS USED        '
                 WRITE(STAT%FHNDL,'("+TRACE...",A)') 'WIND IS COMING FROM WWM - WINDFORMAT', IWINDFORMAT
                 SEWI%TOTL = REAL( (SEWI%EMJD - SEWI%BMJD) * DAY2SEC )
                 SEWI%ISTP = NINT(  SEWI%TOTL / DBLE(SEWI%DELT) ) + 1
                 SEWI%TMJD = SEWI%BMJD
                 WRITE(STAT%FHNDL,'("+TRACE...",A)') 'SPATIAL VARIABLE WIND FIELD IS USED DWD_NETCDF'
                 INQUIRE( FILE = TRIM(WIN%FNAME), EXIST = LFLIVE )
                 IF ( .NOT. LFLIVE ) THEN
                   CALL ERROR_MSG("can't open the wind velocity file", 'PREPARE')
                 ELSE
                   CALL INIT_NETCDF_DWD
                   DTMP = SEWI%BMJD-WIND_TIME_ALL_FILES(1)
                   IF (DTMP .GT. -THR8) THEN
!2do check for the dwd data setup initial condition agian !
                     IFILE = NINT(DTMP/SEWI%DELT/DBLE(NDT_WIND_FILE)) + 1
                     IT    = NINT(DTMP/SEWI%DELT)-(IFILE-1)*NDT_WIND_FILE + 1
                     ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(IFILE), NF90_NOWRITE, WIND_NCID)
                     CALL READ_NETCDF_DWD(IT)
                     ISTAT = NF90_CLOSE(WIND_NCID)
                     CALL INTER_STRUCT_DATA(NDX_WIND,NDY_WIND,DX_WIND,DY_WIND,OFFSET_X_WIND,OFFSET_Y_WIND,WIND_X,WINDXY(:,1))
                     CALL INTER_STRUCT_DATA(NDX_WIND,NDY_WIND,DX_WIND,DY_WIND,OFFSET_X_WIND,OFFSET_Y_WIND,WIND_Y,WINDXY(:,2))
                     CALL INTER_STRUCT_DATA(NDX_WIND,NDY_WIND,DX_WIND,DY_WIND,OFFSET_X_WIND,OFFSET_Y_WIND,ATMO_PRESS,PRESSURE)
                   ELSE
                     WRITE(DBG%FHNDL,*) 'DATA OF STEADY WIND FIELD OUT OF TIME HISTORY OF AVAILABLE WIND DATA'
                   END IF
                 END IF
               ELSE IF (IWINDFORMAT == 3) THEN ! NETCDF created using cdo -f nc copy file.grb file.nc 
 
                 WRITE(STAT%FHNDL,'("+TRACE...",A)') 'NONSTATIONARY WIND FIELD IS USED        '
                 WRITE(STAT%FHNDL,'("+TRACE...",A)') 'WIND IS COMING FROM WWM - WINDFORMAT', IWINDFORMAT
                 SEWI%TOTL = REAL( (SEWI%EMJD - SEWI%BMJD) * DAY2SEC )
                 SEWI%ISTP = NINT(  SEWI%TOTL / DBLE(SEWI%DELT) ) + 1
                 SEWI%TMJD = SEWI%BMJD
                 WRITE(STAT%FHNDL,'("+TRACE...",A)') 'SPATIAL VARIABLE WIND FIELD IS USED NARR_NETCDF'
                 INQUIRE( FILE = TRIM(WIN%FNAME), EXIST = LFLIVE )

                 IF ( .NOT. LFLIVE ) THEN
                   CALL ERROR_MSG("WIND VEL. FILE IS NOT PRESENT", 'PREPARE')
                 ELSE
                   CALL INIT_NETCDF_NARR
                   DTMP = (MAIN%TMJD-WIND_TIME_ALL_FILES(1)) * DAY2SEC
                   IF (DTMP .GT. -THR8) THEN
                     IFILE = INT(DTMP/SEWI%DELT/DBLE(NDT_WIND_FILE)) + 1
                     !IT    = NINT(DTMP/SEWI%DELT)-(IFILE-1)*NDT_WIND_FILE  
                     IT    = NINT(DTMP/SEWI%DELT)-(IFILE-1)*NDT_WIND_FILE + NINT((MAIN%TMJD-WIND_TIME_ALL_FILES(1)) * 4)
                     !IT = IT * 7 + FORECASTHOURS + 1
                     IF (IT .GT. NDT_WIND_FILE) THEN
                       IFILE = IFILE + 1
                       IT    = 1
                     ENDIF
                     IF (WIND_TIME_ALL_FILES(IT+(IFILE-1)*NDT_WIND_FILE) .EQ. WIND_TIME_ALL_FILES(IT+1+(IFILE-1)*NDT_WIND_FILE)) THEN
                       IT = IT + 1 ! Skip the initital condition 
                     END IF
                     WRITE(STAT%FHNDL,'(A36,2I10,4F15.4)') 'CALLING READ_NETCDF_NARR FROM INIT_WIND_CURRENT_WATLEV', IT, FORECASTHOURS, DTMP, SEWI%DELT, MAIN%TMJD, WIND_TIME_ALL_FILES(1)
                     ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(IFILE), NF90_NOWRITE, WIND_NCID)
                     CALL READ_NETCDF_NARR(IT)
                     ISTAT = NF90_CLOSE(WIND_NCID)
                     CALL INTER_STRUCT_DATA(NDX_WIND,NDY_WIND,DX_WIND,DY_WIND,OFFSET_X_WIND,OFFSET_Y_WIND,WIND_X,WINDXY(:,1))
                     CALL INTER_STRUCT_DATA(NDX_WIND,NDY_WIND,DX_WIND,DY_WIND,OFFSET_X_WIND,OFFSET_Y_WIND,WIND_Y,WINDXY(:,2))
                    !CALL INTER_STRUCT_DATA(NDX_WIND,NDY_WIND,DX_WIND,DY_WIND,OFFSET_X_WIND,OFFSET_Y_WIND,ATMO_PRESS,PRESSURE)
                   ELSE
                     WRITE(DBG%FHNDL,*) 'DATA OF WIND FIELD OUT OF TIME HISTORY OF AVAILABLE WIND DATA'
                   END IF
                 END IF

               ELSE IF (IWINDFORMAT == 4) THEN ! NETCDF created using cdo -f nc copy file.grb file.nc 

                 WRITE(STAT%FHNDL,'("+TRACE...",A)') 'NONSTATIONARY WIND FIELD IS USED        '
                 WRITE(STAT%FHNDL,'("+TRACE...",A)') 'WIND IS COMING FROM WWM - WINDFORMAT', IWINDFORMAT
                 SEWI%TOTL = REAL( (SEWI%EMJD - SEWI%BMJD) * DAY2SEC )
                 SEWI%ISTP = NINT(  SEWI%TOTL / DBLE(SEWI%DELT) ) + 1
                 SEWI%TMJD = SEWI%BMJD
                 WRITE(STAT%FHNDL,'("+TRACE...",A)') 'SPATIAL VARIABLE WIND FIELD IS USED NARR_NETCDF'
                 INQUIRE( FILE = TRIM(WIN%FNAME), EXIST = LFLIVE )

                 IF ( .NOT. LFLIVE ) THEN
                   CALL ERROR_MSG("WIND VEL. FILE IS NOT PRESENT", 'PREPARE')
                 ELSE
                   CALL INIT_NETCDF_NARR2
                   DTMP = (MAIN%TMJD-WIND_TIME_ALL_FILES(1)) * DAY2SEC
                   IF (DTMP .GT. -THR8) THEN
                     IFILE = INT(DTMP/SEWI%DELT/DBLE(NDT_WIND_FILE)) + 1
                     IT    = NINT(DTMP/SEWI%DELT)-(IFILE-1)*NDT_WIND_FILE
                     !WRITE(*,*) MAIN%TMJD, WIND_TIME_ALL_FILES(1), SEWI%DELT
                     IF (IT .GT. NDT_WIND_FILE) THEN
                       IFILE = IFILE + 1
                       IT    = 1
                     ENDIF
                     WRITE(STAT%FHNDL,'(A36,2I10,4F15.4)') 'CALLING READ_NETCDF_NARR2 FROM INPUT', IT, DTMP, SEWI%DELT, MAIN%TMJD, WIND_TIME_ALL_FILES(1)
                     ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(IFILE), NF90_NOWRITE, WINDX_NCID)
                     ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(IFILE), NF90_NOWRITE, WINDY_NCID)
                     CALL READ_NETCDF_NARR2(IT)
                     ISTAT = NF90_CLOSE(WIND_NCID)
                     DO IP = 1, MNP
                       !WRITE(*,*) WIND_ELE(IP), INE_WIND(:,WIND_ELE(IP))
                       CALL INTELEMENT(XYPWIND(1,INE_WIND(:,WIND_ELE(IP))),XYPWIND(2,INE_WIND(:,WIND_ELE(IP))),UWND_NARR(INE_WIND(:,WIND_ELE(IP))),XP(IP),YP(IP),Wi,WINDXY(IP,1),.FALSE.) 
                       CALL INTELEMENT(XYPWIND(1,INE_WIND(:,WIND_ELE(IP))),XYPWIND(2,INE_WIND(:,WIND_ELE(IP))),VWND_NARR(INE_WIND(:,WIND_ELE(IP))),XP(IP),YP(IP),Wi,WINDXY(IP,2),.FALSE.)
                     END DO
                   ELSE
                     WRITE(DBG%FHNDL,*) 'DATA OF WIND FIELD OUT OF TIME HISTORY OF AVAILABLE WIND DATA'
                   END IF
                 END IF
#endif
              ENDIF 
            ENDIF
          ENDIF
        END IF ! LWINDFROMWWM
!
!     *** currents ...
!
       IF (.NOT. LCPL) THEN
         CURTXY(:,:) = 0.0
         IF (LSTCU) THEN
            IF (DIMMODE .EQ. 1) THEN
              IF (LCCUR) THEN
                DO IP = 1, MNP
                  CURTXY(IP,1) = CCURTX
                END DO
              ELSE
                INQUIRE( FILE = TRIM(CUR%FNAME), EXIST = LFLIVE )
                IF ( .NOT. LFLIVE ) THEN
                   CALL ERROR_MSG("can't open the current velocity file", 'PREPARE')
                ELSE
                  OPEN(CUR%FHNDL, FILE = TRIM(CUR%FNAME), STATUS = 'OLD')
                  READ(CUR%FHNDL, *, IOSTAT = ISTAT) CURTXY(:,1)
                  IF ( ISTAT > 0 ) CALL ERROR_MSG('error in the current velocity file', 'PREPARE')
                  CLOSE(CUR%FHNDL)
                END IF
              END IF
            ELSE IF (DIMMODE .EQ. 2) THEN
              IF (LCCUR) THEN
                DO IP = 1, MNP
                  CURTXY(IP,1) = CCURTX
                  CURTXY(IP,2) = CCURTY
                END DO
              ELSE
                INQUIRE( FILE = TRIM(CUR%FNAME), EXIST = LFLIVE )
                IF ( .NOT. LFLIVE ) THEN
                   CALL ERROR_MSG("can't open the current velocity file", 'PREPARE')
                ELSE
                  OPEN(CUR%FHNDL, FILE = TRIM(CUR%FNAME), STATUS = 'OLD')
                  READ(CUR%FHNDL, *, IOSTAT = ISTAT) CURTXY(:,1)
                  IF ( ISTAT > 0 ) CALL ERROR_MSG('error in the current velocity file', 'PREPARE')
                  READ(CUR%FHNDL, *, IOSTAT = ISTAT) CURTXY(:,2)
                  IF ( ISTAT > 0 ) CALL ERROR_MSG('error in the current velocity file', 'PREPARE')
                  CLOSE(CUR%FHNDL)
                END IF
              END IF
            END IF
         ELSE IF (LSECU) THEN
            SECU%TOTL = REAL( (SECU%EMJD - SECU%BMJD) * DAY2SEC )
            SECU%ISTP = NINT( SECU%TOTL / DBLE(SECU%DELT) ) + 1
            SECU%TMJD = SECU%BMJD
            LSECN = .FALSE.
            WRITE(STAT%FHNDL,*) 'Serial current Condition -----------'
            WRITE(STAT%FHNDL,*) SECU%BEGT, SECU%ENDT, SECU%ISTP, SECU%TOTL/3600.0, SECU%DELT
            IF (LERGINP) CALL ERG2WWM(SECU%ISTP)
            INQUIRE( FILE = TRIM(CUR%FNAME), EXIST = LFLIVE )
            IF ( .NOT. LFLIVE ) THEN
               CALL ERROR_MSG("can't open the current velocity file", 'PREPARE')
            ELSE
              LSECN = .TRUE.
              OPEN(CUR%FHNDL, FILE = TRIM(CUR%FNAME), STATUS = 'OLD')
              CALL CSEVAL( CUR%FHNDL, TRIM(CUR%FNAME), .FALSE., LCURFILE, MNP, 2, CURTXY, SECU%DELT, MAIN%DELT, DVCURT )
            END IF
         END IF
!
!     *** water levels 
!
         WATLEV    = 0.
         WATLEVOLD = 0.

         IF (LSTWL) THEN

            IF (DIMMODE .EQ. 1) THEN
              IF (LCWLV) THEN
                 WATLEV = CWATLV
                 DEP    = WLDEP + WATLEV
              ELSE
                INQUIRE( FILE = TRIM(WAT%FNAME), EXIST = LFLIVE )
                IF ( .NOT. LFLIVE ) THEN
                  CALL ERROR_MSG("can't open the water level file", 'PREPARE')
                ELSE
                  OPEN(WAT%FHNDL, FILE = TRIM(WAT%FNAME), STATUS = 'OLD')
                  READ(WAT%FHNDL, *, IOSTAT = ISTAT) WATLEV(:)
                  IF ( ISTAT > 0 )  CALL ERROR_MSG('error in the water level file', 'PREPARE')
                  CLOSE(WAT%FHNDL)
                END IF
              END IF
            ELSE IF (DIMMODE .EQ. 2) THEN
              IF (LCWLV) THEN
                 WATLEV = CWATLV
                 DEP    = WLDEP + WATLEV
              ELSE
                INQUIRE( FILE = TRIM(WAT%FNAME), EXIST = LFLIVE )
                IF ( .NOT. LFLIVE ) THEN
                  CALL ERROR_MSG("can't open the water level file", 'PREPARE')
                ELSE
                  OPEN(WAT%FHNDL, FILE = TRIM(WAT%FNAME), STATUS = 'OLD')
                  READ(WAT%FHNDL, *, IOSTAT = ISTAT) WATLEV(:)
                  IF ( ISTAT > 0 )  CALL ERROR_MSG('error in the water level file', 'PREPARE')
                  CLOSE(WAT%FHNDL)
                END IF
              END IF
           END IF
         ELSE IF (LSEWL) THEN
            SEWL%TOTL = REAL( (SEWL%EMJD - SEWL%BMJD) * DAY2SEC )
            SEWL%ISTP = NINT( SEWL%TOTL / DBLE(SEWL%DELT) ) + 1
            SEWL%TMJD = SEWL%BMJD
            IF (LERGINP .AND. .NOT. LSECU) CALL ERG2WWM(SEWL%ISTP)
            LSELN = .FALSE.
            WRITE(STAT%FHNDL,*) 'Serial water level Condition -----------'
            WRITE(STAT%FHNDL,*) SEWL%BEGT, SEWL%ENDT, SEWL%ISTP, SEWL%TOTL/3600.0, SEWL%DELT
            INQUIRE( FILE = TRIM(WAT%FNAME), EXIST = LFLIVE )
            IF ( .NOT. LFLIVE ) THEN
               CALL ERROR_MSG("can't open the water level file", 'PREPARE')
            ELSE
               LSELN = .TRUE.
               OPEN(WAT%FHNDL, FILE = TRIM(WAT%FNAME), STATUS = 'OLD')
               CALL CSEVAL( WAT%FHNDL,TRIM(WAT%FNAME), .FALSE., LWATLFILE, MNP, 1, WATLEV, SEWL%DELT, MAIN%DELT, DVWALV )
            END IF
         END IF
       
       END IF ! LCPL

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_SPATIAL_GRID()

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER :: I, IP, IE, ISTAT, ITMP
         REAL    :: TMP, BNDTMP 
         REAL*8  :: XPDTMP, YPDTMP, ZPDTMP
         CHARACTER (LEN = 10) :: STRNGTMP
         LOGICAL :: LFLIVE
!begin modification by MDS
         REAL*8 DXP1, DXP2, DXP3, DYP1, DYP2, DYP3, DBLTMP
         INTEGER I1, I2, I3
!end modification by MDS
!
!     *** set the mesh
!
!2DO makes subs for reading the mesh file with different formats ... 
!
         INQUIRE( FILE = TRIM(GRD%FNAME), EXIST = LFLIVE )
         IF ( .NOT. LFLIVE ) THEN
            CALL ERROR_MSG("can't open the grid configuration file",'PREPARE')
         ELSE
            SELECT CASE (DIMMODE)
               CASE (1)
                  OPEN(GRD%FHNDL, FILE = GRD%FNAME, STATUS = 'OLD')
                  DO IP = 1, MNP
                     READ(GRD%FHNDL, *, IOSTAT = ISTAT) XP(IP), DEP(IP)
                     IF ( ISTAT /= 0 ) CALL ERROR_MSG('error in the grid configuration file', 'PREPARE')
                  END DO
                  IF (LVAR1D) THEN
                    DX1(0)     = XP(2)- XP(1)
                    DX1(1)     = DX1(0)
                    DX1(MNP)   = XP(MNP) - XP(MNP-1)
                    DX1(MNP+1) = DX1(MNP)
                    DX2(0)     = DX1(0)
                    DX2(MNP+1) = DX1(MNP)
                    DO IP = 2, MNP-1 ! Bandwith at gridpoints
                       DX1(IP) = (XP(IP)-XP(IP-1))/2. + (XP(IP+1)-XP(IP))/2.
                    END DO
                    DO IP = 2, MNP ! Stepwidth between gridpoints K and K-1
                       DX2(IP) = XP(IP) - XP(IP-1)
                    END DO
                    DX2(1) = DX1(0)
                  END IF
                  CLOSE(GRD%FHNDL)
               CASE (2)
                 OPEN(GRD%FHNDL, FILE = GRD%FNAME, STATUS = 'OLD')
                 IF (IGRIDTYPE == 1) THEN ! system.dat ... future now the stupid wwm format ... 
                   DO IP = 1, MNP
                      READ(GRD%FHNDL, *, IOSTAT = ISTAT) XP(IP), YP(IP), DEP(IP)
                      IF ( ISTAT /= 0 ) CALL ERROR_MSG('error in the grid file', 'PREPARE')
                   END DO
                   DO IE = 1, MNE
                      READ(GRD%FHNDL, *, IOSTAT = ISTAT) INE(:,IE)
                      IF ( ISTAT /= 0 )  CALL ERROR_MSG('error in the grid file', 'PREPARE')
                   END DO
                 ELSE IF (IGRIDTYPE == 3) THEN ! selfe gr3
                   READ(GRD%FHNDL,*) 
                   READ(GRD%FHNDL,*) ITMP, ITMP
                   DO IP = 1, MNP
                      READ(GRD%FHNDL, *, IOSTAT = ISTAT) ITMP, XPDTMP, YPDTMP, ZPDTMP
                      XP(IP)  = REAL(XPDTMP)
                      YP(IP)  = REAL(YPDTMP)
                      DEP(IP) = REAL(ZPDTMP)
!2do change code to double precision !!!
                      IF ( ISTAT /= 0 ) CALL ERROR_MSG('error in the grid file', 'PREPARE')
                   END DO
                   DO IE = 1, MNE
                      READ(GRD%FHNDL, *, IOSTAT = ISTAT) ITMP, ITMP, INE(:,IE)
                      IF ( ISTAT /= 0 )  CALL ERROR_MSG('error in the grid file', 'PREPARE')
                   END DO
! begin modification by MDS
                 ELSE IF (IGRIDTYPE == 2) THEN ! periodic grid written by mathieu dutour 
                   DO IP = 1, MNP
                      READ(GRD%FHNDL, *, IOSTAT = ISTAT) DEP(IP)
                      IF ( ISTAT /= 0 ) CALL ERROR_MSG('error in the grid file SELFE 3', 'PREPARE')
                   END DO
                   DO IE = 1, MNE
                      READ(GRD%FHNDL, *, IOSTAT = ISTAT) TRIA(IE)
                      IF ( ISTAT /= 0 )  CALL ERROR_MSG('error in the grid file WWMIDEAL 1', 'PREPARE')
                      READ(GRD%FHNDL, *, IOSTAT = ISTAT) INE(:,IE)
                      IF ( ISTAT /= 0 )  CALL ERROR_MSG('error in the grid file WWMIDEAL 2', 'PREPARE')
                      READ(GRD%FHNDL, *, IOSTAT = ISTAT) DXP1, DXP2, DXP3
                      IF ( ISTAT /= 0 ) CALL ERROR_MSG('error in the grid file WWMIDEAL 3', 'PREPARE')
                      READ(GRD%FHNDL, *, IOSTAT = ISTAT) DYP1, DYP2, DYP3
                      IF ( ISTAT /= 0 ) CALL ERROR_MSG('error in the grid file WWMIDEAL 4', 'PREPARE')
                      IEN(1,IE) = -DYP2
                      IEN(2,IE) = DXP2
                      IEN(3,IE) = -DYP3
                      IEN(4,IE) = DXP3
                      IEN(5,IE) = -DYP1
                      IEN(6,IE) = DXP1
                   END DO
                 END IF
                 CLOSE(GRD%FHNDL)
               CASE DEFAULT
                 STOP 'IGRIDTYPE WRONG'
            END SELECT
         END IF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CHECK_LOGICS()
         USE DATAPOOL
#ifdef SELFE
         use elfe_msgp, only : myrank,parallel_abort
#endif 
         IMPLICIT NONE

         REAL :: TEST

#ifdef SELFE
!        Check MSC,MDC for exchange
         if(MSC<1.or.MDC<1) call parallel_abort('MSC,MDC too small')
#endif

         IF (SMETHOD .GT. 0) THEN

           IF (MESIN .GT. 5) THEN
#ifdef SELFE
             call parallel_abort('CHECK NUMS - MESIN OUT OF RANGE')
#else
             WRITE(DBG%FHNDL,*) 'CHECK NUMS - MESIN OUT OF RANGE', MESIN
             STOP 'CHECK LOGICS'
#endif
           END IF

           IF (MESBR .GT. 1) THEN
#ifdef SELFE
             call parallel_abort('CHECK NUMS - MESBR OUT OF RANGE')
#else
             WRITE(DBG%FHNDL,*) 'CHECK NUMS - MESBR OUT OF RANGE', MESBR
             STOP 'CHECK LOGICS'
#endif
           END IF


           IF (MESBF .GT. 0 .AND. FRICC .LT. 0.) THEN
#ifdef SELFE
             call parallel_abort('CHECK NUMS - FRICTION COEFFICIENT HAS WRONG SIGN')
#else
             WRITE(DBG%FHNDL,*) 'CHECK NUMS - FRICTION COEFFICIENT HAS WRONG SIGN', MESBF
             STOP 'CHECK LOGICS'
#endif
           END IF

           IF (MESTR .GT. 2) THEN
#ifdef SELFE
             call parallel_abort('CHECK NUMS - MESTR OUT OF RANGE')
#else
             WRITE(DBG%FHNDL,*) 'CHECK NUMS - MESTR OUT OF RANGE', MESTR
             STOP 'CHECK LOGICS'
#endif
           END IF

           IF (MESNL .GT. 2) THEN
#ifdef SELFE
             call parallel_abort('CHECK NUMS - MESNL OUT OF RANGE')
#else
             WRITE(DBG%FHNDL,*) 'CHECK NUMS - MESNL OUT OF RANGE', MESNL
             STOP 'CHECK LOGICS'
#endif
           END IF

           IF (MESDS .GT. 5) THEN
#ifdef SELFE
             call parallel_abort('CHECK NUMS - MESDS OUT OF RANGE')
#else
             WRITE(DBG%FHNDL,*) 'CHECK NUMS - MESDS OUT OF RANGE', MESDS
             STOP 'CHECK LOGICS'
#endif
           END IF

           IF (MESNL .GT. 0 .AND. .NOT. LLIMT ) THEN
#ifdef WWMONLY
             STOP 'YOU ARE USING SNL WITHOUT LIMITER CODE WILL STOP NOW'
#elif SELFE
             call parallel_abort('SELFEWWM: error, YOU ARE USING SNL WITHOUT LIMITER CODE WILL STOP NOW ... :(....')
#endif
           END IF 

         ELSE

           MESNL = 0
           MESIN = 0
           MESDS = 0
           MESBR = 0
           MESTR = 0
           MESBR = 0
           MESBF = 0
           LLIMT = .FALSE.

         END IF

         IF (LBCWA .OR. LBCSP) THEN
           IF (PGIVE(7) .LT. THR) THEN
             PGIVE(7) = 0.1
           ELSE IF (PGIVE(8) .LT. THR) THEN
             PGIVE(8) = 3.3
           END IF
         END IF

#ifdef WWMONLY
         IF (LCPL) THEN
           IF (.NOT. LROMS .AND. .NOT. LSHYFEM .AND. .NOT. LTIMOR) THEN
              WRITE(DBG%FHNDL,*) 'LROMS, LSHYFEM or LTIMOR must be true'
              STOP 'LROMS, LSHYFEM or LTIMOR must be true'
           END IF
           IF (MAIN%DTCOUP .LT. MAIN%DELT) THEN
             WRITE(DBG%FHNDL,*) 'COUPLE TIME STEP IS SMALLER AS THE CALCULATION TIME STEP!'
             STOP 'COUPLE TIME STEP IS SMALLER AS THE CALCULATION TIME STEP!'
           END IF
           TEST = MAIN%DTCOUP - NINT(MAIN%DTCOUP/MAIN%DELT)*MAIN%DELT
!2do ... check where else you do some nint stuff ... like that one ...
           IF (ABS(TEST) .GT. THR) THEN
             WRITE(DBG%FHNDL,*) 'MAIN%DTCOUP=', MAIN%DTCOUP
             WRITE(DBG%FHNDL,*) 'MAIN%DELT=', MAIN%DELT
             WRITE(DBG%FHNDL,*) 'TEST=', TEST
             WRITE(DBG%FHNDL,*) 'TIME STEP OF THE WAVEMODELL CANNOT BE DiVIDIED WITHOUT A REST'
             STOP 'TIME STEP OF THE WAVEMODELL CANNOT BE DiVIDIED WITHOUT A REST'
           ELSE
             MAIN%ICPLT = INT(MAIN%DTCOUP/MAIN%DELT)
           END IF
         END IF
#endif 

         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'SWTICHES FOR THE LIMTER'
         WRITE(STAT%FHNDL,*) 'LLIMT', LLIMT
         WRITE(STAT%FHNDL,'("+TRACE...",A)') 'ACTIVATED SOURCE TERMS'
         WRITE(STAT%FHNDL,*) 'MESIN', MESIN
         WRITE(STAT%FHNDL,*) 'MESNL', MESNL
         WRITE(STAT%FHNDL,*) 'MESBR', MESBR
         WRITE(STAT%FHNDL,*) 'MESDS', MESDS
         WRITE(STAT%FHNDL,*) 'MESTR', MESTR

         IF (LSEWD .AND. LSTWD) THEN
           WRITE(DBG%FHNDL,*) 'YOU MUST USE EITHER UNSTEADY OR STEADY WIND'
           WRITE(DBG%FHNDL,*) 'PLEASE CHECK CODE EXITS'
           STOP 'CHECK LSEWL OR LSTDW'
         END IF

         IF (LSTCU .AND. LSECU) THEN
           WRITE(DBG%FHNDL,*) 'YOU MUST USE EITHER UNSTEADY OR STEADY CURRENTS'
           WRITE(DBG%FHNDL,*) 'PLEASE CHECK CODE EXITS'
           STOP 'CHECK LSTCU .AND. LSECU'
         END IF

         IF (LSTCU .AND. LSECU) THEN
           WRITE(DBG%FHNDL,*) 'YOU MUST USE EITHER UNSTEADY OR STEADY CURRENTS'
           WRITE(DBG%FHNDL,*) 'PLEASE CHECK CODE EXITS'
           STOP 'CHECK LSTCU .AND. LSECU'
         END IF

         IF (LSTWL .AND. LSEWL) THEN
           WRITE(DBG%FHNDL,*) 'YOU MUST USE EITHER UNSTEADY OR STEADY CURRENTS'
           WRITE(DBG%FHNDL,*) 'PLEASE CHECK CODE EXITS'
           STOP 'CHECK LSTCU .AND. LSECU'
         END IF

         RETURN

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WRITE_TEST_OUTPUT()
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER :: IS, ID

         IF (LTEST) THEN

         END IF
!
!     *** output format
!

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef NCDF
      SUBROUTINE INIT_NETCDF_DWD
         USE DATAPOOL
         USE NETCDF 
         IMPLICIT NONE

        INTEGER :: ISTAT, IT, IX, IY, IFILE
        INTEGER :: ILON_ID, ILAT_ID, ITIME_ID, I, J, COUNTER
        INTEGER :: D_WIND_X_ID, D_WIND_Y_ID, D_PRESS_ID
        REAL    :: RTMP
        REAL*8  :: DTMP
        REAL*8, ALLOCATABLE :: WIND_TIME(:)
        character ( len = 20 ) chrtmp 
        character ( len = 15 ) chrdate
        character ( len = 40 ) netcfd_fname

        integer, dimension(nf90_max_var_dims) :: dimIDs

        OPEN(WIN%FHNDL,FILE=WIN%FNAME,STATUS='OLD')
!
! count number of netcdf files in list ...
!
        NUM_NETCDF_FILES = 0
        DO 
          READ( WIN%FHNDL, *, IOSTAT = ISTAT ) 
          IF ( ISTAT /= 0 ) EXIT
          NUM_NETCDF_FILES = NUM_NETCDF_FILES + 1
        END DO 
        REWIND (WIN%FHNDL)

        ALLOCATE(NETCDF_FILE_NAMES(NUM_NETCDF_FILES))

        DO IT = 1, NUM_NETCDF_FILES
          READ( WIN%FHNDL, *) NETCDF_FILE_NAMES(IT)
!          WRITE(*,*) IT, NETCDF_FILE_NAMES(IT)
        END DO 
        CLOSE (WIN%FHNDL)
!
! check number of time steps in netcdf file ... it is assumed that all files have the same ammount of time steps ...
!
        ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(1), NF90_NOWRITE, WIND_NCID)
        ISTAT = nf90_inq_varid(WIND_NCID, 'initial_time0_encoded', ITIME_ID)
        ISTAT = NF90_INQUIRE_VARIABLE(WIND_NCID, ITIME_ID, dimids = dimids)
        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(1), len = NDT_WIND_FILE)
!
! check dimensions in the netcdf ... again it is assumed that this is not changing for all files ...
!
        ISTAT = nf90_inq_varid(WIND_NCID, 'g0_lon_2', ILON_ID)
        ISTAT = NF90_INQUIRE_VARIABLE(WIND_NCID, ILON_ID, dimids = dimIDs)
        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(1), len = NDX_WIND)

        ISTAT = nf90_inq_varid(WIND_NCID, 'g0_lat_1', ILAT_ID)
        ISTAT = NF90_INQUIRE_VARIABLE(WIND_NCID, ILAT_ID, dimids = dimIDs)
        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(1), len = NDY_WIND)

        ALLOCATE (COORD_WIND_X(NDX_WIND))
        ALLOCATE (COORD_WIND_Y(NDY_WIND))
!
! read cooridantes from files ....
!
        ISTAT = NF90_GET_VAR(WIND_NCID, ILON_ID, COORD_WIND_X)
        ISTAT = NF90_GET_VAR(WIND_NCID, ILAT_ID, COORD_WIND_Y)
!
! estimate offset ...
!
        OFFSET_X_WIND = MINVAL(COORD_WIND_X)
        OFFSET_Y_WIND = MINVAL(COORD_WIND_Y)
!
! resolution ...
!
        DX_WIND  = ABS(MAXVAL(COORD_WIND_X)-MINVAL(COORD_WIND_X))/(NDX_WIND-1)
        DY_WIND  = ABS(MAXVAL(COORD_WIND_Y)-MINVAL(COORD_WIND_Y))/(NDY_WIND-1)
!
! close netcdf file ...
!
        ISTAT = NF90_CLOSE(WIND_NCID)
!
! total number of time steps ... in all files 
!
        NDT_WIND_ALL_FILES = NDT_WIND_FILE * NUM_NETCDF_FILES

        ALLOCATE (WIND_TIME(NDT_WIND_FILE))
        ALLOCATE (WIND_TIME_ALL_FILES(NDT_WIND_ALL_FILES))
!
! read all time steps in the proper format and transform in wwm time line 
!
        DO IFILE = 1, NUM_NETCDF_FILES 
          ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(IFILE), NF90_NOWRITE, WIND_NCID)
          ISTAT = NF90_GET_VAR(WIND_NCID, ITIME_ID, WIND_TIME)
          DO IT = 1, NDT_WIND_FILE
            WRITE (3001, *) WIND_TIME(IT)
          END DO
          REWIND(3001)
          DO IT = 1, NDT_WIND_FILE
            READ (3001, '(A20)') CHRTMP
            ! YYYYMMDD
            CHRDATE(1:8) = CHRTMP(4:12)
            CHRDATE(9:9) = '.'
            ! HH
            CHRDATE(10:11) = CHRTMP(12:13)
            ! MMSS
            CHRDATE(12:15) = '0000' ! Construct propper format YYYYMMDDHHMMSS character .... len = 15 ... :)
            CALL CT2MJD(CHRDATE,DTMP) ! 
            WIND_TIME(IT) = DTMP    ! Double time with respect to 19000101.000000
          END DO
          CLOSE(3001)
          DO IT = 1, NDT_WIND_FILE
            WIND_TIME_ALL_FILES(IT+(IFILE-1)*NDT_WIND_FILE) = WIND_TIME(IT)
          END DO
        END DO ! IFILE

        SEWI%DELT = (WIND_TIME_ALL_FILES(2) - WIND_TIME_ALL_FILES(1)) * DAY2SEC 

!        DO IFILE = 1, NUM_NETCDF_FILES
!          DO IT = 1, NDT_WIND_FILE
!            WRITE(*,*) IFILE, IT, WIND_TIME_ALL_FILES(IT+(IFILE-1)*NDT_WIND_FILE)
!          END DO
!        END DO

        IF (LWRITE_ORIG_WIND) THEN
          WRITE (3010, '(I10)') 0
          WRITE (3010, '(I10)') NDX_WIND * NDY_WIND 
          COUNTER = 0
          DO I = 1, NDY_WIND 
            DO J = 1, NDX_WIND
              WRITE (3010, '(I10,3F15.4)') COUNTER, OFFSET_X_WIND+(J-1)*DX_WIND ,OFFSET_Y_WIND+(I-1)*DY_WIND , 0.0
              COUNTER = COUNTER + 1
            END DO
          END DO
          WRITE (3010, *) (NDX_WIND-1)*(NDY_WIND-1)*2
          DO J = 0, NDY_WIND-2
            DO I = 0, NDX_WIND-2
              WRITE (3010, '(5I10)')  I+J*NDX_WIND           , NDX_WIND+I+J* NDX_WIND, NDX_WIND+I+1+J*NDX_WIND, 0, 0
              WRITE (3010, '(5I10)')  NDX_WIND+I+1+J*NDX_WIND, I+1+J*NDX_WIND        , I+J*NDX_WIND           , 0, 0
            END DO
          END DO
          OPEN(3011, FILE  = 'ergwindorig.bin', FORM = 'UNFORMATTED')
        END IF
          
     END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_NETCDF_DWD(IT)
         USE DATAPOOL, ONLY : WIND_NCID, WIND_X, WIND_Y, ATMO_PRESS, LINVERTY, NDX_WIND, NDY_WIND, LWRITE_ORIG_WIND
         USE NETCDF
         IMPLICIT NONE
!
!        READS WIND_Y, WIND_X and PRESSURE from a given NCID within one DWD file 
!
         INTEGER, INTENT(IN) :: IT

         INTEGER             :: DWIND_X_ID, DWIND_Y_ID, DPRESS_ID, ISTAT
         INTEGER             :: numLons, numLats, numTime, iy, counter, ip, i, j
         REAL,   ALLOCATABLE :: TMP(:,:)
         REAL, ALLOCATABLE   :: U(:), V(:), H(:)
         REAL, SAVE          :: TIME

         INTEGER, DIMENSION (nf90_max_var_dims) :: dimIDs

         ISTAT = nf90_inq_varid(WIND_NCID, 'V_GDS0_HTGL_13', DWIND_X_ID)
         ISTAT = NF90_INQUIRE_VARIABLE(WIND_NCID, DWIND_X_ID, dimids = dimIDs)
         ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(1), len = numLons)
         ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(2), len = numLats)
         ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(3), len = numTime)
         IF (.NOT. ALLOCATED(WIND_X)) ALLOCATE (WIND_X(NDX_WIND,NDY_WIND))

         ISTAT = nf90_inq_varid(WIND_NCID, 'U_GDS0_HTGL_13', DWIND_Y_ID)
         ISTAT = NF90_INQUIRE_VARIABLE(WIND_NCID, DWIND_Y_ID, dimids = dimIDs)
         ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(1), len = numLons)
         ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(2), len = numLats)
         ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(3), len = numTime)
         IF (.NOT. ALLOCATED(WIND_Y)) ALLOCATE (WIND_Y(NDX_WIND,NDY_WIND))

!         ISTAT = nf90_inq_varid(WIND_NCID, 'PS_MSL_GDS0_MSL_13', DPRESS_ID)
!         ISTAT = NF90_INQUIRE_VARIABLE(WIND_NCID, DPRESS_ID, dimids = dimIDs)
!         ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(1), len = numLons)
!         ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(2), len = numLats)
!         ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(3), len = numTime)
!         IF (.NOT. ALLOCATED(ATMO_PRESS)) ALLOCATE (ATMO_PRESS(NDX_WIND,NDY_WIND))

         ISTAT = NF90_GET_VAR(WIND_NCID, DWIND_X_ID, WIND_X,    start = (/ 1, 1, IT /), count = (/ NDX_WIND, NDY_WIND, 1 /))
         IF (ISTAT .NE. nf90_noerr) STOP 'ERR READING D_WIND_X_ID'
         ISTAT = NF90_GET_VAR(WIND_NCID, DWIND_Y_ID, WIND_Y,    start = (/ 1, 1, IT /), count = (/ NDX_WIND, NDY_WIND, 1 /))
         IF (ISTAT .NE. nf90_noerr) STOP 'ERR READING D_WIND_Y_ID'
!         ISTAT = NF90_GET_VAR(WIND_NCID, DPRESS_ID, ATMO_PRESS, start = (/ 1, 1, IT /), count = (/ NDX_WIND, NDY_WIND, 1 /))
!         IF (ISTAT .NE. nf90_noerr) STOP 'ERR READING D_PRESS_ID'

         IF (LINVERTY) THEN
           ALLOCATE(TMP(NDX_WIND,NDY_WIND))
           DO IY = 1, NDY_WIND 
             tmp(:,NDY_WIND-(IY-1)) = wind_x(:,IY)
           END DO
           wind_x = tmp 
           DO IY = 1, NDY_WIND
             tmp(:,NDY_WIND-(IY-1)) = wind_y(:,IY)
           END DO
           wind_y = tmp
           DO IY = 1, NDY_WIND 
             tmp(:,NDY_WIND-(IY-1)) = atmo_press(:,IY)
           END DO
           atmo_press = tmp
           DEALLOCATE(TMP)
         END IF

         IF (.NOT. ALLOCATED(U)) ALLOCATE(U(NDX_WIND*NDY_WIND))
         IF (.NOT. ALLOCATED(V)) ALLOCATE(V(NDX_WIND*NDY_WIND))
         IF (.NOT. ALLOCATED(H)) ALLOCATE(H(NDX_WIND*NDY_WIND))

         COUNTER = 1
         DO J = 1, NDY_WIND 
           DO I = 1, NDX_WIND 
             U(COUNTER) = WIND_X(I,J)
             V(COUNTER) = WIND_Y(I,J)
             IF (ABS(U(COUNTER)) .GT. 1000.) U(COUNTER) = 0.
             IF (ABS(V(COUNTER)) .GT. 1000.) V(COUNTER) = 0.
             H(COUNTER) = SQRT((U(COUNTER)**2.+V(COUNTER)**2.))
             COUNTER = COUNTER + 1
           END DO
         END DO

        IF (LWRITE_ORIG_WIND) THEN

          IF (.NOT. ALLOCATED(U)) ALLOCATE(U(NDX_WIND*NDY_WIND))
          IF (.NOT. ALLOCATED(V)) ALLOCATE(V(NDX_WIND*NDY_WIND))
          IF (.NOT. ALLOCATED(H)) ALLOCATE(H(NDX_WIND*NDY_WIND))

          COUNTER = 1
          DO J = 1, NDY_WIND
            DO I = 1, NDX_WIND
             U(COUNTER) = WIND_X(I,J)
             V(COUNTER) = WIND_Y(I,J)
             IF (ABS(U(COUNTER)) .GT. 1000.) U(COUNTER) = 0.
             IF (ABS(V(COUNTER)) .GT. 1000.) V(COUNTER) = 0.
             H(COUNTER) = SQRT((U(COUNTER)**2.+V(COUNTER)**2.))
             COUNTER = COUNTER + 1
            END DO
          END DO

          TIME = TIME + 1.
          WRITE(3011) TIME 
          WRITE(3011) (U(IP), V(IP), H(IP), IP = 1, numLons*numLats)

        END IF

      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
     SUBROUTINE INIT_NETCDF_NARR
         USE DATAPOOL
         USE NETCDF 
         IMPLICIT NONE

!for data description consult ftp://nomads.ncdc.noaa.gov/CFSR/HP_time_series/200307/wnd10m.l.gdas.200307.grb2.inv

        INTEGER :: ISTAT, IT, IX, IY, IFILE
        INTEGER :: ILON_ID, ILAT_ID, ITIME_ID, I, J, COUNTER
        INTEGER :: D_WIND_X_ID, D_WIND_Y_ID, D_PRESS_ID
        REAL    :: RTMP
        REAL*8  :: START_TIME 
        REAL*8, ALLOCATABLE :: WIND_TIME(:), WIND_TIME_NETCDF(:)
        character ( len = 20 ) chrtmp 
        character ( len = 15 ) chrdate
        character ( len = 40 ) netcfd_fname, beginn_time

        integer, dimension(nf90_max_var_dims) :: dimIDs

        OPEN(WIN%FHNDL,FILE=WIN%FNAME,STATUS='OLD',IOSTAT = ISTAT)
        IF (ISTAT .LT. 0) THEN
          WRITE(*,*) 'WINDFILES FILE NOT FOUND', WIN%FNAME
          STOP
        END IF
!
! count number of netcdf files in list ...
!
        NUM_NETCDF_FILES = 0
        DO 
          READ( WIN%FHNDL, *, IOSTAT = ISTAT ) 
          IF ( ISTAT /= 0 ) EXIT
          NUM_NETCDF_FILES = NUM_NETCDF_FILES + 1
        END DO 
        REWIND (WIN%FHNDL)

        ALLOCATE(NETCDF_FILE_NAMES(NUM_NETCDF_FILES))

        DO IT = 1, NUM_NETCDF_FILES
          READ( WIN%FHNDL, *) NETCDF_FILE_NAMES(IT)
!          WRITE(*,*) IT, NETCDF_FILE_NAMES(IT)
        END DO 
        CLOSE (WIN%FHNDL)
!
! check number of time steps in netcdf file ... it is assumed that all files have the same ammount of time steps ...
!
        ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(1), NF90_NOWRITE, WIND_NCID)
        ISTAT = nf90_inq_varid(WIND_NCID, 'time', ITIME_ID)
        ISTAT = NF90_INQUIRE_VARIABLE(WIND_NCID, ITIME_ID, dimids = dimids)
        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(1), len = NDT_WIND_FILE)
        ISTAT = nf90_get_att(WIND_NCID, ITIME_ID, 'units', beginn_time)

        CHRDATE(1:4) = beginn_time(13:17)
        CHRDATE(5:6) = beginn_time(18:19)
        CHRDATE(7:8) = beginn_time(21:22)
        CHRDATE(9:9)   = '.'
        CHRDATE(10:15)= '000000'
        CALL CT2MJD(CHRDATE,START_TIME)
!
! check dimensions in the netcdf ... again it is assumed that this is not changing for all files ...
!
        ISTAT = nf90_inq_varid(WIND_NCID, 'lon', ILON_ID)
        ISTAT = NF90_INQUIRE_VARIABLE(WIND_NCID, ILON_ID, dimids = dimIDs)
        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(1), len = NDX_WIND)

        ISTAT = nf90_inq_varid(WIND_NCID, 'lat', ILAT_ID)
        ISTAT = NF90_INQUIRE_VARIABLE(WIND_NCID, ILAT_ID, dimids = dimIDs)
        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(1), len = NDY_WIND)

        ALLOCATE (DCOORD_WIND_X(NDX_WIND))
        ALLOCATE (DCOORD_WIND_Y(NDY_WIND))
!
! read cooridantes from files ....
!
        ISTAT = NF90_GET_VAR(WIND_NCID, ILON_ID, DCOORD_WIND_X)
        ISTAT = NF90_GET_VAR(WIND_NCID, ILAT_ID, DCOORD_WIND_Y)

        DCOORD_WIND_Y = DCOORD_WIND_Y !+ 90.d0

        DO I = 1, NDX_WIND
          DCOORD_WIND_X(I) = DCOORD_WIND_X(I) - 180.d0
        END DO 
!
! estimate offset ...
!
        OFFSET_X_WIND = MINVAL(DCOORD_WIND_X)
        OFFSET_Y_WIND = MINVAL(DCOORD_WIND_Y)
!
! resolution ...
!
        DX_WIND  = ABS(MAXVAL(DCOORD_WIND_X)-MINVAL(DCOORD_WIND_X))/(NDX_WIND-1)
        DY_WIND  = ABS(MAXVAL(DCOORD_WIND_Y)-MINVAL(DCOORD_WIND_Y))/(NDY_WIND-1)
!
! close netcdf file ...
!
        ISTAT = NF90_CLOSE(WIND_NCID)
!
! total number of time steps ... in all files 
!
        NDT_WIND_ALL_FILES = NDT_WIND_FILE * NUM_NETCDF_FILES

        ALLOCATE (WIND_TIME(NDT_WIND_FILE),WIND_TIME_NETCDF(NDT_WIND_FILE))
        ALLOCATE (WIND_TIME_ALL_FILES(NDT_WIND_ALL_FILES))
!
! read all time steps in the proper format and transform in wwm time line 
!
        DO IFILE = 1, NUM_NETCDF_FILES 
          ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(IFILE), NF90_NOWRITE, WIND_NCID)
          ISTAT = NF90_GET_VAR(WIND_NCID, ITIME_ID, WIND_TIME_NETCDF)
!AR: Here is the error CFRS is hourly not 6 hourly ...
          DO IT = 1, NDT_WIND_FILE
            WIND_TIME(IT) = START_TIME + WIND_TIME_NETCDF(IT) * 3600. * SEC2DAY ! We write all the times from the netcdf into the arrayq
                                                                                ! However the 1st entry of each block (7 entries for 6 hours) 
                                                                                ! is the initial condition and NOT the analysis! 
            !WRITE(*,'(2I10,3F15.8)') IT, NDT_WIND_FILE, START_TIME, WIND_TIME_NETCDF(IT) * 3600. * SEC2DAY, WIND_TIME(IT)
          END DO
          DO IT = 1, NDT_WIND_FILE
            WIND_TIME_ALL_FILES(IT+(IFILE-1)*NDT_WIND_FILE) = WIND_TIME(IT)
            !WRITE(*,*) IT, IFILE, WIND_TIME_ALL_FILES(IT+(IFILE-1)*NDT_WIND_FILE)
          END DO
        END DO ! IFILE

        SEWI%DELT = ( WIND_TIME_ALL_FILES(2) - WIND_TIME_ALL_FILES(1) ) * DAY2SEC 

        !WRITE(*,*) WIND_TIME_ALL_FILES(2), WIND_TIME_ALL_FILES(1), DAY2SEC, SEWI%DELT

        IF (LWRITE_ORIG_WIND) THEN
          WRITE (3010, '(I10)') 0
          WRITE (3010, '(I10)') NDX_WIND * NDY_WIND 
          COUNTER = 0
          DO I = 1, NDY_WIND 
            DO J = 1, NDX_WIND
              WRITE (3010, '(I10,3F15.4)') COUNTER, OFFSET_X_WIND+(J-1)*DX_WIND ,OFFSET_Y_WIND+(I-1)*DY_WIND , 0.0
              COUNTER = COUNTER + 1
            END DO
          END DO
          WRITE (3010, *) (NDX_WIND-1)*(NDY_WIND-1)*2
          DO J = 0, NDY_WIND-2
            DO I = 0, NDX_WIND-2
              WRITE (3010, '(5I10)')  I+J*NDX_WIND           , NDX_WIND+I+J* NDX_WIND, NDX_WIND+I+1+J*NDX_WIND, 0, 0
              WRITE (3010, '(5I10)')  NDX_WIND+I+1+J*NDX_WIND, I+1+J*NDX_WIND        , I+J*NDX_WIND           , 0, 0
            END DO
          END DO
          OPEN(3011, FILE  = 'ergwindorig.bin', FORM = 'UNFORMATTED')
        END IF

     END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
     SUBROUTINE INIT_NETCDF_NARR2
         USE DATAPOOL
         USE NETCDF 
         IMPLICIT NONE
!
!2do update ...
!for data description consult ftp://nomads.ncdc.noaa.gov/CFSR/HP_time_series/200307/wnd10m.l.gdas.200307.grb2.inv
!
        INTEGER :: ISTAT, IT, IX, IY, IFILE, II, IP, III
        INTEGER :: ILON_ID, ILAT_ID, ITIME_ID, I, J, COUNTER
        INTEGER :: D_WIND_X_ID, D_WIND_Y_ID, D_PRESS_ID
        INTEGER*4 :: TIMES
        REAL    :: RTMP
        REAL*8  :: START_TIME, OFFSET_TIME
        REAL*8, ALLOCATABLE :: WIND_TIME_NETCDF(:)
        character ( len = 4 ) ch4
        character ( len = 100 ) chrtmp 
        character ( len = 15 ) chrdate
        character ( len = 40 ) netcfd_fname, beginn_time

        integer, dimension(nf90_max_var_dims) :: dimIDs

        OPEN(WIN%FHNDL,FILE=WIN%FNAME,STATUS='OLD',IOSTAT = ISTAT)
        IF (ISTAT .LT. 0) THEN
          WRITE(*,*) 'WINDFILES FILE NOT FOUND', WIN%FNAME
          STOP
        END IF
!
! count number of netcdf files in list ...
!
        NUM_NETCDF_FILES = 0
        DO 
          READ( WIN%FHNDL, *, IOSTAT = ISTAT ) 
          IF ( ISTAT /= 0 ) EXIT
          NUM_NETCDF_FILES = NUM_NETCDF_FILES + 1
        END DO 
        REWIND (WIN%FHNDL)

        ALLOCATE(NETCDF_FILE_NAMES(NUM_NETCDF_FILES))

        DO IT = 1, NUM_NETCDF_FILES
          READ( WIN%FHNDL, *) NETCDF_FILE_NAMES(IT)
!          WRITE(*,*) IT, NETCDF_FILE_NAMES(IT)
        END DO 
        CLOSE (WIN%FHNDL)
!
! check number of time steps in netcdf file ... it is assumed that all files have the same ammount of time steps ...
!
        ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(1), NF90_NOWRITE, WINDX_NCID)
         IF (ISTAT .NE. nf90_noerr) THEN
           CHRTMP = nf90_strerror(ISTAT)
           WRITE(*,*) 'NETCDF ERROR -1-', CHRTMP
           STOP 'NETCDF WW3'
         ENDIF
        ISTAT = nf90_inq_varid(WINDX_NCID, 'time', ITIME_ID)
         IF (ISTAT .NE. nf90_noerr) THEN
           CHRTMP = nf90_strerror(ISTAT)
           WRITE(*,*) 'NETCDF ERROR -2-', CHRTMP
           STOP 'NETCDF WW3'
         ENDIF
        ISTAT = NF90_INQUIRE_VARIABLE(WINDX_NCID, ITIME_ID, dimids = dimIDs)
         IF (ISTAT .NE. nf90_noerr) THEN
           CHRTMP = nf90_strerror(ISTAT)
           WRITE(*,*) 'NETCDF ERROR -3-', CHRTMP
           STOP 'NETCDF WW3'
         ENDIF
        ISTAT = nf90_inquire_dimension(WINDX_NCID, dimIDs(1), len = NDT_WIND_FILE)
         IF (ISTAT .NE. nf90_noerr) THEN
           CHRTMP = nf90_strerror(ISTAT)
           WRITE(*,*) 'NETCDF ERROR -4-', CHRTMP
           STOP 'NETCDF WW3'
         ENDIF
         WRITE(*,*) NDT_WIND_FILE, 'NDT_WIND_FILE'

         ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(2), NF90_NOWRITE, WINDY_NCID)
         IF (ISTAT .NE. nf90_noerr) THEN
           CHRTMP = nf90_strerror(ISTAT)
           WRITE(*,*) 'NETCDF ERROR -5-', CHRTMP
           STOP 'NETCDF WW3'
         ENDIF
!
        !WRITE(*,*) SEBO%BEGT, NDT_WIND_FILE
!
! I make the assumption that the year when the dataset beginns at the year indicated in the bouc section 
! 
        CHRDATE(1:4) = SEBO%BEGT(1:4) 
        CHRDATE(5:6) = '01' 
        CHRDATE(7:8) = '01' 
        CHRDATE(9:9)   = '.'
        CHRDATE(10:15)= '000000'
        CALL CT2MJD(CHRDATE,START_TIME)
!
! the following years are automatically adjusted ...
!

        !WRITE(*,*) CHRDATE, START_TIME
!
! check dimensions in the netcdf ... again it is assumed that this is not changing for all files ...
!
        ISTAT = nf90_inq_varid(WINDX_NCID, 'lon', ILON_ID)
        ISTAT = NF90_INQUIRE_VARIABLE(WINDX_NCID, ILON_ID, dimids = dimIDs)
        ISTAT = nf90_inquire_dimension(WINDX_NCID, dimIDs(1), len = NDX_WIND)
        ISTAT = nf90_inquire_dimension(WINDX_NCID, dimIDs(2), len = NDY_WIND)

        ISTAT = nf90_inq_varid(WINDY_NCID, 'lat', ILAT_ID)
        ISTAT = NF90_INQUIRE_VARIABLE(WINDY_NCID, ILAT_ID, dimids = dimIDs)
        ISTAT = nf90_inquire_dimension(WINDY_NCID, dimIDs(1), len = NDX_WIND)
        ISTAT = nf90_inquire_dimension(WINDY_NCID, dimIDs(2), len = NDY_WIND)

        ALLOCATE (DCOORD_WIND_X2(NDX_WIND,NDY_WIND))
        ALLOCATE (DCOORD_WIND_Y2(NDX_WIND,NDY_WIND))
!
! read cooridantes from files ....
!
        ISTAT = NF90_GET_VAR(WINDX_NCID, ILON_ID, DCOORD_WIND_X2)
        ISTAT = NF90_GET_VAR(WINDY_NCID, ILAT_ID, DCOORD_WIND_Y2)

        !WRITE(*,*) DCOORD_WIND_X2, DCOORD_WIND_Y2

!        DCOORD_WIND_Y = DCOORD_WIND_Y !+ 90.d0

!        DO I = 1, NDX_WIND
!          DCOORD_WIND_X(I) = DCOORD_WIND_X(I) - 180.d0
!        END DO 
!
! estimate offset ...
!
        OFFSET_X_WIND = MINVAL(DCOORD_WIND_X2)
        OFFSET_Y_WIND = MINVAL(DCOORD_WIND_Y2)

        !WRITE(*,*) MINVAL(DCOORD_WIND_X2), MINVAL(DCOORD_WIND_Y2), MAXVAL(DCOORD_WIND_X2), MAXVAL(DCOORD_WIND_Y2) 
!
! close netcdf file ...
!
        ISTAT = NF90_CLOSE(WIND_NCID)
!
! total number of time steps ... in all files 
!
        NDT_WIND_ALL_FILES = NDT_WIND_FILE * NUM_NETCDF_FILES/2

        ALLOCATE (WIND_TIME_NETCDF(NDT_WIND_FILE))
        ALLOCATE (WIND_TIME_ALL_FILES(NDT_WIND_ALL_FILES))
!
! read all time steps in the proper format and transform in wwm time line 
!
        DO IFILE = 1, NUM_NETCDF_FILES, 2
          ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(IFILE), NF90_NOWRITE, WINDX_NCID)
          ISTAT = NF90_GET_VAR(WINDX_NCID, ITIME_ID, WIND_TIME_NETCDF)
          WIND_TIME_NETCDF = WIND_TIME_NETCDF/24 ! Transform to days ...
          OFFSET_TIME = MINVAL(WIND_TIME_NETCDF) - START_TIME  ! in this dataset the time start at 1800 in WWM it start at 1900 
          WIND_TIME_NETCDF = WIND_TIME_NETCDF - OFFSET_TIME ! Now in WWM timeline ...
          DO IT = 1, NDT_WIND_FILE
            WIND_TIME_ALL_FILES(IT+(IFILE-1)*NDT_WIND_FILE) = WIND_TIME_NETCDF(IT)
          END DO
          CH4 = CHRDATE(1:4)
          READ (CH4,'(i4)') I 
          I = I + 1
          WRITE (11111,'(i4)') I 
          REWIND(11111)
          READ(11111,*) ch4
          CHRDATE(1:4) = ch4
          CHRDATE(5:6) = '01'
          CHRDATE(7:8) = '01'
          CHRDATE(9:9)   = '.'
          CHRDATE(10:15)= '000000'
          CALL CT2MJD(CHRDATE,START_TIME)
          !WRITE(*,*) CHRDATE, START_TIME
        END DO ! IFILE

        SEWI%DELT = ( WIND_TIME_ALL_FILES(2) - WIND_TIME_ALL_FILES(1) ) * DAY2SEC 

        !WRITE(*,*) WIND_TIME_ALL_FILES(2), WIND_TIME_ALL_FILES(1), DAY2SEC, SEWI%DELT

        MNE_WIND = (NDX_WIND-1)*(NDY_WIND-1)*2
        MNP_WIND =  NDX_WIND*NDY_WIND

        IF (LWRITE_ORIG_WIND) THEN
          WRITE (3010, '(I10)') 0
          WRITE (3010, '(I10)') MNP_WIND 
        END IF
        COUNTER = 0
        ALLOCATE(XYPWIND(2,NDX_WIND*NDY_WIND))
        DO I = 1, NDX_WIND 
          DO J = 1, NDY_WIND
            IF (DCOORD_WIND_X2(I,J) .GT. 0.) THEN
              DCOORD_WIND_X2(I,J) = -1 * DCOORD_WIND_X2(I,J) - (180.-DCOORD_WIND_X2(I,J)) * 2 ! Transformed to a unified domain extending below -180. 
            END IF
            COUNTER = COUNTER + 1
            XYPWIND(1,COUNTER) = DCOORD_WIND_X2(I,J)
            XYPWIND(2,COUNTER) = DCOORD_WIND_Y2(I,J)
            !WRITE(*,*) COUNTER, XYPWIND(1,COUNTER), XYPWIND(2,COUNTER)
            IF (LWRITE_ORIG_WIND) WRITE (3010, '(I10,3F15.4)') COUNTER, XYPWIND(1,COUNTER), XYPWIND(2,COUNTER), 0.0
          END DO
        END DO

        IF (LWRITE_ORIG_WIND) WRITE (3010, *) MNE_WIND 
        MNE_WIND = (NDX_WIND-1)*(NDY_WIND-1)*2
        ALLOCATE(INE_WIND(3,MNE_WIND)); INE_WIND = 0

        ALLOCATE(UWND_NARR(NDX_WIND*NDY_WIND)); UWND_NARR = 0.
        ALLOCATE(VWND_NARR(NDX_WIND*NDY_WIND)); VWND_NARR = 0.

        II = 0
        III = 0
        DO J = 0, NDX_WIND-2
          DO I = 0, NDY_WIND-2
            IF (MOD(J,2) .EQ. 0) THEN
!            IF (MOD(III,2) .EQ. 1) THEN
              II = II + 1
              INE_WIND(1,II) = I+J*NDY_WIND
              INE_WIND(2,II) = NDY_WIND+I+J* NDY_WIND
              INE_WIND(3,II) = NDY_WIND+I+1+J*NDY_WIND
              !WRITE(*,*) INE_WIND(:,II)
              II = II + 1
              INE_WIND(1,II) = NDY_WIND+I+1+J*NDY_WIND
              INE_WIND(2,II) = I+1+J*NDY_WIND
              INE_WIND(3,II) = I+J*NDY_WIND
!            ELSE
!              II = II + 1
!              INE_WIND(1,II) = I+J*NDY_WIND !I+J*NDY_WIND
!              INE_WIND(2,II) = NDY_WIND+I+J* NDY_WIND !I+1+J*NDY_WIND !NDY_WIND+I+J* NDY_WIND
!              INE_WIND(3,II) = I+1+J*NDY_WIND
!              !WRITE(*,*) INE_WIND(:,II)
!              II = II + 1
!              INE_WIND(1,II) = NDY_WIND+I+1+J*NDY_WIND
!              INE_WIND(2,II) = I+1+J*NDY_WIND
!              INE_WIND(3,II) = NDY_WIND+I+J* NDY_WIND
!              !WRITE(*,*) INE_WIND(:,II)
!            END IF
            ELSE
!            IF (MOD(III,2) .EQ. 0) THEN
!              II = II + 1
!              INE_WIND(1,II) = I+J*NDY_WIND
!              INE_WIND(2,II) = NDY_WIND+I+J* NDY_WIND
!              INE_WIND(3,II) = NDY_WIND+I+1+J*NDY_WIND
!              !WRITE(*,*) INE_WIND(:,II)
!              II = II + 1
!              INE_WIND(1,II) = NDY_WIND+I+1+J*NDY_WIND 
!              INE_WIND(2,II) = I+1+J*NDY_WIND 
!              INE_WIND(3,II) = I+J*NDY_WIND
!            ELSE
              II = II + 1
              INE_WIND(1,II) = I+J*NDY_WIND !I+J*NDY_WIND
              INE_WIND(2,II) = NDY_WIND+I+J* NDY_WIND !I+1+J*NDY_WIND !NDY_WIND+I+J* NDY_WIND
              INE_WIND(3,II) = I+1+J*NDY_WIND 
              !WRITE(*,*) INE_WIND(:,II)
              II = II + 1
              INE_WIND(1,II) = NDY_WIND+I+1+J*NDY_WIND 
              INE_WIND(2,II) = I+1+J*NDY_WIND 
              INE_WIND(3,II) = NDY_WIND+I+J* NDY_WIND 
              !WRITE(*,*) INE_WIND(:,II)
!            END IF
            END IF
            III = III + 1
            IF (LWRITE_ORIG_WIND) THEN
              WRITE (3010, '(5I10)') INE_WIND(:,II-1), 0, 0
              WRITE (3010, '(5I10)') INE_WIND(:,II  ), 0, 0
            ENDIF
          END DO
        END DO

        IF (LWRITE_ORIG_WIND) OPEN(3011, FILE  = 'ergwindorig.bin', FORM = 'UNFORMATTED')

        ALLOCATE(WIND_ELE(MNP)); WIND_ELE = 0
        !WRITE(*,*) MINVAL(XYPWIND(:,1)), MAXVAL(XYPWIND(:,1)), MINVAL(XP), MAXVAL(XP)
        !WRITE(*,*) MINVAL(XYPWIND(:,2)), MAXVAL(XYPWIND(:,2)), MINVAL(YP), MAXVAL(YP)
        INE_WIND = INE_WIND + 1

        DO IP = 1, MNP 
          CALL FIND_ELE_WIND( MNE_WIND, MNP_WIND, INE_WIND, XYPWIND, XP(IP), YP(IP), WIND_ELE(IP)) 
        END DO

     END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_NETCDF_NARR(IT)
         USE DATAPOOL, ONLY : WIND_NCID, WIND_X, WIND_Y, ATMO_PRESS, LINVERTY, NDX_WIND, NDY_WIND, LWRITE_ORIG_WIND
         USE NETCDF
         IMPLICIT NONE
!
!        READS WIND_Y, WIND_X and PRESSURE from a given NCID within one NARR file 
!
         INTEGER, INTENT(IN) :: IT

         INTEGER             :: DWIND_X_ID, DWIND_Y_ID, DPRESS_ID, ISTAT
         INTEGER             :: numLons, numLats, numTime, numHeights, iy, counter, ip, i, j, ix
         REAL,   ALLOCATABLE :: TMP(:,:)
         REAL, ALLOCATABLE   :: U(:), V(:), H(:)
         REAL, SAVE          :: TIME
         CHARACTER (LEN=100)  :: CHRTMP

         INTEGER, DIMENSION (nf90_max_var_dims) :: dimIDs

         ISTAT = nf90_inq_varid(WIND_NCID, '10u', DWIND_X_ID)
         ISTAT = NF90_INQUIRE_VARIABLE(WIND_NCID, DWIND_X_ID, dimids = dimIDs)
         ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(1), len = numLons)
         ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(2), len = numLats)
         ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(3), len = numHeights)
         ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(4), len = numTime)
         IF (.NOT. ALLOCATED(WIND_X)) ALLOCATE (WIND_X(NDX_WIND,NDY_WIND))

         ISTAT = nf90_inq_varid(WIND_NCID, '10v', DWIND_Y_ID)
         ISTAT = NF90_INQUIRE_VARIABLE(WIND_NCID, DWIND_Y_ID, dimids = dimIDs)
         ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(1), len = numLons)
         ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(2), len = numLats)
         ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(3), len = numHeights)
         ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(4), len = numTime)
         IF (.NOT. ALLOCATED(WIND_Y)) ALLOCATE (WIND_Y(NDX_WIND,NDY_WIND))

!        ISTAT = nf90_inq_varid(WIND_NCID, 'PS_MSL_GDS0_MSL_13', DPRESS_ID)
!        ISTAT = NF90_INQUIRE_VARIABLE(WIND_NCID, DPRESS_ID, dimids = dimIDs)
!        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(1), len = numLons)
!        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(2), len = numLats)
!        ISTAT = nf90_inquire_dimension(WIND_NCID, dimIDs(3), len = numTime)
!        IF (.NOT. ALLOCATED(ATMO_PRESS)) ALLOCATE (ATMO_PRESS(NDX_WIND,NDY_WIND))

 

         ISTAT = NF90_GET_VAR(WIND_NCID, DWIND_X_ID, WIND_X, start = (/ 1, 1, 1, IT /), count = (/ NDX_WIND, NDY_WIND,1,1 /))
         IF (ISTAT .NE. nf90_noerr) STOP 'ERR READING D_WIND_X_ID'
         ISTAT = NF90_GET_VAR(WIND_NCID, DWIND_Y_ID, WIND_Y, start = (/ 1, 1, 1, IT /), count = (/ NDX_WIND, NDY_WIND,1,1 /))
         IF (ISTAT .NE. nf90_noerr) STOP 'ERR READING D_WIND_Y_ID'

!        ISTAT = NF90_GET_VAR(WIND_NCID, DPRESS_ID, ATMO_PRESS, start = (/ 1, 1, IT /), count = (/ NDX_WIND, NDY_WIND, 1 /))
!        IF (ISTAT .NE. nf90_noerr) STOP 'ERR READING D_PRESS_ID'

         IF (.TRUE.) THEN
           ALLOCATE(TMP(NDX_WIND,NDY_WIND))
           DO IY = 1, NDY_WIND 
             tmp(:,NDY_WIND-(IY-1)) = wind_x(:,IY)
           END DO
           wind_x = tmp 
           DO IY = 1, NDY_WIND
             tmp(:,NDY_WIND-(IY-1)) = wind_y(:,IY)
           END DO
           wind_y = tmp
!           DO IY = 1, NDY_WIND 
!             tmp(:,NDY_WIND-(IY-1)) = atmo_press(:,IY)
!           END DO
!           atmo_press = tmp
           DEALLOCATE(TMP)
         END IF

         IF (.TRUE.) THEN
           ALLOCATE(TMP(NDX_WIND,NDY_WIND))
           DO IX = 1, NDX_WIND
             IF (IX .GT. NDX_WIND/2) THEN
               tmp(IX-NDX_WIND/2,:) = wind_x(IX,:)
             ELSE
               tmp(IX+NDX_WIND/2,:) = wind_x(IX,:)
             END IF
           END DO
           wind_x = tmp
           DO IX = 1, NDX_WIND
             IF (IX .GT. NDX_WIND/2) THEN
               tmp(IX-NDX_WIND/2,:) = wind_y(IX,:)
             ELSE
               tmp(IX+NDX_WIND/2,:) = wind_y(IX,:)
             END IF
           END DO
           wind_y = tmp
           DEALLOCATE(TMP)
         END IF

!         wind_y = - wind_y 
!         wind_x = - wind_x

        IF (LWRITE_ORIG_WIND) THEN
         IF (.NOT. ALLOCATED(U)) ALLOCATE(U(NDX_WIND*NDY_WIND))
         IF (.NOT. ALLOCATED(V)) ALLOCATE(V(NDX_WIND*NDY_WIND))
!         IF (.NOT. ALLOCATED(H)) ALLOCATE(H(NDX_WIND*NDY_WIND))
         COUNTER = 1 
         DO J = 1, NDY_WIND
           DO I = 1, NDX_WIND
             U(COUNTER) = WIND_X(I,J)
             V(COUNTER) = WIND_Y(I,J)
             IF (ABS(U(COUNTER)) .GT. 1000.) U(COUNTER) = 0.
             IF (ABS(V(COUNTER)) .GT. 1000.) V(COUNTER) = 0.
!             H(COUNTER) = SQRT((U(COUNTER)**2.+V(COUNTER)**2.))
             COUNTER = COUNTER + 1
           END DO
         END DO

          TIME = TIME + 1.
          WRITE(3011) TIME 
          WRITE(3011) (U(IP), V(IP), SQRT(U(IP)**2.+V(IP)**2.), IP = 1, NDX_WIND*NDY_WIND)
          DEALLOCATE(U,V)
        END IF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_NETCDF_NARR2(IT)
         USE DATAPOOL, ONLY : WINDX_NCID, WINDY_NCID, WIND_X, WIND_Y, &
     &                        ATMO_PRESS, LINVERTY, NDX_WIND, NDY_WIND, &
     &                        LWRITE_ORIG_WIND, UWND_NARR, VWND_NARR, &
     &                        NETCDF_FILE_NAMES
         USE NETCDF
         IMPLICIT NONE
!
!        READS WIND_Y, WIND_X and PRESSURE from a given NCID within one NARR file 
!
         INTEGER, INTENT(IN) :: IT

         INTEGER             :: DWIND_X_ID, DWIND_Y_ID, DPRESS_ID, ISTAT
         INTEGER             :: numLons, numLats, numTime, numHeights, iy, counter, ip, i, j, ix
         INTEGER*2, ALLOCATABLE  :: WIND_X4(:,:), WIND_Y4(:,:)
         REAL,   ALLOCATABLE :: TMP(:,:)
         REAL, ALLOCATABLE   :: U(:), V(:), H(:)
         REAL,SAVE           :: TIME
         REAL                :: scale_factor
         CHARACTER (LEN=100)  :: CHRTMP

         INTEGER, DIMENSION (nf90_max_var_dims) :: dimIDs

         ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(1), NF90_NOWRITE, WINDY_NCID) 
         ISTAT = nf90_inq_varid(WINDX_NCID, 'uwnd', DWIND_X_ID)
         ISTAT = NF90_INQUIRE_VARIABLE(WINDX_NCID, DWIND_X_ID, dimids = dimIDs)
         ISTAT = nf90_inquire_dimension(WINDX_NCID, dimIDs(1), len = numLons)
         ISTAT = nf90_inquire_dimension(WINDX_NCID, dimIDs(2), len = numLats)
         ISTAT = nf90_inquire_dimension(WINDX_NCID, dimIDs(3), len = numTime)
         ISTAT = nf90_get_att(WINDY_NCID, DWIND_X_ID, 'scale_factor', scale_factor)

         !WRITE(*,*) scale_factor
         !WRITE(*,*) numLons, numLats, numTime
         !WRITE(*,*) NDX_WIND, NDY_WIND

         IF (.NOT. ALLOCATED(WIND_X4)) ALLOCATE (WIND_X4(NDX_WIND,NDY_WIND))

         ISTAT = NF90_OPEN(NETCDF_FILE_NAMES(2), NF90_NOWRITE, WINDY_NCID)
         ISTAT = nf90_inq_varid(WINDY_NCID, 'vwnd', DWIND_Y_ID)
         ISTAT = NF90_INQUIRE_VARIABLE(WINDY_NCID, DWIND_Y_ID, dimids = dimIDs)
         ISTAT = nf90_inquire_dimension(WINDY_NCID, dimIDs(1), len = numLons)
         ISTAT = nf90_inquire_dimension(WINDY_NCID, dimIDs(2), len = numLats)
         ISTAT = nf90_inquire_dimension(WINDY_NCID, dimIDs(3), len = numTime)

         IF (.NOT. ALLOCATED(WIND_Y4)) ALLOCATE (WIND_Y4(NDX_WIND,NDY_WIND))

         ISTAT = NF90_GET_VAR(WINDX_NCID, DWIND_X_ID, WIND_X4, start = (/ 1, 1, IT /), count = (/ NDX_WIND, NDY_WIND, 1 /))
         IF (ISTAT .NE. nf90_noerr) THEN
           CHRTMP = nf90_strerror(ISTAT)
           WRITE(*,*) 'NETCDF ERROR X', CHRTMP
           STOP 'NETCDF NARR'
         ENDIF

         ISTAT = NF90_GET_VAR(WINDY_NCID, DWIND_Y_ID, WIND_Y4, start = (/ 1, 1, IT /), count = (/ NDX_WIND, NDY_WIND, 1 /))
         IF (ISTAT .NE. nf90_noerr) THEN
           CHRTMP = nf90_strerror(ISTAT)
           WRITE(*,*) 'NETCDF ERROR Y', CHRTMP, WINDX_NCID, WINDY_NCID
           STOP 'NETCDF NARR'
         ENDIF

!        ISTAT = NF90_GET_VAR(WIND_NCID, DPRESS_ID, ATMO_PRESS, start = (/ 1, 1, IT /), count = (/ NDX_WIND, NDY_WIND, 1 /))
!        IF (ISTAT .NE. nf90_noerr) STOP 'ERR READING D_PRESS_ID'

         IF (.FALSE.) THEN
           ALLOCATE(TMP(NDX_WIND,NDY_WIND))
           DO IX = 1, NDX_WIND
             tmp(NDX_WIND-(IX-1),:) = wind_x4(IX,:)
           END DO
           wind_x4 = tmp
           DO IX = 1, NDX_WIND
             tmp(NDX_WIND-(IX-1),:) = wind_y4(IX,:)
           END DO
           wind_y4 = tmp
!           DO IY = 1, NDY_WIND 
!             tmp(:,NDY_WIND-(IY-1)) = atmo_press(:,IY)
!           END DO
!           atmo_press = tmp
           DEALLOCATE(TMP)
         END IF

         IF (.FALSE.) THEN
           ALLOCATE(TMP(NDX_WIND,NDY_WIND))
           DO IX = 1, NDX_WIND
             IF (IX .GT. NDX_WIND/2) THEN
               tmp(IX-NDX_WIND/2,:) = wind_x(IX,:)
             ELSE
               tmp(IX+NDX_WIND/2,:) = wind_x(IX,:)
             END IF
           END DO
           wind_x = tmp
           DO IX = 1, NDX_WIND
             IF (IX .GT. NDX_WIND/2) THEN
               tmp(IX-NDX_WIND/2,:) = wind_y(IX,:)
             ELSE
               tmp(IX+NDX_WIND/2,:) = wind_y(IX,:)
             END IF
           END DO
           wind_y = tmp
           DEALLOCATE(TMP)
         END IF
 
         COUNTER = 1
         DO I = 1, NDX_WIND
           DO J = 1, NDY_WIND
             IF (WIND_X4(I,J) .GT. 32766 .OR. WIND_Y4(I,J) .GT. 32766) THEN
               UWND_NARR(COUNTER) = 0.
               VWND_NARR(COUNTER) = 0. 
               COUNTER = COUNTER + 1
             ELSE
               UWND_NARR(COUNTER) = WIND_X4(I,J) * scale_factor
               VWND_NARR(COUNTER) = WIND_Y4(I,J) * scale_factor 
               COUNTER = COUNTER + 1
             END IF
           END DO
         END DO

         IF (LWRITE_ORIG_WIND) THEN
           TIME = TIME + 1.
           WRITE(3011) TIME 
           WRITE(3011) (UWND_NARR(IP), VWND_NARR(IP), SQRT(UWND_NARR(IP)**2.+VWND_NARR(IP)**2.), IP = 1, NDX_WIND*NDY_WIND)
         END IF

      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_NETCDF_WW3
         USE DATAPOOL
         USE NETCDF
         IMPLICIT NONE

        INTEGER :: ISTAT, IT, IX, IY, IFILE, IVAR, BND_NCID
        INTEGER :: ILON_ID, ILAT_ID, ITIME_ID, I, J, COUNTER
        INTEGER :: D_WIND_X_ID, D_WIND_Y_ID, D_PRESS_ID
        REAL    :: RTMP
        REAL*8  :: DTMP, DTMP1, DTMP2
        REAL*8, ALLOCATABLE :: BND_TIME(:)
        character ( len = 40 ) chrtmp
        character ( len = 15 ) chrdate
        character ( len = 40 ) netcfd_fname
        character ( len = 20 ) dirname

        integer, dimension(nf90_max_var_dims) :: dimIDs

        OPEN(WAV%FHNDL,FILE=WAV%FNAME,STATUS='OLD')
!
! count number of netcdf files in list ...
!
        !WRITE(*,*) WAV%FHNDL, WAV%FNAME, BND%FHNDL, BND%FNAME

        NUM_NETCDF_FILES_BND = 0
        DO
          READ( WAV%FHNDL, *, IOSTAT = ISTAT )
          IF ( ISTAT /= 0 ) EXIT
          NUM_NETCDF_FILES_BND = NUM_NETCDF_FILES_BND + 1
        END DO
        REWIND (WAV%FHNDL)

        !WRITE(*,*) NUM_NETCDF_FILES_BND

        NUM_NETCDF_FILES_BND = NUM_NETCDF_FILES_BND / NUM_NETCDF_VAR_TYPES

        !WRITE(*,*) 'NUM_NETCDF_FILES_BND', NUM_NETCDF_FILES_BND

        ALLOCATE(NETCDF_FILE_NAMES_BND(NUM_NETCDF_FILES_BND,NUM_NETCDF_VAR_TYPES))
 
        DO IT = 1, NUM_NETCDF_FILES_BND
          DO IVAR = 1, NUM_NETCDF_VAR_TYPES
            READ( WAV%FHNDL, *) NETCDF_FILE_NAMES_BND(IT,IVAR)
          END DO
        END DO
        CLOSE (WAV%FHNDL)
!
! three files are read to set up the wave spectra Hs, Tm01, Dir
!
        ALLOCATE(NDT_BND_FILE(NUM_NETCDF_FILES_BND))

        !DO IFILE = 1, NUM_NETCDF_FILES_BND
        !  WRITE(*,'(I10,10X,5A30)') IFILE, NETCDF_FILE_NAMES_BND(IFILE,:)
        !END DO
!
! check number of time steps in netcdf file ... it is assumed that all files have the same ammount of time steps ...
!
        DO IFILE = 1, NUM_NETCDF_FILES_BND
          !WRITE(*,*) ifile, NETCDF_FILE_NAMES_BND(IFILE,1)
          ISTAT = NF90_OPEN(NETCDF_FILE_NAMES_BND(IFILE,1), NF90_NOWRITE, BND_NCID)
          ISTAT = nf90_inq_varid(BND_NCID, 'time', ITIME_ID)
          ISTAT = NF90_INQUIRE_VARIABLE(BND_NCID, ITIME_ID, dimids = dimids)
          ISTAT = nf90_inquire_dimension(BND_NCID, dimIDs(1), len = NDT_BND_FILE(IFILE))
        END DO
!
! check dimensions in the netcdf ... again it is assumed that this is not changing for all files ...
!
        ISTAT = nf90_inq_varid(BND_NCID, 'longitude', ILON_ID)
        ISTAT = NF90_INQUIRE_VARIABLE(BND_NCID, ILON_ID, dimids = dimIDs)
        ISTAT = nf90_inquire_dimension(BND_NCID, dimIDs(1), len = NDX_BND)

        ISTAT = nf90_inq_varid(BND_NCID, 'latitude', ILAT_ID)
        ISTAT = NF90_INQUIRE_VARIABLE(BND_NCID, ILAT_ID, dimids = dimIDs)
        ISTAT = nf90_inquire_dimension(BND_NCID, dimIDs(1), len = NDY_BND)

        ALLOCATE (COORD_BND_X(NDX_BND))
        ALLOCATE (COORD_BND_Y(NDY_BND))
!
! read cooridantes from files ....
!
        ISTAT = NF90_GET_VAR(BND_NCID, ILON_ID, COORD_BND_X)
        ISTAT = NF90_GET_VAR(BND_NCID, ILAT_ID, COORD_BND_Y)
!
! estimate offset ...
!
        OFFSET_X_BND = MINVAL(COORD_BND_X)
        OFFSET_Y_BND = MINVAL(COORD_BND_Y)
!
! resolution ...
!
        DX_BND  = ABS(MAXVAL(COORD_BND_X)-MINVAL(COORD_BND_X))/(NDX_BND-1)
        DY_BND  = ABS(MAXVAL(COORD_BND_Y)-MINVAL(COORD_BND_Y))/(NDY_BND-1)
!
! close netcdf file ...
!
        ISTAT = NF90_CLOSE(BND_NCID)
!
! total number of time steps ... in all files 
!
        NDT_BND_ALL_FILES = 0
        DO IT = 1, NUM_NETCDF_FILES_BND
          NDT_BND_ALL_FILES = NDT_BND_ALL_FILES + NDT_BND_FILE(IT) 
        END DO

        !WRITE(*,*) NDT_BND_ALL_FILES, NDT_BND_FILE

        ALLOCATE (BND_TIME_ALL_FILES(NUM_NETCDF_FILES_BND,MAXVAL(NDT_BND_FILE))); BND_TIME_ALL_FILES = 0.d0
!
! read all time steps in the proper format and transform in wwm time line 
!
        BND_TIME_ALL_FILES = 0.
        DO IFILE = 1, NUM_NETCDF_FILES_BND
          ISTAT = NF90_OPEN(NETCDF_FILE_NAMES_BND(IFILE,1),NF90_NOWRITE,BND_NCID)
          ALLOCATE (BND_TIME(NDT_BND_FILE(IFILE))); BND_TIME = 0.d0
          ISTAT = NF90_GET_VAR(BND_NCID,ITIME_ID,BND_TIME)
          DO IT = 1, NDT_BND_FILE(IFILE)
             BND_TIME_ALL_FILES(IFILE,IT) = BND_TIME(IT)
!             CALL CT2MJD('19000101.000000',DTMP1)
!             CALL CT2MJD('19900101.000000',DTMP2)
!             CALL MJD2CT(DTMP1,chrdate)
!             WRITE(*,*) '19000101.000000', DTMP1, chrdate
!             CALL MJD2CT(DTMP2,chrdate)
!             WRITE(*,*) '19900101.000000', DTMP2, chrdate
!             CALL MJD2CT(0.d0,chrdate)
!             WRITE(*,*) '00000000.000000', 0.d0, chrdate
!             WRITE(*,*) BND_TIME_ALL_FILES(1,1), DT_DIFF_19901900
!             IF (IT == 1 .AND. IFILE ==1) WRITE(*,*) DTMP1, DTMP2, DTMP1+DT_DIFF_19901900
!             IF (IT == 1 .AND. IFILE ==1) WRITE(*,*) IFILE, IT, BND_TIME(IT), chrdate
          END DO
          DEALLOCATE(BND_TIME)
        END DO ! IFILE

        SEBO%DELT = (BND_TIME_ALL_FILES(1,2) - BND_TIME_ALL_FILES(1,1)) * DAY2SEC

        BND_TIME_ALL_FILES = BND_TIME_ALL_FILES + DT_DIFF_19901900 

        IF (LWRITE_ALL_WW3_RESULTS) THEN
          OPEN(3010, FILE  = 'sysglobalboundary.dat', STATUS = 'UNKNOWN')
          WRITE (3010, '(I10)') 0
          WRITE (3010, '(I10)') NDX_BND * NDY_BND
          COUNTER = 0
          DO I = 1, NDY_BND
            DO J = 1, NDX_BND
              WRITE (3010, '(I10,3F15.4)') COUNTER, OFFSET_X_BND+(J-1)*DX_BND,OFFSET_Y_BND+(I-1)*DY_BND, 0.0
              COUNTER = COUNTER + 1
            END DO
          END DO
          WRITE (3010, *) (NDX_BND-1)*(NDY_BND-1)*2
          DO J = 0, NDY_BND-2
            DO I = 0, NDX_BND-2
              WRITE (3010, '(5I10)')  I+J*NDX_BND           , NDX_BND+I+J* NDX_BND, NDX_BND+I+1+J*NDX_BND, 0, 0
              WRITE (3010, '(5I10)')  NDX_BND+I+1+J*NDX_BND, I+1+J*NDX_BND        , I+J*NDX_BND          , 0, 0
            END DO
          END DO
          CLOSE(3010)
        END IF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_NETCDF_WW3(IFILE,IT)

         USE DATAPOOL
         USE NETCDF
         IMPLICIT NONE
!
!        READS WIND_Y, WIND_X and PRESSURE from a given NCID within one DWD file 
!
         INTEGER, INTENT(IN) :: IFILE, IT

         INTEGER              :: HS_WW3_ID, T02_WW3_ID, DIR_WW3_ID, FP_WW3_ID, DSPR_WW3_ID
         INTEGER              :: HS_BND_NCID, T02_BND_NCID, DIR_BND_NCID, FP_BND_NCID, DSPR_BND_NCID
         INTEGER              :: ISTAT
         INTEGER              :: numLons, numLats, numTime, iy, counter, ip, i, j
         INTEGER, ALLOCATABLE :: ITMP(:,:)
         REAL,   ALLOCATABLE  :: TMP(:,:)
         REAL, ALLOCATABLE    :: U(:), V(:), H(:)
         REAL, SAVE           :: TIME, scale_factor

         INTEGER, DIMENSION (nf90_max_var_dims) :: dimIDs
         CHARACTER(LEN=80)    :: CHRTMP

         ALLOCATE (ITMP(NDX_BND,NDY_BND))

         WRITE(DBG%FHNDL,*) IT, IFILE, 'READING GLOBAL DATA'
         ISTAT = NF90_OPEN(NETCDF_FILE_NAMES_BND(IFILE,3),NF90_NOWRITE,HS_BND_NCID)
         IF (ISTAT .NE. nf90_noerr) THEN
           CHRTMP = nf90_strerror(ISTAT)
           WRITE(DBG%FHNDL,*) 'NETCDF ERROR -1-', CHRTMP
           STOP 'NETCDF WW3'
         ENDIF
         ISTAT = nf90_inq_varid(HS_BND_NCID, TRIM(NCDF_HS_NAME), HS_WW3_ID) 
         IF (ISTAT .NE. nf90_noerr) THEN
           CHRTMP = nf90_strerror(ISTAT)
           WRITE(DBG%FHNDL,*) 'NETCDF ERROR -2-', CHRTMP
           STOP 'NETCDF WW3'
         END IF
         ISTAT = nf90_get_att(HS_BND_NCID, HS_WW3_ID, 'scale_factor', scale_factor)
         IF (ISTAT .NE. nf90_noerr) THEN
           CHRTMP = nf90_strerror(ISTAT)
           WRITE(DBG%FHNDL,*) 'NETCDF ERROR -3-', CHRTMP
           STOP 'NETCDF WW3'
         ENDIF
         IF (.NOT. ALLOCATED(HS_WW3)) ALLOCATE (HS_WW3(NDX_BND,NDY_BND)); HS_WW3 = 0.
         ISTAT = NF90_GET_VAR(HS_BND_NCID, HS_WW3_ID, ITMP,  start = (/ 1, 1, IT /), count = (/ NDX_BND, NDY_BND, 1/))
         IF (ISTAT .NE. nf90_noerr) THEN
           CHRTMP = nf90_strerror(ISTAT)
           WRITE(DBG%FHNDL,*) 'NETCDF ERROR -4-', CHRTMP
           STOP 'NETCDF WW3'
         ENDIF
         HS_WW3 = REAL(ITMP) * scale_factor
         ISTAT = nf90_close(HS_BND_NCID)
         IF (ISTAT .NE. nf90_noerr) THEN
           CHRTMP = nf90_strerror(ISTAT)
           WRITE(DBG%FHNDL,*) 'NETCDF ERROR -5-', CHRTMP
           STOP 'NETCDF WW3'
         ENDIF

         ISTAT = NF90_OPEN(NETCDF_FILE_NAMES_BND(IFILE,2),NF90_NOWRITE,FP_BND_NCID)
         IF (ISTAT .NE. nf90_noerr) WRITE(*,*) 'Erorr Reading fp - 1' 
         ISTAT = nf90_inq_varid(FP_BND_NCID, TRIM(NCDF_FP_NAME), FP_WW3_ID)
         IF (ISTAT .NE. nf90_noerr) WRITE(*,*) 'Erorr Reading fp - 2'
         ISTAT = nf90_get_att(FP_BND_NCID, FP_WW3_ID, 'scale_factor', scale_factor)
         IF (ISTAT .NE. nf90_noerr) WRITE(*,*) 'Erorr Reading fp - 3'
         IF (.NOT. ALLOCATED(FP_WW3)) ALLOCATE (FP_WW3(NDX_BND,NDY_BND)); FP_WW3 = 0.
         ISTAT = NF90_GET_VAR(FP_BND_NCID, FP_WW3_ID, ITMP,  start = (/ 1, 1, IT /), count = (/ NDX_BND, NDY_BND, 1/))
         IF (ISTAT .NE. nf90_noerr) WRITE(*,*) 'Erorr Reading fp - 4'
         FP_WW3 = REAL(ITMP) * scale_factor
         ISTAT = nf90_close(FP_BND_NCID)

         ISTAT = NF90_OPEN(NETCDF_FILE_NAMES_BND(IFILE,5),NF90_NOWRITE,T02_BND_NCID)
         IF (ISTAT .NE. nf90_noerr) WRITE(*,*) 'Erorr Reading tm02 - 1'
         ISTAT = nf90_inq_varid(T02_BND_NCID, TRIM(NCDF_F02_NAME), T02_WW3_ID)
         IF (ISTAT .NE. nf90_noerr) WRITE(*,*) 'Erorr Reading tm02 - 2'
         ISTAT = nf90_get_att(T02_BND_NCID, T02_WW3_ID, 'scale_factor', scale_factor)
         IF (ISTAT .NE. nf90_noerr) WRITE(*,*) 'Erorr Reading tm02 - 3'
         IF (.NOT. ALLOCATED(T02_WW3)) ALLOCATE (T02_WW3(NDX_BND,NDY_BND)); T02_WW3 = 0.
         ISTAT = NF90_GET_VAR(T02_BND_NCID, T02_WW3_ID, ITMP,  start = (/ 1, 1, IT /), count = (/ NDX_BND, NDY_BND, 1/))
         IF (ISTAT .NE. nf90_noerr) WRITE(*,*) 'Erorr Reading tm02 - 4'
         T02_WW3 = REAL(ITMP) * scale_factor
         ISTAT = nf90_close(T02_BND_NCID)

         ISTAT = NF90_OPEN(NETCDF_FILE_NAMES_BND(IFILE,4),NF90_NOWRITE,DSPR_BND_NCID)
         IF (ISTAT .NE. nf90_noerr) WRITE(*,*) 'Erorr Reading spr - 1'
         ISTAT = nf90_inq_varid(DSPR_BND_NCID, TRIM(NCDF_SPR_NAME), DSPR_WW3_ID)
         IF (ISTAT .NE. nf90_noerr) WRITE(*,*) 'Erorr Reading spr - 2'
         ISTAT = nf90_get_att(DSPR_BND_NCID, DSPR_WW3_ID, 'scale_factor', scale_factor)
         IF (ISTAT .NE. nf90_noerr) WRITE(*,*) 'Erorr Reading spr - 1'
         IF (.NOT. ALLOCATED(DSPR_WW3)) ALLOCATE (DSPR_WW3(NDX_BND,NDY_BND)); DSPR_WW3 = 0.
         ISTAT = NF90_GET_VAR(DSPR_BND_NCID, DSPR_WW3_ID, ITMP,  start = (/ 1, 1, IT /), count = (/ NDX_BND, NDY_BND, 1/))
         IF (ISTAT .NE. nf90_noerr) WRITE(*,*) 'Erorr Reading spr - 1'
         DSPR_WW3 = REAL(ITMP) * scale_factor
         ISTAT = nf90_close(DSPR_BND_NCID)

         ISTAT = NF90_OPEN(NETCDF_FILE_NAMES_BND(IFILE,1),NF90_NOWRITE,DIR_BND_NCID)
         IF (ISTAT .NE. nf90_noerr) WRITE(*,*) 'Erorr Reading dir - 1'
         ISTAT = nf90_inq_varid(DIR_BND_NCID, TRIM(NCDF_DIR_NAME), DIR_WW3_ID)
         IF (ISTAT .NE. nf90_noerr) WRITE(*,*) 'Erorr Reading dir - 2'
         ISTAT = nf90_get_att(DIR_BND_NCID, DIR_WW3_ID, 'scale_factor', scale_factor)
         IF (ISTAT .NE. nf90_noerr) WRITE(*,*) 'Erorr Reading dir - 3'
         IF (.NOT. ALLOCATED(DIR_WW3)) ALLOCATE (DIR_WW3(NDX_BND,NDY_BND)); DIR_WW3 = 0.
         ISTAT = NF90_GET_VAR(DIR_BND_NCID, DIR_WW3_ID, ITMP,  start = (/ 1, 1, IT /), count = (/ NDX_BND, NDY_BND, 1/))
         IF (ISTAT .NE. nf90_noerr) WRITE(*,*) 'Erorr Reading dir - 4'
         DIR_WW3 = REAL(ITMP) * scale_factor
         ISTAT = nf90_close(DIR_BND_NCID)

         IF (LWRITE_WW3_RESULTS) THEN
           OPEN(3012, FILE  = 'ergwiii.bin', FORM = 'UNFORMATTED')
           IF (.NOT. ALLOCATED(U)) ALLOCATE(U(NDX_BND*NDY_BND))
           IF (.NOT. ALLOCATED(V)) ALLOCATE(V(NDX_BND*NDY_BND))
           IF (.NOT. ALLOCATED(H)) ALLOCATE(H(NDX_BND*NDY_BND))
           COUNTER = 1
           DO J = 1, NDY_BND
             DO I = 1, NDX_BND
               U(COUNTER) = REAL( HS_WW3(I,J)   )
               V(COUNTER) = REAL( DIR_WW3(I,J)  )
               H(COUNTER) = REAL( DSPR_WW3(I,J) )
               COUNTER = COUNTER + 1
             END DO
           END DO
           TIME = TIME + 1.
           WRITE(3012) TIME 
           WRITE(3012) (U(IP), V(IP), H(IP), IP = 1, NDX_BND*NDY_BND)
           DEALLOCATE(U,V,H)
         END IF

      END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE INTER_STRUCT_DATA(NDX,NDY,DX,DY,OFFSET_X,OFFSET_Y,MAT,VAL)
          USE DATAPOOL
          IMPLICIT NONE

          INTEGER, INTENT(IN) :: NDX, NDY
          REAL, INTENT(IN)    :: MAT(NDX,NDY)
          REAL, INTENT(IN)    :: DX, DY, OFFSET_X, OFFSET_Y
          REAL, INTENT(OUT)   :: VAL(MNP)

          INTEGER             :: IP, J_INT, I_INT, INTER_X, INTER_Y
          REAL                :: WX1, WX2, WX3, WX4, HX1, HX2
          REAL                :: DELTA_X, DELTA_Y, LEN_X, LEN_Y

          DO IP = 1, MNP
            LEN_X = XP(IP)-OFFSET_X
            LEN_Y = YP(IP)-OFFSET_Y
            I_INT = INT( LEN_X/DX ) + 1
            J_INT = INT( LEN_Y/DY ) + 1
            DELTA_X = LEN_X - (I_INT - 1) * DX ! Abstand X u. Y
            DELTA_Y = LEN_Y - (J_INT - 1) * DY ! 
            WX1     = MAT(  I_INT   , J_INT  ) ! Unten Links
            WX2     = MAT(  I_INT   , J_INT+1) ! Oben  Links
            WX3     = MAT(  I_INT+1,  J_INT+1) ! Oben  Rechts
            WX4     = MAT(  I_INT+1,  J_INT  ) ! Unten Rechts
            IF (IWINDFORMAT == 2) THEN
              HX1 = WX1 + (WX2-WX1)/DX * DELTA_X
              HX2 = WX4 + (WX3-WX4)/DX * DELTA_X
            ELSE IF (IWINDFORMAT == 3) THEN
              HX1 = WX1 + (WX4-WX1)/DX * DELTA_X
              HX2 = WX2 + (WX3-WX2)/DX * DELTA_X
            END IF
            VAL(IP) = HX1 + (HX2-HX1)/DY * DELTA_Y
         END DO
      END SUBROUTINE        
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE INTER_STRUCT_BOUNDARY(NDX,NDY,DX,DY,OFFSET_X,OFFSET_Y,VAL)
          USE DATAPOOL
          IMPLICIT NONE

          INTEGER, INTENT(IN) :: NDX, NDY
          REAL, INTENT(IN)    :: DX, DY, OFFSET_X, OFFSET_Y
          REAL, INTENT(OUT)   :: VAL(8,IWBMNP)

          INTEGER             :: IP, J_INT, I_INT, INTER_X, INTER_Y
          REAL                :: WX1, WX2, WX3, WX4, HX1, HX2
          REAL                :: DELTA_X, DELTA_Y, LEN_X, LEN_Y

         DO IP = 1, IWBMNP 

            LEN_X = XP(IWBNDLC(IP)) - OFFSET_X
            LEN_Y = YP(IWBNDLC(IP)) - OFFSET_Y
    
            I_INT = INT( LEN_X/DX ) + 1
            J_INT = INT( LEN_Y/DY ) + 1

            DELTA_X   = LEN_X - (I_INT - 1) * DX ! Abstand X u. Y
            DELTA_Y   = LEN_Y - (J_INT - 1) * DY ! 

            !WRITE(*,*) 'XP YP', XP(IWBNDLC(IP)), YP(IWBNDLC(IP)), IWBNDLC(IP)
            !WRITE(*,*) LEN_X, LEN_Y, OFFSET_X, OFFSET_Y, XP(IWBNDLC(IP)), YP(IWBNDLC(IP))

            WX1       = HS_WW3(  I_INT   , J_INT  ) ! Unten Links
            WX2       = HS_WW3(  I_INT   , J_INT+1) ! Oben  Links
            WX3       = HS_WW3(  I_INT+1,  J_INT+1) ! Oben  Rechts
            WX4       = HS_WW3(  I_INT+1,  J_INT  ) ! Unten Rechts

            HX1       = WX1 + (WX4-WX1)/DX * DELTA_X
            HX2       = WX2 + (WX3-WX2)/DX * DELTA_X

            !WRITE(*,*) WX1, WX2, WX3, WX4

            IF (WX1 .LT. 0. .OR. WX2 .LT. 0. .OR. WX3 .LT. 0. .OR. WX4 .LT. 0. ) THEN
              VAL(1,IP) = 0.
            ELSE
              VAL(1,IP) = HX1 + (HX2-HX1)/DY * DELTA_Y
            ENDIF
 
            WX1       = DIR_WW3(  I_INT   , J_INT  ) ! Unten Links
            WX2       = DIR_WW3(  I_INT   , J_INT+1) ! Oben  Links
            WX3       = DIR_WW3(  I_INT+1,  J_INT+1) ! Oben  Rechts
            WX4       = DIR_WW3(  I_INT+1,  J_INT  ) ! Unten Rechts

            HX1       = WX1 + (WX4-WX1)/DX * DELTA_X
            HX2       = WX2 + (WX3-WX2)/DX * DELTA_X

            IF (WX1 .LT. 0. .OR. WX2 .LT. 0. .OR. WX3 .LT. 0. .OR. WX4 .LT. 0. ) THEN
              VAL(3,IP) = 0.
            ELSE
              VAL(3,IP) = HX1 + (HX2-HX1)/DY * DELTA_Y
            ENDIF
          
            WX1       = FP_WW3(  I_INT   , J_INT  ) ! Unten Links
            WX2       = FP_WW3(  I_INT   , J_INT+1) ! Oben  Links
            WX3       = FP_WW3(  I_INT+1,  J_INT+1) ! Oben  Rechts
            WX4       = FP_WW3(  I_INT+1,  J_INT  ) ! Unten Rechts

            HX1       = WX1 + (WX4-WX1)/DX * DELTA_X
            HX2       = WX2 + (WX3-WX2)/DX * DELTA_X

            !IF (WX1 .LT. 0. .OR. WX2 .LT. 0. .OR. WX3 .LT. 0. .OR. WX4 .LT. 0. ) THEN
            !  VAL(2,IP) = 0.
            !ELSE
            !  VAL(2,IP) = HX1 + (HX2-HX1)/DY * DELTA_Y 
            !ENDIF

            !IF (VAL(2,IP) .GT. TINY(1.)) THEN
            !  VAL(2,IP) = 1. / VAL(2,IP) 
            !ELSE
            !  VAL(2,IP) = 0.
            !END IF

            WX1       = T02_WW3(  I_INT   , J_INT  ) ! Unten Links
            WX2       = T02_WW3(  I_INT   , J_INT+1) ! Oben  Links
            WX3       = T02_WW3(  I_INT+1,  J_INT+1) ! Oben  Rechts
            WX4       = T02_WW3(  I_INT+1,  J_INT  ) ! Unten Rechts

            HX1       = WX1 + (WX4-WX1)/DX * DELTA_X
            HX2       = WX2 + (WX3-WX2)/DX * DELTA_X

            IF (WX1 .LT. 0. .OR. WX2 .LT. 0. .OR. WX3 .LT. 0. .OR. WX4 .LT. 0. ) THEN
              VAL(2,IP) = 0.
            ELSE
              VAL(2,IP) = HX1 + (HX2-HX1)/DY * DELTA_Y
            ENDIF

            WX1       = DSPR_WW3(  I_INT   , J_INT  ) ! Unten Links
            WX2       = DSPR_WW3(  I_INT   , J_INT+1) ! Oben  Links
            WX3       = DSPR_WW3(  I_INT+1,  J_INT+1) ! Oben  Rechts
            WX4       = DSPR_WW3(  I_INT+1,  J_INT  ) ! Unten Rechts

            HX1       = WX1 + (WX4-WX1)/DX * DELTA_X
            HX2       = WX2 + (WX3-WX2)/DX * DELTA_X

            IF (WX1 .LT. 0. .OR. WX2 .LT. 0. .OR. WX3 .LT. 0. .OR. WX4 .LT. 0. ) THEN
              VAL(4,IP) = 0.
            ELSE
              VAL(4,IP) = HX1 + (HX2-HX1)/DY * DELTA_Y
            ENDIF

            VAL(5,:)  = -2.
            VAL(6,:)  = 1.
            VAL(7,:)  = 0.1 
            VAL(8,:)  = 3.3 

         END DO

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE INTER_STRUCT_DOMAIN(NDX,NDY,DX,DY,OFFSET_X,OFFSET_Y,VAL)
          USE DATAPOOL
          IMPLICIT NONE

          INTEGER, INTENT(IN) :: NDX, NDY
          REAL, INTENT(IN)    :: DX, DY, OFFSET_X, OFFSET_Y
          REAL, INTENT(OUT)   :: VAL(8,MNP)

          INTEGER             :: IP, J_INT, I_INT, INTER_X, INTER_Y
          REAL                :: WX1, WX2, WX3, WX4, HX1, HX2
          REAL                :: DELTA_X, DELTA_Y, LEN_X, LEN_Y

          DO IP = 1, MNP 

            LEN_X = XP(IP) - OFFSET_X
            LEN_Y = YP(IP) - OFFSET_Y
    
            I_INT = INT( LEN_X/DX ) + 1
            J_INT = INT( LEN_Y/DY ) + 1

            DELTA_X   = LEN_X - (I_INT - 1) * DX ! Abstand X u. Y
            DELTA_Y   = LEN_Y - (J_INT - 1) * DY ! 

            WX1       = HS_WW3(  I_INT   , J_INT  ) ! Unten Links
            WX2       = HS_WW3(  I_INT   , J_INT+1) ! Oben  Links
            WX3       = HS_WW3(  I_INT+1,  J_INT+1) ! Oben  Rechts
            WX4       = HS_WW3(  I_INT+1,  J_INT  ) ! Unten Rechts

            HX1       = WX1 + (WX4-WX1)/DX * DELTA_X
            HX2       = WX2 + (WX3-WX2)/DX * DELTA_X

            IF (WX1 .LT. 0. .OR. WX2 .LT. 0. .OR. WX3 .LT. 0. .OR. WX4 .LT. 0. ) THEN
              VAL(1,IP) = 0.
            ELSE
              VAL(1,IP) = HX1 + (HX2-HX1)/DY * DELTA_Y
            ENDIF
 
            WX1       = DIR_WW3(  I_INT   , J_INT  ) ! Unten Links
            WX2       = DIR_WW3(  I_INT   , J_INT+1) ! Oben  Links
            WX3       = DIR_WW3(  I_INT+1,  J_INT+1) ! Oben  Rechts
            WX4       = DIR_WW3(  I_INT+1,  J_INT  ) ! Unten Rechts

            HX1       = WX1 + (WX4-WX1)/DX * DELTA_X
            HX2       = WX2 + (WX3-WX2)/DX * DELTA_X

            IF (WX1 .LT. 0. .OR. WX2 .LT. 0. .OR. WX3 .LT. 0. .OR. WX4 .LT. 0. ) THEN
              VAL(3,IP) = 0.
            ELSE
              VAL(3,IP) = HX1 + (HX2-HX1)/DY * DELTA_Y
            ENDIF
          
            WX1       = FP_WW3(  I_INT   , J_INT  ) ! Unten Links
            WX2       = FP_WW3(  I_INT   , J_INT+1) ! Oben  Links
            WX3       = FP_WW3(  I_INT+1,  J_INT+1) ! Oben  Rechts
            WX4       = FP_WW3(  I_INT+1,  J_INT  ) ! Unten Rechts

            HX1       = WX1 + (WX4-WX1)/DX * DELTA_X
            HX2       = WX2 + (WX3-WX2)/DX * DELTA_X

!2do: check with fabrice ... peak period looks very strange ...

            IF (WX1 .LT. 0. .OR. WX2 .LT. 0. .OR. WX3 .LT. 0. .OR. WX4 .LT. 0. ) THEN
              VAL(2,IP) = 0.
            ELSE
              VAL(2,IP) = HX1 + (HX2-HX1)/DY * DELTA_Y 
            ENDIF

!            IF (VAL(2,IP) .GT. TINY(1.)) THEN
!              VAL(2,IP) = 1. / VAL(2,IP) 
!            ELSE
!              VAL(2,IP) = 0.
!            END IF

            WX1       = T02_WW3(  I_INT   , J_INT  ) ! Unten Links
            WX2       = T02_WW3(  I_INT   , J_INT+1) ! Oben  Links
            WX3       = T02_WW3(  I_INT+1,  J_INT+1) ! Oben  Rechts
            WX4       = T02_WW3(  I_INT+1,  J_INT  ) ! Unten Rechts

            HX1       = WX1 + (WX4-WX1)/DX * DELTA_X
            HX2       = WX2 + (WX3-WX2)/DX * DELTA_X

            IF (WX1 .LT. 0. .OR. WX2 .LT. 0. .OR. WX3 .LT. 0. .OR. WX4 .LT. 0. ) THEN
              VAL(2,IP) = 0.
            ELSE
              VAL(2,IP) = HX1 + (HX2-HX1)/DY * DELTA_Y
            ENDIF

            WX1       = DSPR_WW3(  I_INT   , J_INT  ) ! Unten Links
            WX2       = DSPR_WW3(  I_INT   , J_INT+1) ! Oben  Links
            WX3       = DSPR_WW3(  I_INT+1,  J_INT+1) ! Oben  Rechts
            WX4       = DSPR_WW3(  I_INT+1,  J_INT  ) ! Unten Rechts

            HX1       = WX1 + (WX4-WX1)/DX * DELTA_X
            HX2       = WX2 + (WX3-WX2)/DX * DELTA_X

            IF (WX1 .LT. 0. .OR. WX2 .LT. 0. .OR. WX3 .LT. 0. .OR. WX4 .LT. 0. ) THEN
              VAL(4,IP) = 0.
            ELSE
              VAL(4,IP) = HX1 + (HX2-HX1)/DY * DELTA_Y
            ENDIF

            VAL(5,IP)  = -2. ! From mean period ...
            VAL(6,IP)  = 1.
            VAL(7,IP)  = 0.1 
            VAL(8,IP)  = 3.3 

            !WRITE(*,'(I10,8F15.4)') IP, VAL(:,IP)
            !PAUSE

         END DO

         IF (LWRITE_INTERPOLATED_WW3_RESULTS) THEN
           OPEN(4013, FILE  = 'erginterwiii.bin', FORM = 'UNFORMATTED')
           WRITE(4013) RTIME
           WRITE(4013) (VAL(3,IP), VAL(2,IP), VAL(4,IP), IP = 1, MNP)
         END IF 

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
