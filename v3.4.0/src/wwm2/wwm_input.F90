#include "wwm_functions.h"
!    25 Mar 2004    4:34 pm
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef MPI_PARALL_GRID
# define wwm_print_namelist(xinp) IF (myrank.eq.0) WRITE(CHK%FHNDL, NML=xinp)
#else
# define wwm_print_namelist(xinp) WRITE(CHK%FHNDL, NML=xinp)
#endif
      SUBROUTINE READ_HISTORY_STATION_NAMELIST()
#ifdef MPI_PARALL_GRID
         use elfe_msgp, only : myrank
#endif
         USE DATAPOOL, only : OUTSTYLE, LOUTITER,                       &
     &        LENERGY, LWXFN, OUT_HISTORY, OUT_STATION, INP, LSP1D,     &
     &        LSP2D, INP, CHK, LOUTS, IOUTS, LSIGMAX, LLOUTS,           &
     &        ILOUTS, OUT, DAY2SEC, FRHIGH, DBG, LINES, VAROUT_HISTORY, &
     &        VAROUT_STATION, GRIDWRITE, RKIND, LVAR_READ,              &
     &        PARAMWRITE_HIS, PARAMWRITE_STAT, wwmerr, LCFL
#ifdef NCDF
         USE NETCDF
         USE DATAPOOL, only : USE_SINGLE_OUT_STAT, USE_SINGLE_OUT_HIS,  &
     &        MULTIPLEOUT_HIS, MULTIPLEOUT_STAT,                        &
     &        NF90_OUTTYPE_STAT, NF90_OUTTYPE_HIS
#endif
         USE DATAPOOL, only : STATION_P => STATION
         USE DATAPOOL, only : MAIN, PRINTMMA, ZERO
         USE DATAPOOL, only : WriteOutputProcess_hot
         USE DATAPOOL, only : WriteOutputProcess_his
         USE DATAPOOL, only : WriteOutputProcess_stat

         IMPLICIT NONE
         CHARACTER(LEN=40)  :: FILEOUT
         INTEGER, PARAMETER :: INUMOUTS = 200 
         CHARACTER(LEN=20)  :: BEGTC, UNITC, ENDTC, NOUTS(INUMOUTS), NLOUTS(INUMOUTS)
         INTEGER            :: NLPOINTS(INUMOUTS)
         REAL(rkind)        :: XOUTS(INUMOUTS), YOUTS(INUMOUTS), CUTOFF(INUMOUTS)
         REAL(rkind)        :: XLOUTS(INUMOUTS), YLOUTS(INUMOUTS)
         REAL(rkind) :: DEFINETC
         INTEGER     :: MULTIPLEOUT
         LOGICAL     :: USE_SINGLE_OUT
         REAL(rkind) :: DELTC
         INTEGER     :: I
         integer istat
         LOGICAL     :: PARAMWRITE
         LOGICAL     :: AC, WK, ACOUT_1D, ACOUT_2D
         LOGICAL     ::   HS, TM01, TM02, TM10, KLM, WLM,               &
     &      ETOTC, ETOTS, DM, DSPR,                                     &
     &      TPPD, CPPD, KPPD, CGPD,                                     &
     &      TPP, CPP, WNPP, CGPP, KPP, LPP, PEAKD, PEAKDSPR,            &
     &      DPEAK, UBOT, ORBITAL, BOTEXPER, TMBOT,                      &
     &      URSELL, UFRIC, Z0, ALPHA_CH, WINDX, WINDY, CD,              &
     &      CURRTX, CURRTY, WATLEV, WATLEVOLD, DEPDT, DEP,              &
     &      WINDMAG, TAUW, TAUWX, TAUWY, TAUHF, TAUTOT,                 &
     &      STOKESBOTTX, STOKESBOTTY,                                   &
     &      STOKESSURFX, STOKESSURFY, STOKESBAROX, STOKESBAROY,         &
     &      RSXX, RSXY, RSYY, CFL1, CFL2, CFL3

         NAMELIST /HISTORY/ BEGTC, DELTC, UNITC, ENDTC, DEFINETC,       &
     &      OUTSTYLE, FILEOUT, LOUTITER,                                &
     &      LENERGY, LWXFN, GRIDWRITE, PARAMWRITE,                      &
     &      MULTIPLEOUT, USE_SINGLE_OUT, PRINTMMA,                      &
     &      HS, TM01, TM02, TM10, KLM, WLM,                             &
     &      ETOTC, ETOTS, DM, DSPR,                                     &
     &      TPPD, CPPD, KPPD, CGPD,                                     &
     &      TPP, CPP, WNPP, CGPP, KPP, LPP, PEAKD, PEAKDSPR,            &
     &      DPEAK, UBOT, ORBITAL, BOTEXPER, TMBOT,                      &
     &      URSELL, UFRIC, Z0, ALPHA_CH, WINDX, WINDY, CD,              &
     &      CURRTX, CURRTY, WATLEV, WATLEVOLD, DEPDT, DEP,              &
     &      WINDMAG, TAUW, TAUWX, TAUWY, TAUHF, TAUTOT,                 &
     &      STOKESBOTTX, STOKESBOTTY,                                   &
     &      STOKESSURFX, STOKESSURFY, STOKESBAROX, STOKESBAROY,         &
     &      RSXX, RSXY, RSYY, CFL1, CFL2, CFL3

         NAMELIST /STATION/ BEGTC, DELTC, UNITC, ENDTC, DEFINETC,       &
     &      OUTSTYLE, USE_SINGLE_OUT, MULTIPLEOUT, PARAMWRITE,          &
     &      FILEOUT, LOUTITER, IOUTS, NOUTS, XOUTS, YOUTS,              &
     &      CUTOFF, LSIGMAX, LSP1D, LSP2D, LLOUTS, ILOUTS, NLOUTS,      &
     &      AC, WK, ACOUT_1D, ACOUT_2D,                                 &
     &      HS, TM01, TM02, TM10, KLM, WLM,                             &
     &      ETOTC, ETOTS, DM, DSPR, TPPD, CPPD, KPPD, CGPD, TPP,        &
     &      CPP, WNPP, CGPP, KPP, LPP, PEAKD, PEAKDSPR, DPEAK,          &
     &      UBOT, ORBITAL, BOTEXPER, TMBOT,                             &
     &      URSELL, UFRIC, Z0, ALPHA_CH, WINDX, WINDY, CD,              &
     &      CURRTX, CURRTY, WATLEV, WATLEVOLD, DEPDT, DEP,              &
     &      WINDMAG, TAUW, TAUWX, TAUWY, TAUHF, TAUTOT,                 &
     &      STOKESSURFX, STOKESSURFY, STOKESBAROX, STOKESBAROY,         &
     &      RSXX, RSXY, RSYY, CFL1, CFL2, CFL3

         XOUTS = 0.
         YOUTS = 0.
         XLOUTS = 0.
         YLOUTS = 0.
         NOUTS = ''
         NLOUTS = ''
         CUTOFF = 0.

!
!     **** HISTORY section
!
#ifdef NCDF
         MULTIPLEOUT=0
         USE_SINGLE_OUT=.TRUE.
         DEFINETC=-1
#endif
         FILEOUT = "zorglub"
         HS=.FALSE.
         TM01=.FALSE.
         TM02=.FALSE.
         TM10=.FALSE.
         KLM=.FALSE.
         WLM=.FALSE.
         ETOTC=.FALSE.
         ETOTS=.FALSE.
         DM=.FALSE.
         DSPR=.FALSE.
         TPPD=.FALSE.
         CPPD=.FALSE.
         KPPD=.FALSE.
         CGPD=.FALSE.
         TPP=.FALSE.
         CPP=.FALSE.
         WNPP=.FALSE.
         CGPP=.FALSE.
         KPP=.FALSE.
         LPP=.FALSE.
         PEAKD=.FALSE.
         PEAKDSPR=.FALSE.
         DPEAK=.FALSE.
         UBOT=.FALSE.
         ORBITAL=.FALSE.
         BOTEXPER=.FALSE.
         TMBOT=.FALSE.
         URSELL=.FALSE.
         UFRIC=.FALSE.
         Z0=.FALSE.
         ALPHA_CH=.FALSE.
         WINDX=.FALSE.
         WINDY=.FALSE.
         CD=.FALSE.
         CURRTX=.FALSE.
         CURRTY=.FALSE.
         WATLEV=.FALSE.
         WATLEVOLD=.FALSE.
         DEPDT=.FALSE.
         DEP=.FALSE.
         WINDMAG=.FALSE.
         TAUW=.FALSE.
         TAUWX=.FALSE.
         TAUWY=.FALSE.
         TAUHF=.FALSE.
         TAUTOT=.FALSE.
         STOKESBOTTX=.FALSE.
         STOKESBOTTY=.FALSE.
         STOKESSURFX=.FALSE.
         STOKESSURFY=.FALSE.
         STOKESBAROX=.FALSE.
         STOKESBAROY=.FALSE.
         RSXX=.FALSE.
         RSXY=.FALSE.
         RSYY=.FALSE.
         CFL1=.FALSE.
         CFL2=.FALSE.
         CFL3=.FALSE.
         BEGTC = MAIN%BEGT
         DELTC = -1
         UNITC = MAIN%UNIT
         ENDTC = MAIN%ENDT
         READ(INP%FHNDL, NML = HISTORY)
         wwm_print_namelist(HISTORY)
         CALL FLUSH(CHK%FHNDL)
         IF (DELTC.lt.MAIN%DELT) THEN
           DELTC=MAIN%DELT
         END IF
#ifdef NCDF
# ifdef MPI_PARALL_GRID
         IF (myrank .gt. 0) THEN
           GRIDWRITE=.FALSE.
         END IF
# endif
         PARAMWRITE_HIS=PARAMWRITE
         USE_SINGLE_OUT_HIS=USE_SINGLE_OUT
         MULTIPLEOUT_HIS=MULTIPLEOUT
         IF (rkind.eq.4) THEN
           NF90_OUTTYPE_HIS=NF90_REAL
         ELSE
           IF (USE_SINGLE_OUT_HIS) THEN
             NF90_OUTTYPE_HIS=NF90_REAL
           ELSE
             NF90_OUTTYPE_HIS=NF90_DOUBLE
           ENDIF
         ENDIF
         OUT_HISTORY % DEFINETC=DEFINETC
         IF (DEFINETC .lt. 0) THEN
           OUT_HISTORY % IDEF = -1
         ELSE
           OUT_HISTORY % IDEF = NINT(DEFINETC/DELTC)
         ENDIF
# ifdef MPI_PARALL_GRID
         IF (MULTIPLEOUT_HIS.eq.1) THEN
           WriteOutputProcess_his=.TRUE.
         ELSE
           IF (MULTIPLEOUT_HIS.eq.0) THEN
             IF (myrank.eq.0) THEN
               WriteOutputProcess_his=.TRUE.
             ELSE
               WriteOutputProcess_his=.FALSE.
             ENDIF
           ELSE
             CALL WWM_ABORT('You must have MULTIPLEOUT=0 or 1')
           ENDIF
         ENDIF
# else
#  ifndef MERGE_OPERATION
         IF (MULTIPLEOUT_HIS.ne.0) THEN
           CALL WWM_ABORT('In Serial for history, you need MULTIPLEOUT=0')
         ENDIF
         WriteOutputProcess_his=.TRUE.
#  endif
# endif
#endif
         OUT_HISTORY%BEGT = BEGTC
         OUT_HISTORY%DELT = DELTC
         OUT_HISTORY%UNIT = UNITC
         OUT_HISTORY%ENDT = ENDTC

         IF (OUT_HISTORY%BEGT .LT. MAIN%BEGT) OUT_HISTORY%BEGT = MAIN%BEGT
         IF (OUT_HISTORY%ENDT .GT. MAIN%ENDT) OUT_HISTORY%ENDT = MAIN%ENDT

         CALL CT2MJD(OUT_HISTORY%BEGT, OUT_HISTORY%BMJD)
         CALL CT2MJD(OUT_HISTORY%ENDT, OUT_HISTORY%EMJD)
         CALL CU2SEC(OUT_HISTORY%UNIT, OUT_HISTORY%DELT)


         OUT_HISTORY%TOTL = (OUT_HISTORY%EMJD - OUT_HISTORY%BMJD) * DAY2SEC
         OUT_HISTORY%ISTP = NINT( OUT_HISTORY%TOTL / OUT_HISTORY%DELT ) + 1
         OUT_HISTORY%TMJD = OUT_HISTORY%BMJD
!
! set the output flag
!
         VAROUT_HISTORY%IOUTP = 1
         IF (     TRIM(OUTSTYLE) == 'NO') THEN
            VAROUT_HISTORY%IOUTP = 0
         ELSE IF (TRIM(OUTSTYLE) == 'XFN') THEN
            VAROUT_HISTORY%IOUTP = 1
         ELSE IF (TRIM(OUTSTYLE) == 'NC') THEN
            VAROUT_HISTORY%IOUTP = 2
         ELSE IF (TRIM(OUTSTYLE) == 'SHP') THEN
            VAROUT_HISTORY%IOUTP = 3
         END IF

         IF (   TRIM(FILEOUT) == 'zorglub') THEN
           IF (     TRIM(OUTSTYLE) == 'XFN') THEN
              FILEOUT='XFNout'
           ELSE IF (TRIM(OUTSTYLE) == 'NC') THEN
              FILEOUT='WWM_output.nc'
           ELSE IF (TRIM(OUTSTYLE) == 'SHP') THEN
              FILEOUT='SHPout'
           END IF
         ENDIF

         OUT_HISTORY%FNAME = FILEOUT

         LVAR_READ( 1)=HS
         LVAR_READ( 2)=TM01
         LVAR_READ( 3)=TM02
         LVAR_READ( 4)=TM10
         LVAR_READ( 5)=KLM
         LVAR_READ( 6)=WLM
         LVAR_READ( 7)=ETOTC
         LVAR_READ( 8)=ETOTS
         LVAR_READ( 9)=DM
         LVAR_READ(10)=DSPR
         LVAR_READ(11)=TPPD
         LVAR_READ(12)=CPPD
         LVAR_READ(13)=KPPD
         LVAR_READ(14)=CGPD
         LVAR_READ(15)=TPP
         LVAR_READ(16)=CPP
         LVAR_READ(17)=WNPP
         LVAR_READ(18)=CGPP
         LVAR_READ(19)=KPP
         LVAR_READ(20)=LPP
         LVAR_READ(21)=PEAKD
         LVAR_READ(22)=PEAKDSPR
         LVAR_READ(23)=DPEAK
         LVAR_READ(24)=UBOT
         LVAR_READ(25)=ORBITAL
         LVAR_READ(26)=BOTEXPER
         LVAR_READ(27)=TMBOT
         LVAR_READ(28)=URSELL
         LVAR_READ(29)=UFRIC
         LVAR_READ(30)=Z0
         LVAR_READ(31)=ALPHA_CH
         LVAR_READ(32)=WINDX
         LVAR_READ(33)=WINDY
         LVAR_READ(34)=CD
         LVAR_READ(35)=CURRTX
         LVAR_READ(36)=CURRTY
         LVAR_READ(37)=WATLEV
         LVAR_READ(38)=WATLEVOLD
         LVAR_READ(39)=DEPDT
         LVAR_READ(40)=DEP
         LVAR_READ(41)=WINDMAG
         LVAR_READ(42)=TAUW
         LVAR_READ(43)=TAUWX
         LVAR_READ(44)=TAUWY
         LVAR_READ(45)=TAUHF
         LVAR_READ(46)=TAUTOT
         LVAR_READ(47)=STOKESBOTTX
         LVAR_READ(48)=STOKESBOTTY
         LVAR_READ(49)=STOKESSURFX
         LVAR_READ(50)=STOKESSURFY
         LVAR_READ(51)=STOKESBAROX
         LVAR_READ(52)=STOKESBAROY
         LVAR_READ(53)=RSXX
         LVAR_READ(54)=RSXY
         LVAR_READ(55)=RSYY
         LVAR_READ(56)=CFL1
         LVAR_READ(57)=CFL2
         LVAR_READ(58)=CFL3
         VAROUT_HISTORY%LVAR=LVAR_READ
         CALL DETERMINE_NEEDED_COMPUTATION(VAROUT_HISTORY)
         IF (.not. LCFL) THEN
           IF (CFL1.or.CFL2.or.CFL3) THEN
             CALL WWM_ABORT('You need to select LCFL=T if asking for LCFLx')
           ENDIF
         ENDIF
!
!     **** STATION section
!
#ifdef NCDF
         MULTIPLEOUT=0
         USE_SINGLE_OUT=.TRUE.
         DEFINETC=-1
#endif
         FILEOUT = "zorglub"
         AC=.FALSE.
         WK=.FALSE.
         ACOUT_1D=.FALSE.
         ACOUT_2D=.FALSE.
         HS=.FALSE.
         TM01=.FALSE.
         TM02=.FALSE.
         TM10=.FALSE.
         KLM=.FALSE.
         WLM=.FALSE.
         ETOTC=.FALSE.
         ETOTS=.FALSE.
         DM=.FALSE.
         DSPR=.FALSE.
         TPPD=.FALSE.
         CPPD=.FALSE.
         KPPD=.FALSE.
         CGPD=.FALSE.
         TPP=.FALSE.
         CPP=.FALSE.
         WNPP=.FALSE.
         CGPP=.FALSE.
         KPP=.FALSE.
         LPP=.FALSE.
         PEAKD=.FALSE.
         PEAKDSPR=.FALSE.
         DPEAK=.FALSE.
         UBOT=.FALSE.
         ORBITAL=.FALSE.
         BOTEXPER=.FALSE.
         TMBOT=.FALSE.
         URSELL=.FALSE.
         UFRIC=.FALSE.
         Z0=.FALSE.
         ALPHA_CH=.FALSE.
         WINDX=.FALSE.
         WINDY=.FALSE.
         CD=.FALSE.
         CURRTX=.FALSE.
         CURRTY=.FALSE.
         WATLEV=.FALSE.
         WATLEVOLD=.FALSE.
         DEPDT=.FALSE.
         DEP=.FALSE.
         WINDMAG=.FALSE.
         TAUW=.FALSE.
         TAUWX=.FALSE.
         TAUWY=.FALSE.
         TAUHF=.FALSE.
         TAUTOT=.FALSE.
         STOKESBOTTX=.FALSE.
         STOKESBOTTY=.FALSE.
         STOKESSURFX=.FALSE.
         STOKESSURFY=.FALSE.
         STOKESBAROX=.FALSE.
         STOKESBAROY=.FALSE.
         RSXX=.FALSE.
         RSXY=.FALSE.
         RSYY=.FALSE.
         CFL1=.FALSE.
         CFL2=.FALSE.
         CFL3=.FALSE.

         BEGTC = MAIN%BEGT
         DELTC = MAIN%DELT
         UNITC = MAIN%UNIT
         ENDTC = MAIN%ENDT
         READ(INP%FHNDL, NML = STATION)
         wwm_print_namelist(STATION)
         CALL FLUSH(CHK%FHNDL)
#ifdef NCDF
         PARAMWRITE_STAT=PARAMWRITE
         USE_SINGLE_OUT_STAT=USE_SINGLE_OUT
         MULTIPLEOUT_STAT=MULTIPLEOUT
         IF (rkind.eq.4) THEN
           NF90_OUTTYPE_STAT=NF90_REAL
         ELSE
           IF (USE_SINGLE_OUT_STAT) THEN
             NF90_OUTTYPE_STAT=NF90_REAL
           ELSE
             NF90_OUTTYPE_STAT=NF90_DOUBLE
           ENDIF
         ENDIF
         OUT_STATION % DEFINETC=DEFINETC
         IF (DEFINETC .lt. 0) THEN
           OUT_STATION % IDEF = -1
         ELSE
           OUT_STATION % IDEF = NINT(DEFINETC/DELTC)
         ENDIF
# ifdef MPI_PARALL_GRID
         IF (MULTIPLEOUT_STAT.eq.1) THEN
           WriteOutputProcess_stat=.TRUE.
         ELSE
           IF (MULTIPLEOUT_STAT.eq.0) THEN
             IF (myrank.eq.0) THEN
               WriteOutputProcess_stat=.TRUE.
             ELSE
               WriteOutputProcess_stat=.FALSE.
             ENDIF
           ELSE
             CALL WWM_ABORT('Station: You must have MULTIPLEOUT=0 or 1')
           ENDIF
         ENDIF
# else
#  ifndef MERGE_OPERATION
         IF (MULTIPLEOUT_STAT.ne.0) THEN
           CALL WWM_ABORT('In serial for station, you need MULTIPLEOUT=0')
         ENDIF
         WriteOutputProcess_stat=.TRUE.
#  endif
# endif
#endif
!         IF (DELTC.lt.0.0_rkind) THEN
!           CALL WWM_ABORT("DELTC is not an optional argument for STATION")
!         END IF
         OUT_STATION%BEGT = BEGTC
         OUT_STATION%DELT = DELTC
         OUT_STATION%UNIT = UNITC
         OUT_STATION%ENDT = ENDTC

         IF (OUT_STATION%BEGT .LT. MAIN%BEGT) OUT_STATION%BEGT = MAIN%BEGT
         IF (OUT_STATION%ENDT .GT. MAIN%ENDT) OUT_STATION%ENDT = MAIN%ENDT

         CALL CT2MJD(OUT_STATION%BEGT, OUT_STATION%BMJD)
         CALL CT2MJD(OUT_STATION%ENDT, OUT_STATION%EMJD)
         CALL CU2SEC(OUT_STATION%UNIT, OUT_STATION%DELT)

         OUT_STATION%TOTL = (OUT_STATION%EMJD - OUT_STATION%BMJD) * DAY2SEC
         OUT_STATION%ISTP = NINT( OUT_STATION%TOTL / OUT_STATION%DELT ) + 1
         OUT_STATION%TMJD = OUT_STATION%BMJD
         VAROUT_STATION%IOUTP = 1
         IF (     TRIM(OUTSTYLE) == 'NO') THEN
            VAROUT_STATION%IOUTP = 0
            LOUTS=.FALSE.
         ELSE IF (TRIM(OUTSTYLE) == 'STE') THEN
            VAROUT_STATION%IOUTP = 1
            LOUTS=.TRUE.
         ELSE IF (TRIM(OUTSTYLE) == 'NC') THEN
            VAROUT_STATION%IOUTP = 2
            LOUTS=.TRUE.
         END IF
         IF (IOUTS.eq.0) THEN
           LOUTS=.FALSE.
         END IF
         IF (   TRIM(FILEOUT) == 'zorglub') THEN
           IF (     TRIM(OUTSTYLE) == 'STE') THEN
              FILEOUT='STEout'
           ELSE IF (TRIM(OUTSTYLE) == 'NC') THEN
              FILEOUT='WWM_stat.nc'
           END IF
         ENDIF
         OUT_STATION%FNAME = FILEOUT
         VAROUT_STATION%AC=AC
         VAROUT_STATION%WK=WK
         VAROUT_STATION%ACOUT_1D=ACOUT_1D
         VAROUT_STATION%ACOUT_2D=ACOUT_2D
         LVAR_READ( 1)=HS
         LVAR_READ( 2)=TM01
         LVAR_READ( 3)=TM02
         LVAR_READ( 4)=TM10
         LVAR_READ( 5)=KLM
         LVAR_READ( 6)=WLM
         LVAR_READ( 7)=ETOTC
         LVAR_READ( 8)=ETOTS
         LVAR_READ( 9)=DM
         LVAR_READ(10)=DSPR
         LVAR_READ(11)=TPPD
         LVAR_READ(12)=CPPD
         LVAR_READ(13)=KPPD
         LVAR_READ(14)=CGPD
         LVAR_READ(15)=TPP
         LVAR_READ(16)=CPP
         LVAR_READ(17)=WNPP
         LVAR_READ(18)=CGPP
         LVAR_READ(19)=KPP
         LVAR_READ(20)=LPP
         LVAR_READ(21)=PEAKD
         LVAR_READ(22)=PEAKDSPR
         LVAR_READ(23)=DPEAK
         LVAR_READ(24)=UBOT
         LVAR_READ(25)=ORBITAL
         LVAR_READ(26)=BOTEXPER
         LVAR_READ(27)=TMBOT
         LVAR_READ(28)=URSELL
         LVAR_READ(29)=UFRIC
         LVAR_READ(30)=Z0
         LVAR_READ(31)=ALPHA_CH
         LVAR_READ(32)=WINDX
         LVAR_READ(33)=WINDY
         LVAR_READ(34)=CD
         LVAR_READ(35)=CURRTX
         LVAR_READ(36)=CURRTY
         LVAR_READ(37)=WATLEV
         LVAR_READ(38)=WATLEVOLD
         LVAR_READ(39)=DEPDT
         LVAR_READ(40)=DEP
         LVAR_READ(41)=WINDMAG
         LVAR_READ(42)=TAUW
         LVAR_READ(43)=TAUWX
         LVAR_READ(44)=TAUWY
         LVAR_READ(45)=TAUHF
         LVAR_READ(46)=TAUTOT
         LVAR_READ(47)=STOKESBOTTX
         LVAR_READ(48)=STOKESBOTTY
         LVAR_READ(49)=STOKESSURFX
         LVAR_READ(50)=STOKESSURFY
         LVAR_READ(51)=STOKESBAROX
         LVAR_READ(52)=STOKESBAROY
         LVAR_READ(53)=RSXX
         LVAR_READ(54)=RSXY
         LVAR_READ(55)=RSYY
         LVAR_READ(56)=CFL1
         LVAR_READ(57)=CFL2
         LVAR_READ(58)=CFL3
         VAROUT_STATION%LVAR=LVAR_READ
         CALL DETERMINE_NEEDED_COMPUTATION(VAROUT_STATION)
         IF (.not. LCFL) THEN
           IF (CFL1.or.CFL2.or.CFL3) THEN
             CALL WWM_ABORT('You need to select LCFL=T if asking for LCFLx')
           ENDIF
         ENDIF
         IF (LOUTS) THEN

           ALLOCATE ( STATION_P(IOUTS), stat=istat)
           IF (istat/=0) CALL WWM_ABORT('wwm_input, allocate error 1')
#ifdef MPI_PARALL_GRID
           STATION_P%IFOUND = 0
           STATION_P%ISUM   = 0
#endif
           DO I = 1, IOUTS
             STATION_P(I)%XELE = 0.
             STATION_P(I)%YELE = 0.
             STATION_P(I)%ZELE = 0.
             STATION_P(I)%WI   = ZERO
             STATION_P(I)%OUTPAR_NODE = 0.
             STATION_P(I)%NAME   = NOUTS(I)
             STATION_P(I)%XCOORD = XOUTS(I)
             STATION_P(I)%YCOORD = YOUTS(I)
             IF (LEN_TRIM(NOUTS(I)).eq.0) THEN
                WRITE(wwmerr,*) 'Station ', I, ' has incorrect name'
                CALL WWM_ABORT(wwmerr)
             ENDIF
           END DO

           IF (LSIGMAX) THEN
             STATION_P(1:IOUTS)%CUTOFF = CUTOFF(1:IOUTS)
           ELSE
             STATION_P(1:IOUTS)%CUTOFF = FRHIGH
           END IF

           WRITE(DBG%FHNDL,*) 'STATION X and Y Coordinates'
           WRITE(DBG%FHNDL,*) STATION_P%XCOORD
           WRITE(DBG%FHNDL,*) STATION_P%YCOORD
           WRITE(DBG%FHNDL,*) 'STATION Names'
           WRITE(DBG%FHNDL,*) STATION_P%NAME
           CALL FLUSH(DBG%FHNDL)

         END IF

         IF (LLOUTS) THEN

           ALLOCATE ( LINES(IOUTS), stat=istat)
           IF (istat/=0) CALL WWM_ABORT('wwm_input, allocate error 2')
#ifdef MPI_PARALL_GRID
           STATION_P%IFOUND = 0
           STATION_P%ISUM   = 0
#endif
           DO I = 1, IOUTS
             STATION_P(I)%XELE = 0.
             STATION_P(I)%YELE = 0.
             STATION_P(I)%ZELE = 0.
             STATION_P(I)%WI   = ZERO
             STATION_P(I)%OUTPAR_NODE = 0.
           END DO
           STATION_P(1:IOUTS)%NAME   = NOUTS(1:IOUTS)
           STATION_P(1:IOUTS)%XCOORD = XOUTS(1:IOUTS)
           STATION_P(1:IOUTS)%YCOORD = YOUTS(1:IOUTS)

           IF (LSIGMAX) THEN
             STATION_P(1:IOUTS)%CUTOFF = CUTOFF(1:IOUTS)
           ELSE
             STATION_P(1:IOUTS)%CUTOFF = FRHIGH
           END IF

           WRITE(DBG%FHNDL,*) 'STATION X and Y Coordinates'
           WRITE(DBG%FHNDL,*) STATION_P%XCOORD
           WRITE(DBG%FHNDL,*) STATION_P%YCOORD
           WRITE(DBG%FHNDL,*) 'STATION Names'
           WRITE(DBG%FHNDL,*) STATION_P%NAME
           CALL FLUSH(DBG%FHNDL)

         END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_MNP_MNE
      USE DATAPOOL, only : MNP, MNE, IGRIDTYPE, FILEGRID, rkind, GRD, STAT, DBG
      implicit none
      integer ISTAT, I, ITMP
      REAL(rkind) ATMP, BTMP, CTMP
      CHARACTER(LEN=100)  :: RHEADER

      OPEN(GRD%FHNDL, FILE=FILEGRID, STATUS = 'OLD')
      IF (IGRIDTYPE == 1) THEN ! system.dat format
        DO I = 1, 2
          READ(GRD%FHNDL, '(A)') RHEADER
        END DO
        READ(GRD%FHNDL, *, IOSTAT = ISTAT) ITMP
        IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=1 error in read mnp/mne')
        READ(GRD%FHNDL, '(A)') RHEADER
        READ(GRD%FHNDL, *, IOSTAT = ISTAT) MNP
        IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=1 error in read mnp/mne')
        DO I = 1, 7
          READ(GRD%FHNDL, '(A)') RHEADER
        END DO
        DO I = 1,MNP
          READ(GRD%FHNDL, *, IOSTAT = ISTAT) ITMP, ATMP, BTMP, CTMP
          IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=1 error in read mnp/mne')
        END DO
        DO I = 1, 2
          READ(GRD%FHNDL, '(A)') RHEADER
        END DO
        READ(GRD%FHNDL, *, IOSTAT = ISTAT) MNE
      ELSE IF (IGRIDTYPE == 2) THEN ! symbolic Mathieu format
        READ(GRD%FHNDL, *, IOSTAT = ISTAT) MNE, MNP
        IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=2 error in read mnp/mne')
      ELSE IF (IGRIDTYPE == 3) THEN ! selfe gr3
        READ(GRD%FHNDL,*)
        READ(GRD%FHNDL,*, IOSTAT = ISTAT) MNE, MNP
        IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=3 error in read mnp/mne')
      ELSE IF (IGRIDTYPE == 4) THEN ! old WWM format
        READ(GRD%FHNDL, *, IOSTAT = ISTAT) MNE 
        READ(GRD%FHNDL, *, IOSTAT = ISTAT) MNP 
        IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=4 error in read mnp/mne')
      ELSE
        CALL WWM_ABORT('IGRIDTYPE WRONG')
      END IF
!      WRITE(STAT%FHNDL,*) 'THE GRIDSIZE IS MNP, MNE:', MNP, MNE

      IF (MNP == 0 .OR. MNE == 0) THEN
        CALL WWM_ABORT('No Nodes and No Elements have been read in')
      ENDIF

      CLOSE(GRD%FHNDL)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_WWMINPUT()
#ifdef NCDF
         USE NETCDF
#endif
         USE DATAPOOL
#ifdef MPI_PARALL_GRID
         use elfe_msgp, only : myrank,parallel_abort
         use elfe_glbl, only : ics
#endif
         IMPLICIT NONE

         LOGICAL           :: LFLIVE
         CHARACTER(LEN=20) :: BEGTC, UNITC, ENDTC
         REAL(rkind)            :: DELTC

         INTEGER           :: I, NI(3)
         REAL(rkind)              :: DEG
         INTEGER :: MULTIPLEIN, MULTIPLEOUT
         NAMELIST /PROC/ PROCNAME, DIMMODE, LSTEA, LQSTEA, LSPHE,       &
     &      LNAUTIN, LNAUTOUT, LMONO_OUT, LMONO_IN,                     &
     &      BEGTC, DELTC, UNITC, ENDTC, DMIN

         NAMELIST /COUPL/ LCPL, LROMS, LTIMOR, LSHYFEM, RADFLAG,        &
     &      LETOT, NLVT, DTCOUP, IMET_DRY

         NAMELIST /GRID/ LCIRD, LSTAG, MINDIR, MAXDIR, MDC, FRLOW,      &
     &      FRHIGH, MSC, FILEGRID, IGRIDTYPE, LSLOP, SLMAX, LVAR1D,     &
     &      LOPTSIG

         NAMELIST /INIT/ LHOTR, LINID, INITSTYLE

         NAMELIST /BOUC/ LBCSE, LBCWA, LBCSP, LINHOM, LBSP1D,           &
     &      LBSP2D, LBINTER, BEGTC, DELTC, UNITC, ENDTC,                &
     &      FILEBOUND, IBOUNDFORMAT, FILEWAVE, LINDSPRDEG, LPARMDIR,    &
     &      WBHS, WBTP, WBDM, WBDS, WBSS, WBDSMS, WBGAUSS, WBPKEN,      &
     &      NCDF_HS_NAME, NCDF_DIR_NAME, NCDF_SPR_NAME, NCDF_FP_NAME,   &
     &      NCDF_F02_NAME

         NAMELIST /WIND/ LSEWD, LSTWD, LCWIN, LWDIR, BEGTC, DELTC,      &
     &      UNITC, ENDTC, LINTERWD, WDIR, WVEL, CWINDX, CWINDY,         &
     &      FILEWIND, WINDFAC, IWINDFORMAT, LWINDFROMWWM

         NAMELIST /CURR/ LSECU, BEGTC, DELTC, UNITC, ENDTC,             &
     &      LINTERCU, LSTCU, LCCUR, CCURTX, CCURTY, FILECUR,            &
     &      LERGINP, CURFAC, ICURRFORMAT

         NAMELIST /WALV/ LSEWL, BEGTC, DELTC, UNITC, ENDTC,             &
     &      LINTERWL, LSTWL, LCWLV, CWATLV, FILEWATL, LERGINP,          &
     &      WALVFAC, IWATLVFORMAT

         NAMELIST /ENGS/ MESNL, MESIN, IFRIC, MESBF, FRICC,             &
     &      MESBR, ICRIT, ALPBJ, BRHD,                                  &
     &      LMAXETOT, MESDS, MESTR, TRICO, TRIRA, TRIURS

         NAMELIST /NUMS/ ICOMP, AMETHOD, SMETHOD, DMETHOD,              &
     &      LITERSPLIT, LFILTERTH, MAXCFLTH, LTHBOUND, FMETHOD,         &
     &      LFILTERCXY, MAXCFLCXY, LFILTERSIG, MAXCFLSIG, LSIGBOUND,    &
     &      LLIMT, LIMFAK, MELIM, LDIFR, IDIFFR, LADVTEST, LSOUBOUND,   &
     &      LCONV, LCFL, LCHKCONV, NQSITER, QSCONV1, QSCONV2,           &
     &      QSCONV3, QSCONV4, QSCONV5, EPSH1, EPSH2, EPSH3, EPSH4,      &
     &      EPSH5, RTHETA, LEXPIMP, FREQEXP, LVECTOR,IVECTOR,           &
     &      DTMIN_DYN, NDYNITER, DTMIN_SIN, DTMIN_SNL4,                 &
     &      DTMIN_SDS, DTMIN_SNL3, DTMIN_SBR, DTMIN_SBF,                &
     &      NDYNITER_SIN, NDYNITER_SNL4, NDYNITER_SDS, NDYNITER_SBR,    &
     &      NDYNITER_SNL3, NDYNITER_SBF, NB_BLOCK, SOLVERTHR, LNANINFCHK

         NAMELIST /HOTFILE/ BEGTC, DELTC, UNITC, ENDTC, LHOTF,          &
     &      LCYCLEHOT, FILEHOT_OUT, HOTSTYLE_IN, HOTSTYLE_OUT,          &
     &      MULTIPLEIN, MULTIPLEOUT, IHOTPOS_IN, FILEHOT_IN

         READ( INP%FHNDL,  NML = PROC)
         wwm_print_namelist(PROC)
         CALL FLUSH(CHK%FHNDL)
#ifdef SELFE
         IF (LSPHE) THEN
           IF (ics /= 2) THEN
             WRITE(DBG%FHNDL) LSPHE, ICS
             CALL FLUSH(DBG%FHNDL)
             CALL WWM_ABORT('You set LSPHE=T but then you need ics=2')
           END IF
         ELSE
           IF (ics /= 1) THEN
             WRITE(DBG%FHNDL) LSPHE, ICS
             CALL FLUSH(DBG%FHNDL)
             CALL WWM_ABORT('You set LSPHE=F but then you need ics=1')
           END IF
         END IF
#endif
         READ( INP%FHNDL,  NML = COUPL)
         wwm_print_namelist(COUPL)
         CALL FLUSH(CHK%FHNDL)
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
         MAIN%TOTL = (MAIN%EMJD - MAIN%BMJD) * DAY2SEC
         MAIN%ISTP = NINT( MAIN%TOTL / MAIN%DELT )
         MAIN%TMJD = MAIN%BMJD

         IF (MAIN%DELT .LT. THR) CALL WWM_ABORT('TIME STEP IS ZERO')

#ifdef SELFE
         IF (DIMMODE .NE. 2 .OR. .NOT. LCPL) THEN
           WRITE(wwmerr,*)'You are running in less than 1d or LCPL = .F.',&
      &     DIMMODE,LCPL
           call parallel_abort(wwmerr)
         endif
!         MAIN%DELT   = DT_WWM   !DELTC
         MAIN%DTCOUP = DT_SELFE !coupling time step
#endif
!
!     *** GRID section
!
         READ (INP%FHNDL,   NML = GRID)
         wwm_print_namelist(GRID)
         CALL FLUSH(CHK%FHNDL)
#ifdef MPI_PARALL_GRID
         IF (TRIM(FILEGRID) /= 'hgrid.gr3') THEN
           CALL WWM_ABORT('In parallel mode you need FILEGRID=hgrid.gr3')
         END IF
#endif
         GRD%FNAME = FILEGRID

         NSPEC=MDC*MSC
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

#ifdef MPI_PARALL_GRID
!AR: I think we can allow that with SELFE we can use different GRIDTYPES even if selfe reads in .gr3 ... let's see ...
         IF (IGRIDTYPE /= 3) CALL WWM_ABORT('In MPI, you need IGRIDTYPE=3')
#endif
#ifndef MPI_PARALL_GRID
         CALL READ_MNP_MNE
         NP_RES=MNP
#endif

         IF (FRLOW > FRHIGH) THEN
           CALL WWM_ABORT('error, the FRHIG must be greater than FRLOW')
         END IF
!
!     *** INIT section
!
         READ(INP%FHNDL,  NML = INIT)
         wwm_print_namelist(INIT)
         CALL FLUSH(CHK%FHNDL)

         IF (LHOTR) THEN
           WRITE(STAT%FHNDL,'("+TRACE...",A)') 'HOTFILE is used as Initital Condition'
         END IF
!
!     *** BOUNDARY CONDITIONS section
!
         READ(INP%FHNDL,  NML = BOUC )
         wwm_print_namelist(BOUC)
         CALL FLUSH(CHK%FHNDL)

         BND%FNAME = FILEBOUND
         WAV%FNAME = FILEWAVE
         IF (LBCWA .and. (.not. LINHOM) .and. (.not. LBCSE)) THEN
           ! Parametric Wave Boundary is prescribed
           ! Inhomogenous in space
           ! Steady in time
           IF (WBTP*FRLOW .gt. 1) THEN
             CALL WWM_ABORT('FRLOW is too high with respect to WBTP')
           END IF
           IF (WBTP*FRHIGH .lt. 1) THEN
             CALL WWM_ABORT('FRHIGH is too low with respect to WBTP')
           END IF
         END IF
         WRITE(STAT%FHNDL,'("+TRACE...",A10,I5)') 'BOUNDARY FILE FORMAT IS', IBOUNDFORMAT

         SEBO%BEGT = BEGTC
         SEBO%DELT = DELTC
         SEBO%UNIT = UNITC
         SEBO%ENDT = ENDTC

!AR: correcting input and output time steps ...
         IF (SEBO%BEGT .LT. MAIN%BEGT) SEBO%BEGT = MAIN%BEGT
         IF (SEBO%ENDT .GT. MAIN%ENDT) SEBO%ENDT = MAIN%ENDT

         CALL CT2MJD(SEBO%BEGT, SEBO%BMJD)
         CALL CT2MJD(SEBO%ENDT, SEBO%EMJD)

         CALL CU2SEC(SEBO%UNIT, SEBO%DELT)

         SEBO%TMJD = SEBO%BMJD
!
!     *** WIND section
!
         READ(INP%FHNDL, NML = WIND)
         wwm_print_namelist(WIND)
         CALL FLUSH(CHK%FHNDL)
!         Print *, 'BEGTC', BEGTC
!         Print *, 'ENDTC', ENDTC


         WIN%FNAME = TRIM(FILEWIND)
         IF (LWINDFROMWWM .and. (LCWIN .eqv. .FALSE.)) THEN
           CALL TEST_FILE_EXIST_DIE("Missing wind file : ", WIN%FNAME)
         END IF

         SEWI%BEGT = BEGTC
         SEWI%DELT = DELTC
         SEWI%UNIT = UNITC
         SEWI%ENDT = ENDTC

         IF (SEWI%BEGT .LT. MAIN%BEGT) SEWI%BEGT = MAIN%BEGT
         IF (SEWI%ENDT .GT. MAIN%ENDT) SEWI%ENDT = MAIN%ENDT
!         Print *, 'SEWI%BEGT=', SEWI%BEGT
!         Print *, 'SEWI%ENDT=', SEWI%ENDT

         CALL CT2MJD(SEWI%BEGT, SEWI%BMJD)
         CALL CT2MJD(SEWI%ENDT, SEWI%EMJD)

         CALL CU2SEC(SEWI%UNIT, SEWI%DELT)

         SEWI%TMJD = 0.0
!
!     *** CURR section
!
         READ(INP%FHNDL, NML = CURR)
         wwm_print_namelist(CURR)
         CALL FLUSH(CHK%FHNDL)

         CUR%FNAME = TRIM(FILECUR)

         IF (LSECU .AND. LSTCU) THEN
           CALL wwm_abort('Error: LSECU and LSTCU are .TRUE. value')
         END IF

         SECU%BEGT = BEGTC
         SECU%DELT = DELTC
         SECU%UNIT = UNITC
         SECU%ENDT = ENDTC

         IF (SECU%BEGT .LT. MAIN%BEGT) SECU%BEGT = MAIN%BEGT
         IF (SECU%ENDT .GT. MAIN%ENDT) SECU%ENDT = MAIN%ENDT

         CALL CT2MJD(SECU%BEGT, SECU%BMJD) ! convert string to internal time ... in double ...
         CALL CT2MJD(SECU%ENDT, SECU%EMJD)

         CALL CU2SEC(SECU%UNIT, SECU%DELT)

         SECU%TMJD = 0.0
!
!     *** water level section
!
         READ(INP%FHNDL, NML = WALV)
         wwm_print_namelist(WALV)
         CALL FLUSH(CHK%FHNDL)

         WAT%FNAME = FILEWATL

         IF (LSEWL .AND. LSTWL) THEN
           CALL wwm_abort('Error: LSEWL and LSTWL have .TRUE. value')
         END IF

         SEWL%BEGT = BEGTC
         SEWL%DELT = DELTC
         SEWL%UNIT = UNITC
         SEWL%ENDT = ENDTC

         IF (SEWL%BEGT .LT. MAIN%BEGT) SEWL%BEGT = MAIN%BEGT
         IF (SEWL%ENDT .GT. MAIN%ENDT) SEWL%ENDT = MAIN%ENDT

         CALL CT2MJD(SEWL%BEGT, SEWL%BMJD)
         CALL CT2MJD(SEWL%ENDT, SEWL%EMJD)

         CALL CU2SEC(SEWL%UNIT, SEWL%DELT)
         SEWL%TMJD = 0.0
!
!     *** ENGS section
!
         READ(INP%FHNDL, NML = ENGS)
         wwm_print_namelist(ENGS)
         CALL FLUSH(CHK%FHNDL)

!
!     *** NUMS section
!
         READ(INP%FHNDL, NML = NUMS)
         wwm_print_namelist(NUMS)
         CALL FLUSH(CHK%FHNDL)
         CALL READ_HISTORY_STATION_NAMELIST()
!
!     **** HOTFILE section

#ifdef NCDF
         IF (rkind == 4) THEN
           NF90_RUNTYPE=NF90_REAL
         ELSE
           NF90_RUNTYPE=NF90_DOUBLE
         ENDIF
#endif
         FILEHOT_OUT='wwm_hotfile_out.nc'
         FILEHOT_IN='wwm_hotfile_in.nc'
         
         BEGTC = MAIN%BEGT
         DELTC = -1
         UNITC = MAIN%UNIT
         ENDTC = MAIN%ENDT
         MULTIPLEOUT=0
         MULTIPLEIN=0
         READ(INP%FHNDL, NML = HOTFILE)
         wwm_print_namelist(HOTFILE)
         CALL FLUSH(CHK%FHNDL)

         MULTIPLEIN_HOT=MULTIPLEIN
         MULTIPLEOUT_HOT=MULTIPLEOUT
#ifdef MPI_PARALL_GRID
         IF (MULTIPLEOUT_HOT.eq.1) THEN
           WriteOutputProcess_hot=.TRUE.
         ELSE
           IF (MULTIPLEOUT_HOT.eq.0) THEN
             IF (myrank.eq.0) THEN
               WriteOutputProcess_hot=.TRUE.
             ELSE
               WriteOutputProcess_hot=.FALSE.
             ENDIF
           ELSE
             CALL WWM_ABORT('Hotfile: You must have MULTIPLEOUT=0 or 1')
           ENDIF
         ENDIF
#else
# ifndef MERGE_OPERATION
         IF (MULTIPLEOUT_HOT.ne.0) THEN
           CALL WWM_ABORT('In Serial for hotfile, you must have MULTIPLEOUT=0')
         ENDIF
         WriteOutputProcess_hot=.TRUE.
# endif
#endif
         HOTOUT%FNAME = FILEHOT_OUT
         HOTIN%FNAME = FILEHOT_IN

         HOTF%BEGT = BEGTC
         HOTF%DELT = DELTC
         HOTF%UNIT = UNITC
         HOTF%ENDT = ENDTC

         IF (HOTF%BEGT .LT. MAIN%BEGT) HOTF%BEGT = MAIN%BEGT
         IF (HOTF%ENDT .GT. MAIN%ENDT) HOTF%ENDT = MAIN%ENDT

         CALL CT2MJD(HOTF%BEGT, HOTF%BMJD)
         CALL CT2MJD(HOTF%ENDT, HOTF%EMJD)

         CALL CU2SEC(HOTF%UNIT, HOTF%DELT)
         IF (HOTF%DELT.lt.MAIN%DELT) THEN
           HOTF%DELT=MAIN%DELT
         END IF

         HOTF%TOTL = (HOTF%EMJD - HOTF%BMJD) * DAY2SEC
         HOTF%ISTP = NINT( HOTF%TOTL / HOTF%DELT ) + 1
         HOTF%TMJD = HOTF%BMJD

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_CURRENT_INPUT
      USE DATAPOOL
#ifdef MPI_PARALL_GRID
      USE ELFE_MSGP, ONLY : myrank, comm, ierr
      USE ELFE_GLBL, ONLY : ipgl, NP_GLOBAL
#endif
      IMPLICIT NONE
      INTEGER :: IP, ISTAT, IT, IFILE, FORECASTHOURS
      LOGICAL :: LFLIVE
      REAL(rkind)    :: WDIRT, wrf_w1, wrf_w2
      CHARACTER(LEN=20) :: TIMESTRING
#ifdef MPI_PARALL_GRID
      INTEGER :: I
      REAL(rkind) :: tmp
      REAL(rkind) :: tmp_arr(np_global)
#endif
      CURTXY(:,:) = 0.0
      IF (LSTCU) THEN
        IF (DIMMODE .EQ. 1) THEN
          IF (LCCUR) THEN
            DO IP = 1, MNP
              CURTXY(IP,1) = CCURTX
            END DO
          ELSE
            CALL TEST_FILE_EXIST_DIE("1: Missing current file : ", CUR%FNAME)
            OPEN(CUR%FHNDL, FILE = TRIM(CUR%FNAME), STATUS = 'OLD')
            READ(CUR%FHNDL, *, IOSTAT = ISTAT) CURTXY(:,1)
            IF ( ISTAT > 0 ) CALL WWM_ABORT('error in the current velocity file')
            CLOSE(CUR%FHNDL)
          END IF
        ELSE IF (DIMMODE .EQ. 2) THEN
          IF (LCCUR) THEN
            DO IP = 1, MNP
              CURTXY(IP,1) = CCURTX
              CURTXY(IP,2) = CCURTY
            END DO
          ELSE
            CALL TEST_FILE_EXIST_DIE("2: Missing current file : ", CUR%FNAME)
            OPEN(CUR%FHNDL, FILE = TRIM(CUR%FNAME), STATUS = 'OLD')
#ifdef MPI_PARALL_GRID
            READ(CUR%FHNDL, *, IOSTAT = ISTAT) tmp_arr
            DO I=1,NP_GLOBAL
              IF (ipgl(I)%rank==myrank) THEN
                IF ( ISTAT > 0 ) CALL WWM_ABORT('error in the wind velocity file')
                CURTXY(ipgl(I)%id,1)=tmp_arr(I)
              END IF
            END DO
            READ(CUR%FHNDL, *, IOSTAT = ISTAT) tmp_arr
            DO I=1,NP_GLOBAL
              IF (ipgl(I)%rank==myrank) THEN
                IF ( ISTAT > 0 ) CALL WWM_ABORT('error in the wind velocity file')
                CURTXY(ipgl(I)%id,2)=tmp_arr(I)
              END IF
            END DO
#else
            READ(CUR%FHNDL, *, IOSTAT = ISTAT) CURTXY(:,1)
            IF ( ISTAT > 0 ) CALL WWM_ABORT('error in the current velocity file')
            READ(CUR%FHNDL, *, IOSTAT = ISTAT) CURTXY(:,2)
            IF ( ISTAT > 0 ) CALL WWM_ABORT('error in the current velocity file')
#endif
            CLOSE(CUR%FHNDL)
          END IF
        END IF
      ELSE IF (LSECU) THEN
        SECU%TOTL = (SECU%EMJD - SECU%BMJD) * DAY2SEC
        SECU%ISTP = NINT( SECU%TOTL / SECU%DELT ) + 1
        SECU%TMJD = SECU%BMJD
        LSECN = .FALSE.
        WRITE(STAT%FHNDL,*) 'Serial current Condition -----------'
        WRITE(STAT%FHNDL,*) SECU%BEGT, SECU%ENDT, SECU%ISTP, SECU%TOTL/3600.0, SECU%DELT
        IF (LERGINP) CALL ERG2WWM(SECU%ISTP)
        CALL TEST_FILE_EXIST_DIE("3: Missing current file : ", CUR%FNAME)
        LSECN = .TRUE.
        OPEN(CUR%FHNDL, FILE = TRIM(CUR%FNAME), STATUS = 'OLD')
        CALL CSEVAL( CUR%FHNDL, TRIM(CUR%FNAME), LCURFILE, 2, CURTXY)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_WATLEV_INPUT
      USE DATAPOOL
#ifdef MPI_PARALL_GRID
      USE ELFE_MSGP, ONLY : myrank, comm, ierr
      USE ELFE_GLBL, ONLY : ipgl, NP_GLOBAL
#endif
      IMPLICIT NONE

      INTEGER :: IP, ISTAT, IT, IFILE, FORECASTHOURS
      LOGICAL :: LFLIVE
      REAL(rkind)    :: WDIRT, wrf_w1, wrf_w2
      CHARACTER(LEN=20) :: TIMESTRING
#ifdef MPI_PARALL_GRID
      INTEGER :: I
      REAL(rkind)    :: tmp
      REAL(rkind)    :: tmp_arr(np_global)
#endif
      WATLEV    = 0.
      WATLEVOLD = 0.
      IF (LSTWL) THEN
        IF (DIMMODE .EQ. 1) THEN
          IF (LCWLV) THEN
            WATLEV = CWATLV
            DEP    = WLDEP + WATLEV
          ELSE
            CALL TEST_FILE_EXIST_DIE("1: Missing watlev file : ", WAT%FNAME)
            OPEN(WAT%FHNDL, FILE = TRIM(WAT%FNAME), STATUS = 'OLD')
#ifdef MPI_PARALL_GRID
            READ(WAT%FHNDL, *, IOSTAT = ISTAT) tmp_arr
            DO I=1,NP_GLOBAL
              IF (ipgl(I)%rank==myrank) THEN
                IF ( ISTAT > 0 ) CALL WWM_ABORT('error in the wind velocity file')
                WATLEV(ipgl(I)%id)=tmp_arr(I)
              END IF
            END DO
#else
            READ(WAT%FHNDL, *, IOSTAT = ISTAT) WATLEV(:)
            IF ( ISTAT > 0 )  CALL WWM_ABORT('error in the water level file')
#endif
            CLOSE(WAT%FHNDL)
          END IF
        ELSE IF (DIMMODE .EQ. 2) THEN
          IF (LCWLV) THEN
            WATLEV = CWATLV
            DEP    = WLDEP + WATLEV
          ELSE
            CALL TEST_FILE_EXIST_DIE("2: Missing watlev file : ", WAT%FNAME)
            OPEN(WAT%FHNDL, FILE = TRIM(WAT%FNAME), STATUS = 'OLD')
            READ(WAT%FHNDL, *, IOSTAT = ISTAT) WATLEV(:)
            IF ( ISTAT > 0 )  CALL WWM_ABORT('error in the water level file')
            CLOSE(WAT%FHNDL)
          END IF
        END IF
      ELSE IF (LSEWL) THEN
        SEWL%TOTL = (SEWL%EMJD - SEWL%BMJD) * DAY2SEC
        SEWL%ISTP = NINT( SEWL%TOTL / SEWL%DELT ) + 1
        SEWL%TMJD = SEWL%BMJD
        IF (LERGINP .AND. .NOT. LSECU) CALL ERG2WWM(SEWL%ISTP)
        LSELN = .FALSE.
        WRITE(STAT%FHNDL,*) 'Serial water level Condition -----------'
        WRITE(STAT%FHNDL,*) SEWL%BEGT, SEWL%ENDT, SEWL%ISTP, SEWL%TOTL/3600.0, SEWL%DELT
        CALL TEST_FILE_EXIST_DIE("LSEWL: Missing watlev file : ", WAT%FNAME)
        LSELN = .TRUE.
        OPEN(WAT%FHNDL, FILE = TRIM(WAT%FNAME), STATUS = 'OLD')
        CALL CSEVAL( WAT%FHNDL,TRIM(WAT%FNAME), LWATLFILE, 1, WATLEV)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_SPATIAL_GRID

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER :: I, IP, IE, ISTAT, ITMP, JTMP
         REAL(rkind)    :: TMP
         REAL(rkind)  :: XPDTMP, YPDTMP, ZPDTMP
         CHARACTER (LEN = 10) :: STRNGTMP
         LOGICAL :: LFLIVE

         REAL(rkind) DXP1, DXP2, DXP3, DYP1, DYP2, DYP3, DBLTMP
         INTEGER KTMP, LTMP, MTMP, NTMP, OTMP
         CHARACTER(LEN=100)              :: RHEADER
!
!2DO makes subs for reading the mesh file with different formats ...
!
         CALL TEST_FILE_EXIST_DIE('Missing grid file : ', GRD%FNAME)

         SELECT CASE (DIMMODE)
           CASE (1)
             OPEN(GRD%FHNDL, FILE = GRD%FNAME, STATUS = 'OLD')
             DO IP = 1, MNP
               READ(GRD%FHNDL, *, IOSTAT = ISTAT) XP(IP), DEP(IP)
               IF ( ISTAT /= 0 ) CALL WWM_ABORT('error in the grid configuration file')
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
             IF ((MNP.eq.0).or.(MNE.eq.0)) THEN
               CALL WWM_ABORT('We have MNP=0 or MNE=0 before reading grid')
             END IF
             IF (IGRIDTYPE == 1) THEN ! system.dat format ... XFN
               DO I = 1, 2
                 READ(GRD%FHNDL, '(A)') RHEADER
               END DO
               READ(GRD%FHNDL, *, IOSTAT = ISTAT) KTMP
               IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=1 error in grid read 1')
               READ(GRD%FHNDL, '(A)') RHEADER
               READ(GRD%FHNDL, *, IOSTAT = ISTAT) JTMP
               IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=1 error in grid read 2')
               DO I = 1, 7
                 READ(GRD%FHNDL, '(A)') RHEADER
               END DO
               DO IP=1,MNP
                 READ(GRD%FHNDL, *, IOSTAT = ISTAT) KTMP, XP(IP), YP(IP), DEP(IP)
                 IF (KTMP+1.ne.IP) THEN
                   CALL WWM_ABORT('IGRIDTYPE=1 error in grid read 3')
                 ENDIF
                 IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=1 error in grid read 4')
               END DO
               DO I = 1, 2
                 READ(GRD%FHNDL, '(A)') RHEADER
               END DO
               READ(GRD%FHNDL, *, IOSTAT = ISTAT) ITMP
               DO I = 1, 3
                 READ(GRD%FHNDL, '(A)') RHEADER
               END DO
               DO IE=1,MNE
                 READ(GRD%FHNDL, *, IOSTAT = ISTAT) KTMP, LTMP, MTMP, NTMP, OTMP
                 IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=1 error in grid read 5')
                 IF (OTMP+1.ne.IE) THEN
                   CALL WWM_ABORT('IGRIDTYPE=1 error in grid read 6')
                 ENDIF
                 INE(1,IE)=KTMP+1
                 INE(2,IE)=LTMP+1
                 INE(3,IE)=MTMP+1
               END DO
             ELSE IF (IGRIDTYPE == 2) THEN ! periodic grid written by mathieu dutour
               READ(GRD%FHNDL,*) ITMP, JTMP
               DO IP = 1, MNP
                 READ(GRD%FHNDL, *, IOSTAT = ISTAT) DEP(IP)
                 IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=2 error in grid read 1')
               END DO
               DO IE = 1, MNE
                 READ(GRD%FHNDL, *, IOSTAT = ISTAT) TRIA(IE)
                 IF ( ISTAT /= 0 )  CALL WWM_ABORT('IGRIDTYPE=2 error in grid read 2')
                 READ(GRD%FHNDL, *, IOSTAT = ISTAT) INE(:,IE)
                 IF ( ISTAT /= 0 )  CALL WWM_ABORT('IGRIDTYPE=2 error in grid read 3')
                 READ(GRD%FHNDL, *, IOSTAT = ISTAT) DXP1, DXP2, DXP3
                 IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=2 error in grid read 4')
                 READ(GRD%FHNDL, *, IOSTAT = ISTAT) DYP1, DYP2, DYP3
                 IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=2 error in grid read 5')
                 IEN(1,IE) = -DYP2
                 IEN(2,IE) = DXP2
                 IEN(3,IE) = -DYP3
                 IEN(4,IE) = DXP3
                 IEN(5,IE) = -DYP1
                 IEN(6,IE) = DXP1
               END DO
             ELSE IF (IGRIDTYPE == 3) THEN ! selfe gr3
               READ(GRD%FHNDL,*)
               READ(GRD%FHNDL,*) ITMP, JTMP
               DO IP = 1, MNP
                 READ(GRD%FHNDL, *, IOSTAT = ISTAT) KTMP, XPDTMP, YPDTMP, ZPDTMP
                 XP(IP)  = XPDTMP
                 YP(IP)  = YPDTMP
                 DEP(IP) = ZPDTMP
                 IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=3 error in grid reading 1')
               END DO
               DO IE = 1, MNE
                 READ(GRD%FHNDL, *, IOSTAT = ISTAT) KTMP, LTMP, INE(:,IE)
                 IF ( ISTAT /= 0 )  CALL WWM_ABORT('IGRIDTYPE=3 error in grid reading 2')
               END DO
             ELSE IF (IGRIDTYPE == 4) THEN ! Old WWM format
               READ(GRD%FHNDL,*) ITMP, JTMP 
               DO IP = 1, MNP
                 READ(GRD%FHNDL, *, IOSTAT = ISTAT) XP(IP), YP(IP), DEP(IP)
                 IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=1 error in grid read 1')
               END DO
               DO IE = 1, MNE
                 READ(GRD%FHNDL, *, IOSTAT = ISTAT) INE(:,IE)
                 IF ( ISTAT /= 0 )  CALL WWM_ABORT('IGRIDTYPE=1 error in grid read 2')
               END DO
             ELSE
               CALL WWM_ABORT('IGRIDTYPE WRONG')
             END IF
             IF ((ITMP.ne.MNE).or.(JTMP.ne.MNP)) THEN
               WRITE(DBG%FHNDL,*) 'ITMP=', ITMP, 'MNE=', MNE
               WRITE(DBG%FHNDL,*) 'JTMP=', JTMP, 'MNP=', MNP
               CALL WWM_ABORT('Inconsistency in reading')
             END IF
             CLOSE(GRD%FHNDL)
           CASE DEFAULT
               CALL WWM_ABORT('WRONG GRID DIMENSION')
         END SELECT
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CHECK_LOGICS()
         USE DATAPOOL
#ifdef MPI_PARALL_GRID
         use elfe_msgp, only : myrank,parallel_abort
#endif
         IMPLICIT NONE

         REAL(rkind) :: TEST
         REAL(rkind) :: DIFFTIME

!        Check timings ...

         IF (OUT_HISTORY%BMJD .GE. OUT_HISTORY%EMJD) CALL WWM_ABORT('CHECK OUTPUT HISTORY TIME STEPS BEGINN TIME STEP IS SMALLER THAN END TIME STEP')
         IF (OUT_STATION%BMJD .GE. OUT_STATION%EMJD) CALL WWM_ABORT('CHECK OUTPUT STATION TIME STEPS BEGINN TIME STEP IS SMALLER THAN END TIME STEP')
         IF (MAIN%BMJD .GE. MAIN%EMJD) CALL WWM_ABORT('CHECK COMPUTATION TIME STEPS BEGINN TIME STEP IS SMALLER THAN END TIME STEP')
         IF (LSEWD) THEN
           IF (SEWI%BMJD .GE. SEWI%EMJD) CALL WWM_ABORT('CHECK WIND TIME STEPS BEGINN TIME STEP IS SMALLER THAN END TIME STEP')
           IF (MAIN%BMJD .LT. SEWI%BMJD) CALL WWM_ABORT('Start time of Wind after begin of run')
           IF (MAIN%EMJD .GT. SEWI%EMJD) CALL WWM_ABORT('End time of Wind before end of run')
         END IF
         IF (LSECU) THEN
           IF (SECU%BMJD .GE. SECU%EMJD) CALL WWM_ABORT('CHECK CURRENT TIME STEPS BEGINN TIME STEP IS SMALLER THAN END TIME STEP')
         END IF
         IF (LHOTF) THEN
           IF (HOTF%BMJD .GE. HOTF%EMJD) CALL WWM_ABORT('CHECK HOTFILE TIME STEPS BEGINN TIME STEP IS SMALLER THAN END TIME STEP')
         END IF
         IF (SEBO%BMJD .GE. SEBO%EMJD) CALL WWM_ABORT('CHECK BOUNDARY TIME STEPS BEGINN TIME STEP IS SMALLER THAN END TIME STEP')
         

!        Check MSC,MDC for exchange
         if(MSC<1.or.MDC<1) call wwm_abort('MSC,MDC too small')

         IF (SMETHOD .GT. 0) THEN

           IF (MESIN .GT. 6) THEN
             call wwm_abort('CHECK NUMS - MESIN OUT OF RANGE')
           END IF

           IF (MESBR .GT. 1) THEN
             call wwm_abort('CHECK NUMS - MESBR OUT OF RANGE')
           END IF
#ifdef PGMCL_COUPLING
           IF (.NOT.LCPL) THEN
             CALL WWM_ABORT('LCPL=T if running with PGMCL')
           ENDIF
           IF (LROMS.or.LTIMOR.or.LSHYFEM) THEN
             CALL WWM_ABORT('LROMS=LTIMOR=LSHYFEM=F if with ROMS_PGMCL')
           ENDIF
#elif SELFE
           IF (.NOT. LCPL) THEN
             CALL WWM_ABORT('LCPL=T if running with SELFE')
           ENDIF
#endif
           IF (MESBF .GT. 0 .AND. FRICC .LT. 0.) THEN
             call wwm_abort('CHECK NUMS - FRICTION COEFFICIENT HAS WRONG SIGN')
           END IF

#ifndef PETSC
           IF (AMETHOD.eq.4) THEN
             call wwm_abort('If AMETHOD=4 then you need PETSC')
           ENDIF
           IF (AMETHOD.eq.5) THEN
             call wwm_abort('If AMETHOD=5 then you need PETSC')
           ENDIF
#endif
#ifndef MPI_PARALL_GRID
           IF (AMETHOD.eq.6) THEN
             call wwm_abort('If AMETHOD=6 then you need MPI')
           ENDIF
#endif
           IF (MESTR .GT. 7) THEN
             call wwm_abort('CHECK NUMS - MESTR OUT OF RANGE')
           END IF

           IF (MESNL .GT. 6) THEN
             call wwm_abort('CHECK NUMS - MESNL OUT OF RANGE')
           END IF

           IF (MESDS .GT. 6) THEN
             call wwm_abort('CHECK NUMS - MESDS OUT OF RANGE')
           END IF

           IF (MESNL .GT. 0 .AND. .NOT. LLIMT ) THEN
!AR: this will be a warning ...
             !call wwm_abort('YOU ARE USING SNL WITHOUT LIMITER CODE WILL STOP NOW')
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

#if !defined PGMCL_COUPLING && !defined SELFE
         IF (LCPL) THEN
           IF (.NOT. LROMS .AND. .NOT. LSHYFEM .AND. .NOT. LTIMOR) THEN
             CALL WWM_ABORT('LROMS, LSHYFEM or LTIMOR must be true')
           END IF
         END IF
#endif
#ifndef PGMCL_COUPLING
         IF (LCPL) THEN
#endif
#ifndef SELFE
           IF (MAIN%DTCOUP .LT. MAIN%DELT) THEN
             CALL WWM_ABORT('COUPLE TIME STEP IS SMALLER AS THE CALCULATION TIME STEP!')
           END IF
           TEST = MAIN%DTCOUP - NINT(MAIN%DTCOUP/MAIN%DELT)*MAIN%DELT
!2do ... check where else you do some nint stuff ... like that one ...
           IF (ABS(TEST) .GT. THR) THEN
             WRITE(DBG%FHNDL,*) 'MAIN%DTCOUP=', MAIN%DTCOUP
             WRITE(DBG%FHNDL,*) 'MAIN%DELT=', MAIN%DELT
             WRITE(DBG%FHNDL,*) 'TEST=', TEST
             WRITE(DBG%FHNDL,*) 'TIME STEP OF THE WAVEMODELL CANNOT BE DiVIDIED WITHOUT A REST'
             CALL WWM_ABORT('TIME STEP OF THE WAVEMODELL CANNOT BE DiVIDIED WITHOUT A REST')
           ELSE
             MAIN%ICPLT = INT(MAIN%DTCOUP/MAIN%DELT)
           END IF
#endif
#ifndef PGMCL_COUPLING
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
           CALL WWM_ABORT('CHECK LSEWL OR LSTDW')
         END IF

         IF (LSTCU .AND. LSECU) THEN
           WRITE(DBG%FHNDL,*) 'YOU MUST USE EITHER UNSTEADY OR STEADY CURRENTS'
           WRITE(DBG%FHNDL,*) 'PLEASE CHECK CODE EXITS'
           CALL WWM_ABORT('CHECK LSTCU .AND. LSECU')
         END IF

         IF (LSTCU .AND. LSECU) THEN
           WRITE(DBG%FHNDL,*) 'YOU MUST USE EITHER UNSTEADY OR STEADY CURRENTS'
           WRITE(DBG%FHNDL,*) 'PLEASE CHECK CODE EXITS'
           CALL WWM_ABORT('CHECK LSTCU .AND. LSECU')
         END IF

         IF (LSTWL .AND. LSEWL) THEN
           WRITE(DBG%FHNDL,*) 'YOU MUST USE EITHER UNSTEADY OR STEADY CURRENTS'
           WRITE(DBG%FHNDL,*) 'PLEASE CHECK CODE EXITS'
           CALL WWM_ABORT('CHECK LSTCU .AND. LSECU')
         END IF

         ! if using PETSc with AMETHOD 4 or 5 ICOMP must be greater equal 1 to init JA IA
         IF (AMETHOD .GE. 4 .AND. ICOMP .LT. 1) THEN
           call wwm_abort('ICMP must be greater equal 1 to use PETSc')
         END IF

         ! if using PETSc with AMETHOD 4 or 5 LVECTOR must be FALSE
         IF (AMETHOD .GE. 4 .AND. LVECTOR .EQV. .TRUE.) THEN
           call wwm_abort('LVECTOR must be FALSE to use PETSc')
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
      SUBROUTINE INIT_NETCDF_WW3_WAVEPARAMETER
         USE DATAPOOL
         USE NETCDF
         IMPLICIT NONE

        INTEGER :: ISTAT, IT, IX, IY, IFILE, IVAR, BND_NCID
        INTEGER :: ILON_ID, ILAT_ID, ITIME_ID, I, J, COUNTER
        INTEGER :: D_WIND_X_ID, D_WIND_Y_ID, D_PRESS_ID
        REAL(rkind)  :: DTMP, DTMP1, DTMP2
        REAL(rkind), ALLOCATABLE :: BND_TIME(:)
        character ( len = 40 ) chrtmp
        character ( len = 15 ) chrdate
        character ( len = 40 ) netcfd_fname
        character ( len = 20 ) dirname
        character ( len = 100) chrerr
        character (len = *), parameter :: CallFct = "INIT_NETCDF_WW3_WAVEPARAMETER"

        integer, dimension(nf90_max_var_dims) :: dimIDs

        OPEN(WAV%FHNDL,FILE=WAV%FNAME,STATUS='OLD')
!
! count number of netcdf files in list ...
!
        WRITE(STAT%FHNDL,*) WAV%FHNDL, WAV%FNAME, BND%FHNDL, BND%FNAME

        NUM_NETCDF_FILES_BND = 0
        DO
          READ( WAV%FHNDL, *, IOSTAT = ISTAT )
          IF ( ISTAT /= 0 ) EXIT
          NUM_NETCDF_FILES_BND = NUM_NETCDF_FILES_BND + 1
        END DO
        REWIND (WAV%FHNDL)

        WRITE(STAT%FHNDL,*) NUM_NETCDF_FILES_BND

        NUM_NETCDF_FILES_BND = NUM_NETCDF_FILES_BND / NUM_NETCDF_VAR_TYPES

        WRITE(STAT%FHNDL,*) 'NUM_NETCDF_FILES_BND', NUM_NETCDF_FILES_BND

        ALLOCATE(NETCDF_FILE_NAMES_BND(NUM_NETCDF_FILES_BND,NUM_NETCDF_VAR_TYPES), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_input, allocate error 3')


        DO IT = 1, NUM_NETCDF_FILES_BND
          DO IVAR = 1, NUM_NETCDF_VAR_TYPES
            READ( WAV%FHNDL, *) NETCDF_FILE_NAMES_BND(IT,IVAR)
          END DO
        END DO
        CLOSE (WAV%FHNDL)
!
! four files are read to set up the wave spectra Hs, Tm01, Dir, Sprd
!
        ALLOCATE(NDT_BND_FILE(NUM_NETCDF_FILES_BND), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_input, allocate error 4')
        NDT_BND_FILE = 0

        DO IFILE = 1, NUM_NETCDF_FILES_BND
          WRITE(STAT%FHNDL,'(I10,10X,5A30)') IFILE, NETCDF_FILE_NAMES_BND(IFILE,:)
        END DO
!
! check number of time steps in netcdf file ... it is assumed that all files have the same ammount of time steps ...
!
        DO IFILE = 1, NUM_NETCDF_FILES_BND
          write(STAT%FHNDL,*) ifile, TRIM(NETCDF_FILE_NAMES_BND(IFILE,1))
          ISTAT = NF90_OPEN(TRIM(NETCDF_FILE_NAMES_BND(IFILE,1)), NF90_NOWRITE, BND_NCID)
          CALL GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)

          ISTAT = nf90_inq_varid(BND_NCID, 'time', ITIME_ID)
          CALL GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)

          ISTAT = NF90_INQUIRE_VARIABLE(BND_NCID, ITIME_ID, dimids = dimids)
          CALL GENERIC_NETCDF_ERROR(CallFct, 3, ISTAT)

          ISTAT = nf90_inquire_dimension(BND_NCID, dimIDs(1), len = NDT_BND_FILE(IFILE))
          CALL GENERIC_NETCDF_ERROR(CallFct, 4, ISTAT)

          write(STAT%FHNDL,*) IFILE, NDT_BND_FILE(IFILE)
        END DO
!
! check dimensions in the netcdf ... again it is assumed that this is not changing for all files ...
!
        ISTAT = nf90_inq_varid(BND_NCID, 'longitude', ILON_ID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 5, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(BND_NCID, ILON_ID, dimids = dimIDs)
        CALL GENERIC_NETCDF_ERROR(CallFct, 6, ISTAT)

        ISTAT = nf90_inquire_dimension(BND_NCID, dimIDs(1), len = NDX_BND)
        CALL GENERIC_NETCDF_ERROR(CallFct, 7, ISTAT)

        ISTAT = nf90_inq_varid(BND_NCID, 'latitude', ILAT_ID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 8, ISTAT)

        ISTAT = NF90_INQUIRE_VARIABLE(BND_NCID, ILAT_ID, dimids = dimIDs)
        CALL GENERIC_NETCDF_ERROR(CallFct, 9, ISTAT)

        ISTAT = nf90_inquire_dimension(BND_NCID, dimIDs(1), len = NDY_BND)
        CALL GENERIC_NETCDF_ERROR(CallFct, 10, ISTAT)

        WRITE(STAT%FHNDL,*) 'Number of Gridpoints', NDX_BND, NDY_BND

        ALLOCATE (COORD_BND_X(NDX_BND), COORD_BND_Y(NDY_BND), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_input, allocate error 5')
!
! read cooridantes from files ....
!
        ISTAT = NF90_GET_VAR(BND_NCID, ILON_ID, COORD_BND_X)
        CALL GENERIC_NETCDF_ERROR(CallFct, 11, ISTAT)

        ISTAT = NF90_GET_VAR(BND_NCID, ILAT_ID, COORD_BND_Y)
        CALL GENERIC_NETCDF_ERROR(CallFct, 12, ISTAT)
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
        CALL GENERIC_NETCDF_ERROR(CallFct, 13, ISTAT)
!
! total number of time steps ... in all files
!
        NDT_BND_ALL_FILES = 0
        write(STAT%FHNDL,*) NUM_NETCDF_FILES_BND
        DO IT = 1, NUM_NETCDF_FILES_BND
          NDT_BND_ALL_FILES = NDT_BND_ALL_FILES + NDT_BND_FILE(IT)
          write(STAT%FHNDL,*) it, NDT_BND_FILE(it)
        END DO

        WRITE(STAT%FHNDL,*) NDT_BND_ALL_FILES, NDT_BND_FILE

        ALLOCATE (BND_TIME_ALL_FILES(NUM_NETCDF_FILES_BND,MAXVAL(NDT_BND_FILE)), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_input, allocate error 6')
        BND_TIME_ALL_FILES = ZERO
!
! read all time steps in the proper format and transform in wwm time line
!
        BND_TIME_ALL_FILES = 0.
        DO IFILE = 1, NUM_NETCDF_FILES_BND
          ISTAT = NF90_OPEN(NETCDF_FILE_NAMES_BND(IFILE,1),NF90_NOWRITE,BND_NCID)
          CALL GENERIC_NETCDF_ERROR(CallFct, 14, ISTAT)

          ALLOCATE (BND_TIME(NDT_BND_FILE(IFILE)), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_input, allocate error 7')
          BND_TIME = ZERO
! MDS: It looks dangerous to use previous id.
          ISTAT = NF90_GET_VAR(BND_NCID,ITIME_ID,BND_TIME)
          CALL GENERIC_NETCDF_ERROR(CallFct, 15, ISTAT)

          DO IT = 1, NDT_BND_FILE(IFILE)
             BND_TIME_ALL_FILES(IFILE,IT) = BND_TIME(IT)
!             CALL CT2MJD('19000101.000000',DTMP1)
!             CALL CT2MJD('19900101.000000',DTMP2)
!             CALL MJD2CT(DTMP1,chrdate)
!             WRITE(*,*) '19000101.000000', DTMP1, chrdate
!             CALL MJD2CT(DTMP2,chrdate)
!             WRITE(*,*) '19900101.000000', DTMP2, chrdate
!             CALL MJD2CT(0.0_rkind,chrdate)
!             WRITE(*,*) '00000000.000000', 0.0_rkind, chrdate
!             WRITE(*,*) BND_TIME_ALL_FILES(1,1), DT_DIFF_19901900
!             IF (IT == 1 .AND. IFILE ==1) WRITE(*,*) DTMP1, DTMP2, DTMP1+DT_DIFF_19901900
!             IF (IT == 1 .AND. IFILE ==1) WRITE(*,*) IFILE, IT, BND_TIME(IT), chrdate
          END DO
          DEALLOCATE(BND_TIME)
        END DO ! IFILE

        SEBO%DELT = (BND_TIME_ALL_FILES(1,2) - BND_TIME_ALL_FILES(1,1)) * DAY2SEC
        write(STAT%FHNDL,*) SEBO%DELT, BND_TIME_ALL_FILES(1,2), BND_TIME_ALL_FILES(1,1)

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
      SUBROUTINE READ_NETCDF_WW3(IFILE,IT,CALLEDFROM)

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
         REAL(rkind),   ALLOCATABLE  :: TMP(:,:)
         REAL(rkind), ALLOCATABLE    :: U(:), V(:), H(:)
         REAL(rkind), SAVE           :: TIME, scale_factor
         character (len = *), parameter :: CallFct = "READ_NETCDF_WW3"
         INTEGER, DIMENSION (nf90_max_var_dims) :: dimIDs
         CHARACTER(LEN=80)    :: CHRTMP
         CHARACTER(LEN=100)   :: CHRERR
         CHARACTER(LEN=25)    :: CALLEDFROM

         ALLOCATE (ITMP(NDX_BND,NDY_BND), stat=istat)
         IF (istat/=0) CALL WWM_ABORT('wwm_input, allocate error 8')

         WRITE(DBG%FHNDL,*) IT, IFILE, 'READING GLOBAL DATA'
         ISTAT = NF90_OPEN(NETCDF_FILE_NAMES_BND(IFILE,3),NF90_NOWRITE,HS_BND_NCID)
         CALL GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)

         ISTAT = nf90_inq_varid(HS_BND_NCID, TRIM(NCDF_HS_NAME), HS_WW3_ID)
         CALL GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)

         ISTAT = nf90_get_att(HS_BND_NCID, HS_WW3_ID, 'scale_factor', scale_factor)
         CALL GENERIC_NETCDF_ERROR(CallFct, 3, ISTAT)

         IF (.NOT. ALLOCATED(HS_WW3)) THEN
           ALLOCATE (HS_WW3(NDX_BND,NDY_BND), stat=istat)
           IF (istat/=0) CALL WWM_ABORT('wwm_input, allocate error 9')
           HS_WW3 = 0.
         END IF
         ISTAT = NF90_GET_VAR(HS_BND_NCID, HS_WW3_ID, ITMP,  start = (/ 1, 1, IT /), count = (/ NDX_BND, NDY_BND, 1/))
         CALL GENERIC_NETCDF_ERROR(CallFct, 4, ISTAT)

         HS_WW3 = MyREAL(ITMP) * scale_factor
         ISTAT = nf90_close(HS_BND_NCID)
         CALL GENERIC_NETCDF_ERROR(CallFct, 5, ISTAT)

         ISTAT = NF90_OPEN(NETCDF_FILE_NAMES_BND(IFILE,2),NF90_NOWRITE,FP_BND_NCID)
         CALL GENERIC_NETCDF_ERROR(CallFct, 6, ISTAT)

         ISTAT = nf90_inq_varid(FP_BND_NCID, TRIM(NCDF_FP_NAME), FP_WW3_ID)
         CALL GENERIC_NETCDF_ERROR(CallFct, 7, ISTAT)

         ISTAT = nf90_get_att(FP_BND_NCID, FP_WW3_ID, 'scale_factor', scale_factor)
         CALL GENERIC_NETCDF_ERROR(CallFct, 8, ISTAT)

         IF (.NOT. ALLOCATED(FP_WW3)) THEN
           ALLOCATE (FP_WW3(NDX_BND,NDY_BND), stat=istat)
           IF (istat/=0) CALL WWM_ABORT('wwm_input, allocate error 10')
           FP_WW3 = 0.
         END IF
         ISTAT = NF90_GET_VAR(FP_BND_NCID, FP_WW3_ID, ITMP,  start = (/ 1, 1, IT /), count = (/ NDX_BND, NDY_BND, 1/))
         CALL GENERIC_NETCDF_ERROR(CallFct, 9, ISTAT)

         FP_WW3 = MyREAL(ITMP) * scale_factor
         ISTAT = nf90_close(FP_BND_NCID)
         CALL GENERIC_NETCDF_ERROR(CallFct, 10, ISTAT)

         ISTAT = NF90_OPEN(NETCDF_FILE_NAMES_BND(IFILE,5),NF90_NOWRITE,T02_BND_NCID)
         CALL GENERIC_NETCDF_ERROR(CallFct, 11, ISTAT)

         ISTAT = nf90_inq_varid(T02_BND_NCID, TRIM(NCDF_F02_NAME), T02_WW3_ID)
         CALL GENERIC_NETCDF_ERROR(CallFct, 12, ISTAT)

         ISTAT = nf90_get_att(T02_BND_NCID, T02_WW3_ID, 'scale_factor', scale_factor)
         CALL GENERIC_NETCDF_ERROR(CallFct, 13, ISTAT)

         IF (.NOT. ALLOCATED(T02_WW3)) THEN
           ALLOCATE (T02_WW3(NDX_BND,NDY_BND), stat=istat)
           IF (istat/=0) CALL WWM_ABORT('wwm_input, allocate error 11')
           T02_WW3 = 0.
         END IF
         ISTAT = NF90_GET_VAR(T02_BND_NCID, T02_WW3_ID, ITMP,  start = (/ 1, 1, IT /), count = (/ NDX_BND, NDY_BND, 1/))
         CALL GENERIC_NETCDF_ERROR(CallFct, 14, ISTAT)

         T02_WW3 = MyREAL(ITMP) * scale_factor
         ISTAT = nf90_close(T02_BND_NCID)
         CALL GENERIC_NETCDF_ERROR(CallFct, 15, ISTAT)

         ISTAT = NF90_OPEN(NETCDF_FILE_NAMES_BND(IFILE,4),NF90_NOWRITE,DSPR_BND_NCID)
         CALL GENERIC_NETCDF_ERROR(CallFct, 16, ISTAT)

         ISTAT = nf90_inq_varid(DSPR_BND_NCID, TRIM(NCDF_SPR_NAME), DSPR_WW3_ID)
         CALL GENERIC_NETCDF_ERROR(CallFct, 17, ISTAT)

         ISTAT = nf90_get_att(DSPR_BND_NCID, DSPR_WW3_ID, 'scale_factor', scale_factor)
         CALL GENERIC_NETCDF_ERROR(CallFct, 18, ISTAT)

         IF (.NOT. ALLOCATED(DSPR_WW3)) THEN
           ALLOCATE (DSPR_WW3(NDX_BND,NDY_BND), stat=istat)
           IF (istat/=0) CALL WWM_ABORT('wwm_input, allocate error 12')
           DSPR_WW3 = 0.
         END IF
         ISTAT = NF90_GET_VAR(DSPR_BND_NCID, DSPR_WW3_ID, ITMP,  start = (/ 1, 1, IT /), count = (/ NDX_BND, NDY_BND, 1/))
         CALL GENERIC_NETCDF_ERROR(CallFct, 19, ISTAT)

         DSPR_WW3 = MyREAL(ITMP) * scale_factor
         ISTAT = nf90_close(DSPR_BND_NCID)
         CALL GENERIC_NETCDF_ERROR(CallFct, 20, ISTAT)

         ISTAT = NF90_OPEN(NETCDF_FILE_NAMES_BND(IFILE,1),NF90_NOWRITE,DIR_BND_NCID)
         CALL GENERIC_NETCDF_ERROR(CallFct, 21, ISTAT)

         ISTAT = nf90_inq_varid(DIR_BND_NCID, TRIM(NCDF_DIR_NAME), DIR_WW3_ID)
         CALL GENERIC_NETCDF_ERROR(CallFct, 22, ISTAT)

         ISTAT = nf90_get_att(DIR_BND_NCID, DIR_WW3_ID, 'scale_factor', scale_factor)
         CALL GENERIC_NETCDF_ERROR(CallFct, 23, ISTAT)

         IF (.NOT. ALLOCATED(DIR_WW3)) THEN
           ALLOCATE (DIR_WW3(NDX_BND,NDY_BND), stat=istat)
           IF (istat/=0) CALL WWM_ABORT('wwm_input, allocate error 13')
           DIR_WW3 = 0.
         END IF
         ISTAT = NF90_GET_VAR(DIR_BND_NCID, DIR_WW3_ID, ITMP,  start = (/ 1, 1, IT /), count = (/ NDX_BND, NDY_BND, 1/))
         CALL GENERIC_NETCDF_ERROR(CallFct, 24, ISTAT)

         DIR_WW3 = MyREAL(ITMP) * scale_factor
         ISTAT = nf90_close(DIR_BND_NCID)
         CALL GENERIC_NETCDF_ERROR(CallFct, 25, ISTAT)

         IF (LWRITE_WW3_RESULTS) THEN
           OPEN(3012, FILE  = 'ergwiii.bin', FORM = 'UNFORMATTED')
           IF (.NOT. ALLOCATED(U)) THEN
             ALLOCATE(U(NDX_BND*NDY_BND), stat=istat)
             IF (istat/=0) CALL WWM_ABORT('wwm_input, allocate error 14')
           END IF
           IF (.NOT. ALLOCATED(V)) THEN
             ALLOCATE(V(NDX_BND*NDY_BND), stat=istat)
             IF (istat/=0) CALL WWM_ABORT('wwm_input, allocate error 15')
           END IF
           IF (.NOT. ALLOCATED(H)) THEN
             ALLOCATE(H(NDX_BND*NDY_BND), stat=istat)
             IF (istat/=0) CALL WWM_ABORT('wwm_input, allocate error 16')
           END IF
           COUNTER = 1
           DO J = 1, NDY_BND
             DO I = 1, NDX_BND
               U(COUNTER) = HS_WW3(I,J)
               V(COUNTER) = DIR_WW3(I,J)
               H(COUNTER) = DSPR_WW3(I,J)
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
          REAL(rkind), INTENT(IN)    :: MAT(NDX,NDY)
          REAL(rkind), INTENT(IN)    :: DX, DY, OFFSET_X, OFFSET_Y
          REAL(rkind), INTENT(OUT)   :: VAL(MNP)

          INTEGER             :: IP, J_INT, I_INT
          REAL(rkind)                :: WX1, WX2, WX3, WX4, HX1, HX2
          REAL(rkind)                :: DELTA_X, DELTA_Y, LEN_X, LEN_Y

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
          REAL(rkind), INTENT(IN)    :: DX, DY, OFFSET_X, OFFSET_Y
          REAL(rkind), INTENT(OUT)   :: VAL(8,IWBMNP)

          INTEGER             :: IP, J_INT, I_INT
          REAL(rkind)                :: WX1, WX2, WX3, WX4, HX1, HX2
          REAL(rkind)                :: DELTA_X, DELTA_Y, LEN_X, LEN_Y

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

            IF (WX1 .LT. 0. .OR. WX2 .LT. 0. .OR. WX3 .LT. 0. .OR. WX4 .LT. 0. ) THEN
              VAL(2,IP) = 0.
            ELSE
              VAL(2,IP) = HX1 + (HX2-HX1)/DY * DELTA_Y
            ENDIF

            IF (VAL(2,IP) .GT. TINY(1.)) THEN
              VAL(2,IP) = 1. / VAL(2,IP)
            ELSE
              VAL(2,IP) = 0.
            END IF

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
          REAL(rkind), INTENT(IN)    :: DX, DY, OFFSET_X, OFFSET_Y
          REAL(rkind), INTENT(OUT)   :: VAL(8,MNP)

          INTEGER             :: IP, J_INT, I_INT
          REAL(rkind)                :: WX1, WX2, WX3, WX4, HX1, HX2
          REAL(rkind)                :: DELTA_X, DELTA_Y, LEN_X, LEN_Y

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

            VAL(5,IP)  = 2. ! From mean period ...
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
