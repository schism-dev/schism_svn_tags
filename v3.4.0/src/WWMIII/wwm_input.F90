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
         USE DATAPOOL, only : OUTSTYLE, LOUTITER,                       &
     &        LENERGY, LWXFN, OUT_HISTORY, OUT_STATION, INP, LSP1D,     &
     &        LSP2D, INP, CHK, LOUTS, IOUTS, LSIGMAX, LLOUTS,           &
     &        ILOUTS, OUT, DAY2SEC, FRHIGH, DBG, LINES, VAROUT_HISTORY, &
     &        VAROUT_STATION, GRIDWRITE, RKIND, LVAR_READ,              &
     &        PARAMWRITE_HIS, PARAMWRITE_STAT, wwmerr, LCFL, myrank,    &
     &        istat
#ifdef NCDF
         USE NETCDF
         USE DATAPOOL, only : USE_SINGLE_OUT_STAT, USE_SINGLE_OUT_HIS,  &
     &        MULTIPLEOUT_HIS, MULTIPLEOUT_STAT,                        &
     &        NF90_OUTTYPE_STAT, NF90_OUTTYPE_HIS
#endif
         USE DATAPOOL, only : STATION_P => STATION
         USE DATAPOOL, only : IOBPD => IOBPD_HISTORY
         USE DATAPOOL, only : MAIN, PRINTMMA, ZERO
         USE DATAPOOL, only : WriteOutputProcess_hot
         USE DATAPOOL, only : WriteOutputProcess_his
         USE DATAPOOL, only : WriteOutputProcess_stat

         IMPLICIT NONE
         CHARACTER(LEN=40)  :: FILEOUT
         INTEGER, PARAMETER :: INUMOUTS = 200 
         CHARACTER(LEN=20)  :: BEGTC, UNITC, ENDTC, NOUTS(INUMOUTS), NLOUTS(INUMOUTS)
         REAL(rkind)        :: XOUTS(INUMOUTS), YOUTS(INUMOUTS), CUTOFF(INUMOUTS)
         REAL(rkind)        :: XLOUTS(INUMOUTS), YLOUTS(INUMOUTS)
         REAL(rkind) :: DEFINETC
         INTEGER     :: MULTIPLEOUT
         LOGICAL     :: USE_SINGLE_OUT
         REAL(rkind) :: DELTC
         INTEGER     :: I
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
     &      RSXX, RSXY, RSYY, CFL1, CFL2, CFL3, ZETA_SETUP

         NAMELIST /HISTORY/ BEGTC, DELTC, UNITC, ENDTC, DEFINETC,       &
     &      OUTSTYLE, FILEOUT, LOUTITER, IOBPD,                         &
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
     &      RSXX, RSXY, RSYY, CFL1, CFL2, CFL3, ZETA_SETUP

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
     &      RSXX, RSXY, RSYY, CFL1, CFL2, CFL3, ZETA_SETUP

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
         ZETA_SETUP=.FALSE.
         BEGTC = MAIN%BEGT
         DELTC = -1
         UNITC = MAIN%UNIT
         ENDTC = MAIN%ENDT
         READ(INP%FHNDL, NML = HISTORY)
         wwm_print_namelist(HISTORY)
         FLUSH(CHK%FHNDL)
         IF (DELTC.lt.MAIN%DELT) THEN
           DELTC=MAIN%DELT
         END IF
#ifdef NCDF
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
         LVAR_READ(59)=ZETA_SETUP
         VAROUT_HISTORY%LVAR=LVAR_READ
         CALL DETERMINE_NEEDED_COMPUTATION(VAROUT_HISTORY)
         IF (.not. LCFL) THEN
           IF (CFL1.or.CFL2.or.CFL3) THEN
             CALL WWM_ABORT('You need to select LCFL=T if asking for CFLx')
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
         ZETA_SETUP=.FALSE.
         BEGTC = MAIN%BEGT
         DELTC = MAIN%DELT
         UNITC = MAIN%UNIT
         ENDTC = MAIN%ENDT
         READ(INP%FHNDL, NML = STATION)
         wwm_print_namelist(STATION)
         FLUSH(CHK%FHNDL)
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
         IF ( TRIM(OUT_STATION%FNAME) == TRIM(OUT_HISTORY%FNAME) ) THEN
           WRITE(DBG%FHNDL,*) 'OUT_STATION%FNAME=', TRIM(OUT_STATION%FNAME)
           WRITE(DBG%FHNDL,*) 'OUT_HISTORY%FNAME=', TRIM(OUT_HISTORY%FNAME)
           CALL WWM_ABORT('You cannot have same name for history and station')
         END IF
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
         LVAR_READ(59)=ZETA_SETUP
         VAROUT_STATION%LVAR=LVAR_READ
         CALL DETERMINE_NEEDED_COMPUTATION(VAROUT_STATION)
         IF (.not. LCFL) THEN
           IF (CFL1.or.CFL2.or.CFL3) THEN
             CALL WWM_ABORT('You need to select LCFL=T if asking for CFLx')
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
           FLUSH(DBG%FHNDL)

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
           FLUSH(DBG%FHNDL)

         END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_WWMINPUT
#ifdef NCDF
         USE NETCDF
#endif
         USE DATAPOOL
#ifdef SELFE
         use elfe_glbl, only : ics
#endif
         IMPLICIT NONE

         CHARACTER(LEN=20) :: BEGTC, UNITC, ENDTC
         CHARACTER(LEN=20) :: BEGTC_OUT, UNITC_OUT, ENDTC_OUT
         CHARACTER(LEN=140) :: NETCDF_OUT_FILE

         REAL(rkind)            :: DELTC, DELTC_OUT

         REAL(rkind)              :: DEG
         INTEGER :: MULTIPLEIN, MULTIPLEOUT
         LOGICAL :: MULTIPLE_IN
         LOGICAL :: NETCDF_OUT_PARAM, NETCDF_OUT_SPECTRA
         REAL(rkind) :: DEFINETC
         LOGICAL     :: USE_SINGLE_OUT
         NAMELIST /PROC/ PROCNAME, DIMMODE, LSTEA, LQSTEA, LSPHE,       &
     &      LNAUTIN, LNAUTOUT, LMONO_OUT, LMONO_IN,                     &
     &      BEGTC, DELTC, UNITC, ENDTC, DMIN, MULTIPLE_OUT_INFO

         NAMELIST /COUPL/ LCPL, LROMS, LTIMOR, LSHYFEM, RADFLAG,        &
     &      LETOT, NLVT, DTCOUP, IMET_DRY

         NAMELIST /GRID/ LCIRD, LSTAG, MINDIR, MAXDIR, MDC, FRLOW,      &
     &      FRHIGH, MSC, FILEGRID, IGRIDTYPE, LSLOP, SLMAX, LVAR1D,     &
     &      LOPTSIG, CART2LATLON, LATLON2CART 

         NAMELIST /INIT/ LHOTR, LINID, INITSTYLE

         NAMELIST /BOUC/ LBCSE, LBCWA, LBCSP, LINHOM, LBSP1D,           &
     &      LBSP2D, LBINTER, BEGTC, DELTC, UNITC, ENDTC,                &
     &      FILEBOUND, IBOUNDFORMAT, FILEWAVE, LINDSPRDEG, LPARMDIR,    &
     &      WBHS, WBTP, WBDM, WBDS, WBSS, WBDSMS, WBGAUSS, WBPKEN,      &
     &      NCDF_HS_NAME, NCDF_DIR_NAME, NCDF_SPR_NAME, NCDF_FP_NAME,   &
     &      NCDF_F02_NAME, MULTIPLE_IN, NETCDF_OUT_PARAM,               &
     &      NETCDF_OUT_SPECTRA, NETCDF_OUT_FILE, USE_SINGLE_OUT,        &
     &      BEGTC_OUT, DELTC_OUT, UNITC_OUT, ENDTC_OUT

         NAMELIST /WIND/ LSEWD, LSTWD, LCWIN, LWDIR, BEGTC, DELTC,      &
     &      UNITC, ENDTC, LINTERWD, WDIR, WVEL, CWINDX, CWINDY,         &
     &      FILEWIND, WINDFAC, IWINDFORMAT, LWINDFROMWWM,               &
     &      SHIFT_WIND_TIME, MULTIPLE_IN

         NAMELIST /CURR/ LSECU, BEGTC, DELTC, UNITC, ENDTC,             &
     &      LINTERCU, LSTCU, LCCUR, CCURTX, CCURTY, FILECUR,            &
     &      LERGINP, CURFAC, ICURRFORMAT, MULTIPLE_IN

         NAMELIST /WALV/ LSEWL, BEGTC, DELTC, UNITC, ENDTC,             &
     &      LINTERWL, LSTWL, LCWLV, CWATLV, FILEWATL, LERGINP,          &
     &      WALVFAC, IWATLVFORMAT, MULTIPLE_IN

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
     &      NDYNITER_SNL3, NDYNITER_SBF, NB_BLOCK, SOLVERTHR, MAXITER,  &
     &      LNANINFCHK, LZETA_SETUP, ZETA_METH, LSOURCESWAM, PMIN,      &
     &      LSOURCESWWIII, BLOCK_GAUSS_SEIDEL, LNONL,                   &
     &      L_SOLVER_NORM, ASPAR_LOCAL_LEVEL

         NAMELIST /HOTFILE/ BEGTC, DELTC, UNITC, ENDTC, LHOTF,          &
     &      LCYCLEHOT, FILEHOT_OUT, HOTSTYLE_IN, HOTSTYLE_OUT,          &
     &      MULTIPLEIN, MULTIPLEOUT, IHOTPOS_IN, FILEHOT_IN

         READ( INP%FHNDL,  NML = PROC)
         wwm_print_namelist(PROC)
         FLUSH(CHK%FHNDL)
#ifdef SELFE
         IF (LSPHE) THEN
           IF (ics /= 2) THEN
             WRITE(DBG%FHNDL) LSPHE, ICS
             FLUSH(DBG%FHNDL)
             CALL WWM_ABORT('You set LSPHE=T but then you need ics=2')
           END IF
         ELSE
           IF (ics /= 1) THEN
             WRITE(DBG%FHNDL) LSPHE, ICS
             FLUSH(DBG%FHNDL)
             CALL WWM_ABORT('You set LSPHE=F but then you need ics=1')
           END IF
         END IF
#endif
         READ( INP%FHNDL,  NML = COUPL)
         wwm_print_namelist(COUPL)
         FLUSH(CHK%FHNDL)
!
!    *** Estimate various timings ...
!
         MAIN%BEGT   = BEGTC
         MAIN%DELT   = DELTC
         MAIN%UNIT   = UNITC
         MAIN%ENDT   = ENDTC
         CALL CT2MJD(MAIN%BEGT, MAIN%BMJD)
         CALL CT2MJD(MAIN%ENDT, MAIN%EMJD)
         CALL CU2SEC(MAIN%UNIT, MAIN%DELT)
         MAIN%TOTL = (MAIN%EMJD - MAIN%BMJD) * DAY2SEC
         MAIN%ISTP = NINT( MAIN%TOTL / MAIN%DELT )
         MAIN%TMJD = MAIN%BMJD

         MAIN%DTCOUP = DTCOUP
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
         FLUSH(CHK%FHNDL)
#if defined MPI_PARALL_GRID && !defined PDLIB
         IF (TRIM(FILEGRID) /= 'hgrid.gr3') THEN
           CALL WWM_ABORT('With SELFE parallelization you need FILEGRID=hgrid.gr3 in wwminput.nml')
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
# ifndef PDLIB
!AR: I think we can allow that with SELFE we can use different GRIDTYPES even if selfe reads in .gr3 ... let's see ...
         IF (IGRIDTYPE /= 3) CALL WWM_ABORT('In MPI, you need IGRIDTYPE=3')
# endif
#endif
         IF (FRLOW > FRHIGH) THEN
           CALL WWM_ABORT('error, the FRHIG must be greater than FRLOW')
         END IF
!
!     *** INIT section
!
         READ(INP%FHNDL,  NML = INIT)
         wwm_print_namelist(INIT)
         FLUSH(CHK%FHNDL)

         IF (LHOTR) THEN
           WRITE(STAT%FHNDL,'("+TRACE...",A)') 'HOTFILE is used as Initital Condition'
         END IF
!
!     *** BOUNDARY CONDITIONS section
!
#ifdef NCDF
         USE_SINGLE_OUT=BOUC_USE_SINGLE_OUT
         DEFINETC=-1
#endif
         MULTIPLE_IN=MULTIPLE_IN_BOUND
         BEGTC=MAIN%BEGT
         ENDTC=MAIN%ENDT
         UNITC=MAIN%UNIT
         DELTC=MAIN%DELT
         BEGTC_OUT=MAIN%BEGT
         ENDTC_OUT=MAIN%ENDT
         UNITC_OUT=MAIN%UNIT
         DELTC_OUT=MAIN%DELT
         NETCDF_OUT_FILE=BOUC_NETCDF_OUT_FILE
#ifdef NCDF
         NETCDF_OUT_PARAM  =BOUC_NETCDF_OUT_PARAM
         NETCDF_OUT_SPECTRA=BOUC_NETCDF_OUT_SPECTRA
         NETCDF_OUT_FILE=BOUC_NETCDF_OUT_FILE
#endif
         READ(INP%FHNDL,  NML = BOUC )
         wwm_print_namelist(BOUC)
         FLUSH(CHK%FHNDL)
         MULTIPLE_IN_BOUND=MULTIPLE_IN
#ifdef NCDF
         BOUC_NETCDF_OUT_PARAM  =NETCDF_OUT_PARAM
         BOUC_NETCDF_OUT_SPECTRA=NETCDF_OUT_SPECTRA
         BOUC_NETCDF_OUT_FILE   =NETCDF_OUT_FILE
         BOUC_USE_SINGLE_OUT=USE_SINGLE_OUT
         IF (rkind.eq.4) THEN
           NF90_OUTTYPE_BOUC=NF90_REAL
         ELSE
           IF (BOUC_USE_SINGLE_OUT) THEN
             NF90_OUTTYPE_BOUC=NF90_REAL
           ELSE
             NF90_OUTTYPE_BOUC=NF90_DOUBLE
           ENDIF
         ENDIF
         OUT_BOUC % BEGT=BEGTC_OUT
         OUT_BOUC % ENDT=ENDTC_OUT
         OUT_BOUC % UNIT=UNITC_OUT
         OUT_BOUC % DELT=DELTC_OUT
         OUT_BOUC % DEFINETC=DEFINETC
         OUT_BOUC % FNAME = NETCDF_OUT_FILE
         CALL CT2MJD(OUT_BOUC % BEGT, OUT_BOUC % BMJD)
         CALL CT2MJD(OUT_BOUC % ENDT, OUT_BOUC % EMJD)
         CALL CU2SEC(OUT_BOUC % UNIT, OUT_BOUC % DELT)
         IF (DEFINETC .lt. 0) THEN
           OUT_BOUC % IDEF = -1
         ELSE
           OUT_BOUC % IDEF = NINT(DEFINETC/DELTC)
         ENDIF
#endif
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
         MULTIPLE_IN=MULTIPLE_IN_WIND
         READ(INP%FHNDL, NML = WIND)
         wwm_print_namelist(WIND)
         FLUSH(CHK%FHNDL)
         MULTIPLE_IN_WIND=MULTIPLE_IN

         WIN%FNAME = TRIM(FILEWIND)
         IF (LWINDFROMWWM .and. (LCWIN .eqv. .FALSE.)) THEN
           CALL TEST_FILE_EXIST_DIE("Missing wind file : ", WIN%FNAME)
         END IF
         IF (IWINDFORMAT .ne. 1) THEN
           BEGTC = MAIN%BEGT
           DELTC = MAIN%DELT
           UNITC = MAIN%UNIT
           ENDTC = MAIN%ENDT
         END IF

         SEWI%BEGT = BEGTC
         SEWI%DELT = DELTC
         SEWI%UNIT = UNITC
         SEWI%ENDT = ENDTC

         IF (SEWI%BEGT .LT. MAIN%BEGT) SEWI%BEGT = MAIN%BEGT
         IF (SEWI%ENDT .GT. MAIN%ENDT) SEWI%ENDT = MAIN%ENDT

         CALL CT2MJD(SEWI%BEGT, SEWI%BMJD)
         CALL CT2MJD(SEWI%ENDT, SEWI%EMJD)

         CALL CU2SEC(SEWI%UNIT, SEWI%DELT)

         SEWI%TMJD = 0.0
!
!     *** CURR section
!
         MULTIPLE_IN=MULTIPLE_IN_CURR
         READ(INP%FHNDL, NML = CURR)
         wwm_print_namelist(CURR)
         FLUSH(CHK%FHNDL)
         MULTIPLE_IN_CURR=MULTIPLE_IN
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
         MULTIPLE_IN_WATLEV=MULTIPLE_IN
         READ(INP%FHNDL, NML = WALV)
         wwm_print_namelist(WALV)
         FLUSH(CHK%FHNDL)
         MULTIPLE_IN_WATLEV=MULTIPLE_IN
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
         FLUSH(CHK%FHNDL)

!
!     *** NUMS section
!
         READ(INP%FHNDL, NML = NUMS)
         wwm_print_namelist(NUMS)
         FLUSH(CHK%FHNDL)
         CALL READ_HISTORY_STATION_NAMELIST()
         IF (ICOMP .eq. 3) THEN
           IF (DMETHOD .GT. 0) THEN
             REFRACTION_IMPL=.TRUE.
           ELSE
             REFRACTION_IMPL=.FALSE.
           END IF
           IF (FMETHOD .GT. 0) THEN
             FREQ_SHIFT_IMPL=.TRUE.
           ELSE
             FREQ_SHIFT_IMPL=.FALSE.
           END IF
           IF (SMETHOD .GT. 0) THEN
             SOURCE_IMPL=.TRUE.
           ELSE
             SOURCE_IMPL=.FALSE.
           END IF
         ELSE
           REFRACTION_IMPL=.FALSE.
           FREQ_SHIFT_IMPL=.FALSE.
         END IF
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
         FLUSH(CHK%FHNDL)

         MULTIPLEIN_HOT=MULTIPLEIN
         MULTIPLEOUT_HOT=MULTIPLEOUT
         IF (DELTC.lt.MAIN%DELT) THEN
           DELTC=MAIN%DELT
         END IF
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
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SINGLE_READ_SPATIAL_GRID_TOTAL
      USE DATAPOOL
#ifdef NCDF
      USE NETCDF
#endif
      IMPLICIT NONE
      INTEGER :: I, IP, IE, ITMP, JTMP
      REAL(rkind)  :: XPDTMP, YPDTMP, ZPDTMP
      REAL(rkind) DXP1, DXP2, DXP3, DYP1, DYP2, DYP3
      INTEGER KTMP, LTMP, MTMP, NTMP, OTMP
      CHARACTER(LEN=100) :: RHEADER
#ifdef NCDF
      INTEGER :: ncid, dimidsB(2), dimidsA(1)
      character (len=20) :: MNEstr, MNPstr
      INTEGER var_id1, var_id2, var_id
      character (len = *), parameter :: CallFct="SINGLE_READ_SPATIAL_GRID_TOTAL"
#endif
      CALL TEST_FILE_EXIST_DIE('Missing grid file : ', GRD%FNAME)
      SELECT CASE (DIMMODE)
        CASE (1)
          OPEN(GRD%FHNDL, FILE = GRD%FNAME, STATUS = 'OLD')
          allocate(XPtotal(np_total), DEPtotal(np_total), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_input, allocate error 3')
          DO IP = 1, NP_TOTAL
            READ(GRD%FHNDL, *, IOSTAT = ISTAT) XPtotal(IP), DEPtotal(IP)
            IF ( ISTAT /= 0 ) CALL WWM_ABORT('error in the grid configuration file')
          END DO
          IF (LVAR1D) THEN
            DX1total(0)     = XPtotal(2)- XPtotal(1)
            DX1total(1)     = DX1total(0)
            DX1total(MNP)   = XPtotal(MNP) - XPtotal(MNP-1)
            DX1total(MNP+1) = DX1total(MNP)
            DX2total(0)     = DX1total(0)
            DX2total(MNP+1) = DX1total(MNP)
            DO IP = 2, NP_TOTAL-1 ! Bandwith at gridpoints
              DX1total(IP) = (XPtotal(IP)-XPtotal(IP-1))/2. + (XPtotal(IP+1)-XPtotal(IP))/2.
            END DO
            DO IP = 2, NP_TOTAL ! Stepwidth between gridpoints K and K-1
              DX2total(IP) = XPtotal(IP) - XPtotal(IP-1)
            END DO
            DX2total(1) = DX1total(0)
          END IF
          CLOSE(GRD%FHNDL)
        CASE (2)
          IF (IGRIDTYPE == 1) THEN ! system.dat format ... XFN
            OPEN(GRD%FHNDL, FILE = GRD%FNAME, STATUS = 'OLD')
            DO I = 1, 2
              READ(GRD%FHNDL, '(A)') RHEADER
            END DO
            READ(GRD%FHNDL, *, IOSTAT = ISTAT) ITMP
            IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=1 error in read mnp/mne')
            READ(GRD%FHNDL, '(A)') RHEADER
            READ(GRD%FHNDL, *, IOSTAT = ISTAT) JTMP 
            IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=1 error in read mnp/mne')
            NP_TOTAL = ITMP + JTMP
            allocate(XPtotal(np_total), YPtotal(np_total), DEPtotal(np_total), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_input, allocate error 4')
            DO I = 1, 7
              READ(GRD%FHNDL, '(A)') RHEADER
            END DO
            DO IP=1,NP_TOTAL
              READ(GRD%FHNDL, *, IOSTAT = ISTAT) KTMP, XPtotal(IP), YPtotal(IP), DEPtotal(IP)
              IF (KTMP+1.ne.IP) THEN
                CALL WWM_ABORT('IGRIDTYPE=1 error in grid read 3')
              ENDIF
              IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=1 error in grid read 4')
            END DO
            DO I = 1, 2
              READ(GRD%FHNDL, '(A)') RHEADER
            END DO
            READ(GRD%FHNDL, *, IOSTAT = ISTAT) NE_TOTAL
            allocate(INEtotal(3, ne_total), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_input, allocate error 5')
            DO I = 1, 3
              READ(GRD%FHNDL, '(A)') RHEADER
            END DO
            DO IE=1,NE_TOTAL
              READ(GRD%FHNDL, *, IOSTAT = ISTAT) KTMP, LTMP, MTMP, NTMP, OTMP
              IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=1 error in grid read 5')
              INEtotal(1,IE)=KTMP+1
              INEtotal(2,IE)=LTMP+1
              INEtotal(3,IE)=MTMP+1
            END DO
            CLOSE(GRD%FHNDL)
          ELSE IF (IGRIDTYPE == 2) THEN ! periodic grid written by mathieu dutour
            OPEN(GRD%FHNDL, FILE = GRD%FNAME, STATUS = 'OLD')
            READ(GRD%FHNDL,*) NE_TOTAL, NP_TOTAL
            allocate(DEPtotal(NP_TOTAL), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_input, allocate error 6')
            DO IP = 1, NP_TOTAL
              READ(GRD%FHNDL, *, IOSTAT = ISTAT) DEPtotal(IP)
              IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=2 error in grid read 1')
            END DO
            allocate(TRIAtotal(NE_TOTAL), INEtotal(3,ne_total), IENtotal(6,ne_total), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_input, allocate error 7')
            DO IE = 1, NE_TOTAL
              READ(GRD%FHNDL, *, IOSTAT = ISTAT) TRIAtotal(IE)
              IF ( ISTAT /= 0 )  CALL WWM_ABORT('IGRIDTYPE=2 error in grid read 2')
              READ(GRD%FHNDL, *, IOSTAT = ISTAT) INEtotal(:,IE)
              IF ( ISTAT /= 0 )  CALL WWM_ABORT('IGRIDTYPE=2 error in grid read 3')
              READ(GRD%FHNDL, *, IOSTAT = ISTAT) DXP1, DXP2, DXP3
              IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=2 error in grid read 4')
              READ(GRD%FHNDL, *, IOSTAT = ISTAT) DYP1, DYP2, DYP3
              IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=2 error in grid read 5')
              IENtotal(1,IE) = -DYP2
              IENtotal(2,IE) = DXP2
              IENtotal(3,IE) = -DYP3
              IENtotal(4,IE) = DXP3
              IENtotal(5,IE) = -DYP1
              IENtotal(6,IE) = DXP1
            END DO
            CLOSE(GRD%FHNDL)
          ELSE IF (IGRIDTYPE == 3) THEN ! selfe gr3
            OPEN(GRD%FHNDL, FILE = GRD%FNAME, STATUS = 'OLD')
            READ(GRD%FHNDL,*)
            READ(GRD%FHNDL,*, IOSTAT = ISTAT) NE_TOTAL, NP_TOTAL
            IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=3 error in read mnp/mne')
            allocate(XPtotal(np_total), YPtotal(np_total), DEPtotal(np_total), INEtotal(3, ne_total), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_input, allocate error 8')
            DO IP=1,NP_TOTAL
              READ(GRD%FHNDL, *, IOSTAT = ISTAT) KTMP, XPDTMP, YPDTMP, ZPDTMP
              XPtotal(IP)  = XPDTMP
              YPtotal(IP)  = YPDTMP
              DEPtotal(IP) = ZPDTMP
              IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=3 error in grid reading 1')
            END DO
            DO IE = 1, NE_TOTAL
              READ(GRD%FHNDL, *, IOSTAT = ISTAT) KTMP, LTMP, INEtotal(:,IE)
              IF ( ISTAT /= 0 )  CALL WWM_ABORT('IGRIDTYPE=3 error in grid reading 2')
            END DO
            CLOSE(GRD%FHNDL)
          ELSE IF (IGRIDTYPE == 4) THEN ! Old WWM format
            OPEN(GRD%FHNDL, FILE = GRD%FNAME, STATUS = 'OLD')
            READ(GRD%FHNDL, *, IOSTAT = ISTAT) NE_TOTAL 
            READ(GRD%FHNDL, *, IOSTAT = ISTAT) NP_TOTAL 
            allocate(XPtotal(np_total), YPtotal(np_total), DEPtotal(np_total), INEtotal(3, ne_total), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('allocate error 9')
            DO IP=1,NP_TOTAL
              READ(GRD%FHNDL, *, IOSTAT = ISTAT) XPtotal(IP), YPtotal(IP), DEPtotal(IP)
              IF ( ISTAT /= 0 ) CALL WWM_ABORT('IGRIDTYPE=4 error in grid read 1')
            END DO
            DO IE=1,NE_TOTAL
              READ(GRD%FHNDL, *, IOSTAT = ISTAT) INEtotal(:,IE)
              IF ( ISTAT /= 0 )  CALL WWM_ABORT('IGRIDTYPE=4 error in grid read 2')
            END DO
            CLOSE(GRD%FHNDL)
#ifdef NCDF
          ELSE IF (IGRIDTYPE == 5) THEN ! Netcdf format
            ISTAT = NF90_OPEN(GRD%FNAME, NF90_NOWRITE, ncid)
            CALL GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)

            ISTAT = nf90_inq_varid(ncid, 'ele', var_id)
            CALL GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)

            ISTAT = nf90_inquire_variable(ncid, var_id, dimids=dimidsB)
            CALL GENERIC_NETCDF_ERROR(CallFct, 3, ISTAT)

            ISTAT = nf90_inquire_dimension(ncid, dimidsB(2), name=MNEstr, len=ne_total)
            CALL GENERIC_NETCDF_ERROR(CallFct, 4, ISTAT)
            WRITE(DBG%FHNDL,*) 'MNEstr=', TRIM(MNEstr)
            FLUSH(DBG%FHNDL)

            ISTAT = nf90_inq_varid(ncid, 'depth', var_id)
            CALL GENERIC_NETCDF_ERROR(CallFct, 5, ISTAT)

            ISTAT = nf90_inquire_variable(ncid, var_id, dimids=dimidsA)
            CALL GENERIC_NETCDF_ERROR(CallFct, 6, ISTAT)

            ISTAT = nf90_inquire_dimension(ncid, dimidsB(1), name=MNPstr, len=np_total)
            CALL GENERIC_NETCDF_ERROR(CallFct, 7, ISTAT)
            WRITE(DBG%FHNDL,*) 'MNPstr=', TRIM(MNPstr)
            FLUSH(DBG%FHNDL)

            allocate(XPtotal(np_total), YPtotal(np_total), DEPtotal(np_total), INEtotal(3, ne_total), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('allocate error 9')

            ISTAT = nf90_inq_varid(ncid, 'depth', var_id)
            CALL GENERIC_NETCDF_ERROR(CallFct, 8, ISTAT)

            ISTAT = nf90_get_var(ncid, var_id, DEPtotal)
            CALL GENERIC_NETCDF_ERROR(CallFct, 9, ISTAT)

            ISTAT = nf90_inq_varid(ncid, 'ele', var_id)
            CALL GENERIC_NETCDF_ERROR(CallFct, 10, ISTAT)

            ISTAT = nf90_get_var(ncid, var_id, INEtotal)
            CALL GENERIC_NETCDF_ERROR(CallFct, 11, ISTAT)

            IF (LSPHE) THEN
              ISTAT = nf90_inq_varid(ncid, 'lon', var_id1)
              CALL GENERIC_NETCDF_ERROR(CallFct, 14, ISTAT)

              ISTAT = nf90_inq_varid(ncid, 'lat', var_id2)
              CALL GENERIC_NETCDF_ERROR(CallFct, 15, ISTAT)
            ELSE
              ISTAT = nf90_inq_varid(ncid, 'x', var_id1)
              CALL GENERIC_NETCDF_ERROR(CallFct, 12, ISTAT)

              ISTAT = nf90_inq_varid(ncid, 'y', var_id2)
              CALL GENERIC_NETCDF_ERROR(CallFct, 13, ISTAT)
            END IF
            ISTAT = nf90_get_var(ncid, var_id1, XPtotal)
            CALL GENERIC_NETCDF_ERROR(CallFct, 16, ISTAT)

            ISTAT = nf90_get_var(ncid, var_id2, YPtotal)
            CALL GENERIC_NETCDF_ERROR(CallFct, 17, ISTAT)

            ISTAT = NF90_CLOSE(ncid)
            CALL GENERIC_NETCDF_ERROR(CallFct, 18, ISTAT)
#endif
          ELSE
            CALL WWM_ABORT('IGRIDTYPE WRONG')
          END IF
        CASE DEFAULT
          CALL WWM_ABORT('WRONG GRID DIMENSION')
      END SELECT
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_SPATIAL_GRID_TOTAL
      USE DATAPOOL
      IMPLICIT NONE
      integer :: rbuf_int(2)
      real(rkind), allocatable :: rbuf_real(:)
      integer iProc, IP, IE, nb_real, idx
#ifdef MPI_PARALL_GRID
      IF (MULTIPLE_IN_GRID) THEN
        CALL SINGLE_READ_SPATIAL_GRID_TOTAL
      ELSE
        IF (DIMMODE .ne. 2) THEN
          CALL WWM_ABORT('Parallel mode only for 2D')
        ENDIF
        IF (myrank .eq. 0) THEN
          CALL SINGLE_READ_SPATIAL_GRID_TOTAL
          rbuf_int(1)=np_total
          rbuf_int(2)=ne_total
          DO iProc=2,nproc
            CALL MPI_SEND(rbuf_int,2,itype, iProc-1, 30, comm, ierr)
          END DO
          DO iProc=2,nproc
            CALL MPI_SEND(INEtotal,3*ne_total,itype, iProc-1, 32, comm, ierr)
          END DO
          IF (IGRIDTYPE .eq. 2) THEN
            nb_real=np_total + 7*ne_total
            allocate(rbuf_real(nb_real), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('allocate error 10')
            idx=0
            DO IP=1,NP_TOTAL
              idx=idx+1
              rbuf_real(idx)=DEPtotal(IP)
            END DO
            DO IE=1,NE_TOTAL
              rbuf_real(idx+1)=TRIAtotal(IE)
              rbuf_real(idx+2)=IENtotal(1,IE)
              rbuf_real(idx+3)=IENtotal(2,IE)
              rbuf_real(idx+4)=IENtotal(3,IE)
              rbuf_real(idx+5)=IENtotal(4,IE)
              rbuf_real(idx+6)=IENtotal(5,IE)
              rbuf_real(idx+7)=IENtotal(6,IE)
              idx=idx+7
            END DO
          ELSE
            nb_real=3*np_total
            allocate(rbuf_real(nb_real), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('allocate error 11')
            idx=0
            DO IP=1,NP_TOTAL
              rbuf_real(idx+1)=XPtotal(IP)
              rbuf_real(idx+2)=YPtotal(IP)
              rbuf_real(idx+3)=DEPtotal(IP)
              idx=idx+3
            END DO
          END IF
          DO iProc=2,nproc
            CALL MPI_SEND(rbuf_real,nb_real,rtype, iProc-1, 34, comm, ierr)
          END DO
          deallocate(rbuf_real)
        ELSE
          CALL MPI_RECV(rbuf_int,2,itype, 0, 30, comm, istatus, ierr)
          np_total=rbuf_int(1)
          ne_total=rbuf_int(2)
          allocate(INEtotal(3,ne_total), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('allocate error 12')
          CALL MPI_RECV(INEtotal,3*ne_total,itype, 0, 32, comm, istatus, ierr)
          IF (IGRIDTYPE .eq. 2) THEN
            nb_real=np_total + 7*ne_total
            allocate(rbuf_real(nb_real), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('allocate error 13')
            CALL MPI_RECV(rbuf_real,nb_real,rtype, 0, 34, comm, istatus, ierr)
            allocate(DEPtotal(np_total), TRIAtotal(ne_total), IENtotal(6,ne_total), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('allocate error 14')
            idx=0
            DO IP=1,NP_TOTAL
              idx=idx+1
              DEPtotal(IP)=rbuf_real(idx)
            END DO
            DO IE=1,NE_TOTAL
              TRIAtotal(IE)=rbuf_real(idx+1)
              IENtotal(1,IE)=rbuf_real(idx+2)
              IENtotal(2,IE)=rbuf_real(idx+3)
              IENtotal(3,IE)=rbuf_real(idx+4)
              IENtotal(4,IE)=rbuf_real(idx+5)
              IENtotal(5,IE)=rbuf_real(idx+6)
              IENtotal(6,IE)=rbuf_real(idx+7)
              idx=idx+7
            END DO
          ELSE
            nb_real=3*np_total
            allocate(rbuf_real(nb_real), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('allocate error 15')
            allocate(DEPtotal(np_total), XPtotal(np_total), YPtotal(np_total), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('allocate error 16')
            CALL MPI_RECV(rbuf_real,nb_real,rtype, 0, 34, comm, istatus, ierr)
            idx=0
            DO IP=1,NP_TOTAL
              XPtotal(IP) = rbuf_real(idx+1)
              YPtotal(IP) = rbuf_real(idx+2)
              DEPtotal(IP)= rbuf_real(idx+3)
              idx=idx+3
            END DO
          END IF
          deallocate(rbuf_real)
        END IF
      END IF
#else
      CALL SINGLE_READ_SPATIAL_GRID_TOTAL
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CHECK_LOGICS()
         USE DATAPOOL
         IMPLICIT NONE

         REAL(rkind) :: TEST

!        Check timings ...

         IF (OUT_HISTORY%BMJD .GE. OUT_HISTORY%EMJD) THEN
           WRITE(STAT%FHNDL,*) 'MAIN%BEGT=',MAIN%BEGT
           WRITE(STAT%FHNDL,*) 'MAIN%ENDT=',MAIN%ENDT
           WRITE(STAT%FHNDL,*) 'MAIN%BMJD=',MAIN%BMJD
           WRITE(STAT%FHNDL,*) 'MAIN%EMJD=',MAIN%EMJD

           WRITE(STAT%FHNDL,*) 'OUT_HISTORY%BEGT=',OUT_HISTORY%BEGT
           WRITE(STAT%FHNDL,*) 'OUT_HISTORY%ENDT=',OUT_HISTORY%ENDT
           WRITE(STAT%FHNDL,*) 'OUT_HISTORY%BMJD=',OUT_HISTORY%BMJD
           WRITE(STAT%FHNDL,*) 'OUT_HISTORY%EMJD=',OUT_HISTORY%EMJD

!           Print *, 'OUT_HISTORY%BMJD=', OUT_HISTORY%BMJD
!           Print *, 'OUT_HISTORY%EMJD=', OUT_HISTORY%EMJD
           CALL WWM_ABORT('CHECK OUTPUT HISTORY TIME STEPS BEGINN TIME STEP IS SMALLER THAN END TIME STEP')
         END IF
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
         
#ifdef MPI_PARALL_GRID
         IF (ICOMP .GT. 0) THEN
           IF ((AMETHOD .eq. 1).or.(AMETHOD .eq. 2).or.(AMETHOD .eq. 3)) THEN
             IF (myrank .gt. 0) CALL WWM_ABORT('The AMETHOD = 1, 2, 3 are not parallelized')
           END IF
         END IF
#endif
         IF (ICOMP .eq. 3) THEN
#if !defined PETSC || !defined MPI_PARALL_GRID
           IF (AMETHOD .eq. 5) THEN
             CALL WWM_ABORT('For ICOMP=3 and AMETHOD=5 we need PETSC')
           END IF
#endif
#ifndef WWM_SOLVER
           IF (AMETHOD .eq. 7) THEN
             CALL WWM_ABORT('For ICOMP=3 and AMETHOD=7 we need WWM_SOLVER')
           END IF
#endif
           IF ((AMETHOD .ne. 5).and.(AMETHOD .ne. 7).and.(AMETHOD .ne. 0)) THEN
             CALL WWM_ABORT('We need AMETHOD=5 or 7 for ICOMP=3')
           END IF
         END IF
!
!        Check MSC,MDC for exchange
!
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

         IF (AMETHOD .eq. 7) THEN
           IF (LNONL .AND. (.NOT. SOURCE_IMPL)) THEN
             CALL WWM_ABORT('SOURCE_IMPL=F and LNONL=T is absurd')
           END IF
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
           call wwm_abort('ICOMP must be greater equal 1 to use PETSc')
         END IF

         ! if using PETSc with AMETHOD 4 or 5 LVECTOR must be FALSE
         IF (AMETHOD .GE. 4 .AND. LVECTOR .EQV. .TRUE.) THEN
           call wwm_abort('LVECTOR must be FALSE to use PETSc')
         END IF

         IF (LSOURCESWAM .AND. MELIM .NE. 3) THEN
           call wwm_abort('FOR WAM YOU NEED MELIM == 3')
         ELSE IF (.NOT. LSOURCESWAM .AND. MELIM .EQ. 3) THEN
           call wwm_abort('FOR WWM SOURCES YOU NEED MELIM .LT. 3') 
         ENDIF

#ifndef GRB
         IF (IWINDFORMAT == 7) CALL wwm_abort('you need to compile with grib')
#endif

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_CURRENT_INPUT
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER :: IP
#ifdef MPI_PARALL_GRID
      INTEGER :: I
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
        CALL CSEVAL( CUR%FHNDL, TRIM(CUR%FNAME), LCURFILE, 2, CURTXY, MULTIPLE_IN_CURR)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_WATLEV_INPUT
      USE DATAPOOL
      IMPLICIT NONE

#ifdef MPI_PARALL_GRID
      INTEGER :: I
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
        CALL CSEVAL( WAT%FHNDL,TRIM(WAT%FNAME), LWATLFILE, 1, WATLEV, MULTIPLE_IN_WATLEV)
      END IF
      END SUBROUTINE
