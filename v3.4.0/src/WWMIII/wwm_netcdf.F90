#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CF_EXTRACT_TIME(eStrUnitTime, ConvertToDay, eTimeStart)
      USE DATAPOOL
      IMPLICIT NONE
      character(len=100), intent(in) :: eStrUnitTime
      real(rkind), intent(out) :: ConvertToDay, eTimeStart
      character (len=100) :: Xname, Yname
      character (len=10) :: YnameYear, YnameMonth, YnameDay
      character (len=10) :: YnameHour, YnameMin, YnameSec
      character (len=50) :: YnameB, YnameD, YnameE
      character (len=50) :: YnameDate, YnameTime, YnameTimeP
      character (len=15) :: eStrTime
      integer alenB, alenC, alenD, alenE, alenTime, alenDate
      integer alen, posBlank
      integer lenHour, lenMin, lenSec, lenMonth, lenDay, posSepDateTime
      alen=LEN_TRIM(eStrUnitTime)
      posBlank=INDEX(eStrUnitTime(1:alen), ' ')
      Xname=eStrUnitTime(1:posBlank-1) ! should be days/hours/seconds
      IF (TRIM(Xname) .eq. 'days') THEN
        ConvertToDay=1
      ELSEIF (TRIM(Xname) .eq. 'hours') THEN
        ConvertToDay=1/24
      ELSEIF (TRIM(Xname) .eq. 'seconds') THEN
        ConvertToDay=1/86400
      ELSE
        CALL WWM_ABORT('Error in the code for conversion')
      END IF
      !
      Yname=eStrUnitTime(posBlank+1:alen)
      alenB=LEN_TRIM(Yname)
      posBlank=INDEX(Yname(1:alenB), ' ')
      YnameB=Yname(posBlank+1:alenB) ! should be 1990-01-01 0:0:0
      !
      alenC=LEN_TRIM(YnameB)
      posSepDateTime=INDEX(YnameB(1:alenC), ' ')
      IF (posSepDateTime .gt. 0) THEN
        YnameDate=YnameB(1:posSepDateTime-1) ! should be 1990-01-01
        YnameTimeP=YnameB(posSepDateTime+1:alenC) ! should be 0:0:0
        alenC=LEN_TRIM(YnameTimeP)
        posBlank=INDEX(YnameTimeP(1:alenC), ' ')
        IF (posBlank .eq. 0) THEN
          YnameTime=YnameTimeP
        ELSE
          YnameTime=YnameTimeP(1:posBlank-1)
        END IF
      ELSE
        YnameDate=YnameB
        eStrTime(10:10)='0'
        eStrTime(11:11)='0'
        eStrTime(12:12)='0'
        eStrTime(13:13)='0'
        eStrTime(14:14)='0'
        eStrTime(15:15)='0'
      END IF
      !
      alenDate=LEN_TRIM(YnameDate)
      posBlank=INDEX(YnameDate(1:alenDate), '-')
      YnameYear=YnameDate(1:posBlank-1) ! should be 1990
      YnameD=YnameDate(posBlank+1:alenDate)
      alenD=LEN_TRIM(YnameD)
      posBlank=INDEX(YnameD(1:alenD), '-')
      YnameMonth=YnameD(1:posBlank-1) ! should be 01
      YnameDay=YnameD(posBlank+1:alenD) ! should be 01
      !
      ! year
      eStrTime( 1: 1)=YnameYear( 1: 1)
      eStrTime( 2: 2)=YnameYear( 2: 2)
      eStrTime( 3: 3)=YnameYear( 3: 3)
      eStrTime( 4: 4)=YnameYear( 4: 4)
      !
      ! month
      WRITE(WINDBG%FHNDL,*) 'YnameMonth=', YnameMonth
      lenMonth=LEN_TRIM(YnameMonth)
      IF (lenMonth .eq. 2) THEN
        eStrTime( 5: 5)=YnameMonth( 1: 1)
        eStrTime( 6: 6)=YnameMonth( 2: 2)
      ELSE
        IF (lenMonth .eq. 1) THEN
          eStrTime( 5: 5)='0'
          eStrTime( 6: 6)=YnameMonth( 1: 1)
        ELSE
          CALL WWM_ABORT('DIE in trying to get the month')
        END IF
      END IF
      !
      ! day
      lenDay=LEN_TRIM(YnameDay)
      IF (lenDay .eq. 2) THEN
        eStrTime( 7: 7)=YnameDay( 1: 1)
        eStrTime( 8: 8)=YnameDay( 2: 2)
      ELSE
        IF (lenDay .eq. 1) THEN
          eStrTime( 7: 7)='0'
          eStrTime( 8: 8)=YnameDay( 1: 1)
        ELSE
          CALL WWM_ABORT('DIE in trying to get the day')
        END IF
      END IF
      !
      eStrTime( 9: 9)='.'
      !
      IF (posSepDateTime .gt. 0) THEN
        !
        alenTime=LEN_TRIM(YnameTime)
        posBlank=INDEX(YnameTime(1:alenTime), ':')
        YnameHour=YnameTime(1:posBlank-1) ! should be 0
        YnameE=YnameTime(posBlank+1:alenTime)
        alenE=LEN_TRIM(YnameE)
        posBlank=INDEX(YnameE(1:alenE), ':')
        YnameMin=YnameE(1:posBlank-1) ! should be 0
        YnameSec=YnameE(posBlank+1:alenE) ! should be 0
        !
        !
        ! Hour
        lenHour=LEN_TRIM(YnameHour)
        IF (lenHour .eq. 2) THEN
          eStrTime(10:10)=YnameHour( 1: 1)
          eStrTime(11:11)=YnameHour( 2: 2)
        ELSE
          IF (lenHour .eq. 1) THEN
            eStrTime(10:10)='0'
            eStrTime(11:11)=YnameHour( 1: 1)
          ELSE
            CALL WWM_ABORT('DIE in trying to get the hour')
          END IF
        END IF
        !
        ! Min
        lenMin=LEN_TRIM(YnameMin)
        IF (lenMin .eq. 2) THEN
          eStrTime(12:12)=YnameMin( 1: 1)
          eStrTime(13:13)=YnameMin( 2: 2)
        ELSE
          IF (lenMin .eq. 1) THEN
            eStrTime(12:12)='0'
            eStrTime(13:13)=YnameMin( 1: 1)
          ELSE
            CALL WWM_ABORT('DIE in trying to get the min')
          END IF
        END IF
        !
        ! Sec
        lenSec=LEN_TRIM(YnameSec)
        IF (lenSec .eq. 2) THEN
          eStrTime(14:14)=YnameSec( 1: 1)
          eStrTime(15:15)=YnameSec( 2: 2)
        ELSE
          IF (lenSec .eq. 1) THEN
            eStrTime(14:14)='0'
            eStrTime(15:15)=YnameSec( 1: 1)
          ELSE
            WRITE(WINDBG%FHNDL,*) 'YnameSec=', TRIM(Ynamesec)
            WRITE(WINDBG%FHNDL,*) 'lenSec=', lenSec
            FLUSH(WINDBG%FHNDL)
            CALL WWM_ABORT('DIE in trying to get the sec')
          END IF
        END IF
      END IF
      WRITE(WINDBG%FHNDL,*) 'eStrTime=', eStrTime
      CALL CT2MJD(eStrTime, eTimeStart)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DETERMINE_NEEDED_COMPUTATION(eVar)
      USE DATAPOOL, only : VAROUT, OUTVARS_COMPLETE
      implicit none
      type(VAROUT), intent(inout) :: eVar
      LOGICAL     ::   HS, TM01, TM02, TM10, KLM, WLM,                  &
     &   ETOTC, ETOTS, DM, DSPR,                                        &
     &   TPPD, CPPD, KPPD, CGPD,                                        &
     &   TPP, CPP, WNPP, CGPP, KPP, LPP, PEAKD, PEAKDSPR,               &
     &   DPEAK, UBOT, ORBITAL, BOTEXPER, TMBOT,                         &
     &   URSELL, UFRIC, Z0, ALPHA_CH, WINDX, WINDY, CD,                 &
     &   CURRTX, CURRTY, WATLEV, WATLEVOLD, DEPDT, DEP,                 &
     &   WINDMAG, TAUW, TAUWX, TAUWY, TAUHF, TAUTOT,                    &
     &   STOKESBOTTX, STOKESBOTTY,                                      &
     &   STOKESSURFX, STOKESSURFY, STOKESBAROX, STOKESBAROY,            &
     &   RSXX, RSXY, RSYY, CFL1, CFL2, CFL3, ZETA_SETUP
      LOGICAL :: ComputeMean, ComputeDirSpread, ComputePeak
      LOGICAL :: ComputeCurr, ComputeUrsell, ComputeStokes
      integer iVar, idx, nbOutVarEff
      integer istat
      HS           = eVar%LVAR( 1)
      TM01         = eVar%LVAR( 2)
      TM02         = eVar%LVAR( 3)
      TM10         = eVar%LVAR( 4)
      KLM          = eVar%LVAR( 5)
      WLM          = eVar%LVAR( 6)
      ETOTC        = eVar%LVAR( 7)
      ETOTS        = eVar%LVAR( 8)
      DM           = eVar%LVAR( 9)
      DSPR         = eVar%LVAR(10)
      TPPD         = eVar%LVAR(11)
      CPPD         = eVar%LVAR(12)
      KPPD         = eVar%LVAR(13)
      CGPD         = eVar%LVAR(14)
      TPP          = eVar%LVAR(15)
      CPP          = eVar%LVAR(16)
      WNPP         = eVar%LVAR(17)
      CGPP         = eVar%LVAR(18)
      KPP          = eVar%LVAR(19)
      LPP          = eVar%LVAR(20)
      PEAKD        = eVar%LVAR(21)
      PEAKDSPR     = eVar%LVAR(22)
      DPEAK        = eVar%LVAR(23)
      UBOT         = eVar%LVAR(24)
      ORBITAL      = eVar%LVAR(25)
      BOTEXPER     = eVar%LVAR(26)
      TMBOT        = eVar%LVAR(27)
      URSELL       = eVar%LVAR(28)
      UFRIC        = eVar%LVAR(29)
      Z0           = eVar%LVAR(30)
      ALPHA_CH     = eVar%LVAR(31)
      WINDX        = eVar%LVAR(32)
      WINDY        = eVar%LVAR(33)
      CD           = eVar%LVAR(34)
      CURRTX       = eVar%LVAR(35)
      CURRTY       = eVar%LVAR(36)
      WATLEV       = eVar%LVAR(37)
      WATLEVOLD    = eVar%LVAR(38)
      DEPDT        = eVar%LVAR(39)
      DEP          = eVar%LVAR(40)
      WINDMAG      = eVar%LVAR(41)
      TAUW         = eVar%LVAR(42)
      TAUWX        = eVar%LVAR(43)
      TAUWY        = eVar%LVAR(44)
      TAUHF        = eVar%LVAR(45)
      TAUTOT       = eVar%LVAR(46)
      STOKESBOTTX  = eVar%LVAR(47)
      STOKESBOTTY  = eVar%LVAR(48)
      STOKESSURFX  = eVar%LVAR(49)
      STOKESSURFY  = eVar%LVAR(50)
      STOKESBAROX  = eVar%LVAR(51)
      STOKESBAROY  = eVar%LVAR(52)
      RSXX         = eVar%LVAR(53)
      RSXY         = eVar%LVAR(54)
      RSYY         = eVar%LVAR(55)
      CFL1         = eVar%LVAR(56)
      CFL2         = eVar%LVAR(57)
      CFL3         = eVar%LVAR(58)
      ZETA_SETUP   = eVar%LVAR(59)
      ComputeMean=.FALSE.
      ComputeDirSpread=.FALSE.
      ComputePeak=.FALSE.
      ComputeCurr=.FALSE.
      ComputeUrsell=.FALSE.
      ComputeStokes=.FALSE.
      IF (HS .or. TM01 .or. TM02 .or. TM10 .or. KLM .or. WLM) THEN
        ComputeMean=.TRUE.
      END IF
      IF (ETOTC .or. ETOTS .or. DM .or. DSPR) THEN
        ComputeDirspread=.TRUE.
      END IF
      IF (TPPD .or. CPPD .or. KPPD .or. CGPD .or. TPP .or. CPP .or. WNPP .or. CGPP .or. KPP .or. LPP .or. PEAKD .or. PEAKDSPR .or. DPEAK) THEN
        ComputePeak=.TRUE.
      END IF
      IF (UBOT .or. ORBITAL .or. BOTEXPER .or. TMBOT) THEN
        ComputeCurr=.TRUE.
      END IF
      IF (URSELL) THEN
        ComputeUrsell=.TRUE.
        ComputeMean=.TRUE.
        ComputePeak=.TRUE.
      END IF
      IF (STOKESSURFX .or. STOKESSURFY .or. STOKESBOTTX .or. STOKESBOTTY .or. STOKESBAROX .or. STOKESBAROY) THEN
        ComputeStokes=.TRUE.
      END IF
      eVar%ComputeMean=ComputeMean
      eVar%ComputeDirspread=ComputeDirspread
      eVar%ComputePeak=ComputePeak
      eVar%ComputeCurr=ComputeCurr
      eVar%ComputeUrsell=ComputeUrsell
      eVar%ComputeStokes=ComputeStokes
      nbOutVarEff=0
      DO iVar=1,OUTVARS_COMPLETE
        IF (eVar%LVAR(iVar)) THEN
          nbOutVarEff=nbOutVarEff+1
        END IF
      END DO
      eVar%nbOutVarEff=nbOutVarEff
      allocate(eVar%ListIdxEff(nbOutVarEff), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_netcdf, allocate error 1')
      idx=0
      DO iVar=1,OUTVARS_COMPLETE
        IF (eVar%LVAR(iVar)) THEN
          idx=idx+1
          eVar%ListIdxEff(idx)=iVar
        END IF
      END DO
      END SUBROUTINE
#ifdef NCDF
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SERIAL_GET_BOUNDARY(np_glob, INEglob, ne_glob, IOBP, NEIGHBOR)
      USE DATAPOOL, ONLY : DBG
      IMPLICIT NONE
      INTEGER, INTENT(IN)             :: np_glob, ne_glob
      INTEGER, INTENT(IN)             :: INEglob(3, ne_glob)
      INTEGER, INTENT(INOUT)          :: IOBP(np_glob)
      INTEGER, INTENT(INOUT)          :: NEIGHBOR(np_glob)

      INTEGER :: STATUS(np_glob), COLLECTED(np_glob)
      INTEGER :: NEXTVERT(np_glob), PREVVERT(np_glob), NBneighbor(np_glob)

      INTEGER :: IE, I, IP, eIdx
      INTEGER :: ISFINISHED, INEXT, IPREV, IPNEXT, IPPREV, ZNEXT
      LOGICAL :: HaveError
      NEIGHBOR=0
      STATUS = 0
      NEXTVERT = 0
      PREVVERT = 0
      DO IE=1,ne_glob
        DO I=1,3
          IF (I.EQ.1) THEN
            IPREV=3
          ELSE
            IPREV=I-1
          END IF
          IF (I.EQ.3) THEN
            INEXT=1
          ELSE
            INEXT=I+1
          END IF
          IP=INEglob(I,IE)
          IPNEXT=INEglob(INEXT,IE)
          IPPREV=INEglob(IPREV,IE)
          IF (STATUS(IP).EQ.0) THEN
            STATUS(IP)=1
            PREVVERT(IP)=IPPREV
            NEXTVERT(IP)=IPNEXT
          END IF
        END DO
      END DO
      STATUS=0
      DO
        COLLECTED=0
        DO IE=1,ne_glob
          DO I=1,3
            IF (I.EQ.1) THEN
              IPREV=3
            ELSE
              IPREV=I-1
            END IF
            IF (I.EQ.3) THEN
              INEXT=1
            ELSE
              INEXT=I+1
            END IF
            IP=INEglob(I,IE)
            IPNEXT=INEglob(INEXT,IE)
            IPPREV=INEglob(IPREV,IE)
            IF (STATUS(IP).eq.0) THEN
              ZNEXT=NEXTVERT(IP)
              IF (ZNEXT.eq.IPPREV) THEN
                COLLECTED(IP)=1
                NEXTVERT(IP)=IPNEXT
                IF (NEXTVERT(IP).eq.PREVVERT(IP)) THEN
                  STATUS(IP)=1
                END IF
              END IF
            END IF
          END DO
        END DO
        ISFINISHED=1
        DO IP=1,np_glob
          IF ((COLLECTED(IP).eq.0).and.(STATUS(IP).eq.0)) THEN
            STATUS(IP)=-1
            NEIGHBOR(IP)=NEXTVERT(IP)     ! new code
          END IF
          IF (STATUS(IP).eq.0) THEN
            ISFINISHED=0
          END IF
        END DO
        IF (ISFINISHED.eq.1) THEN
          EXIT
        END IF
      END DO
      IOBP=0
      DO IP=1,np_glob
        IF (STATUS(IP).eq.-1) THEN
          IOBP(IP)=1
        END IF
      END DO
      NBneighbor=0
      DO IP=1,np_glob
        eIdx=NEIGHBOR(IP)
        IF (eIdx.gt.0) THEN
          NBneighbor(eIdx)=NBneighbor(eIdx)+1
        END IF
      END DO
      HaveError=.FALSE.
      DO IP=1,np_glob
        IF (NBneighbor(IP).gt.1) THEN
          WRITE(DBG%FHNDL,*) 'Inconsistency in the output'
          WRITE(DBG%FHNDL,*) '  Vertex ', IP, ' is ', NBneighbor(IP), ' times neighbor'
          HaveError=.TRUE.
        END IF
        IF ((NBneighbor(IP).eq.1).and.(NEIGHBOR(IP).eq.0)) THEN
          WRITE(DBG%FHNDL,*) 'Inconsistency in the output'
          WRITE(DBG%FHNDL,*) '  Vertex ', IP, ' is a neighbor'
          WRITE(DBG%FHNDL,*) '  but has no neighbor!'
          HaveError=.TRUE.
        END IF
        IF ((NBneighbor(IP).eq.0).and.(NEIGHBOR(IP).gt.0)) THEN
          WRITE(DBG%FHNDL,*) 'Inconsistency in the output'
          WRITE(DBG%FHNDL,*) '  Vertex ', IP, ' has a neighbor'
          WRITE(DBG%FHNDL,*) '  but is not a neighbor!'
          HaveError=.TRUE.
        END IF
      END DO
      IF (HaveError) THEN
        WRITE(DBG%FHNDL,*) 'Find some errors in the output'
        WRITE(DBG%FHNDL,*) 'Please check for node contained in several boundaries'
        CALL WWM_ABORT('Please debug the boundary code')
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SERIAL_WRITE_BOUNDARY(ncid, np_glob, ne_glob, INEglob, Oper)
      USE NETCDF
      USE DATAPOOL
      implicit none
      integer, intent(in) :: ncid, Oper
      integer, intent(in) :: np_glob, ne_glob
      integer, intent(in) :: INEglob(ne_glob,3)
      integer :: STATUS(np_glob)
      integer :: IOBPglob(np_glob)
      integer :: NEIGHBOR(np_glob)
      integer LenBound, IP, NbCycle
      integer IPfirst, IPwork, iCycle, TheLength
      integer, allocatable :: ListSequence(:)
      integer, allocatable :: SequenceNumber(:)
      integer, allocatable :: LengthCycle(:)
      integer iret, lenbound_dims, nbbound_dims, var_id
      integer idx
      character (len = *), parameter :: CallFct="SERIAL_WRITE_BOUNDARY"
      character (len = *), parameter :: FULLNAME = "full-name"
      STATUS = 0
      CALL SERIAL_GET_BOUNDARY(np_glob, INEglob, ne_glob, IOBPglob, NEIGHBOR)

      LenBound=0
      DO IP=1,np_glob
        IF (IOBPglob(IP).eq.1) THEN
          LenBound=LenBound+1
          STATUS(IP)=1
        END IF
      END DO
      ALLOCATE(ListSequence(LenBound), SequenceNumber(LenBound), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_netcdf, allocate error 5')
      NbCycle=0
      DO IP=1,np_glob
        IF (STATUS(IP).eq.1) THEN
          NbCycle=NbCycle+1
          IPfirst=IP
          IPwork=IP
          DO
            IF (IPwork == 0) EXIT
            STATUS(IPwork)=0
            IPwork=NEIGHBOR(IPwork)
            IF (IPwork.eq.IPfirst) THEN
              EXIT
            END IF
          END DO
        END IF
      END DO
      DO IP=1,np_glob
        IF (IOBPglob(IP).eq.1) THEN
          STATUS(IP)=1
        END IF
      END DO
      allocate(LengthCycle(NbCycle), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_netcdf, allocate error 6')
      iCycle=0
      idx=0
      DO IP=1,np_glob
        IF (STATUS(IP).eq.1) THEN
          iCycle=iCycle+1
          IPfirst=IP
          IPwork=IP
          TheLength=0
          DO
            IF (IPwork == 0) EXIT
            STATUS(IPwork)=0
            TheLength=TheLength+1
            idx=idx+1
            ListSequence(idx)=IPwork
            SequenceNumber(idx)=iCycle
            IPwork=NEIGHBOR(IPwork)
            IF (IPwork.eq.IPfirst) THEN
              EXIT
            END IF
          END DO
          LengthCycle(iCycle)=TheLength
        END IF
      END DO
      !
      IF ((Oper == 1).and.(LenBound.gt.0)) THEN
        iret = nf90_def_dim(ncid, 'lenbound', LenBound, lenbound_dims)
        CALL GENERIC_NETCDF_ERROR(CallFct, 1, iret)
        iret = nf90_def_dim(ncid, 'nbbound', NbCycle, nbbound_dims)
        CALL GENERIC_NETCDF_ERROR(CallFct, 2, iret)
        !
        iret=nf90_def_var(ncid,'inode',NF90_INT,(/ lenbound_dims/),var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 3, iret)
        iret=nf90_put_att(ncid,var_id,FULLNAME,'IP of boundary element')
        CALL GENERIC_NETCDF_ERROR(CallFct, 4, iret)
        !
        iret=nf90_def_var(ncid,'icycle',NF90_INT,(/ lenbound_dims/),var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 5, iret)
        iret=nf90_put_att(ncid,var_id,FULLNAME,'index of corresponding cycle')
        CALL GENERIC_NETCDF_ERROR(CallFct, 6, iret)
        !
        iret=nf90_def_var(ncid,'lencycle',NF90_INT,(/ nbbound_dims/),var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 7, iret)
        iret=nf90_put_att(ncid,var_id,FULLNAME,'Length of cycles')
        CALL GENERIC_NETCDF_ERROR(CallFct, 8, iret)
      END IF
      IF ((Oper == 2).and.(LenBound.gt.0)) THEN
        iret=nf90_inq_varid(ncid, "inode", var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 9, iret)
        iret=nf90_put_var(ncid,var_id,ListSequence, start = (/1/), count = (/ LenBound/))
        CALL GENERIC_NETCDF_ERROR(CallFct, 10, iret)
        !
        iret=nf90_inq_varid(ncid, "icycle", var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 11, iret)
        iret=nf90_put_var(ncid,var_id,SequenceNumber, start = (/1/), count = (/ LenBound/))
        CALL GENERIC_NETCDF_ERROR(CallFct, 12, iret)
        !
        iret=nf90_inq_varid(ncid, "lencycle", var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 13, iret)
        iret=nf90_put_var(ncid,var_id,LengthCycle, start = (/1/), count = (/ NbCycle/))
        CALL GENERIC_NETCDF_ERROR(CallFct, 14, iret)
      END IF
      !
      deallocate(LengthCycle, ListSequence, SequenceNumber)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE REPORT_ERROR_INQ(iret, str)
      INTEGER, intent(in) :: iret
      character(*) :: str
      character(256) :: ErrMsg
      IF (iret.ne.0) THEN
        WRITE (ErrMsg,10) 'ERROR while inquiring', str
  10    FORMAT (a,' ',a)
        CALL WWM_ABORT(TRIM(ErrMsg))
      ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE REPORT_ERROR_DEF(iret, str)
      INTEGER, intent(in) :: iret
      character(*) :: str
      character(256) :: ErrMsg
      IF (iret.ne.0) THEN
        WRITE (ErrMsg,10) 'ERROR while defining', str
  10    FORMAT (a,' ',a)
        CALL WWM_ABORT(TRIM(ErrMsg))
      ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GENERIC_NETCDF_ERROR(CallFct, idx, iret)
      USE NETCDF
      USE DATAPOOL, only : wwmerr
      implicit none
      integer, intent(in) :: iret, idx
      character(*), intent(in) :: CallFct
      character(len=500) :: CHRERR
      IF (iret .NE. nf90_noerr) THEN
        CHRERR = nf90_strerror(iret)
        WRITE(wwmerr,*) TRIM(CallFct), ' -', idx, '-', CHRERR
        CALL WWM_ABORT(wwmerr)
      ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE NAMEVARIABLE(idx, eStr, eStrFullName, eStrUnit)
      INTEGER, intent(in) :: idx
      character(len=40), intent(out) :: eStr
      character(len=80), intent(out) :: eStrFullName
      character(len=40), intent(out) :: eStrUnit

      IF (IDX.eq.1) THEN
        eStr="HS"
        eStrFullName="Significant wave height"
        eStrUnit="meter"
      ELSE IF (IDX.eq.2) THEN
        eStr="TM01"
        eStrFullName="mean wave period"
        eStrUnit="second"
      ELSE IF (IDX.eq.3) THEN
        eStr="TM02"
        eStrFullName="Zero-crossing wave period"
        eStrUnit="second"
      ELSE IF (IDX.eq.4) THEN
        eStr="TM10"
        eStrFullName="Mean period of wave over topping/run-up"
        eStrUnit="second"
      ELSE IF (IDX.eq.5) THEN
        eStr="KLM"
        eStrFullName="mean wave number"
        eStrUnit="meter-1"
      ELSE IF (IDX.eq.6) THEN
        eStr="WLM"
        eStrFullName="Mean wave length"
        eStrUnit="meter"
      ELSE IF (IDX.eq.7) THEN
        eStr="ETOTC"
        eStrFullName="model variable"
        eStrUnit="unk"
      ELSE IF (IDX.eq.8) THEN
        eStr="ETOTS"
        eStrFullName="model variable"
        eStrUnit="unk"
      ELSE IF (IDX.eq.9) THEN
        eStr="DM"
        eStrFullName="Mean wave direction"
        eStrUnit="degree"
      ELSE IF (IDX.eq.10) THEN
        eStr="DSPR"
        eStrFullName="Directional spreading"
        eStrUnit="degree"
      ELSE IF (IDX.eq.11) THEN
        eStr="TPPD"
        eStrFullName="Discrete peak wave period"
        eStrUnit="second"
      ELSE IF (IDX.eq.12) THEN
        eStr="CPPD"
        eStrFullName="Discrete peak wave speed"
        eStrUnit="meter second -1"
      ELSE IF (IDX.eq.13) THEN
        eStr="KPPD"
        eStrFullName="discrete peak wave number"
        eStrUnit="meter-1"
      ELSE IF (IDX.eq.14) THEN
        eStr="CGPD"
        eStrFullName="discrete peak group speed"
        eStrUnit="meter second-1"
      ELSE IF (IDX.eq.15) THEN
        eStr="TPP"
        eStrFullName="peak wave period"
        eStrUnit="second"
      ELSE IF (IDX.eq.16) THEN
        eStr="CPP"
        eStrFullName="peak wave speed"
        eStrUnit="meter second-1"
      ELSE IF (IDX.eq.17) THEN
        eStr="WNPP"
        eStrFullName="unk"
        eStrUnit="unk"
      ELSE IF (IDX.eq.18) THEN
        eStr="CGPP"
        eStrFullName="peak group velocity"
        eStrUnit="meter second-1"
      ELSE IF (IDX.eq.19) THEN
        eStr="KPP"
        eStrFullName="peak wave number"
        eStrUnit="meter-1"
      ELSE IF (IDX.eq.20) THEN
        eStr="LPP"
        eStrFullName="Peak wave length"
        eStrUnit="meter"
      ELSE IF (IDX.eq.21) THEN
        eStr="PEAKD"
        eStrFullName="Peak wave direction"
        eStrUnit="degree"
      ELSE IF (IDX.eq.22) THEN
        eStr="PEAKDSPR"
        eStrFullName="Peak directional spreading"
        eStrUnit="degree"
      ELSE IF (IDX.eq.23) THEN
        eStr="DPEAK"
        eStrFullName="discrete peak direction"
        eStrUnit="degree"
      ELSE IF (IDX.eq.24) THEN
        eStr="UBOT"
        eStrFullName="wind-induced bottom orbital velocity"
        eStrUnit="meter second-1"
      ELSE IF (IDX.eq.25) THEN
        eStr="ORBITAL"
        eStrFullName="unk"
        eStrUnit="unk"
      ELSE IF (IDX.eq.26) THEN
        eStr="BOTEXPER"
        eStrFullName="unk"
        eStrUnit="unk"
      ELSE IF (IDX.eq.27) THEN
        eStr="TMBOT"
        eStrFullName="unk"
        eStrUnit="unk"
      ELSE IF (IDX.eq.28) THEN
        eStr="URSELL"
        eStrFullName="ursell number"
        eStrUnit="non-dimensional"
      ELSE IF (IDX.eq.29) THEN
        eStr="UFRIC"
        eStrFullName="air friction velocity"
        eStrUnit="meter second-1"
      ELSE IF (IDX.eq.30) THEN
        eStr="Z0"
        eStrFullName="air roughness length"
        eStrUnit="meter"
      ELSE IF (IDX.eq.31) THEN
        eStr="ALPHA_CH"
        eStrFullName="air Charnock coefficient"
        eStrUnit="non-dimensional"
      ELSE IF (IDX.eq.32) THEN
        eStr="Uwind"
        eStrFullName="wind in X direction"
        eStrUnit="meter second-1"
      ELSE IF (IDX.eq.33) THEN
        eStr="Vwind"
        eStrFullName="wind in Y direction"
        eStrUnit="meter second-1"
      ELSE IF (IDX.eq.34) THEN
        eStr="CD"
        eStrFullName="drag coefficient"
        eStrUnit="non-dimensional"
      ELSE IF (IDX.eq.35) THEN
        eStr="CURTX"
        eStrFullName="current in X direction"
        eStrUnit="meter second-1"
      ELSE IF (IDX.eq.36) THEN
        eStr="CURTY"
        eStrFullName="current in Y direction"
        eStrUnit="meter second-1"
      ELSE IF (IDX.eq.37) THEN
        eStr="WATLEV"
        eStrFullName="water level"
        eStrUnit="meter"
      ELSE IF (IDX.eq.38) THEN
        eStr="WATLEVOLD"
        eStrFullName="water level at previous time step"
        eStrUnit="meter"
      ELSE IF (IDX.eq.39) THEN
        eStr="DEPDT"
        eStrFullName="water level change"
        eStrUnit="meter second-1"
      ELSE IF (IDX.eq.40) THEN
        eStr="DEP"
        eStrFullName="bathymetry"
        eStrUnit="meter"
      ELSE IF (IDX.eq.41) THEN
        eStr="WINDMAG"
        eStrFullName="10-m wind magnitude"
        eStrUnit="meter second-1"
      ELSE IF (IDX.eq.42) THEN
        eStr="TAUW"
        eStrFullName="wave supported surface stress"
        eStrUnit="unk"
      ELSE IF (IDX.eq.43) THEN
        eStr="TAUWX"
        eStrFullName="wave supported surface X-stress"
        eStrUnit="unk"
      ELSE IF (IDX.eq.44) THEN
        eStr="TAUWY"
        eStrFullName="wave supported surface Y-stress"
        eStrUnit="unk"
      ELSE IF (IDX.eq.45) THEN
        eStr="TAUHF"
        eStrFullName="high frequency surface stress"
        eStrUnit="unk"
      ELSE IF (IDX.eq.46) THEN
        eStr="TAUTOT"
        eStrFullName="total surface stress"
        eStrUnit="unk"
      ELSE IF (IDX.eq.47) THEN
        eStr="STOKESBOTTX"
        eStrFullName="bottom Stokes velocity in X direction"
        eStrUnit="meter second -1"
      ELSE IF (IDX.eq.48) THEN
        eStr="STOKESBOTTY"
        eStrFullName="bottom Stokes velocity in Y direction"
        eStrUnit="meter second -1"
      ELSE IF (IDX.eq.49) THEN
        eStr="STOKESSURFX"
        eStrFullName="surface Stokes velocity in X direction"
        eStrUnit="meter second -1"
      ELSE IF (IDX.eq.50) THEN
        eStr="STOKESSURFY"
        eStrFullName="surface Stokes velocity in Y direction"
        eStrUnit="meter second -1"
      ELSE IF (IDX.eq.51) THEN
        eStr="STOKESBAROX"
        eStrFullName="barotropic Stokes velocity in X direction"
        eStrUnit="meter second -1"
      ELSE IF (IDX.eq.52) THEN
        eStr="STOKESBAROY"
        eStrFullName="barotropic Stokes velocity in Y direction"
        eStrUnit="meter second -1"
      ELSE IF (IDX.eq.53) THEN
        eStr="RSXX"
        eStrFullName="barotropic Stress potential Sxx"
        eStrUnit="meter2 second2"
      ELSE IF (IDX.eq.54) THEN
        eStr="RSXY"
        eStrFullName="barotropic Stress potential Sxy"
        eStrUnit="meter2 second2"
      ELSE IF (IDX.eq.55) THEN
        eStr="RSYY"
        eStrFullName="barotropic Stress potential Syy"
        eStrUnit="meter2 second2"
      ELSE IF (IDX.eq.56) THEN
        eStr="CFL1"
        eStrFullName="CFL number 1"
        eStrUnit="non-dimensional"
      ELSE IF (IDX.eq.57) THEN
        eStr="CFL2"
        eStrFullName="CFL number 2"
        eStrUnit="non-dimensional"
      ELSE IF (IDX.eq.58) THEN
        eStr="CFL3"
        eStrFullName="CFL number 3"
        eStrUnit="non-dimensional"
      ELSE IF (IDX.eq.59) THEN
        eStr="ZETA_SETUP"
        eStrFullName="Free-surface elevation induced setup"
        eStrUnit="m"
      ELSE
        CALL WWM_ABORT('Wrong Number')
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WRITE_NETCDF_TIME_HEADER(ncid, nbTime, ntime_dims)
      USE DATAPOOL
      USE NETCDF
      IMPLICIT NONE
      integer, intent(in) :: ncid, nbTime
      integer, intent(inout) :: ntime_dims
      character (len = *), parameter :: UNITS = "units"
      character (len = *), parameter :: CallFct="WRITE_NETCDF_TIME_HEADER"
      integer iret, fifteen_dims, var_id
      iret=nf90_inq_dimid(ncid, 'fifteen', fifteen_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 1, iret)
      IF (nbTime.gt.0) THEN
        iret = nf90_def_dim(ncid, 'ocean_time', nbTime, ntime_dims)
        CALL GENERIC_NETCDF_ERROR(CallFct, 2, iret)
      ELSE
        iret = nf90_def_dim(ncid, 'ocean_time', NF90_UNLIMITED, ntime_dims)
        CALL GENERIC_NETCDF_ERROR(CallFct, 3, iret)
      END IF
      iret=nf90_def_var(ncid,'ocean_time',NF90_RUNTYPE,(/ ntime_dims/), var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 4, iret)
      iret=nf90_put_att(ncid,var_id,UNITS,'seconds since 1858-11-17 00:00:00')
      CALL GENERIC_NETCDF_ERROR(CallFct, 5, iret)
      iret=nf90_put_att(ncid,var_id,"calendar",'gregorian')
      CALL GENERIC_NETCDF_ERROR(CallFct, 6, iret)
      !
      iret=nf90_def_var(ncid,'ocean_time_day',NF90_RUNTYPE,(/ ntime_dims/), var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 7, iret)
      iret=nf90_put_att(ncid,var_id,UNITS,'days since 1858-11-17 00:00:00')
      CALL GENERIC_NETCDF_ERROR(CallFct, 8, iret)
      iret=nf90_put_att(ncid,var_id,"calendar",'gregorian')
      CALL GENERIC_NETCDF_ERROR(CallFct, 9, iret)
      !
      iret=nf90_def_var(ncid,'ocean_time_str',NF90_CHAR,(/ fifteen_dims, ntime_dims/), var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 10, iret)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WRITE_NETCDF_TIME(ncid, idx, eTimeDay)
      USE DATAPOOL, only : DAY2SEC,RKIND, wwmerr
      USE NETCDF
      IMPLICIT NONE
      integer, intent(in) :: ncid, idx
      REAL(rkind), intent(IN) :: eTimeDay
      character (len = *), parameter :: CallFct="WRITE_NETCDF_TIME"
      integer oceantimeday_id, oceantimestr_id, oceantime_id
      integer iret, I
      CHARACTER          :: eChar
      REAL(rkind) eTimeSec
      CHARACTER(LEN=15) :: eTimeStr
      !
      CALL MJD2CT(eTimeDay,eTimeStr)
      eTimeSec=eTimeDay*DAY2SEC
      iret=nf90_inq_varid(ncid, 'ocean_time', oceantime_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 1, iret)
      iret=nf90_put_var(ncid,oceantime_id,eTimeSec,start=(/idx/) )
      CALL GENERIC_NETCDF_ERROR(CallFct, 2, iret)
      !
      iret=nf90_inq_varid(ncid, 'ocean_time_day', oceantimeday_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 3, iret)
      iret=nf90_put_var(ncid,oceantimeday_id,eTimeDay,start=(/idx/) )
      CALL GENERIC_NETCDF_ERROR(CallFct, 4, iret)
      !
      iret=nf90_inq_varid(ncid, 'ocean_time_str', oceantimestr_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 5, iret)
      DO i=1,15
        eChar=eTimeStr(i:i)
        iret=nf90_put_var(ncid,oceantimestr_id,eChar,start=(/i, idx/) )
        CALL GENERIC_NETCDF_ERROR(CallFct, 6, iret)
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GET_IOBPD_OUTPUT(IOBPDoutput, np_write)
      USE DATAPOOL, only : MDC, IOBPD, np_total, MULTIPLEOUT_HIS, MNP
# ifdef MPI_PARALL_GRID
      USE datapool, only : iplg, comm, nproc, istatus, ierr, myrank, itype
# endif
      IMPLICIT NONE
      INTEGER, intent(in)  :: np_write
      INTEGER, INTENT(OUT) :: IOBPDoutput(MDC, np_write)
# ifdef MPI_PARALL_GRID
      integer, allocatable :: rIOBPD(:,:), rStatus(:), Status(:)
      integer iProc, IP, istat
      IF (MULTIPLEOUT_HIS .eq. 1) THEN
        IOBPDoutput=IOBPD
      ELSE
        allocate(rIOBPD(MDC, np_total), rStatus(np_total), Status(np_total), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_netcdf, allocate error 7')
        IOBPDoutput=0
        Status=0
        DO IP=1,MNP
          IOBPDoutput(:, iplg(IP))=IOBPD(:,IP)
          Status(iplg(IP)) = 1
        END DO
        IF (myrank .eq. 0) THEN
          DO iProc=2,nproc
            CALL MPI_RECV(rIOBPD, MDC*np_total, itype, iProc-1, 193, comm, istatus, ierr)
            CALL MPI_RECV(rStatus, np_total, itype, iProc-1, 197, comm, istatus, ierr)
            DO IP=1,np_total
              IF (rStatus(IP) .eq. 1) THEN
                IOBPDoutput(:,IP)=rIOBPD(:,IP)
              END IF
            END DO
          END DO
          DO iProc=2,nproc
            CALL MPI_SEND(IOBPDoutput, MDC*np_total, itype, iProc-1, 199, comm, ierr)
          END DO
        ELSE
          CALL MPI_SEND(IOBPDoutput, MDC*np_total, itype, 0, 193, comm, ierr)
          CALL MPI_SEND(Status, np_total, itype, 0, 197, comm, ierr)
          CALL MPI_RECV(IOBPDoutput, MDC*np_total, itype, 0, 199, comm, istatus, ierr)
        ENDIF
        deallocate(rIOBPD, rStatus, Status)
      END IF
# else
      IOBPDoutput=IOBPD
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WRITE_NETCDF_HEADERS_1(ncid, nbTime, MULTIPLEOUT, np_write, ne_write)
      USE DATAPOOL
      USE NETCDF
      implicit none
      integer, intent(in) :: ncid, nbTime, MULTIPLEOUT
      integer, intent(in) :: np_write, ne_write
      character (len = *), parameter :: UNITS = "units"
      integer one_dims, two_dims, three_dims, fifteen_dims
      integer mnp_dims, mne_dims, msc_dims, mdc_dims
# ifdef MPI_PARALL_GRID
      integer np_global_dims, ne_global_dims
# endif
      integer iret, var_id
      integer ntime_dims
      integer p_dims, e_dims
      integer Oper
      character (len = *), parameter :: CallFct="WRITE_NETCDF_HEADERS_1"
      IF ((np_write.eq.0).or.(ne_write.eq.0)) THEN
        CALL WWM_ABORT('np_write=0 or ne_write=0, not allowed by any mean')
      ENDIF
      iret = nf90_def_dim(ncid, 'one', 1, one_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 1, iret)
      iret = nf90_def_dim(ncid, 'two', 2, two_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 2, iret)
      iret = nf90_def_dim(ncid, 'three', 3, three_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 3, iret)
      iret = nf90_def_dim(ncid, 'fifteen', 15, fifteen_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 4, iret)
      iret = nf90_def_dim(ncid, 'mnp', np_write, mnp_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 5, iret)
      iret = nf90_def_dim(ncid, 'mne', ne_write, mne_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 6, iret)
      iret = nf90_def_dim(ncid, 'msc', MSC, msc_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 7, iret)
      iret = nf90_def_dim(ncid, 'mdc', MDC, mdc_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 8, iret)
# ifdef MPI_PARALL_GRID
      IF (MULTIPLEOUT.eq.1) THEN
        iret = nf90_def_dim(ncid, 'np_global', np_global, np_global_dims)
        CALL GENERIC_NETCDF_ERROR(CallFct, 9, iret)
        iret = nf90_def_dim(ncid, 'ne_global', ne_global, ne_global_dims)
        CALL GENERIC_NETCDF_ERROR(CallFct, 9, iret)
        iret = nf90_def_var(ncid,'iplg',NF90_INT,(/ mnp_dims/),var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 10, iret)
        iret = nf90_put_att(ncid,var_id,'description','local to global indexes')
        CALL GENERIC_NETCDF_ERROR(CallFct, 11, iret)
      END IF
      iret=nf90_def_var(ncid,'nproc',NF90_INT,(/one_dims/),var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 12, iret)
      iret=nf90_put_att(ncid,var_id,UNITS,'integer')
      CALL GENERIC_NETCDF_ERROR(CallFct, 13, iret)
      iret=nf90_put_att(ncid,var_id,'description','number of processors')
      CALL GENERIC_NETCDF_ERROR(CallFct, 14, iret)
      !
      iret=nf90_def_var(ncid,'MULTIPLEOUT',NF90_INT,(/one_dims/),var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 15, iret)
      iret=nf90_put_att(ncid,var_id,UNITS,'integer')
      CALL GENERIC_NETCDF_ERROR(CallFct, 16, iret)
      iret=nf90_put_att(ncid,var_id,'description','multiple status')
      CALL GENERIC_NETCDF_ERROR(CallFct, 17, iret)
# endif
      !
      IF (PARAMWRITE_HIS) THEN
        CALL WRITE_PARAM_1(ncid, one_dims)
      ENDIF
      !
      CALL WRITE_NETCDF_TIME_HEADER(ncid, nbTime, ntime_dims)
      !
      IF (GRIDWRITE) THEN
# ifdef MPI_PARALL_GRID
        IF (MULTIPLEOUT.eq.1) THEN
          e_dims=ne_global_dims
          p_dims=np_global_dims
        ELSE
          e_dims=mne_dims
          p_dims=mnp_dims
        END IF
# else
        e_dims=mne_dims
        p_dims=mnp_dims
# endif
! element
        iret=nf90_def_var(ncid,'ele',NF90_INT,(/three_dims, e_dims/),var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 27, iret)
        iret=nf90_put_att(ncid,var_id,UNITS,'non-dimensional')
        CALL GENERIC_NETCDF_ERROR(CallFct, 28, iret)
! lon
        IF (LSPHE) THEN
          iret=nf90_def_var(ncid,"lon",NF90_RUNTYPE,(/ p_dims/),var_id)
        ELSE
          iret=nf90_def_var(ncid,"x",NF90_RUNTYPE,(/ p_dims/),var_id)
        END IF
        CALL GENERIC_NETCDF_ERROR(CallFct, 29, iret)
        IF (LSPHE) THEN
          iret=nf90_put_att(ncid,var_id,UNITS,'degree')
        ELSE
          iret=nf90_put_att(ncid,var_id,UNITS,'meter')
        END IF
        CALL GENERIC_NETCDF_ERROR(CallFct, 30, iret)
! lat
        IF (LSPHE) THEN
          iret=nf90_def_var(ncid,"lat",NF90_RUNTYPE,(/ p_dims/),var_id)
        ELSE
          iret=nf90_def_var(ncid,"y",NF90_RUNTYPE,(/ p_dims/),var_id)
        END IF
        CALL GENERIC_NETCDF_ERROR(CallFct, 31, iret)
        IF (LSPHE) THEN
          iret=nf90_put_att(ncid,var_id,UNITS,'degree')
        ELSE
          iret=nf90_put_att(ncid,var_id,UNITS,'meter')
        END IF
        CALL GENERIC_NETCDF_ERROR(CallFct, 32, iret)
! IOBP
        iret=nf90_def_var(ncid,"IOBP",NF90_RUNTYPE,(/ p_dims/),var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 33, iret)
! depth
        iret=nf90_def_var(ncid,'depth',NF90_RUNTYPE,(/ p_dims/),var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 34, iret)
        iret=nf90_put_att(ncid,var_id,UNITS,'meters')
        CALL GENERIC_NETCDF_ERROR(CallFct, 35, iret)
! boundary
        Oper=1
        CALL SERIAL_WRITE_BOUNDARY(ncid, np_total, ne_total, INEtotal, Oper)
        !
      END IF
      IF (IOBPD_HISTORY) THEN
# ifdef MPI_PARALL_GRID
        IF (MULTIPLEOUT.eq.1) THEN
          p_dims=np_global_dims
        ELSE
          p_dims=mnp_dims
        END IF
# else
        p_dims=mnp_dims
# endif
        iret=nf90_def_var(ncid,'IOBPD',NF90_INT,(/ mdc_dims, p_dims, ntime_dims/), var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 20, iret)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WRITE_NETCDF_HEADERS_2(ncid, MULTIPLEOUT, WriteOutputProcess, np_write, ne_write)
      USE DATAPOOL
      USE NETCDF
      implicit none
      integer, intent(in) :: ncid, MULTIPLEOUT
      logical, intent(in) :: WriteOutputProcess
      integer, intent(in) :: np_write, ne_write
      integer var_id, iret
      character (len = *), parameter :: CallFct="WRITE_NETCDF_HEADERS_2"
      integer Oper
#ifdef MPI_PARALL_GRID
      integer eInt(1)
#endif
# ifdef MPI_PARALL_GRID
      IF (MULTIPLEOUT.eq.1) THEN
        iret=nf90_inq_varid(ncid, 'iplg', var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 1, iret)
        iret=nf90_put_var(ncid,var_id,iplg,start=(/1/), count = (/ np_write /))
        CALL GENERIC_NETCDF_ERROR(CallFct, 2, iret)
      END IF
      !
      IF (WriteOutputProcess) THEN
        iret=nf90_inq_varid(ncid,'nproc',var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 3, iret)
        eInt(1)=nproc
        iret=nf90_put_var(ncid,var_id,eInt,start = (/1/), count=(/1/))
        CALL GENERIC_NETCDF_ERROR(CallFct, 4, iret)
        !
        iret=nf90_inq_varid(ncid,'MULTIPLEOUT',var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 5, iret)
        eInt(1)=MULTIPLEOUT
        iret=nf90_put_var(ncid,var_id,eInt,start = (/1/), count=(/1/))
        CALL GENERIC_NETCDF_ERROR(CallFct, 6, iret)
      END IF
# endif
      IF (PARAMWRITE_HIS.and.WriteOutputProcess) THEN
        CALL WRITE_PARAM_2(ncid)
      END IF
      IF (GRIDWRITE.and.WriteOutputProcess) THEN
        !
        iret=nf90_inq_varid(ncid, "ele", var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 7, iret)
        iret=nf90_put_var(ncid,var_id,INEtotal, start = (/1,1/), count = (/ 3, ne_total/))
        CALL GENERIC_NETCDF_ERROR(CallFct, 8, iret)
        !
        IF (LSPHE) THEN
          iret=nf90_inq_varid(ncid, "lon", var_id)
        ELSE
          iret=nf90_inq_varid(ncid, "x", var_id)
        ENDIF
        CALL GENERIC_NETCDF_ERROR(CallFct, 9, iret)
        iret=nf90_put_var(ncid,var_id,XPtotal, start = (/1/), count = (/ np_total/))
        CALL GENERIC_NETCDF_ERROR(CallFct, 10, iret)
        !
        IF (LSPHE) THEN
          iret=nf90_inq_varid(ncid, "lat", var_id)
        ELSE
          iret=nf90_inq_varid(ncid, "y", var_id)
        ENDIF
        CALL GENERIC_NETCDF_ERROR(CallFct, 11, iret)
        iret=nf90_put_var(ncid,var_id,YPtotal, start = (/1/), count = (/ np_total/))
        CALL GENERIC_NETCDF_ERROR(CallFct, 12, iret)
        !
        iret=nf90_inq_varid(ncid, "IOBP", var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 13, iret)
        iret=nf90_put_var(ncid,var_id,IOBPtotal, start = (/1/), count = (/ np_total/))
        CALL GENERIC_NETCDF_ERROR(CallFct, 14, iret)
        !
        iret=nf90_inq_varid(ncid, "depth", var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 15, iret)
        iret=nf90_put_var(ncid,var_id,DEPtotal, start = (/1/), count = (/ np_write/))
        CALL GENERIC_NETCDF_ERROR(CallFct, 16, iret)
        !
        Oper=2
        CALL SERIAL_WRITE_BOUNDARY(ncid, np_total, ne_total, INEtotal, Oper)
      ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WRITE_PARAM_1(ncid, one_dims)
      USE NETCDF
      USE DATAPOOL, only : NF90_RUNTYPE
      implicit none
      integer, intent(in) :: ncid, one_dims
      integer :: iret, var_id
      character (len = *), parameter :: UNITS = "units"
      iret = nf90_def_var(ncid,'frlow', NF90_RUNTYPE,(/ one_dims/), var_id)
      CALL REPORT_ERROR_DEF(iret, 'frlow')
      iret = nf90_put_att(ncid,var_id,UNITS,'lower_frequency')
      !
      iret = nf90_def_var(ncid,'frhigh',NF90_RUNTYPE,(/ one_dims/),var_id)
      CALL REPORT_ERROR_DEF(iret, 'frhigh')
      iret = nf90_put_att(ncid,var_id,UNITS,'higher_frequency')
      !
      iret = nf90_def_var(ncid,'MESNL',NF90_INT,(/ one_dims/),var_id)
      CALL REPORT_ERROR_DEF(iret, 'MESNL')
      iret = nf90_put_att(ncid,var_id,'description','nonlinear interaction nl4')
      !
      iret = nf90_def_var(ncid,'MESIN',NF90_INT,(/ one_dims/),var_id)
      CALL REPORT_ERROR_DEF(iret, 'MESIN')
      iret = nf90_put_att(ncid,var_id,'description','wind input source term')
      !
      iret = nf90_def_var(ncid,'MESDS',NF90_INT,(/ one_dims/),var_id)
      CALL REPORT_ERROR_DEF(iret, 'MESDS')
      iret = nf90_put_att(ncid,var_id,'description','dissipation source term')
      !
      iret = nf90_def_var(ncid,'MESBF',NF90_INT,(/ one_dims/),var_id)
      CALL REPORT_ERROR_DEF(iret, 'MESBF')
      iret = nf90_put_att(ncid,var_id,'description','bottom friction')
      !
      iret = nf90_def_var(ncid,'ICOMP',NF90_INT,(/ one_dims/),var_id)
      CALL REPORT_ERROR_DEF(iret, 'ICOMP')
      iret = nf90_put_att(ncid,var_id,'description','implicitness')
      !
      iret = nf90_def_var(ncid,'AMETHOD',NF90_INT,(/ one_dims/),var_id)
      CALL REPORT_ERROR_DEF(iret, 'AMETHOD')
      iret = nf90_put_att(ncid,var_id,'description','advection method')
      !
      iret = nf90_def_var(ncid,'FMETHOD',NF90_INT,(/ one_dims/),var_id)
      CALL REPORT_ERROR_DEF(iret, 'FMETHOD')
      iret = nf90_put_att(ncid,var_id,'description','frequency shifting method')
      !
      iret = nf90_def_var(ncid,'DMETHOD',NF90_INT,(/ one_dims/),var_id)
      CALL REPORT_ERROR_DEF(iret, 'DMETHOD')
      iret = nf90_put_att(ncid,var_id,'description','directional shifting method(refraction)')
      !
      iret = nf90_def_var(ncid,'SMETHOD',NF90_INT,(/ one_dims/),var_id)
      CALL REPORT_ERROR_DEF(iret, 'SMETHOD')
      iret = nf90_put_att(ncid,var_id,'description','source term integration method')
      !
      iret = nf90_def_var(ncid,'LSPHE',NF90_INT,(/ one_dims/),var_id)
      CALL REPORT_ERROR_DEF(iret, 'LSPHE')
      iret = nf90_put_att(ncid,var_id,'description','spherical coordinates')
      END SUBROUTINE WRITE_PARAM_1
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WRITE_PARAM_2(ncid)
      USE DATAPOOL, only : FRLOW, FRHIGH, ICOMP, AMETHOD, FMETHOD,        &
     &    DMETHOD, SMETHOD, MESIN, MESBF, MESDS, MESNL, LSPHE
      USE NETCDF
      IMPLICIT NONE
      integer, intent(in) :: ncid
      integer iret, var_id
      integer LSPHE_INT
      character (len = *), parameter :: CallFct="WRITE_PARAM_2"
      iret=nf90_inq_varid(ncid, "frlow", var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 1, iret)
      iret=nf90_put_var(ncid,var_id,FRLOW,start=(/1/) )
      CALL GENERIC_NETCDF_ERROR(CallFct, 2, iret)
      !
      iret=nf90_inq_varid(ncid, "frhigh", var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 3, iret)
      iret=nf90_put_var(ncid,var_id,FRHIGH, start=(/1/) )
      CALL GENERIC_NETCDF_ERROR(CallFct, 4, iret)
      !
      iret=nf90_inq_varid(ncid, "MESIN", var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 15, iret)
      iret = nf90_put_var(ncid,var_id,MESIN, start=(/1/) )
      CALL GENERIC_NETCDF_ERROR(CallFct, 16, iret)
      !
      iret=nf90_inq_varid(ncid, "MESBF", var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 17, iret)
      iret=nf90_put_var(ncid,var_id,MESBF, start=(/1/) )
      CALL GENERIC_NETCDF_ERROR(CallFct, 18, iret)
      !
      iret=nf90_inq_varid(ncid, "MESDS", var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 19, iret)
      iret=nf90_put_var(ncid,var_id,MESDS, start=(/1/) )
      CALL GENERIC_NETCDF_ERROR(CallFct, 20, iret)
      !
      iret=nf90_inq_varid(ncid, "MESNL", var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 21, iret)
      iret=nf90_put_var(ncid,var_id,MESNL, start=(/1/) )
      CALL GENERIC_NETCDF_ERROR(CallFct, 22, iret)
      !
      iret=nf90_inq_varid(ncid, "ICOMP", var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 5, iret)
      iret=nf90_put_var(ncid,var_id,ICOMP, start=(/1/) )
      CALL GENERIC_NETCDF_ERROR(CallFct, 6, iret)
      !
      iret=nf90_inq_varid(ncid, "AMETHOD", var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 7, iret)
      iret=nf90_put_var(ncid,var_id,AMETHOD, start=(/1/) )
      CALL GENERIC_NETCDF_ERROR(CallFct, 8, iret)
      !
      iret=nf90_inq_varid(ncid, "FMETHOD", var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 9, iret)
      iret=nf90_put_var(ncid,var_id,FMETHOD, start=(/1/) )
      CALL GENERIC_NETCDF_ERROR(CallFct, 10, iret)
      !
      iret=nf90_inq_varid(ncid, "DMETHOD", var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 11, iret)
      iret=nf90_put_var(ncid,var_id,DMETHOD, start=(/1/) )
      CALL GENERIC_NETCDF_ERROR(CallFct, 12, iret)
      !
      iret=nf90_inq_varid(ncid, "SMETHOD", var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 13, iret)
      iret=nf90_put_var(ncid,var_id,SMETHOD, start=(/1/) )
      CALL GENERIC_NETCDF_ERROR(CallFct, 14, iret)
      !
      iret=nf90_inq_varid(ncid, "LSPHE", var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 13, iret)
      IF (LSPHE) THEN
        LSPHE_INT=1
      ELSE
        LSPHE_INT=0
      END IF
      iret=nf90_put_var(ncid,var_id,LSPHE_INT, start=(/1/) )
      CALL GENERIC_NETCDF_ERROR(CallFct, 14, iret)
      END SUBROUTINE WRITE_PARAM_2
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DEFINE_STATION_NC(FILE_NAME, MULTIPLEOUT)
      USE NETCDF
      USE DATAPOOL
      implicit none
      character(len=256), intent(in) :: FILE_NAME
      integer, intent(in) :: MULTIPLEOUT
      character (len = *), parameter :: CallFct="DEFINE_STATION_NC"
      character (len = *), parameter :: UNITS = "units"
      character (len = *), parameter :: FULLNAME = "full-name"
      character(len=40) :: eStr, eStrUnit
      character(len=80) :: eStrFullName
      integer iret, ncid, nbstat_dims, ntime_dims, msc_dims, mdc_dims
      integer one_dims, three_dims, var_id, I
      iret = nf90_create(TRIM(FILE_NAME), NF90_CLOBBER, ncid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 1, iret)

      CALL WRITE_NETCDF_HEADERS_STAT_1(ncid, -1, MULTIPLEOUT)

      iret=nf90_inq_dimid(ncid, 'nbstation', nbstat_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 2, iret)

      iret=nf90_inq_dimid(ncid, 'ocean_time', ntime_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 3, iret)

      iret=nf90_inq_dimid(ncid, 'msc', msc_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 4, iret)

      iret=nf90_inq_dimid(ncid, 'mdc', mdc_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 5, iret)

      iret=nf90_inq_dimid(ncid, 'one', one_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 6, iret)

      iret=nf90_inq_dimid(ncid, 'three', three_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 7, iret)

      IF (VAROUT_STATION%AC) THEN
        iret=nf90_def_var(ncid,'AC',NF90_OUTTYPE_STAT,(/nbstat_dims, msc_dims, mdc_dims, ntime_dims /),var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 11, iret)

        iret=nf90_put_att(ncid,var_id,UNITS,'unk')
        CALL GENERIC_NETCDF_ERROR(CallFct, 12, iret)

        iret=nf90_put_att(ncid,var_id,FULLNAME,'spectral energy density')
        CALL GENERIC_NETCDF_ERROR(CallFct, 13, iret)
      END IF
      IF (VAROUT_STATION%WK) THEN
        iret=nf90_def_var(ncid,'WK',NF90_OUTTYPE_STAT,(/nbstat_dims, msc_dims, ntime_dims /),var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 14, iret)

        iret=nf90_put_att(ncid,var_id,UNITS,'unk')
        CALL GENERIC_NETCDF_ERROR(CallFct, 15, iret)

        iret=nf90_put_att(ncid,var_id,FULLNAME,'wave number by frequency')
        CALL GENERIC_NETCDF_ERROR(CallFct, 16, iret)
      END IF
      IF (VAROUT_STATION%ACOUT_1D) THEN
        iret=nf90_def_var(ncid,'ACOUT_1D',NF90_OUTTYPE_STAT,(/nbstat_dims, msc_dims, three_dims, ntime_dims /),var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 17, iret)

        iret=nf90_put_att(ncid,var_id,UNITS,'unk')
        CALL GENERIC_NETCDF_ERROR(CallFct, 18, iret)

        iret=nf90_put_att(ncid,var_id,FULLNAME,'1-dimensional spectrum')
        CALL GENERIC_NETCDF_ERROR(CallFct, 19, iret)
      END IF
      IF (VAROUT_STATION%ACOUT_2D) THEN
        iret=nf90_def_var(ncid,'ACOUT_2D',NF90_OUTTYPE_STAT,(/nbstat_dims, msc_dims, mdc_dims, ntime_dims /),var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 20, iret)

        iret=nf90_put_att(ncid,var_id,UNITS,'unk')
        CALL GENERIC_NETCDF_ERROR(CallFct, 21, iret)

        iret=nf90_put_att(ncid,var_id,FULLNAME,'2-dimensional spectrum')
        CALL GENERIC_NETCDF_ERROR(CallFct, 22, iret)
      END IF
      DO I=1,OUTVARS_COMPLETE
        IF (VAROUT_STATION%LVAR(I)) THEN
          CALL NAMEVARIABLE(I, eStr, eStrFullName, eStrUnit)
          iret=nf90_def_var(ncid,TRIM(eStr),NF90_OUTTYPE_STAT,(/ nbstat_dims, ntime_dims /),var_id)
          CALL GENERIC_NETCDF_ERROR(CallFct, 23, iret)

          iret=nf90_put_att(ncid,var_id,UNITS,TRIM(eStrUnit))
          CALL GENERIC_NETCDF_ERROR(CallFct, 24, iret)

          iret=nf90_put_att(ncid,var_id,FULLNAME,TRIM(eStrFullName))
          CALL GENERIC_NETCDF_ERROR(CallFct, 25, iret)
        END IF
      END DO
      iret=nf90_close(ncid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 26, iret)

      iret=nf90_open(TRIM(FILE_NAME), NF90_WRITE, ncid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 27, iret)

      CALL WRITE_NETCDF_HEADERS_STAT_2(ncid, MULTIPLEOUT_STAT)
      iret=nf90_close(ncid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 28, iret)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WRITE_NETCDF_HEADERS_STAT_1(ncid, nbTime, MULTIPLEOUT)
      USE DATAPOOL
      USE NETCDF
      implicit none
      integer, intent(in) :: ncid, nbTime, MULTIPLEOUT
      character (len = *), parameter :: CallFct="WRITE_NETCDF_HEADERS_STAT_1"
      character (len = *), parameter :: UNITS = "units"
      integer one_dims, two_dims, three_dims, fifteen_dims
      integer msc_dims, mdc_dims
      integer nbstat_dims
      integer iret, var_id
      integer ntime_dims
      iret = nf90_def_dim(ncid, 'one', 1, one_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 1, iret)
      iret = nf90_def_dim(ncid, 'two', 2, two_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 2, iret)
      iret = nf90_def_dim(ncid, 'three', 3, three_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 3, iret)
      iret = nf90_def_dim(ncid, 'fifteen', 15, fifteen_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 4, iret)
      iret = nf90_def_dim(ncid, 'nbstation', IOUTS, nbstat_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 5, iret)
      iret = nf90_def_dim(ncid, 'msc', MSC, msc_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 6, iret)
      iret = nf90_def_dim(ncid, 'mdc', MDC, mdc_dims)
      CALL GENERIC_NETCDF_ERROR(CallFct, 7, iret)
      !
# ifdef MPI_PARALL_GRID
      iret=nf90_def_var(ncid,'nproc',NF90_INT,(/one_dims/),var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 8, iret)
      iret=nf90_put_att(ncid,var_id,UNITS,'integer')
      CALL GENERIC_NETCDF_ERROR(CallFct, 9, iret)
      iret=nf90_put_att(ncid,var_id,'description','number of processors')
      CALL GENERIC_NETCDF_ERROR(CallFct, 10, iret)
      !
      iret=nf90_def_var(ncid,'MULTIPLEOUT',NF90_INT,(/one_dims/),var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 11, iret)
      iret=nf90_put_att(ncid,var_id,UNITS,'integer')
      CALL GENERIC_NETCDF_ERROR(CallFct, 12, iret)
      iret=nf90_put_att(ncid,var_id,'description','multiple status')
      CALL GENERIC_NETCDF_ERROR(CallFct, 13, iret)
# endif
      !
      IF (PARAMWRITE_STAT) THEN
        CALL WRITE_PARAM_1(ncid, one_dims)
      ENDIF
      !
      CALL WRITE_NETCDF_TIME_HEADER(ncid, nbTime, ntime_dims)
      IF (LSPHE) THEN
        iret=nf90_def_var(ncid,"lon",NF90_RUNTYPE,(/ nbstat_dims/),var_id)
      ELSE
        iret=nf90_def_var(ncid,"x",NF90_RUNTYPE,(/ nbstat_dims/),var_id)
      END IF
      CALL GENERIC_NETCDF_ERROR(CallFct, 23, iret)
      IF (LSPHE) THEN
        iret=nf90_put_att(ncid,var_id,UNITS,'degree')
      ELSE
        iret=nf90_put_att(ncid,var_id,UNITS,'meter')
      END IF
      CALL GENERIC_NETCDF_ERROR(CallFct, 24, iret)
! lat
      IF (LSPHE) THEN
        iret=nf90_def_var(ncid,"lat",NF90_RUNTYPE,(/ nbstat_dims/),var_id)
      ELSE
        iret=nf90_def_var(ncid,"y",NF90_RUNTYPE,(/ nbstat_dims/),var_id)
      END IF
      CALL GENERIC_NETCDF_ERROR(CallFct, 25, iret)
      IF (LSPHE) THEN
        iret=nf90_put_att(ncid,var_id,UNITS,'degree')
      ELSE
        iret=nf90_put_att(ncid,var_id,UNITS,'meter')
      END IF
      CALL GENERIC_NETCDF_ERROR(CallFct, 26, iret)
! cutoff frequency
      iret=nf90_def_var(ncid,'cutoff',NF90_RUNTYPE,(/ nbstat_dims/),var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 27, iret)
      iret=nf90_put_att(ncid,var_id,UNITS,'Hz')
      CALL GENERIC_NETCDF_ERROR(CallFct, 28, iret)
! ismax value
      iret=nf90_def_var(ncid,'ismax',NF90_INT,(/ nbstat_dims/),var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 29, iret)
      iret=nf90_put_att(ncid,var_id,UNITS,'integer')
      CALL GENERIC_NETCDF_ERROR(CallFct, 30, iret)
      IF (MULTIPLEOUT.gt.0) THEN
! Ifound
        iret=nf90_def_var(ncid,'ifound',NF90_INT,(/ nbstat_dims/),var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 31, iret)
        iret=nf90_put_att(ncid,var_id,UNITS,'integer')
        CALL GENERIC_NETCDF_ERROR(CallFct, 32, iret)
      END IF
! Isum
      iret=nf90_def_var(ncid,'isum',NF90_INT,(/ nbstat_dims/),var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 33, iret)
      iret=nf90_put_att(ncid,var_id,UNITS,'integer')
      CALL GENERIC_NETCDF_ERROR(CallFct, 34, iret)
! SPSIG
      iret=nf90_def_var(ncid,'spsig',NF90_RUNTYPE,(/ msc_dims/),var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 35, iret)
      iret=nf90_put_att(ncid,var_id,UNITS,'integer')
      CALL GENERIC_NETCDF_ERROR(CallFct, 36, iret)
! SPDIR
      iret=nf90_def_var(ncid,'spdir',NF90_RUNTYPE,(/ mdc_dims/),var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 31, iret)
      iret=nf90_put_att(ncid,var_id,UNITS,'integer')
      CALL GENERIC_NETCDF_ERROR(CallFct, 37, iret)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WRITE_NETCDF_HEADERS_STAT_2(ncid, MULTIPLEOUT)
      USE DATAPOOL
      USE NETCDF
      implicit none
      integer, intent(in) :: ncid, MULTIPLEOUT
      integer :: eWriteInt(1)
      real(rkind) :: eWriteReal(1)
      integer var_id, iret
      integer I
# ifdef MPI_PARALL_GRID
      integer eInt(1)
# endif
      character (len = *), parameter :: CallFct="WRITE_NETCDF_HEADERS_STAT_2"
      !
# ifdef MPI_PARALL_GRID
      iret=nf90_inq_varid(ncid,'nproc',var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 3, iret)
      eInt(1)=nproc
      iret=nf90_put_var(ncid,var_id,eInt,start = (/1/), count=(/1/))
      CALL GENERIC_NETCDF_ERROR(CallFct, 4, iret)
      !
      iret=nf90_inq_varid(ncid,'MULTIPLEOUT',var_id)
      CALL GENERIC_NETCDF_ERROR(CallFct, 5, iret)
      eInt(1)=MULTIPLEOUT
      iret=nf90_put_var(ncid,var_id,eInt,start = (/1/), count=(/1/))
      CALL GENERIC_NETCDF_ERROR(CallFct, 6, iret)
# endif
      !
      IF (PARAMWRITE_STAT) THEN
        CALL WRITE_PARAM_2(ncid)
      ENDIF
      DO I=1,IOUTS 
        eWriteReal(1)=STATION(I) % XCOORD
        IF (LSPHE) THEN
          iret=nf90_inq_varid(ncid, "lon", var_id)
        ELSE
          iret=nf90_inq_varid(ncid, "x", var_id)
        END IF
        CALL GENERIC_NETCDF_ERROR(CallFct, 1, iret)
        iret=nf90_put_var(ncid,var_id,eWriteReal, start=(/I/), count =(/1/) )
        CALL GENERIC_NETCDF_ERROR(CallFct, 2, iret)
        !
        eWriteReal(1)=STATION(I) % YCOORD
        IF (LSPHE) THEN
          iret=nf90_inq_varid(ncid, "lat", var_id)
        ELSE
          iret=nf90_inq_varid(ncid, "y", var_id)
        END IF
        CALL GENERIC_NETCDF_ERROR(CallFct, 3, iret)
        iret=nf90_put_var(ncid,var_id,eWriteReal, start=(/I/), count =(/1/) )
        CALL GENERIC_NETCDF_ERROR(CallFct, 4, iret)
        !
        eWriteReal(1)=STATION(I) % CUTOFF
        iret=nf90_inq_varid(ncid, "cutoff", var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 5, iret)
        iret=nf90_put_var(ncid,var_id,eWriteReal, start=(/I/), count =(/1/) )
        CALL GENERIC_NETCDF_ERROR(CallFct, 6, iret)
        !
        eWriteInt(1)=STATION(I) % ISMAX
        iret=nf90_inq_varid(ncid, "ismax", var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 7, iret)
        iret=nf90_put_var(ncid,var_id,eWriteInt, start=(/I/), count =(/1/) )
        CALL GENERIC_NETCDF_ERROR(CallFct, 8, iret)
# ifdef MPI_PARALL_GRID
        IF (MULTIPLEOUT.eq.0) THEN
          eWriteInt(1)=STATION(I) % ISUM
          iret=nf90_inq_varid(ncid, "isum", var_id)
          CALL GENERIC_NETCDF_ERROR(CallFct, 9, iret)
          iret=nf90_put_var(ncid,var_id,eWriteInt, start=(/I/), count=(/1/) )
          CALL GENERIC_NETCDF_ERROR(CallFct, 10, iret)
        ELSE
          eWriteInt(1)=STATION(I) % IFOUND
          iret=nf90_inq_varid(ncid, "ifound", var_id)
          CALL GENERIC_NETCDF_ERROR(CallFct, 11, iret)
          iret=nf90_put_var(ncid,var_id,eWriteInt, start=(/I/), count=(/1/) )
          CALL GENERIC_NETCDF_ERROR(CallFct, 12, iret)
          !
          eWriteInt(1)=STATION(I) % ISUM
          iret=nf90_inq_varid(ncid, "isum", var_id)
          CALL GENERIC_NETCDF_ERROR(CallFct, 13, iret)
          iret=nf90_put_var(ncid,var_id,eWriteInt, start=(/I/), count=(/1/) )
          CALL GENERIC_NETCDF_ERROR(CallFct, 14, iret)
        ENDIF
# else
        eWriteInt(1)=STATION(I) % IFOUND
        iret=nf90_inq_varid(ncid, "isum", var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 15, iret)
        iret=nf90_put_var(ncid,var_id,eWriteInt, start=(/I/), count =(/1/) )
        CALL GENERIC_NETCDF_ERROR(CallFct, 16, iret)
# endif
        iret=nf90_inq_varid(ncid, "spsig", var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 17, iret)
        iret=nf90_put_var(ncid,var_id,SPSIG, start=(/1/), count =(/MSC/) )
        CALL GENERIC_NETCDF_ERROR(CallFct, 18, iret)
        !
        iret=nf90_inq_varid(ncid, "spdir", var_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 19, iret)
        iret=nf90_put_var(ncid,var_id,SPDIR, start=(/1/), count =(/MDC/) )
        CALL GENERIC_NETCDF_ERROR(CallFct, 20, iret)
      ENDDO
      END SUBROUTINE
#endif
