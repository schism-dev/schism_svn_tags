      MODULE CF_VARIABLES
      USE NETCDF
      USE schism_glbl, only : npa, np, rkind, pi, np_global
      USE schism_glbl, only : xlon, ylat, iplg
      USE schism_msgp, ONLY : comm, ierr, itype, rtype, myrank
      USE schism_msgp, ONLY : istatus, nproc, parallel_abort
!   Copyright 2014 College of William and Mary
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

!**********************************************************************
!*  Author: Ivica Janeković [ivica.jan@gmail.com]                     *
!*  Author: Mathieu Dutour Sikiric [mathieu.dutour@gmail.com]         *
!**********************************************************************
! This is set of subroutines for handling wind and mslp from netcdf WRF
! to get Sschism forcing fields.
! list of SUBROUTINES:
! INIT_NETCDF_CF 	loads wind_time, compute interp coefs.
! READ_INTERP_NETCDF_CF	reads fields and interpolate them on the FEM
! READ_NETCDF_DIRECT    fileds already interpolated into SCHICM grid
! FIX_COORDS      	tranform coords into Sschism radians
!*********************************************************************
!**********************************************************************
!*  Reading the type of variable                                      *
!*  ---the accessible variable is read from the generic name          *
!**********************************************************************
      TYPE VAR_FILE_ARRAY
        integer nbFile
        character(len=140), allocatable :: ListFile(:)
      END TYPE VAR_FILE_ARRAY
      TYPE VAR_LIST_TIME
        integer nbTime
        REAL(rkind), allocatable :: ListTime(:)
      END TYPE VAR_LIST_TIME
      TYPE VAR_LIST_SERIES_TIME
        integer nbTime
        integer, allocatable :: List_IFILE(:)
        integer, allocatable :: List_ITIME(:)
        REAL(rkind), allocatable :: List_Time(:)
      END TYPE VAR_LIST_SERIES_TIME
      TYPE VAR_NETCDF_CF
        character(len=100) :: eString
        real(rkind) :: cf_scale_factor
        real(rkind) :: cf_add_offset
        integer dimVar
        logical IsDirect
        logical MULTIPLE_IN
        TYPE(VAR_LIST_SERIES_TIME) :: VAR_SERIES
        TYPE(VAR_FILE_ARRAY) :: FileArr
        !
        integer NDX, NDY
        integer, allocatable :: cf_a(:), cf_b(:), cf_c(:), cf_d(:), cf_J(:)
        integer, allocatable :: cf_c11(:,:), cf_c12(:,:), cf_c21(:,:), cf_c22(:,:)
        !
        integer idx1
        integer idx2
        REAL(rkind), allocatable :: RawRead1(:,:)
        REAL(rkind), allocatable :: RawRead2(:,:)
      END TYPE VAR_NETCDF_CF
      TYPE(VAR_NETCDF_CF) :: VAR_wind, VAR_pres
      !
      integer, dimension(:), pointer :: oned_send_rqst
      integer, dimension(:,:), pointer :: oned_send_stat
      integer, dimension(:), pointer :: oned_send_type
      !
      integer, allocatable :: ListNP(:)
      integer, allocatable :: ListIPLG(:)
      integer, allocatable :: ListNPA(:)
      !
      REAL(rkind), allocatable :: outvar1(:,:)
      REAL(rkind), allocatable :: outvar2(:,:)
      REAL(rkind),  PARAMETER  :: DAY2SEC  = 86400.d0
      REAL(rkind),  PARAMETER  :: SEC2DAY  = 1.d0/DAY2SEC
      INTEGER istat
      REAL(rkind) eStartSim
    CONTAINS
      SUBROUTINE GET_NETCDF_VARNAME(string1, string2)
      IMPLICIT NONE
      character(len=140), intent(in) :: string1
      character(len=140), intent(out) :: string2
      string2=TRIM(string1)
      IF (TRIM(string1) == 'WIND10') THEN
        string2='Uwind'
      END IF
      IF (TRIM(string1) == 'SurfCurr') THEN
        string2='UsurfCurr'
      END IF
      END SUBROUTINE
!**********************************************************************
!*  Reading the type of variable                                      *
!*  ---the accessible variable is read from the generic name          *
!**********************************************************************
      SUBROUTINE SCATTER_ONED_ARRAY(Vtotal, Vlocal)
      IMPLICIT NONE
      real(rkind), intent(in) :: Vtotal(np_global)
      real(rkind), intent(out) :: Vlocal(npa)
      integer iProc, IP
      IF (myrank .eq. 0) THEN
        DO iProc=2,nproc
          CALL mpi_isend(Vtotal, 1, oned_send_type(iProc-1), iProc-1, 2030, comm, oned_send_rqst(iProc-1), ierr)
        END DO
        DO IP=1,npa
          Vlocal(IP)=Vtotal(iplg(IP))
        END DO
        IF (nproc > 1) THEN
          CALL MPI_WAITALL(nproc-1, oned_send_rqst, oned_send_stat, ierr)
        END IF
      ELSE
        CALL MPI_RECV(Vlocal, npa, rtype, 0, 2030, comm, istatus, ierr)
      END IF
      END SUBROUTINE
!****************************************************************************
!*                                                                          *
!****************************************************************************
      SUBROUTINE SETUP_ONED_SCATTER_ARRAY
      IMPLICIT NONE
      integer :: ListFirst(nproc)
      integer NPAloc, iProc, IP, IP_glob
      integer, allocatable :: dspl_oned(:), dspl_twod(:)
      integer mpiStatSize
      ListFirst=0
      DO iProc=2,nproc
        ListFirst(iProc)=ListFirst(iProc-1) + ListNPA(iProc-1)
      END DO
      mpiStatSize = size(istatus)
      IF (myrank .eq. 0) THEN
        allocate(oned_send_rqst(nproc-1), oned_send_stat(mpiStatSize,nproc-1), oned_send_type(nproc-1), stat=istat)
        DO iProc=2,nproc
          NPAloc=ListNPA(iProc)
          allocate(dspl_oned(NPAloc), stat=istat)
          DO IP=1,NPAloc
            IP_glob=ListIPLG(IP+ListFirst(iProc))
            dspl_oned(IP)=IP_glob-1
          END DO
          call mpi_type_create_indexed_block(NPAloc,1,dspl_oned,rtype,oned_send_type(iProc-1), ierr)
          call mpi_type_commit(oned_send_type(iProc-1), ierr)
          deallocate(dspl_oned)
        END DO
      END IF
      END SUBROUTINE
!****************************************************************************
!* Specific routines for data output                                        *
!****************************************************************************
      SUBROUTINE COLLECT_ALL_IPLG
      implicit none
      integer, allocatable :: rbuf_int(:)
      integer lenWork, iProc, IP, idx, sumNPA
      allocate(ListNPA(nproc), ListNP(nproc), rbuf_int(2), stat=istat)
      IF (istat/=0) CALL parallel_abort('wwm_aux_parall, allocate error 12')
      IF (myrank == 0) THEN
        ListNPA(1)=npa
        ListNP(1) =np
        DO iProc=2,nproc
          CALL MPI_RECV(rbuf_int,2,itype, iProc-1, 257, comm, istatus, ierr)
          ListNPA(iProc)=rbuf_int(1)
          ListNP(iProc)=rbuf_int(2)
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(ListNPA,nproc,itype, iProc-1, 263, comm, ierr)
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(ListNP,nproc,itype, iProc-1, 571, comm, ierr)
        END DO
      ELSE
        rbuf_int(1)=npa
        rbuf_int(2)=np
        CALL MPI_SEND(rbuf_int,2,itype, 0, 257, comm, ierr)
        CALL MPI_RECV(ListNPA,nproc,itype, 0, 263, comm, istatus, ierr)
        CALL MPI_RECV(ListNP,nproc,itype, 0, 571, comm, istatus, ierr)
      END IF
      deallocate(rbuf_int)
      sumNPA=sum(ListNPA)
      allocate(ListIPLG(sumNPA), stat=istat)
      IF (istat/=0) CALL parallel_abort('wwm_aux_parall, allocate error 13')
      IF (myrank == 0) THEN
        idx=0
        DO IP=1,npa
          idx=idx+1
          ListIPLG(idx)=iplg(IP)
        END DO
        DO iProc=2,nproc
          lenWork=ListNPA(iProc)
          allocate(rbuf_int(lenWork), stat=istat)
          IF (istat/=0) CALL parallel_abort('wwm_aux_parall, allocate error 14')
          CALL MPI_RECV(rbuf_int,lenWork,itype, iProc-1, 269, comm, istatus, ierr)
          DO IP=1,lenWork
            idx=idx+1
            ListIPLG(idx)=rbuf_int(IP)
          END DO
          deallocate(rbuf_int)
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(ListIPLG,sumNPA,itype, iProc-1, 271, comm, ierr)
        END DO
      ELSE
        CALL MPI_SEND(iplg,npa,itype, 0, 269, comm, ierr)
        CALL MPI_RECV(ListIPLG,sumNPA,itype, 0, 271, comm, istatus, ierr)
      END IF
      END SUBROUTINE
!****************************************************************************
!*                                                                          *
!****************************************************************************
      SUBROUTINE DATE2JD(year, month, day, hour, min, sec, eJD)
      IMPLICIT NONE
      integer, intent(in) :: year, month, day, hour, min, sec
      real(rkind), intent(out) :: eJD
      real(rkind) :: eJDbase, eFracDay
      integer a, y, m
      a = floor((DBLE(14) - DBLE(month))/DBLE(12))
      y = year + 4800 - a
      m = month + 12*a - 3
      ! For a date in the Gregorian calendar:                                                                   
      eJDbase = DBLE(day)                                            &
     & + DBLE(floor((DBLE(153)*DBLE(m) + DBLE(2))/DBLE(5)))  &
     & + DBLE(y)*DBLE(365)                                         &
     & + DBLE(floor(DBLE(y)/DBLE(4)))                            &
     & - DBLE(floor(DBLE(y)/DBLE(100)))                          &
     & + DBLE(floor(DBLE(y)/DBLE(400))) - DBLE(32045)
      eFracDay=(DBLE(sec) +                                          &
     &          DBLE(60)*DBLE(min) +                               &
     &          DBLE(3600)*(DBLE(hour) - DBLE(12))               &
     &          )/DBLE(86400)
      eJD=eJDbase + eFracDay
      END SUBROUTINE
!****************************************************************************
!*                                                                          *
!****************************************************************************
      SUBROUTINE DATE_ConvertSix2mjd(year, month, day, hour, min, sec, eMJD)
      IMPLICIT NONE
      integer, intent(in) :: year, month, day, hour, min, sec
      real(rkind), intent(out) :: eMJD
      real(rkind) :: eJD1, eJD2
      CALL DATE2JD(year, month, day, hour, min, sec, eJD1)
      CALL DATE2JD(1858, 11, 17, 0, 0, 0, eJD2)
      eMJD=eJD1-eJD2
      END SUBROUTINE
!****************************************************************************
!*                                                                          *
!****************************************************************************
      SUBROUTINE DATE_ConvertString2six(year, month, day, hour, min, sec, eTimeStr)
      IMPLICIT NONE
      integer, intent(out) :: year, month, day, hour, min, sec
      character(len=15), intent(in) :: eTimeStr
      character(len=4) eYear
      character(len=2) eMonth, eDay, eHour, eMin, eSec
      eYear(1:1)  = eTimeStr(1:1)
      eYear(2:2)  = eTimeStr(2:2)
      eYear(3:3)  = eTimeStr(3:3)
      eYear(4:4)  = eTimeStr(4:4)
      eMonth(1:1) = eTimeStr(5:5)
      eMonth(2:2) = eTimeStr(6:6)
      eDay(1:1)   = eTimeStr(7:7)
      eDay(2:2)   = eTimeStr(8:8)
      eHour(1:1)  = eTimeStr(10:10)
      eHour(2:2)  = eTimeStr(11:11)
      eMin(1:1)   = eTimeStr(12:12)
      eMin(2:2)   = eTimeStr(13:13)
      eSec(1:1)   = eTimeStr(14:14)
      eSec(2:2)   = eTimeStr(15:15)
      read(eYear , '(i10)' ) year
      read(eMonth, '(i10)' ) month
      read(eDay  , '(i10)' ) day
      read(eHour , '(i10)' ) hour
      read(eMin  , '(i10)' ) min
      read(eSec  , '(i10)' ) sec
      END SUBROUTINE
!****************************************************************************
!*                                                                          *
!****************************************************************************
      SUBROUTINE DATE_ConvertSix2string(year, month, day, hour, min, sec, eTimeStr)
      IMPLICIT NONE
      integer, intent(in) :: year, month, day, hour, min, sec
      character(len=15), intent(out) :: eTimeStr
      WRITE(eTimeStr, 20) year, month, day, hour, min, sec
  20  FORMAT (i4.4, i2.2, i2.2, '.', i2.2, i2.2, i2.2)
      END SUBROUTINE
!****************************************************************************
!*                                                                          *
!****************************************************************************
      SUBROUTINE MONTH_LEN(year, month, lenmonth)
      IMPLICIT NONE
      integer, intent(in) :: year, month
      integer, intent(out) :: lenmonth
      IF ((month .eq. 1).or.(month .eq. 3).or.(month .eq. 5).or.(month .eq. 7).or.(month .eq. 8).or.(month .eq. 10).or.(month .eq. 12)) THEN
        lenmonth=31
      END IF
      IF ((month .eq. 4).or.(month .eq. 6).or.(month .eq. 9).or.(month .eq. 11)) THEN
        lenmonth=30
      END IF
      IF (month .eq. 2) THEN
        IF (MOD(year, 4) .ne. 0) THEN
          lenmonth=28
        ELSE
          IF (MOD(year, 100) .ne. 0) THEN
            lenmonth=29
          ELSE
            IF (MOD(year, 400) .ne. 0) THEN
              lenmonth=28
            ELSE
              lenmonth=29
            END IF
          END IF
        END IF
      END IF
      END SUBROUTINE
!****************************************************************************
!*                                                                          *
!****************************************************************************
      SUBROUTINE JD2DATE(year, month, day, hour, min, sec, eJD)
      ! The following algorithm is from the Calendar FAQ.
      IMPLICIT NONE
      integer, intent(out) :: year, month, day, hour, min, sec
      real(rkind), intent(in) :: eJD
      integer ijd, a, b, c, d, e, m
      integer secNear, lenmonth
      real(rkind) :: fjd, second
      ijd = floor(eJD + 0.5_rkind)
      !
      a = ijd + 32044
      b = floor((DBLE(4)*DBLE(a) + DBLE(3)) / DBLE(146097))
      c = a - floor((DBLE(b) * DBLE(146097)) / DBLE(4))
      !
      d = floor((DBLE(4)*DBLE(c) + DBLE(3)) / DBLE(1461))
      e = c - floor((DBLE(1461)*DBLE(d)) / DBLE(4))
      m = floor((DBLE(5) * DBLE(e) + DBLE(2)) / DBLE(153))
      !
      day   = e - floor((DBLE(153) * DBLE(m) + DBLE(2)) / DBLE(5)) + 1
      month = m + 3 - 12 * floor(DBLE(m) / DBLE(10))
      year  = b * 100 + d - 4800 + floor(DBLE(m) / DBLE(10))
      !
      fjd    = eJD - DBLE(ijd) + 0.5_rkind
      second = DBLE(86400) * fjd
      hour   = floor(second/DBLE(3600))
      second = second - DBLE(3600)*DBLE(hour)
      min    = floor(second/DBLE(60))
      sec    = floor(second - DBLE(60)*min)
      !
      ! Now renormalizing
      !
      secNear=NINT(second - DBLE(60)*min)
      IF (secNear .eq. 60) THEN
        sec=0
        min=min+1
      END IF
      IF (min .eq. 60) THEN
        min=0
        hour=hour+1
      END IF
      IF (hour .eq. 24) THEN
        hour=0
        day=day+1
      END IF
      CALL MONTH_LEN(year, month, lenmonth)
      IF (day .eq. lenmonth+1) THEN
        day=1
        month=month+1
      END IF
      IF (month .eq. 13) THEN
        month=1
        year=year+1
      END IF
      END SUBROUTINE
!****************************************************************************
!*                                                                          *
!****************************************************************************
      SUBROUTINE CT2MJD(STIME,XMJD)
      IMPLICIT NONE
      CHARACTER(LEN=15), INTENT(IN) :: STIME
      real(rkind), INTENT(OUT) :: XMJD
      integer year, month, day, hour, min, sec
      CALL DATE_ConvertString2six(year, month, day, hour, min, sec, STIME)
      CALL DATE_ConvertSix2mjd(year, month, day, hour, min, sec, XMJD)
      END SUBROUTINE
!****************************************************************************
!*                                                                          *
!****************************************************************************
      SUBROUTINE MJD2CT(XMJD,STIME)
      IMPLICIT NONE
      CHARACTER(LEN=15), INTENT(OUT) :: STIME
      real(rkind), INTENT(IN) :: XMJD
      integer year, month, day, hour, min, sec
      real(rkind) XMJD_1858, eMJD
      CALL DATE2JD(1858, 11, 17, 0, 0, 0, XMJD_1858)
      eMJD = XMJD + XMJD_1858
      CALL JD2DATE(year, month, day, hour, min, sec, eMJD)
      CALL DATE_ConvertSix2string(year, month, day, hour, min, sec, STIME)
      END SUBROUTINE
!****************************************************************************
!*                                                                          *
!****************************************************************************
      SUBROUTINE READ_LIST_TIME(eVAR_TIME, eFileName, eString)
      IMPLICIT NONE
      TYPE(VAR_LIST_TIME), intent(inout) :: eVAR_TIME
      character(len=140), intent(in) :: eFileName
      character(len=140), intent(in) :: eString
      !
      character(len=140) :: eStringCF
      INTEGER           :: fid, varid
      integer, dimension(nf90_max_var_dims) :: dimidsB
      integer, dimension(nf90_max_var_dims) :: dimids
      character (len=20) :: TimeStr
      character (len=100) :: eStrUnitTime
      real(rkind) :: ConvertToDay
      real(rkind) :: eTimeStart
      real(rkind), allocatable :: PreListTime(:)
      character (len = *), parameter :: CallFct="READ_LIST_TIME"
      integer nbtime_mjd
      integer nbChar
      integer len1, len2, len3
      integer nbDim
      !
      CALL GET_NETCDF_VARNAME(eString, eStringCF)
      !
      ISTAT = nf90_open(TRIM(eFileName), nf90_nowrite, fid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)
      WRITE(16,*) 'eFileName=', TRIM(eFileName)
      FLUSH(16)
      ! Reading wind attributes

      ISTAT = nf90_inq_varid(fid, TRIM(eStringCF), varid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)
      WRITE(16,*) 'eStringCF=', TRIM(eStringCF)
      FLUSH(16)

      ISTAT = nf90_inquire_variable(fid, varid, dimids=dimidsB)
      CALL GENERIC_NETCDF_ERROR(CallFct, 3, ISTAT)

      ISTAT = nf90_inquire_variable(fid, varid, ndims=nbDim)
      CALL GENERIC_NETCDF_ERROR(CallFct, 3, ISTAT)
!      Print *, 'nbDim=', nbDim

!      ISTAT = nf90_inquire_dimension(fid, dimidsB(1), len=len1)
!      CALL GENERIC_NETCDF_ERROR(CallFct, 4, ISTAT)
!      ISTAT = nf90_inquire_dimension(fid, dimidsB(2), len=len2)
!      CALL GENERIC_NETCDF_ERROR(CallFct, 4, ISTAT)
!      ISTAT = nf90_inquire_dimension(fid, dimidsB(3), len=len3)
!      CALL GENERIC_NETCDF_ERROR(CallFct, 4, ISTAT)
!      Print *, 'len123=', len1, len2, len3

      ISTAT = nf90_inquire_dimension(fid, dimidsB(nbDim), name=TimeStr)
      CALL GENERIC_NETCDF_ERROR(CallFct, 4, ISTAT)

      ! Reading time

      ISTAT = nf90_inq_varid(fid, TimeStr, varid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 5, ISTAT)

      ISTAT = nf90_inquire_attribute(fid, varid, "units", len=nbChar)
      CALL GENERIC_NETCDF_ERROR(CallFct, 6, ISTAT)
      ISTAT = nf90_get_att(fid, varid, "units", eStrUnitTime)

      CALL GENERIC_NETCDF_ERROR(CallFct, 7, ISTAT)
      CALL CF_EXTRACT_TIME(eStrUnitTime, ConvertToDay, eTimeStart)

      ISTAT = nf90_inquire_variable(fid, varid, dimids=dimids)
      CALL GENERIC_NETCDF_ERROR(CallFct, 8, ISTAT)

      ISTAT = nf90_inquire_dimension(fid, dimids(1), len=nbtime_mjd)
      CALL GENERIC_NETCDF_ERROR(CallFct, 9, ISTAT)

      eVAR_TIME % nbTime = nbtime_mjd
      allocate(eVAR_TIME % ListTime(nbtime_mjd), PreListTime(nbtime_mjd), stat=istat)

      ISTAT = nf90_get_var(fid, varid, PreListTime)
      CALL GENERIC_NETCDF_ERROR(CallFct, 10, ISTAT)

      ISTAT = nf90_close(fid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 11, ISTAT)

      eVAR_TIME % ListTime = PreListTime(:)*ConvertToDay + eTimeStart
      deallocate(PreListTime)
      END SUBROUTINE
!****************************************************************************
!*                                                                          *
!****************************************************************************
      SUBROUTINE READ_LIST_SERIES_TIME(eVAR_SERIES, VAR_FILE, eString)
      IMPLICIT NONE
      TYPE(VAR_LIST_SERIES_TIME), intent(inout) :: eVAR_SERIES
      TYPE(VAR_FILE_ARRAY), intent(in) :: VAR_FILE
      character(len=140), intent(in) :: eString
      !
      TYPE(VAR_LIST_TIME), allocatable :: ListVar_List(:)
      integer idx, ifile, nbTimeTot, nbTime, iTime, nbFile
      REAL(rkind) DeltaTime, eTime1, eTime2, eTime
      !
      ! Read all the files
      !
      nbFile=VAR_FILE % nbFile
      IF (nbFile .eq. 0) THEN
        Print *, 'We have zero file.'
        Print *, 'Error in READ_LIST_SERIES_TIME'
        CALL parallel_abort('Error here')
      END IF
      allocate(ListVar_List(nbFile))
      CALL READ_LIST_TIME(ListVar_List(1), VAR_FILE % ListFile(1), eString)
      IF (nbFile .ge. 2) THEN
        CALL READ_LIST_TIME(ListVar_List(2), VAR_FILE % ListFile(2), eString)
      END IF
      IF (nbFile .ge. 3) THEN
        CALL READ_LIST_TIME(ListVar_List(3), VAR_FILE % ListFile(3), eString)
      END IF
      IF (nbFile .ge. 4) THEN
        CALL READ_LIST_TIME(ListVar_List(nbFile), VAR_FILE % ListFile(nbFile), eString)
      END IF
      IF (nbFile .gt. 4) THEN
        nbTime=ListVar_List(2) % nbTime
        eTime1 = ListVar_List(2) % ListTime(1)
        eTime2 = ListVar_List(3) % ListTime(1)
        DeltaTime = eTime2 - eTime1
        DO ifile=4,nbFile-1
          ListVar_List(ifile) % nbTime = nbTime
          allocate(ListVar_List(ifile) % ListTime(nbTime))
          DO iTime=1,nbTime
            eTime=ListVar_List(ifile-1) % ListTime(iTime) + DeltaTime
            ListVar_List(ifile) % ListTime(iTime) = eTime
          END DO
        END DO
      END IF
      !
      ! Read all the files
      !
      nbTimeTot=0
      DO ifile=1,nbFile
        nbTimeTot = nbTimeTot + ListVar_List(ifile) % nbTime
      END DO
      eVAR_SERIES % nbTime = nbTimeTot
      allocate(eVAR_SERIES % List_IFILE(nbTimeTot))
      allocate(eVAR_SERIES % List_ITIME(nbTimeTot))
      allocate(eVAR_SERIES % List_Time(nbTimeTot))
      idx=0
      DO ifile=1,nbFile
        nbTime=ListVar_List(ifile) % nbTime
        DO iTime=1,nbTime
          idx=idx+1
          eVAR_SERIES % List_IFILE(idx) = ifile
          eVAR_SERIES % List_ITIME(idx) = iTime
          eVAR_SERIES % List_Time (idx) = ListVar_List(ifile) % ListTime(iTime)
        END DO
      END DO
      END SUBROUTINE
!****************************************************************************
!*                                                                          *
!****************************************************************************
      SUBROUTINE READ_OFFSET_SCALE_FACTOR(eFileName, cf_scale_factor, cf_add_offset, eString)
      IMPLICIT NONE
      character(len=140), intent(in) :: eFileName
      REAL(rkind), intent(inout) :: cf_scale_factor
      REAL(rkind), intent(inout) :: cf_add_offset
      character(len=140), intent(in) :: eString
      !
      character(len=140) :: eStringCF
      character (len = *), parameter :: CallFct="READ_OFFSET_SCALE_FACTOR"
      character(len=100) :: CHRERR
      integer fid, varid
      character (len=20) :: TimeStr
      integer, dimension(nf90_max_var_dims) :: dimidsB
      !
      CALL GET_NETCDF_VARNAME(eString, eStringCF)
      !
      ISTAT = nf90_open(TRIM(eFileName), nf90_nowrite, fid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)

      ! Reading wind attributes

      ISTAT = nf90_inq_varid(fid, TRIM(eStringCF), varid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)

      ISTAT = nf90_inquire_variable(fid, varid, dimids=dimidsB)
      CALL GENERIC_NETCDF_ERROR(CallFct, 3, ISTAT)

      ISTAT = nf90_inquire_dimension(fid, dimidsB(2), name=TimeStr)
      CALL GENERIC_NETCDF_ERROR(CallFct, 4, ISTAT)

      ISTAT = nf90_get_att(fid, varid, "scale_factor", cf_scale_factor)
      IF (ISTAT /= 0) THEN
        CHRERR = nf90_strerror(ISTAT)
        cf_scale_factor=1.0_rkind
      ENDIF
      ISTAT = nf90_get_att(fid, varid, "add_offset", cf_add_offset)
      IF (ISTAT /= 0) THEN
        CHRERR = nf90_strerror(ISTAT)
        cf_add_offset=0.0_rkind
      ENDIF

      ISTAT = nf90_close(fid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 5, ISTAT)
      END SUBROUTINE
!****************************************************************************
!*                                                                          *
!****************************************************************************
      SUBROUTINE GET_INTERPOLATION_ARRAY(eVAR, eFileName, eString)
      IMPLICIT NONE
      TYPE(VAR_NETCDF_CF), intent(inout) :: eVAR
      character(len=140), intent(in) :: eString
      character(len=140), intent(in) :: eFileName
      !
      character(len=100) :: StrCoord
      character(len=100) :: Yname1, Yname2
      integer len1, posBlank
      character (len = *), parameter :: CallFct="GET_INTERPOLATION_ARRAY"
      REAL(rkind), allocatable :: CF_LON(:,:), CF_LAT(:,:)
      REAL(rkind), allocatable :: dist(:,:)
      character(len=140) :: eStringCF
      integer NDX, NDY
      REAL(rkind) d_lat, d_lon, closest(2)
      integer i11, i12, i21, i22, j11, j12, j21
      integer I
      integer fid, varid
      integer, dimension(nf90_max_var_dims) :: dimids
      !
      CALL GET_NETCDF_VARNAME(eString, eStringCF)
      !
      ISTAT = nf90_open(TRIM(eFileName), nf90_nowrite, fid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)
      !
      ! Reading variable
      !
      ISTAT = nf90_inq_varid(fid, TRIM(eStringCF), varid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)
      !
      ISTAT = nf90_get_att(fid, varid, "coordinates", StrCoord)
      CALL GENERIC_NETCDF_ERROR(CallFct, 3, ISTAT)
      !
      len1=LEN_TRIM(StrCoord)
      posBlank=INDEX(StrCoord(1:len1), ' ')
      Yname1 = StrCoord(1:posBlank-1)
      Yname2 = StrCoord(posBlank+1:len1)
      !
      ISTAT = nf90_inq_varid(fid, TRIM(Yname1), varid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 4, ISTAT)
      !
      ISTAT = nf90_inquire_variable(fid, varid, dimids=dimids)
      CALL GENERIC_NETCDF_ERROR(CallFct, 5, ISTAT)
      !
      ISTAT = nf90_inquire_dimension(fid, dimids(1), len=NDX)
      CALL GENERIC_NETCDF_ERROR(CallFct, 6, ISTAT)
      !
      ISTAT = nf90_inquire_dimension(fid, dimids(2), len=NDY)
      CALL GENERIC_NETCDF_ERROR(CallFct, 7, ISTAT)
      eVAR%NDX = NDX
      eVAR%NDY = NDY
      !
      allocate(CF_LON(NDX, NDY), CF_LAT(NDX, NDY), stat=istat)
      !
      ! loading LON
      !
      ISTAT = nf90_get_var(fid, varid, CF_LON)
      CALL GENERIC_NETCDF_ERROR(CallFct, 8, ISTAT)
      !
      ! loading LAT
      !
      ISTAT = nf90_inq_varid(fid, TRIM(Yname2), varid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 9, ISTAT)
      !
      ISTAT = nf90_get_var(fid, varid, CF_LAT)
      CALL GENERIC_NETCDF_ERROR(CallFct, 10, ISTAT)
      !
      ISTAT = nf90_close(fid)
      CALL GENERIC_NETCDF_ERROR(CallFct, 11, ISTAT)

      CF_LON(:,:)=CF_LON(:,:)*pi/180.0_rkind
      CF_LAT(:,:)=CF_LAT(:,:)*pi/180.0_rkind

      if(myrank==0) then
        WRITE(16,*) 'NDX_WIND_FD=', NDX
        WRITE(16,*) 'NYX_WIND_FD=', NDY
        WRITE(16,*) 'Min/Max CF_LON',minval(CF_LON),maxval(CF_LON)
        WRITE(16,*) 'Min/Max CF_LAT',minval(CF_LAT),maxval(CF_LAT)
        WRITE(16,*) 'Min/Max xlon',minval(xlon),maxval(xlon)
        WRITE(16,*) 'Min/Max ylat',minval(ylat),maxval(ylat)
        WRITE(16,*) 'Done with UVP.nc init phase, calculating interp coefs'
        FLUSH(16)
      endif
      !
! compute nodes and coefs for bilinear interpolation for whole grid
      ALLOCATE(eVAR%cf_c11(npa,2), eVAR%cf_c12(npa,2), eVAR%cf_c21(npa,2), eVAR%cf_c22(npa,2), stat=istat)
      IF (istat/=0) WRITE(16,*) 'Problem with allocate cf11,12,21,22'
      ALLOCATE(eVAR%cf_a(npa), eVAR%cf_b(npa), eVAR%cf_c(npa), eVAR%cf_d(npa), eVAR%cf_J(npa), stat=istat)
      IF (istat/=0) WRITE(16,*) 'Problem with allocate cf_a,b,c,f,J'
      ALLOCATE(dist(NDX, NDY), stat=istat)
      IF (istat/=0) WRITE(16,*) 'Problem with allocate dist'

      DO I = 1, npa
        dist(:,:) = ABS( CMPLX(xlon(I)-CF_LON(:,:), ylat(I)-CF_LAT(:,:)) )
        closest(1:2) = MINLOC(dist)
        d_lon = xlon(I)-CF_LON(closest(1),closest(2)) 
        d_lat = ylat(I)-CF_LAT(closest(1),closest(2))
        IF ((d_lon.ge.0).and.(d_lat.ge.0)) THEN ! point is in the I kvadrant
            eVAR%cf_c11(I,:) = closest(:)
            eVAR%cf_c21(I,1) = closest(1) + 1
            eVAR%cf_c22(I,1) = closest(1) + 1
            eVAR%cf_c12(I,1) = closest(1)
            eVAR%cf_c21(I,2) = closest(2)
            eVAR%cf_c22(I,2) = closest(2) + 1
            eVAR%cf_c12(I,2) = closest(2) + 1
        END IF
        IF ((d_lon.ge.0).and.(d_lat.le.0)) THEN ! point is in the IV kvadrant
            eVAR%cf_c11(I,1) = closest(1)
            eVAR%cf_c21(I,1) = closest(1) + 1
            eVAR%cf_c22(I,1) = closest(1) + 1
            eVAR%cf_c12(I,:) = closest(:)
            eVAR%cf_c11(I,2) = closest(2) - 1
            eVAR%cf_c21(I,2) = closest(2) - 1
            eVAR%cf_c22(I,2) = closest(2) 
        END IF
        IF ((d_lon.le.0).and.(d_lat.ge.0)) THEN ! point is in the II kvadrant
            eVAR%cf_c11(I,1) = closest(1) - 1 
            eVAR%cf_c21(I,:) = closest(:)
            eVAR%cf_c22(I,1) = closest(1)
            eVAR%cf_c12(I,1) = closest(1) - 1
            eVAR%cf_c11(I,2) = closest(2)
            eVAR%cf_c22(I,2) = closest(2) + 1
            eVAR%cf_c12(I,2) = closest(2) + 1 
        END IF
        IF ((d_lon.le.0).and.(d_lat.le.0)) THEN ! point is in the III kvadrant
            eVAR%cf_c11(I,1) = closest(1) - 1
            eVAR%cf_c21(I,1) = closest(1)
            eVAR%cf_c22(I,:) = closest(:)
            eVAR%cf_c12(I,1) = closest(1) - 1
            eVAR%cf_c11(I,2) = closest(2) - 1
            eVAR%cf_c21(I,2) = closest(2) - 1
            eVAR%cf_c12(I,2) = closest(2) 
        END IF
        ! J =1/((x2-x1)*(y2-y1))
        i11=eVAR%cf_c11(I,1)
        j11=eVAR%cf_c11(I,2)
        i12=eVAR%cf_c12(I,1)
        j12=eVAR%cf_c12(I,2)
        i21=eVAR%cf_c21(I,1)
        j21=eVAR%cf_c21(I,2)
        eVAR%cf_J(I)=1.0/( (CF_LON(i21,j21)-CF_LON(i11,j11))*(CF_LAT(i12,j12)-CF_LAT(i11,j11)) )
        eVAR%cf_a(I) = CF_LON(i21,j21) - xlon(I) ! x2-x
        eVAR%cf_b(I) = xlon(I) - CF_LON(i11,j11) ! x-x1
        eVAR%cf_c(I) = CF_LAT(i12,j12) - ylat(I) ! y2-y
        eVAR%cf_d(I) = ylat(I) - CF_LAT(i11,j11) ! y-y1
      END DO
      DEALLOCATE(dist, CF_LON, CF_LAT)
      if(myrank==0) then
        WRITE(16,*) 'Done with UVP.nc interp coefs, all done in INIT_NETCDF_CF'
        FLUSH(16)
      endif
      END SUBROUTINE
!****************************************************************************
!*                                                                          *
!****************************************************************************
      SUBROUTINE READ_FILE_ARRAY(eVAR_FILE, ePrefix_File)
      IMPLICIT NONE
      TYPE(VAR_FILE_ARRAY), intent(inout) :: eVAR_FILE
      CHARACTER(LEN=140), intent(in) :: ePrefix_File
      !
      LOGICAL LFLIVE
      character(len=140) :: eFILE, eSuffix
      integer ifile, nbFile
      INQUIRE( FILE = TRIM(ePrefix_File), EXIST = LFLIVE )
      IF (LFLIVE) THEN
        eVAR_FILE % nbFile = 1
        allocate(eVAR_FILE % ListFile(1))
        eVAR_FILE % ListFile(1) = ePrefix_File
        RETURN
      END IF
      ifile=1
      nbFile=0
      DO
!        Print *, 'ifile=', ifile
!        Print *, 'ePrefix_File=', TRIM(ePrefix_File)
        WRITE (eSuffix,10) ifile
!        Print *, 'eSuffix=', TRIM(eSuffix)
        eFILE=TRIM(ePrefix_File) // eSuffix
!        Print *, 'eFILE=', TRIM(eFILE)
        INQUIRE( FILE = TRIM(eFILE), EXIST = LFLIVE )
        IF (.NOT. LFLIVE) THEN
          EXIT
        END IF
        nbFile=ifile
        ifile = ifile + 1
      END DO
      IF (nbFile .eq. 0) THEN
        Print *, 'We have zero file error'
        Print *, 'ePrefix_File=', TRIM(ePrefix_File)
        CALL parallel_abort('Error')
      END IF
      eVAR_FILE % nbFile = nbFile
      allocate(eVAR_FILE % ListFile(nbFile))
      DO ifile=1,nbFile
        WRITE (eSuffix,10) ifile
        eFILE=TRIM(ePrefix_File) // eSuffix
        INQUIRE( FILE = TRIM(eFILE), EXIST = LFLIVE )
        Print *, 'ifile=', ifile, ' eFILE=', TRIM(eFILE)
        eVAR_FILE % ListFile(ifile) = eFILE
      END DO
      RETURN
10    FORMAT (i4.4,'.nc')
      END SUBROUTINE
!****************************************************************************
!*                                                                          *
!****************************************************************************
      SUBROUTINE INIT_NETCDF_CF(eVAR, MULTIPLE_IN, ePrefix_File, eString, nws)
      IMPLICIT NONE
      TYPE(VAR_NETCDF_CF), intent(inout) :: eVAR
      logical, intent(in) :: MULTIPLE_IN
      character(len=140), intent(in) :: ePrefix_File
      character(len=140), intent(in) :: eString
      integer, intent(in) :: nws
      !
      real(rkind) :: cf_scale_factor, cf_add_offset
      integer eInt(1)
      real(rkind) :: eReal(2)
      integer IPROC
      integer dimVar
      LOGICAL IsDirect
      !
      CALL READ_FILE_ARRAY(eVAR % FileArr, ePrefix_File)
      eVAR % MULTIPLE_IN = MULTIPLE_IN
      IF (MULTIPLE_IN .or. (myrank .eq. 0)) THEN
        CALL READ_LIST_SERIES_TIME(eVAR % VAR_SERIES, eVAR % FileArr, eString)
        CALL READ_OFFSET_SCALE_FACTOR(eVAR % FileArr % ListFile(1), cf_scale_factor, cf_add_offset, eString)
      END IF
      !
      IF (nws == 5) THEN
        IsDirect=.FALSE.
      ELSE
        IsDirect=.TRUE.
      END IF
      IF (IsDirect .eqv. .FALSE.) THEN
        CALL GET_INTERPOLATION_ARRAY(eVAR, eVAR % FileArr % ListFile(1), eString)
      END IF
      eVAR % IsDirect = IsDirect
      !
      eVAR % cf_scale_factor = cf_scale_factor
      eVAR % cf_add_offset = cf_add_offset
      IF (.NOT. MULTIPLE_IN) THEN
        IF (myrank .eq. 0) THEN
          eInt(1)=eVAR % VAR_SERIES % nbTime
          DO IPROC=2,nproc
            CALL MPI_SEND(eInt,1,itype, iProc-1, 811, comm, ierr)
          END DO
          DO IPROC=2,nproc
            CALL MPI_SEND(eVAR % VAR_SERIES % List_IFILE,eInt(1),itype, iProc-1, 811, comm, ierr)
            CALL MPI_SEND(eVAR % VAR_SERIES % List_ITIME,eInt(1),itype, iProc-1, 812, comm, ierr)
            CALL MPI_SEND(eVAR % VAR_SERIES % List_Time,eInt(1),rtype, iProc-1, 813, comm, ierr)
          END DO
          eReal(1)=cf_scale_factor
          eReal(2)=cf_add_offset
          DO IPROC=2,nproc
            CALL MPI_SEND(eReal,2,rtype, iProc-1, 814, comm, ierr)
          END DO
        ELSE
          CALL MPI_RECV(eInt,1,itype, 0, 811, comm, istatus, ierr)
          eVAR % VAR_SERIES % nbTime=eInt(1)
          allocate(eVAR % VAR_SERIES % List_IFILE(eInt(1)))
          allocate(eVAR % VAR_SERIES % List_ITIME(eInt(1)))
          allocate(eVAR % VAR_SERIES % List_Time (eInt(1)))
          CALL MPI_RECV(EVAR % VAR_SERIES % List_IFILE, eInt(1),itype, 0, 811, comm, istatus, ierr)
          CALL MPI_RECV(EVAR % VAR_SERIES % List_ITIME, eInt(1),itype, 0, 812, comm, istatus, ierr)
          CALL MPI_RECV(EVAR % VAR_SERIES % List_Time , eInt(1),rtype, 0, 813, comm, istatus, ierr)
          !
          CALL MPI_RECV(eReal,2,rtype, 0, 814, comm, istatus, ierr)
          cf_scale_factor=eReal(1)
          cf_add_offset=eReal(2)
        END IF
      END IF
      dimVar=1
      IF (TRIM(eString) == 'WIND10') THEN
        dimVar = 2
      END IF
      IF (TRIM(eString) == 'SurfCurr') THEN
        dimVar = 2
      END IF
      eVAR % cf_scale_factor = cf_scale_factor
      eVAR % cf_add_offset = cf_add_offset
      eVAR % eString = TRIM(eString)
      eVAR % dimVar = dimVar
      !
      eVAR % idx1 = 0
      eVAR % idx2 = 0
      allocate(eVAR % RawRead1(npa,dimVar), eVAR % RawRead2(npa,dimVar), stat=istat)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_INTERP_NETCDF_CF(eVAR, idx, eFileName, outvar)
      IMPLICIT NONE
      TYPE(VAR_NETCDF_CF)                :: eVAR
      INTEGER, INTENT(in)                :: idx
      character(len=140), intent(in)     :: eFileName
      REAL(rkind), INTENT(out)           :: outvar(npa,eVAR % dimVar)
      !
      INTEGER                             :: FID, ID, ISTAT, I, ntimes
      REAL(rkind)			  :: DirRead(eVAR%NDX, eVAR%NDY)
      REAL(rkind)                         :: RawRead(eVAR%NDX, eVAR%NDY,eVAR%dimVar)
      character (len = *), parameter :: CallFct="READ_INTERP_NETCDF_CF"
      REAL(rkind) cf_add_offset, cf_scale_factor
      character(len=10) :: eStr1, eStr2
      !
      ISTAT = NF90_OPEN(TRIM(eFileName), NF90_NOWRITE, FID)
      CALL GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)
      cf_add_offset = eVAR % cf_add_offset
      cf_scale_factor = eVAR % cf_scale_factor
      IF (eVAR % dimVar == 2) THEN
        IF (TRIM(eVAR % eString) == 'WIND10') THEN
          eStr1='Uwind'
          eStr2='Vwind'
        END IF
        IF (TRIM(eVAR % eString) == 'SurfCurr') THEN
          eStr1='UsurfCurr'
          eStr2='VsurfCurr'
        END IF
      ELSE
        eStr1=TRIM(eVAR % eString)
      END IF

      ISTAT = NF90_inq_varid(FID, TRIM(eStr1), ID)
      CALL GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)

      ISTAT = NF90_GET_VAR(FID, ID, DirRead, start = (/ 1, 1, idx /), count = (/ eVAR%NDX, eVAR%NDY, 1 /))
      CALL GENERIC_NETCDF_ERROR(CallFct, 3, ISTAT)
      RawRead(:,:,1)=cf_add_offset + cf_scale_factor*DirRead(:,:)

      IF (eVAR % dimVar == 2) THEN
        ISTAT = NF90_inq_varid(FID, TRIM(eStr2), ID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 4, ISTAT)

        ISTAT = NF90_GET_VAR(FID, ID, DirRead, start = (/ 1, 1, idx /), count = (/ eVAR%NDX, eVAR%NDY, 1 /))
        CALL GENERIC_NETCDF_ERROR(CallFct, 5, ISTAT)
        RawRead(:,:,2)=cf_add_offset + cf_scale_factor*DirRead(:,:)
      END IF

      ISTAT = NF90_CLOSE(FID)
      CALL GENERIC_NETCDF_ERROR(CallFct, 6, ISTAT)

      DO I = 1, npa
        outvar(I,:) =  eVAR%cf_J(I)*(                                &
     &     RawRead(eVAR%cf_c11(I,1),eVAR%cf_c11(I,2),:)*eVAR%cf_a(I)*eVAR%cf_c(I)+        &
     &     RawRead(eVAR%cf_c21(I,1),eVAR%cf_c21(I,2),:)*eVAR%cf_b(I)*eVAR%cf_c(I)+        &
     &     RawRead(eVAR%cf_c12(I,1),eVAR%cf_c12(I,2),:)*eVAR%cf_a(I)*eVAR%cf_d(I)+        &
     &     RawRead(eVAR%cf_c22(I,1),eVAR%cf_c22(I,2),:)*eVAR%cf_b(I)*eVAR%cf_d(I) )
      END DO
      END SUBROUTINE READ_INTERP_NETCDF_CF
!****************************************************************************
!*                                                                          *
!****************************************************************************
      SUBROUTINE READ_DIRECT_NETCDF_CF(eVAR, itime, eFileName, outvar)
      IMPLICIT NONE
      TYPE(VAR_NETCDF_CF), intent(in)    :: eVAR
      INTEGER, INTENT(in)                :: itime
      character(len=140), intent(in)     :: eFileName
      REAL(rkind), INTENT(out)           :: outvar(npa,eVAR % dimVar)
      !
      character (len = *), parameter :: CallFct="READ_DIRECT_NETCDF_CF"
      INTEGER                            :: FID, varid
      real(rkind) :: DirectRead(np_global)
      real(rkind) :: Vtotal(np_global,eVAR%dimVar)
      real(rkind) :: cf_scale_factor, cf_add_offset
      integer idim
      integer IP_glob, IP
      character(len=10) :: eStr1, eStr2
      cf_scale_factor = eVAR % cf_scale_factor
      cf_add_offset = eVAR % cf_add_offset
      IF (eVAR % dimVar == 2) THEN
        IF (TRIM(eVAR % eString) == 'WIND10') THEN
          eStr1='Uwind'
          eStr2='Vwind'
        END IF
        IF (TRIM(eVAR % eString) == 'SurfCurr') THEN
          eStr1='UsurfCurr'
          eStr2='VsurfCurr'
        END IF
      ELSE
        eStr1=TRIM(eVAR % eString)
      END IF
      cf_add_offset = eVAR % cf_add_offset
      cf_scale_factor = eVAR % cf_scale_factor
      IF (eVAR % MULTIPLE_IN .or. (myrank .eq. 0)) THEN
        CALL TEST_FILE_EXIST_DIE("Missing file : ", TRIM(eFileName))
        ISTAT = NF90_OPEN(TRIM(eFileName), NF90_NOWRITE, FID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 1, ISTAT)

        ISTAT = NF90_inq_varid(FID, TRIM(eStr1), varid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 2, ISTAT)

        ISTAT = NF90_GET_VAR(FID, varid, DirectRead, start = (/ 1, itime /), count = (/ np_global, 1 /))
        CALL GENERIC_NETCDF_ERROR(CallFct, 3, ISTAT)
        Vtotal(:,1) = cf_add_offset + cf_scale_factor*DirectRead
        IF (eVAR % dimVar == 2) THEN
          ISTAT = NF90_inq_varid(FID, TRIM(eStr2), varid)
          CALL GENERIC_NETCDF_ERROR(CallFct, 4, ISTAT)

          ISTAT = NF90_GET_VAR(FID, varid, DirectRead, start = (/ 1, itime /), count = (/ np_global, 1 /))
          CALL GENERIC_NETCDF_ERROR(CallFct, 5, ISTAT)
          Vtotal(:,2) = cf_add_offset + cf_scale_factor*DirectRead
        END IF
        ISTAT = NF90_CLOSE(FID)
        CALL GENERIC_NETCDF_ERROR(CallFct, 6, ISTAT)
      END IF
      IF (eVAR % MULTIPLE_IN) THEN
        DO IP=1,npa
          IP_glob=iplg(IP)
          outvar(IP,:)=Vtotal(IP_glob,:)
        END DO
      ELSE
        DO idim=1,eVAR % dimVar
          CALL SCATTER_ONED_ARRAY(Vtotal(:,idim), outvar(:,idim))
        END DO
      END IF
      END SUBROUTINE
!****************************************************************************
!*                                                                          *
!****************************************************************************
      SUBROUTINE READ_NETCDF_CF(eVAR, idx, outvar)
      IMPLICIT NONE
      TYPE(VAR_NETCDF_CF), intent(in)    :: eVAR
      INTEGER, INTENT(in)                :: idx
      REAL(rkind), INTENT(out)           :: outvar(npa,eVAR % dimVar)
      !
      INTEGER iFile, iTime
      character(len=140) eFileName
      !
      iFile = eVAR % VAR_SERIES % List_IFILE(idx)
      iTime = eVAR % VAR_SERIES % List_ITIME(idx)
      eFileName = eVAR % FileArr % ListFile(iFile)
      IF (eVAR % IsDirect) THEN
        CALL READ_DIRECT_NETCDF_CF(eVAR, iTime, eFileName, outvar)
      ELSE
        CALL READ_INTERP_NETCDF_CF(eVAR, iTime, eFileName, outvar)
      END IF
      END SUBROUTINE
!****************************************************************************
!*                                                                          *
!****************************************************************************
      SUBROUTINE INIT_NETCDF_CF_ALLVAR(nws, eYear, eMonth, eDay, eHour, eMin, eSec)
      IMPLICIT NONE
      integer, intent(in) :: nws
      REAL(rkind), intent(in) :: eYear, eMonth, eDay, eHour, eMin, eSec
      !
      character (len = 140), parameter :: WindFile = "Wind.nc"
      character (len = 140), parameter :: BulkFile = "Bulk.nc"
      character(len=140) :: strWind = "WIND10"
      character(len=140) :: strPair = "Pair"
      LOGICAL :: MULTIPLE_IN = .FALSE.
      integer year_i, mon_i, day_i, hour_i, min_i, sec_i
      !
      year_i = INT(eYear)
      mon_i  = INT(eMonth)
      day_i  = INT(eDay)
      hour_i = INT(eHour)
      min_i  = INT(eMin)
      sec_i  = INT(eSec)
      Print *, 'eYear  =', eYear,  ' year_i =', year_i
      Print *, 'eMonth =', eMonth, ' mon_i  =', mon_i
      Print *, 'eDay   =', eDay,   ' day_i  =', day_i
      Print *, 'eHour  =', eHour,  ' hour_i =', hour_i
      Print *, 'eMin   =', eMin,   ' min_i  =', min_i
      Print *, 'eSec   =', eSec,   ' sec_i  =', sec_i
      CALL DATE_ConvertSix2mjd(year_i, mon_i, day_i, hour_i, min_i, sec_i, eStartSim)

      allocate(outvar1(npa,1), outvar2(npa,2))
      CALL INIT_NETCDF_CF(VAR_wind, MULTIPLE_IN, WindFile, strWind, nws)
      CALL INIT_NETCDF_CF(VAR_pres, MULTIPLE_IN, BulkFile, strPair, nws)

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_NETCDF_CF_ALLVAR(time, INP_windx, INP_windy, INP_pr)
      IMPLICIT NONE
      REAL(rkind), intent(out) :: INP_windx(npa), INP_windy(npa), INP_pr(npa)
      REAL(rkind), intent(in) :: time
      !
      REAL(rkind) eTimeDay
      eTimeDay = time * SEC2DAY + eStartSim
      CALL READ_VAR_NETCDF_CF_SPECTIME(VAR_wind, eTimeDay, outvar2)
      INP_windx = outvar2(:,1)
      INP_windy = outvar2(:,2)
      CALL READ_VAR_NETCDF_CF_SPECTIME(VAR_pres, eTimeDay, outvar1)
      INP_pr = 100*outvar1(:,1)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_VAR_NETCDF_CF_SPECTIME(eVAR, eTimeDay, outvar)
      IMPLICIT NONE
      TYPE(VAR_NETCDF_CF), intent(inout)    :: eVAR
      REAL(rkind), intent(in) :: eTimeDay
      REAL(rkind), intent(out) :: outvar(npa, eVAR % dimVar)
      !
      integer record1, record2
      REAL(rkind) w1, w2
      !
      CALL GET_FRC_REC(eVAR % VAR_SERIES, eTimeDay, record1, record2, w1, w2)
      IF (record1 .ne. eVAR % idx1) THEN
        eVAR % idx1 = record1
        CALL READ_NETCDF_CF(eVAR, record1, eVAR % RawRead1)
      END IF
      IF (record2 .ne. eVAR % idx2) THEN
        eVAR % idx2 = record2
        CALL READ_NETCDF_CF(eVAR, record2, eVAR % RawRead2)
      END IF
      !
      outvar = w1*eVAR % RawRead1 + w2*eVAR % RawRead2
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************      
      SUBROUTINE GENERIC_NETCDF_ERROR(CallFct, idx, iret)
      IMPLICIT NONE
      INTEGER, intent(in) :: iret, idx
      character(*), intent(in) :: CallFct
      character(len=500) :: CHRERR
      IF (iret .NE. nf90_noerr) THEN
        CHRERR = nf90_strerror(iret)
        WRITE(16,*) TRIM(CallFct), ' -', idx, '-', TRIM(CHRERR)
        Print *, TRIM(CallFct), ' -', idx, '-', TRIM(CHRERR)
        CALL parallel_abort('Netcdf type of error')
      ENDIF
      END SUBROUTINE GENERIC_NETCDF_ERROR
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CF_EXTRACT_TIME(eStrUnitTime, ConvertToDay, eTimeStart)
      IMPLICIT NONE
      character(len=100), intent(in) :: eStrUnitTime
      REAL(rkind), intent(out) :: ConvertToDay
      REAL(rkind), intent(out) :: eTimeStart
      !
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
        ConvertToDay=DBLE(1) / DBLE(24)
      ELSEIF (TRIM(Xname) .eq. 'seconds') THEN
        ConvertToDay=DBLE(1) / DBLE(86400)
      ELSE
        CALL parallel_abort('Error in reading the unit of time')
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
      lenMonth=LEN_TRIM(YnameMonth)
      IF (lenMonth .eq. 2) THEN
        eStrTime( 5: 5)=YnameMonth( 1: 1)
        eStrTime( 6: 6)=YnameMonth( 2: 2)
      ELSE
        IF (lenMonth .eq. 1) THEN
          eStrTime( 5: 5)='0'
          eStrTime( 6: 6)=YnameMonth( 1: 1)
        ELSE
          CALL parallel_abort('Error in lenMonth')
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
          call parallel_abort('DIE in trying to get the day')
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
            CALL parallel_abort('DIE in trying to get the hour')
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
            call parallel_abort('DIE in trying to get the min')
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
            CALL parallel_abort('DIE in trying to get the sec')
          END IF
        END IF
      END IF
      CALL CT2MJD(eStrTime, eTimeStart)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GET_FRC_REC(eVAR, eTimeDay, record1, record2, w1, w2)
      ! For given WTIME return records to get and weights for time
      ! interpolation F(wwm_time)=F(rec1)*w1 + F(rec2)*w2
      ! wind_time is in the same units (s) as model_time 
      IMPLICIT NONE
      TYPE(VAR_LIST_SERIES_TIME), intent(in) :: eVAR
      REAL(rkind), INTENT(IN)             :: eTimeDay
      REAL(rkind), INTENT(OUT)            :: w1, w2
      INTEGER, INTENT(OUT)                :: record1, record2
      REAL(rkind) :: eTime1, eTime2
      INTEGER  :: iTime, I
      INTEGER  :: nbTime
      nbTime = eVAR % nbTime
 
      DO iTime=2,nbTime
        eTime1=eVAR % List_Time(iTime-1)
        eTime2=eVAR % List_Time(iTime)
        IF ((eTime1 .le. eTimeDay).and.(eTimeDay .le. eTime2)) THEN
          record2=iTime
          record1=iTime-1
          w2=(eTimeDay - eTime1)/(eTime2-eTime1)
          w1=(eTime2 - eTimeDay)/(eTime2-eTime1)
          RETURN
        END IF
      END DO
      Print *, 'eTimeDay=', eTimeDay
      Print *, 'nbTime=', nbTime
      Print *, 'List_Time(   1  )=', eVAR % List_Time(1)
      Print *, 'List_Time(nbTime)=', eVAR % List_Time(nbTime)
      CALL parallel_abort('Error in the GET_FRC_REC')
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE TEST_FILE_EXIST_DIE(string1, string2)
       CHARACTER(LEN=*) :: string1
       CHARACTER(LEN=*) :: string2
       CHARACTER(LEN=512) :: ErrMsg
       LOGICAL :: LFLIVE
       !      Print *, 'string1=', TRIM(string1)
       !      Print *, 'string2=', TRIM(string2)
       INQUIRE( FILE = TRIM(string2), EXIST = LFLIVE )
       IF ( .NOT. LFLIVE ) THEN
         WRITE(ErrMsg,10) TRIM(string1), TRIM(string2)
10       FORMAT(a, ' ', a)
         CALL parallel_abort(TRIM(ErrMsg))
       END IF
       END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      END MODULE
