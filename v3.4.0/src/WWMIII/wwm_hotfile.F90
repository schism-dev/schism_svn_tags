#include "wwm_functions.h"
MODULE wwm_hotfile_mod
      TYPE Subset
         integer eRankProc
         integer nbNeedEntries
         integer NPLOC
         integer, dimension(:), pointer :: ListNeedIndexFile
         integer, dimension(:), pointer :: ListNeedIndexMemory
      END TYPE Subset
      TYPE ReconstructInfo
         LOGICAL IsEasy
         integer nbNeedProc
         type(Subset), dimension(:), pointer :: ListSubset
      END TYPE ReconstructInfo
      CONTAINS
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PRE_CREATE_LOCAL_HOTNAME(FILEHOT, FILERET,            &
     &  MULTIPLE, HOTSTYLE, eRank)
      IMPLICIT NONE
      character(len=*), intent(in) :: FILEHOT
      character(len=*), intent(out) :: FILERET
      integer, intent(in) :: MULTIPLE, HOTSTYLE, eRank
!
      character(len=1), parameter :: ePoint = '.'
      integer, parameter :: powerproc = 6
      character(len=6) :: eStrProc
      integer :: LPOS
      integer POSITION_BEFORE_POINT
      IF (MULTIPLE.eq.0) THEN
        LPOS=POSITION_BEFORE_POINT(FILEHOT)
        IF (HOTSTYLE == 1) THEN
          WRITE (FILERET,10) TRIM(FILEHOT(1:LPOS))
  10      FORMAT (a,'.dat')
        ELSE
          WRITE (FILERET,20) TRIM(FILEHOT(1:LPOS))
  20      FORMAT (a,'.nc')
        ENDIF
      ELSE
        LPOS=POSITION_BEFORE_POINT(FILEHOT)
        CALL GETSTRING(powerproc, eRank, eStrProc)
        IF (HOTSTYLE == 1) THEN
          WRITE (FILERET,40) TRIM(FILEHOT(1:LPOS)),eStrProc
  40      FORMAT (a,'_',a,'.dat')
        ELSE
          WRITE (FILERET,50) TRIM(FILEHOT(1:LPOS)),eStrProc
  50      FORMAT (a,'_',a,'.nc')
        ENDIF
      ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CREATE_LOCAL_HOTNAME(FILEHOT, FILERET, MULTIPLE, HOTSTYLE)
      USE datapool, only : myrank
      IMPLICIT NONE
      character(LEN=*), intent(in) :: FILEHOT
      character(LEN=140), intent(OUT) :: FILERET
      integer, intent(in) :: MULTIPLE, HOTSTYLE
#ifdef MPI_PARALL_GRID
      CALL PRE_CREATE_LOCAL_HOTNAME(FILEHOT, FILERET,                  &
     &  MULTIPLE, HOTSTYLE, myrank+1)
#else
      FILERET=FILEHOT
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_AC_SIMPLE(FILEHOT, NPCALL, ACread)
      USE DATAPOOL, only : MSC, MDC, rkind, HOTIN
      USE DATAPOOL, only : NP_TOTAL, NE_TOTAL, FRLOW, FRHIGH
      IMPLICIT NONE
      character(len=*), intent(in) :: FILEHOT
      integer, intent(in) :: NPCALL
      real(rkind), intent(inout) :: ACread(MSC, MDC, NPCALL)
      integer :: NPLOC
      integer istat
      integer, allocatable :: IPLGloc(:)
      INTEGER :: HMNP, HMNE
      INTEGER :: HMSC, HMDC
      REAL(rkind) :: HFRLOW, HFRHIGH
      OPEN(HOTIN%FHNDL, FILE = TRIM(FILEHOT), STATUS = 'OLD', FORM = 'UNFORMATTED')
      READ(HOTIN%FHNDL) HMNP, HMNE
      READ(HOTIN%FHNDL) HMSC, HMDC, HFRLOW, HFRHIGH
      IF ( HMNP .NE. NP_TOTAL .OR. HMNE .NE. NE_TOTAL .OR.          &
     &     HMSC .NE. MSC      .OR. HFRLOW .NE. FRLOW .OR.           &
     &    HFRHIGH .NE. FRHIGH ) THEN
        CALL WWM_ABORT('THE HOTFILE GEOMETRY DOES NOT FIT THE INPUT FILE')
      ENDIF
      READ(HOTIN%FHNDL) NPLOC
      allocate(IPLGloc(NPLOC), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_hotfile, allocate error 1')
      READ(HOTIN%FHNDL) IPLGloc
      deallocate(IPLGloc)
!todo ordering      
      READ(HOTIN%FHNDL) ACread
      CLOSE(HOTIN%FHNDL)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DETERMINE_NUMBER_PROC(FILEHOT, HOTSTYLE, nbProc)
      IMPLICIT NONE
      character(len=*), intent(in) :: FILEHOT
      integer, intent(in) :: HOTSTYLE
      integer, intent(out) :: nbProc
      INTEGER :: MULTIPLE = 1
      character(len=140) :: FILERET
      integer :: iRankTest
      logical :: test
      iRankTest=1
      DO
        CALL PRE_CREATE_LOCAL_HOTNAME(FILEHOT, FILERET, MULTIPLE, HOTSTYLE, iRankTest)
        INQUIRE(FILE=TRIM(FILERET), EXIST=test)
        IF (test.eqv..false.) THEN
          nbProc=iRankTest-1
          EXIT
        END IF
        iRankTest=iRankTest+1
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE READ_IPLG(HOTSTYLE, FILEHOT, eRank, NPLOC, IPLG)
      USE DATAPOOL, only : HOTIN, rkind, np_total, ne_total, MSC, MDC
      USE DATAPOOL, only : FRLOW, FRHIGH
#ifdef NCDF
      USE NETCDF
#endif
      IMPLICIT NONE
      character (len = *), parameter :: CallFct="READ_IPLG"
      character(len=*), intent(in) :: FILEHOT
      integer, intent(IN) :: HOTSTYLE, eRank
      integer, intent(out) :: NPLOC, IPLG(np_total)
#ifdef NCDF
      integer :: iret, ncid, mnp_dims
#endif
      INTEGER :: HMNP, HMNE
      INTEGER :: HMSC, HMDC
      INTEGER, allocatable :: IPLGin(:)
#ifdef NCDF
      INTEGER :: iplg_id
#endif
      REAL(rkind) :: HFRLOW, HFRHIGH
      integer :: MULTIPLE, istat
      character(len=140) :: FILERET
      MULTIPLE=1
      CALL PRE_CREATE_LOCAL_HOTNAME(FILEHOT, FILERET, MULTIPLE, HOTSTYLE, eRank)
      IF (HOTSTYLE.eq.1) THEN
        OPEN(HOTIN%FHNDL, FILE = TRIM(FILERET), STATUS = 'OLD', FORM = 'UNFORMATTED')
        READ(HOTIN%FHNDL) HMNP, HMNE
        READ(HOTIN%FHNDL) HMSC, HMDC, HFRLOW, HFRHIGH
        IF ( HMNP .NE. NP_TOTAL .OR. HMNE .NE. NE_TOTAL .OR.            &
     &       HMSC .NE. MSC      .OR. HFRLOW .NE. FRLOW .OR.             &
     &       HFRHIGH .NE. FRHIGH ) THEN
          CALL WWM_ABORT('THE HOTFILE GEOMETRY DOES NOT FIT THE INPUT FILE')
        ENDIF
        READ(HOTIN%FHNDL) NPLOC
        allocate(IPLGin(NPLOC), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_hotfile, allocate error 2')
        READ(HOTIN%FHNDL) IPLGin
        IPLG(1:NPLOC)=IPLGin
        deallocate(IPLGin)
        CLOSE(HOTIN%FHNDL)
      ELSE
#ifdef NCDF
        iret=nf90_open(TRIM(FILERET), nf90_nowrite, ncid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 1, iret)
        iret=nf90_inq_dimid(ncid,"mnp",mnp_dims)
        CALL GENERIC_NETCDF_ERROR(CallFct, 2, iret)
        iret=nf90_inquire_dimension(ncid, mnp_dims, len=nploc)
        CALL GENERIC_NETCDF_ERROR(CallFct, 3, iret)
        allocate(IPLGin(NPLOC), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_hotfile, allocate error 3')
        iret=nf90_inq_varid(ncid, "iplg", iplg_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 4, iret)
        iret=nf90_get_var(ncid, iplg_id, IPLGin, start=(/1/), count=(/NPLOC/))
        CALL GENERIC_NETCDF_ERROR(CallFct, 5, iret)
        IPLG(1:NPLOC)=IPLGin
        deallocate(IPLGin)
#endif
      ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DEALLOCATE_ReconsArr(eRecons)
      IMPLICIT NONE
      type(ReconstructInfo), intent(inout) :: eRecons
      integer :: nbNeedProc, iProc
      IF (eRecons % IsEasy.eqv..false.) THEN
        nbNeedProc=eRecons % nbNeedProc
        DO iProc=1,nbNeedProc
          DEALLOCATE(eRecons % ListSubset(iProc) % ListNeedIndexFile)
          DEALLOCATE(eRecons % ListSubset(iProc) % ListNeedIndexMemory)
        END DO
        DEALLOCATE(eRecons % ListSubset)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DETERMINE_NEEDED_HOTFILES(HOTSTYLE, FILEHOT, eRecons)
#ifdef MPI_PARALL_GRID
      USE DATAPOOL, only : iplg, myrank, nproc 
#endif
      USE DATAPOOL, only : MNP, np_total
      IMPLICIT NONE
      INTEGER, intent(in) :: HOTSTYLE
      character(len=*), intent(in) :: FILEHOT
      type(ReconstructInfo), intent(inout) :: eRecons
      character(len=140) :: errmsg
      integer, allocatable :: IPLGtot(:)
      integer, allocatable :: eStatus(:)
      integer, allocatable :: ListAttained(:)
      integer :: eDiff, I, iProc, idx, eIdx, nbProc, IP
      integer :: nbNeedProc, nbF, NPLOC, nbZero, idxB
      integer istat
#ifndef MPI_PARALL_GRID
      integer, allocatable :: iplg(:)
      integer :: nproc, myrank
      nproc=1
      myrank=0
      allocate(iplg(MNP), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_hotfile, allocate error 4')
      DO IP=1,MNP
        iplg(IP)=IP
      END DO
#endif
      allocate(IPLGtot(np_total), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_hotfile, allocate error 5')
      CALL DETERMINE_NUMBER_PROC(FILEHOT, HOTSTYLE, nbProc)
      IF (nbProc.eq.nproc) THEN
        CALL READ_IPLG(HOTSTYLE, FILEHOT, myrank+1, NPLOC, IPLGtot)
        IF (NPLOC.eq.MNP) THEN
          eDiff=0
          DO IP=1,MNP
            eDiff=eDiff + abs(IPLG(IP) - IPLGtot(IP))
          END DO
          IF (eDiff.eq.0) THEN
            eRecons % IsEasy = .TRUE.
            deallocate(IPLGtot)
            RETURN
          END IF
        END IF
      END IF
      eRecons % IsEasy = .FALSE.
      allocate(eStatus(np_total), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_hotfile, allocate error 6')
      eStatus=0
      DO IP=1,MNP
        eStatus(IPLG(IP))=IP
      END DO
      nbNeedProc=0
      IF (myrank.eq.0) THEN
        allocate(ListAttained(np_total), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_hotfile, allocate error 7')
        ListAttained=0
      ENDIF
      DO iProc=1,nbProc
        CALL READ_IPLG(HOTSTYLE, FILEHOT, iProc, NPLOC, IPLGtot)
        nbF=0
        DO I=1,NPLOC
          IF (eStatus(IPLGtot(I)).gt.0) THEN
            nbF=nbF+1
          END IF
          IF (myrank.eq.0) THEN
            ListAttained(IPLGtot(I))=1
          ENDIF
        END DO
        IF (nbF.gt.0) THEN
          nbNeedProc=nbNeedProc+1
        END IF
      END DO
      IF (myrank.eq.0) THEN
        nbZero=0
        DO IP=1,np_total
          IF (ListAttained(IP).eq.0) THEN
            nbZero=nbZero+1
          ENDIF
        END DO
        IF (nbZero.gt.0) THEN
          WRITE(errmsg, *) 'Not enough data nbProc=', nbProc, ' nbMissedMode=', nbZero
          CALL WWM_ABORT(errmsg)
        ENDIF
        DEALLOCATE(ListAttained)
      ENDIF
      eRecons % nbNeedProc=nbNeedProc
      allocate(eRecons % ListSubset(nbNeedProc), stat=istat)
      IF (istat/=0) CALL WWM_ABORT('wwm_hotfile, allocate error 8')
      idx=0
      DO iProc=1,nbProc
        CALL READ_IPLG(HOTSTYLE, FILEHOT, iProc, NPLOC, IPLGtot)
        nbF=0
        DO IP=1,NPLOC
          IF (eStatus(IPLGtot(IP)).gt.0) THEN
            nbF=nbF+1
          END IF
        END DO
        IF (nbF.gt.0) THEN
          idx=idx+1
          eRecons % ListSubset(idx) % eRankProc=iProc
          eRecons % ListSubset(idx) % nbNeedEntries=nbF
          eRecons % ListSubset(idx) % NPLOC=NPLOC
          ALLOCATE(eRecons % ListSubset(idx) % ListNeedIndexFile(nbF), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_hotfile, allocate error 9')
          ALLOCATE(eRecons % ListSubset(idx) % ListNeedIndexMemory(nbF), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_hotfile, allocate error 10')
          idxB=0
          DO I=1,NPLOC
            eIdx=IPLGtot(I)
            IP=eStatus(eIdx)
            IF (IP.gt.0) THEN
              idxB=idxB+1
              eRecons % ListSubset(idx) % ListNeedIndexFile(idxB)=I
              eRecons % ListSubset(idx) % ListNeedIndexMemory(idxB)=IP
            END IF
          END DO
        END IF
      END DO
      deallocate(eStatus, IPLGtot)
#ifndef MPI_PARALL_GRID
      deallocate(iplg)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INPUT_HOTFILE
        USE DATAPOOL
        IMPLICIT NONE
        IF (HOTSTYLE_IN == 1) THEN
          CALL INPUT_HOTFILE_BINARY
#ifdef NCDF
        ELSE IF (HOTSTYLE_IN == 2) THEN
          CALL INPUT_HOTFILE_NETCDF
#endif
        ELSE
          CALL WWM_ABORT('Wrong choice of HOTSTYLE_IN')
        END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE OUTPUT_HOTFILE
        USE DATAPOOL
        IMPLICIT NONE

        IF (HOTSTYLE_OUT == 1) THEN
          CALL OUTPUT_HOTFILE_BINARY
#ifdef NCDF
        ELSE IF (HOTSTYLE_OUT == 2) THEN
          CALL OUTPUT_HOTFILE_NETCDF
#endif
        ELSE
          CALL WWM_ABORT('Wrong choice of HOTSTYLE_OUT')
        END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INPUT_HOTFILE_BINARY
#ifdef MPI_PARALL_GRID
        USE DATAPOOL, ONLY: iplg, np_global, rkind, MULTIPLEIN_HOT, MDC, MSC, MNP, HOTIN
#else
        USE DATAPOOL, ONLY: rkind, MULTIPLEIN_HOT, MDC, MSC, MNP, HOTIN
#endif
        USE DATAPOOL, ONLY: AC2, HOTSTYLE_IN
        IMPLICIT NONE
        INTEGER idxFil, idxMem, iProc, istat
        REAL(rkind), ALLOCATABLE :: ACinB(:,:,:)
        character(len=140) :: FILERET
        type(ReconstructInfo) :: eRecons
        integer :: nbF, NPLOC, eRank, I
#ifdef MPI_PARALL_GRID
        integer IP
#endif
        IF (MULTIPLEIN_HOT.eq.0) THEN
#ifdef MPI_PARALL_GRID
          ALLOCATE(ACinB(MSC,MDC,NP_GLOBAL), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_hotfile, allocate error 11')
!todo ordering of ACinB to (IS,ID,IP) ?
          CALL READ_AC_SIMPLE(HOTIN%FNAME, NP_GLOBAL, ACinB)
          DO IP=1,MNP
            AC2(:,:,IP)=ACinB(:,:,iplg(IP))
          ENDDO
          DEALLOCATE(ACinB)
#else
          CALL READ_AC_SIMPLE(HOTIN%FNAME, MNP, AC2)
#endif
        ELSE
          CALL DETERMINE_NEEDED_HOTFILES(HOTSTYLE_IN, TRIM(HOTIN%FNAME), eRecons)
          IF (eRecons % IsEasy) THEN
            CALL CREATE_LOCAL_HOTNAME(HOTIN%FNAME, FILERET, MULTIPLEIN_HOT, HOTSTYLE_IN)
            CALL READ_AC_SIMPLE(FILERET, MNP, AC2)
          ELSE
            DO iProc=1,eRecons % nbNeedProc
              eRank=eRecons % ListSubset(iProc) % eRankProc
              nbF=eRecons % ListSubset(iProc) % nbNeedEntries
              NPLOC=eRecons % ListSubset(iProc) % NPLOC
              CALL PRE_CREATE_LOCAL_HOTNAME(HOTIN%FNAME, FILERET, MULTIPLEIN_HOT, HOTSTYLE_IN, eRank)
              allocate(ACinB(MSC, MDC,NPLOC), stat=istat)
              IF (istat/=0) CALL WWM_ABORT('wwm_hotfile, allocate error 12')
              CALL READ_AC_SIMPLE(FILERET, NPLOC, ACinB)
              DO I=1,nbF
                idxFil=eRecons % ListSubset(iProc) % ListNeedIndexFile(I)
                idxMem=eRecons % ListSubset(iProc) % ListNeedIndexMemory(I)
                AC2(:,:,idxMem)=ACinB(:,:,idxFil)
              END DO
              deallocate(ACinB)
            END DO
          END IF
          CALL DEALLOCATE_ReconsArr(eRecons)
        ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE OUTPUT_HOTFILE_BINARY
#ifdef MPI_PARALL_GRID
      USE DATAPOOL, ONLY: iplg, np_global, myrank, comm, ierr, rtype
      USE DATAPOOL, ONLY: nwild_gb
#endif
      USE DATAPOOL, ONLY: NP_TOTAL, NE_TOTAL, FRLOW, FRHIGH
        USE DATAPOOL, ONLY: MSC, MDC, MNP, AC2, rkind, HOTOUT, HOTSTYLE_OUT, MULTIPLEOUT_HOT
        IMPLICIT NONE
#ifdef MPI_PARALL_GRID
        include 'mpif.h'
#endif
        CHARACTER(len=140) :: FILERET
#ifdef MPI_PARALL_GRID
        REAL(rkind), allocatable :: ACreturn(:,:,:)
        integer IP, IS, ID, istat
        REAL(rkind), allocatable :: VALB(:), VALB_SUM(:)
#endif
        CALL CREATE_LOCAL_HOTNAME(HOTOUT%FNAME, FILERET, MULTIPLEOUT_HOT, HOTSTYLE_OUT)
        OPEN(HOTOUT%FHNDL, FILE = TRIM(FILERET), STATUS = 'UNKNOWN',  FORM = 'UNFORMATTED')
        WRITE(HOTOUT%FHNDL) NP_TOTAL, NE_TOTAL
        WRITE(HOTOUT%FHNDL) MSC, MDC, FRLOW, FRHIGH
#ifndef MPI_PARALL_GRID
        WRITE(HOTOUT%FHNDL) AC2
#else
        IF (MULTIPLEOUT_HOT.eq.0) THEN
          IF (myrank.eq.0) THEN
            allocate(ACreturn(MSC, MDC,np_global), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_hotfile, allocate error 13')
          END IF
          allocate(VALB(np_global), VALB_SUM(np_global), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_hotfile, allocate error 13')
          DO ID=1,MDC
            DO IS=1,MSC
              VALB=0
              DO IP=1,MNP
                VALB(iplg(IP))=AC2(IS,ID,IP)
              END DO
              call mpi_reduce(VALB,VALB_SUM,NP_GLOBAL,rtype, MPI_SUM,0,comm,ierr)
              IF (myrank.eq.0) THEN
                ACreturn(IS,ID,:)=VALB_SUM
              END IF
            END DO
          END DO
          deallocate(VALB, VALB_SUM)
          IF (myrank.eq.0) THEN
            DO IP=1,MNP
              ACreturn(:,:,IP)=ACreturn(:,:,IP)*nwild_gb(IP)
            END DO
            WRITE(HOTOUT%FHNDL) ACreturn
            deallocate(ACreturn)
          END IF
        ELSE
          WRITE(HOTOUT%FHNDL) MNP
          WRITE(HOTOUT%FHNDL) IPLG
          WRITE(HOTOUT%FHNDL) AC2
        ENDIF
#endif
        CLOSE(HOTOUT%FHNDL)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef NCDF
      SUBROUTINE INPUT_HOTFILE_NETCDF
# ifdef MPI_PARALL_GRID
        USE DATAPOOL, ONLY : iplg, np_global, rkind, MULTIPLEIN_HOT, MDC, AC2
# else
        USE DATAPOOL, ONLY : rkind, MULTIPLEIN_HOT, MDC, AC2
# endif
        USE DATAPOOL, ONLY: HOTIN, MSC, MNP, IHOTPOS_IN, HOTSTYLE_IN
        USE NETCDF
        IMPLICIT NONE
# ifdef MPI_PARALL_GRID
        INTEGER :: IP, ID
        REAL(rkind) :: ACLOC(MSC,MDC)
# endif
        INTEGER :: NPLOC, eRank, I
        character (len = *), parameter :: CallFct="INPUT_HOTFILE_NETCDF"
        REAL(rkind), ALLOCATABLE :: ACinB(:,:,:)
        type(ReconstructInfo) :: eRecons
        character(len=140) :: FILERET
        INTEGER :: iret, ncid, ac_id, iProc, idxFil, idxMem
        INTEGER :: nbF
        integer istat
        IF (MULTIPLEIN_HOT.eq.0) THEN
# ifdef MPI_PARALL_GRID
          iret=nf90_open(HOTIN%FNAME, nf90_nowrite, ncid)
          CALL GENERIC_NETCDF_ERROR(CallFct, 1, iret)
          iret=nf90_inq_varid(ncid, "ac", ac_id)
          CALL GENERIC_NETCDF_ERROR(CallFct, 2, iret)
          DO IP=1,MNP
            iret=nf90_get_var(ncid,ac_id,ACLOC, start=(/1,1,iplg(IP),IHOTPOS_IN/), count = (/MSC, MDC, 1, 1 /))
            IF (iret /= 0) THEN
              Print *, 'This time send direcly your bug to Mathieu.Dutour@gmail.com'
              CALL GENERIC_NETCDF_ERROR(CallFct, 3, iret)
            END IF
            AC2(:,:,IP)=ACLOC
          END DO
          iret=nf90_close(ncid)
          CALL GENERIC_NETCDF_ERROR(CallFct, 4, iret)
# else
          iret=nf90_open(HOTIN%FNAME, nf90_nowrite, ncid)
          CALL GENERIC_NETCDF_ERROR(CallFct, 5, iret)
          iret=nf90_inq_varid(ncid, "ac", ac_id)
          CALL GENERIC_NETCDF_ERROR(CallFct, 6, iret)
          iret=nf90_get_var(ncid,ac_id,AC2, start=(/1,1,1,IHOTPOS_IN/), count=(/MSC,MDC,MNP,1/))
          CALL GENERIC_NETCDF_ERROR(CallFct, 7, iret)
          iret=nf90_close(ncid)
          CALL GENERIC_NETCDF_ERROR(CallFct, 8, iret)
# endif
        ELSE
          CALL DETERMINE_NEEDED_HOTFILES(HOTSTYLE_IN, TRIM(HOTIN%FNAME), eRecons)
          IF (eRecons % IsEasy) THEN
            CALL CREATE_LOCAL_HOTNAME(HOTIN%FNAME, FILERET, MULTIPLEIN_HOT, HOTSTYLE_IN)
            iret=nf90_open(FILERET, nf90_nowrite, ncid)
            CALL GENERIC_NETCDF_ERROR(CallFct, 9, iret)
            iret=nf90_inq_varid(ncid, "ac", ac_id)
            CALL GENERIC_NETCDF_ERROR(CallFct, 10, iret)
            iret=nf90_get_var(ncid,ac_id,AC2, start=(/1,1,1,IHOTPOS_IN/),  count = (/MSC, MDC, MNP, 1 /))
            CALL GENERIC_NETCDF_ERROR(CallFct, 11, iret)
            iret=nf90_close(ncid)
            CALL GENERIC_NETCDF_ERROR(CallFct, 12, iret)
          ELSE
            DO iProc=1,eRecons % nbNeedProc
              eRank=eRecons % ListSubset(iProc) % eRankProc
              nbF=eRecons % ListSubset(iProc) % nbNeedEntries
              NPLOC=eRecons % ListSubset(iProc) % NPLOC
              CALL PRE_CREATE_LOCAL_HOTNAME(HOTIN%FNAME, FILERET, MULTIPLEIN_HOT, HOTSTYLE_IN, eRank)
              iret=nf90_open(TRIM(FILERET), nf90_nowrite, ncid)
              CALL GENERIC_NETCDF_ERROR(CallFct, 13, iret)
              allocate(ACinB(MSC,MDC,NPLOC), stat=istat)
              IF (istat/=0) CALL WWM_ABORT('wwm_hotfile, allocate error 15')
              iret=nf90_inq_varid(ncid, "ac", ac_id)
              CALL GENERIC_NETCDF_ERROR(CallFct, 14, iret)
              iret=nf90_get_var(ncid,ac_id,ACinB, start=(/1,1,1,IHOTPOS_IN/), count=(/MSC, MDC, NPLOC, 1 /))
              CALL GENERIC_NETCDF_ERROR(CallFct, 15, iret)
              iret=nf90_close(ncid)
              CALL GENERIC_NETCDF_ERROR(CallFct, 16, iret)
              DO I=1,nbF
                idxFil=eRecons % ListSubset(iProc) % ListNeedIndexFile(I)
                idxMem=eRecons % ListSubset(iProc) % ListNeedIndexMemory(I)
                AC2(:,:,idxMem)=ACinB(:,:,idxFil)
              END DO
              deallocate(ACinB)
            END DO
          ENDIF
          CALL DEALLOCATE_ReconsArr(eRecons)
        ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE OUTPUT_HOTFILE_NETCDF
# ifdef MPI_PARALL_GRID
      USE DATAPOOL, ONLY: ierr, comm, rtype, myrank, iplg, np_global, ne_global, nwild_gb
# endif
      use datapool, only: MULTIPLEOUT_HOT, MULTIPLEIN_HOT, HOTSTYLE_OUT, MDC, MSC, MNP, MNE, RKIND
      use datapool, only: HOTOUT, IDXHOTOUT, LCYCLEHOT, WriteOutputProcess_hot, NF90_RUNTYPE
      use datapool, only: AC2, MAIN
      USE NETCDF
      IMPLICIT NONE
# ifdef MPI_PARALL_GRID
      include 'mpif.h'
# endif
      character (len = *), parameter :: CallFct="OUTPUT_HOTFILE_NETCDF"
      INTEGER :: POS
      integer :: iret, ncid, ntime_dims, mnp_dims, msc_dims, mdc_dims
      integer :: ac_id
      integer :: nbTime
      REAL(rkind)  :: eTimeDay
      character (len = *), parameter :: UNITS = "units"
      character(len=140) :: FILERET
      integer np_write, ne_write
# ifdef MPI_PARALL_GRID
      integer ID, IS, IP, ISTAT
      REAL(rkind), ALLOCATABLE :: ACreturn(:,:,:)
      REAL(rkind), ALLOCATABLE :: VALB(:), VALB_SUM(:)
# endif
# ifdef MPI_PARALL_GRID
      IF (MULTIPLEOUT_HOT.eq.0) THEN
        np_write=np_global
        ne_write=ne_global
      ELSE
        np_write=MNP
        ne_write=MNE
      ENDIF
# else
      np_write=MNP
      ne_write=MNE
# endif
      CALL CREATE_LOCAL_HOTNAME(HOTOUT%FNAME, FILERET, MULTIPLEOUT_HOT, HOTSTYLE_OUT)
      IF (IDXHOTOUT.eq.0) THEN
!$OMP MASTER
        IF (LCYCLEHOT) THEN
          nbTime=2
        ELSE
          nbTime=-1
        END IF
        IF (WriteOutputProcess_hot) THEN
          iret = nf90_create(FILERET, NF90_CLOBBER, ncid)
          CALL GENERIC_NETCDF_ERROR(CallFct, 1, iret)
          CALL WRITE_NETCDF_HEADERS_1(ncid, nbTime, MULTIPLEOUT_HOT, np_write, ne_write)
          iret=nf90_inq_dimid(ncid, "mnp", mnp_dims)
          CALL REPORT_ERROR_INQ(iret, "mnp")

          iret=nf90_inq_dimid(ncid, "msc", msc_dims)
          CALL REPORT_ERROR_INQ(iret, "msc")

          iret=nf90_inq_dimid(ncid, "mdc", mdc_dims)
          CALL REPORT_ERROR_INQ(iret, 'mdc')

          iret=nf90_inq_dimid(ncid, 'ocean_time', ntime_dims)
          CALL REPORT_ERROR_INQ(iret, 'ocean_time')

          iret=nf90_def_var(ncid,"ac",NF90_RUNTYPE,(/ msc_dims, mdc_dims, mnp_dims, ntime_dims/),ac_id)
          CALL GENERIC_NETCDF_ERROR(CallFct, 2, iret)

          iret=nf90_put_att(ncid,ac_id,UNITS,'unknown')
          CALL GENERIC_NETCDF_ERROR(CallFct, 3, iret)

          iret=nf90_close(ncid)
          CALL GENERIC_NETCDF_ERROR(CallFct, 4, iret)
          !
          iret = nf90_open(FILERET, nf90_write, ncid)
          CALL GENERIC_NETCDF_ERROR(CallFct, 5, iret)
        END IF
        CALL WRITE_NETCDF_HEADERS_2(ncid, MULTIPLEOUT_HOT, WriteOutputProcess_hot, np_write, ne_write)
        IF (WriteOutputProcess_hot) THEN
          iret = nf90_close(ncid)
          CALL GENERIC_NETCDF_ERROR(CallFct, 6, iret)
        END IF
!$OMP END MASTER
      END IF
# ifdef MPI_PARALL_GRID
      IF (MULTIPLEOUT_HOT.eq.0) THEN
        IF (myrank.eq.0) THEN
          allocate(ACreturn(MSC,MDC,np_global), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_hotfile, allocate error 16')
        END IF
        allocate(VALB(np_global), VALB_SUM(np_global), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_hotfile, allocate error 17')
        DO ID=1,MDC
          DO IS=1,MSC
            VALB=0
            DO IP=1,MNP
              VALB(iplg(IP))=AC2(IS,ID,IP)
            END DO
            call mpi_reduce(VALB,VALB_SUM,NP_GLOBAL,rtype, MPI_SUM,0,comm,ierr)
            IF (myrank.eq.0) THEN
              ACreturn(IS,ID,:)=VALB_SUM
            END IF
          END DO
        END DO
        IF (myrank.eq.0) THEN
          DO IP=1,np_global
            ACreturn(:,:,IP)=ACreturn(:,:,IP)*nwild_gb(IP)
          END DO
        END IF
        deallocate(VALB, VALB_SUM)
      ENDIF
# endif
      IF (WriteOutputProcess_hot) THEN
        IF (LCYCLEHOT) THEN
          POS=mod(IDXHOTOUT,2)+1
        ELSE
          POS=IDXHOTOUT+1
        END IF
        eTimeDay=MAIN%TMJD
        !
        iret=nf90_open(FILERET, nf90_write, ncid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 7, iret)
        CALL WRITE_NETCDF_TIME(ncid, POS, eTimeDay)
        iret=nf90_inq_varid(ncid, "ac", ac_id)
        CALL GENERIC_NETCDF_ERROR(CallFct, 8, iret)
# ifdef MPI_PARALL_GRID
        IF (MULTIPLEOUT_HOT.eq.0) THEN
          iret=nf90_put_var(ncid,ac_id,ACreturn,start=(/1, 1, 1, POS/), count=(/ MSC, MDC, np_global, 1 /))
          CALL GENERIC_NETCDF_ERROR(CallFct, 9, iret)
          deallocate(ACreturn);
        ELSE
          iret=nf90_put_var(ncid,ac_id,AC2,start=(/1, 1, 1, POS/), count=(/ MSC, MDC, MNP, 1 /))
          CALL GENERIC_NETCDF_ERROR(CallFct, 10, iret)
        ENDIF
# else
        iret=nf90_put_var(ncid,ac_id,AC2,start=(/1, 1, 1, POS/), count=(/ MSC, MDC, MNP, 1 /))
        CALL GENERIC_NETCDF_ERROR(CallFct, 11, iret)
# endif
        iret=nf90_close(ncid)
        CALL GENERIC_NETCDF_ERROR(CallFct, 12, iret)
      ENDIF
      IDXHOTOUT=IDXHOTOUT+1
      END SUBROUTINE
#endif
END MODULE wwm_hotfile_mod
