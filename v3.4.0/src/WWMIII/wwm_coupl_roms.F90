#include "wwm_functions.h"
#define DEBUG_WWM
#undef DEBUG_WWM
#if !defined PGMCL_COUPLING
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_PIPES_ROMS()
      USE DATAPOOL
      IMPLICIT NONE
!
! open pipe data files for coupling
!
      LSEWL = .TRUE.
      LSECU = .TRUE.
      WRITE(DBG%FHNDL,'("+TRACE...",A)') 'OPEN PIPE ROMS'
      FLUSH(DBG%FHNDL)
!     Pipes that are read by the wave model
      OPEN(1000,file='pipe/ExchRW'  ,form='unformatted', action='read')
      WRITE(DBG%FHNDL,*) 'WWM: open pipe ExchImport'
      FLUSH(DBG%FHNDL)
!     Pipes that are written by the wave modell
      OPEN(101 ,file='pipe/ExchWR' ,form='unformatted', action='write')
      WRITE(DBG%FHNDL,*) 'WWM: open pipe ExchExport'
      FLUSH(DBG%FHNDL)
      WRITE(DBG%FHNDL,'("+TRACE...",A)') 'END OPEN PIPE ROMS'
      FLUSH(DBG%FHNDL)

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE TERMINATE_PIPES_ROMS()
      USE DATAPOOL
      IMPLICIT NONE
      close(1000)
      close(101)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PIPE_ROMS_IN(K,IFILE,IT)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: K,IFILE,IT
      INTEGER              :: IP
# ifdef WWM_MPI
      REAL(rkind), allocatable :: WINDXY_TOT(:,:), CURTXY_TOT(:,:), WATLEV_TOT(:)
      real(rkind), allocatable :: rbuf_real(:)
      integer idx, iProc
# endif
        LCALC=.TRUE.
      IF ( K-INT(K/MAIN%ICPLT)*MAIN%ICPLT .EQ. 0 ) THEN
        WATLEVOLD=WATLEV
        DELTAT_WATLEV = MAIN%DTCOUP
        LCALC=.TRUE.
        WRITE(DBG%FHNDL,'("+TRACE...",A)') 'READING PIPE'
        FLUSH(DBG%FHNDL)
# ifndef WWM_MPI
        DO IP = 1, MNP
          READ(1000) WINDXY(IP,1), WINDXY(IP,2), CURTXY(IP,1), CURTXY(IP,2), WATLEV(IP)
        END DO
# else
        allocate(WINDXY_TOT(np_global,2), CURTXY_TOT(np_global,2), WATLEV_TOT(np_global), rbuf_real(np_global*5), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_coupl_roms, allocate err')
        IF (myrank.eq.0) THEN
          DO IP = 1, np_global
            READ(1000) WINDXY_TOT(IP,1), WINDXY_TOT(IP,2), CURTXY_TOT(IP,1), CURTXY_TOT(IP,2), WATLEV_TOT(IP)
          END DO
          DO IP=1,np_global
            idx=idx+1
            rbuf_real(idx)=WINDXY_TOT(IP,1)
            idx=idx+1
            rbuf_real(idx)=WINDXY_TOT(IP,2)
            idx=idx+1
            rbuf_real(idx)=CURTXY_TOT(IP,1)
            idx=idx+1
            rbuf_real(idx)=CURTXY_TOT(IP,2)
            idx=idx+1
            rbuf_real(idx)=WATLEV_TOT(IP)
          END DO
          DO iProc=2,nproc
            CALL MPI_SEND(rbuf_real,np_global*5,MPI_REAL8, iProc-1, 196, comm, ierr)
          END DO
        ELSE
          CALL MPI_RECV(rbuf_real,np_global*5,MPI_REAL8, 0, 196, comm, istatus, ierr)
          idx=0
          DO IP=1,np_global
            idx=idx+1
            WINDXY_TOT(IP,1)=rbuf_real(idx)
            idx=idx+1
            WINDXY_TOT(IP,2)=rbuf_real(idx)
            idx=idx+1
            CURTXY_TOT(IP,1)=rbuf_real(idx)
            idx=idx+1
            CURTXY_TOT(IP,2)=rbuf_real(idx)
            idx=idx+1
            WATLEV_TOT(IP)=rbuf_real(idx)
          END DO
        END IF
        DO IP = 1, MNP
          WINDXY(IP,:)=WINDXY_TOT(iplg(IP),:)
          CURTXY(IP,:)=CURTXY_TOT(iplg(IP),:)
          WATLEV(IP)=WATLEV_TOT(iplg(IP))
        END DO
        deallocate(rbuf_real, WINDXY_TOT, CURTXY_TOT, WATLEV_TOT)
# endif
        WRITE(DBG%FHNDL,'("+TRACE...",A)') 'END READING PIPE'
        FLUSH(DBG%FHNDL)
      END IF
      IF (K == 1) CALL INITIAL_CONDITION(IFILE,IT)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PIPE_ROMS_OUT(K)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(IN)  :: K
      INTEGER              :: IP
      REAL(rkind)          :: ACLOC(MSC,MDC)
      REAL(rkind)          :: HS,WLM,LPP,FPP,CPP,BOTEXPER
      REAL(rkind)          :: UBOT,TM01,TM10
      REAL(rkind)          :: TMBOT, KPP,DM,DSPR,ORBITAL,ETOTS,ETOTC,WNPP,TPP,CGPP
      REAL(rkind)          :: PEAKDSPR, PEAKDM, HSWE, HSLIM, TM02, KLM, DPEAK
      REAL(rkind)          :: TPPD,KPPD,CGPD,CPPD
# ifdef WWM_MPI
      REAL(rkind), allocatable :: OUTT(:,:), OUTT_TOT(:,:)
      REAL(rkind)    :: TP
# endif
      IF ( K-INT(K/MAIN%ICPLT)*MAIN%ICPLT .EQ. 0 ) THEN
# ifndef WWM_MPI
        DO IP = 1, MNP
          ACLOC = AC2(:,:,IP)
          CALL MEAN_PARAMETER(IP,ACLOC,MSC,HS,TM01,TM02,TM10,KLM,WLM)
          CALL WAVE_CURRENT_PARAMETER(IP,ACLOC,UBOT,ORBITAL,BOTEXPER,TMBOT,'PIPE_ROMS_OUT 1')
          CALL MEAN_DIRECTION_AND_SPREAD(IP,ACLOC,MSC,ETOTS,ETOTC,DM,DSPR)
          CALL PEAK_PARAMETER(IP,ACLOC,MSC,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD)
! HS, HSWE, HSLIM  ! - Significant wave height (m) -- HS
! DM  ! - Wave direction (degrees)
! TPP ! - Surface wave relative peak period (s) -- TP
! WLM, KME, SME ! - Average Wave Length [m] - LME
! ORBITAL(IP) ! - Wave bottom orbital velocity (m/s)
! TMBOT ! - Bottom wave period (s)
! DISSIPATION(IP) ! - Wave energy dissipation (W/m2)
! QBLOCAL(IP) ! - Percent of breakig waves (nondimensional)
! DSPR ! - directional spreading
! PEAKDSPR ! - peak directional spreading
! PEAKDM ! - Peak direction
!
!AR: what for you need this HSWE and HSLIM and so on ... what is exactly SME in your definition ...
!AR: I have deleted them ...
!
          HSWE=0
          HSLIM=0
          WRITE(101)  HS, HSWE,                                        &
     &                   HSLIM, DM,                                    &
     &                   TPP, WLM,                                     &
     &                   KLM, TM01,                                    &
     &                   ORBITAL, TMBOT,                               &
     &                   DISSIPATION(IP), QBLOCAL(IP),                 &
     &                   DSPR, PEAKDSPR,                               &
     &                   PEAKDM, TM02
          FLUSH(101)
        END DO
# else
        allocate(OUTT(np_global,16), OUTT_TOT(np_global,16), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_coupl_roms, allocate err')
        OUTT=0
        DO IP = 1, MNP
          ACLOC = AC2(:,:,IP)
          CALL MEAN_PARAMETER(IP,ACLOC,MSC,HS,TM01,TM02,TM10,KLM,WLM)
          CALL WAVE_CURRENT_PARAMETER(IP,ACLOC,UBOT,ORBITAL,BOTEXPER,TMBOT,'PIPE_ROMS_OUT 2')
          CALL MEAN_DIRECTION_AND_SPREAD(IP,ACLOC,MSC,ETOTS,ETOTC,DM,DSPR)
          CALL PEAK_PARAMETER(IP,ACLOC,MSC,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD)
          HSWE=0
          HSLIM=0
          OUTT(iplg(IP), 1)=HS
          OUTT(iplg(IP), 2)=HSWE
          OUTT(iplg(IP), 3)=HSLIM
          OUTT(iplg(IP), 4)=DM
          OUTT(iplg(IP), 5)=TPP
          OUTT(iplg(IP), 6)=WLM
          OUTT(iplg(IP), 7)=KLM
          OUTT(iplg(IP), 8)=TM01
          OUTT(iplg(IP), 9)=ORBITAL
          OUTT(iplg(IP),10)=TMBOT
          OUTT(iplg(IP),11)=DISSIPATION(IP)
          OUTT(iplg(IP),12)=QBLOCAL(IP)
          OUTT(iplg(IP),13)=DSPR
          OUTT(iplg(IP),14)=PEAKDSPR
          OUTT(iplg(IP),15)=PEAKDM
          OUTT(iplg(IP),16)=TM02
        END DO
        call mpi_reduce(OUTT,OUTT_TOT,NP_GLOBAL*16,rtype,MPI_SUM,0,comm,ierr)
        IF (myrank.eq.0) THEN
          DO IP=1,NP_GLOBAL
            OUTT_TOT(IP,:)=OUTT_TOT(IP,:)/nwild_gb(IP)
            WRITE(101) OUTT_TOT(IP, 1), OUTT_TOT(IP, 2),                 &
     &                    OUTT_TOT(IP, 3), OUTT_TOT(IP, 4),              &
     &                    OUTT_TOT(IP, 5), OUTT_TOT(IP, 6),              &
     &                    OUTT_TOT(IP, 7), OUTT_TOT(IP, 8),              &
     &                    OUTT_TOT(IP, 9), OUTT_TOT(IP,10),              &
     &                    OUTT_TOT(IP,11), OUTT_TOT(IP,12),              &
     &                    OUTT_TOT(IP,13), OUTT_TOT(IP,14),              &
     &                    OUTT_TOT(IP,15), OUTT_TOT(IP,16)
          END DO
        END IF
        deallocate(OUTT, OUTT_TOT)
# endif
      END IF
      WRITE(DBG%FHNDL,*) 'export WWM: ending of writing data'
      FLUSH(DBG%FHNDL)
      END SUBROUTINE
#endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef PGMCL_COUPLING
# ifdef WWM_MPI
      SUBROUTINE WWM_CreateMatrixPartition
        USE DATAPOOL
        USE mod_coupler
        implicit none
        integer, allocatable :: rbuf_int(:)
        integer, allocatable :: TheIndex(:)
        integer, allocatable :: NumberNode(:), NumberTrig(:)
        integer, allocatable :: All_LocalToGlobal(:,:)
        integer i, eIdx, iProc, MNPloc, MNEloc, idx
        integer IPc, IP
#  ifdef DEBUG_WWM
        integer MinValIndex, MinValIndexInv, eVal
#  endif
        allocate(MatrixBelongingWAV(np_global, NnodesWAV), NumberNode(NnodesWAV), NumberTrig(NnodesWAV), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_coupl_roms, allocate error 3')
        IF (myrank.ne.MyRankLocal) THEN
          CALL WWM_ABORT('die from ignominious death')
        END IF
        IF (MyRankLocal.eq.0) THEN
          allocate(All_LocalToGlobal(np_global, NnodesWAV), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_coupl_roms, allo error 4')
          All_LocalToGlobal=0
          MatrixBelongingWAV=0
          DO i=1,MNP
            eIdx=iplg(i)
            MatrixBelongingWAV(eIdx,1)=i
            All_LocalToGlobal(i,1)=eIdx
          ENDDO
          NumberNode(1)=MNP
          DO iProc=2,NnodesWAV
            allocate(rbuf_int(1), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_coupl_roms, alloc error 5')
            CALL MPI_RECV(rbuf_int,1,itype, iProc-1, 194, WAV_COMM_WORLD, istatus, ierr)
            MNPloc=rbuf_int(1)
            NumberNode(iProc)=MNPloc
            deallocate(rbuf_int)
            !
            allocate(rbuf_int(MNPloc), stat=istat)
            IF (istat/=0) CALL WWM_ABORT('wwm_coupl_roms, alloc error 6')
            CALL MPI_RECV(rbuf_int,MNPloc,itype, iProc-1, 195, WAV_COMM_WORLD, istatus, ierr)
            DO IP=1,MNPloc
              eIdx=rbuf_int(IP)
              MatrixBelongingWAV(eIdx,iProc)=IP
              All_LocalToGlobal(IP,iProc)=eIdx
            END DO
            deallocate(rbuf_int)
          END DO
          !
          allocate(rbuf_int(np_global*NnodesWAV), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_coupl_roms, alloc error 7')
          idx=0
          DO iProc=1,NnodesWAV
            DO IP=1,np_global
              idx=idx+1
              rbuf_int(idx)=MatrixBelongingWAV(IP,iProc)
            END DO
          END DO
          DO iProc=2,NnodesWAV
            CALL MPI_SEND(rbuf_int,np_global*NnodesWAV,itype, iProc-1, 196, WAV_COMM_WORLD, ierr)
          END DO
          deallocate(rbuf_int)
        ELSE
          allocate(rbuf_int(1), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_coupl_roms, alloc error 8')
          rbuf_int(1)=MNP
          CALL MPI_SEND(rbuf_int,1,itype, 0, 194, WAV_COMM_WORLD, ierr)
          deallocate(rbuf_int)

          CALL MPI_SEND(iplg,MNP,itype, 0, 195, WAV_COMM_WORLD, ierr)
!
          allocate(rbuf_int(np_global*NnodesWAV), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_coupl_roms, alloc error 10')
          CALL MPI_RECV(rbuf_int,np_global*NnodesWAV,itype, 0, 196, WAV_COMM_WORLD, istatus, ierr)
          idx=0
          DO iProc=1,NnodesWAV
            DO i=1,np_global
              idx=idx+1
              MatrixBelongingWAV(i,iProc)=rbuf_int(idx)
            END DO
          END DO
          deallocate(rbuf_int)
        ENDIF
        deallocate(NumberNode, NumberTrig)
      END SUBROUTINE
# else
      SUBROUTINE WWM_CreateMatrixPartition
        USE DATAPOOL
        USE mod_coupler
        IMPLICIT NONE
        integer IP, istat
        allocate(MatrixBelongingWAV(MNP, 1), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_coupl_roms, allocate error 12')
        DO IP=1,MNP
          MatrixBelongingWAV(IP,1)=IP
        ENDDO
      END SUBROUTINE
# endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WWM_DeallocMatrixPartition
      USE DATAPOOL
      USE mod_coupler
      IMPLICIT NONE
      deallocate(MatrixBelongingWAV)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WWM_CreateGlobalLON_LAT_ListTrig
      USE mod_coupler
      USE DATAPOOL
      implicit none
      CALL GET_GRID_ARRAY_FE_r8(NP_TOTAL, NE_TOTAL, XPtotal, YPtotal, INEtotal, eGrid_wav)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ROMS_COUPL_INITIALIZE
        USE DATAPOOL, only : MNP, rkind, DEP, XP, YP, np_total, ne_total
        USE DATAPOOL, only : DBG, istatus, ierr, itype, myrank
        USE mod_coupler
        USE PGMCL_LIBRARY
        USE pgmcl_interp
        implicit none
        logical DoNearest
        integer rbuf_int(1)
        integer IP, iNodeSel, idx, eRankRecv
        integer istat
        real(rkind) eDiff, AbsDiff, SumDep1, SumDep2, SumDiff
        real(rkind) minBathy, maxBathy
!        character(len=40) :: FileSave1, FileSave2
!        character(len=3) :: eStrFi
        real(rkind) SumDepReceive
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 1, rnk=', myrank
        FLUSH(DBG%FHNDL)
# endif
        CALL SetComputationalNodes(ArrLocal, NnodesWAV, OCNid)
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 1.2, rnk=', myrank
        FLUSH(DBG%FHNDL)
# endif
        CALL WWM_CreateMatrixPartition
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 1.3, rnk=', myrank
        FLUSH(DBG%FHNDL)
# endif
        CALL WWM_CreateGlobalLON_LAT_ListTrig
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 2, rnk=', myrank
        FLUSH(DBG%FHNDL)
# endif
        IF (MyRankLocal.eq.0) THEN
          CALL M2M_send_fem(ArrLocal, OCNid, eGrid_wav)
          CALL M2M_send_node_partition(ArrLocal, OCNid,                 &
     &        np_total, NnodesWAV, MatrixBelongingWAV)
        ENDIF
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 3, rnk=', myrank
        FLUSH(DBG%FHNDL)
# endif
        CALL M2M_recv_grid(ArrLocal, OCNid, eGrid_ocn_rho)
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 4, rnk=', myrank
        FLUSH(DBG%FHNDL)
# endif
        CALL M2M_recv_grid(ArrLocal, OCNid, eGrid_ocn_u)
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 5, rnk=', myrank
        FLUSH(DBG%FHNDL)
# endif
        CALL M2M_recv_grid(ArrLocal, OCNid, eGrid_ocn_v)
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 6, rnk=', myrank
        FLUSH(DBG%FHNDL)
# endif
        CALL M2M_recv_node_partition(ArrLocal, OCNid,                   &
     &   NnodeRho, NnodesOCN, MatrixBelongingOCN_rho)
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 7, rnk=', myrank
        FLUSH(DBG%FHNDL)
# endif
        CALL M2M_recv_node_partition(ArrLocal, OCNid,                   &
     &   NnodeU, NnodesOCN, MatrixBelongingOCN_u)
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 8, rnk=', myrank
        FLUSH(DBG%FHNDL)
# endif
        CALL M2M_recv_node_partition(ArrLocal, OCNid,                   &
     &   NnodeV, NnodesOCN, MatrixBelongingOCN_v)
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 9, rnk=', myrank
        FLUSH(DBG%FHNDL)
# endif
        eRankRecv=ArrLocal % ListFirstRank(OCNid)
        CALL MPI_RECV(rbuf_int,1,itype, eRankRecv, 103, MPI_COMM_WORLD, istatus, ierr)
        Nlevel=rbuf_int(1)
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 10, rnk=', myrank
        FLUSH(DBG%FHNDL)
# endif
        ALLOCATE(z_w_loc(0:Nlevel), eUSTOKES_loc(Nlevel), eVSTOKES_loc(Nlevel), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_coupl_roms, allocate error 16')
        DoNearest=.TRUE.
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 11, rnk=', myrank
        FLUSH(DBG%FHNDL)
# endif
        CALL GetString(MyRankGlobal, eStr)
        FileSave_OCNtoWAV_rho='InterpSave_OCNtoWAV_rho'
        CALL SAVE_CreateInterpolationSparseMatrix(                      &
     &    FileSave_OCNtoWAV_rho, eStr, mMat_OCNtoWAV_rho, DoNearest,    &
     &    eGrid_ocn_rho, eGrid_wav)
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 12, rnk=', myrank
        FLUSH(DBG%FHNDL)
# endif
        FileSave_OCNtoWAV_u='InterpSave_OCNtoWAV_u'
        CALL SAVE_CreateInterpolationSparseMatrix(                      &
     &    FileSave_OCNtoWAV_u, eStr, mMat_OCNtoWAV_u, DoNearest,        &
     &    eGrid_ocn_u, eGrid_wav)
        FileSave_OCNtoWAV_v='InterpSave_OCNtoWAV_v'
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 13, rnk=', myrank
        FLUSH(DBG%FHNDL)
# endif
        CALL SAVE_CreateInterpolationSparseMatrix(                      &
     &    FileSave_OCNtoWAV_v, eStr, mMat_OCNtoWAV_v, DoNearest,        &
     &    eGrid_ocn_v, eGrid_wav)
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 14', myrank
        FLUSH(DBG%FHNDL)
# endif
        CALL MPI_INTERP_GetSystemOutputSide(ArrLocal, OCNid, WAVid,     &
     &    MatrixBelongingOCN_rho, MatrixBelongingWAV,                   &
     &    mMat_OCNtoWAV_rho, TheArr_OCNtoWAV_rho)
# ifndef NO_ASYNC
        CALL MPI_INTERP_GetAsyncOutput_r8(TheArr_OCNtoWAV_rho,          &
     &    3, TheAsync_OCNtoWAV_uvz)
        CALL MPI_INTERP_GetAsyncOutput_r8(TheArr_OCNtoWAV_rho,          &
     &    Nlevel+1, TheAsync_OCNtoWAV_rho)
# endif
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'nbNeedTot=', TheArr_OCNtoWAV_rho % nbNeedTot
        WRITE(DBG%FHNDL,*) 'nbProc=', TheArr_OCNtoWAV_rho % nbProc
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, WAV, step 14'
        FLUSH(DBG%FHNDL)
# endif
        CALL MPI_INTERP_GetSystemOutputSide(ArrLocal, OCNid, WAVid,     &
     &    MatrixBelongingOCN_u, MatrixBelongingWAV,                     &
     &    mMat_OCNtoWAV_u, TheArr_OCNtoWAV_u)
# ifndef NO_ASYNC
#  ifdef FIRST_ORDER_ARDHUIN
        CALL MPI_INTERP_GetAsyncOutput_r8(TheArr_OCNtoWAV_u,            &
     &    2, TheAsync_OCNtoWAV_u)
#  else
        CALL MPI_INTERP_GetAsyncOutput_r8(TheArr_OCNtoWAV_u,            &
     &    Nlevel, TheAsync_OCNtoWAV_u)
#  endif
# endif
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 15, rnk=', myrank
        FLUSH(DBG%FHNDL)
# endif
        CALL MPI_INTERP_GetSystemOutputSide(ArrLocal, OCNid, WAVid,     &
     &    MatrixBelongingOCN_v, MatrixBelongingWAV,                     &
     &    mMat_OCNtoWAV_v, TheArr_OCNtoWAV_v)
# ifndef NO_ASYNC
#  ifdef FIRST_ORDER_ARDHUIN
        CALL MPI_INTERP_GetAsyncOutput_r8(TheArr_OCNtoWAV_v,            &
     &    2, TheAsync_OCNtoWAV_v)
#  else
        CALL MPI_INTERP_GetAsyncOutput_r8(TheArr_OCNtoWAV_v,            &
     &    Nlevel, TheAsync_OCNtoWAV_v)
#  endif
# endif
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 16, rnk=', myrank
        FLUSH(DBG%FHNDL)
# endif
        CALL DeallocSparseMatrix(mMat_OCNtoWAV_rho)
        CALL DeallocSparseMatrix(mMat_OCNtoWAV_u)
        CALL DeallocSparseMatrix(mMat_OCNtoWAV_v)
        FileSave_WAVtoOCN_rho='InterpSave_WAVtoOCN_rho'
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 17, rnk=', myrank
        FLUSH(DBG%FHNDL)
# endif
        CALL SAVE_CreateInterpolationSparseMatrix(                      &
     &    FileSave_WAVtoOCN_rho, eStr, mMat_WAVtoOCN_rho, DoNearest,    &
     &    eGrid_wav, eGrid_ocn_rho)
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 18, rnk=', myrank
        FLUSH(DBG%FHNDL)
# endif
        FileSave_WAVtoOCN_u='InterpSave_WAVtoOCN_u'
        CALL SAVE_CreateInterpolationSparseMatrix(                      &
     &    FileSave_WAVtoOCN_u, eStr, mMat_WAVtoOCN_u, DoNearest,        &
     &    eGrid_wav, eGrid_ocn_u)
        FileSave_WAVtoOCN_v='InterpSave_WAVtoOCN_v'
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 19, rnk=', myrank
        FLUSH(DBG%FHNDL)
# endif
        CALL SAVE_CreateInterpolationSparseMatrix(                      &
     &    FileSave_WAVtoOCN_v, eStr, mMat_WAVtoOCN_v, DoNearest,        &
     &    eGrid_wav, eGrid_ocn_v)
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 20, rnk=', myrank
        FLUSH(DBG%FHNDL)
# endif
        CALL MPI_INTERP_GetSystemInputSide(ArrLocal, WAVid, OCNid,      &
     &    MatrixBelongingWAV, MatrixBelongingOCN_rho,                   &
     &    mMat_WAVtoOCN_rho, TheArr_WAVtoOCN_rho)
# ifndef NO_ASYNC
        CALL MPI_INTERP_GetAsyncInput_r8(TheArr_WAVtoOCN_rho,           &
     &    19, TheAsync_WAVtoOCN_stat)
# endif
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 21, rnk=', myrank
        FLUSH(DBG%FHNDL)
# endif
        CALL MPI_INTERP_GetSystemInputSide(ArrLocal, WAVid, OCNid,      &
     &    MatrixBelongingWAV, MatrixBelongingOCN_u,                     &
     &    mMat_WAVtoOCN_u, TheArr_WAVtoOCN_u)
# ifndef NO_ASYNC
#  ifdef STOKES_DRIFT_USING_INTEGRAL
        CALL MPI_INTERP_GetAsyncInput_r8(TheArr_WAVtoOCN_u,             &
     &    Nlevel+1, TheAsync_WAVtoOCN_u)
#  else
        CALL MPI_INTERP_GetAsyncInput_r8(TheArr_WAVtoOCN_u,             &
     &    1, TheAsync_WAVtoOCN_u)
#  endif
# endif
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 22, rnk=', myrank
        FLUSH(DBG%FHNDL)
# endif
        CALL MPI_INTERP_GetSystemInputSide(ArrLocal, WAVid, OCNid,      &
     &    MatrixBelongingWAV, MatrixBelongingOCN_v,                     &
     &    mMat_WAVtoOCN_v, TheArr_WAVtoOCN_v)
# ifndef NO_ASYNC
#  ifdef STOKES_DRIFT_USING_INTEGRAL
        CALL MPI_INTERP_GetAsyncInput_r8(TheArr_WAVtoOCN_v,             &
     &    Nlevel+1, TheAsync_WAVtoOCN_v)
#  else
        CALL MPI_INTERP_GetAsyncInput_r8(TheArr_WAVtoOCN_v,             &
     &    1, TheAsync_WAVtoOCN_v)
#  endif
# endif
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 23, rnk=', myrank
        FLUSH(DBG%FHNDL)
# endif
        CALL DeallocSparseMatrix(mMat_WAVtoOCN_rho)
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 23.1'
        FLUSH(DBG%FHNDL)
# endif
        CALL DeallocSparseMatrix(mMat_WAVtoOCN_u)
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 23.2'
        FLUSH(DBG%FHNDL)
# endif
        CALL DeallocSparseMatrix(mMat_WAVtoOCN_v)
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 23.3'
        FLUSH(DBG%FHNDL)
# endif
# ifdef FIRST_ORDER_ARDHUIN
        allocate(A_wav_ur_3D(2,MNP), A_wav_vr_3D(2,MNP), U_wav(2,MNP), V_wav(2,MNP), stat=istat)
# else
        allocate(A_wav_ur_3D(Nlevel,MNP), A_wav_vr_3D(Nlevel,MNP), U_wav(Nlevel, MNP), V_wav(Nlevel,MNP), stat=istat)
# endif
        IF (istat/=0) CALL WWM_ABORT('wwm_coupl_roms, allocate error 23.4')
        allocate(CosAng(MNP), SinAng(MNP), dep_rho(MNP), A_wav_rho_3D(Nlevel+1,MNP), A_wav_stat(19,MNP), A_wav_uvz(3,MNP), A_wav_rho(MNP), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_coupl_roms, allocate error 16.1')
# ifdef STOKES_DRIFT_USING_INTEGRAL
        allocate(A_wav_u_3D(Nlevel+1,MNP), A_wav_v_3D(Nlevel+1,MNP), stat=istat)
# else
        allocate(A_wav_u_3D(1,MNP), A_wav_v_3D(1,MNP), stat=istat)
# endif
        IF (istat/=0) CALL WWM_ABORT('wwm_coupl_roms, allocate error 17')
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 24, rnk=', myrank
        FLUSH(DBG%FHNDL)
# endif
        CALL MPI_INTERP_RECV_r8(TheArr_OCNtoWAV_rho, 23, A_wav_rho)
        DO IP=1,MNP
          CosAng(IP)=COS(A_wav_rho(IP))
          SinAng(IP)=SIN(A_wav_rho(IP))
        END DO
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, step 25, rnk=', myrank
        WRITE(DBG%FHNDL,*) 'MyRankGlobal=', MyRankGlobal
        FLUSH(DBG%FHNDL)
# endif
        CALL MPI_INTERP_RECV_r8(TheArr_OCNtoWAV_rho, 217, A_wav_rho)
# ifdef DEBUG_WWM
        SumDepReceive=0
# endif
        DO IP=1,MNP
          dep_rho(IP)=A_wav_rho(IP)
# ifdef DEBUG_WWM
          SumDepReceive=SumDepReceive + abs(A_wav_rho(idx))
# endif
        END DO
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'SumDepReceive=', SumDepReceive
        WRITE(DBG%FHNDL,*) 'WAV, ROMS_COUPL_INITIALIZE, WAV, step 33'
        WRITE(DBG%FHNDL,*) 'WAV, rnk=', myrank
        FLUSH(DBG%FHNDL)
        AbsDiff=0
        SumDep1=0
        SumDep2=0
        SumDiff=0
        minBathy=140000
        maxBathy=0
        DO IP=1,MNP
          IF (dep_rho(IP).lt.minBathy) THEN
            minBathy=dep_rho(IP)
          END IF
          IF (dep_rho(IP).gt.maxBathy) THEN
            maxBathy=dep_rho(IP)
          END IF
        END DO
        WRITE(DBG%FHNDL,*) 'dep_rho, min=', minBathy, ' max=', maxBathy
        FLUSH(DBG%FHNDL)
        minBathy=140000
        maxBathy=0
        DO IP=1,MNP
          IF (DEP(IP).lt.minBathy) THEN
            minBathy=DEP(IP)
          END IF
          IF (DEP(IP).gt.maxBathy) THEN
            maxBathy=DEP(IP)
          END IF
        END DO
        WRITE(DBG%FHNDL,*) 'DEP, min=', minBathy, ' max=', maxBathy
        FLUSH(DBG%FHNDL)


        iNodeSel=-1
!        CALL MyGetString(MyRankGlobal, eStrFi)
!        FileSave1='DEP_infos' // eStrFi
!        FileSave2='Lookup_infos' // eStrFi
!        open(745, FILE=TRIM(FileSave1))
!        open(746, FILE=TRIM(FileSave2))
        DO IP=1,MNP
          DEP(IP)=dep_rho(IP)
          eDiff=abs(dep_rho(IP) - DEP(IP))
          SumDiff=SumDiff + eDiff
          SumDep1=SumDep1 + dep_rho(IP)
          SumDep2=SumDep2 + DEP(IP)
          IF (eDiff.gt.AbsDiff) THEN
            AbsDiff=eDiff
            iNodeSel=IP
          END IF
          IF ((DEP(IP).ge.200).and.(eDiff.ge.10)) THEN
            WRITE(DBG%FHNDL,*) 'AD, IP=', IP, dep_rho(IP), DEP(IP)
            WRITE(DBG%FHNDL,*) 'AD, xp, yp=', XP(IP), YP(IP)
            FLUSH(DBG%FHNDL)
          END IF
        END DO
!        close(745)
!        close(746)
!        WRITE(DBG%FHNDL,*) 'AD, AbsDiff=', AbsDiff
!        WRITE(DBG%FHNDL,*) 'AD, IP=', iNodeSel, dep_rho(iNodeSel), DEP(iNodeSel)
!        WRITE(DBG%FHNDL,*) 'AD, xp, yp=', XP(iNodeSel), YP(iNodeSel)
        WRITE(DBG%FHNDL,*) 'AD, SumDep1=', SumDep1, ' SumDep2=', SumDep2
        WRITE(DBG%FHNDL,*) 'AD, SumDiff=', SumDiff
        FLUSH(DBG%FHNDL)
# endif
        allocate(z_r(Nlevel), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_coupl_roms, allocate error 18')
# ifdef FIRST_ORDER_ARDHUIN
        allocate(PartialU1(1), PartialV1(1), stat=istat)
# else
        allocate(PartialU1(Nlevel), PartialV1(Nlevel), PartialU2(Nlevel), PartialV2(Nlevel), stat=istat)
# endif
        IF (istat/=0) CALL WWM_ABORT('wwm_coupl_roms, allocate error 18.1')
        allocate(z_w_wav(0:Nlevel, MNP), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_coupl_roms, allocate error 19')
        allocate(USTOKES_wav(Nlevel,MNP), VSTOKES_wav(Nlevel,MNP), ZETA_CORR(MNP), J_PRESSURE(MNP), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_coupl_roms, allocate error 20')
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'End ROMS_COUPL_INITIALIZE'
        FLUSH(DBG%FHNDL)
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ROMS_COUPL_DEALLOCATE
      USE DATAPOOL, only : MNP, rkind, DEP, XP, YP, np_total, ne_total, istatus, ierr, itype, myrank
      USE mod_coupler
      USE PGMCL_LIBRARY
      USE pgmcl_interp
      implicit none
      logical DoNearest
      integer IP, iNodeSel, idx, eRankRecv
      real(rkind) eDiff, AbsDiff, SumDep1, SumDep2, SumDiff
      real(rkind) minBathy, maxBathy
      real(rkind) SumDepReceive
      CALL WWM_DeallocMatrixPartition
      deallocate(LONtrig_wav, LATtrig_wav, ListTrig_wav)
      deallocate(LON_rho_ocn, LAT_rho_ocn, MSK_rho_ocn)
      deallocate(LON_u_ocn, LAT_u_ocn, MSK_u_ocn)
      deallocate(LON_v_ocn, LAT_v_ocn, MSK_v_ocn)
      deallocate(MatrixBelongingOCN_rho)
      deallocate(MatrixBelongingOCN_u)
      deallocate(MatrixBelongingOCN_v)
      DEALLOCATE(z_w_loc, eUSTOKES_loc, eVSTOKES_loc)
      CALL DEALLOCATE_Arr(TheArr_OCNtoWAV_rho)
      CALL DEALLOCATE_Arr(TheArr_OCNtoWAV_u)
      CALL DEALLOCATE_Arr(TheArr_OCNtoWAV_v)
      CALL DEALLOCATE_Arr(TheArr_WAVtoOCN_rho)
      CALL DEALLOCATE_Arr(TheArr_WAVtoOCN_u)
      CALL DEALLOCATE_Arr(TheArr_WAVtoOCN_v)
      deallocate(CosAng, SinAng, dep_rho)
      deallocate(A_wav_rho_3D, A_wav_stat, A_wav_uvz, A_wav_rho)
      deallocate(A_wav_u_3D, A_wav_v_3D, A_wav_ur_3D, A_wav_vr_3D)
      deallocate(z_r, z_w_wav, U_wav, V_wav)
      deallocate(PartialU1, PartialV1, PartialU2, PartialV2)
      deallocate(USTOKES_wav, VSTOKES_wav)
      deallocate(ZETA_CORR, J_PRESSURE)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE STOKES_STRESS_INTEGRAL_ROMS
        USE mod_coupler
        USE DATAPOOL
        implicit none
        integer IP, k, ID, IS
        real(rkind) eF1, eF2, eDelta, TheInt, eDep, eHeight
        real(rkind) eFrac, eFracB, eQuot, TheIntChk
        real(rkind) eFct, eQuot1, eQuot2, eQuot3, eScal, eZeta
        real(rkind) eOmega, eMult, MFACT, kD, eSinc
        real(rkind) USTOKES1, USTOKES2, USTOKES3
        real(rkind) VSTOKES1, VSTOKES2, VSTOKES3
        real(rkind) USTOKESpart, VSTOKESpart, eJPress
        real(rkind) ACLOC, eWk, eSigma, eLoc, eSinhkd, eSinh2kd, eSinhkd2
        real(rkind) zMid, HS, ETOT, MinVal_MFACT, MaxVal_MFACT, MaxVal_eQuot1
        real(rkind) MinVal_MFACT_gl, MaxVal_MFACT_gl, MaxHS, SumHS, AvgHS
        real(rkind) WLM, KLM, AvgStokesNormA, AvgStokesNormB
        real(rkind) PPTAIL, CETAIL, CKTAIL
        real(rkind) ETOT1, EKTOT
        real(rkind) eQuotDispersion, eMaxAC, TotSumAC, eQuotAC, eQuotK
        real(rkind) StokesNormA, StokesNormB, cPhase
        integer IDsel, ISsel, SelectedK
        logical DoTail
        real(rkind) SumNormStokesA(Nlevel), SumNormStokesB(Nlevel)
        real(rkind) SumZetaCorr, MaxZetaCorr, AvgZetaCorr
        real(rkind) eMinMfact, eMaxMfact, SelectedHS
        real(rkind) MaxStokesNorm, MaxValSinc, StokesNorm, SelectedDEP
        real(rkind) CritError, USTOKES_bar, VSTOKES_bar
        real(rkind) USTOKES_bar_int, VSTOKES_bar_int
        real(rkind) eSum_tot, eSum_tot_int, eWkReal
        real(rkind) eSum_totA, eSum_totA_int
        real(rkind) eSum_totB, eSum_totB_int
        real(rkind) TotalBarotropicErrorUstokes, TotalBarotropicErrorVstokes
        real(rkind) TotalSumUstokes, TotalSumVstokes
        real(rkind) SumHeight
        real(rkind) eJPress_loc, eZetaCorr_loc, eProd, eUint, eVint
# ifndef FIRST_ORDER_ARDHUIN
        eMinMfact=-3
        eMaxMfact=5
# endif
        DO IP=1,MNP
          DO k=1,Nlevel
            z_r(k)=(z_w_wav(k,IP)+z_w_wav(k-1,IP))/2
          END DO
          z_w_loc=z_w_wav(:,IP)
          eDep=z_w_loc(Nlevel)-z_w_loc(0)
# ifdef FIRST_ORDER_ARDHUIN
          PartialU1(1)=(U_wav(2,IP) - U_wav(1,IP))/(z_r(Nlevel)-z_r(Nlevel-1))
          PartialV1(1)=(V_wav(2,IP) - V_wav(1,IP))/(z_r(Nlevel)-z_r(Nlevel-1))
# else
          DO k=2,Nlevel
            PartialU1(k)=(U_wav(k,IP) - U_wav(k-1,IP))/(z_r(k)-z_r(k-1))
            PartialV1(k)=(V_wav(k,IP) - V_wav(k-1,IP))/(z_r(k)-z_r(k-1))
          END DO
          PartialU1(1)=PartialU1(2)
          PartialV1(1)=PartialV1(2)
          DO k=2,Nlevel-1
            !we compute second differential with three values.
            !We have classic formula
            ! d2f/dx2 = (f(x+h) + f(x-h) -2f(x))/h^2
            ! and this is extended to three arbitrary positions
            ! but only first order accuracy.
            eF1=(z_r(k)-z_r(k-1))/(z_r(k+1)-z_r(k-1))
            eF2=(z_r(k+1)-z_r(k))/(z_r(k+1)-z_r(k-1))
            eDelta=(z_r(k) - z_r(k+1))*(z_r(k-1) - z_r(k))
            PartialU2(k)=(U_wav(k+1,IP)*eF1 + U_wav(k-1,IP)*eF2 - U_wav(k,IP))/eDelta
            PartialV2(k)=(V_wav(k+1,IP)*eF1 + V_wav(k-1,IP)*eF2 - V_wav(k,IP))/eDelta
          END DO
          PartialU2(1)=PartialU2(2)
          PartialV2(1)=PartialV2(2)
          PartialU2(Nlevel)=PartialU2(Nlevel-1)
          PartialV2(Nlevel)=PartialV2(Nlevel-1)
# endif
          eUSTOKES_loc=0
          eVSTOKES_loc=0
          eJpress_loc=0
          eZetaCorr_loc=0
# ifdef FIRST_ORDER_ARDHUIN
!todo IS ID ordering
          DO IS=1,MSC
            eMult=SPSIG(IS)*DDIR*DS_INCR(IS)
            eWk=WK(IS,IP)
            kD=MIN(KDMAX, eWk*eDep)
            eWkReal=kD/eDep
            eSinh2kd=MySINH(2*kD)
            eSinhkd=MySINH(kD)
            eSinhkd2=eSinhkd**2
            eSigma=SPSIG(IS)
#  ifdef STOKES_DRIFT_USING_INTEGRAL
            eUint=0
            eVint=0
#  endif
            DO ID=1,MDC
              eLoc=AC2(IS,ID,IP)*eMult
              eScal=COSTH(ID)*PartialU1(1)+SINTH(ID)*PartialV1(1)
              eZeta=eWk/eSinhkd + (eWk/eSigma)*eScal
              eZetaCorr_loc=eZetaCorr_loc + eLoc*eZeta
              eJPress=G9*(kD/eSinh2kd)*(1/eDep) * eLoc
              eJPress_loc=eJPress_loc + eJPress
#  ifdef STOKES_DRIFT_USING_INTEGRAL
              eUint=eUint + eLoc*COSTH(ID)
              eVint=eVint + eLoc*SINTH(ID)
#  endif
            END DO
#  ifdef STOKES_DRIFT_USING_INTEGRAL
            DO k=1,Nlevel
              eFrac=(z_r(k) - z_w_loc(0))/eDep
              eHeight=z_w_loc(k)-z_w_loc(k-1)
              eFracB=eHeight/eDep
              eSinc=SINH(kD*eFracB)/(kD*eFracB)
              eQuot1=eSinc*MyCOSH(2*kD*eFrac)/eSinhkd2
              eProd=eSigma*eWkReal*eQuot1
              eUSTOKES_loc(k)=eUSTOKES_loc(k) + eUint*eProd
              eVSTOKES_loc(k)=eVSTOKES_loc(k) + eVint*eProd
            ENDDO
#  endif
          END DO
# else
!todo ID IS ordering
          DO IS=1,MSC
            eMult=SPSIG(IS)*DDIR*DS_INCR(IS)
            eWk=WK(IS,IP)
            kD=MIN(KDMAX, eWk*eDep)
            eWkReal=kD/eDep
            eSinh2kd=MySINH(2*kD)
            eSinhkd=MySINH(kD)
            eSinhkd2=eSinhkd**2
            eSigma=SPSIG(IS)
            DO ID=1,MDC
              eLoc=AC2(IS,ID,IP)*eMult
              TheInt=0
              DO k=1,Nlevel
                eHeight=z_w_loc(k)-z_w_loc(k-1)
                zMid=0.5*(z_w_loc(k)+z_w_loc(k-1))
                eFrac=(zMid - z_w_loc(0))/eDep
                eFracB=eHeight/eDep
                eSinc=MySINH(eFracB*kD)/(eFracB*kD)
                eQuot=eWkReal*2*MyCOSH(2*kD*eFrac)/eSinh2kd
                eFct=U_wav(k,IP)*COSTH(ID)+V_wav(k,IP)*SINTH(ID)
                TheInt=TheInt+eHeight*eFct*eQuot*eSinc
              END DO
              eOmega=eSigma + TheInt*eWkReal
#  ifdef STOKES_DRIFT_USING_INTEGRAL
              DO k=1,Nlevel
                MFACT=eSigma/(eOmega - (U_wav(k,IP)*COSTH(ID)+V_wav(k,IP)*SINTH(ID))*eWkReal)
                MFACT=MAX(MFACT, eMinMfact)
                MFACT=MIN(MFACT, eMaxMfact)
                eFrac=(z_r(k) - z_w_loc(0))/eDep
                eHeight=z_w_loc(k)-z_w_loc(k-1)
                eFracB=eHeight/eDep
                eSinc=SINH(kD*eFracB)/(kD*eFracB)
                eQuot1=eSinc*MyCOSH(2*kD*eFrac)/eSinhkd2
                USTOKES1=MFACT*eSigma*COSTH(ID)*eWkReal*eQuot1
                VSTOKES1=MFACT*eSigma*SINTH(ID)*eWkReal*eQuot1
                eQuot2=eSinc*MySINH(2*kD*eFrac)/eSinhkd2
                eQuot3=eSinc*(MySINH(kD*eFrac)/eSinhkd)**2
                eScal=PartialU1(k)*COSTH(ID) + PartialV1(k)*SINTH(ID)
                USTOKES2=0.5*(MFACT**2)*COSTH(ID)*eWkReal*eQuot2*eScal
                USTOKES3=0.5*MFACT*PartialU2(k)*eQuot3
                VSTOKES2=0.5*(MFACT**2)*SINTH(ID)*eWkReal*eQuot2*eScal
                VSTOKES3=0.5*MFACT*PartialV2(k)*eQuot3
                USTOKESpart=eLoc*(USTOKES1+USTOKES2+USTOKES3)
                VSTOKESpart=eLoc*(VSTOKES1+VSTOKES2+VSTOKES3)
                eUSTOKES_loc(k)=eUSTOKES_loc(k) + USTOKESpart
                eVSTOKES_loc(k)=eVSTOKES_loc(k) + VSTOKESpart
              ENDDO
#  else
              MFACT=eSigma/(eOmega - (U_wav(Nlevel,IP)*COSTH(ID)+V_wav(Nlevel,IP)*SINTH(ID))*eWkReal)
#  endif
              eScal=COSTH(ID)*PartialU1(Nlevel)+SINTH(ID)*PartialV1(Nlevel)
              eZeta=eWk/eSinhkd + (MFACT*eWk/eSigma)*eScal
              eZetaCorr_loc=eZetaCorr_loc + MFACT*eLoc*eZeta
              eJPress=G9*(kD/eSinh2kd)*(1/eDep) * eLoc
              eJPress_loc=eJPress_loc + eJPress
            END DO
          END DO
# endif
          USTOKES_wav(:,IP)=eUSTOKES_loc
          VSTOKES_wav(:,IP)=eVSTOKES_loc
          ZETA_CORR(IP)=eZetaCorr_loc
          J_PRESSURE(IP)=eJPress_loc
        ENDDO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PGMCL_ROMS_IN(K,IFILE,IT)
        USE pgmcl_library
        USE datapool
        USE mod_coupler
        implicit none
        INTEGER, INTENT(IN) :: K,IFILE,IT
        integer IP, kLev, i, idx
        real(rkind) u1, v1, u2, v2, z1
# ifdef DEBUG_WWM
        real(rkind) :: MaxUwind, SumUwind, avgUwind
        real(rkind) :: MaxVwind, SumVwind, avgVwind
        real(rkind) :: NbPoint
# endif
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WWM: Begin PGMCL_ROMS_IN'
        FLUSH(DBG%FHNDL)
# endif
# ifdef DUMMY_COUPLING
#  ifdef DUMMYB
        CALL MPI_INTERP_DummyB_Recv(ArrLocal, 201, OCNid)
#  else
        CALL MPI_INTERP_Dummy_Recv(ArrLocal, 201, OCNid)
#  endif
# else
#  ifdef NO_ASYNC
        CALL MPI_INTERP_RECV_3D_r8(TheArr_OCNtoWAV_rho, 201, 3, A_wav_uvz)
#  else
        CALL MPI_INTERP_ARECV_3D_r8(TheAsync_OCNtoWAV_uvz, 201, A_wav_uvz)
#  endif
# endif
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_IN, After Data receive'
        FLUSH(DBG%FHNDL)
# endif
# ifdef DEBUG_WWM
        MaxUwind=0.0_r8
        SumUwind=0.0_r8
        MaxVwind=0.0_r8
        SumVwind=0.0_r8
        NbPoint=0
# endif
        WATLEVOLD=WATLEV
        DELTAT_WATLEV = MAIN%DTCOUP
        LCALC=.TRUE.
        DO IP=1,MNP
          u1=A_wav_uvz(1,IP)
          v1=A_wav_uvz(2,IP)
# ifdef DEBUG_WWM
          IF (abs(u1).gt.MaxUwind) THEN
            MaxUwind=abs(u1)
          ENDIF
          IF (abs(v1).gt.MaxVwind) THEN
            MaxVwind=abs(v1)
          ENDIF
          SumUwind=SumUwind + abs(u1)
          SumVwind=SumVwind + abs(v1)
          NbPoint=NbPoint+1
# endif
          u2=u1*CosAng(IP)-v1*SinAng(IP)
          v2=v1*CosAng(IP)+u1*SinAng(IP)
          z1=A_wav_uvz(3,IP)
          WINDXY(IP,1)=u2
          WINDXY(IP,2)=v2
          WATLEV(IP)=z1
        END DO
# ifdef DEBUG_WWM
        avgUwind=SumUwind/NbPoint
        avgVwind=SumVwind/NbPoint
        WRITE(DBG%FHNDL,*) 'WAV, MaxUwind=', MaxUwind, ' avgUwind=', avgUwind
        WRITE(DBG%FHNDL,*) 'WAV, MaxVwind=', MaxVwind, ' avgVwind=', avgVwind
        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_IN, step 2'
        FLUSH(DBG%FHNDL)
# endif
# ifdef DUMMY_COUPLING
#  ifdef DUMMYB
        CALL MPI_INTERP_DummyB_Recv(ArrLocal, 203, OCNid)
#  else
        CALL MPI_INTERP_Dummy_Recv(ArrLocal, 203, OCNid)
#  endif
# else
#  ifdef NO_ASYNC
        CALL MPI_INTERP_RECV_3D_r8(TheArr_OCNtoWAV_rho, 203, Nlevel+1, A_wav_rho_3D)
#  else
        CALL MPI_INTERP_ARECV_3D_r8(TheAsync_OCNtoWAV_rho, 203, A_wav_rho_3D)
#  endif
# endif
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_IN, step 3'
        FLUSH(DBG%FHNDL)
# endif
        DO IP=1,MNP
          DO kLev=0,Nlevel
            z_w_wav(kLev,IP)=A_wav_rho_3D(kLev+1,IP)
          END DO
        END DO
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_IN, step 4'
        FLUSH(DBG%FHNDL)
# endif
# ifdef DUMMY_COUPLING
#  ifdef DUMMYB
        CALL MPI_INTERP_DummyB_Recv(ArrLocal, 204, OCNid)
#  else
        CALL MPI_INTERP_Dummy_Recv(ArrLocal, 204, OCNid)
#  endif
# else
#  ifdef NO_ASYNC
#   ifdef FIRST_ORDER_ARDHUIN
        CALL MPI_INTERP_RECV_3D_r8(TheArr_OCNtoWAV_u, 204, 2, A_wav_ur_3D)
#   else
        CALL MPI_INTERP_RECV_3D_r8(TheArr_OCNtoWAV_u, 204, Nlevel, A_wav_ur_3D)
#   endif
#  else
        CALL MPI_INTERP_ARECV_3D_r8(TheAsync_OCNtoWAV_u, 204, A_wav_ur_3D)
#  endif
# endif
        DO IP=1,MNP
          U_wav(:,IP)=A_wav_ur_3D(:,IP)
        END DO
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_IN, step 5'
        FLUSH(DBG%FHNDL)
# endif
# ifdef DUMMY_COUPLING
#  ifdef DUMMYB
        CALL MPI_INTERP_DummyB_Recv(ArrLocal, 205, OCNid)
#  else
        CALL MPI_INTERP_Dummy_Recv(ArrLocal, 205, OCNid)
#  endif
# else
#  ifdef NO_ASYNC
#   ifdef FIRST_ORDER_ARDHUIN
        CALL MPI_INTERP_RECV_3D_r8(TheArr_OCNtoWAV_v, 205, 2, A_wav_vr_3D)
#   else
        CALL MPI_INTERP_RECV_3D_r8(TheArr_OCNtoWAV_v, 205, Nlevel, A_wav_vr_3D)
#   endif
#  else
        CALL MPI_INTERP_ARECV_3D_r8(TheAsync_OCNtoWAV_v, 205, A_wav_vr_3D)
#  endif
# endif
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'After the receive'
        FLUSH(DBG%FHNDL)
# endif
        DO IP=1,MNP
          V_wav(:,IP)=A_wav_vr_3D(:,IP)
        END DO
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_IN, step 6'
        FLUSH(DBG%FHNDL)
# endif
        DO IP=1,MNP
#ifdef FIRST_ORDER_ARDHUIN
          DO kLev=1,2
#else
          DO kLev=1,Nlevel
#endif
            u1=U_wav(kLev,IP)
            v1=V_wav(kLev,IP)
            u2=u1*CosAng(IP)-v1*SinAng(IP)
            v2=v1*CosAng(IP)+u1*SinAng(IP)
            U_wav(kLev,IP)=u2
            V_wav(kLev,IP)=v2
          END DO
          CURTXY(IP,1)=u2
          CURTXY(IP,2)=v2
        END DO
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_IN, step 7'
        FLUSH(DBG%FHNDL)
# endif
      END SUBROUTINE PGMCL_ROMS_IN
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PGMCL_ROMS_OUT(K)
        USE DATAPOOL
        USE pgmcl_library
        USE mod_coupler
# ifdef ST41
        USE W3SRC4MD_OLD, only : Z0, TAUWX, TAUWY, UFRIC, CD
# endif
# ifdef ST42
        USE W3SRC4MD_NEW, only : Z0, TAUWX, TAUWY, UFRIC, CD
# endif
        implicit none
        INTEGER, INTENT(IN)  :: K
        integer IP, kLev, idx
        real(rkind) u1, v1, u2, v2
        real(rkind) HS, TM01, TM02, KLM, WLM, TM10
        real(rkind) UBOT, ORBITAL, BOTEXPER, TMBOT
        real(rkind) FPP, TPP, CPP, WNPP, CGPP, KPP, LPP, PEAKDSPR, PEAKDM, DPEAK
        real(rkind) ETOTS,ETOTC,DM,DSPR
        REAL(RKIND) :: ACLOC(MSC,MDC)
        real(rkind) cPhase, eStokesNorm
        real(rkind) kD
        real(rkind) :: TPPD, KPPD, CGPD, CPPD
# ifdef DEBUG_WWM
        real(rkind) SumNormTau, MaxNormTau, AvgNormTau, eNorm
        real(rkind) AvgUFRICsqr, SumUFRICsqr
        real(rkind) AvgCd, SumCd
        real(rkind) AvgStressCd, SumStressCd, eStressCd, eMag
        real(rkind) AvgAlpha, SumAlpha, eAlpha, NbAlpha
        real(rkind) SumWind, AvgWind
        real(rkind) :: MaxHwave, SumHwave, avgHwave, NbPoint
        real(rkind) :: MaxLwave, SumLwave, avgLwave
        real(rkind) :: MaxTM02, SumTM02, AvgTM02
        real(rkind) :: MaxStokesNorm, SumStokesNorm, avgStokesNorm
        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_OUT, step 1'
        FLUSH(DBG%FHNDL)
        SumNormTau=0
        MaxNormTau=0
# endif
        CALL STOKES_STRESS_INTEGRAL_ROMS
        DO IP=1,MNP
# ifdef STOKES_DRIFT_USING_INTEGRAL
          DO kLev=1,Nlevel
            u1=USTOKES_wav(kLev,IP)
            v1=VSTOKES_wav(kLev,IP)
            u2=u1*CosAng(IP)+v1*SinAng(IP)
            v2=v1*CosAng(IP)-u1*SinAng(IP)
            A_wav_u_3D(kLev,IP)=u2
            A_wav_v_3D(kLev,IP)=v2
          END DO
# endif
          u1=TAUWX(IP)
          v1=TAUWY(IP)
          u2=u1*CosAng(IP)+v1*SinAng(IP)
          v2=v1*CosAng(IP)-u1*SinAng(IP)
# ifdef STOKES_DRIFT_USING_INTEGRAL
          A_wav_u_3D(Nlevel+1,IP)=u2
          A_wav_v_3D(Nlevel+1,IP)=v2
# else
          A_wav_u_3D(1,IP)=u2
          A_wav_v_3D(1,IP)=v2
# endif
# ifdef DEBUG_WWM
          eNorm=SQRT(u2*u2 + v2*v2)
          IF (eNorm.gt.MaxNormTau) THEN
            MaxNormTau=eNorm
          END IF
          SumNormTau=SumNormTau + eNorm
# endif
        END DO
# ifdef DEBUG_WWM
        AvgNormTau=SumNormTau / MNP
        WRITE(DBG%FHNDL,*) 'AvgNormTau=', AvgNormTau, 'MaxNormTau=', MaxNormTau
        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_OUT, step 5.1'
        FLUSH(DBG%FHNDL)
# endif
# ifdef DUMMY_COUPLING
#  ifdef DUMMYB
        CALL MPI_INTERP_DummyB_Send(ArrLocal, 208, OCNid)
        CALL MPI_INTERP_DummyB_Send(ArrLocal, 210, OCNid)
        CALL MPI_INTERP_DummyB_Send(ArrLocal, 209, OCNid)
        CALL MPI_INTERP_DummyB_Send(ArrLocal, 211, OCNid)
#  else
        CALL MPI_INTERP_Dummy_Send(ArrLocal, 208, OCNid)
        CALL MPI_INTERP_Dummy_Send(ArrLocal, 210, OCNid)
        CALL MPI_INTERP_Dummy_Send(ArrLocal, 209, OCNid)
        CALL MPI_INTERP_Dummy_Send(ArrLocal, 211, OCNid)
#  endif
# else
#  ifdef STOKES_DRIFT_USING_INTEGRAL
#   ifdef NO_ASYNC
        CALL MPI_INTERP_SEND_3D_r8(TheArr_WAVtoOCN_u, 208, Nlevel+1, A_wav_u_3D)
        CALL MPI_INTERP_SEND_3D_r8(TheArr_WAVtoOCN_v, 210, Nlevel+1, A_wav_u_3D)
        CALL MPI_INTERP_SEND_3D_r8(TheArr_WAVtoOCN_u, 209, Nlevel+1, A_wav_v_3D)
        CALL MPI_INTERP_SEND_3D_r8(TheArr_WAVtoOCN_v, 211, Nlevel+1, A_wav_v_3D)
#   else
        CALL MPI_INTERP_ASEND_3D_r8(TheAsync_WAVtoOCN_u, 208, A_wav_u_3D)
        CALL MPI_INTERP_ASEND_3D_r8(TheAsync_WAVtoOCN_v, 210, A_wav_u_3D)
        CALL MPI_INTERP_ASEND_3D_r8(TheAsync_WAVtoOCN_u, 209, A_wav_v_3D)
        CALL MPI_INTERP_ASEND_3D_r8(TheAsync_WAVtoOCN_v, 211, A_wav_v_3D)
#   endif
#  else
#   ifdef NO_ASYNC
        CALL MPI_INTERP_SEND_3D_r8(TheArr_WAVtoOCN_u, 208, 1, A_wav_u_3D)
        CALL MPI_INTERP_SEND_3D_r8(TheArr_WAVtoOCN_v, 211, 1, A_wav_v_3D)
#   else
        CALL MPI_INTERP_ASEND_3D_r8(TheAsync_WAVtoOCN_u, 208, A_wav_u_3D)
        CALL MPI_INTERP_ASEND_3D_r8(TheAsync_WAVtoOCN_v, 211, A_wav_v_3D)
#   endif
#  endif
# endif
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_OUT, step 5.3'
        FLUSH(DBG%FHNDL)
# endif
# ifdef DEBUG_WWM
        MaxHwave=0.0
        SumHwave=0.0
        MaxTM02=0
        SumTM02=0
        MaxLwave=0.0
        SumLwave=0.0
        SumStokesNorm=0
        MaxStokesNorm=0
        NbPoint=0.0
        SumUFRICsqr=0
        SumCd=0
        SumWind=0
        SumStressCd=0
        SumAlpha=0
        NbAlpha=0
# endif
        DO IP = 1, MNP
          ACLOC = AC2(:,:,IP)
          CALL MEAN_PARAMETER(IP,ACLOC,MSC,HS,TM01,TM02,TM10,KLM,WLM)
          CALL WAVE_CURRENT_PARAMETER(IP,ACLOC,UBOT,ORBITAL,BOTEXPER,TMBOT,'PGMCL_ROMS_OUT')
          CALL MEAN_DIRECTION_AND_SPREAD(IP,ACLOC,MSC,ETOTS,ETOTC,DM,DSPR)
          CALL PEAK_PARAMETER(IP,ACLOC,MSC,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD)
# ifdef DEBUG_WWM
          SumUFRICsqr=SumUFRICsqr + UFRIC(IP)*UFRIC(IP)
          eMag=SQRT(WINDXY(IP,1)**2 + WINDXY(IP,2)**2)
          eStressCd=CD(IP)*eMag*eMag
          SumWind=SumWind + eMag
          SumStressCd=SumStressCd + eStressCd
          IF (UFRIC(IP).gt.0) THEN
            eAlpha=G9*Z0(IP)/(UFRIC(IP) * UFRIC(IP))
            SumAlpha=SumAlpha + eAlpha
            NbAlpha=NbAlpha+1
          END IF
          IF (HS.gt.MaxHwave) THEN
            MaxHwave=HS
          ENDIF
          SumHwave=SumHwave + HS
          IF (TM02.gt.MaxTM02) THEN
            MaxTM02=TM02
          ENDIF
          SumTM02=SumTM02 + TM02
          IF (WLM.gt.MaxLwave) THEN
            MaxLwave=WLM
          ENDIF
          SumLwave=SumLwave + WLM
          kD=MIN(KDMAX, KLM*DEP(IP))
          cPhase=SQRT((G9/KLM)*REAL(MySINH(kD)/MyCOSH(kD)) )
          eStokesNorm=(G9*HS*HS/REAL(16))*2*(KLM/cPhase)
          IF (eStokesNorm.ne.eStokesNorm) THEN
            WRITE(DBG%FHNDL,*) 'eStokesNorm=', eStokesNorm
            WRITE(DBG%FHNDL,*) 'KLM=', KLM, 'WLM=', WLM
            WRITE(DBG%FHNDL,*) 'cPhase=', cPhase, 'kD=', kD
            WRITE(DBG%FHNDL,*) 'HS=', HS, ' DEP=', DEP(IP)
            FLUSH(DBG%FHNDL)
          END IF
          IF (eStokesNorm.gt.MaxStokesNorm) THEN
            MaxStokesNorm=eStokesNorm
          END IF
          SumStokesNorm=SumStokesNorm + eStokesNorm
          NbPoint=NbPoint + 1
# endif
          A_wav_stat(1, IP)=HS
          A_wav_stat(2, IP)=TM01
          A_wav_stat(3, IP)=TM02
          A_wav_stat(4, IP)=KLM
          A_wav_stat(5, IP)=WLM
          A_wav_stat(6, IP)=ORBITAL
          A_wav_stat(7, IP)=TMBOT
          A_wav_stat(8, IP)=DISSIPATION(IP)
          A_wav_stat(9, IP)=QBLOCAL(IP)
          A_wav_stat(10,IP)=DM
          A_wav_stat(11,IP)=TPP
          A_wav_stat(12,IP)=DSPR
          A_wav_stat(13,IP)=PEAKDSPR
          A_wav_stat(14,IP)=PEAKDM
          A_wav_stat(15,IP)=UFRIC(IP)
          A_wav_stat(16,IP)=Z0(IP)
          A_wav_stat(17,IP)=CD(IP)
          A_wav_stat(18,IP)=J_PRESSURE(IP)
          A_wav_stat(19,IP)=ZETA_CORR(IP)
        END DO
# ifdef DEBUG_WWM
        avgHwave=SumHwave/NbPoint
        avgLwave=SumLwave/NbPoint
        AvgTM02=SumTM02/NbPoint
        avgStokesNorm=SumStokesNorm/NbPoint
        WRITE(DBG%FHNDL,*) 'WAV, MaxHwave=', MaxHwave, ' avgHwave=', avgHwave
        WRITE(DBG%FHNDL,*) 'WAV, MaxLwave=', MaxLwave, ' avgLwave=', avgLwave
        WRITE(DBG%FHNDL,*) 'WAV, MaxStokesNorm=', MaxStokesNorm
        WRITE(DBG%FHNDL,*) 'WAV, avgStokesNorm=', avgStokesNorm
        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_OUT, step 6'
        FLUSH(DBG%FHNDL)
        AvgUFRICsqr=SumUFRICsqr/MNP
        AvgStressCd=SumStressCd/MNP
        AvgAlpha=SumAlpha/NbAlpha
        AvgWind=SumWind/MNP
        AvgCd=SumCd/MNP
        WRITE(DBG%FHNDL,*) 'AvgNormTau=', AvgNormTau, 'AvgNormFV2=', AvgUFRICsqr
        WRITE(DBG%FHNDL,*) 'AvgNormTau=', AvgNormTau, 'AvgCdU2=', AvgStressCd
        WRITE(DBG%FHNDL,*) 'AvgCd=', AvgCd, ' AvgAlpha=', AvgAlpha
        WRITE(DBG%FHNDL,*) 'AvgWind=', AvgWind
        FLUSH(DBG%FHNDL)
# endif
# ifdef DUMMY_COUPLING
#  ifdef DUMMYB
        CALL MPI_INTERP_DummyB_Send(ArrLocal, 212, OCNid)
#  else
        CALL MPI_INTERP_Dummy_Send(ArrLocal, 212, OCNid)
#  endif
# else
#  ifdef NO_ASYNC
        CALL MPI_INTERP_SEND_3D_r8(TheArr_WAVtoOCN_rho, 212, 19, A_wav_stat)
#  else
        CALL MPI_INTERP_ASEND_3D_r8(TheAsync_WAVtoOCN_stat, 212, A_wav_stat)
#  endif
# endif
# ifdef DEBUG_WWM
        WRITE(DBG%FHNDL,*) 'WWM: PGMCL_ROMS_OUT, step 11'
        FLUSH(DBG%FHNDL)
# endif
      END SUBROUTINE
#endif
