#include "wwm_functions.h"
! I5 is under the assumptions that we have a lot of memory
!    and tries to minimize the number of MPI exchanges and
!    to have the processors be as busy as possible
!**********************************************************************
!* We have to think on how the system is solved. Many questions are   *
!* mixed: the ordering of the nodes, the ghost nodes, the aspar array *
!* Here is a repository of the conclusions that have been reached     *
!*                                                                    *
!* Ordering 1> should be that way: We have two global nodes i and j.  *
!* -- If i and j belong to a common local grid, then we select the    *
!*    grid of lowest color and decide whether ipgl(i) < ipgl(j)       *
!* -- If i and j belong to two different grid then                    *
!*     ---If Color(i) < Color(j) or reverse we decide by that         *
!*     ---If Color(i) = Color(j) we decide by i<j or not (but it      *
!*        does not matter to the solution)                            *
!* The functions WRITE_EXPLICIT_ORDERING does exactly that and        *
!* provides an ordering that can be used. That is we start with the   *
!* nodes of lowest color and index until we arrive at highest color   *
!*                                                                    *
!* The ASPAR is computed correctly only on 1:NP_RES but this can be   *
!* extended by exchange routines.                                     *
!* We Compute on the resident nodes only. This means loops over       *
!* IP=1,NP_RES and backwards. This means that we do not have to do    *
!* exchanges of ASPAR values. Only the resident nodes are sent.       *
!* This is smaller and this is all that we ever need.                 *
!*                                                                    *
!* WRONG APPROACHES:                                                  *
!* to use all nodes 1:MNP may look simpler but it forces to have the  *
!* following property of the NNZ, IA, JA arrays. If two vertices      *
!* i and j are adjacent in a grid G, then they are adjacent in ANY    *
!* of the grid in which they are both contained.                      *
!* This property is actually not satisfied in general.                *
!* We could extend the IA, JA arrays by                               *
!* adding some vertices but that looks quite hazardous idea and it    *
!* it is actually not needed by the ILU0 preconditioner and other     *
!*                                                                    *
!* PROBLEM:                                                           *
!* There is an asymmetry in the construction of the ordering.         *
!* We start from low colors and upwards. If we had started with       *
!* high colors and gone downwards, then we get a different ordering   *
!* (even if we take the opposite, because of the interface nodes)     *
!* This requires the construction of many mappings.                   *
!* Our approach is actually to rebuild separate node sets.            *
!*                                                                    *
!* CHECKS:                                                            *
!* ---The sum of number of non-zero entries in Jstatus_L over all     *
!*    nodes should be equal to the sum of number of non-zero entries  *
!*    of Jstatus_U over all nodes.                                    *
!*    This is because number of upper diagonal entries should be      *
!*    equal to number of lower diagonal entries.                      *
!* ---We CANNOT have Jstatus_L(J)=1 and Jstatus_U(J)=1, i.e. a matrix *
!*    entry cannot be both lower and upper.                           *
!* ---We have sum of Jstatus_L + sum J_status_U + np_global should    *
!*    be equal to NNZ_global                                          *
!*                                                                    *
!* So, procedure is as follows:                                       *
!* ---compute ASPAR on 1,NP_RES nodes and no synchronization          *
!* ---compute on IP=1,NP_RES for L solving                            *
!* ---export to grids of higher rank, the values on nodes 1,NP_RES    *
!*    only. The other ghost points have invalid values or are         *
!*    resident of other grids of lower rank. (at this stage, some     *
!*    ghost values are wrong but are not exported)                    *
!* ---export to grid of lower rank in order to correct their ghost    *
!*    values and get the value                                        *
!* ---compute on IP=NP_RES,1,-1 for U solving                         *
!* ---export to grid of lower rank.                                   *
!*    Do everything similarly to L solve.                             *
!*                                                                    *
!* The basic approach is that we compute at a node S if and only if   *
!* it is a resident node. The twist come because some nodes are       *
!* resident for TWO domains. This is why we have the CovLower         *
!* We need to create disjoint domains for each node, so that          *
!* we have sum   sum(CovLower) = np_global                            *
!*                                                                    *
!* At the end of the resolution of the system, we need to do the      *
!* synchronization with respect to the unused nodes.                  *
!* Mystery?: When we apply the function, we do it on 1:NP_RES and     *
!* then call synchronizer. So, this means we need to do the           *
!* synchronization after the call to the preconditioner.              *
!* But there may be space for improvements here.                      *
!*                                                                    *
!* Description of specific exchange arrays:                           *
!* ---wwm_p2dsend_type/wwm_p2drecv_type                               *
!*    The points of 1:NP_RES are sent to nodes that contained them    *
!*    length=1                                                        *
!* ---wwmtot_p2dsend_type/wwmtot_p2drecv_type                         *
!*    same as above but length=MSC*MDC                                *
!* ---blk_p2dsend_type/blk_p2drecv_type                               *
!*    same as above but length=maxBlockLength for matrix exchanges    *
!* ---wwmmat_p2dsend_type/wwmmat_p2drecv_type                         *
!*    exchange of correct matrix elements, i.e. elements A(I,J)       *
!*    with I<=NP_RES                                                  *
!*    length is 1.                                                    *
!* ---u2l_p2dsend_type/u2l_p2drecv_type                               *
!*    upper 2 lower exchange arrays, depends on CovLower, so on       *
!*    the coloring chosen. length=maxBlockLength                      *
!*    exchange are from upper to lower.                               *
!* ---sync_p2dsend_type/sync_p2drecv_type                             *
!*    synchronize value, i.e. the CovLower=1 values  are send to all  *
!*    nodes. Length is MSC*MDC                                        *
!**********************************************************************
!* Mystery of the mpi_exchange routines.                              *
!* One way to have derived data types is to do                        *
!* call mpi_type_create_indexed_block(nbCommon,1,dspl_send,           *
!*                                    rtype,eType,ierr)               *
!* call mpi_type_commit(eType, ierr)                                  *
!* here 1 refers to the block length.                                 *
!* The problem is how to make exchanges. The standard method          *
!* mpi_isend(eMes, 1, eType, ....)                                    *
!* and so to send only 1 vector. A priori, it seems impossible        *
!* to send several vectors together because when we created the       *
!* derived data type, the length was not precised in this creation    *
!* one way could be to use mpi_type_resized, but this seems quite     *
!* inoperative or create seg-fault. So, we do not know how to make    *
!* exchanges of the form mpi_isend(AC, MSC*MDC, eType, ....)          *
!* and this means that we cannot send AC(MNP, MSC, MDC) simply.       *
!* Instead, we have to use U(MSC,MDC,MNP) in order to use the above   *
!**********************************************************************
MODULE WWM_PARALL_SOLVER
#if defined WWM_SOLVER && defined MPI_PARALL_GRID
      TYPE Graph
         integer nbVert
         integer MaxDeg
         integer nbEdge
         integer, dimension(:), pointer :: ListDegree
         integer, dimension(:,:), pointer :: ListEdge
      END TYPE Graph
      CONTAINS
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE EXCHANGE_P4D_WWM_TR(AC)
      USE DATAPOOL, only : MNP, MSC, MDC, rkind
      USE elfe_msgp, only : exchange_p4d_wwm
      implicit none
      real(rkind), intent(inout) :: AC(MNP,MSC,MDC)
      real(rkind) :: U(MSC,MDC,MNP)
      INTEGER :: IP, IS, ID
      DO IP = 1, MNP
        DO IS = 1, MSC
          DO ID = 1, MDC
            U(IS,ID,IP) = AC(IP,IS,ID)
          END DO
        END DO
      END DO
      CALL EXCHANGE_P4D_WWM(U)
      DO IP = 1, MNP
        DO IS = 1, MSC
          DO ID = 1, MDC
            AC(IP,IS,ID) = U(IS,ID,IP)
          END DO
        END DO
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE EXCHANGE_ASPAR(ASPAR_bl)
      USE DATAPOOL
      USE elfe_msgp, only : comm, ierr, myrank
      implicit none
      integer I, iProc
      real(rkind), intent(inout) :: ASPAR_bl(NNZ,MSC,MDC)
      real(rkind) :: U(MSC,MDC,NNZ)
      INTEGER :: IZ, IS, ID
      DO IZ = 1, NNZ
        DO IS = 1, MSC
          DO ID = 1, MDC
            U(IS,ID,IZ) = ASPAR_bl(IZ,IS,ID)
          END DO
        END DO
      END DO
# ifdef DEBUG
      WRITE(myrank+740,*) 'wwm_nnbr_m_send=', wwm_nnbr_m_send
# endif
      do I=1,wwm_nnbr_m_send
        iProc=wwm_ListNeigh_m_send(I)
# ifdef DEBUG
        WRITE(myrank+740,*) 'I=', I, 'iProc=', iProc
# endif
        call mpi_isend(U,1,wwmmat_p2dsend_type(I),iProc-1,991,comm,wwmmat_p2dsend_rqst(i),ierr)
      enddo
# ifdef DEBUG
      WRITE(myrank+740,*) 'wwm_nnbr_m_recv=', wwm_nnbr_m_recv
# endif
      do I=1,wwm_nnbr_m_recv
        iProc=wwm_ListNeigh_m_recv(I)
# ifdef DEBUG
        WRITE(myrank+740,*) 'I=', I, 'iProc=', iProc
# endif
        call mpi_irecv(U,1,wwmmat_p2drecv_type(I),iProc-1,991,comm,wwmmat_p2drecv_rqst(i),ierr)
      enddo
# ifdef DEBUG
      WRITE(myrank+740,*) 'wwm_nnbr_m_recv=', wwm_nnbr_m_recv
# endif
      IF (wwm_nnbr_m_recv .gt. 0) THEN
        call mpi_waitall(wwm_nnbr_m_recv,wwmmat_p2drecv_rqst,wwmmat_p2drecv_stat,ierr)
      END IF
# ifdef DEBUG
      WRITE(myrank+740,*) 'wwm_nnbr_m_send=', wwm_nnbr_m_send
# endif
      IF (wwm_nnbr_m_send .gt. 0) THEN
        call mpi_waitall(wwm_nnbr_m_send,wwmmat_p2dsend_rqst,wwmmat_p2dsend_stat,ierr)
# ifdef DEBUG
        WRITE(myrank+740,*) 'ierr=', ierr
# endif
      END IF
      DO IZ = 1, NNZ
        DO IS = 1, MSC
          DO ID = 1, MDC
            ASPAR_bl(IZ,IS,ID) = U(IS,ID,IZ)
          END DO
        END DO
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SYMM_GRAPH_BUILD_PROCESSOR_ADJACENCY(AdjGraph)
      USE DATAPOOL, only : wwm_nnbr, wwm_ListNeigh
      USE elfe_msgp, only : myrank
      implicit none
      type(Graph), intent(inout) :: AdjGraph
# ifdef DEBUG
      WRITE(myrank+740,*) 'wwm_nnbr=', wwm_nnbr
# endif
      CALL KERNEL_GRAPH_BUILD_PROCESSOR_ADJACENCY(AdjGraph, wwm_nnbr, wwm_ListNeigh)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GRAPH_BUILD_PROCESSOR_ADJACENCY(AdjGraph)
      USE elfe_msgp, only : nnbr_p, nbrrank_p
      implicit none
      type(Graph), intent(inout) :: AdjGraph
      integer :: ListNe(nnbr_p)
      integer I
      DO I=1,nnbr_p
        ListNe(I)=nbrrank_p(I)+1
      END DO
      CALL KERNEL_GRAPH_BUILD_PROCESSOR_ADJACENCY(AdjGraph, nnbr_p, ListNe)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE KERNEL_GRAPH_BUILD_PROCESSOR_ADJACENCY(AdjGraph, nb, ListNe)
      USE elfe_msgp, only : myrank, nproc, ierr, comm, istatus, itype
      USE DATAPOOL
      implicit none
      integer, intent(in) :: nb
      integer, intent(in) :: ListNe(nb)
      type(Graph), intent(inout) :: AdjGraph
      integer, allocatable :: rbuf_int(:)
      integer ierror, I, iProc
      integer idx, eDeg, nbEdge, iEdge
      AdjGraph % nbVert=nproc
# ifdef DEBUG
      WRITE(myrank+740,*) 'myrank=', myrank, ' KERNEL_GRAPH_BUILD_PROC...'
# endif
      IF (myrank.eq.0) THEN
        allocate(AdjGraph % ListDegree(nproc))
        AdjGraph % ListDegree(1)=nb
        allocate(rbuf_int(1))
        DO iProc=2,nproc
          CALL MPI_RECV(rbuf_int,1,itype, iProc-1, 19, comm, istatus, ierror)
          AdjGraph % ListDegree(iProc)=rbuf_int(1)
        END DO
        deallocate(rbuf_int)
        nbEdge=0
        DO iProc=1,nproc
          nbEdge=nbEdge + AdjGraph % ListDegree(iProc)
        END DO
        AdjGraph % nbEdge=nbEdge
# ifdef DEBUG
        WRITE(myrank+740,*) 'nbEdge=', nbEdge
# endif
        allocate(AdjGraph % ListEdge(nbEdge,2))
        idx=0
        eDeg=AdjGraph % ListDegree(1)
        DO I=1,eDeg
          idx=idx+1
          AdjGraph % ListEdge(idx,1)=1
          AdjGraph % ListEdge(idx,2)=ListNe(I)
        END DO
        DO iProc=2,nproc
          eDeg=AdjGraph % ListDegree(iProc)
          allocate(rbuf_int(eDeg))
          CALL MPI_RECV(rbuf_int,eDeg,itype, iProc-1, 24, comm, istatus, ierr)
          DO I=1,eDeg
            idx=idx+1
            AdjGraph % ListEdge(idx,1)=iProc
            AdjGraph % ListEdge(idx,2)=rbuf_int(I)
          END DO
          deallocate(rbuf_int)
        END DO
        allocate(rbuf_int(1))
        rbuf_int(1)=nbEdge
        DO iProc=2,nproc
          CALL MPI_SEND(rbuf_int,1,itype, iProc-1, 30, comm, ierr)
        END DO
        deallocate(rbuf_int)
        !
        allocate(rbuf_int(nproc))
        DO iProc=1,nproc
          rbuf_int(iProc)=AdjGraph % ListDegree(iProc)
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(rbuf_int,nproc,itype, iProc-1, 32, comm, ierr)
        END DO
        deallocate(rbuf_int)
        !
        allocate(rbuf_int(nbEdge*2))
        DO iEdge=1,nbEdge
          rbuf_int(2*(iEdge-1)+1)=AdjGraph % ListEdge(iEdge,1)
          rbuf_int(2*(iEdge-1)+2)=AdjGraph % ListEdge(iEdge,2)
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(rbuf_int,2*nbEdge,itype, iProc-1, 34, comm, ierr)
        END DO
# ifdef DEBUG
        WRITE(myrank+740,*) 'ListEdge exported'
# endif
      ELSE
        allocate(rbuf_int(1))
        rbuf_int(1)=nb
        CALL MPI_SEND(rbuf_int,1,itype, 0, 19, comm, ierr)
        deallocate(rbuf_int)
        !
        allocate(rbuf_int(nb))
        DO I=1,nb
          rbuf_int(I)=ListNe(I)
        END DO
        CALL MPI_SEND(rbuf_int,nb,itype, 0, 24, comm, ierr)
        deallocate(rbuf_int)
        !
        allocate(rbuf_int(1))
        CALL MPI_RECV(rbuf_int,1,itype, 0, 30, comm, istatus, ierr)
        nbEdge=rbuf_int(1)
        deallocate(rbuf_int)
        AdjGraph % nbEdge=nbEdge
# ifdef DEBUG
        WRITE(myrank+740,*) 'nbEdge=', nbEdge
# endif
        !
        allocate(rbuf_int(nproc))
        CALL MPI_RECV(rbuf_int,nproc,itype, 0, 32, comm, istatus, ierr)
        allocate(AdjGraph % ListDegree(nproc))
        AdjGraph % ListDegree=rbuf_int
        deallocate(rbuf_int)
# ifdef DEBUG
        WRITE(myrank+740,*) 'ListDegree assigned'
# endif
        !
        allocate(rbuf_int(2*nbEdge))
        CALL MPI_RECV(rbuf_int,2*nbEdge,itype, 0, 34, comm, istatus, ierr)
        allocate(AdjGraph % ListEdge(nbEdge,2))
        DO iEdge=1,nbEdge
          AdjGraph % ListEdge(iEdge,1)=rbuf_int(2*(iEdge-1)+1)
          AdjGraph % ListEdge(iEdge,2)=rbuf_int(2*(iEdge-1)+2)
        END DO
        deallocate(rbuf_int)
# ifdef DEBUG
        WRITE(myrank+740,*) 'ListEdge assigned'
# endif
      ENDIF
      AdjGraph % MaxDeg=maxval(AdjGraph % ListDegree)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GRAPH_TEST_CONNECT(AdjGraph, result)
      implicit none
      type(Graph), intent(in) :: AdjGraph
      integer, intent(out) :: result
      integer, allocatable :: ListStatus(:)
      integer, allocatable :: ListPosFirst(:)
      integer idx, iVert, nbVert, nbVertIsFinished, eAdj
      integer eDeg, sizConn, I, IsFinished
      idx=0
      nbVert=AdjGraph%nbVert
      allocate(ListPosFirst(nbVert))
      allocate(ListStatus(nbVert))
      ListStatus=0
      DO iVert=1,nbVert
        ListPosFirst(iVert)=idx
        idx=idx+AdjGraph % ListDegree(iVert)
      END DO
      ListStatus(1)=1
      DO
        IsFinished=1
        DO iVert=1,nbVert
          IF (ListStatus(iVert) == 1) THEN
            eDeg=AdjGraph%ListDegree(iVert)
            idx=ListPosFirst(iVert)
            DO I=1,eDeg
              eAdj=AdjGraph%ListEdge(idx+I,1)
              IF (ListStatus(eAdj) == 0) THEN
                IsFinished=0
              END IF
              ListStatus(eAdj)=1
            END DO
          END IF
        END DO
        IF (IsFinished == 1) THEN
          EXIT
        END IF
      END DO
      sizConn=0
      DO iVert=1,nbVert
        IF (ListStatus(iVert) == 1) THEN
          sizConn=sizConn+1
        END IF
      END DO
      IF (sizConn == nbVert) THEN
        result=1
      END IF
      result=0
      deallocate(ListStatus)
      deallocate(ListPosFirst)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GRAPH_TEST_UNIDIRECT(AdjGraph)
      implicit none
      type(Graph), intent(in) :: AdjGraph
      integer nbEdge, iEdge, jEdge, jEdgeF
      integer eVert1, eVert2, fVert1, fVert2
      nbEdge=AdjGraph%nbEdge
      DO iEdge=1,nbEdge
        eVert1=AdjGraph%ListEdge(iEdge,1)
        eVert2=AdjGraph%ListEdge(iEdge,2)
        jEdgeF=-1
        DO jEdge=1,nbEdge
          fVert1=AdjGraph%ListEdge(jEdge,1)
          fVert2=AdjGraph%ListEdge(jEdge,2)
          IF ((eVert1.eq.fVert2).and.(eVert2.eq.fVert1)) THEN
            jEdgeF=jEdge
          END IF
        END DO
        IF (jEdgeF.eq.-1) THEN
          CALL WWM_ABORT('Error at test of symmetry')
        END IF
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COLLECT_ALL_IPLG
      USE DATAPOOL
      USE elfe_msgp
      USE elfe_glbl
      implicit none
      integer, allocatable :: rbuf_int(:)
      integer len, iProc, IP, idx, sumMNP
      allocate(ListMNP(nproc))
      allocate(ListNP_RES(nproc))
      allocate(rbuf_int(2))
      IF (myrank == 0) THEN
        ListMNP(1)=MNP
        ListNP_RES(1)=NP_RES
        DO iProc=2,nproc
          CALL MPI_RECV(rbuf_int,2,itype, iProc-1, 257, comm, istatus, ierr)
          ListMNP(iProc)=rbuf_int(1)
          ListNP_RES(iProc)=rbuf_int(2)
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(ListMNP,nproc,itype, iProc-1, 263, comm, ierr)
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(ListNP_RES,nproc,itype, iProc-1, 571, comm, ierr)
        END DO
      ELSE
        rbuf_int(1)=MNP
        rbuf_int(2)=NP_RES
        CALL MPI_SEND(rbuf_int,2,itype, 0, 257, comm, ierr)
        CALL MPI_RECV(ListMNP,nproc,itype, 0, 263, comm, istatus, ierr)
        CALL MPI_RECV(ListNP_RES,nproc,itype, 0, 571, comm, istatus, ierr)
      END IF
      deallocate(rbuf_int)
      sumMNP=sum(ListMNP)
# ifdef DEBUG
      DO iProc=1,nproc
        WRITE(myrank+740,*) 'iProc=', iProc, ' np_res=', ListNP_RES(iProc)
      END DO
      WRITE(myrank+740,*) 'max(np_res)=', maxval(ListNP_RES)
      WRITE(myrank+740,*) 'nnproc=', nproc
      WRITE(myrank+740,*) 'sumMNP=', sumMNP
# endif
      allocate(ListIPLG(sumMNP))
      IF (myrank == 0) THEN
        idx=0
        DO IP=1,MNP
          idx=idx+1
          ListIPLG(idx)=iplg(IP)
        END DO
        DO iProc=2,nproc
          len=ListMNP(iProc)
          allocate(rbuf_int(len))
          CALL MPI_RECV(rbuf_int,len,itype, iProc-1, 269, comm, istatus, ierr)
          DO IP=1,len
            idx=idx+1
            ListIPLG(idx)=rbuf_int(IP)
          END DO
          deallocate(rbuf_int)
        END DO
        IF (idx /= sumMNP) THEN
          CALL WWM_ABORT('Inconsistency in IPLG creation')
        END IF
        DO iProc=2,nproc
          CALL MPI_SEND(ListIPLG,sumMNP,itype, iProc-1, 271, comm, ierr)
        END DO
      ELSE
        CALL MPI_SEND(iplg,MNP,itype, 0, 269, comm, ierr)
        CALL MPI_RECV(ListIPLG,sumMNP,itype, 0, 271, comm, istatus, ierr)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COLLECT_ALL_COVLOWER(LocalColor)
      USE DATAPOOL
      USE elfe_msgp
      USE elfe_glbl
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      integer, allocatable :: rbuf_int(:)
      integer len, iProc, IP, idx, sumMNP
      sumMNP=sum(ListMNP)
      allocate(LocalColor % ListCovLower(sumMNP))
      IF (myrank == 0) THEN
        idx=0
        DO IP=1,MNP
          idx=idx+1
          LocalColor % ListCovLower(idx)=LocalColor % CovLower(IP)
        END DO
        DO iProc=2,nproc
          len=ListMNP(iProc)
          allocate(rbuf_int(len))
          CALL MPI_RECV(rbuf_int,len,itype, iProc-1, 809, comm, istatus, ierr)
          DO IP=1,len
            idx=idx+1
            LocalColor % ListCovLower(idx)=rbuf_int(IP)
          END DO
          deallocate(rbuf_int)
        END DO
        IF (idx /= sumMNP) THEN
          CALL WWM_ABORT('Inconsistency in CovLower creation')
        END IF
        DO iProc=2,nproc
          CALL MPI_SEND(LocalColor % ListCovLower,sumMNP,itype, iProc-1, 811, comm, ierr)
        END DO
      ELSE
        CALL MPI_SEND(LocalColor % CovLower,MNP,itype, 0, 809, comm, ierr)
        CALL MPI_RECV(LocalColor % ListCovLower,sumMNP,itype, 0, 811, comm, istatus, ierr)
      END IF
# ifdef DEBUG
      WRITE(740+myrank,*) 'COLLECT_ALL_COVLOWER'
      WRITE(740+myrank,*) 'sum(ListCovLower)=', sum(LocalColor % ListCovLower)
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WWM_SOLVER_INIT
      USE DATAPOOL
      implicit none
      NblockFreqDir = NB_BLOCK 
      CALL SYMM_INIT_COLORING(MainLocalColor, NblockFreqDir)
      CALL WWM_SOLVER_ALLOCATE(SolDat)
      IF (PCmethod .eq. 2) THEN
!        CALL CREATE_ASPAR_EXCHANGE_ARRAY(LocalColor)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COLLECT_ALL_IA_JA
      USE DATAPOOL
      USE elfe_msgp
      USE elfe_glbl
      implicit none
      integer, allocatable :: rbuf_int(:)
      integer len, iProc, IP, idx, sumMNP
      integer sumIAsiz, sumNNZ
      allocate(ListNNZ(nproc))
      !
      ! Collecting NNZ
      !
      IF (myrank == 0) THEN
        ListNNZ(1)=NNZ
        allocate(rbuf_int(1))
        DO iProc=2,nproc
          CALL MPI_RECV(rbuf_int,1,itype, iProc-1, 257, comm, istatus, ierr)
          ListNNZ(iProc)=rbuf_int(1)
        END DO
        deallocate(rbuf_int)
        DO iProc=2,nproc
          CALL MPI_SEND(ListNNZ,nproc,itype, iProc-1, 263, comm, ierr)
        END DO
      ELSE
        allocate(rbuf_int(1))
        rbuf_int(1)=NNZ
        CALL MPI_SEND(rbuf_int,1,itype, 0, 257, comm, ierr)
        deallocate(rbuf_int)
        CALL MPI_RECV(ListNNZ,nproc,itype, 0, 263, comm, istatus, ierr)
      END IF
      !
      ! Collecting IA
      !
      sumIAsiz=sum(ListMNP) + nproc
# ifdef DEBUG
      DO iProc=1,nproc
        WRITE(myrank+740,*) 'iProc=', iProc, 'nnz=', ListNNZ(iProc)
      END DO
      WRITE(myrank+740,*) 'sumIAsiz=', sumIAsiz
# endif
      allocate(ListIA(sumIAsiz))
      IF (myrank == 0) THEN
        idx=0
        DO IP=1,MNP+1
          idx=idx+1
          ListIA(idx)=IA(IP)
        END DO
        DO iProc=2,nproc
          len=ListMNP(iProc)+1
          allocate(rbuf_int(len))
          CALL MPI_RECV(rbuf_int,len,itype, iProc-1, 269, comm, istatus, ierr)
          DO IP=1,len
            idx=idx+1
            ListIA(idx)=rbuf_int(IP)
          END DO
          deallocate(rbuf_int)
        END DO
        IF (idx /= sumIAsiz) THEN
          CALL WWM_ABORT('Inconsistency in sumIAsiz')
        END IF
        DO iProc=2,nproc
          CALL MPI_SEND(ListIA,sumIAsiz,itype, iProc-1, 271, comm, ierr)
        END DO
      ELSE
        CALL MPI_SEND(IA,MNP+1,itype, 0, 269, comm, ierr)
        CALL MPI_RECV(ListIA,sumIAsiz,itype, 0, 271, comm, istatus, ierr)
      END IF
      !
      ! Collecting JA
      !
      sumNNZ=sum(ListNNZ)
      allocate(ListJA(sumNNZ))
      IF (myrank == 0) THEN
        idx=0
        DO IP=1,NNZ
          idx=idx+1
          ListJA(idx)=JA(IP)
        END DO
        DO iProc=2,nproc
          len=ListNNZ(iProc)
          allocate(rbuf_int(len))
          CALL MPI_RECV(rbuf_int,len,itype, iProc-1, 569, comm, istatus, ierr)
          DO IP=1,len
            idx=idx+1
            ListJA(idx)=rbuf_int(IP)
          END DO
          deallocate(rbuf_int)
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(ListJA,sumNNZ,itype, iProc-1, 467, comm, ierr)
        END DO
      ELSE
        CALL MPI_SEND(JA,NNZ,itype, 0, 569, comm, ierr)
        CALL MPI_RECV(ListJA,sumNNZ,itype, 0, 467, comm, istatus, ierr)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
# ifdef DEBUG
      SUBROUTINE TEST_ASPAR_SYMMETRY
      USE DATAPOOL, only : MNP, IA, JA
      USE elfe_msgp, only : myrank
      implicit none
      integer :: IsSymm
      integer IP, JP, IPB, J, J2, JFOUND
      IsSymm=1
      DO IP=1,MNP
        DO J=IA(IP),IA(IP+1)-1
          JP=JA(J)
          JFOUND=-1
          DO J2=IA(JP),IA(JP+1)-1
            IPb=JA(J2)
            IF (IPb == IP) THEN
              JFOUND=J2
            END IF
          END DO
          IF (JFOUND == -1) THEN
            IsSymm=0
          END IF
        END DO
      END DO
      IF (IsSymm == 0) THEN
        WRITE(740+myrank,*) 'NNZ_IA_JA is NOT symmetric'
      ELSE
        WRITE(740+myrank,*) 'NNZ_IA_JA IS symmetric'
      END IF
      END SUBROUTINE
# endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CREATE_WWM_P2D_EXCH
      USE DATAPOOL
      USE elfe_msgp
      USE elfe_glbl
      implicit none
      include 'mpif.h'
      integer :: ListFirst(nproc)
      integer :: ListNeigh01(nproc)
      integer :: ListCommon_recv(nproc)
      integer :: ListCommon_send(nproc)
      integer :: ListMapped(np_global)
      integer :: ListMappedB(np_global)
      integer, allocatable :: dspl_send(:), dspl_recv(:)
      integer, allocatable :: dspl_send_tot(:), dspl_recv_tot(:)
      integer IP, IP_glob, iProc, WeMatch, MNPloc, idx, NP_RESloc
      integer iNeigh, IPmap, eSize, eSizeRed, nbCommon
      integer nbCommon_send, nbCommon_recv, idx_send, idx_recv
      integer sumNbCommon_send, sumNbCommon_recv
      integer idxDspl_send, idxDspl_recv
      integer eExtent, eExtentRed, NewExtent, eLB, sizRType
      integer eType1, eType2
      ListFirst=0
      DO iProc=2,nproc
        ListFirst(iProc)=ListFirst(iProc-1) + ListMNP(iProc-1)
      END DO
      ListMapped=0
      DO IP=1,MNP
        IP_glob=iplg(IP)
        IF (ListMapped(IP_glob) .gt. 0) THEN
          CALL WWM_ABORT('Clear error in ListMapped');
        ENDIF
        ListMapped(IP_glob)=IP
      END DO
      wwm_nnbr_send=0
      wwm_nnbr_recv=0
      ListCommon_send=0
      ListCommon_recv=0
      DO iProc=1,nproc
        IF (iPROC .ne. myrank+1) THEN
          MNPloc=ListMNP(iProc)
          NP_RESloc=ListNP_RES(iProc)
          ListMappedB=0
          DO IP=1,MNPloc
            IP_glob=ListIPLG(IP+ListFirst(iProc))
            IF (ListMappedB(IP_glob) .gt. 0) THEN
              CALL WWM_ABORT('Clear error in ListMappedB');
            ENDIF
            ListMappedB(IP_glob)=IP
          END DO
          !
          nbCommon_recv=0
          DO IP=1,NP_RESloc
            IP_glob=ListIPLG(IP+ListFirst(iProc))
            IF (ListMapped(IP_glob).gt.0) THEN
              nbCommon_recv=nbCommon_recv+1
            END IF
          END DO
          IF (nbCommon_recv .gt. 0) THEN
            wwm_nnbr_recv=wwm_nnbr_recv+1
            ListCommon_recv(iProc)=nbCommon_recv
          END IF
          !
          nbCommon_send=0
          DO IP=1,NP_RES
            IP_glob=iplg(IP)
            IF (ListMappedB(IP_glob).gt.0) THEN
              nbCommon_send=nbCommon_send+1
            END IF
          END DO
          IF (nbCommon_send .gt. 0) THEN
            wwm_nnbr_send=wwm_nnbr_send+1
            ListCommon_send(iProc)=nbCommon_send
          END IF
        END IF
      END DO
# ifdef DEBUG
      WRITE(740+myrank,*) 'WWM_P2D: wwm_nnbr_send=', wwm_nnbr_send
      WRITE(740+myrank,*) 'WWM_P2D: wwm_nnbr_recv=', wwm_nnbr_recv
# endif
      allocate(wwm_ListNbCommon_send(wwm_nnbr_send))
      allocate(wwm_ListNbCommon_recv(wwm_nnbr_recv))
      allocate(wwm_ListNeigh_send(wwm_nnbr_send))
      allocate(wwm_ListNeigh_recv(wwm_nnbr_recv))
      idx_send=0
      idx_recv=0
      sumNbCommon_send=0
      sumNbCommon_recv=0
      DO iProc=1,nproc
        IF (ListCommon_send(iProc) .gt. 0) THEN
          idx_send=idx_send+1
          wwm_ListNeigh_send(idx_send)=iProc
          nbCommon=ListCommon_send(iProc)
          wwm_ListNbCommon_send(idx_send)=nbCommon
          sumNbCommon_send=sumNbCommon_send+nbCommon
        END IF
        IF (ListCommon_recv(iProc) .gt. 0) THEN
          idx_recv=idx_recv+1
          wwm_ListNeigh_recv(idx_recv)=iProc
          nbCommon=ListCommon_recv(iProc)
          wwm_ListNbCommon_recv(idx_recv)=nbCommon
          sumNbCommon_recv=sumNbCommon_recv+nbCommon
        END IF
      END DO
      allocate(wwm_ListDspl_send(sumNbCommon_send))
      allocate(wwm_ListDspl_recv(sumNbCommon_recv))
      !
      ! Now the symmetric exchanges for color computations
      ! 
      wwm_nnbr=0
      DO iProc=1,nproc
        IF ((ListCommon_send(iProc).gt.0).or.(ListCommon_recv(iProc).gt.0)) THEN
          wwm_nnbr=wwm_nnbr+1
        END IF
      END DO
      allocate(wwm_ListNeigh(wwm_nnbr))
      idx=0
      DO iProc=1,nproc
        IF ((ListCommon_send(iProc).gt.0).or.(ListCommon_recv(iProc).gt.0)) THEN
          idx=idx+1
          wwm_ListNeigh(idx)=iProc
        END IF
      END DO
# ifdef DEBUG
      WRITE(740+myrank,*) 'WWM_P2D: wwm_nnbr=', wwm_nnbr
      WRITE(740+myrank,*) 'WWM_P2D: wwm_ListNeigh built'
# endif
      allocate(wwm_p2dsend_rqst(wwm_nnbr_send))
      allocate(wwm_p2drecv_rqst(wwm_nnbr_recv))
      allocate(wwm_p2dsend_stat(MPI_STATUS_SIZE,wwm_nnbr_send))
      allocate(wwm_p2drecv_stat(MPI_STATUS_SIZE,wwm_nnbr_recv))
      allocate(wwm_p2dsend_type(wwm_nnbr_send))
      allocate(wwm_p2drecv_type(wwm_nnbr_recv))
      allocate(wwmtot_p2dsend_type(wwm_nnbr_send))
      allocate(wwmtot_p2drecv_type(wwm_nnbr_recv))
# ifdef DEBUG
      WRITE(740+myrank,*) 'WWM_P2D: alloc done'
# endif
      idxDspl_send=0
      DO iNeigh=1,wwm_nnbr_send
        iProc=wwm_ListNeigh_send(iNeigh)
        MNPloc=ListMNP(iProc)
        NP_RESloc=ListNP_RES(iProc)
        ListMappedB=0
        DO IP=1,MNPloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          IF (ListMappedB(IP_glob) .gt. 0) THEN
            CALL WWM_ABORT('Clear error in ListMappedB');
          ENDIF
          ListMappedB(IP_glob)=IP
        END DO
        nbCommon=wwm_ListNbCommon_send(iNeigh)
        allocate(dspl_send(nbCommon))
        allocate(dspl_send_tot(nbCommon))
        idx=0
        DO IP=1,NP_RES
          IP_glob=iplg(IP)
          IF (ListMappedB(IP_glob).gt.0) THEN
            IPmap=ListMappedB(IP_glob)
            idx=idx+1
            dspl_send(idx)=IP-1
            dspl_send_tot(idx)=MSC*(IP-1)
            idxDspl_send=idxDspl_send+1
            wwm_ListDspl_send(idxDspl_send)=IP
          END IF
        END DO
        call mpi_type_create_indexed_block(nbCommon,1,dspl_send,rtype,wwm_p2dsend_type(iNeigh), ierr)
        call mpi_type_commit(wwm_p2dsend_type(iNeigh), ierr)
        call mpi_type_create_indexed_block(nbCommon,MSC,dspl_send_tot,rtype,wwmtot_p2dsend_type(iNeigh), ierr)
        call mpi_type_commit(wwmtot_p2dsend_type(iNeigh), ierr)
        deallocate(dspl_send)
        deallocate(dspl_send_tot)
      END DO
      idxDspl_recv=0
      DO iNeigh=1,wwm_nnbr_recv
        iProc=wwm_ListNeigh_recv(iNeigh)
        MNPloc=ListMNP(iProc)
        NP_RESloc=ListNP_RES(iProc)
        nbCommon=wwm_ListNbCommon_recv(iNeigh)
        allocate(dspl_recv(nbCommon))
        allocate(dspl_recv_tot(nbCommon))
        idx=0
        DO IP=1,NP_RESloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          IF (ListMapped(IP_glob).gt.0) THEN
            IPmap=ListMapped(IP_glob)
            idx=idx+1
            dspl_recv(idx)=IPmap-1
            dspl_recv_tot(idx)=MSC*(IPmap-1)
            idxDspl_recv=idxDspl_recv+1
            wwm_ListDspl_recv(idxDspl_recv)=IPmap
          END IF
        END DO
        !
        call mpi_type_create_indexed_block(nbCommon,1,dspl_recv,rtype,wwm_p2drecv_type(iNeigh),ierr)
        call mpi_type_commit(wwm_p2drecv_type(iNeigh), ierr)
        call mpi_type_create_indexed_block(nbCommon,MSC,dspl_recv_tot,rtype,wwmtot_p2drecv_type(iNeigh),ierr)
        call mpi_type_commit(wwmtot_p2drecv_type(iNeigh), ierr)
        deallocate(dspl_recv)
        deallocate(dspl_recv_tot)
      END DO
# ifdef DEBUG
      WRITE(740+myrank,*) 'Leaving CREATE_WWM_P2D_EXCH'
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CREATE_WWM_MAT_P2D_EXCH
      USE DATAPOOL
      USE elfe_msgp
      USE elfe_glbl
      implicit none
      include 'mpif.h'
      integer :: ListFirstMNP(nproc), ListFirstNNZ(nproc)
      integer :: ListCommon_send(nproc), ListCommon_recv(nproc)
      integer :: ListMapped(np_global)
      integer :: ListMappedB(np_global)
      integer, allocatable :: dspl_send(:), dspl_recv(:)
      integer nbCommon_send, nbCommon_recv
      integer IAfirst
      integer IP, JP, I, J, J2, IP_glob, JP_glob, iProc
      integer WeMatch, MNPloc, NP_RESloc, JP_j
      integer IPloc, JPloc, Jfound, idx
      integer iNeigh, IPmap, nbCommon, nbCommonB, eSize, eSizeRed
      integer sumNbCommon_send, sumNbCommon_recv
      integer idxDspl_send, idxDspl_recv
      ListFirstNNZ=0
      ListFirstMNP=0
      DO iProc=2,nproc
        ListFirstMNP(iProc)=ListFirstMNP(iProc-1) + ListMNP(iProc-1)
        ListFirstNNZ(iProc)=ListFirstNNZ(iProc-1) + ListNNZ(iProc-1)
      END DO
      ListMapped=0
      DO IP=1,MNP
        IP_glob=iplg(IP)
        IF (ListMapped(IP_glob) .gt. 0) THEN
          CALL WWM_ABORT('Clear error in ListMapped');
        ENDIF
        ListMapped(IP_glob)=IP
      END DO
      wwm_nnbr_m_recv=0
      wwm_nnbr_m_send=0
      ListCommon_recv=0
      DO I=1,wwm_nnbr_recv
        iProc=wwm_ListNeigh_recv(I)
        NP_RESloc=ListNP_RES(iProc)
        MNPloc=ListMNP(iProc)
        ListMappedB=0
        DO IP=1,MNPloc
          IP_glob=ListIPLG(IP+ListFirstMNP(iProc))
          IF (ListMappedB(IP_glob) .gt. 0) THEN
            CALL WWM_ABORT('Clear error in ListMappedB');
          ENDIF
          ListMappedB(IP_glob)=IP
        END DO
        nbCommon_recv=0
        DO IP=1,NP_RESloc
          IP_glob=ListIPLG(IP+ListFirstMNP(iProc))
          IPloc=ListMapped(IP_glob)
          IF (IPloc.gt.0) THEN
            IAfirst=ListFirstMNP(iProc) + iProc-1
            DO J=ListIA(IP+IAfirst),ListIA(IP+IAfirst+1)-1
              JP=ListJA(J+ListFirstNNZ(iProc))
              JP_glob=ListIPLG(JP+ListFirstMNP(iProc))
              JPloc=ListMapped(JP_glob)
              IF (JPloc.gt.0) THEN
                JFOUND=-1
                DO J2=IA(IPloc),IA(IPloc+1)-1
                  IF (JA(J2) == JPloc) THEN
                    JFOUND=J2
                  END IF
                END DO
                IF (JFOUND .gt. 0) THEN
                  nbCommon_recv=nbCommon_recv+1
                END IF
              END IF
            END DO
          END IF
        END DO
        IF (nbCommon_recv .gt. 0) THEN
          wwm_nnbr_m_recv=wwm_nnbr_m_recv+1
          ListCommon_recv(iProc)=nbCommon_recv
        END IF
      END DO
      ListCommon_send=0
      DO I=1,wwm_nnbr_send
        iProc=wwm_ListNeigh_send(I)
        NP_RESloc=ListNP_RES(iProc)
        MNPloc=ListMNP(iProc)
        ListMappedB=0
        DO IP=1,MNPloc
          IP_glob=ListIPLG(IP+ListFirstMNP(iProc))
          IF (ListMappedB(IP_glob) .gt. 0) THEN
            CALL WWM_ABORT('Clear error in ListMappedB');
          ENDIF
          ListMappedB(IP_glob)=IP
        END DO
        nbCommon_send=0
        DO IP=1,NP_RES
          IP_glob=iplg(IP)
          IPloc=ListMappedB(IP_glob)
          IF (IPloc.gt.0) THEN
            DO J=IA(IP),IA(IP+1)-1
              JP=JA(J)
              JP_glob=iplg(JP)
              JPloc=ListMappedB(JP_glob)
              IF (JPloc.gt.0) THEN
                IAfirst=ListFirstMNP(iProc) + iProc-1
                JFOUND=-1
                DO J2=ListIA(IPloc+IAfirst),ListIA(IPloc+IAfirst+1)-1
                  JP_j=ListJA(J2+ListFirstNNZ(iProc))
                  IF (JP_j == JPloc) THEN
                    JFOUND=J2
                  END IF
                END DO
                IF (JFOUND .gt. 0) THEN
                  nbCommon_send=nbCommon_send+1
                END IF
              END IF
            END DO
          END IF
        END DO
        IF (nbCommon_send .gt. 0) THEN
          wwm_nnbr_m_send=wwm_nnbr_m_send+1
          ListCommon_send(iProc)=nbCommon_send
        END IF
      END DO
# ifdef DEBUG
      WRITE(740+myrank,*) 'WWM_MAT_P2D: wwm_nnbr_m_send=', wwm_nnbr_m_send
      WRITE(740+myrank,*) 'WWM_MAT_P2D: wwm_nnbr_m_recv=', wwm_nnbr_m_recv
# endif
      allocate(wwmmat_p2dsend_rqst(wwm_nnbr_m_send))
      allocate(wwmmat_p2drecv_rqst(wwm_nnbr_m_recv))
      allocate(wwmmat_p2dsend_stat(MPI_STATUS_SIZE,wwm_nnbr_m_send))
      allocate(wwmmat_p2drecv_stat(MPI_STATUS_SIZE,wwm_nnbr_m_recv))
      allocate(wwmmat_p2dsend_type(wwm_nnbr_m_send))
      allocate(wwmmat_p2drecv_type(wwm_nnbr_m_recv))
      allocate(wwm_ListNbCommon_m_send(wwm_nnbr_m_send))
      allocate(wwm_ListNbCommon_m_recv(wwm_nnbr_m_recv))
      allocate(wwm_ListNeigh_m_recv(wwm_nnbr_m_recv))
      allocate(wwm_ListNeigh_m_send(wwm_nnbr_m_send))
# ifdef DEBUG
      WRITE(740+myrank,*) 'WWM_MAT_P2D: alloc done'
# endif
      idx=0
      sumNbCommon_send=0
      DO iProc=1,nproc
        IF (ListCommon_send(iProc) .gt. 0) THEN
          idx=idx+1
          wwm_ListNeigh_m_send(idx)=iProc
          nbCommon=ListCommon_send(iProc)
          wwm_ListNbCommon_m_send(idx)=nbCommon
          sumNbCommon_send=sumNbCommon_send+nbCommon
        END IF
      END DO
# ifdef DEBUG
      WRITE(740+myrank,*) 'WWM_MAT_P2D: wwm_ListNeigh_m_send built'
# endif
      idx=0
      sumNbCommon_recv=0
      DO iProc=1,nproc
        IF (ListCommon_recv(iProc) .gt. 0) THEN
          idx=idx+1
          wwm_ListNeigh_m_recv(idx)=iProc
          nbCommon=ListCommon_recv(iProc)
          wwm_ListNbCommon_m_recv(idx)=nbCommon
          sumNbCommon_recv=sumNbCommon_recv+nbCommon
        END IF
      END DO
# ifdef DEBUG
      WRITE(740+myrank,*) 'WWM_MAT_P2D: wwm_ListNeigh_m_recv built'
# endif
      allocate(wwm_ListDspl_m_send(sumNbCommon_send))
      idxDspl_send=0
      DO iNeigh=1,wwm_nnbr_m_send
        iProc=wwm_ListNeigh_m_send(iNeigh)
        MNPloc=ListMNP(iProc)
        ListMappedB=0
        DO IP=1,MNPloc
          IP_glob=ListIPLG(IP+ListFirstMNP(iProc))
# ifdef DEBUG
          IF (ListMappedB(IP_glob) .gt. 0) THEN
            CALL WWM_ABORT('Clear error in ListMappedB');
          ENDIF
# endif
          ListMappedB(IP_glob)=IP
        END DO
        nbCommon=wwm_ListNbCommon_m_send(iNeigh)
        allocate(dspl_send(nbCommon))
        idx=0
        DO IP=1,NP_RES
          IP_glob=iplg(IP)
          IPloc=ListMappedB(IP_glob)
          IF (IPloc.gt.0) THEN
            DO J=IA(IP),IA(IP+1)-1
              JP=JA(J)
              JP_glob=iplg(JP)
              JPloc=ListMappedB(JP_glob)
              IF (JPloc .gt. 0) THEN
                IAfirst=ListFirstMNP(iProc) + iProc-1
                JFOUND=-1
                DO J2=ListIA(IPloc+IAfirst),ListIA(IPloc+IAfirst+1)-1
                  JP_j=ListJA(J2+ListFirstNNZ(iProc))
                  IF (JP_j == JPloc) THEN
                    JFOUND=J2
                  END IF
                END DO
                IF (JFOUND .gt. 0) THEN
                  idxDspl_send=idxDspl_send+1
                  wwm_ListDspl_m_send(idxDspl_send)=J
                  idx=idx+1
                  dspl_send(idx)=J-1
                END IF
              END IF
            END DO
          END IF
        END DO
        call mpi_type_create_indexed_block(nbCommon,1,dspl_send,rtype,wwmmat_p2dsend_type(iNeigh), ierr)
        call mpi_type_commit(wwmmat_p2dsend_type(iNeigh), ierr)
        deallocate(dspl_send)
      END DO
# ifdef DEBUG
      WRITE(740+myrank,*) 'sumNbCommon_recv=', sumNbCommon_recv
# endif
      allocate(wwm_ListDspl_m_recv(sumNbCommon_recv))
      idxDspl_recv=0
      DO iNeigh=1,wwm_nnbr_m_recv
        iProc=wwm_ListNeigh_m_recv(iNeigh)
        nbCommon=wwm_ListNbCommon_m_recv(iNeigh)
# ifdef DEBUG
        WRITE(740+myrank,*) 'nbCommon=', nbCommon
# endif
        allocate(dspl_recv(nbCommon))
        NP_RESloc=ListNP_RES(iProc)
        idx=0
        DO IP=1,NP_RESloc
          IP_glob=ListIPLG(IP+ListFirstMNP(iProc))
          IPloc=ListMapped(IP_glob)
          IF (IPloc.gt.0) THEN
            IAfirst=ListFirstMNP(iProc) + iProc-1
            DO J=ListIA(IP+IAfirst),ListIA(IP+IAfirst+1)-1
              JP=ListJA(J+ListFirstNNZ(iProc))
              JP_glob=ListIPLG(JP+ListFirstMNP(iProc))
              JPloc=ListMapped(JP_glob)
              IF (JPloc.gt.0) THEN
                JFOUND=-1
                DO J2=IA(IPloc),IA(IPloc+1)-1
                  IF (JA(J2) == JPloc) THEN
                    JFOUND=J2
                  END IF
                END DO
                IF (JFOUND /= -1) THEN
                  idxDspl_recv=idxDspl_recv+1
                  wwm_ListDspl_m_recv(idxDspl_recv)=Jfound
                  idx=idx+1
                  dspl_recv(idx)=Jfound-1
                END IF
              END IF
            END DO
          END IF
        END DO
        call mpi_type_create_indexed_block(nbCommon,1,dspl_recv,rtype,wwmmat_p2drecv_type(iNeigh), ierr)
        call mpi_type_commit(wwmmat_p2drecv_type(iNeigh), ierr)
        deallocate(dspl_recv)
      END DO
# ifdef DEBUG
      WRITE(740+myrank,*) 'Leaving CREATE_WWM_MAT_P2D_EXCH'
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
# ifdef DEBUG
      SUBROUTINE CHECK_STANDARD_SELFE_EXCH
      USE elfe_msgp, only : exchange_p2d, myrank
      USE DATAPOOL
      REAL(rkind) :: XPcopy(MNP)
      REAL(rkind) :: SumErr
      XPcopy=XP
      CALL exchange_p2d(XP)
      SumErr=sum(abs(XPcopy - XP))
!      WRITE(740+myrank,*) 'SELFE_EXCH SumErr=', SumErr
      END SUBROUTINE
# endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
# ifdef DEBUG
      SUBROUTINE CHECK_EXCHANGE_P4D_WWM
      USE DATAPOOL, only : rkind, MDC, MSC, MNP, XP, YP
      USE elfe_msgp, only : myrank, EXCHANGE_P4D_WWM
      USE elfe_glbl, only : iplg
      implicit none
      integer :: IS, ID, idx
      REAL(rkind) :: eField(MSC, MDC, MNP)
      REAL(rkind) :: eFieldR(MSC, MDC, MNP)
      REAL(rkind) :: MaxErr
      idx=0
      DO IS=1,MSC
        DO ID=1,MDC
          idx=idx+1
          IF (idx == 4) THEN
            idx=1
          END IF
          IF (idx == 1) THEN
            eField(IS,ID,:)=XP
          END IF
          IF (idx == 2) THEN
            eField(IS,ID,:)=YP
          END IF
          IF (idx == 3) THEN
            eField(IS,ID,:)=MyREAL(iplg)
          END IF
        END DO
      END DO
      eFieldR=eField
      CALL EXCHANGE_P4D_WWM(eField)
      MaxErr=maxval(abs(eFieldR - eField))
      WRITE(740+myrank,*) 'Test EXCHANGE_P4D_WWM MaxErr=', MaxErr
      END SUBROUTINE
# endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
# ifdef DEBUG
      SUBROUTINE CHECK_P2D_EXCH_BLOCK(LocalColor, eField, siz)
      USE elfe_msgp
      USE DATAPOOL
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      integer, intent(in) :: siz
      real(rkind), intent(in) :: eField(MNP, siz)
      integer iLow, iUpp, nbLow_send, nbUpp_send, nbLow_recv, nbUpp_recv
      integer I, iRank, IP, iSiz
      real(rkind) :: MaxErr, eDiff
      real(rkind), allocatable :: TheRecv(:,:), p2d_data_send(:,:), eFieldC(:,:)
      integer nbDiff
      IF (siz .gt. LocalColor % maxBlockLength) THEN
        Print *, 'siz=', siz
        Print *, 'maxBlockLength=', LocalColor % maxBlockLength
        CALL WWM_ABORT('Error between siz and maxBlockLength')
      ENDIF
      allocate(eFieldC(LocalColor % maxBlockLength, MNP))
      allocate(TheRecv(LocalColor % maxBlockLength, MNP))
      allocate(p2d_data_send(LocalColor % maxBlockLength, MNP))
      eFieldC=0
      DO I=1,siz
        eFieldC(I,:)=eField(:,I)
      END DO
      nbLow_recv=LocalColor%nbLow_recv
      nbUpp_send=LocalColor%nbUpp_send
      p2d_data_send=eFieldC
      TheRecv=eFieldC
      DO iUpp=1,nbUpp_send
        I=LocalColor % ListIdxUpper_send(iUpp)
        iRank=wwm_ListNeigh_send(i)
        call mpi_isend(p2d_data_send,1,LocalColor % blk_p2dsend_type(i),iRank-1,13,comm,LocalColor % Upp_s_rq(iUpp),ierr)
      END DO
      DO iLow=1,nbLow_recv
        I=LocalColor % ListIdxLower_recv(iLow)
        iRank=wwm_ListNeigh_recv(i)
        call mpi_irecv(TheRecv,1,LocalColor % blk_p2drecv_type(i),iRank-1,13,comm,LocalColor % Low_r_rq(iLow),ierr)
      END DO
      IF (nbUpp_send > 0) THEN
        call mpi_waitall(nbUpp_send, LocalColor%Upp_s_rq, LocalColor%Upp_s_stat,ierr)
      END IF
      IF (nbLow_recv > 0) THEN
        call mpi_waitall(nbLow_recv, LocalColor%Low_r_rq, LocalColor%Low_r_stat,ierr)
      END IF
      MaxErr=maxval(abs(eFieldC - TheRecv))
      WRITE(740+myrank,*) 'CHECK_P2D_EXCH_BLOCK 1: MaxErr=', MaxErr
      !
      !
      nbLow_send=LocalColor%nbLow_send
      nbUpp_recv=LocalColor%nbUpp_recv
      DO iLow=1,nbLow_send
        I=LocalColor % ListIdxLower_send(iLow)
        iRank=wwm_ListNeigh_send(i)
        call mpi_isend(p2d_data_send,1,LocalColor % blk_p2dsend_type(i),iRank-1,17,comm,LocalColor%Low_s_rq(iLow),ierr)
      END DO
      DO iUpp=1,nbUpp_recv
        I=LocalColor % ListIdxUpper_recv(iUpp)
        iRank=wwm_ListNeigh_recv(i)
        call mpi_irecv(TheRecv,1,LocalColor % blk_p2drecv_type(i),iRank-1,17,comm,LocalColor%Upp_r_rq(iUpp),ierr)
      END DO
      IF (nbLow_send > 0) THEN
        call mpi_waitall(nbLow_send, LocalColor%Low_s_rq, LocalColor%Low_s_stat,ierr)
      END IF
      IF (nbUpp_recv > 0) THEN
        call mpi_waitall(nbUpp_recv, LocalColor%Upp_r_rq, LocalColor%Upp_r_stat,ierr)
      END IF
      MaxErr=maxval(abs(eFieldC - TheRecv))
      WRITE(740+myrank,*) 'CHECK_P2D_EXCH_BLOCK 2: MaxErr=', MaxErr
      DO iSiz=1,siz
        MaxErr=maxval(abs(eFieldC(iSiz,:) - TheRecv(iSiz,:)))
        if (MaxErr .gt. 0) THEN
          WRITE(740+myrank,*) 'iSiz=', iSiz, '/', siz
          nbDiff=0
          DO IP=1,MNP
            eDiff=abs(eFieldC(iSiz,IP) - TheRecv(iSiz,IP))
            IF (eDiff .gt. 0) THEN
              WRITE(740+myrank,*) IP, LocalColor%CovLower(IP), eFieldC(iSiz,IP), TheRecv(iSiz,IP)
              nbDiff=nbDiff+1
            END IF
          END DO
          WRITE(740+myrank,*) 'nbDiff=', nbDiff
        END IF
      END DO
      END SUBROUTINE
# endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
# ifdef DEBUG
      SUBROUTINE CHECK_P2D_EXCH_BLOCK_REV(LocalColor, eField, siz)
      USE elfe_msgp
      USE DATAPOOL
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      integer, intent(in) :: siz
      real(rkind), intent(in) :: eField(MNP, siz)
      integer iLow, iUpp, nbLow_send, nbUpp_send, nbLow_recv, nbUpp_recv
      integer I, iRank, IP, iSiz
      real(rkind) :: MaxErr, eDiff
      real(rkind), allocatable :: TheRecv(:,:), p2d_data_send(:,:), eFieldC(:,:)
      integer nbDiff
      IF (siz .gt. LocalColor % maxBlockLength) THEN
        Print *, 'siz=', siz
        Print *, 'maxBlockLength=', LocalColor % maxBlockLength
        CALL WWM_ABORT('Error between siz and maxBlockLength')
      ENDIF
      allocate(eFieldC(MNP, LocalColor % maxBlockLength))
      allocate(TheRecv(MNP, LocalColor % maxBlockLength))
      allocate(p2d_data_send(MNP, LocalColor % maxBlockLength))
      eFieldC=0
      DO I=1,siz
        eFieldC(:,I)=eField(:,I)
      END DO
      nbLow_recv=LocalColor%nbLow_recv
      nbUpp_send=LocalColor%nbUpp_send
      p2d_data_send=eFieldC
      TheRecv=eFieldC
      DO iUpp=1,nbUpp_send
        I=LocalColor % ListIdxUpper_send(iUpp)
        iRank=wwm_ListNeigh_send(i)
        call mpi_isend(p2d_data_send,1,LocalColor % blk_p2dsend_type(i),iRank-1,13,comm,LocalColor % Upp_s_rq(iUpp),ierr)
      END DO
      DO iLow=1,nbLow_recv
        I=LocalColor % ListIdxLower_recv(iLow)
        iRank=wwm_ListNeigh_recv(i)
        call mpi_irecv(TheRecv,1,LocalColor % blk_p2drecv_type(i),iRank-1,13,comm,LocalColor % Low_r_rq(iLow),ierr)
      END DO
      IF (nbUpp_send > 0) THEN
        call mpi_waitall(nbUpp_send, LocalColor%Upp_s_rq, LocalColor%Upp_s_stat,ierr)
      END IF
      IF (nbLow_recv > 0) THEN
        call mpi_waitall(nbLow_recv, LocalColor%Low_r_rq, LocalColor%Low_r_stat,ierr)
      END IF
      MaxErr=maxval(abs(eFieldC - TheRecv))
      WRITE(740+myrank,*) 'CHECK_P2D_EXCH_BLOCK_REV 1: MaxErr=', MaxErr
      !
      !
      nbLow_send=LocalColor%nbLow_send
      nbUpp_recv=LocalColor%nbUpp_recv
      DO iLow=1,nbLow_send
        I=LocalColor % ListIdxLower_send(iLow)
        iRank=wwm_ListNeigh_send(i)
        call mpi_isend(p2d_data_send,1,LocalColor % blk_p2dsend_type(i),iRank-1,17,comm,LocalColor%Low_s_rq(iLow),ierr)
      END DO
      DO iUpp=1,nbUpp_recv
        I=LocalColor % ListIdxUpper_recv(iUpp)
        iRank=wwm_ListNeigh_recv(i)
        call mpi_irecv(TheRecv,1,LocalColor % blk_p2drecv_type(i),iRank-1,17,comm,LocalColor%Upp_r_rq(iUpp),ierr)
      END DO
      IF (nbLow_send > 0) THEN
        call mpi_waitall(nbLow_send, LocalColor%Low_s_rq, LocalColor%Low_s_stat,ierr)
      END IF
      IF (nbUpp_recv > 0) THEN
        call mpi_waitall(nbUpp_recv, LocalColor%Upp_r_rq, LocalColor%Upp_r_stat,ierr)
      END IF
      MaxErr=maxval(abs(eFieldC - TheRecv))
      WRITE(740+myrank,*) 'CHECK_P2D_EXCH_BLOCK_REV 2: MaxErr=', MaxErr
      DO iSiz=1,siz
        MaxErr=maxval(abs(eFieldC(iSiz,:) - TheRecv(iSiz,:)))
        if (MaxErr .gt. 0) THEN
          WRITE(740+myrank,*) 'iSiz=', iSiz, '/', siz
          nbDiff=0
          DO IP=1,MNP
            eDiff=abs(eFieldC(IP,iSiz) - TheRecv(IP,iSiz))
            IF (eDiff .gt. 0) THEN
              WRITE(740+myrank,*) IP, LocalColor%CovLower(IP), eFieldC(IP,iSiz), TheRecv(IP,iSiz)
              nbDiff=nbDiff+1
            END IF
          END DO
          WRITE(740+myrank,*) 'nbDiff=', nbDiff
        END IF
      END DO
      END SUBROUTINE
# endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
# ifdef DEBUG
      SUBROUTINE CHECK_P2D_EXCH_STACKED(LocalColor, eField, siz)
      USE elfe_msgp
      USE DATAPOOL
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      integer, intent(in) :: siz
      real(rkind), intent(in) :: eField(MNP, siz)
      integer iLow, iUpp, nbLow_send, nbUpp_send, nbLow_recv, nbUpp_recv
      integer I, iRank, IP, iSiz
      real(rkind) :: TheRecv(MNP,siz), p2d_data_send(MNP,siz)
      real(rkind) :: MaxErr, eDiff
      integer nbDiff
      nbLow_recv=LocalColor%nbLow_recv
      nbUpp_send=LocalColor%nbUpp_send
      p2d_data_send=eField
      TheRecv=eField
      DO iUpp=1,nbUpp_send
        I=LocalColor % ListIdxUpper_send(iUpp)
        iRank=wwm_ListNeigh_send(i)
        call mpi_isend(p2d_data_send,siz,wwm_p2dsend_type(i),iRank-1,13,comm,LocalColor % Upp_r_rq(iUpp),ierr)
      END DO
      DO iLow=1,nbLow_recv
        I=LocalColor % ListIdxLower_recv(iLow)
        iRank=wwm_ListNeigh_recv(i)
        call mpi_irecv(TheRecv,siz,wwm_p2drecv_type(i),iRank-1,13,comm,LocalColor % Low_r_rq(iLow),ierr)
      END DO
      IF (nbUpp_send > 0) THEN
        call mpi_waitall(nbUpp_send, LocalColor%Upp_s_rq, LocalColor%Upp_s_stat,ierr)
      END IF
      IF (nbLow_recv > 0) THEN
        call mpi_waitall(nbLow_recv, LocalColor%Low_r_rq, LocalColor%Low_r_stat,ierr)
      END IF
      MaxErr=maxval(abs(eField - TheRecv))
      WRITE(740+myrank,*) 'CHECK_P2D_EXCH_STACKED 1: MaxErr=', MaxErr
      DO iSiz=1,siz
        MaxErr=maxval(abs(eField(:,iSiz) - TheRecv(:,iSiz)))
        if (MaxErr .gt. 0) THEN
          WRITE(740+myrank,*) 'iSiz=', iSiz, '/', siz
          nbDiff=0
          DO IP=1,MNP
            eDiff=abs(eField(IP,iSiz) - TheRecv(IP,iSiz))
            IF (eDiff .gt. 0) THEN
              WRITE(740+myrank,*) IP, LocalColor%CovLower(IP), eField(IP,iSiz), TheRecv(IP,iSiz), 'diff'
              nbDiff=nbDiff+1
            ELSE
              WRITE(740+myrank,*) IP, LocalColor%CovLower(IP), eField(IP,iSiz), TheRecv(IP,iSiz)
            ENDIF
          END DO
          WRITE(740+myrank,*) 'nbDiff=', nbDiff
        END IF
      END DO
      !
      !
      nbLow_send=LocalColor%nbLow_send
      nbUpp_recv=LocalColor%nbUpp_recv
      p2d_data_send=eField
      TheRecv=eField
      DO iLow=1,nbLow_send
        I=LocalColor % ListIdxLower_send(iLow)
        iRank=wwm_ListNeigh_send(i)
        call mpi_isend(p2d_data_send,siz,wwm_p2dsend_type(i),iRank-1,17,comm,LocalColor%Low_s_rq(iLow),ierr)
      END DO
      DO iUpp=1,nbUpp_recv
        I=LocalColor % ListIdxUpper_recv(iUpp)
        iRank=wwm_ListNeigh_recv(i)
        call mpi_irecv(TheRecv,siz,wwm_p2drecv_type(i),iRank-1,17,comm,LocalColor%Upp_r_rq(iUpp),ierr)
      END DO
      IF (nbLow_send > 0) THEN
        call mpi_waitall(nbLow_send, LocalColor%Low_s_rq, LocalColor%Low_s_stat,ierr)
      END IF
      IF (nbUpp_recv > 0) THEN
        call mpi_waitall(nbUpp_recv, LocalColor%Upp_r_rq, LocalColor%Upp_r_stat,ierr)
      END IF
      MaxErr=maxval(abs(eField - TheRecv))
      WRITE(740+myrank,*) 'CHECK_P2D_EXCH_STACKED 2: MaxErr=', MaxErr
      DO iSiz=1,siz
        MaxErr=maxval(abs(eField(:,iSiz) - TheRecv(:,iSiz)))
        if (MaxErr .gt. 0) THEN
          WRITE(740+myrank,*) 'iSiz=', iSiz, '/', siz
          nbDiff=0
          DO IP=1,MNP
            eDiff=abs(eField(IP,iSiz) - TheRecv(IP,iSiz))
            IF (eDiff .gt. 0) THEN
              WRITE(740+myrank,*) IP, LocalColor%CovLower(IP), eField(IP,iSiz), TheRecv(IP,iSiz)
              nbDiff=nbDiff+1
            END IF
          END DO
          WRITE(740+myrank,*) 'nbDiff=', nbDiff
        END IF
      END DO
      END SUBROUTINE
# endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
# ifdef DEBUG
      SUBROUTINE CHECK_P2D_EXCH(LocalColor, eField)
      USE elfe_msgp
      USE DATAPOOL
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      real(rkind), intent(in) :: eField(MNP)
      integer iLow, iUpp, nbLow_send, nbUpp_send, nbLow_recv, nbUpp_recv
      integer I, iRank
      real(rkind) :: TheRecv(MNP), p2d_data_send(MNP)
      real(rkind) :: MaxErr
      nbLow_recv=LocalColor%nbLow_recv
      nbUpp_send=LocalColor%nbUpp_send
      p2d_data_send=eField
      DO iUpp=1,nbUpp_send
        I=LocalColor % ListIdxUpper_send(iUpp)
        iRank=wwm_ListNeigh_send(i)
        call mpi_isend(p2d_data_send,1,wwm_p2dsend_type(i),iRank-1,13,comm,LocalColor % Upp_s_rq(iUpp),ierr)
      END DO
      TheRecv=eField
      DO iLow=1,nbLow_recv
        I=LocalColor % ListIdxLower_recv(iLow)
        iRank=wwm_ListNeigh_recv(i)
        call mpi_irecv(TheRecv,1,wwm_p2drecv_type(i),iRank-1,13,comm,LocalColor % Low_r_rq(iLow),ierr)
      END DO
      IF (nbUpp_send > 0) THEN
        call mpi_waitall(nbUpp_send, LocalColor%Upp_s_rq, LocalColor%Upp_s_stat,ierr)
      END IF
      IF (nbLow_recv > 0) THEN
        call mpi_waitall(nbLow_recv, LocalColor%Low_r_rq, LocalColor%Low_r_stat,ierr)
      END IF
      MaxErr=maxval(abs(eField - TheRecv))
      WRITE(740+myrank,*) 'CHECK_P2D_EXCH 1: MaxErr=', MaxErr
      !
      !
      nbLow_send=LocalColor%nbLow_send
      nbUpp_recv=LocalColor%nbUpp_recv
      DO iLow=1,nbLow_send
        I=LocalColor % ListIdxLower_send(iLow)
        iRank=wwm_ListNeigh_send(i)
        call mpi_isend(p2d_data_send,1,wwm_p2dsend_type(i),iRank-1,17,comm,LocalColor%Low_s_rq(iLow),ierr)
      END DO
      DO iUpp=1,nbUpp_recv
        I=LocalColor % ListIdxUpper_recv(iUpp)
        iRank=wwm_ListNeigh_recv(i)
        call mpi_irecv(TheRecv,1,wwm_p2drecv_type(i),iRank-1,17,comm,LocalColor%Upp_r_rq(iUpp),ierr)
      END DO
      IF (nbLow_send > 0) THEN
        call mpi_waitall(nbLow_send, LocalColor%Low_s_rq, LocalColor%Low_s_stat,ierr)
      END IF
      IF (nbUpp_recv > 0) THEN
        call mpi_waitall(nbUpp_recv, LocalColor%Upp_r_rq, LocalColor%Upp_r_stat,ierr)
      END IF
      MaxErr=maxval(abs(eField - TheRecv))
      WRITE(740+myrank,*) 'CHECK_P2D_EXCH 2: MaxErr=', MaxErr
      END SUBROUTINE
# endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
# ifdef DEBUG
      SUBROUTINE WRITE_EXPLICIT_ORDERING(ListPos, ListPosRev, ListColor)
      USE DATAPOOL, only : MNP, ListMNP, ListNP_RES, ListIPLG
      USE elfe_msgp, only : nproc, myrank
      USE elfe_glbl, only : np_global
      IMPLICIT NONE
      integer, intent(inout) :: ListPos(np_global)
      integer, intent(inout) :: ListPosRev(np_global)
      integer, intent(in) :: ListColor(nproc)
      integer :: ListFirst(nproc)
      integer :: ListTotal(np_global)
      integer minColor, maxColor, eColor
      integer iProc, idx, IP, IPglob
      minColor=minval(ListColor)
      maxColor=maxval(ListColor)
      ListTotal=1
      ListFirst=0
      DO iProc=2,nproc
        ListFirst(iProc)=ListFirst(iProc-1) + ListMNP(iProc-1)
      END DO
      ListPos=0
      ListPosRev=0
      idx=0
      DO eColor=minColor,maxColor
        DO iProc=1,nproc
          IF (ListColor(iProc) .eq. eColor) THEN
            DO IP=1,ListNP_RES(iProc)
              IPglob=ListIPLG(IP+ListFirst(iProc))
              IF (ListTotal(IPglob) .eq. 1) THEN
                idx=idx+1
                ListPos(idx)=IPglob
                ListPosRev(IPglob)=idx
                ListTotal(IPglob)=0
              END IF
            END DO
          END IF
        END DO
      END DO
      IF (minval(ListPosRev) .eq. 0) THEN
        CALL WWM_ABORT('Please correct ')
      END IF
      IF (idx .ne. np_global) THEN
        DO IP=1,np_global
          WRITE(myrank+540,*) 'IP/Pos/PosR=', IP, ListPos(IP), ListPosRev(IP)
        END DO
        CALL WWM_ABORT('One more bug to solve')
      END IF
      END SUBROUTINE
# endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE BUILD_MULTICOLORING(AdjGraph, ListColor)
      USE elfe_msgp, only : myrank
      implicit none
      type(Graph), intent(in) :: AdjGraph
      integer, intent(out) :: ListColor(AdjGraph%nbVert)
      integer, allocatable :: CurrColor(:)
      integer MaxDeg, iVert, eVert, eColor, eDeg
      integer idx, I, ChromaticNr, nbVert
      integer, allocatable :: ListPosFirst(:)
      integer, allocatable :: TheOrdering(:)
      integer eColorF, iVertFound, eAdjColor
      integer nbUndef, MinDeg, eAdj, MinUndef, PosMin
      MaxDeg=AdjGraph % MaxDeg
      nbVert=AdjGraph % nbVert
      allocate(CurrColor(MaxDeg+1))
      ListColor=0
      idx=0
      allocate(ListPosFirst(nbVert))
      DO iVert=1,nbVert
        ListPosFirst(iVert)=idx
        idx=idx+AdjGraph % ListDegree(iVert)
      END DO
      MinDeg=MaxDeg+3
      PosMin=-1
      DO iVert=1,nbVert
        eDeg=AdjGraph % ListDegree(iVert)
        IF (eDeg .lt. MinDeg) THEN
          MinDeg=eDeg
          PosMin=iVert
        END IF
      END DO
      idx=ListPosFirst(PosMin)
      DO I=0,MinDeg
        IF (I.eq.0) THEN
          eVert=PosMin
        ELSE
          eVert=AdjGraph % ListEdge(idx+I,2)
        END IF
        ListColor(eVert)=I+1
      END DO
      DO
        MinUndef=nbVert
        iVertFound=0
        DO iVert=1,nbVert
          IF (ListColor(iVert) == 0) THEN
            idx=ListPosFirst(iVert)
            eDeg=AdjGraph % ListDegree(iVert)
            nbUndef=0
            DO I=1,eDeg
              eAdj=AdjGraph % ListEdge(idx+I,2)
              eAdjColor=ListColor(eAdj)
              IF (eAdjColor == 0) THEN
                nbUndef=nbUndef+1
              END IF
            END DO
            IF (nbUndef .lt. MinUndef) THEN
              MinUndef=nbUndef
              iVertFound=iVert
            END IF
          END IF
        END DO
        IF (iVertFound == 0) THEN
          EXIT
        END IF
        eDeg=AdjGraph % ListDegree(iVertFound)
        idx=ListPosFirst(iVertFound)
        CurrColor=0
        DO I=1,eDeg
          eVert=AdjGraph % ListEdge(idx+I,2)
          eColor=ListColor(eVert)
          IF (eColor.gt.0) THEN
            CurrColor(eColor)=1
          END IF
        END DO
        eColorF=-1
        DO I=1,MaxDeg+1
          IF (eColorF == -1) THEN
            IF (CurrColor(I) == 0) THEN
              eColorF=I
            END IF
          END IF
        END DO
        ListColor(iVertFound)=eColorF
      END DO
      deallocate(ListPosFirst)
      deallocate(CurrColor)
      ChromaticNr=maxval(ListColor)
# ifdef DEBUG
      WRITE(740+myrank,*) 'ChromaticNr=', ChromaticNr
      DO iVert=1,nbVert
        WRITE(740+myrank,*) 'iVert=', iVert, 'eColor=', ListColor(iVert)
      END DO
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DeallocateGraph(TheGraph)
      implicit none
      type(Graph), intent(inout) :: TheGraph
      deallocate(TheGraph % ListDegree)
      deallocate(TheGraph % ListEdge)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_BLOCK_FREQDIR(LocalColor, Nblock)
      USE DATAPOOL, only : MNP, MSC, MDC, LocalColorInfo, rkind
      USE DATAPOOL, only : wwm_nnbr_send, wwm_nnbr_recv
      USE DATAPOOL, only : wwm_ListNbCommon_send, wwm_ListNbCommon_recv
      USE elfe_msgp, only : rtype, ierr, myrank
      USE elfe_glbl, only : iplg
      USE DATAPOOL, only : XP, YP
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      integer, intent(in) :: Nblock
      integer Ntot, Hlen, Delta, iBlock, idx, ID, IS
      integer lenBlock, maxBlockLength
      integer IC, nbCommon, eFirst, I, idxSend, idxRecv
      REAL(rkind) :: eFieldStackA(MNP,3)
      Ntot=MyREAL(MSC*MDC)
      Hlen=INT(Ntot/Nblock)
      Delta=Ntot - Hlen*Nblock
      iBlock=1
      idx=1
      LocalColor % Nblock=Nblock
      IF (Delta == 0) THEN
        maxBlockLength=Hlen
      ELSE
        maxBlockLength=Hlen+1
      ENDIF
      allocate(LocalColor % ISindex(Nblock, maxBlockLength))
      allocate(LocalColor % IDindex(Nblock, maxBlockLength))
      DO IS=1,MSC
        DO ID=1,MDC
          LocalColor % ISindex(iBlock, idx)=IS
          LocalColor % IDindex(iBlock, idx)=ID
          IF (iBlock <= Delta) THEN
            lenBlock=Hlen+1
          ELSE
            lenBlock=Hlen
          END IF
          idx=idx+1
          IF (idx > lenBlock) THEN
            iBlock=iBlock+1
            idx=1
          ENDIF
        END DO
      END DO
      allocate(LocalColor % BlockLength(Nblock))
      DO iBlock=1,Nblock
        IF (iBlock <= Delta) THEN
          lenBlock=Hlen+1
        ELSE
          lenBlock=Hlen
        END IF
        LocalColor % BlockLength(iBlock)=lenBlock
      END DO
      LocalColor % maxBlockLength = maxBlockLength
      allocate(LocalColor % ACexch(maxBlockLength, MNP))
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_BLK_ARRAY(LocalColor)
      USE DATAPOOL, only : MNP, MSC, MDC, LocalColorInfo, rkind
      USE DATAPOOL, only : wwm_nnbr_send, wwm_nnbr_recv
      USE DATAPOOL, only : wwm_ListNbCommon_send, wwm_ListNbCommon_recv
      USE DATAPOOL, only : wwm_ListDspl_send, wwm_ListDspl_recv
      USE elfe_msgp, only : rtype, ierr, myrank
      USE elfe_glbl, only : iplg
      USE DATAPOOL, only : XP, YP
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      integer maxBlockLength, idxSend, idxRecv
      integer I, IC, nbCommon, eFirst
      integer ListFirstCommon_send(wwm_nnbr_send)
      integer ListFirstCommon_recv(wwm_nnbr_recv)
      integer, allocatable :: dspl_send(:), dspl_recv(:)
      maxBlockLength=LocalColor % maxBlockLength
      ListFirstCommon_send=0
      DO I=2,wwm_nnbr_send
        ListFirstCommon_send(I)=ListFirstCommon_send(I-1)+wwm_ListNbCommon_send(I-1)
      END DO
      ListFirstCommon_recv=0
      DO I=2,wwm_nnbr_recv
        ListFirstCommon_recv(I)=ListFirstCommon_recv(I-1)+wwm_ListNbCommon_recv(I-1)
      END DO
      allocate(LocalColor % blk_p2dsend_type(wwm_nnbr_send))
      allocate(LocalColor % blk_p2drecv_type(wwm_nnbr_recv))
# ifdef DEBUG
      WRITE(740+myrank,*) 'maxBlockLength=', maxBlockLength
# endif
      DO I=1,wwm_nnbr_send
        nbCommon=wwm_ListNbCommon_send(I)
        eFirst=ListFirstCommon_send(I)
        ALLOCATE(dspl_send(nbCommon))
        DO IC=1,nbCommon
          idxSend=wwm_ListDspl_send(eFirst+IC)
          dspl_send(IC)=idxSend-1
        END DO
        call mpi_type_create_indexed_block(nbCommon,maxBlockLength,dspl_send,rtype,LocalColor % blk_p2dsend_type(I),ierr)
        call mpi_type_commit(LocalColor % blk_p2dsend_type(I), ierr)
        DEALLOCATE(dspl_send)
      END DO
      DO I=1,wwm_nnbr_recv
        nbCommon=wwm_ListNbCommon_recv(I)
        eFirst=ListFirstCommon_recv(I)
        ALLOCATE(dspl_recv(nbCommon))
        DO IC=1,nbCommon
          idxRecv=wwm_ListDspl_recv(eFirst+IC)
          dspl_recv(IC)=idxRecv-1
        END DO
        call mpi_type_create_indexed_block(nbCommon,maxBlockLength,dspl_recv,rtype,LocalColor % blk_p2drecv_type(I),ierr)
        call mpi_type_commit(LocalColor % blk_p2drecv_type(I), ierr)
        DEALLOCATE(dspl_recv)
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SYMM_INIT_COLORING(LocalColor, NbBlock)
      USE DATAPOOL, only : LocalColorInfo, MNP, rkind, XP, YP
      USE DATAPOOL, only : DO_SOLVE_L, DO_SOLVE_U
      USE DATAPOOL, only : DO_SYNC_UPP_2_LOW, DO_SYNC_LOW_2_UPP, DO_SYNC_FINAL
      USE elfe_msgp, only : myrank, nproc
      USE elfe_glbl, only : iplg
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      integer, intent(in) :: NbBlock
      type(Graph) :: AdjGraph
      integer :: ListColor(nproc)
      integer :: ListColorWork(nproc)
      real(rkind) :: eFieldStackA(MNP,3)
      real(rkind) :: eFieldStackB(MNP,2)
      real(rkind) :: eFieldStackC(MNP,1)
      real(rkind) :: eFieldStackRevA(3,MNP)
      real(rkind) :: eFieldStackRevB(2,MNP)
      real(rkind) :: eFieldStackRevC(1,MNP)
      integer TheRes
# ifdef DEBUG
      CALL TEST_ASPAR_SYMMETRY
      WRITE(740+myrank,*) 'After TEST_ASPAR_SYMMETRY'
# endif
# ifdef DEBUG
      CALL COMPUTE_TOTAL_INDEX_SHIFT(TheRes)
      WRITE(740+myrank,*) 'Total residual shift=', TheRes
# endif
      CALL COLLECT_ALL_IPLG
# ifdef DEBUG
      WRITE(740+myrank,*) 'After COLLECT_ALL_IPLG'
# endif
      CALL COLLECT_ALL_IA_JA
# ifdef DEBUG
      WRITE(740+myrank,*) 'After COLLECT_ALL_IA_JA'
# endif
      CALL CREATE_WWM_P2D_EXCH
# ifdef DEBUG
      WRITE(740+myrank,*) 'After CREATE_WWM_P2D_EXCH'
# endif
      CALL CREATE_WWM_MAT_P2D_EXCH
# ifdef DEBUG
      WRITE(740+myrank,*) 'After CREATE_WWM_MAT_P2D_EXCH'
# endif
      CALL SYMM_GRAPH_BUILD_PROCESSOR_ADJACENCY(AdjGraph)
# ifdef DEBUG
      WRITE(740+myrank,*) 'After SYMM_GRAPH_BUILD_PROCESSOR_ADJACENCY'
# endif
      CALL GRAPH_TEST_UNIDIRECT(AdjGraph)
# ifdef DEBUG
      WRITE(740+myrank,*) 'After GRAPH_TEST_UNIDIRECT'
# endif
      CALL BUILD_MULTICOLORING(AdjGraph, ListColor)
# ifdef DEBUG
      WRITE(740+myrank,*) 'After BUILD_MULTICOLORING'
# endif
      CALL DeallocateGraph(AdjGraph)
      ListColorWork=-ListColor
      allocate(LocalColor % ListColor(nproc))
      LocalColor % ListColor=ListColorWork
# ifdef DEBUG
      WRITE(740+myrank,*) 'Before INIT_LOW_2_UPP_ARRAYS'
# endif
      CALL INIT_LOW_2_UPP_ARRAYS(LocalColor, ListColorWork)
      !
!# ifdef DEBUG
!      WRITE(740+myrank,*) 'Before INIT_UPP_2_LOW_ARRAYS (not needed)'
!# endif
!      CALL INIT_UPP_2_LOW_ARRAYS(LocalColor, ListColorWork)
      !
# ifdef DEBUG
      WRITE(740+myrank,*) 'Before CALL_BLOCK_FREQDIR'
# endif
      CALL INIT_BLOCK_FREQDIR(LocalColor, NbBlock)
      !
# ifdef DEBUG
      WRITE(740+myrank,*) 'Before INIT_BLK_ARRAY'
# endif
      CALL INIT_BLK_ARRAY(LocalColor)
      !
# ifdef DEBUG
      WRITE(740+myrank,*) 'Before COLLECT_ALL_COVLOWER'
# endif
      CALL COLLECT_ALL_COVLOWER(LocalColor)
      !
# ifdef DEBUG
      WRITE(740+myrank,*) 'Before INIT_COVLOWER_ARRAY'
# endif
      CALL INIT_COVLOWER_ARRAY(LocalColor)
      !
# ifdef DEBUG
      WRITE(740+myrank,*) 'Before DETERMINE_JSTATUS_L_U'
# endif
      CALL DETERMINE_JSTATUS_L_U(LocalColor)
      !
      DO_SOLVE_L=.TRUE.
      DO_SOLVE_U=.TRUE.
      DO_SYNC_UPP_2_LOW=.TRUE.
      DO_SYNC_LOW_2_UPP=.TRUE.
      DO_SYNC_FINAL=.TRUE.
# ifdef DEBUG
      WRITE(740+myrank,*) 'Checking XP exchanges'
      CALL CHECK_P2D_EXCH(LocalColor, XP)
      WRITE(740+myrank,*) 'Checking YP exchanges'
      CALL CHECK_P2D_EXCH(LocalColor, YP)
      WRITE(740+myrank,*) 'Checking IPLG exchanges'
      CALL CHECK_P2D_EXCH(LocalColor, MyREAL(iplg))
      !
      WRITE(740+myrank,*) 'Check Standard SELFE'
      CALL CHECK_STANDARD_SELFE_EXCH
      WRITE(740+myrank,*) 'Checking stacked exchanges A'
      eFieldStackA(:,1)=XP
      eFieldStackA(:,2)=YP
      eFieldStackA(:,3)=MyREAL(iplg)
!      CALL CHECK_P2D_EXCH_STACKED(LocalColor, eFieldStackA, 3)
      WRITE(740+myrank,*) 'Checking stacked exchanges B'
      eFieldStackB(:,1)=XP
      eFieldStackB(:,2)=YP
!      CALL CHECK_P2D_EXCH_STACKED(LocalColor, eFieldStackB, 2)
      WRITE(740+myrank,*) 'Checking stacked exchanges C'
      eFieldStackC(:,1)=XP
!      CALL CHECK_P2D_EXCH_STACKED(LocalColor, eFieldStackC, 1)
      WRITE(740+myrank,*) 'Checking stacked exchanges D'
      eFieldStackB(:,1)=XP
      eFieldStackB(:,2)=XP
!      CALL CHECK_P2D_EXCH_STACKED(LocalColor, eFieldStackB, 2)
      WRITE(740+myrank,*) 'Checking stacked exchanges E'
      eFieldStackB(:,1)=XP
      eFieldStackB(:,2)=0
!      CALL CHECK_P2D_EXCH_STACKED(LocalColor, eFieldStackB, 2)
      WRITE(740+myrank,*) 'Checking stacked exchanges F'
      eFieldStackB(:,1)=0
      eFieldStackB(:,2)=XP
!      CALL CHECK_P2D_EXCH_STACKED(LocalColor, eFieldStackB, 2)
      !
      WRITE(740+myrank,*) 'Checking stacked exchanges RevA'
      eFieldStackA(1,:)=XP
      eFieldStackA(2,:)=YP
      eFieldStackA(3,:)=MyREAL(iplg)
!      CALL CHECK_P2D_EXCH_STACKED_REV(LocalColor, eFieldStackRevA, 3)
      WRITE(740+myrank,*) 'Checking stacked exchanges RevB'
      eFieldStackB(1,:)=XP
      eFieldStackB(2,:)=YP
!      CALL CHECK_P2D_EXCH_STACKED_REV(LocalColor, eFieldStackRevB, 2)
      WRITE(740+myrank,*) 'Checking stacked exchanges RevC'
      eFieldStackC(1,:)=XP
!      CALL CHECK_P2D_EXCH_STACKED_REV(LocalColor, eFieldStackRevC, 1)
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_LOW_2_UPP_ARRAYS(LocalColor, ListColor)
      USE DATAPOOL, only : LocalColorInfo, MNP, rkind
      USE DATAPOOL, only : NNZ, IA, JA, NP_RES
      USE DATAPOOL, only : wwm_ListNbCommon_send, wwm_ListNbCommon_recv
      USE DATAPOOL, only : wwm_nnbr_send, wwm_nnbr_recv
      USE DATAPOOL, only : wwm_p2drecv_type, wwm_p2dsend_type
      USE DATAPOOL, only : wwm_ListNeigh_send, wwm_ListNeigh_recv
      USE DATAPOOL, only : wwm_ListDspl_recv, wwm_ListDspl_send
      USE elfe_msgp, only : myrank, nproc, comm, ierr, nbrrank_p
      implicit none
      include 'mpif.h'
      type(LocalColorInfo), intent(inout) :: LocalColor
      real(rkind) :: p2d_data_send(MNP)
      real(rkind) :: CovLower(MNP), CovLower_meth2(MNP), CovLower_meth3(MNP)
      real(rkind) :: SumErr, SumDiff
      integer, intent(in) :: ListColor(nproc)
      integer eColor, fColor, I, iRank, J
      integer nbLow_send, nbUpp_send, nbLow_recv, nbUpp_recv
      integer iProc, stat, eSize, iLow, iUpp, DoOper
      integer IC, eFirst, nbCommon, IP, IPloc, JP
      integer ListFirstCommon_send(wwm_nnbr_send)
      integer ListFirstCommon_recv(wwm_nnbr_recv)
      eColor=ListColor(myrank+1)
      nbUpp_send=0
      DO I=1,wwm_nnbr_send
        iRank=wwm_ListNeigh_send(I)
        fColor=ListColor(iRank)
# ifdef DEBUG
        WRITE(740+myrank,*) 'I=', I, 'iRank=', iRank, 'fColor=', fColor
        IF (fColor.eq.eColor) THEN
          call wwm_abort('Major error in the code')
        END IF
# endif
        IF (fColor.gt.eColor) THEN
          nbUpp_send=nbUpp_send + 1
        ENDIF
      END DO
      nbLow_recv=0
      DO I=1,wwm_nnbr_recv
        iRank=wwm_ListNeigh_recv(I)
        fColor=ListColor(iRank)
# ifdef DEBUG
        IF (fColor.eq.eColor) THEN
          call wwm_abort('Major error in the code')
        END IF
# endif
        IF (fColor.lt.eColor) THEN
          nbLow_recv=nbLow_recv + 1
        ENDIF
      END DO
      ListFirstCommon_send=0
      DO I=2,wwm_nnbr_send
        ListFirstCommon_send(I)=ListFirstCommon_send(I-1)+wwm_ListNbCommon_send(I-1)
      END DO
      ListFirstCommon_recv=0
      DO I=2,wwm_nnbr_recv
        ListFirstCommon_recv(I)=ListFirstCommon_recv(I-1)+wwm_ListNbCommon_recv(I-1)
      END DO
# ifdef DEBUG
      WRITE(740+myrank,*) 'SIC: nbLow_recv=', nbLow_recv, ' nbUpp_send=', nbUpp_send
# endif
      LocalColor % nbUpp_send=nbUpp_send
      LocalColor % nbLow_recv=nbLow_recv
      allocate(LocalColor % ListIdxUpper_send(nbUpp_send))
      allocate(LocalColor % ListIdxLower_recv(nbLow_recv))
      iUpp=0
      DO I=1,wwm_nnbr_send
        iRank=wwm_ListNeigh_send(I)
        fColor=ListColor(iRank)
        IF (fColor.gt.eColor) THEN
          iUpp=iUpp + 1
          LocalColor % ListIdxUpper_send(iUpp)=I
        ENDIF
      END DO
      iLow=0
      DO I=1,wwm_nnbr_recv
        iRank=wwm_ListNeigh_recv(I)
        fColor=ListColor(iRank)
        IF (fColor.lt.eColor) THEN
          iLow=iLow + 1
          LocalColor % ListIdxLower_recv(iLow)=I
        ENDIF
      END DO
      allocate(LocalColor % Upp_s_rq(nbUpp_send))
      allocate(LocalColor % Upp_s_stat(MPI_STATUS_SIZE, nbUpp_send))
      allocate(LocalColor % Low_r_rq(nbLow_recv))
      allocate(LocalColor % Low_r_stat(MPI_STATUS_SIZE, nbLow_recv))
      p2d_data_send=0
      CovLower=1
      CovLower_meth2=1
      DO iUpp=1,nbUpp_send
        I=LocalColor % ListIdxUpper_send(iUpp)
        iRank=wwm_ListNeigh_send(i)
# ifdef DEBUG
        WRITE(740+myrank,*) 'ISEND: iUpp=', iUpp, 'I=', I, 'iRank=', iRank
# endif
        call mpi_isend(p2d_data_send,1,wwm_p2dsend_type(i),iRank-1,13,comm,LocalColor % Upp_s_rq(iUpp),ierr)
      END DO
      DO iLow=1,nbLow_recv
        I=LocalColor % ListIdxLower_recv(iLow)
        iRank=wwm_ListNeigh_recv(i)
# ifdef DEBUG
        WRITE(740+myrank,*) 'IRECV: iLow=', iLow, 'I=', I, 'iRank=', iRank
# endif
        nbCommon=wwm_ListNbCommon_recv(I)
        eFirst=ListFirstCommon_recv(I)
        DO IC=1,nbCommon
          IPloc=wwm_ListDspl_recv(eFirst+IC)
          CovLower_meth2(IPloc)=0
        END DO
        call mpi_irecv(CovLower,1,wwm_p2drecv_type(i),iRank-1,13,comm,LocalColor % Low_r_rq(iLow),ierr)
      END DO
      IF (nbUpp_send > 0) THEN
        call mpi_waitall(nbUpp_send, LocalColor%Upp_s_rq, LocalColor%Upp_s_stat,ierr)
      END IF
      IF (nbLow_recv > 0) THEN
        call mpi_waitall(nbLow_recv, LocalColor%Low_r_rq, LocalColor%Low_r_stat,ierr)
      END IF
      allocate(LocalColor % CovLower(MNP))
      LocalColor % CovLower=INT(CovLower)
      SumErr=sum(abs(CovLower-CovLower_meth2))
# ifdef DEBUG
      WRITE(740+myrank,*) 'SumErr(meth1/meth2) CovLower=', SumErr
      WRITE(740+myrank,*) 'MNP=', MNP, ' sum(CovLower)=', sum(CovLower)
# endif
      IF (SumErr .gt. 0) THEN
        DO IP=1,MNP
          WRITE(740+myrank,*) IP, CovLower(IP), CovLower_meth2(IP)
        END DO
      ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_UPP_2_LOW_ARRAYS(LocalColor, ListColor)
      USE DATAPOOL, only : LocalColorInfo, MNP, rkind
      USE DATAPOOL, only : NNZ, IA, JA, NP_RES
      USE DATAPOOL, only : wwm_ListNbCommon_send, wwm_ListNbCommon_recv
      USE DATAPOOL, only : wwm_nnbr_send, wwm_nnbr_recv
      USE DATAPOOL, only : wwm_p2drecv_type, wwm_p2dsend_type
      USE DATAPOOL, only : wwm_ListNeigh_send, wwm_ListNeigh_recv
      USE DATAPOOL, only : wwm_ListDspl_recv, wwm_ListDspl_send
      USE elfe_msgp, only : myrank, nproc, comm, ierr, nbrrank_p
      implicit none
      include 'mpif.h'
      type(LocalColorInfo), intent(inout) :: LocalColor
      real(rkind) :: p2d_data_send(MNP)
      real(rkind) :: CovUpper(MNP), CovUpper_meth2(MNP), CovUpper_meth3(MNP)
      real(rkind) :: SumErr, SumDiff
      integer, intent(in) :: ListColor(nproc)
      integer eColor, fColor, I, iRank, J
      integer nbLow_send, nbUpp_send, nbLow_recv, nbUpp_recv
      integer iProc, stat, eSize, iLow, iUpp, DoOper
      integer IC, eFirst, nbCommon, IP, IPloc, JP
      integer ListFirstCommon_send(wwm_nnbr_send)
      integer ListFirstCommon_recv(wwm_nnbr_recv)
      eColor=ListColor(myrank+1)
# ifdef DEBUG
      WRITE(740+myrank,*) 'eColor=', eColor
      allocate(LocalColor % ListColor(nproc))
      LocalColor % ListColor=ListColor
# endif
      nbLow_send=0
      DO I=1,wwm_nnbr_send
        iRank=wwm_ListNeigh_send(I)
        fColor=ListColor(iRank)
# ifdef DEBUG
        WRITE(740+myrank,*) 'I=', I, 'iRank=', iRank, 'fColor=', fColor
        IF (fColor.eq.eColor) THEN
          call wwm_abort('Major error in the code')
        END IF
# endif
        IF (fColor.lt.eColor) THEN
          nbLow_send=nbLow_send + 1
        ENDIF
      END DO
      nbUpp_recv=0
      DO I=1,wwm_nnbr_recv
        iRank=wwm_ListNeigh_recv(I)
        fColor=ListColor(iRank)
        IF (fColor.eq.eColor) THEN
          call wwm_abort('Major error in the code')
        END IF
        IF (fColor.gt.eColor) THEN
          nbUpp_recv=nbUpp_recv + 1
        ENDIF
      END DO
      ListFirstCommon_send=0
      DO I=2,wwm_nnbr_send
        ListFirstCommon_send(I)=ListFirstCommon_send(I-1)+wwm_ListNbCommon_send(I-1)
      END DO
      ListFirstCommon_recv=0
      DO I=2,wwm_nnbr_recv
        ListFirstCommon_recv(I)=ListFirstCommon_recv(I-1)+wwm_ListNbCommon_recv(I-1)
      END DO
# ifdef DEBUG
      WRITE(740+myrank,*) 'SIC: nbLow_send=', nbLow_send, ' nbUpp_send=', nbUpp_send
      WRITE(740+myrank,*) 'SIC: nbLow_recv=', nbLow_recv, ' nbUpp_recv=', nbUpp_recv
# endif
      LocalColor % nbLow_send=nbLow_send
      LocalColor % nbUpp_recv=nbUpp_recv
      allocate(LocalColor % ListIdxLower_send(nbLow_send))
      allocate(LocalColor % ListIdxUpper_recv(nbUpp_recv))
      iLow=0
      DO I=1,wwm_nnbr_send
        iRank=wwm_ListNeigh_send(I)
        fColor=ListColor(iRank)
        IF (fColor.lt.eColor) THEN
          iLow=iLow + 1
          LocalColor % ListIdxLower_send(iLow)=I
        ENDIF
      END DO
      iUpp=0
      DO I=1,wwm_nnbr_recv
        iRank=wwm_ListNeigh_recv(I)
        fColor=ListColor(iRank)
        IF (fColor.gt.eColor) THEN
          iUpp=iUpp + 1
          LocalColor % ListIdxUpper_recv(iUpp)=I
        ENDIF
      END DO
      allocate(LocalColor % Low_s_rq(nbLow_send))
      allocate(LocalColor % Low_s_stat(MPI_STATUS_SIZE, nbLow_send))
      allocate(LocalColor % Upp_r_rq(nbUpp_recv))
      allocate(LocalColor % Upp_r_stat(MPI_STATUS_SIZE, nbUpp_recv))
      !
      !
      p2d_data_send=0
      CovUpper=1
      DO iLow=1,nbLow_send
        I=LocalColor % ListIdxLower_send(iLow)
        iRank=wwm_ListNeigh_send(i)
        call mpi_isend(p2d_data_send,1,wwm_p2dsend_type(i),iRank-1,17,comm,LocalColor%Low_s_rq(iLow),ierr)
      END DO
      DO iUpp=1,nbUpp_recv
        I=LocalColor % ListIdxUpper_recv(iUpp)
        iRank=wwm_ListNeigh_recv(i)
        nbCommon=wwm_ListNbCommon_recv(I)
        eFirst=ListFirstCommon_recv(I)
        DO IC=1,nbCommon
          IPloc=wwm_ListDspl_recv(eFirst+IC)
          CovUpper_meth2(IPloc)=0
        END DO
        call mpi_irecv(CovUpper,1,wwm_p2drecv_type(i),iRank-1,17,comm,LocalColor%Upp_r_rq(iUpp),ierr)
      END DO
      IF (nbLow_send > 0) THEN
        call mpi_waitall(nbLow_send, LocalColor%Low_s_rq, LocalColor%Low_s_stat,ierr)
      END IF
      IF (nbUpp_recv > 0) THEN
        call mpi_waitall(nbUpp_recv, LocalColor%Upp_r_rq, LocalColor%Upp_r_stat,ierr)
      END IF
      allocate(LocalColor % CovUpper(MNP), stat=stat)
      if(stat/=0) call wwm_abort('CovUpper: error in allocation')
      LocalColor % CovUpper=INT(CovUpper)
# ifdef DEBUG
      WRITE(740+myrank,*) 'MNP=', MNP, ' sum(CovUpper)=', sum(CovUpper)
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE INIT_COVLOWER_ARRAY(LocalColor)
      USE DATAPOOL
      USE elfe_msgp
      USE elfe_glbl
      implicit none
      include 'mpif.h'
      type(LocalColorInfo), intent(inout) :: LocalColor
      integer :: ListFirst(nproc)
      integer :: ListNeigh01(nproc)
      integer :: ListCommon_recv(nproc)
      integer :: ListCommon_send(nproc)
      integer :: ListMapped0(np_global)
      integer :: ListMapped1(np_global)
      integer :: ListMapped0_B(np_global)
      integer :: ListMapped1_B(np_global)
      integer, allocatable :: dspl_send(:), dspl_recv(:)
      integer IP, IP_glob, iProc, WeMatch, MNPloc, idx, NP_RESloc
      integer iNeigh, IPmap, eSize, eSizeRed, nbCommon
      integer nbCommon_send, nbCommon_recv, idx_send, idx_recv
      integer sumNbCommon_send, sumNbCommon_recv
      integer eExtent, eExtentRed, NewExtent, eLB, sizRType
      integer eType1, eType2
      integer u2l_nnbr_send, u2l_nnbr_recv
      integer sync_nnbr_send, sync_nnbr_recv
      integer eCov, eColor, fColor, iSync
      integer maxBlockLength
      integer nbCase1, nbCase2
      integer nbMap0, nbMap1
      integer DoOper
      ListFirst=0
      DO iProc=2,nproc
        ListFirst(iProc)=ListFirst(iProc-1) + ListMNP(iProc-1)
      END DO
      ListMapped0=0
      ListMapped1=0
      nbMap0=0
      nbMap1=0
      DO IP=1,MNP
        IP_glob=iplg(IP)
        eCov=LocalColor % CovLower(IP)
        IF (eCov == 0) THEN
          ListMapped0(IP_glob)=IP
          nbMap0=nbMap0+1
        END IF
        IF (eCov == 1) THEN
          ListMapped1(IP_glob)=IP
          nbMap1=nbMap1+1
        END IF
      END DO
# ifdef DEBUG
      WRITE(740+myrank,*) 'nbMap0=', nbMap0, ' nbMap1=', nbMap1
# endif
      !
      ! First the Upper to lower (u2l) block arrays 
      !
      u2l_nnbr_send=0
      u2l_nnbr_recv=0
      ListCommon_send=0
      ListCommon_recv=0
      eColor=LocalColor % ListColor(myrank+1)
# ifdef DEBUG
      WRITE(740+myrank,*) 'U2L eColor=', eColor
# endif
      DO iNeigh=1,wwm_nnbr
        iProc=wwm_ListNeigh(iNeigh)
        fColor=LocalColor % ListColor(iProc)
# ifdef DEBUG
        WRITE(740+myrank,*) 'U2L iNeigh=', iNeigh, ' iProc=', iProc, 'fColor=', fColor
# endif
        MNPloc=ListMNP(iProc)
        NP_RESloc=ListNP_RES(iProc)
        IF (fColor .ge. eColor) THEN
          nbCommon_recv=0
          DO IP=1,NP_RESloc
            IP_glob=ListIPLG(IP+ListFirst(iProc))
            eCov=LocalColor % ListCovLower(IP+ListFirst(iProc))
            IF (eCov .eq. 1) THEN
              IPmap=ListMapped0(IP_glob)
              IF (IPmap .gt. 0) THEN
                nbCommon_recv=nbCommon_recv+1
              ELSE
                IPmap=ListMapped1(IP_glob)
                IF ((IPmap .gt. 0).and.(IPmap .gt. NP_RES)) THEN
                  nbCommon_recv=nbCommon_recv+1
                END IF
              END IF
            END IF
          END DO
          IF (nbCommon_recv .gt. 0) THEN
            u2l_nnbr_recv=u2l_nnbr_recv+1
            ListCommon_recv(iProc)=nbCommon_recv
          END IF
# ifdef DEBUG
          WRITE(740+myrank,*) '   U2L nbCommon_recv=', nbCommon_recv
# endif
        END IF
        IF (fColor .le. eColor) THEN
          nbCommon_send=0
          DO IP=1,MNPloc
            IP_glob=ListIPLG(IP+ListFirst(iProc))
            eCov=LocalColor % ListCovLower(IP+ListFirst(iProc))
            IPmap=ListMapped1(IP_glob)
            IF ((IPmap .gt. 0).and.(IPmap .le. NP_RES)) THEN
              IF ((eCov .eq. 0).or.(IP.gt.NP_RESloc)) THEN
                nbCommon_send=nbCommon_send+1
              END IF
            END IF
          END DO
          IF (nbCommon_send .gt. 0) THEN
            u2l_nnbr_send=u2l_nnbr_send+1
            ListCommon_send(iProc)=nbCommon_send
          END IF
# ifdef DEBUG
          WRITE(740+myrank,*) '   U2L nbCommon_send=', nbCommon_send
# endif
        END IF
      END DO
# ifdef DEBUG
      WRITE(740+myrank,*) 'WWM_P2D: u2l_nnbr_send=', u2l_nnbr_send
      WRITE(740+myrank,*) 'WWM_P2D: u2l_nnbr_recv=', u2l_nnbr_recv
# endif
      LocalColor % u2l_nnbr_send=u2l_nnbr_send
      LocalColor % u2l_nnbr_recv=u2l_nnbr_recv
      allocate(LocalColor % u2l_ListNbCommon_send(u2l_nnbr_send))
      allocate(LocalColor % u2l_ListNbCommon_recv(u2l_nnbr_recv))
      allocate(LocalColor % u2l_ListNeigh_send(u2l_nnbr_send))
      allocate(LocalColor % u2l_ListNeigh_recv(u2l_nnbr_recv))
      idx_send=0
      idx_recv=0
      DO iProc=1,nproc
        IF (ListCommon_send(iProc) .gt. 0) THEN
          idx_send=idx_send+1
          LocalColor % u2l_ListNeigh_send(idx_send)=iProc
          nbCommon=ListCommon_send(iProc)
          LocalColor % u2l_ListNbCommon_send(idx_send)=nbCommon
        END IF
        IF (ListCommon_recv(iProc) .gt. 0) THEN
          idx_recv=idx_recv+1
          LocalColor % u2l_ListNeigh_recv(idx_recv)=iProc
          nbCommon=ListCommon_recv(iProc)
          LocalColor % u2l_ListNbCommon_recv(idx_recv)=nbCommon
        END IF
      END DO
      !
      ! Now creating the u2l exchange
      ! 
# ifdef DEBUG
      WRITE(740+myrank,*) 'WWM_P2D: wwm_nnbr=', wwm_nnbr
      WRITE(740+myrank,*) 'WWM_P2D: wwm_ListNeigh built'
# endif
      allocate(LocalColor % u2l_p2dsend_rqst(u2l_nnbr_send))
      allocate(LocalColor % u2l_p2drecv_rqst(u2l_nnbr_recv))
      allocate(LocalColor % u2l_p2dsend_stat(MPI_STATUS_SIZE,u2l_nnbr_send))
      allocate(LocalColor % u2l_p2drecv_stat(MPI_STATUS_SIZE,u2l_nnbr_recv))
      allocate(LocalColor % u2l_p2dsend_type(u2l_nnbr_send))
      allocate(LocalColor % u2l_p2drecv_type(u2l_nnbr_recv))
# ifdef DEBUG
      WRITE(740+myrank,*) 'WWM_P2D: alloc done'
# endif
      maxBlockLength=LocalColor % maxBlockLength
      DO iNeigh=1,u2l_nnbr_send
        iProc=LocalColor % u2l_ListNeigh_send(iNeigh)
        MNPloc=ListMNP(iProc)
        NP_RESloc=ListNP_RES(iProc)
        ListMapped0_B=0
        ListMapped1_B=0
        DO IP=1,MNPloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          eCov=LocalColor % ListCovLower(IP+ListFirst(iProc))
          IF (eCov == 0) THEN
            ListMapped0_B(IP_glob)=IP
          END IF
          IF (eCov == 1) THEN
            ListMapped1_B(IP_glob)=IP
          END IF
        END DO
        nbCommon=LocalColor % u2l_ListNbCommon_send(iNeigh)
        allocate(dspl_send(nbCommon))
        idx=0
        DO IP=1,NP_RES
          IF (LocalColor % CovLower(IP) .eq. 1) THEN
            IP_glob=iplg(IP)
            IPmap=ListMapped0_B(IP_glob)
            IF (IPmap .gt. 0) THEN
              idx=idx+1
              dspl_send(idx)=IP-1
            ELSE
              IPmap=ListMapped1_B(IP_glob)
              IF ((IPmap .gt. 0).and.(IPmap .gt. NP_RESloc)) THEN
                idx=idx+1
                dspl_send(idx)=IP-1
              END IF
            END IF
          END IF
        END DO
# ifdef DEBUG
        IF (idx .ne. nbCommon) THEN
          CALL WWM_ABORT('error in u2l_p2dsend')
        END IF
        WRITE(740+myrank,*) '   U2L idx=', idx, ' nbCommon=', nbCommon
# endif
        call mpi_type_create_indexed_block(nbCommon,maxBlockLength,dspl_send,rtype,LocalColor % u2l_p2dsend_type(iNeigh), ierr)
        call mpi_type_commit(LocalColor % u2l_p2dsend_type(iNeigh), ierr)
        deallocate(dspl_send)
      END DO
      DO iNeigh=1,u2l_nnbr_recv
        iProc=LocalColor % u2l_ListNeigh_recv(iNeigh)
        MNPloc=ListMNP(iProc)
        NP_RESloc=ListNP_RES(iProc)
        ListMapped0_B=0
        ListMapped1_B=0
        DO IP=1,MNPloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          eCov=LocalColor % ListCovLower(IP+ListFirst(iProc))
          IF (eCov == 0) THEN
            ListMapped0_B(IP_glob)=IP
          END IF
          IF (eCov == 1) THEN
            ListMapped1_B(IP_glob)=IP
          END IF
        END DO
        nbCommon=LocalColor % u2l_ListNbCommon_recv(iNeigh)
        allocate(dspl_recv(nbCommon))
        idx=0
        DO IP=1,NP_RESloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          eCov=LocalColor % ListCovLower(IP+ListFirst(iProc))
          IF (eCov == 1) THEN
            IPmap=ListMapped0(IP_glob)
            IF (IPmap .gt. 0) THEN
              idx=idx+1
              dspl_recv(idx)=IPmap-1
            ELSE
              IPmap=ListMapped1(IP_glob)
              IF ((IPmap .gt. 0).and.(IPmap .gt. NP_RES)) THEN
                idx=idx+1
                dspl_recv(idx)=IPmap-1
              END IF
            END IF
          END IF
        END DO
# ifdef DEBUG
        IF (idx .ne. nbCommon) THEN
          CALL WWM_ABORT('error in u2l_p2drecv')
        END IF
        WRITE(740+myrank,*) '   U2L idx=', idx, ' nbCommon=', nbCommon
# endif
        call mpi_type_create_indexed_block(nbCommon,maxBlockLength,dspl_recv,rtype,LocalColor % u2l_p2drecv_type(iNeigh),ierr)
        call mpi_type_commit(LocalColor % u2l_p2drecv_type(iNeigh), ierr)
        deallocate(dspl_recv)
      END DO
      !
      ! Now the synchronization arrays
      !
      sync_nnbr_send=0
      sync_nnbr_recv=0
      ListCommon_send=0
      ListCommon_recv=0
# ifdef DEBUG
      WRITE(740+myrank,*) 'wwm_nnbr=', wwm_nnbr
      WRITE(740+myrank,*) 'sum(ListCovLower)=', sum(LocalColor % ListCovLower)
# endif
      DO iNeigh=1,wwm_nnbr
        iProc=wwm_ListNeigh(iNeigh)
# ifdef DEBUG
        WRITE(740+myrank,*) 'iNeigh=', iNeigh, 'iProc=', iProc
# endif
        fColor=LocalColor % ListColor(iProc)
        MNPloc=ListMNP(iProc)
        NP_RESloc=ListNP_RES(iProc)
        ListMapped0_B=0
        ListMapped1_B=0
        DO IP=1,MNPloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          eCov=LocalColor % ListCovLower(IP+ListFirst(iProc))
          IF (eCov == 0) THEN
            ListMapped0_B(IP_glob)=IP
          END IF
          IF (eCov == 1) THEN
            ListMapped1_B(IP_glob)=IP
          END IF
        END DO
        nbCommon_recv=0
        nbCase1=0
        nbCase2=0
        DO IP=1,NP_RESloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          eCov=LocalColor % ListCovLower(IP+ListFirst(iProc))
          IF (eCov .eq. 1) THEN
            IF (ListMapped0(IP_glob) .gt. 0) THEN
              nbCommon_recv=nbCommon_recv+1
              nbCase1=nbCase1+1
            ELSE
              IPmap=ListMapped1(IP_glob)
              IF ((IPmap .gt. 0).and.(IPmap .gt. NP_RES)) THEN
                nbCommon_recv=nbCommon_recv+1
                nbCase2=nbCase2+1
              ENDIF
            END IF
          END IF
        END DO
# ifdef DEBUG
        WRITE(740+myrank,*) 'i=', iProc-1,  ' RnbCase12=', nbCase1, nbCase2
# endif
        IF (nbCommon_recv .gt. 0) THEN
          sync_nnbr_recv=sync_nnbr_recv+1
          ListCommon_recv(iProc)=nbCommon_recv
        END IF
        nbCommon_send=0
        DO IP=1,MNPloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          eCov=LocalColor % ListCovLower(IP+ListFirst(iProc))
          IPmap=ListMapped1(IP_glob)
          IF ((IPmap .gt. 0).and.(IPmap .le. NP_RES)) THEN
            IF ((eCov .eq. 0).or.(IP.gt.NP_RESloc)) THEN
              nbCommon_send=nbCommon_send+1
            END IF
          END IF
        END DO
        IF (nbCommon_send .gt. 0) THEN
          sync_nnbr_send=sync_nnbr_send+1
          ListCommon_send(iProc)=nbCommon_send
        END IF
# ifdef DEBUG
        WRITE(740+myrank,*) '   nbCommon(send/recv)=', nbCommon_send, nbCommon_recv
# endif
      END DO
# ifdef DEBUG
      WRITE(740+myrank,*) 'WWM_P2D: sync_nnbr_send=', sync_nnbr_send
      WRITE(740+myrank,*) 'WWM_P2D: sync_nnbr_recv=', sync_nnbr_recv
# endif
      LocalColor % sync_nnbr_send=sync_nnbr_send
      LocalColor % sync_nnbr_recv=sync_nnbr_recv
      allocate(LocalColor % sync_ListNbCommon_send(sync_nnbr_send))
      allocate(LocalColor % sync_ListNbCommon_recv(sync_nnbr_recv))
      allocate(LocalColor % sync_ListNeigh_send(sync_nnbr_send))
      allocate(LocalColor % sync_ListNeigh_recv(sync_nnbr_recv))
      idx_send=0
      idx_recv=0
      DO iProc=1,nproc
        IF (ListCommon_send(iProc) .gt. 0) THEN
          idx_send=idx_send+1
          LocalColor % sync_ListNeigh_send(idx_send)=iProc
          nbCommon=ListCommon_send(iProc)
          LocalColor % sync_ListNbCommon_send(idx_send)=nbCommon
        END IF
        IF (ListCommon_recv(iProc) .gt. 0) THEN
          idx_recv=idx_recv+1
          LocalColor % sync_ListNeigh_recv(idx_recv)=iProc
          nbCommon=ListCommon_recv(iProc)
          LocalColor % sync_ListNbCommon_recv(idx_recv)=nbCommon
        END IF
      END DO
      !
      ! Now creating the sync exchange
      !
      allocate(LocalColor % sync_p2dsend_rqst(sync_nnbr_send))
      allocate(LocalColor % sync_p2drecv_rqst(sync_nnbr_recv))
      allocate(LocalColor % sync_p2dsend_stat(MPI_STATUS_SIZE,sync_nnbr_send))
      allocate(LocalColor % sync_p2drecv_stat(MPI_STATUS_SIZE,sync_nnbr_recv))
      allocate(LocalColor % sync_p2dsend_type(sync_nnbr_send))
      allocate(LocalColor % sync_p2drecv_type(sync_nnbr_recv))
# ifdef DEBUG
      WRITE(740+myrank,*) 'SYNC sync_nnbr_send=', sync_nnbr_send
# endif
      DO iNeigh=1,sync_nnbr_send
        iProc=LocalColor % sync_ListNeigh_send(iNeigh)
        MNPloc=ListMNP(iProc)
        NP_RESloc=ListNP_RES(iProc)
        ListMapped0_B=0
        ListMapped1_B=0
        DO IP=1,MNPloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          eCov=LocalColor % ListCovLower(IP+ListFirst(iProc))
          IF (eCov == 0) THEN
            ListMapped0_B(IP_glob)=IP
          END IF
          IF (eCov == 1) THEN
            ListMapped1_B(IP_glob)=IP
          END IF
        END DO
        nbCommon=LocalColor % sync_ListNbCommon_send(iNeigh)
# ifdef DEBUG
        WRITE(740+myrank,*) '   SYNC iNeigh=', iNeigh, ' nbCommon=', nbCommon
# endif
        allocate(dspl_send(nbCommon))
        idx=0
        DO IP=1,NP_RES
          IF (LocalColor % CovLower(IP) .eq. 1) THEN
            DoOper=0
            IP_glob=iplg(IP)
            IPmap=ListMapped0_B(IP_glob)
            IF (IPmap .gt. 0) THEN
              DoOper=1
            ELSE
              IPmap=ListMapped1_B(IP_glob)
              IF ((IPmap .gt. 0).and.(IPmap.gt.NP_RESloc)) THEN
                DoOper=1
              END IF
            END IF
            IF (DoOper == 1) THEN
              idx=idx+1
              dspl_send(idx)=IP-1
# ifdef DEBUG
              WRITE(740+myrank,*) 'idx=', idx, 'IP=', IP
              WRITE(740+myrank,*) '  IP_glob=', IP_glob, 'IPmap=', IPmap
# endif
            END IF
          END IF
        END DO
# ifdef DEBUG
        IF (idx .ne. nbCommon) THEN
          CALL WWM_ABORT('error in SYNC_p2dsend')
        END IF
        WRITE(740+myrank,*) '   SYNC_p2dsend iProc=', iProc-1, ' nbCommon=', nbCommon
# endif
        call mpi_type_create_indexed_block(nbCommon,MSC*MDC,dspl_send,rtype,LocalColor % sync_p2dsend_type(iNeigh), ierr)
        call mpi_type_commit(LocalColor % sync_p2dsend_type(iNeigh), ierr)
        deallocate(dspl_send)
      END DO
# ifdef DEBUG
      WRITE(740+myrank,*) 'SYNC sync_nnbr_recv=', sync_nnbr_recv
# endif
      DO iNeigh=1,sync_nnbr_recv
        iProc=LocalColor % sync_ListNeigh_recv(iNeigh)
        MNPloc=ListMNP(iProc)
        NP_RESloc=ListNP_RES(iProc)
        ListMapped0_B=0
        ListMapped1_B=0
        DO IP=1,MNPloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          eCov=LocalColor % ListCovLower(IP+ListFirst(iProc))
          IF (eCov == 0) THEN
            ListMapped0_B(IP_glob)=IP
          END IF
          IF (eCov == 1) THEN
            ListMapped1_B(IP_glob)=IP
          END IF
        END DO
        nbCommon=LocalColor%sync_ListNbCommon_recv(iNeigh)
# ifdef DEBUG
        WRITE(740+myrank,*) '   SYNC iNeigh=', iNeigh, ' nbCommon=', nbCommon
# endif
        allocate(dspl_recv(nbCommon))
        idx=0
        DO IP=1,NP_RESloc
          IP_glob=ListIPLG(IP+ListFirst(iProc))
          eCov=LocalColor % ListCovLower(IP+ListFirst(iProc))
          IF (eCov == 1) THEN
            IPmap=ListMapped0(IP_glob)
            DoOper=0
            IF (IPmap .gt. 0) THEN
              DoOper=1
            ELSE
              IPmap=ListMapped1(IP_glob)
              IF ((IPmap .gt. 0).and.(IPmap.gt.NP_RES)) THEN
                DoOper=1
              END IF
            END IF
            IF (DoOper == 1) THEN
              idx=idx+1
              dspl_recv(idx)=IPmap-1
# ifdef DEBUG
              WRITE(740+myrank,*) 'idx=', idx, 'IP=', IP
              WRITE(740+myrank,*) '  IP_glob=', IP_glob, 'IPmap=', IPmap
# endif
            END IF
          END IF
        END DO
# ifdef DEBUG
        IF (idx .ne. nbCommon) THEN
          CALL WWM_ABORT('error in SYNC_p2drecv')
        END IF
        WRITE(740+myrank,*) '   SYNC_p2drecv iProc=', iProc-1, ' nbCommon=', nbCommon
# endif
# ifdef DEBUG
        WRITE(740+myrank,*) '   nbCase1=', nbCase1, ' nbCase2=', nbCase2
# endif
        call mpi_type_create_indexed_block(nbCommon,MSC*MDC,dspl_recv,rtype,LocalColor % sync_p2drecv_type(iNeigh),ierr)
        call mpi_type_commit(LocalColor % sync_p2drecv_type(iNeigh), ierr)
        deallocate(dspl_recv)
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DETERMINE_JSTATUS_L_U(LocalColor)
      USE DATAPOOL, only : LocalColorInfo, MNP, rkind
      USE DATAPOOL, only : NNZ, IA, JA, NP_RES
      USE DATAPOOL, only : wwm_nnbr_send, wwm_nnbr_recv
      USE DATAPOOL, only : wwm_p2drecv_type, wwm_p2dsend_type
      USE DATAPOOL, only : wwm_ListNeigh_send, wwm_ListNeigh_recv
      USE elfe_msgp, only : myrank, nproc, comm, ierr, nbrrank_p
      implicit none
      include 'mpif.h'
      type(LocalColorInfo), intent(inout) :: LocalColor
      integer Jstatus_L(NNZ), Jstatus_U(NNZ)
      integer IP, J, JP, DoOper
      Jstatus_L=0
      Jstatus_U=0
      DO IP=1,NP_RES
        IF (LocalColor%CovLower(IP) == 1) THEN
          DO J=IA(IP),IA(IP+1)-1
            JP=JA(J)
            IF (JP /= IP) THEN
              IF (JP .lt. IP) THEN
                DoOper=0
              ELSE
                IF (LocalColor % CovLower(JP) == 0) THEN
                  DoOper=0
                ELSE
                  DoOper=1
                END IF
              END IF
            ELSE
              DoOper=0
            END IF
            Jstatus_L(J)=DoOper
          END DO
          DO J=IA(IP),IA(IP+1)-1
            JP=JA(J)
            IF (JP /= IP) THEN
              IF (JP .lt. IP) THEN
                DoOper=1
              ELSE
                IF (LocalColor % CovLower(JP) == 0) THEN
                  DoOper=1
                ELSE
                  DoOper=0
                END IF
              END IF
            ELSE
              DoOper=0
            END IF
            Jstatus_U(J)=DoOper
          END DO
        END IF
      END DO
! Those are passive checks that we may pass
! and still be incorrect
# ifdef DEBUG
      DO IP=1,NP_RES
        DO J=IA(IP),IA(IP+1)-1
          JP=JA(J)
          IF ((Jstatus_L(J).eq.1).and.(Jstatus_U(J).eq.1)) THEN
            WRITE(myrank+919,*) 'MNP=', MNP, 'NP_RES=', NP_RES
            WRITE(myrank+919,*) 'IP=', IP, ' JP=', JP
            WRITE(myrank+919,*) 'IPcovLower=', LocalColor%CovLower(IP)
            WRITE(myrank+919,*) 'JPcovLower=', LocalColor%CovLower(JP)
            WRITE(myrank+919,*) 'We have major error'
            CALL WWM_ABORT('Please panic and debug')
          END IF
        END DO
      END DO
      WRITE(myrank+740,*) 'sum(Jstatus_L)=', sum(Jstatus_L)
      WRITE(myrank+740,*) 'sum(Jstatus_U)=', sum(Jstatus_U)
# endif
      allocate(LocalColor % Jstatus_L(NNZ))
      allocate(LocalColor % Jstatus_U(NNZ))
      LocalColor % Jstatus_L=Jstatus_L
      LocalColor % Jstatus_U=Jstatus_U
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5_RECV_ASPAR_PC(LocalColor, SolDat)
      USE DATAPOOL, only : LocalColorInfo, I5_SolutionData
      USE DATAPOOL, only : MNP, MSC, MDC, rkind
      USE elfe_msgp, only : ierr, comm, rtype, istatus, nbrrank_p
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      type(I5_SolutionData), intent(inout) :: SolDat
      real(rkind), allocatable :: ASPAR_rs(:)
      integer idx, iNNZ, jNNZ, IS, ID, NNZ_l, siz
      integer iProc, i, iRank
      DO iProc=1,LocalColor % nbLow_send
        i=LocalColor % ListIdxUpper_send(iProc)
        iRank=nbrrank_p(i)
        NNZ_l=LocalColor % NNZ_len_r(iProc)
        siz=NNZ_l*MSC*MDC
        allocate(ASPAR_rs(NNZ_l*MSC*MDC))
        CALL mpi_recv(ASPAR_rs,siz,rtype,iRank,45,comm,istatus,ierr)
        idx=0
        DO iNNZ=1,NNZ_l
          jNNZ=LocalColor % NNZ_index_r(iProc,iNNZ)
          DO IS=1,MSC
            DO ID=1,MDC
              idx=idx+1
              SolDat % ASPAR_pc(jNNZ,IS,ID)=ASPAR_rs(idx)
            END DO
          END DO
        END DO
        deallocate(ASPAR_rs)
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5_SEND_ASPAR_PC(LocalColor, SolDat)
      USE DATAPOOL, only : LocalColorInfo, I5_SolutionData
      USE DATAPOOL, only : MSC, MDC, rkind
      USE elfe_msgp, only : ierr, comm, rtype, nbrrank_p
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      type(I5_SolutionData), intent(inout) :: SolDat
      real(rkind), allocatable :: ASPAR_rs(:)
      integer idx, iNNZ, jNNZ, IS, ID, NNZ_l, siz
      integer iProc, iRank, i
      DO iProc=1,LocalColor % nbUpp_send
        i=LocalColor % ListIdxUpper_send(iProc)
        iRank=nbrrank_p(i)
        NNZ_l=LocalColor % NNZ_len_s(iProc)
        siz=NNZ_l*MSC*MDC
        allocate(ASPAR_rs(NNZ_l*MSC*MDC))
        idx=0
        DO iNNZ=1,NNZ_l
          jNNZ=LocalColor % NNZ_index_s(iProc,iNNZ)
          DO IS=1,MSC
            DO ID=1,MDC
              idx=idx+1
              ASPAR_rs(idx)=SolDat % ASPAR_pc(jNNZ,IS,ID)
            END DO
          END DO
        END DO
        CALL mpi_isend(ASPAR_rs,siz,rtype,iRank,45,comm,LocalColor%Upp_s_rq(iProc),ierr)
        deallocate(ASPAR_rs)
      END DO
      IF (LocalColor % nbUpp_send > 0) THEN
        call mpi_waitall(LocalColor %nbUpp_send, LocalColor % Upp_s_rq, LocalColor % Upp_s_stat,ierr)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5_CREATE_PRECOND_ILU0(LocalColor, SolDat)
      USE DATAPOOL, only : MNP, NP_RES, LocalColorInfo, I5_SolutionData
      USE DATAPOOL, only : MSC, MDC, IA, JA, I_DIAG, rkind
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      type(I5_SolutionData), intent(inout) :: SolDat
      integer IP, JP, J, JP2, J2, J_FOUND, IS, ID
      integer, allocatable :: ListJ(:)
      real(rkind) tl
      SolDat%ASPAR_pc=SolDat%ASPAR_block
      CALL I5_RECV_ASPAR_PC(LocalColor, SolDat)
      allocate(ListJ(MNP))
      DO IP=1,NP_RES
        IF (LocalColor % CovLower(IP) == 1) THEN
          DO J=IA(IP),IA(IP+1)-1
            JP=JA(J)
            ListJ(JP)=J
          END DO
          DO J=IA(IP),I_DIAG(IP)-1
            JP=JA(J)
            DO IS=1,MSC
              DO ID=1,MDC
                tl=SolDat%ASPAR_pc(J,IS,ID)*SolDat%ASPAR_pc(I_DIAG(JP),IS,ID)
                DO J2=IA(JP),IA(JP+1)-1
                  JP2=JA(J2)
                  J_FOUND=ListJ(JP2)
                  IF (J_FOUND.gt.0) THEN ! Here is ILU0 approximation
                    SolDat%ASPAR_pc(J_FOUND,IS,ID)=SolDat%ASPAR_pc(J_FOUND,IS,ID) - tl*SolDat%ASPAR_pc(J2,IS,ID)
                  END IF
                END DO
                SolDat%ASPAR_pc(J,IS,ID)=tl
              END DO
            END DO
          END DO
          DO J=IA(IP),IA(IP+1)-1
            JP=JA(J)
            ListJ(JP)=0
          END DO
          J=I_DIAG(IP)
          DO IS=1,MSC
            DO ID=1,MDC
              SolDat%ASPAR_pc(J,IS,ID)=1.0_rkind/SolDat%ASPAR_pc(J,IS,ID)
            END DO
          END DO
        END IF
      END DO
      deallocate(ListJ)
      CALL I5_SEND_ASPAR_PC(LocalColor, SolDat)
      END SUBROUTINE
!**********************************************************************
!* We assign the values only for CovLower(IP)=1                       *
!* We could with some effort assign values for all with some effort   *
!* but the values would not be used                                   *
!**********************************************************************
      SUBROUTINE I5_CREATE_PRECOND_SOR(LocalColor, SolDat)
      USE DATAPOOL, only : MNP, MSC, MDC, IA, JA, I_DIAG, NP_RES, rkind, ONE
      USE DATAPOOL, only : LocalColorInfo, I5_SolutionData
# ifdef DEBUG
      USE elfe_msgp, only : myrank
      USE elfe_glbl, only : iplg
# endif
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      type(I5_SolutionData), intent(inout) :: SolDat
      integer IP, ID, IS, JP, J1, J, IPglob, JPglob
      real(rkind) eVal
# ifdef DEBUG
      integer ISsel, IDsel
      integer eCov
      real(rkind) eCoeff
# endif
      DO IP=1,NP_RES
        J=I_DIAG(IP)
        DO IS=1,MSC
          DO ID=1,MDC
            eVal=ONE/SolDat % ASPAR_block(J,IS,ID)
            SolDat%AC2(IP,IS,ID)=eVal
          END DO
        END DO
      END DO
      CALL EXCHANGE_P4D_WWM_TR(SolDat%AC2)
# ifdef DEBUG
      ISsel=2
      IDsel=2
      DO IP=1,MNP
        IPglob=iplg(IP)
        eCoeff=SolDat%AC2(IP,ISsel,IDsel)
        J=I_DIAG(IP)
        WRITE(myrank+267,*) 'IP/IPglob/dia=', IP, IPglob, eCoeff, J
      END DO
      WRITE(myrank+267,*) 'ENDING'
      WRITE(myrank+277,*) 'NP_RES=', NP_RES
      DO IP=1,MNP
        IPglob=iplg(IP)
        eCoeff=SolDat%AC2(IP,ISsel,IDsel)
        eCov=LocalColor%CovLower(IP)
        WRITE(myrank+277,*) IPglob, eCoeff, eCov
      END DO
      WRITE(myrank+277,*) 'ENDING'
# endif
# ifdef DEBUG
      DO IP=1,MNP
        DO J=IA(IP),IA(IP+1)-1
          JP=JA(J)
          IPglob=iplg(IP)
          JPglob=iplg(JP)
          eCoeff=SolDat % ASPAR_block(J,ISsel,IDsel)
          WRITE(myrank+167,*) IPglob, JPglob, J, eCoeff
        END DO
      END DO
# endif
      DO IP=1,NP_RES
        IF (LocalColor%CovLower(IP) == 1) THEN
          DO J=IA(IP),IA(IP+1)-1
            IF (LocalColor%Jstatus_L(J) == 1) THEN
              JP=JA(J)
# ifdef DEBUG
              IPglob=iplg(IP)
              JPglob=iplg(JP)
              eCoeff=SolDat % ASPAR_block(J,ISsel,IDsel)
              WRITE(myrank+167,*) IPglob, JPglob, J, eCoeff
# endif
              DO IS=1,MSC
                DO ID=1,MDC
                  eVal=SolDat%AC2(JP,IS,ID)
                  SolDat % ASPAR_pc(J,IS,ID)=SolDat % ASPAR_block(J,IS,ID)*eVal
# ifdef DEBUG
                  IF ((IS.eq.2).and.(ID.eq.2)) THEN
                    WRITE(myrank+177,*) 'iplg(IP/JP), J=', iplg(IP), iplg(JP), J
                    WRITE(myrank+177,*) 'eVal/pc/block=', eVal, SolDat % ASPAR_pc(J,IS,ID), SolDat % ASPAR_block(J,IS,ID)
                  END IF
# endif
                END DO
              END DO
            END IF
          ENDDO
          J=I_DIAG(IP)
          DO IS=1,MSC
            DO ID=1,MDC
              SolDat % ASPAR_pc(J,IS,ID)=SolDat%AC2(IP,IS,ID)
            END DO
          END DO
          DO J=IA(IP),IA(IP+1)-1
            IF (LocalColor%Jstatus_U(J) == 1) THEN
              JP=JA(J)
              DO IS=1,MSC
                DO ID=1,MDC
                  SolDat % ASPAR_pc(J,IS,ID)=SolDat % ASPAR_block(J,IS,ID)
                END DO
              END DO
            END IF
          END DO
        END IF
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5_CREATE_PRECOND(LocalColor, SolDat, TheMethod)
      USE DATAPOOL, only : LocalColorInfo, I5_SolutionData
# ifdef DEBUG
      USE elfe_msgp, only : myrank
# endif
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      type(I5_SolutionData), intent(inout) :: SolDat
      integer, intent(in) :: TheMethod
      IF (TheMethod == 1) THEN ! SOR 
        CALL I5_CREATE_PRECOND_SOR(LocalColor, SolDat)
      ELSE IF (TheMethod == 2) THEN ! ILU0
        CALL I5_CREATE_PRECOND_ILU0(LocalColor, SolDat)
      ELSE
        CALL WWM_ABORT('Wrong choice of preconditioner')
      ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5_SYNC_UPP_2_LOW_SENDRECV(LocalColor, AC)
      USE DATAPOOL, only : LocalColorInfo, MNP, MSC, MDC, rkind
      USE DATAPOOL, only : wwm_ListNeigh_send, wwmtot_p2dsend_type
      USE DATAPOOL, only : wwm_ListNeigh_recv, wwmtot_p2drecv_type
      USE elfe_msgp, only : comm, ierr, myrank
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      real(rkind), intent(inout) :: AC(MNP, MSC, MDC)
      real(rkind) :: U(MSC,MDC,MNP)
      INTEGER :: IP, IS, ID, i, iRank
      integer iUpp, iLow
      integer nbUpp_recv, nbLow_send
      nbLow_send=LocalColor % nbLow_send
      nbUpp_recv=LocalColor % nbUpp_recv
# ifdef DEBUG
      write(myrank+740,*) 'nbLow_send=', nbLow_send, 'nbUpp_recv=', nbUpp_recv
# endif
      DO IP = 1, MNP
        DO IS = 1, MSC
          DO ID = 1, MDC
            U(IS,ID,IP) = AC(IP,IS,ID)
          END DO
        END DO
      END DO
      DO iLow=1,nbLow_send
        i=LocalColor % ListIdxLower_send(iLow)
        iRank=wwm_ListNeigh_send(i)-1
        CALL mpi_isend(U, 1, wwmtot_p2dsend_type(i), iRank, 47, comm, LocalColor%Low_s_rq(iLow), ierr)
      END DO
      DO iUpp=1,nbUpp_recv
        I=LocalColor%ListIdxUpper_recv(iUpp)
        iRank=wwm_ListNeigh_recv(i)
        call mpi_irecv(U,1,wwmtot_p2drecv_type(i),iRank-1,47,comm,LocalColor%Upp_r_rq(iUpp),ierr)
      END DO
      IF (nbUpp_recv > 0) THEN
        call mpi_waitall(nbUpp_recv, LocalColor%Upp_r_rq, LocalColor%Upp_r_stat,ierr)
      END IF
      IF (nbLow_send > 0) THEN
        call mpi_waitall(nbLow_send, LocalColor % Low_s_rq, LocalColor % Low_s_stat,ierr)
      END IF
      DO IP = 1, MNP
        DO IS = 1, MSC
          DO ID = 1, MDC
            AC(IP,IS,ID) = U(IS,ID,IP)
          END DO
        END DO
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5_SYNC_LOW_2_UPP_SENDRECV(LocalColor, AC)
      USE DATAPOOL, only : LocalColorInfo, MNP, MSC, MDC, rkind
      USE DATAPOOL, only : wwm_ListNeigh_send, wwmtot_p2dsend_type
      USE DATAPOOL, only : wwm_ListNeigh_recv, wwmtot_p2drecv_type
      USE elfe_msgp, only : comm, ierr, myrank
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      real(rkind), intent(inout) :: AC(MNP, MSC, MDC)
      real(rkind) :: U(MSC,MDC,MNP)
      INTEGER :: IP, IS, ID, i, iRank
      integer iUpp, iLow
      integer nbUpp_send, nbLow_recv
      nbLow_recv=LocalColor%nbLow_recv
      nbUpp_send=LocalColor%nbUpp_send
# ifdef DEBUG
      write(myrank+740,*) 'nbLow_recv=', nbLow_recv, ' nbUpp_send=', nbUpp_send
# endif
      DO IP = 1, MNP
        DO IS = 1, MSC
          DO ID = 1, MDC
            U(IS,ID,IP) = AC(IP,IS,ID)
          END DO
        END DO
      END DO
      DO iUpp=1,nbUpp_send
        i=LocalColor % ListIdxUpper_send(iUpp)
        iRank=wwm_ListNeigh_send(i)
        write(myrank+740,*) 'iUpp=', iUpp, 'I=', I, 'iRank=', iRank
        CALL mpi_isend(U, 1, wwmtot_p2dsend_type(i), iRank-1, 59, comm, LocalColor%Upp_s_rq(iUpp), ierr)
      END DO
      DO iLow=1,nbLow_recv
        I=LocalColor%ListIdxLower_recv(iLow)
        iRank=wwm_ListNeigh_recv(i)
# ifdef DEBUG
        write(myrank+740,*) 'iLow=', iLow, 'I=', I, 'iRank=', iRank
# endif
        call mpi_irecv(U,1,wwmtot_p2drecv_type(i),iRank-1,59,comm,LocalColor%Low_r_rq(iLow),ierr)
      END DO
      IF (nbUpp_send > 0) THEN
        call mpi_waitall(nbUpp_send, LocalColor%Upp_s_rq, LocalColor%Upp_s_stat,ierr)
      END IF
      IF (nbLow_recv > 0) THEN
        call mpi_waitall(nbLow_recv, LocalColor % Low_r_rq, LocalColor % Low_r_stat,ierr)
      END IF
      DO IP = 1, MNP
        DO IS = 1, MSC
          DO ID = 1, MDC
            AC(IP,IS,ID) = U(IS,ID,IP)
          END DO
        END DO
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5_EXCHANGE_P3_LOW_2_UPP_Send(LocalColor, AC, iBlock)
      USE DATAPOOL, only : LocalColorInfo, MNP, MSC, MDC, rkind
      USE DATAPOOL, only : wwm_ListNeigh_send, wwm_p2dsend_type
      USE elfe_msgp, only : comm, ierr, myrank
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
# ifdef DEBUG
      real(rkind) :: Hfield(MNP)
# endif
      REAL(rkind), intent(in) :: AC(MNP, MSC, MDC)
      INTEGER, intent(in) :: iBlock
      integer iUpp, i, iRank, idx, lenBlock, maxBlockLength, IS, ID, nbUpp_send
      lenBlock=LocalColor % BlockLength(iBlock)
      maxBlockLength=LocalColor % maxBlockLength
      DO idx=1,lenBlock
        IS=LocalColor % ISindex(iBlock, idx)
        ID=LocalColor % IDindex(iBlock, idx)
        LocalColor % ACexch(idx,:)=AC(:,IS,ID)
      END DO
# ifdef DEBUG
      write(myrank+740,*) 'I5_P3_LOW_2_UPP_Send(I), iBlock=', iBlock
      DO idx=1,lenBlock
        IS=LocalColor % ISindex(iBlock, idx)
        ID=LocalColor % IDindex(iBlock, idx)
        Hfield(:)=AC(:,IS,ID)
        write(myrank+740,*) IS, ID, sum(Hfield)
      END DO
# endif
      nbUpp_send=LocalColor % nbUpp_send
      DO iUpp=1,nbUpp_send
        i=LocalColor % ListIdxUpper_send(iUpp)
        iRank=wwm_ListNeigh_send(i)-1
        CALL mpi_isend(LocalColor % ACexch, 1, LocalColor % blk_p2dsend_type(i), iRank, 7, comm, LocalColor%Upp_s_rq(iUpp), ierr)
      END DO
      IF (nbUpp_send > 0) THEN
        call mpi_waitall(nbUpp_send, LocalColor % Upp_s_rq, LocalColor % Upp_s_stat,ierr)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5_EXCHANGE_P3_LOW_2_UPP_Recv(LocalColor, AC, iBlock)
      USE DATAPOOL, only : LocalColorInfo, MNP, MSC, MDC, rkind
      USE DATAPOOL, only : wwm_ListNeigh_recv, wwm_p2drecv_type
      USE elfe_msgp, only : comm, ierr, myrank
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      REAL(rkind), intent(inout) :: AC(MNP, MSC, MDC)
# ifdef DEBUG
      real(rkind) :: Hfield(MNP)
# endif
      INTEGER, intent(in) :: iBlock
      integer iProc, i, iRank, idx, lenBlock, maxBlockLength, IS, ID, nbLow_recv
      lenBlock=LocalColor % BlockLength(iBlock)
      maxBlockLength=LocalColor % maxBlockLength
      DO idx=1,lenBlock
        IS=LocalColor % ISindex(iBlock, idx)
        ID=LocalColor % IDindex(iBlock, idx)
        LocalColor % ACexch(idx,:)=AC(:,IS,ID)
      END DO
# ifdef DEBUG
      write(myrank+740,*) 'I5_P3_LOW_2_UPP_Recv(I), iBlock=', iBlock
      DO idx=1,lenBlock
        IS=LocalColor % ISindex(iBlock, idx)
        ID=LocalColor % IDindex(iBlock, idx)
        Hfield(:)=AC(:,IS,ID)
        write(myrank+740,*) IS, ID, sum(Hfield)
      END DO
# endif
      nbLow_recv=LocalColor % nbLow_recv
      DO iProc=1,nbLow_recv
        i=LocalColor % ListIdxLower_recv(iProc)
        iRank=wwm_ListNeigh_recv(i)-1
        call mpi_irecv(LocalColor % ACexch,1,LocalColor % blk_p2drecv_type(i),iRank,7,comm,LocalColor % Low_r_rq(iProc),ierr)
      END DO
      IF (nbLow_recv > 0) THEN
        call mpi_waitall(nbLow_recv, LocalColor%Low_r_rq, LocalColor%Low_r_stat,ierr)
      END IF
      DO idx=1,lenBlock
        IS=LocalColor % ISindex(iBlock, idx)
        ID=LocalColor % IDindex(iBlock, idx)
        AC(:,IS,ID)=LocalColor % ACexch(idx,:)
      END DO
# ifdef DEBUG
      write(myrank+740,*) 'I5_P3_LOW_2_UPP_Recv(II), iBlock=', iBlock
      DO idx=1,lenBlock
        IS=LocalColor % ISindex(iBlock, idx)
        ID=LocalColor % IDindex(iBlock, idx)
        Hfield(:)=AC(:,IS,ID)
        write(myrank+740,*) IS, ID, sum(Hfield)
      END DO
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5_EXCHANGE_P3_UPP_2_LOW_Send_V1(LocalColor, AC, iBlock)
      USE DATAPOOL, only : LocalColorInfo, MNP, MSC, MDC, rkind
      USE DATAPOOL, only : wwm_ListNeigh_send, wwm_p2dsend_type
      USE elfe_msgp, only : comm, ierr, myrank
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      REAL(rkind), intent(in) :: AC(MNP, MSC, MDC)
      INTEGER, intent(in) :: iBlock
# ifdef DEBUG
      real(rkind) :: Hfield(MNP)
# endif
      integer iProc, i, iRank, idx, lenBlock, maxBlockLength, IS, ID, nbLow_send
      lenBlock=LocalColor % BlockLength(iBlock)
      maxBlockLength=LocalColor % maxBlockLength
      DO idx=1,lenBlock
        IS=LocalColor % ISindex(iBlock, idx)
        ID=LocalColor % IDindex(iBlock, idx)
        LocalColor % ACexch(idx,:)=AC(:,IS,ID)
      END DO
# ifdef DEBUG
      write(myrank+740,*) 'I5_P3_UPP_2_LOW_Send(I), iBlock=', iBlock
      DO idx=1,lenBlock
        IS=LocalColor % ISindex(iBlock, idx)
        ID=LocalColor % IDindex(iBlock, idx)
        Hfield(:)=AC(:,IS,ID)
        write(myrank+740,*) IS, ID, sum(Hfield)
      END DO
# endif
      nbLow_send=LocalColor % nbLow_send
      DO iProc=1,nbLow_send
        i=LocalColor % ListIdxLower_send(iProc)
        iRank=wwm_ListNeigh_send(i)-1
        call mpi_isend(LocalColor % ACexch,1,LocalColor%blk_p2dsend_type(i),iRank,3,comm,LocalColor%Low_s_rq(iProc),ierr)
      END DO
      IF (nbLow_send > 0) THEN
        call mpi_waitall(nbLow_send, LocalColor%Low_s_rq, LocalColor%Low_s_stat,ierr)
      END IF
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5_EXCHANGE_P3_UPP_2_LOW_Recv_V1(LocalColor, AC, iBlock)
      USE DATAPOOL, only : LocalColorInfo, MNP, MSC, MDC, rkind
      USE DATAPOOL, only : wwm_ListNeigh_recv, wwm_p2drecv_type
      USE elfe_msgp, only : comm, ierr, myrank
      implicit none
      include 'mpif.h'
# ifdef DEBUG
      real(rkind) :: Hfield(MNP)
# endif
      type(LocalColorInfo), intent(inout) :: LocalColor
      REAL(rkind), intent(inout) :: AC(MNP, MSC, MDC)
      INTEGER, intent(in) :: iBlock
      integer iProc, i, iRank, idx, lenBlock, maxBlockLength
      integer nbUpp_recv, IS, ID
      lenBlock=LocalColor % BlockLength(iBlock)
      maxBlockLength=LocalColor % maxBlockLength
      DO idx=1,lenBlock
        IS=LocalColor % ISindex(iBlock, idx)
        ID=LocalColor % IDindex(iBlock, idx)
        LocalColor % ACexch(idx,:)=AC(:,IS,ID)
      END DO
# ifdef DEBUG
      write(myrank+740,*) 'I5_P3_UPP_2_LOW_Recv(I), iBlock=', iBlock
      DO idx=1,lenBlock
        IS=LocalColor % ISindex(iBlock, idx)
        ID=LocalColor % IDindex(iBlock, idx)
        Hfield(:)=AC(:,IS,ID)
        write(myrank+740,*) IS, ID, sum(Hfield)
      END DO
# endif
      nbUpp_recv=LocalColor % nbUpp_recv
      DO iProc=1,nbUpp_recv
        i=LocalColor % ListIdxUpper_recv(iProc)
        iRank=wwm_ListNeigh_recv(i)-1
        call mpi_irecv(LocalColor % ACexch,1,LocalColor % blk_p2drecv_type(i),iRank,3,comm,LocalColor % Upp_r_rq(iProc),ierr)
      END DO
      IF (nbUpp_recv > 0) THEN
        call mpi_waitall(nbUpp_recv, LocalColor%Upp_r_rq, LocalColor%Upp_r_stat,ierr)
      END IF
      DO idx=1,lenBlock
        IS=LocalColor % ISindex(iBlock, idx)
        ID=LocalColor % IDindex(iBlock, idx)
        AC(:,IS,ID)=LocalColor % ACexch(idx,:)
      END DO
# ifdef DEBUG
      write(myrank+740,*) 'I5_P3_UPP_2_LOW_Recv(II), iBlock=', iBlock
      DO idx=1,lenBlock
        IS=LocalColor % ISindex(iBlock, idx)
        ID=LocalColor % IDindex(iBlock, idx)
        Hfield(:)=AC(:,IS,ID)
        write(myrank+740,*) IS, ID, sum(Hfield)
      END DO
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5_EXCHANGE_P3_UPP_2_LOW_Send(LocalColor, AC, iBlock)
      USE DATAPOOL, only : LocalColorInfo, MNP, MSC, MDC, rkind
      USE DATAPOOL, only : wwm_ListNeigh_send, wwm_p2dsend_type
      USE elfe_msgp, only : comm, ierr, myrank
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      REAL(rkind), intent(in) :: AC(MNP, MSC, MDC)
      INTEGER, intent(in) :: iBlock
# ifdef DEBUG
      real(rkind) :: Hfield(MNP)
# endif
      integer iProc, iRank, idx, lenBlock, maxBlockLength, IS, ID, nbLow_send
      lenBlock=LocalColor % BlockLength(iBlock)
      DO idx=1,lenBlock
        IS=LocalColor % ISindex(iBlock, idx)
        ID=LocalColor % IDindex(iBlock, idx)
        LocalColor % ACexch(idx,:)=AC(:,IS,ID)
      END DO
# ifdef DEBUG
      write(myrank+740,*) 'I5_P3_UPP_2_LOW_Send(I), iBlock=', iBlock
      DO idx=1,lenBlock
        IS=LocalColor % ISindex(iBlock, idx)
        ID=LocalColor % IDindex(iBlock, idx)
        Hfield(:)=AC(:,IS,ID)
        write(myrank+740,*) IS, ID, sum(Hfield)
      END DO
# endif
      nbLow_send=LocalColor % u2l_nnbr_send
      DO iProc=1,nbLow_send
        iRank=LocalColor % u2l_ListNeigh_send(iProc)
        call mpi_isend(LocalColor % ACexch,1,LocalColor%u2l_p2dsend_type(iProc),iRank-1,1151,comm,LocalColor%u2l_p2dsend_rqst(iProc),ierr)
      END DO
      IF (nbLow_send > 0) THEN
        call mpi_waitall(nbLow_send, LocalColor%u2l_p2dsend_rqst, LocalColor%u2l_p2dsend_stat,ierr)
      END IF
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5_EXCHANGE_P3_UPP_2_LOW_Recv(LocalColor, AC, iBlock)
      USE DATAPOOL, only : LocalColorInfo, MNP, MSC, MDC, rkind
      USE DATAPOOL, only : wwm_ListNeigh_recv, wwm_p2drecv_type
      USE elfe_msgp, only : comm, ierr, myrank
      implicit none
      include 'mpif.h'
# ifdef DEBUG
      real(rkind) :: Hfield(MNP)
# endif
      type(LocalColorInfo), intent(inout) :: LocalColor
      REAL(rkind), intent(inout) :: AC(MNP, MSC, MDC)
      INTEGER, intent(in) :: iBlock
      integer iProc, i, iRank, idx, lenBlock, maxBlockLength
      integer nbUpp_recv, IS, ID
      lenBlock=LocalColor % BlockLength(iBlock)
      DO idx=1,lenBlock
        IS=LocalColor % ISindex(iBlock, idx)
        ID=LocalColor % IDindex(iBlock, idx)
        LocalColor % ACexch(idx,:)=AC(:,IS,ID)
      END DO
# ifdef DEBUG
      write(myrank+740,*) 'I5_P3_UPP_2_LOW_Recv(I), iBlock=', iBlock
      DO idx=1,lenBlock
        IS=LocalColor % ISindex(iBlock, idx)
        ID=LocalColor % IDindex(iBlock, idx)
        Hfield(:)=AC(:,IS,ID)
        write(myrank+740,*) IS, ID, sum(Hfield)
      END DO
# endif
      nbUpp_recv=LocalColor % u2l_nnbr_recv
      DO iProc=1,nbUpp_recv
        iRank=LocalColor % u2l_ListNeigh_recv(iProc)
        call mpi_irecv(LocalColor % ACexch,1,LocalColor%u2l_p2drecv_type(iProc),iRank-1,1151,comm,LocalColor % u2l_p2drecv_rqst(iProc),ierr)
      END DO
      IF (nbUpp_recv > 0) THEN
        call mpi_waitall(nbUpp_recv, LocalColor%u2l_p2drecv_rqst, LocalColor%u2l_p2drecv_stat,ierr)
      END IF
      DO idx=1,lenBlock
        IS=LocalColor % ISindex(iBlock, idx)
        ID=LocalColor % IDindex(iBlock, idx)
        AC(:,IS,ID)=LocalColor % ACexch(idx,:)
      END DO
# ifdef DEBUG
      write(myrank+740,*) 'I5_P3_UPP_2_LOW_Recv(II), iBlock=', iBlock
      DO idx=1,lenBlock
        IS=LocalColor % ISindex(iBlock, idx)
        ID=LocalColor % IDindex(iBlock, idx)
        Hfield(:)=AC(:,IS,ID)
        write(myrank+740,*) IS, ID, sum(Hfield)
      END DO
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
# ifdef DEBUG
      SUBROUTINE COMPUTE_TOTAL_INDEX_SHIFT(TheRes)
      USE DATAPOOL, only : NP_RES, IA, JA
      implicit none
      integer, intent(out) :: TheRes
      integer :: IP, JP, J
      TheRes=0
      DO IP=1,NP_RES
        DO J=IA(IP),IA(IP+1)-1
          JP=JA(J)
          TheRes=TheRes+abs(IP-JP)
        END DO
      END DO
      END SUBROUTINE
# endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5_PARTIAL_SOLVE_L(LocalColor, SolDat, iBlock, ACret)
      USE DATAPOOL, only : LocalColorInfo, I5_SolutionData
      USE DATAPOOL, only : IA, JA, I_DIAG, MSC, MDC, MNP, rkind, NP_RES
      USE DATAPOOL, only : DO_SOLVE_L, DO_SOLVE_U
      USE elfe_msgp, only : myrank
      USE elfe_glbl, only : iplg
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      type(I5_SolutionData), intent(in) :: SolDat
      integer, intent(in) :: iBlock
# ifdef DEBUG
      real(rkind) :: Hfield(MNP)
# endif
      real(rkind) :: eCoeff
      real(rkind), intent(inout) :: ACret(MNP,MSC,MDC)
      integer IP, idx, ID, IS, J, IP_glob
      integer lenBlock, JP
      lenBlock=LocalColor % BlockLength(iBlock)
# ifdef DEBUG
      write(myrank+740,*) 'I5_PARTIAL_SOLVE_L, MNP=', MNP
      write(myrank+740,*) 'I5_PARTIAL_SOLVE_L, sum(CovLower)=', sum(LocalColor%CovLower)
      WRITE(myrank+790,*) 'eColor=', LocalColor%ListColor(myrank+1)
# endif
      DO IP=1,NP_RES
        IF (LocalColor % CovLower(IP) == 1) THEN
          DO J=IA(IP),IA(IP+1)-1
            IF (LocalColor % Jstatus_L(J) .eq. 1) THEN
              JP=JA(J)
# ifdef DEBUG
              WRITE(myrank+790,*) 'iplg IP/JP/J=', iplg(IP), iplg(JP), J
              WRITE(myrank+490,*) iplg(IP), iplg(JP)
# endif
              DO idx=1,lenBlock
                IS=LocalColor % ISindex(iBlock, idx)
                ID=LocalColor % IDindex(iBlock, idx)
                eCoeff=SolDat % ASPAR_pc(J,IS,ID)
# ifdef DEBUG
                IF ((IS.eq.2).and.(ID.eq.2)) THEN
                  WRITE(myrank+1490,*) iplg(IP), iplg(JP), eCoeff
                  WRITE(myrank+790,*) 'eCoeff, AC12=', eCoeff, ACret(IP,IS,ID), ACret(JP,IS,ID)
                END IF
# endif
                IF (DO_SOLVE_L) THEN
                  ACret(IP,IS,ID)=ACret(IP,IS,ID) - eCoeff*ACret(JP,IS,ID)
                END IF
              END DO
            END IF
          END DO
        ENDIF
      END DO
# ifdef DEBUG
      WRITE(myrank+790,*) 'ENDING'
      WRITE(myrank+490,*) 'ENDING'
      WRITE(myrank+1490,*) 'ENDING'
      write(myrank+740,*) 'I5_PARTIAL_SOLVE_L, iBlock=', iBlock
      write(myrank+740,*) 'MNP=', MNP, 'NP_RES=', NP_RES
      DO IP=1,MNP
        IP_glob=iplg(IP)
        write(myrank+740,*) 'IP(l/g/c)=', IP, IP_glob, LocalColor%CovLower(IP), 'AC=', ACret(IP,2,2)
      END DO
      DO idx=1,lenBlock
        IS=LocalColor % ISindex(iBlock, idx)
        ID=LocalColor % IDindex(iBlock, idx)
        Hfield(:)=ACret(:,IS,ID)
        write(myrank+740,*) IS, ID, sum(Hfield)
      END DO
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5_PARTIAL_SOLVE_U(LocalColor, SolDat, iBlock, ACret)
      USE DATAPOOL, only : LocalColorInfo, I5_SolutionData
      USE DATAPOOL, only : MNP, IA, JA, I_DIAG, MSC, MDC, rkind, NP_RES
      USE DATAPOOL, only : DO_SOLVE_L, DO_SOLVE_U
      USE elfe_msgp, only : myrank
      USE elfe_glbl, only : iplg
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      type(I5_SolutionData), intent(in) :: SolDat
      integer, intent(in) :: iBlock
# ifdef DEBUG
      real(rkind) :: Hfield(MNP)
# endif
      real(rkind) :: eCoeff
      integer lenBlock, IP, JP, idx, J, IS, ID, IP_glob
      integer DoOper
      real(rkind), intent(inout) :: ACret(MNP,MSC,MDC)
      lenBlock=LocalColor % BlockLength(iBlock)
      DO IP=NP_RES,1,-1
        IF (LocalColor % CovLower(IP) == 1) THEN
          DO J=IA(IP),IA(IP+1)-1
            IF (LocalColor % Jstatus_U(J) .eq. 1) THEN
              JP=JA(J)
#ifdef DEBUG
              WRITE(myrank+800,*) 'iplg(IP/JP)/J=', iplg(IP), iplg(JP), J
              WRITE(myrank+500,*) iplg(IP), iplg(JP)
#endif
              DO idx=1,lenBlock
                IS=LocalColor % ISindex(iBlock, idx)
                ID=LocalColor % IDindex(iBlock, idx)
                eCoeff=SolDat % ASPAR_pc(J,IS,ID)
#ifdef DEBUG
                IF ((IS.eq.2).and.(ID.eq.2)) THEN
                  WRITE(myrank+1500,*) iplg(IP), iplg(JP), eCoeff
                  WRITE(myrank+800,*) 'eCoeff, AC12=', eCoeff, ACret(IP,IS,ID), ACret(JP,IS,ID)
                END IF
#endif
                IF (DO_SOLVE_U) THEN
                  ACret(IP,IS,ID)=ACret(IP,IS,ID) - eCoeff*ACret(JP,IS,ID)
                END IF
              END DO
            END IF
          END DO
          J=I_DIAG(IP)
          DO idx=1,lenBlock
            IS=LocalColor % ISindex(iBlock, idx)
            ID=LocalColor % IDindex(iBlock, idx)
            IF (DO_SOLVE_U) THEN
              ACret(IP,IS,ID)=ACret(IP,IS,ID)*SolDat % ASPAR_pc(J,IS,ID)
            END IF
          END DO
        ENDIF
      END DO
# ifdef DEBUG
      WRITE(myrank+500,*) 'ENDING'
      WRITE(myrank+1500,*) 'ENDING'
      write(myrank+740,*) 'I5_PARTIAL_SOLVE_U, iBlock=', iBlock
      write(myrank+740,*) 'MNP=', MNP, 'NP_RES=', NP_RES
      DO IP=1,MNP
        IP_glob=iplg(IP)
        write(myrank+740,*) 'IP(l/g/c)=', IP, IP_glob, LocalColor % CovLower(IP), 'AC=', ACret(IP,2,2)
      END DO
      DO idx=1,lenBlock
        IS=LocalColor % ISindex(iBlock, idx)
        ID=LocalColor % IDindex(iBlock, idx)
        Hfield(:)=ACret(:,IS,ID)
        write(myrank+740,*) IS, ID, sum(Hfield)
      END DO
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5_SYNC_SENDRECV(LocalColor, AC)
      USE DATAPOOL, only : LocalColorInfo, MNP, MSC, MDC, rkind
      USE DATAPOOL, only : wwm_ListNeigh_send, wwm_ListNeigh_recv
      USE elfe_msgp, only : comm, ierr, myrank
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      real(rkind), intent(inout) :: AC(MNP, MSC, MDC)
      real(rkind) :: U(MSC,MDC,MNP)
      INTEGER :: IP, IS, ID, i, iRank
      integer iSync
      integer nbSync_send, nbSync_recv
      nbSync_recv=LocalColor%sync_nnbr_recv
      nbSync_send=LocalColor%sync_nnbr_send
# ifdef DEBUG
      write(myrank+740,*) 'nbSync_recv=', nbSync_recv
      write(myrank+740,*) 'nbSync_send=', nbSync_send
      CALL FLUSH(myrank+740)
# endif
      DO IP = 1, MNP
        DO IS = 1, MSC
          DO ID = 1, MDC
            U(IS,ID,IP) = AC(IP,IS,ID)
          END DO
        END DO
      END DO
      DO iSync=1,nbSync_send
        iRank=LocalColor % sync_ListNeigh_send(iSync)
# ifdef DEBUG
        write(myrank+740,*) 'SEND iSync=', iSync, ' iRank=', iRank
        CALL FLUSH(myrank+740)
# endif
        CALL mpi_isend(U, 1, LocalColor%sync_p2dsend_type(iSync), iRank-1, 1009, comm, LocalColor%sync_p2dsend_rqst(iSync), ierr)
      END DO
      DO iSync=1,nbSync_recv
        iRank=LocalColor % sync_ListNeigh_recv(iSync)
# ifdef DEBUG
        write(myrank+740,*) 'RECV iSync=', iSync, ' iRank=', iRank
        CALL FLUSH(myrank+740)
# endif
        call mpi_irecv(U,1,LocalColor%sync_p2drecv_type(iSync),iRank-1,1009,comm,LocalColor%sync_p2drecv_rqst(iSync),ierr)
      END DO
      IF (nbSync_send > 0) THEN
        call mpi_waitall(nbSync_send, LocalColor%sync_p2dsend_rqst, LocalColor%sync_p2dsend_stat,ierr)
      END IF
      IF (nbSync_recv > 0) THEN
        call mpi_waitall(nbSync_recv, LocalColor%sync_p2drecv_rqst, LocalColor%sync_p2drecv_stat,ierr)
      END IF
      DO IP = 1, MNP
        DO IS = 1, MSC
          DO ID = 1, MDC
            AC(IP,IS,ID) = U(IS,ID,IP)
          END DO
        END DO
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5_APPLY_PRECOND(LocalColor, SolDat, ACret)
      USE DATAPOOL, only : LocalColorInfo, I5_SolutionData
      USE DATAPOOL, only : MSC, MDC, MNP, rkind
      USE elfe_msgp, only : myrank
      USE DATAPOOL, only : DO_SYNC_UPP_2_LOW, DO_SYNC_LOW_2_UPP, DO_SYNC_FINAL
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      type(I5_SolutionData), intent(in) :: SolDat
      real(rkind), intent(inout) :: ACret(MNP, MSC, MDC)
      integer iBlock, lenBlock, idx, IS, ID
      integer maxBlockLength
      maxBlockLength=LocalColor % maxBlockLength
# ifdef DEBUG
      write(myrank+740,*) 'maxBlockLength=', maxBlockLength
      CALL FLUSH(myrank+740)
# endif
      DO iBlock=1,LocalColor % Nblock
# ifdef DEBUG
        write(myrank+740,*) 'L1: ACr max=', maxval(ACret), 'sum=', sum(ACret)
        CALL FLUSH(myrank+740)
        write(myrank+740,*) 'iBlock=', iBlock
        CALL FLUSH(myrank+740)
# endif
        IF (DO_SYNC_LOW_2_UPP) THEN
          CALL I5_EXCHANGE_P3_LOW_2_UPP_Recv(LocalColor, ACret, iBlock)
        END IF
# ifdef DEBUG
        write(myrank+740,*) 'L2: ACr max=', maxval(ACret), 'sum=', sum(ACret)
        CALL FLUSH(myrank+740)
        write(myrank+740,*) 'I5_EXCHANGE_P3_LOW_2_UOOER_Recv'
        CALL FLUSH(myrank+740)
# endif
        CALL I5_PARTIAL_SOLVE_L(LocalColor, SolDat, iBlock, ACret)
# ifdef DEBUG
        write(myrank+740,*) 'L3: ACr max=', maxval(ACret), 'sum=', sum(ACret)
        CALL FLUSH(myrank+740)
        write(myrank+740,*) 'I5_PARTIAL_SOLVE_L'
        CALL FLUSH(myrank+740)
# endif
        IF (DO_SYNC_LOW_2_UPP) THEN
          CALL I5_EXCHANGE_P3_LOW_2_UPP_Send(LocalColor, ACret, iBlock)
        END IF
# ifdef DEBUG
        write(myrank+740,*) 'L4: ACr max=', maxval(ACret), 'sum=', sum(ACret)
        CALL FLUSH(myrank+740)
        write(myrank+740,*) 'I5_EXCHANGE_P3_LOW_2_UPP_Send'
        CALL FLUSH(myrank+740)
# endif
      END DO
      DO iBlock=1,LocalColor%Nblock
# ifdef DEBUG
        write(myrank+740,*) 'U1: ACr max=', maxval(ACret), 'sum=', sum(ACret)
        CALL FLUSH(myrank+740)
        write(myrank+740,*) 'iBlock=', iBlock
        CALL FLUSH(myrank+740)
# endif
        IF (DO_SYNC_UPP_2_LOW) THEN
          CALL I5_EXCHANGE_P3_UPP_2_LOW_Recv(LocalColor, ACret, iBlock)
        END IF
# ifdef DEBUG
        write(myrank+740,*) 'U2: ACr max=', maxval(ACret), 'sum=', sum(ACret)
        CALL FLUSH(myrank+740)
        write(myrank+740,*) 'I5_EXCHANGE_P3_UPP_2_LOW_Recv'
        CALL FLUSH(myrank+740)
# endif
        CALL I5_PARTIAL_SOLVE_U(LocalColor, SolDat, iBlock, ACret)
# ifdef DEBUG
        write(myrank+740,*) 'U3: ACr max=', maxval(ACret), 'sum=', sum(ACret)
        CALL FLUSH(myrank+740)
        write(myrank+740,*) 'I5_PARTIAL_SOLVE_U'
        CALL FLUSH(myrank+740)
# endif
        IF (DO_SYNC_UPP_2_LOW) THEN
          CALL I5_EXCHANGE_P3_UPP_2_LOW_Send(LocalColor, ACret, iBlock)
        END IF
# ifdef DEBUG
        write(myrank+740,*) 'U4: ACr max=', maxval(ACret), 'sum=', sum(ACret)
        write(myrank+740,*) 'I5_EXCHANGE_P3_UPP_2_LOW_Send'
        CALL FLUSH(myrank+740)
# endif
      END DO
# ifdef DEBUG
      write(myrank+740,*) 'Before call to I5_SYNC_SENDRECV'
      CALL FLUSH(myrank+740)
# endif
      IF (DO_SYNC_FINAL) THEN
        CALL I5_SYNC_SENDRECV(LocalColor, ACret)
      END IF
# ifdef DEBUG
      write(myrank+740,*) 'End of APPLY_PRECOND'
      CALL FLUSH(myrank+740)
# endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5_APPLY_FCT(SolDat,  ACin, ACret)
      USE DATAPOOL, only : I5_SolutionData, IA, JA, NP_RES, MSC, MDC, MNP, rkind
      implicit none
      integer IP, J, idx
      type(I5_SolutionData), intent(inout) :: SolDat
      REAL(rkind) :: eSum(MSC,MDC)
      REAL(rkind), intent(in) :: ACin(MNP, MSC, MDC)
      REAL(rkind), intent(inout) :: ACret(MNP, MSC, MDC)
      DO IP=1,NP_RES
        eSum=0
        DO J=IA(IP),IA(IP+1)-1
          idx=JA(J)
          eSum=eSum + SolDat % ASPAR_block(J,:,:)*ACin(idx,:,:)
        END DO
        ACret(IP,:,:)=eSum
      END DO
      CALL EXCHANGE_P4D_WWM_TR(ACret)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE REPLACE_NAN_ZERO(LScal)
      USE DATAPOOL, only : rkind, MSC, MDC
      implicit none
      real(rkind), intent(inout) :: LScal(MSC,MDC)
      integer IS, ID
      DO IS=1,MSC
        DO ID=1,MDC
          IF (LScal(IS,ID) .ne. LScal(IS,ID)) THEN
            LScal(IS,ID)=0
          END IF
        END DO
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5_SCALAR(ACw1, ACw2, LScal)
      USE DATAPOOL, only : rkind, MNP, MSC, MDC
      USE DATAPOOL, only : nwild_loc_res, NP_RES
      USE elfe_msgp, only : myrank, comm, ierr, nproc, istatus, rtype
      implicit none
      real(rkind), intent(in) :: ACw1(MNP, MSC, MDC)
      real(rkind), intent(in) :: ACw2(MNP, MSC, MDC)
      real(rkind), intent(inout) :: LScal(MSC, MDC)
      real(rkind) :: RScal(MSC, MDC)
      integer IP, iProc
      LScal=0
      DO IP=1,NP_RES
        LScal=LScal + nwild_loc_res(IP)*ACw1(IP,:,:)*ACw2(IP,:,:)
      END DO
      IF (myrank == 0) THEN
        DO iProc=2,nproc
          CALL MPI_RECV(RScal,MSC*MDC,rtype, iProc-1, 19, comm, istatus, ierr)
          LScal = LScal + RScal
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(LScal,MSC*MDC,rtype, iProc-1, 23, comm, ierr)
        END DO
      ELSE
        CALL MPI_SEND(LScal,MSC*MDC,rtype, 0, 19, comm, ierr)
        CALL MPI_RECV(LScal,MSC*MDC,rtype, 0, 23, comm, istatus, ierr)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5_TOTAL_COHERENCY_ERROR(ACw, Lerror)
      USE DATAPOOL, only : MNP, MSC, MDC, rkind
      USE DATAPOOL, only : ListIPLG, ListMNP
      USE elfe_msgp, only : istatus, ierr, comm, rtype, myrank, nproc
      USE elfe_glbl, only : iplg, np_global
      implicit none
      real(rkind), intent(in) :: ACw(MNP, MSC, MDC)
      real(rkind), intent(out) :: Lerror
      real(rkind), allocatable :: ACtotal(:,:,:)
      real(rkind), allocatable :: ACloc(:,:,:)
      real(rkind), allocatable :: rbuf_real(:)
      integer, allocatable :: ListFirstMNP(:)
      integer, allocatable :: eStatus(:)
      integer IP, iProc, IPglob, IS, ID
      integer MNPloc
      IF (myrank == 0) THEN
        Lerror=0
        allocate(ListFirstMNP(nproc))
        ListFirstMNP=0
        DO iProc=2,nproc
          ListFirstMNP(iProc)=ListFirstMNP(iProc-1) + ListMNP(iProc-1)
        END DO
        allocate(eStatus(np_global))
        allocate(ACtotal(np_global, MSC, MDC))
        eStatus=0
        DO IP=1,MNP
          IPglob=iplg(IP)
          ACtotal(IPglob,:,:)=ACw(IP,:,:)
          eStatus(IPglob)=1
        END DO
        DO iProc=2,nproc
          MNPloc=ListMNP(iProc)
          allocate(ACloc(MNPloc, MSC, MDC))
          CALL MPI_RECV(ACloc,MNPloc*MSC*MDC,rtype, iProc-1, 53, comm, istatus, ierr)
          DO IP=1,MNPloc
            IPglob=ListIPLG(IP+ListFirstMNP(iProc))
            IF (eStatus(IPglob) == 1) THEN
              DO IS=1,MSC
                DO ID=1,MDC
                  Lerror=Lerror+abs(ACtotal(IPglob,IS,ID)-ACloc(IP,IS,ID))
                END DO
              END DO
            ELSE
              eStatus(IPglob)=1
              ACtotal(IPglob,:,:)=ACloc(IP,:,:)
            END IF
          END DO
          deallocate(ACloc)
        END DO
        deallocate(ListFirstMNP)
        deallocate(ACtotal)
        deallocate(eStatus)
        allocate(rbuf_real(1))
        rbuf_real(1)=Lerror
        DO iProc=2,nproc
          CALL MPI_SEND(rbuf_real,1,rtype, iProc-1, 23, comm, ierr)
        END DO
        deallocate(rbuf_real)
      ELSE
        CALL MPI_SEND(ACw,MNP*MSC*MDC,rtype, 0, 53, comm, ierr)
        allocate(rbuf_real(1))
        CALL MPI_RECV(rbuf_real,1,rtype, 0, 23, comm, istatus, ierr)
        Lerror=rbuf_real(1)
        deallocate(rbuf_real)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5_SUM_MAX(ACw, LSum, LMax)
      USE DATAPOOL, only : rkind, MNP, MSC, MDC
      USE DATAPOOL, only : nwild_loc_res, NP_RES
      USE elfe_msgp, only : myrank, comm, ierr, nproc, istatus, rtype
      implicit none
      real(rkind), intent(in) :: ACw(MNP, MSC, MDC)
      real(rkind), intent(inout) :: LSum(MSC, MDC)
      real(rkind), intent(inout) :: LMax(MSC, MDC)
      real(rkind) :: RScal(MSC, MDC)
      integer IP, iProc, IS, ID
      LSum=0
      DO IP=1,NP_RES
        LSum=LSum + nwild_loc_res(IP)*ACw(IP,:,:)
      END DO
      DO IS=1,MSC
        DO ID=1,MDC
          LMax(IS,ID)=maxval(ACw(:,IS,ID))
        END DO
      END DO
      IF (myrank == 0) THEN
        DO iProc=2,nproc
          CALL MPI_RECV(RScal,MSC*MDC,rtype, iProc-1, 53, comm, istatus, ierr)
          LSum = LSum + RScal
        END DO
        DO iProc=2,nproc
          CALL MPI_RECV(RScal,MSC*MDC,rtype, iProc-1, 59, comm, istatus, ierr)
          DO IS=1,MSC
            DO ID=1,MDC
              IF (RScal(IS,ID) .gt. LMax(IS,ID)) THEN
                LMax(IS,ID)=RScal(IS,ID)
              END IF
            END DO
          END DO
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(LSum,MSC*MDC,rtype, iProc-1, 197, comm, ierr)
        END DO
        DO iProc=2,nproc
          CALL MPI_SEND(LMax,MSC*MDC,rtype, iProc-1, 199, comm, ierr)
        END DO
      ELSE
        CALL MPI_SEND(LSum,MSC*MDC,rtype, 0, 53, comm, ierr)
        CALL MPI_SEND(LMax,MSC*MDC,rtype, 0, 59, comm, ierr)
        CALL MPI_RECV(LSum,MSC*MDC,rtype, 0, 197, comm, istatus, ierr)
        CALL MPI_RECV(LMax,MSC*MDC,rtype, 0, 199, comm, istatus, ierr)
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
! We use the notations of
! http://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method
! and we use K1=Id
! In this algorithm, the use of v_{i-1}, v_i can be replace to just "v"
! The same for x, r
! 
      SUBROUTINE I5_BCGS_SOLVER(LocalColor, SolDat)
      USE DATAPOOL, only : MSC, MDC, MNP, NP_RES, NNZ, AC2, SOLVERTHR
      USE DATAPOOL, only : LocalColorInfo, I5_SolutionData, rkind
      USE DATAPOOL, only : PCmethod
# ifdef DEBUG
      USE DATAPOOL, only : MNE, XP, YP, INE
      USE DATAPOOL, only : IA, JA, PCmethod
      USE elfe_msgp, only : myrank, nproc
      USE elfe_glbl, only : iplg, np_global
#  ifdef NCDF
      USE NETCDF
#  endif
# endif
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      type(I5_SolutionData), intent(inout) :: SolDat
      REAL(rkind) :: Rho(MSC,MDC)
      REAL(rkind) :: Prov(MSC,MDC)
      REAL(rkind) :: Alpha(MSC,MDC)
      REAL(rkind) :: Beta(MSC,MDC)
      REAL(rkind) :: Omega(MSC,MDC)
      REAL(rkind) :: MaxError, CritVal
      REAL(rkind) :: eSum1, eSum2
      REAL(rkind) :: TheTol
# ifdef DEBUG
      REAL(rkind) :: LSum(MSC,MDC), LMax(MSC,MDC)
# endif
      integer :: MaxIter = 30
      integer IP, IS, ID, nbIter
# ifdef DEBUG
      integer, SAVE :: iSystem = 1
      integer nbDiff
      real(rkind) :: Lerror
      integer idAC1, idAC2, idAC3, idAC4, idAC5, idAC6, idAC7, idAC8, idAC9
      integer idProv, idOmega, idRho, idAlpha, idBeta, idCritVal
      integer iret, ncid, var_id
      integer iter_dims, three_dims
      integer npgl_dims, mnp_dims, mnpp_dims, np_res_dims
      integer nnz_dims, msc_dims, mdc_dims, mne_dims
      character(len =256) :: FILE_NAME, PRE_FILE_NAME
      integer :: ListPos(np_global)
      integer ISsel, IDsel
      integer J, JP, IPglob, JPglob, IPos, JPos
      integer :: ListPosB(np_global)
      integer :: ListPosBRev(np_global)
      ISsel=2
      IDsel=2
# endif
# ifdef DEBUG
      DO IS=1,MSC
        DO ID=1,MDC
          IF ((IS.ne.ISsel).or.(ID.ne.IDsel)) THEN
            AC2(:,IS,ID)=0.0_rkind
          END IF
        END DO
      END DO
      WRITE(myrank+240,*) 'MNP=', MNP, 'NP_RES=', NP_RES
# endif
      MaxError=SOLVERTHR
      CALL I5_APPLY_FCT(SolDat,  AC2, SolDat % AC3)
      SolDat % AC1=0                               ! y
      SolDat % AC2=AC2                             ! x solution
# ifdef DEBUG
      CALL I5_SUM(SolDat%AC2, eSum1)
      CALL I5_SUM(AC2, eSum2)
      CALL I5_LOCATE_MAX(AC2)
      WRITE(myrank+240,*) 'AC2 max=', maxval(SolDat%AC2), 'sum=', sum(SolDat % AC2)
      WRITE(myrank+240,*) 'eSum1=', eSum1, ' eSum2=', eSum2
      WRITE(myrank+240,*) 'AC3 max=', maxval(SolDat%AC3), 'sum=', sum(SolDat % AC3)
# endif
      SolDat % AC3=SolDat % B_block - SolDat % AC3 ! r residual
      SolDat % AC4=SolDat % AC3                    ! hat{r_0} term
      SolDat % AC5=0                               ! v
      SolDat % AC6=0                               ! p
      SolDat % AC7=0                               ! s
      SolDat % AC8=0                               ! z
      SolDat % AC9=0                               ! t
      Rho=1
      Alpha=1
      Omega=1
      nbIter=0
# if defined DEBUG && defined NCDF
      PRE_FILE_NAME='DebugAC'
      WRITE (FILE_NAME,10) TRIM(PRE_FILE_NAME),nproc, iSystem, myrank
  10  FORMAT (a,'_np',i2.2,'_',i3.3,'_',i4.4, '.nc')
      iret = nf90_create(TRIM(FILE_NAME), NF90_CLOBBER, ncid)
      iret = nf90_def_dim(ncid, 'iter', NF90_UNLIMITED, iter_dims)
      iret = nf90_def_dim(ncid, 'three', 3, three_dims)
      iret = nf90_def_dim(ncid, 'msc', MSC, msc_dims)
      iret = nf90_def_dim(ncid, 'mdc', MDC, mdc_dims)
      iret = nf90_def_dim(ncid, 'mnp', MNP, mnp_dims)
      iret = nf90_def_dim(ncid, 'mnpp', MNP+1, mnpp_dims)
      iret = nf90_def_dim(ncid, 'np_global', np_global, npgl_dims)
      iret = nf90_def_dim(ncid, 'np_res', NP_RES, np_res_dims)
      iret = nf90_def_dim(ncid, 'mne', MNE, mne_dims)
      iret = nf90_def_dim(ncid, 'nnz', NNZ, nnz_dims)
      iret=nf90_def_var(ncid,'ASPAR',NF90_DOUBLE,(/nnz_dims,msc_dims, mdc_dims/),var_id)
      iret=nf90_def_var(ncid,'IA',NF90_INT,(/mnpp_dims/),var_id)
      iret=nf90_def_var(ncid,'JA',NF90_INT,(/nnz_dims/),var_id)
      iret=nf90_def_var(ncid,'ListPos',NF90_INT,(/npgl_dims/),var_id)
      iret=nf90_def_var(ncid,'AC1',NF90_DOUBLE,(/mnp_dims,msc_dims, mdc_dims, iter_dims /),var_id)
      iret=nf90_def_var(ncid,'AC2',NF90_DOUBLE,(/mnp_dims,msc_dims, mdc_dims, iter_dims /),var_id)
      iret=nf90_def_var(ncid,'AC3',NF90_DOUBLE,(/mnp_dims,msc_dims, mdc_dims, iter_dims /),var_id)
      iret=nf90_def_var(ncid,'AC4',NF90_DOUBLE,(/mnp_dims,msc_dims, mdc_dims, iter_dims /),var_id)
      iret=nf90_def_var(ncid,'AC5',NF90_DOUBLE,(/mnp_dims,msc_dims, mdc_dims, iter_dims /),var_id)
      iret=nf90_def_var(ncid,'AC6',NF90_DOUBLE,(/mnp_dims,msc_dims, mdc_dims, iter_dims /),var_id)
      iret=nf90_def_var(ncid,'AC7',NF90_DOUBLE,(/mnp_dims,msc_dims, mdc_dims, iter_dims /),var_id)
      iret=nf90_def_var(ncid,'AC8',NF90_DOUBLE,(/mnp_dims,msc_dims, mdc_dims, iter_dims /),var_id)
      iret=nf90_def_var(ncid,'AC9',NF90_DOUBLE,(/mnp_dims,msc_dims, mdc_dims, iter_dims /),var_id)
      !
      iret=nf90_def_var(ncid,'Prov',NF90_DOUBLE,(/msc_dims, mdc_dims, iter_dims /),var_id)
      iret=nf90_def_var(ncid,'Alpha',NF90_DOUBLE,(/msc_dims, mdc_dims, iter_dims /),var_id)
      iret=nf90_def_var(ncid,'Beta',NF90_DOUBLE,(/msc_dims, mdc_dims, iter_dims /),var_id)
      iret=nf90_def_var(ncid,'Rho',NF90_DOUBLE,(/msc_dims, mdc_dims, iter_dims /),var_id)
      iret=nf90_def_var(ncid,'Omega',NF90_DOUBLE,(/msc_dims, mdc_dims, iter_dims /),var_id)
      iret=nf90_def_var(ncid,'CritVal',NF90_DOUBLE,(/msc_dims, mdc_dims, iter_dims /),var_id)
      iret = nf90_def_var(ncid,'iplg',NF90_INT,(/ mnp_dims/),var_id)
      iret=nf90_def_var(ncid,'XP',NF90_DOUBLE,(/ mnp_dims/),var_id)
      iret=nf90_def_var(ncid,'YP',NF90_DOUBLE,(/ mnp_dims/),var_id)
      iret=nf90_def_var(ncid,'ine',NF90_INT,(/ three_dims, mne_dims/),var_id)
      iret=nf90_close(ncid)
      idAC1=0
      idAC2=0
      idAC3=0
      idAC4=0
      idAC5=0
      idAC6=0
      idAC7=0
      idAC8=0
      idAC9=0
      idProv=0
      idAlpha=0
      idBeta=0
      idRho=0
      idOmega=0
      idCritVal=0
      iret = nf90_open(TRIM(FILE_NAME), NF90_WRITE, ncid)
      iret=nf90_inq_varid(ncid, 'iplg', var_id)
      iret=nf90_put_var(ncid,var_id,iplg,start=(/1/), count = (/ MNP /))
      !
      iret=nf90_inq_varid(ncid, 'XP', var_id)
      iret=nf90_put_var(ncid,var_id,XP,start=(/1/), count = (/ MNP /))
      !
      iret=nf90_inq_varid(ncid, 'YP', var_id)
      iret=nf90_put_var(ncid,var_id,YP,start=(/1/), count = (/ MNP /))
      !
      iret=nf90_inq_varid(ncid, 'ine', var_id)
      iret=nf90_put_var(ncid,var_id,ine,start=(/1,1/), count = (/ 3, MNE /))
      !
      idAC2=idAC2+1
      iret=nf90_inq_varid(ncid, 'AC2', var_id)
      iret=nf90_put_var(ncid,var_id,SolDat%AC2,start=(/1,1,1,idAC2/), count = (/ MNP, MSC, MDC, 1 /))
      idAC3=idAC3+1
      iret=nf90_inq_varid(ncid, 'AC3', var_id)
      iret=nf90_put_var(ncid,var_id,SolDat%AC3,start=(/1,1,1,idAC3/), count = (/ MNP, MSC, MDC, 1 /))
      idAC4=idAC4+1
      iret=nf90_inq_varid(ncid, 'AC4', var_id)
      iret=nf90_put_var(ncid,var_id,SolDat%AC4,start=(/1,1,1,idAC4/), count = (/ MNP, MSC, MDC, 1 /))
      iret=nf90_inq_varid(ncid, 'ASPAR', var_id)
      iret=nf90_put_var(ncid,var_id,SolDat%aspar_block,start=(/1,1,1/), count = (/ NNZ, MSC, MDC/))

      iret=nf90_inq_varid(ncid, 'IA', var_id)
      iret=nf90_put_var(ncid,var_id,IA,start=(/1/), count = (/ MNP+1/))
      iret=nf90_inq_varid(ncid, 'JA', var_id)
      iret=nf90_put_var(ncid,var_id,JA,start=(/1/), count = (/ NNZ/))
      iret=nf90_inq_varid(ncid, 'ListPos', var_id)
      iret=nf90_put_var(ncid,var_id,ListPos,start=(/1/), count = (/ np_global/))
# endif
      DO
        nbIter=nbIter+1
# ifdef DEBUG
        IF (myrank .eq. 0) THEN
          Print *, 'nbIter=', nbIter
        END IF
        WRITE(myrank+240,*) 'nbIter=', nbIter
# endif

        ! L1: Rhoi =(\hat{r}_0, r_{i-1}
        CALL I5_SCALAR(SolDat % AC4, SolDat % AC3, Prov)
# if defined DEBUG && defined NCDF
        WRITE(myrank+240,*) 'Prov max=', maxval(Prov), 'sum=', sum(Prov)
        idProv=idProv+1
        iret=nf90_inq_varid(ncid, 'Prov', var_id)
        iret=nf90_put_var(ncid,var_id,Prov,start=(/1,1,idProv/), count = (/ MSC, MDC, 1 /))
# endif

        ! L2: Beta=(RhoI/Rho(I-1))  *  (Alpha/Omega(i-1))
        Beta=(Prov/Rho)*(Alpha/Omega)
# if defined DEBUG && defined NCDF
        WRITE(myrank+240,*) 'Beta max=', maxval(Beta), 'sum=', sum(Beta)
        idBeta=idBeta+1
        iret=nf90_inq_varid(ncid, 'Beta', var_id)
        iret=nf90_put_var(ncid,var_id,Beta,start=(/1,1,idBeta/), count = (/ MSC, MDC, 1 /))
# endif
        CALL REPLACE_NAN_ZERO(Beta)
        Rho=Prov
# if defined DEBUG && defined NCDF
        idRho=idRho+1
        iret=nf90_inq_varid(ncid, 'Rho', var_id)
        iret=nf90_put_var(ncid,var_id,Rho,start=(/1,1,idBeta/), count = (/ MSC, MDC, 1 /))
# endif

        ! L3: Pi = r(i-1) + Beta*(p(i-1) -omega(i-1)*v(i-1))
        DO IP=1,MNP
          SolDat%AC6(IP,:,:)=SolDat%AC3(IP,:,:)                        &
     &      + Beta(:,:)*SolDat%AC6(IP,:,:)                            &
     &      - Beta(:,:)*Omega(:,:)*SolDat%AC5(IP,:,:)
        END DO
# if defined DEBUG && defined NCDF
        WRITE(myrank+240,*) 'AC6(p) max=', maxval(SolDat % AC6), 'sum=', sum(SolDat % AC6)
        CALL I5_TOTAL_COHERENCY_ERROR(SolDat%AC6, Lerror)
        CALL I5_SUM_MAX(SolDat%AC6, LSum, LMax)
        IF (myrank .eq. 0) THEN
          Print *, 'AC6(coherr)=', Lerror
          Print *, 'AC6(sum/max)=', LSum(ISsel, IDsel), LMax(ISsel, IDsel)
        END IF
        idAC6=idAC6+1
        iret=nf90_inq_varid(ncid, 'AC6', var_id)
        iret=nf90_put_var(ncid,var_id,SolDat%AC6,start=(/1,1,1,idAC6/), count = (/ MNP, MSC, MDC, 1 /))
# endif

        ! L4 y=K^(-1) Pi
        SolDat%AC1=SolDat%AC6
# ifdef DEBUG
        WRITE(myrank+240,*) 'Before apply precond'
# endif
        IF (PCmethod .gt. 0) THEN
          CALL I5_APPLY_PRECOND(LocalColor, SolDat, SolDat%AC1)
        ENDIF
# if defined DEBUG && defined NCDF
        WRITE(myrank+240,*) 'AC1(y) max=', maxval(SolDat % AC1), 'sum=', sum(SolDat % AC1)
        CALL I5_TOTAL_COHERENCY_ERROR(SolDat%AC1, Lerror)
        CALL I5_SUM_MAX(SolDat%AC1, LSum, LMax)
        IF (myrank .eq. 0) THEN
          Print *, 'AC1(coherr)=', Lerror
          Print *, 'AC1(sum/max)=', LSum(ISsel, IDsel), LMax(ISsel, IDsel)
        END IF
        idAC1=idAC1+1
        iret=nf90_inq_varid(ncid, 'AC1', var_id)
        iret=nf90_put_var(ncid,var_id,SolDat%AC1,start=(/1,1,1,idAC1/), count = (/ MNP, MSC, MDC, 1 /))
        WRITE(myrank+240,*) 'idAC1=', idAC1, ' idAC6=', idAC6
# endif

        ! L5 vi=Ay
        CALL I5_APPLY_FCT(SolDat,  SolDat%AC1, SolDat%AC5)
# if defined DEBUG && defined NCDF
        WRITE(myrank+240,*) 'AC5(v) max=', maxval(SolDat % AC5), 'sum=', sum(SolDat % AC5)
        idAC5=idAC5+1
        iret=nf90_inq_varid(ncid, 'AC5', var_id)
        iret=nf90_put_var(ncid,var_id,SolDat%AC5,start=(/1,1,1,idAC5/), count = (/ MNP, MSC, MDC, 1 /))
# endif

        ! L6 Alpha=Rho/(hat(r)_0, v_i)
        CALL I5_SCALAR(SolDat % AC4, SolDat % AC5, Prov)
        Alpha(:,:)=Rho(:,:)/Prov(:,:)
        CALL REPLACE_NAN_ZERO(Alpha)
# if defined DEBUG && defined NCDF
        WRITE(myrank+240,*) 'Alpha max=', maxval(Alpha), 'sum=', sum(Alpha)
        idAlpha=idAlpha+1
        iret=nf90_inq_varid(ncid, 'Alpha', var_id)
        iret=nf90_put_var(ncid,var_id,Alpha,start=(/1,1,idAlpha/), count = (/ MSC, MDC, 1 /))
# endif

        ! L7 s=r(i-1) - alpha v(i)
        DO IP=1,MNP
          SolDat%AC7(IP,:,:)=SolDat%AC3(IP,:,:)                        &
     &      - Alpha(:,:)*SolDat%AC5(IP,:,:)
        END DO
# if defined DEBUG && defined NCDF
        WRITE(myrank+240,*) 'AC7(s) max=', maxval(SolDat%AC7), 'sum=', sum(SolDat%AC7)
        idAC7=idAC7+1
        iret=nf90_inq_varid(ncid, 'AC7', var_id)
        iret=nf90_put_var(ncid,var_id,SolDat%AC7,start=(/1,1,1,idAC7/), count = (/ MNP, MSC, MDC, 1 /))
# endif

        ! L8 z=K^(-1) s
        SolDat%AC8=SolDat%AC7
        IF (PCmethod .gt. 0) THEN
          CALL I5_APPLY_PRECOND(LocalColor, SolDat, SolDat%AC8)
        END IF
# ifdef DEBUG
        WRITE(myrank+240,*) 'max(AC8-7)=', maxval(SolDat%AC8 - SolDat%AC7)
        WRITE(myrank+240,*) 'min(AC8-7)=', minval(SolDat%AC8 - SolDat%AC7)
        WRITE(myrank+240,*) 'SYNCERR'
        TheTol=0.1_rkind
        nbDiff=0
        DO IP=1,MNP
          IF (abs(SolDat%AC8(IP,2,2) - SolDat%AC7(IP,2,2)) .gt. TheTol) THEN
            WRITE(240+myrank,*) 'IPlg=', IP, iplg(IP), '87=', SolDat%AC8(IP,2,2), SolDat%AC7(IP,2,2)
            nbDiff=nbDiff+1
          END IF
        END DO
        WRITE(240+myrank,*) 'nbDiff=', nbDiff
        WRITE(240+myrank,*) 'AC8(z) max=', maxval(SolDat%AC8), 'sum=', sum(SolDat%AC8)
        idAC8=idAC8+1
        iret=nf90_inq_varid(ncid, 'AC8', var_id)
        iret=nf90_put_var(ncid,var_id,SolDat%AC8,start=(/1,1,1,idAC8/), count = (/ MNP, MSC, MDC, 1 /))
# endif

        ! L9 t=Az
        CALL I5_APPLY_FCT(SolDat,  SolDat%AC8, SolDat%AC9)
# if defined DEBUG && defined NCDF
        WRITE(myrank+240,*) 'AC9(t) max=', maxval(SolDat%AC9), 'sum=', sum(SolDat%AC9)
        idAC9=idAC9+1
        iret=nf90_inq_varid(ncid, 'AC9', var_id)
        iret=nf90_put_var(ncid,var_id,SolDat%AC9,start=(/1,1,1,idAC9/), count = (/ MNP, MSC, MDC, 1 /))
# endif

        ! L10 omega=(t,s)/(t,t)
        CALL I5_SCALAR(SolDat % AC9, SolDat % AC7, Omega)
        CALL I5_SCALAR(SolDat % AC9, SolDat % AC9, Prov)
        Omega(:,:)=Omega(:,:)/Prov(:,:)
        CALL REPLACE_NAN_ZERO(Omega)
# if defined DEBUG && defined NCDF
        WRITE(myrank+240,*) 'Omega max=', maxval(Omega), 'sum=', sum(Omega)
        idOmega=idOmega+1
        iret=nf90_inq_varid(ncid, 'Omega', var_id)
        iret=nf90_put_var(ncid,var_id,Omega,start=(/1,1,idOmega/), count = (/ MSC, MDC, 1 /))
# endif

        ! L11 x(i)=x(i-1) + Alpha y + Omega z
        DO IP=1,MNP
          SolDat%AC2(IP,:,:)=SolDat%AC2(IP,:,:)                        &
     &      + Alpha(:,:)*SolDat%AC1(IP,:,:)                            &
     &      + Omega(:,:)*SolDat%AC8(IP,:,:)
        END DO
# if defined DEBUG && defined NCDF
        WRITE(myrank+240,*) 'AC2(x) max=', maxval(SolDat%AC2), 'sum=', sum(SolDat%AC2)
        idAC2=idAC2+1
        iret=nf90_inq_varid(ncid, 'AC2', var_id)
        iret=nf90_put_var(ncid,var_id,SolDat%AC2,start=(/1,1,1,idAC2/), count = (/ MNP, MSC, MDC, 1 /))
# endif

        ! L12 If x is accurate enough finish
        CALL I5_APPLY_FCT(SolDat,  SolDat%AC2, SolDat%AC1)
        SolDat%AC8=SolDat%AC1 - SolDat%B_block
        CALL I5_SCALAR(SolDat % AC8, SolDat % AC8, Prov)
        CritVal=maxval(Prov)
# if defined DEBUG && defined NCDF
        WRITE(myrank+240,*) 'CritVal=', CritVal
        idCritVal=idCritVal+1
        iret=nf90_inq_varid(ncid, 'CritVal', var_id)
        iret=nf90_put_var(ncid,var_id,Prov,start=(/1,1,idCritVal/), count = (/ MSC, MDC, 1 /))
# endif

        IF (maxval(Prov) .lt. MaxError) THEN
          EXIT
        ENDIF
        IF (nbIter .gt. MaxIter) THEN
          EXIT
        ENDIF

        ! L13 r=s-omega t
        DO IP=1,MNP
          SolDat%AC3(IP,:,:)=SolDat%AC7(IP,:,:)                        &
     &      - Omega(:,:)*SolDat%AC9(IP,:,:)
        END DO
# if defined DEBUG && defined NCDF
        WRITE(myrank+240,*) 'AC3(r) max=', maxval(SolDat%AC3), 'sum=', sum(SolDat%AC3)
        idAC3=idAC3+1
        iret=nf90_inq_varid(ncid, 'AC3', var_id)
        iret=nf90_put_var(ncid,var_id,SolDat%AC3,start=(/1,1,1,idAC3/), count = (/ MNP, MSC, MDC, 1 /))
# endif
      END DO
# if defined DEBUG
      CALL WRITE_EXPLICIT_ORDERING(ListPosB, ListPosBRev, LocalColor%ListColor)
      DO IP=1,NP_RES
        IF (LocalColor % CovLower(IP) == 1) THEN
          DO J=IA(IP),IA(IP+1)-1
            JP=JA(J)
            IPglob=iplg(IP)
            JPglob=iplg(JP)
            IPos=ListPosBRev(IPglob)
            JPos=ListPosBRev(JPglob)
            IF (JPos .lt. IPos) THEN
              IF (LocalColor%Jstatus_L(J) .ne. 1) THEN
                WRITE(240+myrank,*) 'Begin Pair'
                WRITE(240+myrank,*) 'IP=', IP, ' IPglob=', IPglob, 'IPos=', IPos
                WRITE(240+myrank,*) 'JP=', JP, ' JPglob=', JPglob, 'JPos=', JPos
                WRITE(240+myrank,*) 'End Pair'
                CALL FLUSH(240+myrank)
                CALL WWM_ABORT('Major ordering inconsistency here 1')
              END IF
              IF (LocalColor%Jstatus_U(J) .ne. 0) THEN
                CALL WWM_ABORT('Major ordering inconsistency here 2')
              END IF
            END IF
            IF (JPos .gt. IPos) THEN
              IF (LocalColor%Jstatus_U(J) .ne. 1) THEN
                CALL WWM_ABORT('Major ordering inconsistency here 1')
              END IF
              IF (LocalColor%Jstatus_L(J) .ne. 0) THEN
                CALL WWM_ABORT('Major ordering inconsistency here 2')
              END IF
            END IF
          END DO
        END IF
      END DO
# endif
# if defined DEBUG && defined NCDF
      iret=nf90_close(ncid)
      IF (myrank == 0) THEN
        CALL SINGLE_PROC_SOLVE_CHECK(LocalColor, iSystem, nbiter-1, ISsel, IDsel)
      END IF
      iSystem=iSystem+1
# endif
      AC2=SolDat%AC2
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
# if defined DEBUG && defined NCDF
      SUBROUTINE SINGLE_PROC_SOLVE_CHECK(LocalColor, iSystem, nbEnt, IS, ID)
      USE DATAPOOL, only : LocalColorInfo
      USE DATAPOOL, only : MNP, ListIA, ListJA, ListIPLG, rkind
      USE DATAPOOL, only : ListMNP, ListNP_RES, ListNNZ
      USE DATAPOOL, only : NNZ, IA, JA, I_DIAG
      USE DATAPOOL, only : DO_SOLVE_L, DO_SOLVE_U
      USE elfe_msgp, only : nproc, myrank
      USE elfe_glbl, only : np_global
      USE NETCDF
      implicit none
      type(LocalColorInfo), intent(in) :: LocalColor
      integer, intent(in) :: iSystem, nbEnt, IS, ID
      integer idAC
      real(rkind) :: AC1glob(np_global)
      real(rkind) :: AC6glob(np_global)
      real(rkind), allocatable :: ASPARglob(:)
      integer,     allocatable :: JAglob(:)
      integer :: ListNBglob(np_global)
      integer :: IAglob(np_global+1)
      integer :: I_DIAGglob(np_global)
      integer :: ListPos(np_global)
      integer :: ListPosRev(np_global)
      integer IPglob, ShiftIA, ShiftMNP
      integer IP, IP_b, JP, JP_b
      integer IA1, IA2, len, J, idx, JPglob
      integer NNZglob, NNZloc, MNPloc, NP_RESloc
      character(len =256) :: FILE_NAME, PRE_FILE_NAME
      integer iret, ncid, var_id, iProc
      integer ListFirstMNP(nproc), ListFirstNNZ(nproc)
      real(rkind) :: eCoeff, eCoeffInv
      real(rkind) :: ACw(np_global)
      real(rkind) :: TotalErr, Sum1, Sum2
      real(rkind), allocatable :: AC6loc(:), AC1loc(:), ASPARloc(:)
      real(rkind) :: TheTol
      ListFirstNNZ=0
      ListFirstMNP=0
      DO iProc=2,nproc
        ListFirstMNP(iProc)=ListFirstMNP(iProc-1) + ListMNP(iProc-1)
        ListFirstNNZ(iProc)=ListFirstNNZ(iProc-1) + ListNNZ(iProc-1)
      END DO
      CALL WRITE_EXPLICIT_ORDERING(ListPos, ListPosRev, LocalColor % ListColor)
      DO IP=1,NP_GLOBAL
        Write(640,*) 'IP/Pos/PosR=', IP, ListPos(IP), ListPosRev(IP)
      END DO
      ListNBglob=-1
      DO iProc=1,nproc
        MNPloc=ListMNP(iProc)
        NP_RESloc=ListNP_RES(iProc)
        DO IP=1,NP_RESloc
          ShiftMNP=ListFirstMNP(iProc)
          ShiftIA=ListFirstMNP(iProc) + iProc - 1
          IPglob=ListIPLG(IP + ShiftMNP)
          IA1=ListIA(IP+ShiftIA)
          IA2=ListIA(IP+1+ShiftIA)
          len=IA2-IA1
          IF (ListNBglob(IPglob) == -1) THEN
            ListNBglob(IPglob)=len
          ELSE
            IF (ListNBglob(IPglob) .ne. len) THEN
              CALL WWM_ABORT('Error in length assignment')
            END IF
          END IF
        ENDDO
      END DO
      NNZglob=sum(ListNBglob)
      Write(640,*) 'NNZ=', NNZ
      Write(640,*) 'NNZglob=', NNZglob
      allocate(JAglob(NNZglob))
      allocate(ASPARglob(NNZglob))
      IAglob(1)=1
      DO IP=1,np_global
        IAglob(IP+1)=IAglob(IP) + ListNBglob(IP)
      END DO
      DO iProc=1,nproc
        NP_RESloc=ListNP_RES(iProc)
        DO IP=1,NP_RESloc
          ShiftMNP=ListFirstMNP(iProc)
          ShiftIA=ListFirstMNP(iProc) + iProc - 1
          IPglob=ListIPLG(IP + ShiftMNP)
          IA1=ListIA(IP+ShiftIA)
          IA2=ListIA(IP+1+ShiftIA)
          DO J=IA1,IA2-1
            idx=J-IA1
            JP=ListJA(J+ListFirstNNZ(iProc))
            JPglob=ListIPLG(JP+ShiftMNP)
            JAglob(idx+IAglob(IPglob))=JPglob
            Write(640,*) 'IP/JPglob, J=', IPglob, JPglob, idx+IAglob(IPglob)
          END DO
        END DO
      END DO
      I_DIAGglob=0
      DO IP=1,NP_GLOBAL
        Write(640,*) 'IP/IA1/IA2=', IP, IAglob(IP),IAglob(IP+1)
        DO J=IAglob(IP),IAglob(IP+1)-1
          JP=JAglob(J)
          Write(640,*) 'IP/JP, J=', IP, JP, J
          IF (IP .eq. JP) THEN
            I_DIAGglob(IP)=J
          END IF
        END DO
        IF (I_DIAGglob(IP) .eq. 0) THEN
          Print *, 'IP=', IP
          CALL WWM_ABORT('We found value 0 for I_DIAGglob')
        ENDIF
        Write(640,*) 'IP, I_DIAG=', IP, I_DIAGglob(IP)
      END DO
      IF (nproc .eq. 1) THEN
        WRITE(290,*) 'diff(I_DIAG)=', sum(abs(I_DIAG - I_DIAGglob))
        WRITE(290,*) 'diff(IA)=', sum(abs(IA - IAglob))
        WRITE(290,*) 'diff(JA)=', sum(abs(JA - JAglob))
      END IF
      DO idAC=1,nbEnt
        WRITE(2407,*) 'idAC=', idAC
        DO iProc=1,nproc
          PRE_FILE_NAME='DebugAC'
          WRITE (FILE_NAME,10) TRIM(PRE_FILE_NAME),nproc, iSystem, iProc-1
  10      FORMAT (a,'_np',i2.2,'_',i3.3,'_',i4.4, '.nc')
          NNZloc=ListNNZ(iProc)
          MNPloc=ListMNP(iProc)
          NP_RESloc=ListNP_RES(iProc)
          allocate(ASPARloc(NNZloc))
          allocate(AC1loc(MNPloc))
          allocate(AC6loc(MNPloc))
          iret=nf90_open(TRIM(FILE_NAME), NF90_NOWRITE, ncid)
          iret=nf90_inq_varid(ncid, "ASPAR", var_id)
          iret=nf90_get_var(ncid, var_id, ASPARloc, start=(/1,IS,ID/), count=(/NNZloc, 1, 1/) )
          iret=nf90_inq_varid(ncid, "AC1", var_id)
          iret=nf90_get_var(ncid, var_id, AC1loc, start=(/1,IS,ID,idAC/), count=(/MNPloc, 1, 1, 1/) )
          iret=nf90_inq_varid(ncid, "AC6", var_id)
          iret=nf90_get_var(ncid, var_id, AC6loc, start=(/1,IS,ID,idAC/), count=(/MNPloc, 1, 1, 1/) )
          DO IP=1,MNPloc
            IPglob=ListIPLG(IP+ListFirstMNP(iProc))
            AC1glob(IPglob)=AC1loc(IP)
            AC6glob(IPglob)=AC6loc(IP)
          END DO
          deallocate(AC1loc)
          deallocate(AC6loc)
          DO IP=1,NP_RESloc
            ShiftMNP=ListFirstMNP(iProc)
            ShiftIA=ListFirstMNP(iProc) + iProc - 1
            IPglob=ListIPLG(IP + ShiftMNP)
            IA1=ListIA(IP+ShiftIA)
            IA2=ListIA(IP+1+ShiftIA)
            DO J=IA1,IA2-1
              idx=J-IA1
              ASPARglob(idx+IAglob(IPglob))=ASPARloc(J)
            END DO
          END DO
          deallocate(ASPARloc)
          iret=nf90_close(ncid)
        END DO
        DO IP=1,np_global
          DO J=IAglob(IP),IAglob(IP+1)-1
            JP=JAglob(J)
            eCoeff=ASPARglob(J)
            WRITE(199,*) IP, JP, J, eCoeff
          END DO
        END DO
        DO IP=1,np_global
          J=I_DIAGglob(IP)
          eCoeff=ASPARglob(J)
          eCoeffInv=1.0_rkind/eCoeff
          WRITE(299,*) IP, eCoeff, eCoeffInv
        END DO
        ACw=AC6glob
        WRITE(2407,*) 'Before solver'
        WRITE(2407,*) 'ACw(sum/max)=', sum(ACw), maxval(ACw)
        Print *, 'Before solver'
        Print *, 'ACw(sum/max)=', sum(ACw), maxval(ACw)
        DO IP=1,NP_GLOBAL
          IP_b=ListPos(IP)
          DO J=IAglob(IP_b),IAglob(IP_b+1)-1
            JP_b=JAglob(J)
            JP=ListPosRev(JP_b)
            IF (JP .lt. IP) THEN
              WRITE(640,*) 'J=', J
              WRITE(640,*) 'JP_b=', JP_b
              WRITE(640,*) 'I_DIAGglob=', I_DIAGglob(JP_b)
              eCoeff=ASPARglob(J)
              WRITE(199,*) IP_b, JP_b, eCoeff
              eCoeff=ASPARglob(J)/ASPARglob(I_DIAGglob(JP_b))
!              WRITE(890,*) 'IP/JP=', IP, JP
              WRITE(890,*) 'IP_b/JP_b/J=', IP_b, JP_b, J
              WRITE(390,*) IP_b, JP_b
              WRITE(1390,*) IP_b, JP_b, eCoeff
              WRITE(890,*) 'AS(gl/di/co)=', ASPARglob(J), ASPARglob(I_DIAGglob(JP_b)), eCoeff
              WRITE(890,*) 'eCoeff, AC12=', eCoeff, ACw(IP_b), ACw(JP_b)
              IF (DO_SOLVE_L) THEN
                ACw(IP_b)=ACw(IP_b) - eCoeff*ACw(JP_b)
              END IF
            END IF
          END DO
        END DO
        WRITE(2407,*) 'After solve_L'
        WRITE(2407,*) 'ACw(sum/max)=', sum(ACw), maxval(ACw)
        Print *, 'After solve_L'
        Print *, 'ACw(sum/max)=', sum(ACw), maxval(ACw)
        WRITE(890,*) 'ENDING'
        WRITE(390,*) 'ENDING'
        WRITE(1390,*) 'ENDING'
        DO IP=NP_GLOBAL,1,-1
          IP_b=ListPos(IP)
          DO J=IAglob(IP_b),IAglob(IP_b+1)-1
            JP_b=JAglob(J)
            JP=ListPosRev(JP_b)
            IF (JP .gt. IP) THEN
              eCoeff=ASPARglob(J)
!              WRITE(900,*) 'IP/JP=', IP, JP
              WRITE(900,*) 'IP_b/JP_b/J=', IP_b, JP_b, J
              WRITE(400,*) IP_b, JP_b
              WRITE(1400,*) IP_b, JP_b, eCoeff
              WRITE(900,*) 'eCoeff, AC12=', eCoeff, ACw(IP_b), ACw(JP_b)
              IF (DO_SOLVE_U) THEN
                ACw(IP_b)=ACw(IP_b) - eCoeff*ACw(JP_b)
              END IF
            END IF
          END DO
          J=I_DIAGglob(IP)
          IF (DO_SOLVE_U) THEN
            ACw(IP_b)=ACw(IP_b)/ASPARglob(I_DIAGglob(IP_b))
          END IF
        END DO
        WRITE(2407,*) 'After solve_U'
        WRITE(2407,*) 'ACw(sum/max)=', sum(ACw), maxval(ACw)
        Print *, 'After solve_U'
        Print *, 'ACw(sum/max)=', sum(ACw), maxval(ACw)
        WRITE(400,*) 'ENDING'
        WRITE(1400,*) 'ENDING'
        TotalErr=sum(abs(ACw - AC1glob))
        TheTol=0.1_rkind
        Sum1=sum(abs(ACw))
        Sum2=sum(abs(AC1glob))
        WRITE(257,*) 'idAC=', idAC, '/', nbEnt
        WRITE(257,*) 'Err/1/2=', TotalErr, Sum1, Sum2
        DO IP=1,np_global
          IF (abs(ACw(IP) - AC1glob(IP)) .gt. TheTol) THEN
            WRITE(257,*) 'IP, AC12=', IP, ACw(IP), AC1glob(IP)
          END IF
        END DO
        Print *, 'Err/1/2=', TotalErr, Sum1, Sum2
      END DO
      END SUBROUTINE
# endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WWM_SOLVER_ALLOCATE(SolDat)
      USE DATAPOOL, only : I5_SolutionData, MNP, MSC, MDC, NNZ
      implicit none
      type(I5_SolutionData), intent(inout) :: SolDat
      allocate(SolDat % AC1(MNP, MSC,MDC))
      allocate(SolDat % AC2(MNP, MSC,MDC))
      allocate(SolDat % AC3(MNP, MSC,MDC))
      allocate(SolDat % AC4(MNP, MSC,MDC))
      allocate(SolDat % AC5(MNP, MSC,MDC))
      allocate(SolDat % AC6(MNP, MSC,MDC))
      allocate(SolDat % AC7(MNP, MSC,MDC))
      allocate(SolDat % AC8(MNP, MSC,MDC))
      allocate(SolDat % AC9(MNP, MSC,MDC))
      allocate(SolDat % ASPAR_block(NNZ, MSC, MDC))
      allocate(SolDat % B_block(MNP, MSC, MDC))
      allocate(SolDat % ASPAR_pc(NNZ, MSC, MDC))
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5_FREE(SolDat)
      USE DATAPOOL, only : I5_SolutionData
      implicit none
      type(I5_SolutionData), intent(inout) :: SolDat
      deallocate(SolDat % AC1)
      deallocate(SolDat % AC2)
      deallocate(SolDat % AC3)
      deallocate(SolDat % AC4)
      deallocate(SolDat % AC5)
      deallocate(SolDat % AC6)
      deallocate(SolDat % AC7)
      deallocate(SolDat % AC8)
      deallocate(SolDat % AC9)
      deallocate(SolDat % ASPAR_block)
      deallocate(SolDat % B_block)
      deallocate(SolDat % ASPAR_pc)
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE I5_SUM(AC, eSum)
      USE DATAPOOL, only : MNP, MSC, MDC, rkind, ZERO
      implicit none
      real(rkind), intent(in) :: AC(MNP,MSC,MDC)
      real(rkind), intent(out) :: eSum
      integer IP,IS,ID
      eSum=ZERO
      DO IP=1,MNP
        DO IS=1,MSC
          DO ID=1,MDC
            eSum=eSum + AC(IP,IS,ID)
          END DO
        END DO
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
# ifdef DEBUG
      SUBROUTINE SYNCHRONIZE_DEATH
      USE elfe_msgp, only : myrank, nproc, itype, comm, istatus, ierr
      implicit none
      integer :: TheVal(nproc)
      integer :: eVal(1)
      integer iProc, I
      IF (myrank == 0) THEN
        TheVal(1)=1
        DO iProc=2,nproc
          CALL MPI_RECV(eVal,1,itype, iProc-1, 997, comm, istatus, ierr)
          TheVal(iProc)=eVal(1)
        END DO
      ELSE
        eVal(1)=myrank+1
        CALL MPI_SEND(eVal,1,itype, 0, 997, comm, ierr)
      END IF
      DO I=1,10000
        WRITE(740+myrank,*) 'GARBAGE to flush the file'
      END DO
      CALL WWM_ABORT('SYNCHRONIZE DEATH FOR YOU')
      END SUBROUTINE
# endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
# ifdef DEBUG
      SUBROUTINE I5_LOCATE_MAX(AC)
      USE DATAPOOL, only : MNP, MSC, MDC, rkind, ZERO
      implicit none
      real(rkind), intent(in) :: AC(MNP,MSC,MDC)
      real(rkind) :: TheMax
      integer IP,IS,ID
      integer IPf,ISf,IDf
      TheMax=ZERO
      IPf=0
      ISf=0
      IDf=0
      DO IP=1,MNP
        DO IS=1,MSC
          DO ID=1,MDC
            IF (AC(IP,IS,ID) .gt. TheMax) THEN
              TheMax=AC(IP,IS,ID)
              IPf=IP
              ISf=IS
              IDf=ID
            END IF
          END DO
        END DO
      END DO
      Print *, 'f IP=', IPf, 'IS=', ISf, 'ID=', IDf
      END SUBROUTINE
# endif
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WWM_SOLVER_EIMPS(LocalColor, SolDat)
      USE DATAPOOL, only : LocalColorInfo, I5_SolutionData
      USE DATAPOOL, only : rkind, MSC, MDC, AC2, MNP, NNZ
      USE DATAPOOL, only : PCmethod, IOBPD, ZERO
      USE elfe_msgp, only : myrank
      implicit none
      type(LocalColorInfo), intent(inout) :: LocalColor
      type(I5_SolutionData), intent(inout) :: SolDat
      real(rkind) :: U(MNP), ASPAR(NNZ), B(MNP)
      integer IS, ID, IP
      DO IS=1,MSC
        DO ID=1,MDC
          U=AC2(:,IS,ID)
          CALL EIMPS_ASPAR_B(IS, ID, ASPAR, B, U)
          SolDat % ASPAR_block(:,IS,ID)=ASPAR
          SolDat % B_block(:,IS,ID)=B
        END DO
      END DO
# ifdef DEBUG
      WRITE(myrank+740,*) 'Before EXCHANGE_P4D_WWM_TR'
# endif
      CALL EXCHANGE_P4D_WWM_TR(SolDat % B_block)
# ifdef DEBUG
      WRITE(myrank+740,*) 'Before EXCHANGE_ASPAR'
# endif
      CALL EXCHANGE_ASPAR(SolDat % ASPAR_block)
# ifdef DEBUG
      WRITE(myrank+740,*) 'maxval ASPAR_block=', maxval(SolDat % ASPAR_block)
# endif
      CALL I5_CREATE_PRECOND(LocalColor, SolDat, PCmethod)
# ifdef DEBUG
      WRITE(myrank+740,*) 'maxval ASPAR_pc=', maxval(SolDat % ASPAR_pc)
# endif
      CALL I5_BCGS_SOLVER(LocalColor, SolDat)
      DO IP=1,MNP
        DO IS=1,MSC
          DO ID=1,MDC
            AC2(IP,IS,ID)=MAX(ZERO, AC2(IP,IS,ID))*MyREAL(IOBPD(ID,IP))
          END DO
        END DO
      END DO
# ifdef DEBUG
      IF (myrank == 1) THEN
        Write(myrank+591,*) 'Clearing ENDING'
      END IF
# endif
      END SUBROUTINE
#endif
END MODULE WWM_PARALL_SOLVER
