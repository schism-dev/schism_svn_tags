#ifdef PETSC
      MODULE PETSCPOOL

#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
#include "finclude/petscao.h"
#include "finclude/petscis.h"
#include "finclude/petscdraw.h"
#include "finclude/petscviewer.h"
#include "finclude/petscviewer.h90"
#include "finclude/petscviewerdef.h"
#include "finclude/petscpc.h"

      ! petsc stuff
      Mat matrix;
      Vec :: myB, myX
      PetscScalar, pointer :: myBtemp(:), myXtemp(:)

      KSP Solver;
      PC Prec;
      KSPConvergedReason reason;
      PetscInt iterationen;
      PetscErrorCode ierr;

      Vec     :: r;                         ! solution, residual vectors
      PetscViewer :: viewer ! draw matrix
      PetscDraw   :: draw   ! extended viewer. e.g. pause
      real*8 :: matInfo(MAT_INFO_SIZE)

      ! petsc parallel stuff

      ! rank of a process, number of processors
      PetscMPIInt :: rank,nProcs
      integer :: nNodesWithoutInterfaceGhosts
      ! number of ghost nodes (interface + ghost)
      integer :: nghost

      ! Mappings
      ! the "old mapping" has resident and interface nodes
      ! the others have only resident nodes
      PetscInt, allocatable :: ALO2PGO(:), AGO2PGO(:)
      PetscInt, allocatable :: PLO2PGO(:), PGO2AGO(:)
      PetscInt, allocatable :: AGO2ALO(:), ALO2AGO(:)
      PetscInt, allocatable :: ALOold2ALO(:), ALO2ALOold(:)
      PetscInt, allocatable :: PGO2PLO(:)

      ! Application Ordering
      AO :: appOrdering

      integer, allocatable :: onlyNodes(:), onlyGhosts(:), onlyGhostsOldLocalMapping(:)

      ! holds the range of vector X that every processor own
      PetscInt:: ranges

      ! Debug stuff
      integer :: timePetscCalled = 0
      ! Same as first matrix but filled with app order instead of petsc order. easyer to debug
      Mat matrixAppOrder;
      Vec myBAppOrder;

      contains

      ! set the sover tolerances. arguments are optional
      subroutine setSolverTolerance(rtol, abstol, dtol, maxits)
        Real*8, optional :: rtol    ! the relative convergence tolerance
        Real*8, optional :: abstol  ! the absolute convergence tolerance
        Real*8, optional :: dtol    ! the divergence tolerance
        Integer, optional ::  maxits  ! maximum number of iterations

        PetscReal oldrtol
        PetscReal oldabstol
        PetscReal olddtol
        PetscInt  oldmaxits

        PetscReal newrtol
        PetscReal newabstol
        PetscReal newdtol
        PetscInt  newmaxits

        call KSPGetTolerances(solver, oldrtol, oldabstol, olddtol, oldmaxits, ierr);CHKERRQ(ierr)

        if(present(rtol)) then
          newrtol = rtol
        else
          newrtol = oldrtol;
        endif

        if(present(abstol)) then
          newabstol = abstol
        else
          newabstol = oldabstol
        endif

        if(present(dtol)) then
          newdtol = dtol
        else
          newdtol = olddtol
        endif

        if(present(maxits)) then
          newmaxits = maxits
        else
          newmaxits = oldmaxits
        endif

        call KSPSetTolerances(solver, rtol, abstol, dtol, maxits, ierr);CHKERRQ(ierr)

      end subroutine

      ! print out the solver tolerances
      subroutine printSolverTolerance
        use datapool, only : stat
        PetscReal rtol    ! the relative convergence tolerance
        PetscReal abstol  ! the absolute convergence tolerance
        PetscReal dtol    ! the divergence tolerance
        PetscInt  maxits  ! maximum number of iterations

        call KSPGetTolerances(solver, rtol, abstol, dtol, maxits, ierr);CHKERRQ(ierr)
        write(stat%fhndl,*) "relative convergence tolerance", rtol
        write(stat%fhndl,*) "absolute convergence tolerance", abstol
        write(stat%fhndl,*) "divergence tolerance", dtol
        write(stat%fhndl,*) "maximum number of iterations", maxits
      end subroutine

      ! print the matrix, x and RHs vector
#ifdef SELFE
      subroutine printPDE(ISS, IDD)
         use elfe_msgp, only: myrank
         use datapool, only : stat
         IMPLICIT NONE

          integer, intent(in) :: ISS, IDD
          PetscViewer :: viewer

          ! create a viewer for vectors that print the index number
          call PetscViewerCreate(PETSC_COMM_WORLD, viewer, ierr);CHKERRQ(ierr)
          call PetscViewerSetType(viewer, PETSCVIEWERASCII, ierr);CHKERRQ(ierr)
          call PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_INDEX, ierr);CHKERRQ(ierr)

          if(myrank == 0) write(stat%fhndl,*) "IS ID petscCalls", ISS, IDD, timePetscCalled

          if(myrank == 0) write(stat%fhndl,*) "matrix: in Application Order"
          ! draw the matrix
          call MatView(matrixAppOrder, viewer, ierr );CHKERRQ(ierr)

          ! draw x
          if(myrank == 0) write(stat%fhndl,*) "X in PETSC Order"
          call VecView(myX, viewer, ierr);CHKERRQ(ierr)

          ! draw rhs
          if(myrank == 0) write(stat%fhndl,*) "rhs in Application Order"
          call VecView(myBAppOrder, viewer, ierr);CHKERRQ(ierr)
      end subroutine printPDE

      subroutine blah(ISS, IDD, ASPAR, B, X, problemNode)
        USE DATAPOOL
         use elfe_glbl, only: np_global, np, npg, npa, nnp, inp, iplg
         ! iplg1 points to elfe_glbl::ipgl because ipgl exist allreay as integer in this function
         use elfe_glbl, only: ipgl1=> ipgl
         use elfe_msgp
         IMPLICIT NONE
         include 'mpif.h'

        integer, intent(in) :: ISS, IDD
        REAL*8, intent(in) ::  ASPAR(:)
        REAL*8, intent(in)  :: X(:)
        REAL*8, intent(in)  :: B(:)
        integer, intent(in) :: problemNode

        integer :: i
        PetscScalar :: eEntry
        PetscInt :: eCol
!         PetscInt :: eRow

        if(myrank == 1 .or. nProcs == 1) then
          if(myrank == 1) then
            write(stat%fhndl,*) "Infos from myrank 1"
          else
            write(stat%fhndl,*) "Infos from myrank 0"
          endif
          write(stat%fhndl,*) "infos for global node ", problemNode, "+1"
          write(stat%fhndl,*) "IS", ISS, " ID", IDD

          ! print only one row
          write(stat%fhndl,*) "matrix koeff"
          do i = 1, np_global
            eCol = AGO2PGO(i-1)
!              matrix, nRows, global row numner, number cols, global col number, value
            CALL MatGetValues(matrix, 1, AGO2PGO(problemNode), 1, eCol, eEntry, ierr); CHKERRQ(ierr);
            if(dabs(eEntry) > 1e-16)  write(stat%fhndl,*) "AGO value", PGO2AGO(ecol), eEntry
!           if(dabs(eEntry) > 1e-16)  write(*,*) PGO2AGO(ecol)
          END DO

        endif
      end subroutine
#endif
      ! print some matrix information like number of nonzeros, memory allocated
      !subroutine printMatrixInformation(Mat matrix)
      !   call MatGetInfo(matrix, MAT_LOCAL, matInfo, ierr);CHKERRQ(ierr);
      !   write(*,*) rank, "block size", matInfo(MAT_INFO_BLOCK_SIZE)
      !   write(*,*) rank, "number of nonzeros allocated", matInfo(MAT_INFO_NZ_ALLOCATED)
      !   write(*,*) rank, "number of nonzeros used", matInfo(MAT_INFO_NZ_USED)
      !   write(*,*) rank, "number of nonzeros uneeded", matInfo(MAT_INFO_NZ_UNNEEDED)
      !   write(*,*) rank, "memory allocated", matInfo(MAT_INFO_MEMORY)
      !   write(*,*) rank, "number of matrix assemblies called ", matInfo(MAT_INFO_ASSEMBLIES)
      !   write(*,*) rank, "number of mallocs during MatSetValues()", matInfo(MAT_INFO_MALLOCS)
      !   write(*,*) rank, "fill ratio for LU/ILU given", matInfo(MAT_INFO_FILL_RATIO_GIVEN)
      !   write(*,*) rank, "fill ratio for LU/ILU needed", matInfo(MAT_INFO_FILL_RATIO_NEEDED)
      !   write(*,*) rank, "number of mallocs during factorization", matInfo(MAT_INFO_FACTOR_MALLOCS)
      !end subroutine

      SUBROUTINE PETSC_INIT()

        USE DATAPOOL

        IMPLICIT NONE

        ! petsc wants indices startet from 0
        IA = IA -1
        JA = JA -1


         call PetscInitialize(PETSC_NULL_CHARACTER,ierr);CHKERRQ(ierr)

! create sparse matrix
!        call MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, MNP, MNP, IA, JA, ASPAR, matrix, ierr);CHKERRQ(ierr)
        call MatCreate(PETSC_COMM_SELF, matrix, ierr);CHKERRQ(ierr)
        call MatSetType(matrix, MATSEQAIJ , ierr);CHKERRQ(ierr)
        call MatSetSizes(matrix, PETSC_DECIDE, PETSC_DECIDE, MNP, MNP, ierr);CHKERRQ(ierr)

!create x vector
         call VecCreate(PETSC_COMM_SELF, myX, ierr);CHKERRQ(ierr)
         call VecSetSizes(myX, PETSC_DECIDE, MNP, ierr);CHKERRQ(ierr)
         call VecSetType(myX,VECSEQ, ierr);CHKERRQ(ierr)

! create vec myB
         call VecCreate(PETSC_COMM_SELF, myB, ierr);CHKERRQ(ierr)
         call VecSetSizes(myB, PETSC_DECIDE, MNP, ierr);CHKERRQ(ierr)
         call VecSetType(myB,VECSEQ, ierr);CHKERRQ(ierr)

! create solver
         call KSPCreate(PETSC_COMM_SELF, Solver, ierr);CHKERRQ(ierr)
!          call KSPSetOperators(Solver, matrix, matrix, SAME_NONZERO_PATTERN, ierr); CHKERRQ(ierr)
         call KSPSetOperators(Solver, matrix, matrix, 0, ierr); CHKERRQ(ierr)

        call KSPSetType(Solver,KSPBCGS, ierr);CHKERRQ(ierr)
!          call KSPSetType(Solver, KSPIBCGS, ierr);CHKERRQ(ierr) ! this solver create a segfault

! Create preconditioner
         call KSPGetPC(Solver, Prec, ierr);CHKERRQ(ierr)
         call PCSetType(Prec, PCILU, ierr);CHKERRQ(ierr)


      END SUBROUTINE

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE  EIMPS_PETSC(ISS, IDD)

         USE DATAPOOL
         !use elfe_glbl, only: iplg
         !use elfe_msgp
         IMPLICIT NONE

         !include 'mpif.h'

         INTEGER   :: ISS,IDD

         INTEGER :: I, J
!          INTEGER :: KK

         INTEGER :: IP, IPGL, IE, POS

         INTEGER :: I1, I2, I3

         REAL*8 :: DTK, TMP3

         REAL*8 :: LAMBDA(2)
         REAL*8 :: FL11, FL12, FL21, FL22, FL31, FL32
         REAL*8 :: CRFS(3), K1, KM(3), K(3), TRIA03

         REAL*8 :: DELTAL(3,MNE)
         REAL*8 :: KP(3,MNE), NM(MNE)
         REAL*8 :: U(MNP), C(2,MNP)

         REAL*8 :: X(MNP)
         REAL*8 :: B(MNP)

         REAL*8 ::  ASPAR(NNZ)

         REAL    :: TIME1, TIME2, TIME3, TIME4, TIME5

         INTEGER :: POS_TRICK(3,2)


         ! petsc stuff
         integer :: counter
         integer :: nnznew
         integer :: temp
         PetscInt :: ncols
         PetscInt :: eCol, eRow
         PetscScalar :: eValue

         CALL CPU_TIME(TIME1)

         POS_TRICK(1,1) = 2
         POS_TRICK(1,2) = 3
         POS_TRICK(2,1) = 3
         POS_TRICK(2,2) = 1
         POS_TRICK(3,1) = 1
         POS_TRICK(3,2) = 2

         CALL CADVXY(ISS,IDD,C)
!
!        Calculate countour integral quantities ...
!
         DO IE = 1, MNE
           I1 = INE(1,IE)
           I2 = INE(2,IE)
           I3 = INE(3,IE)
           LAMBDA(1) = ONESIXTH * (C(1,I1)+C(1,I2)+C(1,I3))
           LAMBDA(2) = ONESIXTH * (C(2,I1)+C(2,I2)+C(2,I3))
           K(1)  = LAMBDA(1) * IEN(1,IE) + LAMBDA(2) * IEN(2,IE)
           K(2)  = LAMBDA(1) * IEN(3,IE) + LAMBDA(2) * IEN(4,IE)
           K(3)  = LAMBDA(1) * IEN(5,IE) + LAMBDA(2) * IEN(6,IE)
           KP(1,IE) = MAX(0.d0,K(1))
           KP(2,IE) = MAX(0.d0,K(2))
           KP(3,IE) = MAX(0.d0,K(3))
           KM(1) = MIN(0.d0,K(1))
           KM(2) = MIN(0.d0,K(2))
           KM(3) = MIN(0.d0,K(3))
           FL11 = C(1,I2)*IEN(1,IE)+C(2,I2)*IEN(2,IE)
           FL12 = C(1,I3)*IEN(1,IE)+C(2,I3)*IEN(2,IE)
           FL21 = C(1,I3)*IEN(3,IE)+C(2,I3)*IEN(4,IE)
           FL22 = C(1,I1)*IEN(3,IE)+C(2,I1)*IEN(4,IE)
           FL31 = C(1,I1)*IEN(5,IE)+C(2,I1)*IEN(6,IE)
           FL32 = C(1,I2)*IEN(5,IE)+C(2,I2)*IEN(6,IE)
           CRFS(1) =  - ONESIXTH *  (2.0d0 *FL31 + FL32 + FL21 + 2.0d0 * FL22 )
           CRFS(2) =  - ONESIXTH *  (2.0d0 *FL32 + 2.0d0 * FL11 + FL12 + FL31 )
           CRFS(3) =  - ONESIXTH *  (2.0d0 *FL12 + 2.0d0 * FL21 + FL22 + FL11 )
           DELTAL(:,IE) = CRFS(:)- KP(:,IE)
           NM(IE)       = 1.d0/MIN(DBLE(THR),SUM(KM(:)))
         END DO

         CALL CPU_TIME(TIME2)

         U = DBLE(AC2(:,ISS,IDD))

         J     = 0    ! Counter ...
         ASPAR = 0.d0 ! Mass matrix ...
         B     = 0.d0 ! Right hand side ...
!
! ... assembling the linear equation system ....
!
         DO IP = 1, MNP
           DO I = 1, CCON(IP)
             J = J + 1
             IE    =  IE_CELL(J)
             POS   =  POS_CELL(J)
             K1    =  KP(POS,IE) ! Flux Jacobian
             TRIA03 = ONETHIRD * TRIA(IE)
             DTK   =  K1 * DT4A
             TMP3  =  DTK * NM(IE)
             I1    =  POSI(1,J) ! Position of the recent entry in the ASPAR matrix ... ASPAR is shown in fig. 42, p.122
             I2    =  POSI(2,J)
             I3    =  POSI(3,J)
!             IF (IOBP(IP) .NE. 2 .AND. IOBPD(IDD,IP) .EQ. 1) THEN
               ASPAR(I1) =  TRIA03 + DTK - TMP3 * DELTAL(POS             ,IE) + ASPAR(I1)  ! Diagonal entry
               ASPAR(I2) =               - TMP3 * DELTAL(POS_TRICK(POS,1),IE) + ASPAR(I2)  ! off diagonal entries ...
               ASPAR(I3) =               - TMP3 * DELTAL(POS_TRICK(POS,2),IE) + ASPAR(I3)
               B(IP)     =  B(IP) + TRIA03 * U(IP)
!            END IF
           END DO !I: loop over connected elements ...
         END DO !IP

         DO IP = 1, MNP
           IF (IOBPD(IDD,IP) .EQ. 0) THEN
             ASPAR(I_DIAG(IP)) = SI(IP)
             B(IP)             = 0.d0
           END IF
         END DO

         IF (LBCWA .OR. LBCSP) THEN
           IF (LINHOM) THEN
             DO IP = 1, IWBMNP
               IPGL = IWBNDLC(IP)
               ASPAR(I_DIAG(IPGL)) = SI(IPGL) ! Add source term to the diagonal
               B(IPGL)             = SI(IPGL) * WBAC(ISS,IDD,IP)
             END DO
           ELSE
             DO IP = 1, IWBMNP
               IPGL = IWBNDLC(IP)
                ASPAR(I_DIAG(IPGL)) = SI(IPGL) ! Add source term to IDthe diagonal
               B(IPGL)              = SI(IPGL) * WBAC(ISS,IDD,1)
             END DO
           ENDIF
         END IF

         IF (ICOMP .GE. 2) THEN
           DO IP = 1, MNP
             IF (IOBP(IP) .EQ. 0) THEN
               ASPAR(I_DIAG(IP)) = ASPAR(I_DIAG(IP)) + IMATDAA(IP,ISS,IDD) * DT4A * SI(IP) * IOBPD(IDD,IP)  ! Add source term to the diagonal
               B(IP)             = B(IP) + IMATRAA(IP,ISS,IDD) * DT4A * SI(IP) * IOBPD(IDD,IP) ! Add source term to the right hand side
             END IF
           END DO
         END IF

         CALL CPU_TIME(TIME3)

!> \todo for efficiency reason this will be moved out of this part when the code is working, AR.
! fill the matrix
         call MatZeroEntries(matrix, ierr);CHKERRQ(ierr)
         counter = 1
         nnznew=0

         do i = 1, MNP
           ncols = IA(i+1) - IA(i)
           eRow = i - 1

           ! insert col by col into matrix
           do j = 1, ncols
             ! the value we want to insert
             eValue=ASPAR(counter)
             ! get the col index in old counting
             temp = JA(counter)
             counter = counter + 1
             nnznew=nnznew + 1

             ! increase temp because iplg counts from one and JA from zero
             eCol = temp
             call MatSetValue(matrix, eRow, eCol, eValue, ADD_VALUES, ierr);CHKERRQ(ierr)
           end do
         end do ! np
         call MatAssemblyBegin(matrix,MAT_FINAL_ASSEMBLY, ierr);CHKERRQ(ierr);
         call MatAssemblyEnd(matrix,MAT_FINAL_ASSEMBLY, ierr);CHKERRQ(ierr);


!fill RHS vector
         call VecGetArrayF90(myB, myBtemp, ierr);CHKERRQ(ierr)
         myBtemp = B
         call VecRestoreArrayF90(myB, myBtemp, ierr);CHKERRQ(ierr)
         call VecAssemblyBegin(myB, ierr);CHKERRQ(ierr);
         call VecAssemblyEnd(myB, ierr);CHKERRQ(ierr);

         CALL CPU_TIME(TIME4)
! Solve
         call KSPSolve(Solver, myB, myX, ierr);CHKERRQ(ierr)

         CALL CPU_TIME(TIME5)

         !WRITE(*,*) TIME2-TIME1, TIME3-TIME2, TIME4-TIME3, TIME5-TIME4

         call KSPGetConvergedReason(Solver, reason, ierr);CHKERRQ(ierr)
         if (reason .LT. 0) then
           write(stat%fhndl,*) "Failure to converge\n"
           stop
         else
         call KSPGetIterationNumber(Solver, iterationen, ierr);CHKERRQ(ierr)
            if(iterationen /= 0)  write(stat%fhndl,*) "Number of iterations", iss,idd,iterationen
!            write(*,*) "Number of iterations", iterationen
         endif

         ! write the soluten from vec into fortran array
         X    = 0.d0
         call VecGetArrayF90(myX, myXtemp, ierr);CHKERRQ(ierr)
         X = myXtemp
         call VecRestoreArrayF90(myX, myXtemp, ierr);CHKERRQ(ierr)

         IF (SUM(X) .NE. SUM(X)) STOP 'NaN in X'
!          AC2(:, ISS, IDD) = SNGL(MAX(0.d0,X))

         AC2(:,ISS,IDD) = X

!          CALL CPU_TIME(TIME4)

         RETURN
        END SUBROUTINE
!
! - Initialize Petsc
! - seperate interface nodes from resident nodes
! - create the mappings
! - create the matrix and vectors
! - create the solver
!
#ifdef SELFE
      SUBROUTINE PETSC_INIT_PARALLEL()
        USE DATAPOOL, only: MNP, CCON, IA, JA, NNZ, STAT
        ! np_global - # nodes gloabl
        ! np        - # nodes local non augmented
        ! npg       - # ghost
        ! npa       - # nodes aufmented
        ! nea       - # elemnts augmented
        ! int::iplg(ip)           global element index. local to global LUT
        ! llsit_type::ipgl(ipgb)  ipgb is a global node. global to local LUT
        ! int::nnp(ip)            total # of surrounding nodes for node ip
        ! int::inp(ip, 1:nnp(ip)) list of surrounding nodes
        USE elfe_glbl, only: np_global, np, npg, npa, ne_global, nea, iplg, ipgl, nnp, inp, llist_type
        IMPLICIT NONE
        integer :: i

        call PetscInitialize(PETSC_NULL_CHARACTER,ierr);CHKERRQ(ierr);
        !call PetscPrintf(PETSC_COMM_WORLD, "PETSC_INIT_PARALLEL\n", ierr)

        call MPI_Comm_rank(MPI_COMM_WORLD,rank, ierr);CHKERRQ(ierr);
        call MPI_Comm_size(MPI_COMM_WORLD,nProcs, ierr);CHKERRQ(ierr);

        call PetscViewerDrawOpen(PETSC_COMM_WORLD,PETSC_NULL_CHARACTER, PETSC_NULL_CHARACTER, 0, 0, 600, 600,viewer,ierr);CHKERRQ(ierr);

        ! petsc wants indices startet from 0
        IA = IA -1
        JA = JA -1

! exclude interface nodes from the local mapping,list,iplg whatever
! a interface node is not a ghost node. so one can find them in the range from 1 to np.
! a interface node is owned by two (or more threads?)
! to the detect an interface node i use this method:

! a) first geht the global id
! b) with the gloabl id, get the local rank
! c) check if there is a linked list associated with the local node
! d) if the rank from the next node in the list ist less then the curren rank, ignore this node
! f) else copy this node it into the new mappings

! there are two runs. in the first, count the number of non interface nodes
! in the second run, allocate and fill the mappings
        nNodesWithoutInterfaceGhosts = 0
        do i = 1, np
            if(ASSOCIATED(ipgl( iplg(i) )%next) ) then
              if(( ipgl( iplg(i) )%next%rank .le. rank )) then
                cycle
              end if
            end if
            nNodesWithoutInterfaceGhosts = nNodesWithoutInterfaceGhosts + 1
        end do

! create App Local <-> App Global mappings
! create App Local Old <-> App Local new mapping (without interface nodes)
        allocate(ALO2AGO(0:nNodesWithoutInterfaceGhosts-1))
        allocate(AGO2ALO(0:np_global-1))
        AGO2ALO = -1

        allocate(ALO2ALOold(0:nNodesWithoutInterfaceGhosts-1))
        allocate(ALOold2ALO(0:np + npg -1))
        ! nodes in the old App local mapping become a -990 to indicate, that this node is an interface node
        ALOold2ALO = -999

        nNodesWithoutInterfaceGhosts=0
        do i = 1, np
            if(ASSOCIATED(ipgl( iplg(i) )%next) ) then
              if(( ipgl( iplg(i) )%next%rank .le. rank )) then
                cycle
              end if
            end if

            ALO2AGO(nNodesWithoutInterfaceGhosts) = iplg(i) -1
            AGO2ALO(iplg(i)-1) = nNodesWithoutInterfaceGhosts

            ALOold2ALO(i-1) = nNodesWithoutInterfaceGhosts
            ALO2ALOold(nNodesWithoutInterfaceGhosts) = i-1

            nNodesWithoutInterfaceGhosts = nNodesWithoutInterfaceGhosts + 1
        end do

! create PETsc Local -> Global mapping
        allocate(PLO2PGO(0:nNodesWithoutInterfaceGhosts))
        allocate(PGO2PLO(0:np_global))
        PGO2PLO = -999
        call MPI_Scan(nNodesWithoutInterfaceGhosts, ranges, 1, MPIU_INTEGER, MPI_SUM, PETSC_COMM_WORLD, ierr);CHKERRQ(ierr)

        do i = 1, nNodesWithoutInterfaceGhosts
          PLO2PGO(i-1) = ranges - nNodesWithoutInterfaceGhosts + i -1
!> \todo warum issen das auskommentiert? PGO2PLO wird mom auch nicht gebraucht, daher nicht so schlimm
!> \todo dont need PGO2PLO atm
!         PGO2PLO(ranges - nNodesWithoutInterfaceGhosts + i -1 ) = i - 1
        end do

! create App Global <-> PETsc Global mapping
        allocate(AGO2PGO(0:np_global-1))
        allocate(PGO2AGO(0:np_global-1))
        call AOCreateBasic(MPI_COMM_WORLD, nNodesWithoutInterfaceGhosts, ALO2AGO, PLO2PGO, appOrdering, ierr);CHKERRQ(ierr)
!         call AOView(appOrdering, PETSC_VIEWER_STDOUT_WORLD, ierr) ;CHKERRQ(ierr)

        do i = 1, np_global
          AGO2PGO(i-1) = i-1
          PGO2AGO(i-1) = i-1
        end do

        call AOApplicationToPetsc(appOrdering, np_global, AGO2PGO, ierr);CHKERRQ(ierr)
        call AOPetscToApplication(appOrdering, np_global, PGO2AGO, ierr);CHKERRQ(ierr)
        call AODestroy(appOrdering, ierr);CHKERRQ(ierr)

! create onlyGhosts array
        nghost=np+npg-nNodesWithoutInterfaceGhosts
        allocate(onlyGhosts(0:nghost-1))
        allocate(onlyGhostsOldLocalMapping(0:nghost-1))
        nghost = 0
        do i = 1, np
            if(ASSOCIATED(ipgl( iplg(i) )%next) ) then
              if(( ipgl( iplg(i) )%next%rank .le. rank )) then
                onlyGhosts(nghost)= AGO2PGO( iplg(i)-1 )
                onlyGhostsOldLocalMapping(nghost) = i
                nghost = nghost + 1
                cycle
              end if
            end if
        end do

        DO i=1,npg
          onlyGhosts(nghost) = AGO2PGO( iplg(np+i)-1)
          onlyGhostsOldLocalMapping(nghost) = np+i
          nghost = nghost + 1
        END DO

        write(stat%fhndl,*) rank, "Global Number of Nodes" , np_global
        write(stat%fhndl,*) rank, "Local Number of resident nodes", np
        write(stat%fhndl,*) rank, "Local Number of ghost nodes", npg
        write(stat%fhndl,*) rank, "local Number of nodes in augmented subdomain (np+npg)", npa
        write(stat%fhndl,*) rank, "Local Number of nodes without interface and ghost nodes", nNodesWithoutInterfaceGhosts
        write(stat%fhndl,*) rank, "Local Number of ghost + interface nodes", nghost

! create matrix
         call MatCreate(MPI_COMM_WORLD, matrix, ierr);CHKERRQ(ierr)
         call MatSetType(matrix, MATMPIAIJ, ierr);CHKERRQ(ierr)
         call MatSetSizes(matrix, nNodesWithoutInterfaceGhosts, nNodesWithoutInterfaceGhosts, np_global, np_global, ierr);CHKERRQ(ierr)

         call MatCreate(MPI_COMM_WORLD, matrixAppOrder, ierr);CHKERRQ(ierr)
         call MatSetType(matrixAppOrder, MATMPIAIJ, ierr);CHKERRQ(ierr)
         call MatSetSizes(matrixAppOrder, nNodesWithoutInterfaceGhosts, nNodesWithoutInterfaceGhosts, np_global, np_global, ierr);CHKERRQ(ierr)

! create X vector
         call VecCreateGhost(MPI_COMM_WORLD, nNodesWithoutInterfaceGhosts, np_global, nghost, onlyGhosts, myX, ierr);CHKERRQ(ierr)

!          call VecCreate(MPI_COMM_WORLD, myX, ierr); CHKERRQ(ierr)
!          call VecSetSizes(myX, nNodesWithoutInterfaceGhosts, np_global, ierr);CHKERRQ(ierr)
!          call VecSetType(myX, VECMPI, ierr); CHKERRQ(ierr)


! create B vector
          call VecCreateGhost(MPI_COMM_WORLD, nNodesWithoutInterfaceGhosts, np_global, nghost, onlyGhosts, myB, ierr);CHKERRQ(ierr)
          call VecCreateGhost(MPI_COMM_WORLD, nNodesWithoutInterfaceGhosts, np_global, nghost, onlyGhosts, myBAppOrder, ierr);CHKERRQ(ierr)

!          call VecCreate(MPI_COMM_WORLD, myB, ierr); CHKERRQ(ierr)
!          call VecSetSizes(myB, nNodesWithoutInterfaceGhosts, np_global, ierr);CHKERRQ(ierr)
!          call VecSetType(myB, VECMPI, ierr); CHKERRQ(ierr)

! create solver
         call KSPCreate(MPI_COMM_WORLD,Solver,ierr);CHKERRQ(ierr)
!          call KSPSetOperators(Solver, matrix, matrix, SAME_NONZERO_PATTERN, ierr); CHKERRQ(ierr)
         call KSPSetOperators(Solver,matrix,matrix,0, ierr);CHKERRQ(ierr)
         call KSPSetType(Solver,KSPLGMRES, ierr); CHKERRQ(ierr) ! one good solution ...
         !call KSPSetType(Solver,KSPDGMRES, ierr); CHKERRQ(ierr) ! one good solution ...
         !call KSPSetType(Solver,KSPBCGSL,ierr); CHKERRQ(ierr)
         ! Use the old solution vom AC2 to make the solver faster
         call KSPSetInitialGuessNonzero(Solver,PETSC_TRUE,ierr);CHKERRQ(ierr)
         call setSolverTolerance(1.D-8,1.D-8,10000.D0,10000)
         call printSolverTolerance

! Create preconditioner
         call KSPGetPC(Solver, Prec, ierr);CHKERRQ(ierr);
         call PCSetType(Prec, PCHYPRE, ierr);CHKERRQ(ierr);
         call PCHYPRESetType(Prec,"euclid",ierr); CHKERRQ(ierr)


      END SUBROUTINE
#endif

#ifdef SELFE
      SUBROUTINE  EIMPS_PETSC_PARALLEL(ISS, IDD)
         USE DATAPOOL
         use elfe_glbl, only: np_global, np, npg, npa, nnp, inp, iplg
         ! iplg1 points to elfe_glbl::ipgl because ipgl exist allreay as integer in this function
         use elfe_glbl, only: ipgl1=> ipgl
         use elfe_msgp
         IMPLICIT NONE
         include 'mpif.h'
         INTEGER, intent(in) :: ISS, IDD

         INTEGER :: I, J
         INTEGER :: IP, IPGL, IE, POS
         INTEGER :: I1, I2, I3
         INTEGER :: POS_TRICK(3,2)

         REAL*8  :: DTK, TMP3
         REAL*8  :: LAMBDA(2)
         REAL*8  :: FL11, FL12, FL21, FL22, FL31, FL32
         REAL*8  :: CRFS(3), K1, KM(3), K(3), TRIA03
         REAL*8  :: DELTAL(3,MNE)
         REAL*8  :: KP(3,MNE), NM(MNE)
         REAL*8  :: U(MNP), C(2,MNP)
         REAL*8  :: X(MNP)
         REAL*8  :: B(MNP)
         REAL*8  ::  ASPAR(NNZ)

         REAL    ::  TIME2

!
! Petsc stuff
!
         PetscInt :: ncols
         PetscInt :: eCol, eRow
         PetscScalar :: eEntry
         PetscScalar :: eValue

         integer :: counter
         integer :: nnznew

         INTEGER :: DoPrint
         INTEGER :: temp

        VecScatter     ctx;
        Vec            vout;


         ! Norms for myB and myX
         PetscScalar :: eNormB1, eNormX1
         PetscScalar :: eNormB2, eNormX2
         PetscScalar :: eNormBi, eNormXi

!         timePetscCalled = timePetscCalled +1

         POS_TRICK(1,1) = 2
         POS_TRICK(1,2) = 3
         POS_TRICK(2,1) = 3
         POS_TRICK(2,2) = 1
         POS_TRICK(3,1) = 1
         POS_TRICK(3,2) = 2

         CALL CADVXY(ISS,IDD,C)
!
!        Calculate countour integral quantities ...
!
         DO IE = 1, MNE
           I1 = INE(1,IE)
           I2 = INE(2,IE)
           I3 = INE(3,IE)
           LAMBDA(1) = ONESIXTH * (C(1,I1)+C(1,I2)+C(1,I3))
           LAMBDA(2) = ONESIXTH * (C(2,I1)+C(2,I2)+C(2,I3))
           K(1)  = LAMBDA(1) * IEN(1,IE) + LAMBDA(2) * IEN(2,IE)
           K(2)  = LAMBDA(1) * IEN(3,IE) + LAMBDA(2) * IEN(4,IE)
           K(3)  = LAMBDA(1) * IEN(5,IE) + LAMBDA(2) * IEN(6,IE)
           KP(1,IE) = MAX(0.d0,K(1))
           KP(2,IE) = MAX(0.d0,K(2))
           KP(3,IE) = MAX(0.d0,K(3))
           KM(1) = MIN(0.d0,K(1))
           KM(2) = MIN(0.d0,K(2))
           KM(3) = MIN(0.d0,K(3))
           FL11 = C(1,I2)*IEN(1,IE)+C(2,I2)*IEN(2,IE)
           FL12 = C(1,I3)*IEN(1,IE)+C(2,I3)*IEN(2,IE)
           FL21 = C(1,I3)*IEN(3,IE)+C(2,I3)*IEN(4,IE)
           FL22 = C(1,I1)*IEN(3,IE)+C(2,I1)*IEN(4,IE)
           FL31 = C(1,I1)*IEN(5,IE)+C(2,I1)*IEN(6,IE)
           FL32 = C(1,I2)*IEN(5,IE)+C(2,I2)*IEN(6,IE)
           CRFS(1) =  - ONESIXTH *  (2.0d0 *FL31 + FL32 + FL21 + 2.0d0 * FL22 )
           CRFS(2) =  - ONESIXTH *  (2.0d0 *FL32 + 2.0d0 * FL11 + FL12 + FL31 )
           CRFS(3) =  - ONESIXTH *  (2.0d0 *FL12 + 2.0d0 * FL21 + FL22 + FL11 )
           DELTAL(:,IE) = CRFS(:)- KP(:,IE)
           NM(IE)       = 1.d0/MIN(DBLE(THR),SUM(KM(:)))
         END DO

         CALL CPU_TIME(TIME2)

         U = DBLE(AC2(:,ISS,IDD))
         call exchange_p2d(U)

         J     = 0    ! Counter ...
         ASPAR = 0.d0 ! Mass matrix ...
         B     = 0.d0 ! Right hand side ...
!
! ... assembling the linear equation system ....
!
         DO IP = 1, MNP
           DO I = 1, CCON(IP)
             J = J + 1
             IE    =  IE_CELL(J)
             POS   =  POS_CELL(J)
             K1    =  KP(POS,IE) ! Flux Jacobian
             TRIA03 = ONETHIRD * TRIA(IE)
             DTK   =  K1 * DT4A
             TMP3  =  DTK * NM(IE)
             I1    =  POSI(1,J) ! Position of the recent entry in the ASPAR matrix ... ASPAR is shown in fig. 42, p.122
             I2    =  POSI(2,J)
             I3    =  POSI(3,J)
             IF (IOBP(IP) .NE. 2 .AND. IOBPD(IDD,IP) .EQ. 1) THEN
               ASPAR(I1) =  TRIA03 + DTK - TMP3 * DELTAL(POS             ,IE) + ASPAR(I1)  ! Diagonal entry
               ASPAR(I2) =               - TMP3 * DELTAL(POS_TRICK(POS,1),IE) + ASPAR(I2)  ! off diagonal entries ...
               ASPAR(I3) =               - TMP3 * DELTAL(POS_TRICK(POS,2),IE) + ASPAR(I3)
               B(IP)     =  B(IP) + TRIA03 * U(IP)
             END IF
           END DO !I: loop over connected elements ...
         END DO !IP

         DO IP = 1, MNP
           IF (IOBPD(IDD,IP) .EQ. 0) THEN
             ASPAR(I_DIAG(IP)) = SI(IP)
             B(IP)             = 0.d0
           END IF
         END DO

         IF (LBCWA .OR. LBCSP) THEN
           IF (LINHOM) THEN
             DO IP = 1, IWBMNP
               IPGL = IWBNDLC(IP)
               ASPAR(I_DIAG(IPGL)) = SI(IPGL) ! Add source term to the diagonal
               B(IPGL)             = SI(IPGL) * WBAC(ISS,IDD,IP)
             END DO
           ELSE
             DO IP = 1, IWBMNP
               IPGL = IWBNDLC(IP)
               ASPAR(I_DIAG(IPGL)) = SI(IPGL)
               B(IPGL)              = SI(IPGL) * WBAC(ISS,IDD,1)
             END DO
           ENDIF
         END IF

         IF (ICOMP .GE. 2 .AND. SMETHOD .GT. 0) THEN
           DO IP = 1, MNP
             IF (IOBP(IP) .EQ. 0) THEN
               ASPAR(I_DIAG(IP)) = ASPAR(I_DIAG(IP)) + IMATDAA(IP,ISS,IDD) * DT4A * SI(IP) * IOBPD(IDD,IP)  ! Add source term to the diagonal
               B(IP)             = B(IP)             + IMATRAA(IP,ISS,IDD) * DT4A * SI(IP) * IOBPD(IDD,IP) ! Add source term to the right hand side
             END IF
           END DO
         END IF


!> \todo for efficiency reason this will be moved out of this part when the code is working, AR.
! fill the matrix
         call MatZeroEntries(matrix, ierr);CHKERRQ(ierr)
!          call MatZeroEntries(matrixAppOrder, ierr);CHKERRQ(ierr)
         counter = 1
         nnznew=0

! iterate over all resident (and interface) nodes
! use the mapping to determine if the node is an interface node. if so, ignore it
! if it is an resident node, map it to petsc global ordering (this is also the row of the matrix)
! then get the connected nodes
! map the connected nodes to petsc global ordering (this is also the col of the matrix)
! and insert the value vom ASPAR into the matrix
         do i = 1, NP_RES
           ncols = IA(i+1) - IA(i)
           ! this is a interface node (row). ignore it. just increase counter
           if(ALOold2ALO(i-1) .eq. -999) then
             counter = counter + ncols
             cycle
           end if
           ! continue if this is a resident node
           ! map to petsc global order
           eRow = AGO2PGO(iplg(i) - 1)
           ! insert col by col into matrix
           do j = 1, ncols
             ! the value we want to insert
             eValue=ASPAR(counter)
             ! get the col index in old counting
             temp = JA(counter)
             counter = counter + 1
             nnznew=nnznew+1
             ! increase temp because iplg counts from one and JA from zero
             eCol = AGO2PGO(iplg((temp+1)) - 1 )
             call MatSetValue(matrix, eRow, eCol, eValue, ADD_VALUES, ierr);CHKERRQ(ierr)
             ! insert values in App global order
!              call MatSetValue(matrixAppOrder, iplg(i) - 1, iplg((temp+1)) - 1 , eValue, ADD_VALUES, ierr);CHKERRQ(ierr)
           end do
         end do ! np
         call MatAssemblyBegin(matrix,MAT_FINAL_ASSEMBLY, ierr);CHKERRQ(ierr);
         call MatAssemblyEnd(matrix,MAT_FINAL_ASSEMBLY, ierr);CHKERRQ(ierr);

!          call MatAssemblyBegin(matrixAppOrder,MAT_FINAL_ASSEMBLY, ierr);CHKERRQ(ierr);
!          call MatAssemblyEnd(matrixAppOrder,MAT_FINAL_ASSEMBLY, ierr);CHKERRQ(ierr);

! print some matrix values, just 4 debug. ignore it
!          DoPrint = 0
!           ! draw the matrix
!          IF (DoPrint.eq.1) THEN
!            call MatView(matrix, viewer, ierr ); CHKERRQ(ierr);CHKERRQ(ierr)
!            call PetscViewerDrawGetDraw(viewer, 0, draw, ierr);CHKERRQ(ierr)
!            call PetscDrawSetPause(draw,0, ierr);CHKERRQ(ierr);CHKERRQ(ierr)
!            call PetscDrawPause(draw, ierr);CHKERRQ(ierr)
!            ! print only one row
!          else if (DoPrint .eq. 2) then
!            if(rank == 1) then
!
! !              DO i = 1, nNodesWithoutInterfaceGhosts
! !                eCol = PLO2PGO(i-1)
!               do i = 1, np_global
!                 eCol = AGO2PGO(i-1)
! !              matrix, nRows, global row numner, number cols, global col number, value
!                CALL MatGetValues(matrix, 1, AGO2PGO(509), 1, eCol, eEntry, ierr); CHKERRQ(ierr);
! !                if(dabs(eEntry) > 1e-16)  write(*,*) "AGO value", PGO2AGO(ecol), eEntry
!                if(dabs(eEntry) > 1e-16)  write(*,*) PGO2AGO(ecol)
!              END DO
!            end if
!          end if

!print some matrix information

!fill RHS vector
!iterate over all resident (and interface) nodes
!map it to petsc global ordering
!and insert the value from B into RHS vector
         eEntry = 0;
         call VecSet(myB, eEntry, ierr);CHKERRQ(ierr)
!          call VecSet(myBAppOrder, eEntry, ierr);CHKERRQ(ierr)

         do i= 1, np
           ! this is a interface node (row). ignore it. just increase counter
           if(ALOold2ALO(i-1) .eq. -999) cycle 
           ! map to petsc global order
           eCol = AGO2PGO(iplg(i) - 1 )
           eEntry = B(i)
           call VecSetValue(myB, eCol, eEntry, ADD_VALUES, ierr);CHKERRQ(ierr)
!          call VecSetValue(myBAppOrder, iplg(i) - 1, eEntry, ADD_VALUES, ierr);CHKERRQ(ierr)
         end do

!> \todo it seems the ghost values must not be set by hand
!         do i = 1, nghost
!           eCol = i + nNodesWithoutInterfaceGhosts
!           eEntry = B(onlyGhostsOldLocalMapping(i-1))
!           call VecSetValue(myB, eCol, eEntry, ADD_VALUES, ierr);CHKERRQ(ierr)
!         end do

         call VecAssemblyBegin(myB, ierr);CHKERRQ(ierr)
         call VecAssemblyEnd(myB, ierr);CHKERRQ(ierr)

!          call VecAssemblyBegin(myBAppOrder, ierr);CHKERRQ(ierr);
!          call VecAssemblyEnd(myBAppOrder, ierr);CHKERRQ(ierr);
!          call VecGhostUpdateBegin(myB,INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr);
!          call VecGhostUpdateEnd(myB,INSERT_VALUES,SCATTER_FORWARD,ierr);CHKERRQ(ierr);



!Calc norms for myB
!          call VecNorm(myB, NORM_1, eNormB1, ierr);CHKERRQ(ierr);
!          call VecNorm(myB, NORM_2, eNormB2, ierr);CHKERRQ(ierr);
!          call VecNorm(myB, NORM_INFINITY, eNormBi, ierr);CHKERRQ(ierr);
!          Print *, rank, ' L1_B=', eNormB1
!          Print *, rank, ' L2_B=', eNormB2
!          Print *, rank, ' Li_B=', eNormBi



! print matrix myX and myB
!         if(SUM(U) > 10e-2) then
!         if(ISS == 2 .and. IDD == 2) then
!           call printPDE(ISS, IDD)
!           if(timePetscCalled == 5) stop
!         endif


!Solve
         ! Copy the old solution from AC2 to myX to make the solver faster
         do i = 1, np
          eCol = AGO2PGO(iplg(i)-1)
          eEntry = DBLE(AC2(i, ISS, IDD))
          call VecSetValue(myX, eCol, eEntry, INSERT_VALUES, ierr);CHKERRQ(ierr)
         end do
         call VecAssemblyBegin(myX, ierr);CHKERRQ(ierr)
         call VecAssemblyEnd(myX, ierr);CHKERRQ(ierr)

         ! Solve!
         call KSPSolve(Solver, myB, myX, ierr);CHKERRQ(ierr)
         call KSPGetConvergedReason(Solver, reason, ierr);CHKERRQ(ierr)
         if (reason .LT. 0) then
           write(stat%fhndl,*) "Failure to converge\n"
         else
          call KSPGetIterationNumber(Solver, iterationen, ierr);CHKERRQ(ierr)
         !   !if(iterationen /= 0)  write(*,*) "Number of iterations", iss,idd,iterationen
         !   write(stat%fhndl,*) "Number of iterations", iterationen
         endif

! print matrix myX and myB
!         if(SUM(U) > 10e-2) then
!         call printPDE(ISS, IDD)
!         endif

!Calc norms for myX
!          call VecNorm(myX, NORM_1, eNormX1, ierr);CHKERRQ(ierr)
!          call VecNorm(myX, NORM_2, eNormX2, ierr);CHKERRQ(ierr)
!          call VecNorm(myX, NORM_INFINITY, eNormXi, ierr);CHKERRQ(ierr)
!          Print *, rank, ' L1_X=', eNormX1
!          Print *, rank, ' L2_X=', eNormX2
!          Print *, rank, ' Li_X=', eNormXi

!          call VecAssemblyBegin(myX, ierr);CHKERRQ(ierr);
!          call VecAssemblyEnd(myX, ierr);CHKERRQ(ierr);
!          call VecGhostUpdateBegin(myX, INSERT_VALUES, SCATTER_FORWARD, ierr);CHKERRQ(ierr);
!          call VecGhostUpdateEnd(myX, INSERT_VALUES, SCATTER_FORWARD, ierr);CHKERRQ(ierr);
      X = 0.d0!DBLE(AC2(:, ISS, IDD))
      ! get the solution back to fortran. version 2
      ! create a scatter that collects from all threads
      ! i think this is the same as call EXCHANGE_P2D(X)
!        call VecScatterCreateToAll(myX, ctx, vout, ierr); CHKERRQ(ierr)
!        call VecScatterBegin(ctx, myX, vout, INSERT_VALUES, SCATTER_FORWARD, ierr); CHKERRQ(ierr)
!        call VecScatterEnd(ctx, myX, vout, INSERT_VALUES, SCATTER_FORWARD,ierr); CHKERRQ(ierr)
!        call VecScatterDestroy(ctx, ierr); CHKERRQ(ierr)

       ! copy the solution back to the fortran vector
!       call VecGetArrayF90(vout, myXtemp, IERR); CHKERRQ(IERR)
!       do i=1, npa
!         X(i) = myXtemp(AGO2PGO( iplg(i) -1) +1)
!       end do
!       call VecRestoreArrayF90(vout, myXtemp, ierr); CHKERRQ(ierr);
!       call VecDestroy(vout, ierr);; CHKERRQ(ierr)

       !get the solution back to fortran. version 1
       !iterate over all resident nodes (without interface and ghost nodes)
       !map the solution from petsc local ordering back to app old local ordering
       !(the app old ordering contains interface nodes)
         call VecGetArrayF90(myX, myXtemp, IERR); CHKERRQ(IERR)
         do i = 1, nNodesWithoutInterfaceGhosts
           X(ipgl1((PGO2AGO(PLO2PGO(i-1)))+1)%id) = myXtemp(i)
         end do
         call VecRestoreArrayF90(myX, myXtemp, ierr); CHKERRQ(ierr);
         !IF (SUM(X) .NE. SUM(X)) STOP 'NaN in X'
         AC2(:, ISS, IDD) = SNGL(MAX(0.d0,X))

      END SUBROUTINE
#endif
    END MODULE
#endif
