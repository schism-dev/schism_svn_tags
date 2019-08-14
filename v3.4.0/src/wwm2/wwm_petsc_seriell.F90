#include "wwm_functions.h"
#ifdef PETSC
      MODULE PETSC_SERIELL

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


      KSPConvergedReason reason;
      PetscInt iterationen;

      contains

      SUBROUTINE PETSC_INIT_SERIELL()

        USE DATAPOOL
        use petscpool
        IMPLICIT NONE

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
      SUBROUTINE  EIMPS_PETSC_SERIELL(ISS, IDD)

         USE DATAPOOL
         use petscpool
         !use elfe_msgp
         IMPLICIT NONE

         include 'mpif.h'

         INTEGER, INTENT(IN) :: ISS,IDD

         INTEGER :: I, J
!          INTEGER :: KK

         INTEGER :: IP, IPGL, IE, POS

         INTEGER :: I1, I2, I3, IPrel

         real(rkind) :: DTK, TMP3

         real(rkind) :: LAMBDA(2)
         real(rkind) :: FL11, FL12, FL21, FL22, FL31, FL32
         real(rkind) :: CRFS(3), K1, KM(3), K(3), TRIA03

         real(rkind) :: DELTAL(3,MNE)
         real(rkind) :: KP(3,MNE), NM(MNE)
         real(rkind) :: U(MNP), C(2,MNP)

         real(rkind) :: X(MNP)
         real(rkind) :: B(MNP)

         real(rkind) ::  ASPAR(NNZ)

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
           KP(1,IE) = MAX(ZERO,K(1))
           KP(2,IE) = MAX(ZERO,K(2))
           KP(3,IE) = MAX(ZERO,K(3))
           KM(1) = MIN(ZERO,K(1))
           KM(2) = MIN(ZERO,K(2))
           KM(3) = MIN(ZERO,K(3))
           FL11 = C(1,I2)*IEN(1,IE)+C(2,I2)*IEN(2,IE)
           FL12 = C(1,I3)*IEN(1,IE)+C(2,I3)*IEN(2,IE)
           FL21 = C(1,I3)*IEN(3,IE)+C(2,I3)*IEN(4,IE)
           FL22 = C(1,I1)*IEN(3,IE)+C(2,I1)*IEN(4,IE)
           FL31 = C(1,I1)*IEN(5,IE)+C(2,I1)*IEN(6,IE)
           FL32 = C(1,I2)*IEN(5,IE)+C(2,I2)*IEN(6,IE)
           CRFS(1) =  - ONESIXTH *  (TWO *FL31 + FL32 + FL21 + TWO * FL22 )
           CRFS(2) =  - ONESIXTH *  (TWO *FL32 + TWO * FL11 + FL12 + FL31 )
           CRFS(3) =  - ONESIXTH *  (TWO *FL12 + TWO * FL21 + FL22 + FL11 )
           DELTAL(:,IE) = CRFS(:)- KP(:,IE)
           NM(IE)       = 1.0_rkind/MIN(-THR,SUM(KM(:)))
         END DO

         CALL CPU_TIME(TIME2)

         U = AC2(:,ISS,IDD)

         J     = 0    ! Counter ...
         ASPAR = ZERO ! Mass matrix ...
         B     = ZERO ! Right hand side ...
!
! ... assembling the linear equation system ....
!
         DO IP = 1, NP_RES
           IF (IOBPD(IDD,IP) .EQ. 1 .AND. IOBWB(IP) .EQ. 1 .AND. DEP(IP) .GT. DMIN) THEN
             DO I = 1, CCON(IP)
               J = J + 1
               IE    =  IE_CELL(J)
               POS   =  POS_CELL(J)
               K1    =  KP(POS,IE) ! Flux Jacobian
               TRIA03 = ONETHIRD * TRIA(IE)
               DTK   =  K1 * DT4A * IOBPD(IDD,IP)
               TMP3  =  DTK * NM(IE)
               I1    =  POSI(1,J) ! Position of the recent entry in the ASPAR matrix ... ASPAR is shown in fig. 42, p.122
               I2    =  POSI(2,J)
               I3    =  POSI(3,J)
               ASPAR(I1) =  TRIA03 + DTK - TMP3 * DELTAL(POS             ,IE) + ASPAR(I1)  ! Diagonal entry
               ASPAR(I2) =               - TMP3 * DELTAL(POS_TRICK(POS,1),IE) + ASPAR(I2)  ! off diagonal entries ...
               ASPAR(I3) =               - TMP3 * DELTAL(POS_TRICK(POS,2),IE) + ASPAR(I3)
               B(IP)     =  B(IP) + TRIA03 * U(IP)
             END DO !I: loop over connected elements ...
           ELSE
             DO I = 1, CCON(IP)
               J = J + 1
               IE    =  IE_CELL(J)
               TRIA03 = ONETHIRD * TRIA(IE)
               I1    =  POSI(1,J) ! Position of the recent entry in the ASPAR matrix ... ASPAR is shown in fig. 42, p.122
               ASPAR(I1) =  TRIA03 + ASPAR(I1)  ! Diagonal entry
               B(IP)     =  0.!B(IP)  + TRIA03 * 0.
             END DO !I: loop over connected elements ...
           END IF
         END DO !IP

         IF (LBCWA .OR. LBCSP) THEN
           DO IP = 1, IWBMNP
             IF (LINHOM) THEN
               IPrel=IP
             ELSE
               IPrel=1
             ENDIF
             IPGL = IWBNDLC(IP)
             ASPAR(I_DIAG(IPGL)) = SI(IPGL) ! Set boundary on the diagonal
             B(IPGL)             = SI(IPGL) *  WBAC(ISS,IDD,IPrel)
           END DO
         END IF

         IF (ICOMP .GE. 2 .AND. SMETHOD .GT. 0) THEN
           DO IP = 1, NP_RES
             IF (IOBWB(IP) .EQ. 1) THEN
               ASPAR(I_DIAG(IP)) = ASPAR(I_DIAG(IP)) + IMATDAA(IP,ISS,IDD) * DT4A * SI(IP) ! Add source term to the diagonal
               B(IP)             = B(IP) + IMATRAA(IP,ISS,IDD) * DT4A * SI(IP) ! Add source term to the right hand side
             ENDIF
           END DO
         ENDIF

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
           write(*,*) "Failure to converge\n"
          stop 'wwm_petsc_seriell l.260'
         else
         call KSPGetIterationNumber(Solver, iterationen, ierr);CHKERRQ(ierr)
            if(iterationen /= 0)  write(*,*) "Number of iterations", iss,idd,iterationen
!            write(*,*) "Number of iterations", iterationen
         endif

         ! write the soluten from vec into fortran array
         X    = ZERO
         call VecGetArrayF90(myX, myXtemp, ierr);CHKERRQ(ierr)
         X = myXtemp
         call VecRestoreArrayF90(myX, myXtemp, ierr);CHKERRQ(ierr)

         IF (SUM(X) .NE. SUM(X)) CALL WWM_ABORT('NaN in X')
!          AC2(:, ISS, IDD) = MAX(ZERO,X)

         AC2(:,ISS,IDD) = X

!          CALL CPU_TIME(TIME4)

         RETURN
        END SUBROUTINE

      !> cleanup memory. You never need to call this function by hand. It will automaticly called by PETSC_FINALIZE()
      subroutine PETSC_FINALIZE_SERIELL()
        implicit none

      end subroutine
end module
#endif
