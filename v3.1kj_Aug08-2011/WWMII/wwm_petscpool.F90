#ifdef PETSC
      MODULE PETSCPOOL
      
#include "finclude/petscsys.h"
#include "finclude/petscvec.h"
#include "finclude/petscvec.h90"
#include "finclude/petscmat.h"
#include "finclude/petscksp.h"
               
         ! petsc stuff

         Vec :: myB;
         PetscScalar, pointer :: myBtemp(:)
         Vec :: myX
         PetscScalar, pointer :: myXtemp(:)
             
         Mat matrix;
         PetscErrorCode ierr;
         KSP Solver;
         PC Prec;
         KSPConvergedReason reason;
         PetscInt iterationen;
         
         contains

      SUBROUTINE PETSC_INIT()
      
        USE DATAPOOL
       
        IMPLICIT NONE
        
        ! petsc wants indices startet from 0
        IA = IA -1
        JA = JA -1
         
         call PetscInitialize(PETSC_NULL_CHARACTER,ierr) 
         CHKERRQ(ierr);
        
         ! set up the solver
         call KSPCreate(PETSC_COMM_SELF, Solver, ierr)
         CHKERRQ(ierr);

         !Create preconditioner
         call KSPGetPC(Solver, Prec, ierr)
         CHKERRQ(ierr);
         call PCSetType(Prec, PCILU, ierr)
         CHKERRQ(ierr);
        
         !create x vector  
         call VecCreate(PETSC_COMM_SELF, myX, ierr);
         CHKERRQ(ierr);
         call VecSetSizes(myX, PETSC_DECIDE, MNP, ierr);
         CHKERRQ(ierr);
         call VecSetType(myX,VECSEQ, ierr);
         CHKERRQ(ierr);
        
         ! create vec myB
         call VecCreate(PETSC_COMM_SELF, myB, ierr)
         CHKERRQ(ierr);
         call VecSetSizes(myB, PETSC_DECIDE, MNP, ierr);
         CHKERRQ(ierr);
         call VecSetType(myB,VECSEQ, ierr);
         CHKERRQ(ierr);
        
        
      END SUBROUTINE
            
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE  EIMPS_PETSC(ISS, IDD)

         USE DATAPOOL
#ifdef SELFE
         use elfe_glbl, only: iplg
         use elfe_msgp
#endif
         IMPLICIT NONE

#ifdef SELFE
         include 'mpif.h'
#endif

         INTEGER   :: ISS,IDD

         INTEGER :: I, J, KK

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

   
         real*8, ALLOCATABLE :: temp(:);


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
 
         ! create sparse matrix
         call MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, MNP, MNP, IA, JA, ASPAR, matrix, ierr);
         CHKERRQ(ierr)

         ! fill myB
         call VecGetArrayF90(myB, myBtemp, ierr)
         CHKERRQ(ierr)

         myBtemp = B

         call VecRestoreArrayF90(myB, myBtemp, ierr)
         CHKERRQ(ierr)
  
         call KSPSetOperators(Solver, matrix, matrix, SAME_NONZERO_PATTERN, ierr)
         CHKERRQ(ierr)

!         call KSPSetType(Solver,KSPBCGS, ierr)
         call KSPSetType(Solver, KSPIBCGS, ierr)
         CHKERRQ(ierr)

         CALL CPU_TIME(TIME4)
                 
         call KSPSolve(Solver, myB, myX, ierr)
         CHKERRQ(ierr)

         CALL CPU_TIME(TIME5)

         !WRITE(*,*) TIME2-TIME1, TIME3-TIME2, TIME4-TIME3, TIME5-TIME4
      
         call KSPGetConvergedReason(Solver, reason, ierr)
         CHKERRQ(ierr)

         if (reason .LT. 0) then
           write(*,*) "Failure to converge\n"
           stop
         else 
         call KSPGetIterationNumber(Solver, iterationen, ierr)
           CHKERRQ(ierr)
           write(*,*) "Number of iterations", iterationen
         endif
  
         ! write the soluten from vec into fortran array
         X    = 0.d0
         call VecGetArrayF90(myX, myXtemp, ierr)
         CHKERRQ(ierr)
         X = myXtemp
         call VecRestoreArrayF90(myX, myXtemp, ierr)
         CHKERRQ(ierr)
       
         AC2(:,ISS,IDD) = X

!          CALL CPU_TIME(TIME4)

         RETURN
        END SUBROUTINE
             
      END MODULE
#endif PETSC
