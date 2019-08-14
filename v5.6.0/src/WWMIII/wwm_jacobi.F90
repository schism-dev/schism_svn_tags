#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
!AR:todo: code duplication ... and buggy since boundary pointer are not taken into account 
      SUBROUTINE COMPUTE_CFL_N_SCHEME_EXPLICIT(CFLadvgeoOutI)
      USE DATAPOOL
      IMPLICIT NONE
      integer,intent(out) :: CFLadvgeoOutI(MNP)
      REAL(rkind) C(2,MNP), KELEM(3,MNE)
      REAL(rkind) KKSUM(MNP)
      REAL(rkind) LAMBDA(2), KTMP(3)
      integer I1, I2, I3
      real(rkind) DTMAX_EXP
      integer IE, IP, I, J, POS
      integer IS, ID
      REAL(rkind) :: CFLadvgeoOut(MNP)
      CFLadvgeoOut=0
      DO IS=1,NUMSIG
        DO ID=1,NUMDIR
          CALL CADVXY(IS,ID,C)
          DO IE = 1, MNE
            I1 = INE(1,IE)
            I2 = INE(2,IE)
            I3 = INE(3,IE)
            LAMBDA(1) = ONESIXTH *(C(1,I1)+C(1,I2)+C(1,I3))
            LAMBDA(2) = ONESIXTH *(C(2,I1)+C(2,I2)+C(2,I3))
            KELEM(1,IE) = LAMBDA(1) * IEN(1,IE) + LAMBDA(2) * IEN(2,IE)
            KELEM(2,IE) = LAMBDA(1) * IEN(3,IE) + LAMBDA(2) * IEN(4,IE)
            KELEM(3,IE) = LAMBDA(1) * IEN(5,IE) + LAMBDA(2) * IEN(6,IE)
            KTMP  = KELEM(:,IE)
            KELEM(:,IE) = MAX(ZERO,KTMP)
          END DO
          KKSUM = ZERO
          J    = 0
          DO IP = 1, MNP
            DO I = 1, CCON(IP)
              J = J + 1
              IE    = IE_CELL(J)
              POS   = POS_CELL(J)
              KKSUM(IP)  = KKSUM(IP) + MAX(KELEM(POS,IE),ZERO)
            END DO
          END DO
          DO IP=1,MNP
            DTMAX_EXP = SI(IP)/MAX(THR,KKSUM(IP))
            CFLadvgeoOut(IP) = MAX(CFLadvgeoOut(IP), DT4A / DTMAX_EXP)
          END DO
        END DO
      END DO
#ifdef MPI_PARALL_GRID     
      CALL EXCHANGE_P2D(CFLadvgeoOut)
#endif
      DO IP=1,MNP
        CFLadvgeoOutI(IP) = NINT(CFLadvgeoOut(IP))
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE EIMPS_ASPAR_BLOCK(ASPAR)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), intent(out) :: ASPAR(NUMSIG, NUMDIR, NNZ)
      INTEGER :: POS_TRICK(3,2)
      REAL(rkind) :: FL11(NUMSIG,NUMDIR), FL12(NUMSIG,NUMDIR), FL21(NUMSIG,NUMDIR), FL22(NUMSIG,NUMDIR), FL31(NUMSIG,NUMDIR), FL32(NUMSIG,NUMDIR)
      REAL(rkind) :: CRFS(NUMSIG,NUMDIR,3), K1(NUMSIG,NUMDIR), KM(NUMSIG,NUMDIR,3), K(NUMSIG,NUMDIR,3), TRIA03
      REAL(rkind) :: CXY(2,NUMSIG,NUMDIR,3)
      REAL(rkind) :: DIFRU, USOC, WVC
      REAL(rkind) :: DELTAL(NUMSIG,NUMDIR,3)
      REAL(rkind) :: KP(NUMSIG,NUMDIR,3), NM(NUMSIG,NUMDIR)
      INTEGER     :: I1, I2, I3
      INTEGER     :: IP, ID, IS, IE
      INTEGER     :: I, IPGL1

      REAL(rkind) :: DTK(NUMSIG,NUMDIR), TMP3(NUMSIG,NUMDIR)
      REAL(rkind) :: LAMBDA(2,NUMSIG,NUMDIR)
      REAL(rkind) :: CXnorm
      REAL(rkind) sumTMP3, sumCXY, sumCG, sumASPARdiag
      REAL(rkind) DTeffect
      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2
!
!     Calculate countour integral quantities ...
!
#ifdef DEBUG
      WRITE(STAT%FHNDL,*) 'MNP=', MNP
      WRITE(STAT%FHNDL,*) 'sum(IOBWB)=', sum(IOBWB)
      WRITE(STAT%FHNDL,*) 'sum(IOBPD)=', sum(IOBPD)
      WRITE(STAT%FHNDL,*) 'sum(IOBDP)=', sum(IOBDP)
      WRITE(STAT%FHNDL,*) 'LSPHE=', LSPHE
      sumTMP3=0
      sumCXY=0
      sumCG=0
#endif
      IF (IMPL_GEOADVECT) THEN
        DTeffect = DT4A
      ELSE
        DTeffect = ZERO
      END IF
      IF (LCFL) THEN
        CFLCXY(1,:) = ZERO
        CFLCXY(2,:) = ZERO
        CFLCXY(3,:) = LARGE
      END IF
      ASPAR = 0.0_rkind ! Mass matrix ...
      DO IE = 1, MNE
        DO I=1,3
          IP = INE(I,IE)
          DO IS=1,NUMSIG
            DO ID=1,NUMDIR
              IF (LSECU .OR. LSTCU) THEN
                CXY(1,IS,ID,I) = CG(IS,IP)*COSTH(ID)+CURTXY(IP,1)
                CXY(2,IS,ID,I) = CG(IS,IP)*SINTH(ID)+CURTXY(IP,2)
              ELSE
                CXY(1,IS,ID,I) = CG(IS,IP)*COSTH(ID)
                CXY(2,IS,ID,I) = CG(IS,IP)*SINTH(ID)
              END IF
              IF (LSPHE) THEN
                CXY(1,IS,ID,I) = CXY(1,IS,ID,I)*INVSPHTRANS(IP,1)
                CXY(2,IS,ID,I) = CXY(2,IS,ID,I)*INVSPHTRANS(IP,2)
              END IF
              IF (LDIFR) THEN
                CXY(1,IS,ID,I) = CXY(1,IS,ID,I)*DIFRM(IP)
                CXY(2,IS,ID,I) = CXY(2,IS,ID,I)*DIFRM(IP)
                IF (LSECU .OR. LSTCU) THEN
                  IF (IDIFFR .GT. 1) THEN
                    WVC = SPSIG(IS)/WK(IS,IP)
                    USOC = (COSTH(ID)*CURTXY(IP,1) + SINTH(ID)*CURTXY(IP,2))/WVC
                    DIFRU = ONE + USOC * (ONE - DIFRM(IP))
                  ELSE
                    DIFRU = DIFRM(IP)
                  END IF
                  CXY(1,IS,ID,I) = CXY(1,IS,ID,I) + DIFRU*CURTXY(IP,1)
                  CXY(2,IS,ID,I) = CXY(2,IS,ID,I) + DIFRU*CURTXY(IP,2)
                END IF
              END IF
              IF (LCFL) THEN
                CXnorm=SQRT(CXY(1,IS,ID,I)**2 + CXY(2,IS,ID,I)**2)
                CFLCXY(1,IP) = MAX(CFLCXY(1,IP), CXnorm)
              END IF
            END DO
          END DO
        END DO
#ifdef DEBUG
        sumCXY = sumCXY + sum(abs(CXY))
        sumCG  = sumCG  + sum(abs(CG ))
#endif
        LAMBDA(:,:,:) = ONESIXTH * (CXY(:,:,:,1) + CXY(:,:,:,2) + CXY(:,:,:,3))
        K(:,:,1)  = LAMBDA(1,:,:) * IEN(1,IE) + LAMBDA(2,:,:) * IEN(2,IE)
        K(:,:,2)  = LAMBDA(1,:,:) * IEN(3,IE) + LAMBDA(2,:,:) * IEN(4,IE)
        K(:,:,3)  = LAMBDA(1,:,:) * IEN(5,IE) + LAMBDA(2,:,:) * IEN(6,IE)
        FL11(:,:) = CXY(1,:,:,2)*IEN(1,IE)+CXY(2,:,:,2)*IEN(2,IE)
        FL12(:,:) = CXY(1,:,:,3)*IEN(1,IE)+CXY(2,:,:,3)*IEN(2,IE)
        FL21(:,:) = CXY(1,:,:,3)*IEN(3,IE)+CXY(2,:,:,3)*IEN(4,IE)
        FL22(:,:) = CXY(1,:,:,1)*IEN(3,IE)+CXY(2,:,:,1)*IEN(4,IE)
        FL31(:,:) = CXY(1,:,:,1)*IEN(5,IE)+CXY(2,:,:,1)*IEN(6,IE)
        FL32(:,:) = CXY(1,:,:,2)*IEN(5,IE)+CXY(2,:,:,2)*IEN(6,IE)
        CRFS(:,:,1) = - ONESIXTH *  (TWO *FL31(:,:) + FL32(:,:) + FL21(:,:) + TWO * FL22(:,:) )
        CRFS(:,:,2) = - ONESIXTH *  (TWO *FL32(:,:) + TWO * FL11(:,:) + FL12(:,:) + FL31(:,:) )
        CRFS(:,:,3) = - ONESIXTH *  (TWO *FL12(:,:) + TWO * FL21(:,:) + FL22(:,:) + FL11(:,:) )
        KM = MIN(ZERO,K)
        KP(:,:,:) = MAX(ZERO,K)
        DELTAL(:,:,:) = CRFS(:,:,:)- KP(:,:,:)
        NM(:,:)=ONE/MIN(-THR,KM(:,:,1) + KM(:,:,2) + KM(:,:,3))
        TRIA03 = ONETHIRD * TRIA(IE)
        DO I=1,3
          IP=INE(I,IE)
          I1=JA_IE(I,1,IE)
          I2=JA_IE(I,2,IE)
          I3=JA_IE(I,3,IE)
          K1(:,:) =  KP(:,:,I)
          DO ID=1,NUMDIR
            DTK(:,ID) =  K1(:,ID) * DTeffect * IOBPD(ID,IP) * IOBWB(IP) * IOBDP(IP)
          END DO
          TMP3(:,:)  =  DTK(:,:) * NM(:,:)
#ifdef DEBUG
          sumTMP3 = sumTMP3 + sum(abs(TMP3))
#endif
          ASPAR(:,:,I1) =  TRIA03+DTK(:,:)- TMP3(:,:) * DELTAL(:,:,I             ) + ASPAR(:,:,I1)
          ASPAR(:,:,I2) =                 - TMP3(:,:) * DELTAL(:,:,POS_TRICK(I,1)) + ASPAR(:,:,I2)
          ASPAR(:,:,I3) =                 - TMP3(:,:) * DELTAL(:,:,POS_TRICK(I,2)) + ASPAR(:,:,I3)
        END DO
      END DO
#ifdef DEBUG
      WRITE(STAT%FHNDL,*) 'sumTMP3=', sumTMP3
      WRITE(STAT%FHNDL,*) 'sumCXY=', sumCXY
      WRITE(STAT%FHNDL,*) 'sumCG=', sumCG
      WRITE(STAT%FHNDL,*) 'LBCWA=', LBCWA
      WRITE(STAT%FHNDL,*) 'LBCSP=', LBCSP
      WRITE(STAT%FHNDL,*) 'IWBMNP=', IWBMNP
#endif
      IF (LBCWA .OR. LBCSP) THEN
        DO IP = 1, IWBMNP
          IPGL1 = IWBNDLC(IP)
          ASPAR(:,:,I_DIAG(IPGL1)) = SI(IPGL1) ! Set boundary on the diagonal
        END DO
      END IF
#ifdef DEBUG
      sumASPARdiag=0
      DO IP=1,MNP
        sumASPARdiag = sumASPARdiag + sum(abs(ASPAR(:,:,I_DIAG(IP))))
      END DO
      WRITE(STAT%FHNDL,*) 'sumASPARdiag=', sumASPARdiag
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ADD_FREQ_DIR_TO_ASPAR_COMP_CADS(ASPAR)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), intent(inout) :: ASPAR(NUMSIG,NUMDIR,NNZ)
      REAL(rkind) :: TheVal, eFact
      REAL(rkind) :: CASS(0:NUMSIG+1), CP_SIG(0:NUMSIG+1), CM_SIG(0:NUMSIG+1)
      REAL(rkind) :: CAD(NUMSIG,NUMDIR), CAS(NUMSIG,NUMDIR)
      REAL(rkind) :: CP_THE(NUMSIG,NUMDIR), CM_THE(NUMSIG,NUMDIR)
      REAL(rkind) :: B_SIG(NUMSIG)
      INTEGER     :: IS, ID, IP
      IF (REFRACTION_IMPL) THEN
        DO IP=1,NP_RES
          TheVal=1
          IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LTHBOUND) TheVal=0
          IF (DEP(IP) .LT. DMIN) TheVal=0
          IF (IOBP(IP) .EQ. 2) TheVal=0
          IF (TheVal .eq. 1) THEN
            CALL PROPTHETA(IP,CAD)
          ELSE
            CAD=ZERO
          END IF
          CP_THE = MAX(ZERO,CAD)
          CM_THE = MIN(ZERO,CAD)
          eFact=(DT4D/DDIR)*SI(IP)
          CAD_THE(:,:,IP)=CAD
          ASPAR(:,:,I_DIAG(IP)) = ASPAR(:,:,I_DIAG(IP)) + eFact * (CP_THE(:,:) - CM_THE(:,:))
        END DO
      END IF
      IF (FREQ_SHIFT_IMPL) THEN
        DO IP=1,NP_RES
          TheVal=1
          IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LSIGBOUND) TheVal=0
          IF (DEP(IP) .LT. DMIN) TheVal=0
          IF (IOBP(IP) .EQ. 2) TheVal=0
          IF (TheVal .eq. 1) THEN
            CALL PROPSIGMA(IP,CAS)
          ELSE
            CAS=ZERO
          END IF
          CAS_SIG(:,:,IP)=CAS
          eFact=DT4F*SI(IP)
          DO ID = 1, NUMDIR
            CASS(1:NUMSIG) = CAS(:,ID)
            CASS(0)     = 0.
            CASS(NUMSIG+1) = CASS(NUMSIG)
            CP_SIG = MAX(ZERO,CASS)
            CM_SIG = MIN(ZERO,CASS)
            ! Now forming the tridiagonal system
            DO IS=1,NUMSIG
              B_SIG(IS)=eFact*(CP_SIG(IS)/DS_INCR(IS-1) - CM_SIG(IS) /DS_INCR(IS))
            END DO
            B_SIG(NUMSIG) = B_SIG(NUMSIG) + eFact*CM_SIG(NUMSIG+1)/DS_INCR(NUMSIG) * TAIL_ARR(5)
            ASPAR(:,ID,I_DIAG(IP))=ASPAR(:,ID,I_DIAG(IP)) + B_SIG
          END DO
        END DO
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!* For the refraction, we use the Upwind implicit scheme
!* N^{n+1} = N^n + f(N^(n+1))
!* This solves the differential equation N'=f(N)
!*
!* Courant, R., Isaacson, E., and Rees, M. (1952). "On the Solution of
!* Nonlinear Hyperbolic Differential Equations by Finite Differences",
!* Comm. Pure Appl. Math., 5, 243–255.
!*
!* This gives
!*
!* (1) N_i^(n+1) = N_i^n + (delta t/delta x)
!*            (u_(i+1)^n N_(i+1)^(n+1) - u_i^n N_i^(n+1)) if u_i^n > 0
!* (2) N_i^(n+1) = N_i^n + (delta t/delta x)
!*            (u_i^n N_i^(n+1) - u_(i-1)^n N_(i-1)^(n+1)) if u_i^n < 0
!* 
!* The notations for tridiagonal system are available from
!* http://en.wikipedia.org/wiki/Tridiagonal_matrix
!*
!* For frequency shifting, things are complicated:
!* 
!* Boundary condition: For low frequency, energy disappear. For
!* high frequency, we prolongate the energy by using a parametrization
!* of the tail: TAIL_ARR(5).
!* 
!* Grid: the gridsize is variable. DS_INCR(IS) is essentially defined
!* as  DS_INCR(IS) = SPSIG(IS) - SPSIG(IS-1)
!* Therefore the system that needs to be resolved is for i=1,NUMSIG
!*
!* We write f_{n,+} = 1 if u_i^n > 0
!*                    0 otherwise
!* We write f_{n,-} = 0 if u_i^n > 0
!*                    1 otherwise
!* We write u_{i,n,+} = u_i f_{n,+} and similarly for other variables.
!* 
!* If we continue like that then we eventually get a non-conservative
!* scheme. See below for details.
!* 
!* N_i^(n+1) = N_i^n + (delta t) [
!* + { u_(i+1,n,+) N_(i+1)^(n+1) - u_(i,n,+)   N_i^(n+1)     }/DS_INCR_i+1
!*   { u_(i,n,-)   N_i^(n+1)     - u_(i-1,n,-) N_(i-1)^(n+1) }/DS_INCR_i
!* 
!* which after rewrites give us
!* N_i^n = N_i^(n+1)     [1 + Delta t { u_(i,n,+)/DS_i+1  
!*                                  -   u_(i,n,-)/DS_i          }    ]
!*       + N_(i-1)^(n+1) [    Delta t {  u_(i-1,n,-)/DS_i       }    ]
!*       + N_(i+1)^(n+1) [    Delta t { -u_(i+1,n,+)/DS_(i+1)   }    ]
!*
!* Instead, we set u_{i,n,+} = u_{i,+} = max(u_i, 0)
!*                 u_{i,n,-} = u_{i,-} = min(u_i, 0)
!* and the equations become simpler:
!* N_i^n = N_i^(n+1)     [1 + Delta t { u_(i,+)/DS_i+1  
!*                                  -   u_(i,-)/DS_i          }    ]
!*       + N_(i-1)^(n+1) [    Delta t {  u_(i-1,-)/DS_i       }    ]
!*       + N_(i+1)^(n+1) [    Delta t { -u_(i+1,+)/DS_(i+1)   }    ]
!* 
!* The boundary conditions are expressed as
!* N_0^{n+1}=0 and N_{NUMSIG+1}^{n+1} = N_{NUMSIG}^{n+1} TAIL_ARR(5)
!* 
!* 
!**********************************************************************
      SUBROUTINE GET_FREQ_DIR_CONTRIBUTION(IP, ASPAR_DIAG, A_THE, C_THE, A_SIG, C_SIG)
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind), intent(inout) :: ASPAR_DIAG(NUMSIG,NUMDIR)
      REAL(rkind), intent(out) :: A_THE(NUMSIG,NUMDIR), C_THE(NUMSIG,NUMDIR)
      REAL(rkind), intent(out) :: A_SIG(NUMSIG,NUMDIR), C_SIG(NUMSIG,NUMDIR)

      REAL(rkind) :: TheVal, eFact
      REAL(rkind) :: CP_SIG, CM_SIG
      REAL(rkind) :: CP_SIG_ip1, CM_SIG_im1
      REAL(rkind) :: FP, FM
      REAL(rkind) :: CAD(NUMSIG,NUMDIR), CAS(NUMSIG,NUMDIR)
      REAL(rkind) :: CP_THE(NUMSIG,NUMDIR), CM_THE(NUMSIG,NUMDIR)
      REAL(rkind) :: B_SIG(NUMSIG)
      INTEGER     :: ID1, ID2, IS, ID, IP
      IF (REFRACTION_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LTHBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPTHETA(IP,CAD)
        ELSE
          CAD=ZERO
        END IF
        CP_THE = MAX(ZERO,CAD)
        CM_THE = MIN(ZERO,CAD)
        eFact=(DT4D/DDIR)*SI(IP)
        DO ID=1,NUMDIR
          ID1 = ID_PREV(ID)
          ID2 = ID_NEXT(ID)
          A_THE(:,ID) = - eFact *  CP_THE(:,ID1)
          C_THE(:,ID) =   eFact *  CM_THE(:,ID2)
        END DO
        ASPAR_DIAG=ASPAR_DIAG + eFact * (CP_THE(:,:) - CM_THE(:,:))
      END IF
      IF (FREQ_SHIFT_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LSIGBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPSIGMA(IP,CAS)
        ELSE
          CAS=ZERO
        END IF
        eFact=DT4F*SI(IP)
        DO ID = 1, NUMDIR
          DO IS=1,NUMSIG
            IF (CAS(IS,ID) .gt. 0) THEN
              FP=ONE
              FM=ZERO
            ELSE
              FP=ZERO
              FM=ONE
            END IF
            CP_SIG=CAS(IS,ID) * FP
            CM_SIG=CAS(IS,ID) * FM
            B_SIG(IS)=eFact*(CP_SIG/DS_INCR(IS+1) - CM_SIG/DS_INCR(IS))
            IF (IS .eq. NUMSIG) THEN
              CP_SIG_ip1=CAS(NUMSIG,ID)*FP*TAIL_ARR(5)
              B_SIG(NUMSIG)=B_SIG(NUMSIG) - eFact*CP_SIG_ip1/DS_INCR(IS)
            END IF
            IF (IS .gt. 1) THEN
              CM_SIG_im1=CAS(IS-1,ID)*FM
              A_SIG(IS,ID)=eFact*CM_SIG_im1/DS_INCR(IS)
            END IF
            IF (IS .lt. NUMSIG) THEN
              CP_SIG_ip1=CAS(IS+1,ID)*FP
              C_SIG(IS,ID)=-eFact*CP_SIG_ip1/DS_INCR(IS+1)
            END IF
          END DO
          ASPAR_DIAG(:,ID)=ASPAR_DIAG(:,ID) + B_SIG
        END DO
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef NCDF      
      SUBROUTINE DEBUG_EIMPS_TOTAL_JACOBI(iPass, iIter, FieldOut1)
      USE DATAPOOL
      USE NETCDF  
      IMPLICIT NONE
      INTEGER, intent(in) :: iPass, iIter
      REAL(rkind), intent(in) :: FieldOut1(MNP)
      character (len = *), parameter :: CallFct="DEBUG_EIMPS_TOTAL_JACOBI"
      REAL(rkind) :: FieldOutTotal1(np_total)
      REAL(rkind), allocatable :: ARRAY_loc(:)
      character(len=256) :: FileSave
      REAL(rkind) eTimeDay
      integer ncid, iret, nbTime, mnp_dims, ntime_dims, var_id
      integer fifteen_dims
      integer IP, IPloc, IPglob, NP_RESloc
      integer iProc
      integer, allocatable :: ListFirstMNP(:)
      WRITE(FileSave, 10) 'DebugJacobi', iPass
10    FORMAT(a, '_', i4.4,'.nc')
# ifdef MPI_PARALL_GRID
      IF (myrank .eq. 0) THEN
        allocate(ListFirstMNP(nproc), stat=istat)
        IF (istat/=0) CALL WWM_ABORT('wwm_jacobi, allocate error 1')
        ListFirstMNP=0
        DO iProc=2,nproc
          ListFirstMNP(iProc)=ListFirstMNP(iProc-1) + ListMNP(iProc-1)
        END DO
        DO IP=1,NP_RES
          IPglob=iplg(IP)
          FieldOutTotal1(IPglob)=FieldOut1(IP)
        END DO
        DO iPROC=2,nproc
          NP_RESloc=ListNP_RES(iPROC)
          allocate(ARRAY_loc(NP_RESloc), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_jacobi, allocate error 2')
          !
          CALL MPI_RECV(ARRAY_loc, NP_RESloc, rtype, iProc-1, 511, comm, istatus, ierr)
          DO IPloc=1,NP_RESloc
            IPglob=ListIPLG(IPloc + ListFirstMNP(iProc))
            FieldOutTotal1(IPglob)=ARRAY_loc(IPloc)
          END DO
          deallocate(ARRAY_loc)
        END DO
        deallocate(ListFirstMNP)
      ELSE
        CALL MPI_SEND(FieldOut1, NP_RES, rtype, 0, 511, comm, ierr)
      END IF
# else
      FieldOutTotal1 = FieldOut1
# endif
      !
      ! Now writing to netcdf file
      ! 
# ifdef MPI_PARALL_GRID
      IF (myrank .eq. 0) THEN
# endif
        IF (iIter .eq. 1) THEN
          iret = nf90_create(TRIM(FileSave), NF90_CLOBBER, ncid)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 1, iret)
          !
          iret = nf90_def_dim(ncid, 'fifteen', 15, fifteen_dims)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 2, iret)
          !
          nbTime=0
          CALL WRITE_NETCDF_TIME_HEADER(ncid, nbTime, ntime_dims)
          !
          iret = nf90_def_dim(ncid, 'mnp', np_total, mnp_dims)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 3, iret)
          !
          iret=nf90_def_var(ncid,"FieldOut1",NF90_RUNTYPE,(/ mnp_dims, ntime_dims/),var_id)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 4, iret)
          !
          iret = nf90_close(ncid)
          CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 5, iret)
        END IF
        !
        ! Writing data
        !
        iret = nf90_open(TRIM(FileSave), NF90_WRITE, ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 6, iret)
        !
        eTimeDay = MAIN%BMJD + MyREAL(iIter-1)*MyREAL(3600)/MyREAL(86400)
        CALL WRITE_NETCDF_TIME(ncid, iIter, eTimeDay)
        !
        iret=nf90_inq_varid(ncid, "FieldOut1", var_id)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 7, iret)
        !
        iret=nf90_put_var(ncid,var_id,FieldOutTotal1,start=(/1, iIter/), count=(/ np_total, 1 /))
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 8, iret)
        !
        iret = nf90_close(ncid)
        CALL GENERIC_NETCDF_ERROR_WWM(CallFct, 9, iret)
# ifdef MPI_PARALL_GRID
      END IF
# endif
      END SUBROUTINE DEBUG_EIMPS_TOTAL_JACOBI
#endif    
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE EIMPS_TOTAL_JACOBI_ITERATION
      USE DATAPOOL
      IMPLICIT NONE
      REAL(rkind) :: MaxNorm, SumNorm, p_is_converged
      REAL(rkind) :: eSum(NUMSIG,NUMDIR)
      REAL(rkind) :: BSIDE(NUMSIG,NUMDIR), DIAG(NUMSIG,NUMDIR)
      REAL(rkind) :: WALOC(NUMSIG,NUMDIR)
      REAL(rkind) :: CAD(NUMSIG,NUMDIR), CAS(NUMSIG,NUMDIR)
      REAL(rkind) :: BLOC(NUMSIG,NUMDIR)
      REAL(rkind) :: ASPAR_DIAG(NUMSIG,NUMDIR)
      REAL(rkind) :: SSLIM(NUMSIG,NUMDIR), SSBRL(NUMSIG,NUMDIR)
      LOGICAL     :: test=.true., LCONVERGED(MNP)
      REAL(rkind) :: ASPAR_LOC(NUMSIG,NUMDIR,MAX_DEG), MAXDAC(NUMSIG)
#ifdef DEBUG_ITERATION_LOOP
      integer iIter
      integer, save :: iPass = 0
      REAL(rkind) :: FieldOut1(MNP)
#endif
#ifdef TIMINGS
      REAL(rkind) :: TIME1, TIME2, TIME3, TIME4, TIME5
#endif
      REAL(rkind) :: eFact
      REAL(rkind) :: Sum_new, DiffNew
      INTEGER :: IP, J, idx, nbIter, is_converged(1), itmp(1)
      INTEGER :: JDX
      LOGICAL, SAVE :: InitCFLadvgeo = .FALSE.
      integer nbPassive
#ifdef DEBUG
      REAL(rkind) sumESUM
#endif
      WRITE(STAT%FHNDL,*) 'LCALC=', LCALC
      WRITE(STAT%FHNDL,*) 'SOURCE_IMPL=', SOURCE_IMPL
      WRITE(STAT%FHNDL,*) 'LNONL=', LNONL
      WRITE(STAT%FHNDL,*) 'REFRACTION_IMPL=', REFRACTION_IMPL
      WRITE(STAT%FHNDL,*) 'FREQ_SHIFT_IMPL=', FREQ_SHIFT_IMPL

      IF (WAE_JGS_CFL_LIM) THEN
        IF (InitCFLadvgeo .eqv. .FALSE.) THEN
          allocate(CFLadvgeoI(MNP), NumberOperationJGS(MNP), stat=istat)
          IF (istat/=0) CALL WWM_ABORT('wwm_jacobi, allocate error 3')
        END IF
        InitCFLadvgeo=.TRUE.
        NumberOperationJGS = 0
        IF (LCALC) THEN
          CALL COMPUTE_CFL_N_SCHEME_EXPLICIT(CFLadvgeoI)
        END IF
      END IF
#ifdef DEBUG
      CALL LOCAL_NODE_PRINT(20506, "Before Jacobi iteration")
#endif

#ifdef TIMINGS
      CALL WAV_MY_WTIME(TIME1)
#endif
      p_is_converged=0
      IF (ASPAR_LOCAL_LEVEL .le. 1) THEN
        CALL EIMPS_ASPAR_BLOCK(ASPAR_JAC)
      END IF
#ifdef DEBUG
      WRITE(STAT%FHNDL,*) 'sum(abs(ASPAR_JAC))=', sum(abs(ASPAR_JAC))
      WRITE(STAT%FHNDL,*) 'sum(    ASPAR_JAC )=', sum(ASPAR_JAC)
#endif
      IF ((ASPAR_LOCAL_LEVEL .ge. 5).and.(ASPAR_LOCAL_LEVEL .le. 7)) THEN
        CALL COMPUTE_K_CRFS_XYU
      END IF
#ifdef TIMINGS
      CALL WAV_MY_WTIME(TIME2)
#endif
      !
      IF (ASPAR_LOCAL_LEVEL .eq. 0) THEN
        CALL ADD_FREQ_DIR_TO_ASPAR_COMP_CADS(ASPAR_JAC)
#ifdef DEBUG
        WRITE(STAT%FHNDL,*) 'Aft/Refr/Freq sum(abs(ASPAR_JAC))=', sum(abs(ASPAR_JAC))
        WRITE(STAT%FHNDL,*) 'Aft/Refr/Freq sum(    ASPAR_JAC )=', sum(ASPAR_JAC)
#endif
      END IF

      IF (ASPAR_LOCAL_LEVEL .le. 1) THEN
        IF ((.NOT. LNONL) .AND. SOURCE_IMPL) THEN 
          DO IP = 1, NP_RES
            CALL GET_BSIDE_DIAG(IP, AC2, AC2, BSIDE, DIAG, BLOC)
            ASPAR_JAC(:,:,I_DIAG(IP)) = ASPAR_JAC(:,:,I_DIAG(IP)) + DIAG
            B_JAC(:,:,IP)             = BSIDE + BLOC
          ENDDO
        END IF
      END IF

#ifdef TIMINGS
      CALL WAV_MY_WTIME(TIME3)
#endif
      NumberIterationSolver = 0
      nbIter=0

      LCONVERGED = .FALSE. 

      DO
        is_converged(1) = 0
        JDX=0
#ifdef DEBUG_ITERATION_LOOP
        FieldOut1 = 0
#endif
#ifdef DEBUG
        WRITE(STAT%FHNDL,*) 'Before iteration sum(AC2)=', sum(abs(AC2))
        sumESUM=0
#endif
        nbPassive = 0
       
        DO IP=1,NP_RES

          IF (IOBDP(IP) .EQ. 0 .OR. LCONVERGED(IP)) THEN
            is_converged(1) = is_converged(1) + 1
            cycle
          END IF

          WALOC = AC2(:,:,IP)

          IF (WAE_JGS_CFL_LIM) THEN
            IF (NumberOperationJGS(IP) .lt. CFLadvgeoI(IP)) THEN
              test=.TRUE.
            ELSE
              test=.FALSE.
            END IF
          END IF

          IF (.true.) THEN
!            WRITE(STAT%FHNDL,*) 'IP=', IP
            NumberIterationSolver(IP) = NumberIterationSolver(IP) + 1
            CALL SINGLE_VERTEX_COMPUTATION(JDX, WALOC, eSum, ASPAR_DIAG)
#ifdef DEBUG
            sumESUM = sumESUM + sum(abs(eSum))
#endif
            eSum=eSum/ASPAR_DIAG

            IF (MELIM .EQ. 1) THEN
              CALL GET_MAXDAC(IP,MAXDAC)
              CALL ACTION_LIMITER_LOCAL(MAXDAC,AC1(:,:,IP),ESUM,SSLIM)
            ENDIF
            IF (LMAXETOT) CALL BREAKING_LIMITER_LOCAL(IP,ESUM,SSBRL)

            IF (BLOCK_GAUSS_SEIDEL) THEN
              AC2(:,:,IP)=eSum
            ELSE
              U_JACOBI(:,:,IP)=eSum
            END IF

            IF (JGS_CHKCONV) THEN
              Sum_new = sum(eSum)
              if (Sum_new .gt. thr8) then
                DiffNew=sum(abs(WALOC - eSum))
                p_is_converged = DiffNew/Sum_new
              else
                p_is_converged = zero
              endif
#ifdef DEBUG_ITERATION_LOOP
              FieldOut1(IP)=p_is_converged
#endif
!              WRITE(STAT%FHNDL,*) 'p_is_converged=', p_is_converged
              IF (IPstatus(IP) .eq. 1) THEN
                IF (p_is_converged .LT. THR) THEN ! not real ... mostly never touched point or whatorever ...
                  LCONVERGED(IP) = .FALSE. 
                ELSE IF (.NOT. p_is_converged .LT. THR .AND. p_is_converged .LT. jgs_diff_solverthr) THEN
                  !write(*,*) ip, p_is_converged, jgs_diff_solverthr
                  LCONVERGED(IP) = .TRUE.
                  is_converged(1) = is_converged(1) + 1
                  IF (WAE_JGS_CFL_LIM) THEN
                    NumberOperationJGS(IP) = NumberOperationJGS(IP) +1
                  END IF
                ENDIF
              END IF
            END IF!JGS_CHKCONV
          ELSE!test
            nbPassive = nbPassive + 1
            IF (JGS_CHKCONV .and. (IPstatus(IP) .eq. 1)) THEN
              is_converged(1) = is_converged(1) + 1
            END IF
          END IF!test
        END DO!IP

!        WRITE(*,*) SIZE(LCONVERGED), COUNT(LCONVERGED .eqv. .TRUE.)
#ifdef DEBUG
        WRITE(STAT%FHNDL,*) 'sumESUM=', sumESUM
#endif
!        WRITE(STAT%FHNDL,*) 'is_converged(1)=', is_converged(1)
!        WRITE(STAT%FHNDL,*) 'NP_RES=', NP_RES
!        WRITE(STAT%FHNDL,*) 'nbPassive=', nbPassive
!        WRITE(STAT%FHNDL,*) 'diffconv=', NP_RES - is_converged(1)
        IF (JGS_CHKCONV) THEN
#ifdef MPI_PARALL_GRID
          CALL MPI_ALLREDUCE(is_converged, itmp(1), 1, itype, MPI_SUM, COMM, ierr)
          is_converged = itmp
#endif
          p_is_converged = (real(np_total) - real(is_converged(1)))/real(np_total) * 100.
        ENDIF 

#ifdef MPI_PARALL_GRID
        IF (BLOCK_GAUSS_SEIDEL) THEN
          CALL EXCHANGE_P4D_WWM(AC2)
        ELSE
          CALL EXCHANGE_P4D_WWM(U_JACOBI)
        END IF
#endif
        IF (.NOT. BLOCK_GAUSS_SEIDEL) THEN
          AC2 = U_JACOBI
        ENDIF
#ifdef DEBUG
        WRITE(STAT%FHNDL,*) ' After iteration sum(AC2)=', sum(abs(AC2))
#endif
#ifdef DEBUG
        iIter=nbIter + 1
# ifdef NCDF        
        CALL DEBUG_EIMPS_TOTAL_JACOBI(iPass, iIter, FieldOut1)
# endif        
#endif
!
! The termination criterions several can be chosen
!
        WRITE(STAT%FHNDL,'(A10,4I10,E30.20,F10.5)') 'solver', nbiter, nbPassive, is_converged(1), np_total-is_converged(1), p_is_converged, pmin
        !
        ! Number of iterations. If too large the exit.
        !
        nbIter=nbIter+1
        IF (nbiter .eq. maxiter) THEN
          EXIT
        ENDIF
        !
        ! Check via number of converged points
        !
        IF (JGS_CHKCONV) THEN
          IF (p_is_converged .le. pmin) EXIT
        ENDIF
        !
        ! Check via the norm
        !
        IF (L_SOLVER_NORM) THEN
          CALL COMPUTE_JACOBI_SOLVER_ERROR(MaxNorm, SumNorm)
          IF (sqrt(SumNorm) .le. WAE_SOLVERTHR) THEN
            EXIT
          END IF
        END IF
      END DO
      WRITE(STAT%FHNDL,*) 'nbIter=', nbIter
#ifdef DEBUG
      CALL LOCAL_NODE_PRINT(20506, "After Jacobi Iteration")
#endif
#ifdef TIMINGS
      CALL WAV_MY_WTIME(TIME4)
#endif
      AC2 = MAX(ZERO, AC2) ! Make sure there is no negative energy left ... 

#ifdef TIMINGS
      CALL WAV_MY_WTIME(TIME5)
#endif

#ifdef TIMINGS
# ifdef MPI_PARALL_GRID
      IF (myrank == 0) THEN
# endif
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'PREPROCESSING SOURCES AND ADVECTION  ', TIME2-TIME1
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'PREPROCESSING REFRACTION             ', TIME3-TIME2
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'ITERATION                            ', TIME4-TIME3
        WRITE(STAT%FHNDL,'("+TRACE...",A,F15.6)') 'STORE RESULT                         ', TIME5-TIME4
        FLUSH(STAT%FHNDL)
# ifdef MPI_PARALL_GRID
      ENDIF
# endif
#endif
#ifdef DEBUG_ITERATION_LOOP
      iPass=iPass+1
#endif
      CONTAINS
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GET_BSIDE_DIAG(IP, ACin1, ACin2, BSIDE, DIAG, BLOC)
      IMPLICIT NONE
      INTEGER, intent(in) :: IP
      REAL(rkind), intent(in)  :: ACin1(NUMSIG,NUMDIR,MNP)
      REAL(rkind), intent(in)  :: ACin2(NUMSIG,NUMDIR,MNP)
      REAL(rkind), intent(out) :: BSIDE(NUMSIG,NUMDIR)
      REAL(rkind), intent(out) :: DIAG (NUMSIG,NUMDIR)
      REAL(rkind), intent(out) :: BLOC(NUMSIG,NUMDIR)
      REAL(rkind)              :: PHI(NUMSIG,NUMDIR)
      REAL(rkind)              :: DPHIDN(NUMSIG,NUMDIR)
      REAL(rkind)              :: eVal
      REAL(rkind)              :: ePHI, TheFactor, DVS1, DVS2
      REAL(rkind)              :: DAM(NUMSIG), MAXDAC, eDAM
      REAL(rkind)              :: DELFL(NUMSIG), eFric, USFM
      INTEGER IS, ID

      PHI=ZERO
      DPHIDN=ZERO

      IF (LNONL) THEN
         CALL SOURCES_IMPLICIT 
      ELSE
        DPHIDN = DPHIDNA(:,:,IP)
        PHI    = PHIA(:,:,IP)
      END IF

      eVal = SI(IP) * DT4A

      CALL GET_BLOCAL(IP, ACin1, BLOC)
!
      BSIDE =     eVal * (PHI - MIN(ZERO,DPHIDN) * Acin2(:,:,IP))
      DIAG  =   - eVal * MIN(ZERO,DPHIDN) ! AR: The minus put the DHPIDN on the left side of the equation as diagonal contributions with the right sign ... it inverts the sign ... however this is wrong now for IBREAK = 2 the SWAN stuff 
!
#ifdef DEBUG_SOURCE_TERM
      WRITE(STAT%FHNDL,'(I10,10G20.10,A40)') IP, SUM(ACin1), SUM(ACin2), SUM(PHI), SUM(DPHIDN), SUM(BSIDE), SUM(DIAG), SUM(BLOC), eval, 'GET_BSIDE_DIAG'
#endif

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE GET_BLOCAL(IP, Ac1in, BLOC)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IP
      REAL(rkind), intent(in)  :: AC1in(NUMSIG,NUMDIR,MNP)
      REAL(rkind), INTENT(OUT) :: BLOC(NUMSIG,NUMDIR)
      INTEGER ID, idx
      idx=IWBNDLC_REV(IP)
      IF ((LBCWA .OR. LBCSP).and.(idx.gt.0)) THEN
        BLOC = WBAC(:,:,idx)  * SI(IP)
      ELSE
        DO ID=1,NUMDIR
          BLOC(:,ID) = Ac1in(:,ID,IP) * IOBPD(ID,IP)*IOBWB(IP)*IOBDP(IP)*SI(IP)
        ENDDO
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE LINEAR_ASPAR_LOCAL(IP, ASPAR_LOC, ASPAR_DIAG, A_THE, C_THE, A_SIG, C_SIG)
      IMPLICIT NONE
      INTEGER, intent(in) :: IP
      REAL(rkind), intent(out) :: ASPAR_LOC(NUMSIG,NUMDIR,MAX_DEG)
      REAL(rkind), intent(out) :: ASPAR_DIAG(NUMSIG,NUMDIR)
      REAL(rkind), intent(out) :: A_THE(NUMSIG,NUMDIR), C_THE(NUMSIG,NUMDIR)
      REAL(rkind), intent(out) :: A_SIG(NUMSIG,NUMDIR), C_SIG(NUMSIG,NUMDIR)
      INTEGER :: POS_TRICK(3,2)
      REAL(rkind) :: FL11(NUMSIG,NUMDIR), FL12(NUMSIG,NUMDIR), FL21(NUMSIG,NUMDIR), FL22(NUMSIG,NUMDIR), FL31(NUMSIG,NUMDIR), FL32(NUMSIG,NUMDIR)
      REAL(rkind) :: CRFS(NUMSIG,NUMDIR,3), K1(NUMSIG,NUMDIR), KM(NUMSIG,NUMDIR,3), K(NUMSIG,NUMDIR,3), TRIA03
      REAL(rkind) :: CXY(2,NUMSIG,NUMDIR,3)
      REAL(rkind) :: DIFRU, USOC, WVC
      REAL(rkind) :: DELTAL(NUMSIG,NUMDIR,3)
      REAL(rkind) :: KP(NUMSIG,NUMDIR,3), NM(NUMSIG,NUMDIR)
      REAL(rkind) :: DTK(NUMSIG,NUMDIR), TMP3(NUMSIG,NUMDIR)
      REAL(rkind) :: LAMBDA(2,NUMSIG,NUMDIR)
      INTEGER     :: I1, I2, I3
      INTEGER     :: ID, IS, IE, IPOS
      INTEGER     :: I, ICON
      INTEGER     :: IP_fall, IPie, TheVal
      INTEGER     :: ID1, ID2, POS1, POS2
      REAL(rkind) :: CAD(NUMSIG,NUMDIR)
      REAL(rkind) :: CAS(NUMSIG,NUMDIR)
      REAL(rkind) :: CP_THE(NUMSIG,NUMDIR), CM_THE(NUMSIG,NUMDIR)
      REAL(rkind) :: CASS(0:NUMSIG+1), B_SIG(NUMSIG)
      REAL(rkind) :: CP_SIG(0:NUMSIG+1), CM_SIG(0:NUMSIG+1)
      REAL(rkind) :: eFact
      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2

      ASPAR_LOC=ZERO
      ASPAR_DIAG=ZERO
      DO ICON = 1, CCON(IP)
        IE     =  IE_CELL2(IP,ICON)
        IPOS   = POS_CELL2(IP,ICON)
        I1 = INE(1,IE)
        I2 = INE(2,IE)
        I3 = INE(3,IE)

        DO I=1,3
          IPie = INE(I,IE)
          DO ID=1,NUMDIR
            DO IS=1,NUMSIG
              IF (LSECU .OR. LSTCU) THEN
                CXY(1,IS,ID,I) = CG(IS,IPie)*COSTH(ID)+CURTXY(IPie,1)
                CXY(2,IS,ID,I) = CG(IS,IPie)*SINTH(ID)+CURTXY(IPie,2)
              ELSE
                CXY(1,IS,ID,I) = CG(IS,IPie)*COSTH(ID)
                CXY(2,IS,ID,I) = CG(IS,IPie)*SINTH(ID)
              END IF
              IF (LSPHE) THEN
                CXY(1,IS,ID,I) = CXY(1,IS,ID,I)*INVSPHTRANS(IPie,1)
                CXY(2,IS,ID,I) = CXY(2,IS,ID,I)*INVSPHTRANS(IPie,2)
              END IF
              IF (LDIFR) THEN
                CXY(1,IS,ID,I) = CXY(1,IS,ID,I)*DIFRM(IPie)
                CXY(2,IS,ID,I) = CXY(2,IS,ID,I)*DIFRM(IPie)
                IF (LSECU .OR. LSTCU) THEN
                  IF (IDIFFR .GT. 1) THEN
                    WVC = SPSIG(IS)/WK(IS,IPie)
                    USOC = (COSTH(ID)*CURTXY(IPie,1) + SINTH(ID)*CURTXY(IPie,2))/WVC
                    DIFRU = ONE + USOC * (ONE - DIFRM(IPie))
                  ELSE
                    DIFRU = DIFRM(IPie)
                  END IF
                  CXY(1,IS,ID,I) = CXY(1,IS,ID,I) + DIFRU*CURTXY(IPie,1)
                  CXY(2,IS,ID,I) = CXY(2,IS,ID,I) + DIFRU*CURTXY(IPie,2)
                END IF
              END IF
            END DO
          END DO
        END DO

        LAMBDA(:,:,:) = ONESIXTH * (CXY(:,:,:,1) + CXY(:,:,:,2) + CXY(:,:,:,3))
        K(:,:,1)  = LAMBDA(1,:,:) * IEN(1,IE) + LAMBDA(2,:,:) * IEN(2,IE)
        K(:,:,2)  = LAMBDA(1,:,:) * IEN(3,IE) + LAMBDA(2,:,:) * IEN(4,IE)
        K(:,:,3)  = LAMBDA(1,:,:) * IEN(5,IE) + LAMBDA(2,:,:) * IEN(6,IE)
        FL11(:,:) = CXY(1,:,:,2)*IEN(1,IE)+CXY(2,:,:,2)*IEN(2,IE)
        FL12(:,:) = CXY(1,:,:,3)*IEN(1,IE)+CXY(2,:,:,3)*IEN(2,IE)
        FL21(:,:) = CXY(1,:,:,3)*IEN(3,IE)+CXY(2,:,:,3)*IEN(4,IE)
        FL22(:,:) = CXY(1,:,:,1)*IEN(3,IE)+CXY(2,:,:,1)*IEN(4,IE)
        FL31(:,:) = CXY(1,:,:,1)*IEN(5,IE)+CXY(2,:,:,1)*IEN(6,IE)
        FL32(:,:) = CXY(1,:,:,2)*IEN(5,IE)+CXY(2,:,:,2)*IEN(6,IE)
        CRFS(:,:,1) = - ONESIXTH *  (TWO *FL31(:,:) + FL32(:,:) + FL21(:,:) + TWO * FL22(:,:) )
        CRFS(:,:,2) = - ONESIXTH *  (TWO *FL32(:,:) + TWO * FL11(:,:) + FL12(:,:) + FL31(:,:) )
        CRFS(:,:,3) = - ONESIXTH *  (TWO *FL12(:,:) + TWO * FL21(:,:) + FL22(:,:) + FL11(:,:) )
        KM = MIN(ZERO,K)
        KP(:,:,:) = MAX(ZERO,K)
        DELTAL(:,:,:) = CRFS(:,:,:) - KP(:,:,:)
        NM(:,:)=ONE/MIN(-THR,KM(:,:,1) + KM(:,:,2) + KM(:,:,3))
        TRIA03 = ONETHIRD * TRIA(IE)
        !
        IP_fall=INE(IPOS,IE)
        IF (IP_fall .ne. IP) THEN
          CALL WWM_ABORT('Bugs and many more bugs')
        END IF
        POS1=POS_IP_ADJ(1,IPOS,IE)
        POS2=POS_IP_ADJ(2,IPOS,IE)
!        I1=JA_IE(IPOS,1,IE)
!        I2=JA_IE(IPOS,2,IE)
!        I3=JA_IE(IPOS,3,IE)
        K1(:,:) =  KP(:,:,IPOS)
        DO ID=1,NUMDIR
          DTK(:,ID) =  K1(:,ID) * DT4A * IOBPD(ID,IP) * IOBWB(IP) * IOBDP(IP)
        END DO
        TMP3(:,:)  =  DTK(:,:) * NM(:,:)
        ASPAR_DIAG=ASPAR_DIAG + TRIA03+DTK(:,:)- TMP3(:,:) * DELTAL(:,:,IPOS)
        ASPAR_LOC(:,:,POS1)=ASPAR_LOC(:,:,POS1)-TMP3(:,:)*DELTAL(:,:,POS_TRICK(IPOS,1))
        ASPAR_LOC(:,:,POS2)=ASPAR_LOC(:,:,POS2)-TMP3(:,:)*DELTAL(:,:,POS_TRICK(IPOS,2))
      END DO
      IF (REFRACTION_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LTHBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPTHETA(IP,CAD)
        ELSE
          CAD=ZERO
        END IF
        CP_THE = MAX(ZERO,CAD)
        CM_THE = MIN(ZERO,CAD)
        eFact=(DT4D/DDIR)*SI(IP)
        DO ID=1,NUMDIR
          ID1 = ID_PREV(ID)
          ID2 = ID_NEXT(ID)
          A_THE(:,ID) = - eFact *  CP_THE(:,ID1)
          C_THE(:,ID) =   eFact *  CM_THE(:,ID2)
        END DO
        ASPAR_DIAG = ASPAR_DIAG + eFact * (CP_THE(:,:) - CM_THE(:,:))
      ELSE
        A_THE=ZERO
        C_THE=ZERO
      END IF
      IF (FREQ_SHIFT_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LSIGBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPSIGMA(IP,CAS)
        ELSE
          CAS=ZERO
        END IF
        eFact=DT4F*SI(IP)
        DO ID = 1, NUMDIR
          CASS(1:NUMSIG) = CAS(:,ID)
          CASS(0)     = 0.
          CASS(NUMSIG+1) = CASS(NUMSIG)
          CP_SIG = MAX(ZERO,CASS)
          CM_SIG = MIN(ZERO,CASS)
          DO IS=1,NUMSIG
            B_SIG(IS)=eFact*(CP_SIG(IS)/DS_INCR(IS-1) - CM_SIG(IS) /DS_INCR(IS))
          END DO
          DO IS=2,NUMSIG
            A_SIG(IS,ID) = - eFact*CP_SIG(IS-1)/DS_INCR(IS-1)
          END DO
          DO IS=1,NUMSIG-1
            C_SIG(IS,ID) = eFact*CM_SIG(IS+1)/DS_INCR(IS)
          END DO
          B_SIG(NUMSIG) = B_SIG(NUMSIG) + eFact*CM_SIG(NUMSIG+1)/DS_INCR(NUMSIG) * TAIL_ARR(5)
          ASPAR_DIAG(:,ID)=ASPAR_DIAG(:,ID) + B_SIG
        END DO
      ELSE
        A_SIG=ZERO
        C_SIG=ZERO
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SINGLE_VERTEX_COMPUTATION(JDX, WALOC, eSum, ASPAR_DIAG)
      IMPLICIT NONE
      integer, intent(inout) :: JDX
      REAL(rkind), intent(in) :: WALOC(NUMSIG,NUMDIR)
      REAL(rkind), intent(out) :: eSum(NUMSIG,NUMDIR), ASPAR_DIAG(NUMSIG,NUMDIR)
      integer ID1, ID2, ID, IS, IP_ADJ, IADJ
      REAL(rkind) :: NEG_P(NUMSIG,NUMDIR)
      REAL(rkind) :: CP_THE(NUMSIG,NUMDIR), CM_THE(NUMSIG,NUMDIR)
      REAL(rkind) :: CP_SIG(NUMSIG,NUMDIR), CM_SIG(NUMSIG,NUMDIR)
      REAL(rkind) :: A_THE(NUMSIG,NUMDIR), C_THE(NUMSIG,NUMDIR)
      REAL(rkind) :: A_SIG(NUMSIG,NUMDIR), C_SIG(NUMSIG,NUMDIR)

      IF (ASPAR_LOCAL_LEVEL .eq. 0) THEN
        ASPAR_DIAG=ASPAR_JAC(:,:,I_DIAG(IP))
        IF (SOURCE_IMPL) THEN
          IF (LNONL) THEN
            CALL GET_BSIDE_DIAG(IP, AC1, AC2, BSIDE, DIAG, BLOC)
            ASPAR_DIAG = ASPAR_DIAG + DIAG
            eSum = BLOC + BSIDE
          ELSE
            eSum = B_JAC(:,:,IP)
          END IF
        ELSE
          CALL GET_BLOCAL(IP, AC1, eSum)
        END IF
        DO J=IA(IP),IA(IP+1)-1 
          IF (J .ne. I_DIAG(IP)) eSum = eSum - ASPAR_JAC(:,:,J) * AC2(:,:,JA(J))
        END DO
        IF (REFRACTION_IMPL) THEN
          CAD=CAD_THE(:,:,IP)
          CP_THE = MAX(ZERO,CAD)
          CM_THE = MIN(ZERO,CAD)
          eFact=(DT4D/DDIR)*SI(IP)
          DO ID=1,NUMDIR
            ID1 = ID_PREV(ID)
            ID2 = ID_NEXT(ID)
            eSum(:,ID) = eSum(:,ID) + eFact*CP_THE(:,ID1)*WALOC(:,ID1)
            eSum(:,ID) = eSum(:,ID) - eFact*CM_THE(:,ID2)*WALOC(:,ID2)
          END DO
        END IF
        IF (FREQ_SHIFT_IMPL) THEN
          CAS=CAS_SIG(:,:,IP)
          CP_SIG = MAX(ZERO,CAS)
          CM_SIG = MIN(ZERO,CAS)
          eFact=DT4F*SI(IP)
          DO ID=1,NUMDIR
            DO IS=2,NUMSIG
              eSum(IS,ID)=eSum(IS,ID) + eFact*(CP_SIG(IS-1,ID)/DS_INCR(IS-1))*WALOC(IS-1,ID)
            END DO
            DO IS=1,NUMSIG-1
              eSum(IS,ID)=eSum(IS,ID) - eFact*(CM_SIG(IS+1,ID)/DS_INCR(IS))*WALOC(IS+1,ID)
            END DO
          END DO
        END IF
      ELSE IF (ASPAR_LOCAL_LEVEL .eq. 1) THEN
        ASPAR_DIAG=ASPAR_JAC(:,:,I_DIAG(IP))
        IF (SOURCE_IMPL) THEN
          IF (LNONL) THEN
            CALL GET_BSIDE_DIAG(IP, AC1, AC2, BSIDE, DIAG, BLOC)
            ASPAR_DIAG = ASPAR_DIAG + DIAG
            eSum = BLOC + BSIDE
          ELSE
            eSum = B_JAC(:,:,IP)
          END IF
        ELSE
          CALL GET_BLOCAL(IP, AC1, eSum)
        END IF
        CALL GET_FREQ_DIR_CONTRIBUTION(IP, ASPAR_DIAG, A_THE, C_THE, A_SIG, C_SIG)
        DO J=IA(IP),IA(IP+1)-1 
          IF (J .ne. I_DIAG(IP)) eSum = eSum - ASPAR_JAC(:,:,J) * AC2(:,:,JA(J))
        END DO
        IF (REFRACTION_IMPL) THEN
          DO ID=1,NUMDIR
            ID1 = ID_PREV(ID)
            ID2 = ID_NEXT(ID)
            eSum(:,ID) = eSum(:,ID) - A_THE(:,ID)*WALOC(:,ID1)
            eSum(:,ID) = eSum(:,ID) - C_THE(:,ID)*WALOC(:,ID2)
          END DO
        END IF
        IF (FREQ_SHIFT_IMPL) THEN
          DO ID=1,NUMDIR
            DO IS=2,NUMSIG
              eSum(IS,ID)=eSum(IS,ID) - A_SIG(IS,ID)*WALOC(IS-1,ID)
            END DO
            DO IS=1,NUMSIG-1
              eSum(IS,ID)=eSum(IS,ID) - C_SIG(IS,ID)*WALOC(IS+1,ID)
            END DO
          END DO
        END IF
      ELSE IF (ASPAR_LOCAL_LEVEL .eq. 2) THEN
        CALL LINEAR_ASPAR_LOCAL(IP, ASPAR_LOC, ASPAR_DIAG, A_THE, C_THE, A_SIG, C_SIG)
        IF (SOURCE_IMPL) THEN
          IF (LNONL) THEN
            CALL GET_BSIDE_DIAG(IP, AC1, AC2, BSIDE, DIAG, eSum)
          ELSE
            CALL GET_BSIDE_DIAG(IP, AC1, AC1, BSIDE, DIAG, eSum)
          END IF
          ASPAR_DIAG = ASPAR_DIAG + DIAG
          eSum = eSum + BSIDE
        END IF
        DO IADJ=1,VERT_DEG(IP)
          IP_ADJ=LIST_ADJ_VERT(IADJ,IP)
          eSum=eSum - ASPAR_LOC(:,:,IADJ)*AC2(:,:,IP_ADJ)
        END DO
        IF (REFRACTION_IMPL) THEN
          DO ID=1,NUMDIR
            ID1 = ID_PREV(ID)
            ID2 = ID_NEXT(ID)
            eSum(:,ID) = eSum(:,ID) - A_THE(:,ID)*WALOC(:,ID1)
            eSum(:,ID) = eSum(:,ID) - C_THE(:,ID)*WALOC(:,ID2)
          END DO
        END IF
        IF (FREQ_SHIFT_IMPL) THEN
          DO ID=1,NUMDIR
            DO IS=2,NUMSIG
              eSum(IS,ID)=eSum(IS,ID) - A_SIG(IS,ID)*WALOC(IS-1,ID)
            END DO
            DO IS=1,NUMSIG-1
              eSum(IS,ID)=eSum(IS,ID) - C_SIG(IS,ID)*WALOC(IS+1,ID)
            END DO
          END DO
        END IF
      ELSE IF (ASPAR_LOCAL_LEVEL .eq. 3) THEN
        CALL NEGATIVE_PART(IP, NEG_P, ASPAR_DIAG)
        IF (SOURCE_IMPL) THEN
          IF (LNONL) THEN
            CALL GET_BSIDE_DIAG(IP ,AC1, AC2, BSIDE, DIAG, eSum)
          ELSE
            CALL GET_BSIDE_DIAG(IP, AC1, AC1, BSIDE, DIAG, eSum)
          END IF
          ASPAR_DIAG = ASPAR_DIAG + DIAG
          eSum = eSum + BSIDE
        END IF
        eSum=eSum - NEG_P
      ELSE IF (ASPAR_LOCAL_LEVEL .eq. 4) THEN
        CALL NEGATIVE_PART_B(IP, NEG_P, ASPAR_DIAG)
        IF (SOURCE_IMPL) THEN
          IF (LNONL) THEN
            CALL GET_BSIDE_DIAG(IP, AC1, AC2, BSIDE, DIAG, eSum)
          ELSE
            CALL GET_BSIDE_DIAG(IP, AC1, AC1, BSIDE, DIAG, eSum)
          END IF
          ASPAR_DIAG = ASPAR_DIAG + DIAG
          eSum = eSum + BSIDE
        END IF
        eSum = eSum - NEG_P
      ELSE IF (ASPAR_LOCAL_LEVEL .eq. 5) THEN
        CALL NEGATIVE_PART_C(JDX, IP, NEG_P, ASPAR_DIAG)
        IF (SOURCE_IMPL) THEN
          IF (LNONL) THEN
            CALL GET_BSIDE_DIAG(IP,AC1, AC2, BSIDE, DIAG, eSum)
          ELSE
            CALL GET_BSIDE_DIAG(IP, AC1, AC1, BSIDE, DIAG, eSum)
          END IF
          ASPAR_DIAG = ASPAR_DIAG + DIAG
          eSum = eSum + BSIDE
        END IF
        eSum=eSum - NEG_P
      ELSE IF (ASPAR_LOCAL_LEVEL .eq. 6) THEN
        CALL NEGATIVE_PART_D(JDX, IP, NEG_P, ASPAR_DIAG)
        IF (SOURCE_IMPL) THEN
          IF (LNONL) THEN
            CALL GET_BSIDE_DIAG(IP, AC1, AC2, BSIDE, DIAG, eSum)
          ELSE
            CALL GET_BSIDE_DIAG(IP, AC1, AC1, BSIDE, DIAG, eSum)
          END IF
          ASPAR_DIAG = ASPAR_DIAG + DIAG
          eSum = eSum + BSIDE
        END IF
        eSum=eSum - NEG_P
      ELSE IF (ASPAR_LOCAL_LEVEL .eq. 7) THEN
        CALL NEGATIVE_PART_E(JDX, IP, NEG_P, ASPAR_DIAG)
        IF (SOURCE_IMPL) THEN
          IF (LNONL) THEN
            CALL GET_BSIDE_DIAG(IP, AC1, AC2, BSIDE, DIAG, eSum)
          ELSE
            CALL GET_BSIDE_DIAG(IP, AC1, AC1, BSIDE, DIAG, eSum)
          END IF
          ASPAR_DIAG = ASPAR_DIAG + DIAG
          eSum = eSum + BSIDE
        END IF
        eSum=eSum - NEG_P
      ELSE IF (ASPAR_LOCAL_LEVEL .eq. 8) THEN
        CALL NEGATIVE_PART_F(IP, NEG_P, ASPAR_DIAG)
        IF (SOURCE_IMPL) THEN
          IF (LNONL) THEN
            CALL GET_BSIDE_DIAG(IP, AC1, AC2, BSIDE, DIAG, eSum)
          ELSE
            CALL GET_BSIDE_DIAG(IP, AC1, AC1, BSIDE, DIAG, eSum)
          END IF
          ASPAR_DIAG = ASPAR_DIAG + DIAG
          eSum = eSum + BSIDE
        END IF
        eSum=eSum - NEG_P
      ELSE IF (ASPAR_LOCAL_LEVEL .eq. 9) THEN
        CALL NEGATIVE_PART_G(IP, NEG_P, ASPAR_DIAG)
        IF (SOURCE_IMPL) THEN
          IF (LNONL) THEN
            CALL GET_BSIDE_DIAG(IP, AC1, AC2, BSIDE, DIAG, eSum)
          ELSE
            CALL GET_BSIDE_DIAG(IP, AC1, AC1, BSIDE, DIAG, eSum)
          END IF
          ASPAR_DIAG = ASPAR_DIAG + DIAG
          eSum = eSum + BSIDE
        END IF
        eSum=eSum - NEG_P
      ELSE IF (ASPAR_LOCAL_LEVEL .eq. 10) THEN
        CALL NEGATIVE_PART_H(IP, NEG_P, ASPAR_DIAG)
        IF (SOURCE_IMPL) THEN
          IF (LNONL) THEN
            CALL GET_BSIDE_DIAG(IP, AC1, AC2, BSIDE, DIAG, eSum)
          ELSE
            CALL GET_BSIDE_DIAG(IP, AC1, AC1, BSIDE, DIAG, eSum)
          END IF
          ASPAR_DIAG = ASPAR_DIAG + DIAG
          eSum = eSum + BSIDE
        END IF
        eSum=eSum - NEG_P
      ELSE
        CALL WWM_ABORT('Not defined')
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_JACOBI_SOLVER_ERROR(MaxNorm, SumNorm)
      IMPLICIT NONE
      real(rkind), intent(out) :: MaxNorm, SumNorm
      integer IP
      integer ID1, ID2, ID, IS, IP_ADJ, IADJ
#ifdef MPI_PARALL_GRID
      REAL(rkind) :: Norm_L2_gl(NUMSIG,NUMDIR), Norm_LINF_gl(NUMSIG,NUMDIR)
#endif
      REAL(rkind) :: Norm_L2(NUMSIG,NUMDIR), Norm_LINF(NUMSIG,NUMDIR)
      REAL(rkind) :: ASPAR_DIAG(NUMSIG,NUMDIR), WALOC(NUMSIG,NUMDIR)
      REAL(rkind) :: CP_THE(NUMSIG,NUMDIR), CM_THE(NUMSIG,NUMDIR)
      REAL(rkind) :: CP_SIG(NUMSIG,NUMDIR), CM_SIG(NUMSIG,NUMDIR)
      REAL(rkind) :: A_THE(NUMSIG,NUMDIR), C_THE(NUMSIG,NUMDIR)
      REAL(rkind) :: A_SIG(NUMSIG,NUMDIR), C_SIG(NUMSIG,NUMDIR)
      REAL(rkind) :: BLOC(NUMSIG,NUMDIR)
      REAL(rkind) :: NEG_P(NUMSIG,NUMDIR)
      Norm_L2=0
      Norm_LINF=0
      DO IP=1,NP_RES
        WALOC = AC2(:,:,IP)
        IF (ASPAR_LOCAL_LEVEL .eq. 0) THEN
          ASPAR_DIAG=ASPAR_JAC(:,:,I_DIAG(IP))
          IF (SOURCE_IMPL) THEN
            IF (LNONL) THEN
              CALL GET_BSIDE_DIAG(IP, AC1, AC2, BSIDE, DIAG, BLOC)
              ASPAR_DIAG = ASPAR_DIAG + DIAG
              eSum = BLOC + BSIDE
            ELSE
              eSum = B_JAC(:,:,IP)
            END IF
          ELSE
            CALL GET_BLOCAL(IP, AC1, eSum)
          END IF
          DO J=IA(IP),IA(IP+1)-1
            idx=JA(J)
            IF (J .eq. I_DIAG(IP)) THEN
              eSum=eSum - ASPAR_DIAG*AC2(:,:,idx)
            ELSE
              eSum=eSum - ASPAR_JAC(:,:,J)*AC2(:,:,idx)
            END IF
          END DO
          IF (REFRACTION_IMPL) THEN
            CAD=CAD_THE(:,:,IP)
            CP_THE = MAX(ZERO,CAD)
            CM_THE = MIN(ZERO,CAD)
            eFact=(DT4D/DDIR)*SI(IP)
            DO ID=1,NUMDIR
              ID1 = ID_PREV(ID)
              ID2 = ID_NEXT(ID)
              eSum(:,ID) = eSum(:,ID) + eFact*CP_THE(:,ID1)*AC2(:,ID1,IP)
              eSum(:,ID) = eSum(:,ID) - eFact*CM_THE(:,ID2)*AC2(:,ID2,IP)
            END DO
          END IF
          IF (FREQ_SHIFT_IMPL) THEN
            CAS=CAS_SIG(:,:,IP)
            CP_SIG = MAX(ZERO,CAS)
            CM_SIG = MIN(ZERO,CAS)
            eFact=DT4F*SI(IP)
            DO ID=1,NUMDIR
              DO IS=2,NUMSIG
                eSum(IS,ID)=eSum(IS,ID) + eFact*(CP_SIG(IS-1,ID)/DS_INCR(IS-1))*AC2(IS-1,ID,IP)
              END DO
              DO IS=1,NUMSIG-1
                eSum(IS,ID)=eSum(IS,ID) - eFact*(CM_SIG(IS+1,ID)/DS_INCR(IS))*AC2(IS+1,ID,IP)
              END DO
            END DO
          END IF
        ELSE IF (ASPAR_LOCAL_LEVEL .eq. 1) THEN
          ASPAR_DIAG=ASPAR_JAC(:,:,I_DIAG(IP))
          IF (SOURCE_IMPL) THEN
            IF (LNONL) THEN
              CALL GET_BSIDE_DIAG(IP, AC1, AC2, BSIDE, DIAG, BLOC)
              ASPAR_DIAG = ASPAR_DIAG + DIAG
              eSum = BLOC + BSIDE
            ELSE
              eSum = B_JAC(:,:,IP)
            END IF
          ELSE
            CALL GET_BLOCAL(IP, AC1, eSum)
          END IF
          CALL GET_FREQ_DIR_CONTRIBUTION(IP, ASPAR_DIAG, A_THE, C_THE, A_SIG, C_SIG)
          DO J=IA(IP),IA(IP+1)-1
            idx=JA(J)
            IF (J .eq. I_DIAG(IP)) THEN
              eSum=eSum - ASPAR_DIAG*AC2(:,:,idx)
            ELSE
              eSum=eSum - ASPAR_JAC(:,:,J)*AC2(:,:,idx)
            END IF
          END DO
          IF (REFRACTION_IMPL) THEN
            DO ID=1,NUMDIR
              ID1 = ID_PREV(ID)
              ID2 = ID_NEXT(ID)
              eSum(:,ID) = eSum(:,ID) - A_THE(:,ID)*WALOC(:,ID1)
              eSum(:,ID) = eSum(:,ID) - C_THE(:,ID)*WALOC(:,ID2)
            END DO
          END IF
          IF (FREQ_SHIFT_IMPL) THEN
            DO ID=1,NUMDIR
              DO IS=2,NUMSIG
                eSum(IS,ID)=eSum(IS,ID) - A_SIG(IS,ID)*WALOC(IS-1,ID)
              END DO
              DO IS=1,NUMSIG-1
                eSum(IS,ID)=eSum(IS,ID) - C_SIG(IS,ID)*WALOC(IS+1,ID)
              END DO
            END DO
          END IF
        ELSE IF (ASPAR_LOCAL_LEVEL .eq. 2) THEN
          CALL LINEAR_ASPAR_LOCAL(IP, ASPAR_LOC, ASPAR_DIAG, A_THE, C_THE, A_SIG, C_SIG)
          IF (SOURCE_IMPL) THEN
            IF (LNONL) THEN
              CALL GET_BSIDE_DIAG(IP, AC1, AC2, BSIDE, DIAG, eSum)
            ELSE
              CALL GET_BSIDE_DIAG(IP, AC1, AC1, BSIDE, DIAG, eSum)
            END IF
            ASPAR_DIAG = ASPAR_DIAG + DIAG
            eSum = eSum + BSIDE
          END IF
          DO IADJ=1,VERT_DEG(IP)
            IP_ADJ=LIST_ADJ_VERT(IADJ,IP)
            eSum=eSum - ASPAR_LOC(:,:,IADJ)*AC2(:,:,IP_ADJ)
          END DO
          eSum=eSum - ASPAR_DIAG*AC2(:,:,IP)
          IF (REFRACTION_IMPL) THEN
            DO ID=1,NUMDIR
              ID1 = ID_PREV(ID)
              ID2 = ID_NEXT(ID)
              eSum(:,ID) = eSum(:,ID) - A_THE(:,ID)*WALOC(:,ID1)
              eSum(:,ID) = eSum(:,ID) - C_THE(:,ID)*WALOC(:,ID2)
            END DO
          END IF
          IF (FREQ_SHIFT_IMPL) THEN
            DO ID=1,NUMDIR
              DO IS=2,NUMSIG
                eSum(IS,ID)=eSum(IS,ID) - A_SIG(IS,ID)*WALOC(IS-1,ID)
              END DO
              DO IS=1,NUMSIG-1
                eSum(IS,ID)=eSum(IS,ID) - C_SIG(IS,ID)*WALOC(IS+1,ID)
              END DO
            END DO
          END IF
        ELSE IF (ASPAR_LOCAL_LEVEL .eq. 3) THEN
          CALL NEGATIVE_PART(IP, NEG_P, ASPAR_DIAG)
          IF (SOURCE_IMPL) THEN
            IF (LNONL) THEN
              CALL GET_BSIDE_DIAG(IP, AC1, AC2, BSIDE, DIAG, eSum)
            ELSE
              CALL GET_BSIDE_DIAG(IP, AC1, AC2, BSIDE, DIAG, eSum)
            END IF
            ASPAR_DIAG = ASPAR_DIAG + DIAG
            eSum = eSum + BSIDE
          END IF
          eSum = eSum - NEG_P - ASPAR_DIAG*AC2(:,:,IP)
        ELSE IF (ASPAR_LOCAL_LEVEL .eq. 4) THEN
          CALL NEGATIVE_PART_B(IP, NEG_P, ASPAR_DIAG)
          IF (SOURCE_IMPL) THEN
            IF (LNONL) THEN
              CALL GET_BSIDE_DIAG(IP, AC1, AC2, BSIDE, DIAG, eSum)
            ELSE
              CALL GET_BSIDE_DIAG(IP, AC1, AC2, BSIDE, DIAG, eSum)
            END IF
            ASPAR_DIAG = ASPAR_DIAG + DIAG
            eSum = eSum + BSIDE
          END IF
          eSum = eSum - NEG_P - ASPAR_DIAG*AC2(:,:,IP)
        ELSE
          CALL WWM_ABORT('Wrong selection')
        END IF
        IF (IPstatus(IP) .eq. 1) THEN
          Norm_L2 = Norm_L2 + (eSum**2)
        END IF
        Norm_LINF = max(Norm_LINF, abs(eSum))
      END DO
#ifdef MPI_PARALL_GRID
      CALL MPI_ALLREDUCE(Norm_LINF, Norm_LINF_gl, NUMSIG*NUMDIR,rtype,MPI_MAX,comm,ierr)
      CALL MPI_ALLREDUCE(Norm_L2, Norm_L2_gl, NUMSIG*NUMDIR, rtype,MPI_SUM,comm,ierr)
      MaxNorm = maxval(Norm_L2_gl)
      SumNorm = sum(Norm_L2_gl)
#else
      MaxNorm = maxval(Norm_L2)
      SumNorm = sum(Norm_L2)
#endif
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE NEGATIVE_PART_B(IP, NEG_P, ASPAR_DIAG)
      IMPLICIT NONE
      INTEGER, intent(in) :: IP
      REAL(rkind), intent(out) :: NEG_P(NUMSIG,NUMDIR)
      REAL(rkind), intent(out) :: ASPAR_DIAG(NUMSIG,NUMDIR)
      INTEGER :: POS_TRICK(3,2)
      REAL(rkind) :: FL11(NUMSIG,NUMDIR), FL12(NUMSIG,NUMDIR), FL21(NUMSIG,NUMDIR), FL22(NUMSIG,NUMDIR), FL31(NUMSIG,NUMDIR), FL32(NUMSIG,NUMDIR)
      REAL(rkind) :: CRFS(NUMSIG,NUMDIR,3), K1(NUMSIG,NUMDIR), KM(NUMSIG,NUMDIR,3), K(NUMSIG,NUMDIR,3), TRIA03
      REAL(rkind) :: CXY(2,NUMSIG,NUMDIR,3)
      REAL(rkind) :: DIFRU, USOC, WVC
      REAL(rkind) :: DELTAL(NUMSIG,NUMDIR,3)
      REAL(rkind) :: KP(NUMSIG,NUMDIR,3), NM(NUMSIG,NUMDIR)
      REAL(rkind) :: DTK(NUMSIG,NUMDIR), TMP3(NUMSIG,NUMDIR)
      REAL(rkind) :: LAMBDA(2,NUMSIG,NUMDIR)
      REAL(rkind) :: eF(NUMSIG,NUMDIR)
      INTEGER     :: I1, I2, I3
      INTEGER     :: ID, IS, IE, IPOS
      INTEGER     :: I, ICON
      INTEGER     :: IPie, TheVal
      INTEGER     :: ID1, ID2, IP1, IP2
      REAL(rkind) :: CAD(NUMSIG,NUMDIR)
      REAL(rkind) :: CAS(NUMSIG,NUMDIR)
      REAL(rkind) :: CP_THE(NUMSIG,NUMDIR), CM_THE(NUMSIG,NUMDIR)
      REAL(rkind) :: CASS(0:NUMSIG+1), B_SIG(NUMSIG)
      REAL(rkind) :: CP_SIG(0:NUMSIG+1), CM_SIG(0:NUMSIG+1)
      REAL(rkind) :: eFact
      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2

      NEG_P=ZERO
      ASPAR_DIAG=ZERO
      DO ICON = 1, CCON(IP)
        IE     =  IE_CELL2(IP,ICON)
        IPOS   = POS_CELL2(IP,ICON)
        I1 = INE(1,IE)
        I2 = INE(2,IE)
        I3 = INE(3,IE)
        DO I=1,3
          IPie = INE(I,IE)
          DO ID=1,NUMDIR
            DO IS=1,NUMSIG
              IF (LSECU .OR. LSTCU) THEN
                CXY(1,IS,ID,I) = CG(IS,IPie)*COSTH(ID)+CURTXY(IPie,1)
                CXY(2,IS,ID,I) = CG(IS,IPie)*SINTH(ID)+CURTXY(IPie,2)
              ELSE
                CXY(1,IS,ID,I) = CG(IS,IPie)*COSTH(ID)
                CXY(2,IS,ID,I) = CG(IS,IPie)*SINTH(ID)
              END IF
              IF (LSPHE) THEN
                CXY(1,IS,ID,I) = CXY(1,IS,ID,I)*INVSPHTRANS(IPie,1)
                CXY(2,IS,ID,I) = CXY(2,IS,ID,I)*INVSPHTRANS(IPie,2)
              END IF
              IF (LDIFR) THEN
                CXY(1,IS,ID,I) = CXY(1,IS,ID,I)*DIFRM(IPie)
                CXY(2,IS,ID,I) = CXY(2,IS,ID,I)*DIFRM(IPie)
                IF (LSECU .OR. LSTCU) THEN
                  IF (IDIFFR .GT. 1) THEN
                    WVC = SPSIG(IS)/WK(IS,IPie)
                    USOC = (COSTH(ID)*CURTXY(IPie,1) + SINTH(ID)*CURTXY(IPie,2))/WVC
                    DIFRU = ONE + USOC * (ONE - DIFRM(IPie))
                  ELSE
                    DIFRU = DIFRM(IPie)
                  END IF
                  CXY(1,IS,ID,I) = CXY(1,IS,ID,I) + DIFRU*CURTXY(IPie,1)
                  CXY(2,IS,ID,I) = CXY(2,IS,ID,I) + DIFRU*CURTXY(IPie,2)
                END IF
              END IF
            END DO
          END DO
        END DO
        LAMBDA(:,:,:) = ONESIXTH * (CXY(:,:,:,1) + CXY(:,:,:,2) + CXY(:,:,:,3))
        K(:,:,1)  = LAMBDA(1,:,:) * IEN(1,IE) + LAMBDA(2,:,:) * IEN(2,IE)
        K(:,:,2)  = LAMBDA(1,:,:) * IEN(3,IE) + LAMBDA(2,:,:) * IEN(4,IE)
        K(:,:,3)  = LAMBDA(1,:,:) * IEN(5,IE) + LAMBDA(2,:,:) * IEN(6,IE)
        FL11(:,:) = CXY(1,:,:,2)*IEN(1,IE)+CXY(2,:,:,2)*IEN(2,IE)
        FL12(:,:) = CXY(1,:,:,3)*IEN(1,IE)+CXY(2,:,:,3)*IEN(2,IE)
        FL21(:,:) = CXY(1,:,:,3)*IEN(3,IE)+CXY(2,:,:,3)*IEN(4,IE)
        FL22(:,:) = CXY(1,:,:,1)*IEN(3,IE)+CXY(2,:,:,1)*IEN(4,IE)
        FL31(:,:) = CXY(1,:,:,1)*IEN(5,IE)+CXY(2,:,:,1)*IEN(6,IE)
        FL32(:,:) = CXY(1,:,:,2)*IEN(5,IE)+CXY(2,:,:,2)*IEN(6,IE)
        CRFS(:,:,1) = - ONESIXTH *  (TWO *FL31(:,:) + FL32(:,:) + FL21(:,:) + TWO * FL22(:,:) )
        CRFS(:,:,2) = - ONESIXTH *  (TWO *FL32(:,:) + TWO * FL11(:,:) + FL12(:,:) + FL31(:,:) )
        CRFS(:,:,3) = - ONESIXTH *  (TWO *FL12(:,:) + TWO * FL21(:,:) + FL22(:,:) + FL11(:,:) )
        KM = MIN(ZERO,K)
        KP(:,:,:) = MAX(ZERO,K)
        DELTAL(:,:,:) = CRFS(:,:,:) - KP(:,:,:)
        NM(:,:)=ONE/MIN(-THR,KM(:,:,1) + KM(:,:,2) + KM(:,:,3))
        TRIA03 = ONETHIRD * TRIA(IE)
        !
        IP1=INE(POS_TRICK(IPOS,1),IE)
        IP2=INE(POS_TRICK(IPOS,2),IE)
        K1(:,:) =  KP(:,:,IPOS)
        DO ID=1,NUMDIR
          DTK(:,ID) =  K1(:,ID) * DT4A * IOBPD(ID,IP) * IOBWB(IP) * IOBDP(IP)
        END DO
        TMP3(:,:)  =  DTK(:,:) * NM(:,:)
        ASPAR_DIAG=ASPAR_DIAG + TRIA03+DTK(:,:)- TMP3(:,:) * DELTAL(:,:,IPOS)
        eF(:,:) = -TMP3(:,:)*DELTAL(:,:,POS_TRICK(IPOS,1))
        NEG_P=NEG_P  + eF(:,:)*AC2(:,:,IP1)
        eF(:,:) = -TMP3(:,:)*DELTAL(:,:,POS_TRICK(IPOS,2))
        NEG_P=NEG_P  + eF(:,:)*AC2(:,:,IP2)
      END DO
      IF (REFRACTION_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LTHBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPTHETA(IP,CAD)
          CP_THE = MAX(ZERO,CAD)
          CM_THE = MIN(ZERO,CAD)
          eFact=(DT4D/DDIR)*SI(IP)
          DO ID=1,NUMDIR
            ID1 = ID-1
            ID2 = ID+1
            IF (ID == 1) ID1 = NUMDIR
            IF (ID == NUMDIR) ID2 = 1
            NEG_P(:,ID)=NEG_P(:,ID) - eFact*CP_THE(:,ID1)*AC2(:,ID1,IP)
            NEG_P(:,ID)=NEG_P(:,ID) + eFact*CM_THE(:,ID2)*AC2(:,ID2,IP)
          END DO
          ASPAR_DIAG = ASPAR_DIAG + eFact * (CP_THE(:,:) - CM_THE(:,:))
        END IF
      END IF
      IF (FREQ_SHIFT_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LSIGBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPSIGMA(IP,CAS)
          eFact=DT4F*SI(IP)
          DO ID = 1, NUMDIR
            CASS(1:NUMSIG) = CAS(:,ID)
            CASS(0)     = 0.
            CASS(NUMSIG+1) = CASS(NUMSIG)
            CP_SIG = MAX(ZERO,CASS)
            CM_SIG = MIN(ZERO,CASS)
            DO IS=1,NUMSIG
              B_SIG(IS)=eFact*(CP_SIG(IS)/DS_INCR(IS-1) - CM_SIG(IS) /DS_INCR(IS))
            END DO
            DO IS=2,NUMSIG
              NEG_P(IS,ID)=NEG_P(IS,ID) - eFact*CP_SIG(IS-1)/DS_INCR(IS-1)*AC2(IS-1,ID,IP)
            END DO
            DO IS=1,NUMSIG-1
              NEG_P(IS,ID)=NEG_P(IS,ID) + eFact*CM_SIG(IS+1)/DS_INCR(IS)*AC2(IS+1,ID,IP)
            END DO
            B_SIG(NUMSIG) = B_SIG(NUMSIG) + eFact*CM_SIG(NUMSIG+1)/DS_INCR(NUMSIG) * TAIL_ARR(5)
            ASPAR_DIAG(:,ID)=ASPAR_DIAG(:,ID) + B_SIG
          END DO
        END IF
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE NEGATIVE_PART(IP, NEG_P, ASPAR_DIAG)
      IMPLICIT NONE
      INTEGER, intent(in) :: IP
      REAL(rkind), intent(out) :: NEG_P(NUMSIG,NUMDIR)
      REAL(rkind), intent(out) :: ASPAR_DIAG(NUMSIG,NUMDIR)
      REAL(rkind) :: FL11_X, FL12_X, FL21_X, FL22_X, FL31_X, FL32_X
      REAL(rkind) :: FL11_Y, FL12_Y, FL21_Y, FL22_Y, FL31_Y, FL32_Y
      REAL(rkind) :: FL11_U, FL12_U, FL21_U, FL22_U, FL31_U, FL32_U
      REAL(rkind) :: CRFS(3), KM(3), K(3), TRIA03
      REAL(rkind) :: DELTAL(3)
      REAL(rkind) :: KP(3), NM, val1, val2
      REAL(rkind) :: K_X(3), K_Y(3), CRFS_X(3), CRFS_Y(3)
      REAL(rkind) :: CSX(3), CSY(3)
      REAL(rkind) :: CRFS_U(3), K_U(3), LAMBDA_UX, LAMBDA_UY
      REAL(rkind) :: UV_CUR(3,2)
      REAL(rkind) :: DTK, TMP3
      REAL(rkind) :: LAMBDA_X, LAMBDA_Y
      INTEGER     :: I1, I2, I3
      INTEGER     :: ID, IS, IE, IPOS
      INTEGER     :: I, ICON
      INTEGER     :: IPie, TheVal, IP1, IP2
      INTEGER     :: ID1, ID2
      REAL(rkind) :: CAD(NUMSIG,NUMDIR)
      REAL(rkind) :: CAS(NUMSIG,NUMDIR)
      REAL(rkind) :: CP_THE(NUMSIG,NUMDIR), CM_THE(NUMSIG,NUMDIR)
      REAL(rkind) :: CASS(0:NUMSIG+1), B_SIG(NUMSIG)
      REAL(rkind) :: CP_SIG(0:NUMSIG+1), CM_SIG(0:NUMSIG+1)
      REAL(rkind) :: eFact
      INTEGER :: POS_TRICK(3,2)
      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2

      NEG_P=ZERO
      ASPAR_DIAG=ZERO
      DO ICON = 1, CCON(IP)
        IE     =  IE_CELL2(IP,ICON)
        IPOS   = POS_CELL2(IP,ICON)
        I1 = INE(1,IE)
        I2 = INE(2,IE)
        I3 = INE(3,IE)
        IF (LSECU .OR. LSTCU) THEN
          IF (LSPHE) THEN
            DO I=1,3
              IPie=INE(I,IE)
              UV_CUR(I,:)=CURTXY(IPie,:)*INVSPHTRANS(IPie,:)
            END DO
          ELSE
            DO I=1,3
              IPie=INE(I,IE)
              UV_CUR(I,:)=CURTXY(IPie,:)
            END DO
          END IF
          LAMBDA_UX=ONESIXTH*(UV_CUR(1,1)+UV_CUR(2,1)+UV_CUR(3,1))
          LAMBDA_UY=ONESIXTH*(UV_CUR(1,2)+UV_CUR(2,2)+UV_CUR(3,2))
          K_U(1)  = LAMBDA_UX * IEN(1,IE) + LAMBDA_UY * IEN(2,IE)
          K_U(2)  = LAMBDA_UX * IEN(3,IE) + LAMBDA_UY * IEN(4,IE)
          K_U(3)  = LAMBDA_UX * IEN(5,IE) + LAMBDA_UY * IEN(6,IE)
          FL11_U = UV_CUR(2,1)*IEN(1,IE)+UV_CUR(2,2)*IEN(2,IE)
          FL12_U = UV_CUR(3,1)*IEN(1,IE)+UV_CUR(3,2)*IEN(2,IE)
          FL21_U = UV_CUR(3,1)*IEN(3,IE)+UV_CUR(3,2)*IEN(4,IE)
          FL22_U = UV_CUR(1,1)*IEN(3,IE)+UV_CUR(1,2)*IEN(4,IE)
          FL31_U = UV_CUR(1,1)*IEN(5,IE)+UV_CUR(1,2)*IEN(6,IE)
          FL32_U = UV_CUR(2,1)*IEN(5,IE)+UV_CUR(2,2)*IEN(6,IE)
          CRFS_U(1) = - ONESIXTH*(TWO *FL31_U + FL32_U + FL21_U + TWO * FL22_U)
          CRFS_U(2) = - ONESIXTH*(TWO *FL32_U + TWO * FL11_U + FL12_U + FL31_U)
          CRFS_U(3) = - ONESIXTH*(TWO *FL12_U + TWO * FL21_U + FL22_U + FL11_U)
        ELSE
          K_U=ZERO
          CRFS_U=ZERO
        END IF
        IP1=INE(POS_TRICK(IPOS,1),IE)
        IP2=INE(POS_TRICK(IPOS,2),IE)
        DO IS=1,NUMSIG
          IF (LSPHE) THEN
            DO I=1,3
              IPie=INE(I,IE)
              CSX(I)=CG(IS,IPie)*INVSPHTRANS(IPie,1)
              CSY(I)=CG(IS,IPie)*INVSPHTRANS(IPie,2)
            END DO
          ELSE
            DO I=1,3
              IPie=INE(I,IE)
              CSX(I)=CG(IS,IPie)
              CSY(I)=CG(IS,IPie)
            END DO
          END IF
          LAMBDA_X=ONESIXTH * (CSX(1) + CSX(2) + CSX(3))
          LAMBDA_Y=ONESIXTH * (CSY(1) + CSY(2) + CSY(3))
          K_X(1)  = LAMBDA_X * IEN(1,IE)
          K_X(2)  = LAMBDA_X * IEN(3,IE)
          K_X(3)  = LAMBDA_X * IEN(5,IE)
          K_Y(1)  = LAMBDA_Y * IEN(2,IE)
          K_Y(2)  = LAMBDA_Y * IEN(4,IE)
          K_Y(3)  = LAMBDA_Y * IEN(6,IE)

          FL11_X = CSX(2)*IEN(1,IE)
          FL12_X = CSX(3)*IEN(1,IE)
          FL21_X = CSX(3)*IEN(3,IE)
          FL22_X = CSX(1)*IEN(3,IE)
          FL31_X = CSX(1)*IEN(5,IE)
          FL32_X = CSX(2)*IEN(5,IE)
          FL11_Y = CSY(2)*IEN(2,IE)
          FL12_Y = CSY(3)*IEN(2,IE)
          FL21_Y = CSY(3)*IEN(4,IE)
          FL22_Y = CSY(1)*IEN(4,IE)
          FL31_Y = CSY(1)*IEN(6,IE)
          FL32_Y = CSY(2)*IEN(6,IE)

          CRFS_X(1)= - ONESIXTH*(TWO *FL31_X + FL32_X + FL21_X + TWO * FL22_X )
          CRFS_X(2)= - ONESIXTH*(TWO *FL32_X + TWO * FL11_X + FL12_X + FL31_X )
          CRFS_X(3)= - ONESIXTH*(TWO *FL12_X + TWO * FL21_X + FL22_X + FL11_X )
          CRFS_Y(1)= - ONESIXTH*(TWO *FL31_Y + FL32_Y + FL21_Y + TWO * FL22_Y )
          CRFS_Y(2)= - ONESIXTH*(TWO *FL32_Y + TWO * FL11_Y + FL12_Y + FL31_Y )
          CRFS_Y(3)= - ONESIXTH*(TWO *FL12_Y + TWO * FL21_Y + FL22_Y + FL11_Y )

          DO ID=1,NUMDIR
            DO I=1,3
              K(I)=K_X(I)*COSTH(ID) + K_Y(I)*SINTH(ID) + K_U(I)
              CRFS(I)=CRFS_X(I)*COSTH(ID) + CRFS_Y(I)*SINTH(ID) + CRFS_U(I)
            END DO
            KM = MIN(ZERO,K)
            KP = MAX(ZERO,K)
            DELTAL = CRFS - KP
            NM=ONE/MIN(-THR,SUM(KM))
            TRIA03 = ONETHIRD * TRIA(IE)
            !
            DTK =  KP(IPOS) * DT4A * IOBPD(ID,IP) * IOBWB(IP) * IOBDP(IP)
            TMP3  =  DTK * NM
            ASPAR_DIAG(IS,ID)=ASPAR_DIAG(IS,ID) + TRIA03+DTK- TMP3 * DELTAL(IPOS)
            val1=-TMP3*DELTAL(POS_TRICK(IPOS,1))
            val2=-TMP3*DELTAL(POS_TRICK(IPOS,2))
            NEG_P(IS,ID)=NEG_P(IS,ID) + val1*AC2(IS,ID,IP1)
            NEG_P(IS,ID)=NEG_P(IS,ID) + val2*AC2(IS,ID,IP2)
          END DO
        END DO
      END DO
      IF (REFRACTION_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LTHBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPTHETA(IP,CAD)
          CP_THE = MAX(ZERO,CAD)
          CM_THE = MIN(ZERO,CAD)
          eFact=(DT4D/DDIR)*SI(IP)
          DO ID=1,NUMDIR
            ID1 = ID-1
            ID2 = ID+1
            IF (ID == 1) ID1 = NUMDIR
            IF (ID == NUMDIR) ID2 = 1
            NEG_P(:,ID)=NEG_P(:,ID) - eFact*CP_THE(:,ID1)*AC2(:,ID1,IP)
            NEG_P(:,ID)=NEG_P(:,ID) + eFact*CM_THE(:,ID2)*AC2(:,ID2,IP)
          END DO
          ASPAR_DIAG = ASPAR_DIAG + eFact * (CP_THE(:,:) - CM_THE(:,:))
        END IF
      END IF
      IF (FREQ_SHIFT_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LSIGBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPSIGMA(IP,CAS)
          eFact=DT4F*SI(IP)
          DO ID = 1, NUMDIR
            CASS(1:NUMSIG) = CAS(:,ID)
            CASS(0)     = 0.
            CASS(NUMSIG+1) = CASS(NUMSIG)
            CP_SIG = MAX(ZERO,CASS)
            CM_SIG = MIN(ZERO,CASS)
            DO IS=1,NUMSIG
              B_SIG(IS)=eFact*(CP_SIG(IS)/DS_INCR(IS-1) - CM_SIG(IS) /DS_INCR(IS))
            END DO
            DO IS=2,NUMSIG
              NEG_P(IS,ID)=NEG_P(IS,ID) - eFact*CP_SIG(IS-1)/DS_INCR(IS-1)*AC2(IS-1,ID,IP)
            END DO
            DO IS=1,NUMSIG-1
              NEG_P(IS,ID)=NEG_P(IS,ID) + eFact*CM_SIG(IS+1)/DS_INCR(IS)*AC2(IS+1,ID,IP)
            END DO
            B_SIG(NUMSIG) = B_SIG(NUMSIG) + eFact*CM_SIG(NUMSIG+1)/DS_INCR(NUMSIG) * TAIL_ARR(5)
            ASPAR_DIAG(:,ID)=ASPAR_DIAG(:,ID) + B_SIG
          END DO
        END IF
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE COMPUTE_K_CRFS_XYU
      IMPLICIT NONE
      INTEGER :: IP, J, ICON, IPie
      REAL(rkind) :: FL11_X, FL12_X, FL21_X, FL22_X, FL31_X, FL32_X
      REAL(rkind) :: FL11_Y, FL12_Y, FL21_Y, FL22_Y, FL31_Y, FL32_Y
      REAL(rkind) :: FL11_U, FL12_U, FL21_U, FL22_U, FL31_U, FL32_U
      REAL(rkind) :: CSX(3), CSY(3)
      REAL(rkind) :: K_X(3), K_Y(3), CRFS_X(3), CRFS_Y(3)
      REAL(rkind) :: CRFS_U(3), K_U(3), LAMBDA_UX, LAMBDA_UY
      REAL(rkind) :: UV_CUR(3,2)
      INTEGER :: IE, IPOS, I1, I2, I3, I, IP1, IP2, IS
      REAL(rkind) :: LAMBDA_X, LAMBDA_Y
      INTEGER :: POS_TRICK(3,2)
      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2
      J=0
      DO IP=1,NP_RES
        DO ICON = 1, CCON(IP)
          J=J+1
          IE     =  IE_CELL2(IP,ICON)
          IPOS   = POS_CELL2(IP,ICON)
          I1 = INE(1,IE)
          I2 = INE(2,IE)
          I3 = INE(3,IE)
          IF (LSECU .OR. LSTCU) THEN
            IF (LSPHE) THEN
              DO I=1,3
                IPie=INE(I,IE)
                UV_CUR(I,:)=CURTXY(IPie,:)*INVSPHTRANS(IPie,:)
              END DO
            ELSE
              DO I=1,3
                IPie=INE(I,IE)
                UV_CUR(I,:)=CURTXY(IPie,:)
              END DO
            END IF
            LAMBDA_UX=ONESIXTH*(UV_CUR(1,1)+UV_CUR(2,1)+UV_CUR(3,1))
            LAMBDA_UY=ONESIXTH*(UV_CUR(1,2)+UV_CUR(2,2)+UV_CUR(3,2))
            K_U(1)  = LAMBDA_UX * IEN(1,IE) + LAMBDA_UY * IEN(2,IE)
            K_U(2)  = LAMBDA_UX * IEN(3,IE) + LAMBDA_UY * IEN(4,IE)
            K_U(3)  = LAMBDA_UX * IEN(5,IE) + LAMBDA_UY * IEN(6,IE)
            FL11_U = UV_CUR(2,1)*IEN(1,IE)+UV_CUR(2,2)*IEN(2,IE)
            FL12_U = UV_CUR(3,1)*IEN(1,IE)+UV_CUR(3,2)*IEN(2,IE)
            FL21_U = UV_CUR(3,1)*IEN(3,IE)+UV_CUR(3,2)*IEN(4,IE)
            FL22_U = UV_CUR(1,1)*IEN(3,IE)+UV_CUR(1,2)*IEN(4,IE)
            FL31_U = UV_CUR(1,1)*IEN(5,IE)+UV_CUR(1,2)*IEN(6,IE)
            FL32_U = UV_CUR(2,1)*IEN(5,IE)+UV_CUR(2,2)*IEN(6,IE)
            CRFS_U(1) = - ONESIXTH*(TWO *FL31_U + FL32_U + FL21_U + TWO * FL22_U)
            CRFS_U(2) = - ONESIXTH*(TWO *FL32_U + TWO * FL11_U + FL12_U + FL31_U)
            CRFS_U(3) = - ONESIXTH*(TWO *FL12_U + TWO * FL21_U + FL22_U + FL11_U)
          ELSE
            K_U=ZERO
            CRFS_U=ZERO
          END IF
          K_CRFS_U(1:3,J)=K_U
          K_CRFS_U(4:6,J)=CRFS_U
          IP1=INE(POS_TRICK(IPOS,1),IE)
          IP2=INE(POS_TRICK(IPOS,2),IE)
          DO IS=1,NUMSIG
            IF (LSPHE) THEN
              DO I=1,3
                IPie=INE(I,IE)
                CSX(I)=CG(IS,IPie)*INVSPHTRANS(IPie,1)
                CSY(I)=CG(IS,IPie)*INVSPHTRANS(IPie,2)
              END DO
            ELSE
              DO I=1,3
                IPie=INE(I,IE)
                CSX(I)=CG(IS,IPie)
                CSY(I)=CG(IS,IPie)
              END DO
            END IF
            LAMBDA_X=ONESIXTH * (CSX(1) + CSX(2) + CSX(3))
            LAMBDA_Y=ONESIXTH * (CSY(1) + CSY(2) + CSY(3))
            K_X(1)  = LAMBDA_X * IEN(1,IE)
            K_X(2)  = LAMBDA_X * IEN(3,IE)
            K_X(3)  = LAMBDA_X * IEN(5,IE)
            K_Y(1)  = LAMBDA_Y * IEN(2,IE)
            K_Y(2)  = LAMBDA_Y * IEN(4,IE)
            K_Y(3)  = LAMBDA_Y * IEN(6,IE)

            FL11_X = CSX(2)*IEN(1,IE)
            FL12_X = CSX(3)*IEN(1,IE)
            FL21_X = CSX(3)*IEN(3,IE)
            FL22_X = CSX(1)*IEN(3,IE)
            FL31_X = CSX(1)*IEN(5,IE)
            FL32_X = CSX(2)*IEN(5,IE)
            FL11_Y = CSY(2)*IEN(2,IE)
            FL12_Y = CSY(3)*IEN(2,IE)
            FL21_Y = CSY(3)*IEN(4,IE)
            FL22_Y = CSY(1)*IEN(4,IE)
            FL31_Y = CSY(1)*IEN(6,IE)
            FL32_Y = CSY(2)*IEN(6,IE)

            CRFS_X(1)= - ONESIXTH*(TWO *FL31_X + FL32_X + FL21_X + TWO * FL22_X )
            CRFS_X(2)= - ONESIXTH*(TWO *FL32_X + TWO * FL11_X + FL12_X + FL31_X )
            CRFS_X(3)= - ONESIXTH*(TWO *FL12_X + TWO * FL21_X + FL22_X + FL11_X )
            CRFS_Y(1)= - ONESIXTH*(TWO *FL31_Y + FL32_Y + FL21_Y + TWO * FL22_Y )
            CRFS_Y(2)= - ONESIXTH*(TWO *FL32_Y + TWO * FL11_Y + FL12_Y + FL31_Y )
            CRFS_Y(3)= - ONESIXTH*(TWO *FL12_Y + TWO * FL21_Y + FL22_Y + FL11_Y )
            K_CRFS_NUMSIG(1:3,IS,J)=K_X
            K_CRFS_NUMSIG(4:6,IS,J)=K_Y
            K_CRFS_NUMSIG(7:9,IS,J)=CRFS_X
            K_CRFS_NUMSIG(10:12,IS,J)=CRFS_Y
          END DO
        END DO
      END DO
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE NEGATIVE_PART_C(J, IP, NEG_P, ASPAR_DIAG)
      IMPLICIT NONE
      INTEGER, intent(in) :: IP
      INTEGER, intent(inout) :: J
      REAL(rkind), intent(out) :: NEG_P(NUMSIG,NUMDIR)
      REAL(rkind), intent(out) :: ASPAR_DIAG(NUMSIG,NUMDIR)
      REAL(rkind) :: K(NUMSIG,NUMDIR,3), CRFS(NUMSIG,NUMDIR,3)
      REAL(rkind) :: DELTAL(NUMSIG,NUMDIR,3)
      REAL(rkind) :: KM(NUMSIG,NUMDIR,3), KP(NUMSIG,NUMDIR,3)
      REAL(rkind) :: NM(NUMSIG,NUMDIR), K1(NUMSIG,NUMDIR)
      REAL(rkind) :: DTK(NUMSIG,NUMDIR), TMP3(NUMSIG,NUMDIR)
      REAL(rkind) :: eF(NUMSIG,NUMDIR)
      REAL(rkind) :: CAD(NUMSIG,NUMDIR)
      REAL(rkind) :: CAS(NUMSIG,NUMDIR)
      REAL(rkind) :: CP_THE(NUMSIG,NUMDIR), CM_THE(NUMSIG,NUMDIR)
      REAL(rkind) :: CASS(0:NUMSIG+1), B_SIG(NUMSIG)
      REAL(rkind) :: CP_SIG(0:NUMSIG+1), CM_SIG(0:NUMSIG+1)
      INTEGER ICON, ID, IS, idx, IE, IPOS, IP1, IP2, TheVal
      INTEGER ID1, ID2
      INTEGER :: POS_TRICK(3,2)
      REAL(rkind) :: TRIA03
      REAL(rkind) :: eFact
      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2
      NEG_P=ZERO
      ASPAR_DIAG=ZERO
      DO ICON = 1, CCON(IP)
        J=J+1
        IE     =  IE_CELL2(IP,ICON)
        IPOS   = POS_CELL2(IP,ICON)
        DO ID=1,NUMDIR
          DO idx=1,3
            K(:,ID,idx)=K_CRFS_NUMSIG(idx,:,J)*COSTH(ID) + K_CRFS_NUMSIG(idx+3,:,J)*SINTH(ID) + K_CRFS_U(idx,J)
            CRFS(:,ID,idx)=K_CRFS_NUMSIG(idx+6,:,J)*COSTH(ID) + K_CRFS_NUMSIG(idx+9,:,J)*SINTH(ID) + K_CRFS_U(idx+3,J)
          END DO
        END DO
        KM = MIN(ZERO,K)
        KP = MAX(ZERO,K)
        DELTAL = CRFS - KP
        NM(:,:)=ONE/MIN(-THR,KM(:,:,1) + KM(:,:,2) + KM(:,:,3))
        TRIA03 = ONETHIRD * TRIA(IE)
        !
        IP1=INE(POS_TRICK(IPOS,1),IE)
        IP2=INE(POS_TRICK(IPOS,2),IE)
        K1(:,:) =  KP(:,:,IPOS)
        DO ID=1,NUMDIR
          DTK(:,ID) =  K1(:,ID) * DT4A * IOBPD(ID,IP) * IOBWB(IP) * IOBDP(IP)
        END DO
        TMP3(:,:)  =  DTK(:,:) * NM(:,:)
        ASPAR_DIAG=ASPAR_DIAG + TRIA03+DTK(:,:)- TMP3(:,:) * DELTAL(:,:,IPOS)
        eF(:,:) = -TMP3(:,:)*DELTAL(:,:,POS_TRICK(IPOS,1))
        NEG_P=NEG_P  + eF(:,:)*AC2(:,:,IP1)
        eF(:,:) = -TMP3(:,:)*DELTAL(:,:,POS_TRICK(IPOS,2))
        NEG_P=NEG_P  + eF(:,:)*AC2(:,:,IP2)
      END DO
      IF (REFRACTION_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LTHBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPTHETA(IP,CAD)
          CP_THE = MAX(ZERO,CAD)
          CM_THE = MIN(ZERO,CAD)
          eFact=(DT4D/DDIR)*SI(IP)
          DO ID=1,NUMDIR
            ID1 = ID-1
            ID2 = ID+1
            IF (ID == 1) ID1 = NUMDIR
            IF (ID == NUMDIR) ID2 = 1
            NEG_P(:,ID)=NEG_P(:,ID) - eFact*CP_THE(:,ID1)*AC2(:,ID1,IP)
            NEG_P(:,ID)=NEG_P(:,ID) + eFact*CM_THE(:,ID2)*AC2(:,ID2,IP)
          END DO
          ASPAR_DIAG = ASPAR_DIAG + eFact * (CP_THE(:,:) - CM_THE(:,:))
        END IF
      END IF
      IF (FREQ_SHIFT_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LSIGBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPSIGMA(IP,CAS)
          eFact=DT4F*SI(IP)
          DO ID = 1, NUMDIR
            CASS(1:NUMSIG) = CAS(:,ID)
            CASS(0)     = 0.
            CASS(NUMSIG+1) = CASS(NUMSIG)
            CP_SIG = MAX(ZERO,CASS)
            CM_SIG = MIN(ZERO,CASS)
            DO IS=1,NUMSIG
              B_SIG(IS)=eFact*(CP_SIG(IS)/DS_INCR(IS-1) - CM_SIG(IS) /DS_INCR(IS))
            END DO
            DO IS=2,NUMSIG
              NEG_P(IS,ID)=NEG_P(IS,ID) - eFact*CP_SIG(IS-1)/DS_INCR(IS-1)*AC2(IS-1,ID,IP)
            END DO
            DO IS=1,NUMSIG-1
              NEG_P(IS,ID)=NEG_P(IS,ID) + eFact*CM_SIG(IS+1)/DS_INCR(IS)*AC2(IS+1,ID,IP)
            END DO
            B_SIG(NUMSIG) = B_SIG(NUMSIG) + eFact*CM_SIG(NUMSIG+1)/DS_INCR(NUMSIG) * TAIL_ARR(5)
            ASPAR_DIAG(:,ID)=ASPAR_DIAG(:,ID) + B_SIG
          END DO
        END IF
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE NEGATIVE_PART_D(J, IP, NEG_P, ASPAR_DIAG)
      IMPLICIT NONE
      INTEGER, intent(in) :: IP
      INTEGER, intent(inout) :: J
      REAL(rkind), intent(out) :: NEG_P(NUMSIG,NUMDIR)
      REAL(rkind), intent(out) :: ASPAR_DIAG(NUMSIG,NUMDIR)
      REAL(rkind) :: K(NUMSIG,3), CRFS(NUMSIG,3)
      REAL(rkind) :: DELTAL(NUMSIG,3)
      REAL(rkind) :: KM(NUMSIG,3), KP(NUMSIG,3)
      REAL(rkind) :: NM(NUMSIG)
      REAL(rkind) :: DTK(NUMSIG), TMP3(NUMSIG)
      REAL(rkind) :: eF(NUMSIG)
      REAL(rkind) :: CAD(NUMSIG,NUMDIR)
      REAL(rkind) :: CAS(NUMSIG,NUMDIR)
      REAL(rkind) :: CP_THE(NUMSIG,NUMDIR), CM_THE(NUMSIG,NUMDIR)
      REAL(rkind) :: CASS(0:NUMSIG+1), B_SIG(NUMSIG)
      REAL(rkind) :: CP_SIG(0:NUMSIG+1), CM_SIG(0:NUMSIG+1)
      INTEGER ICON, ID, IS, idx, IE, IPOS, IP1, IP2, TheVal
      INTEGER ID1, ID2
      INTEGER :: POS_TRICK(3,2)
      REAL(rkind) :: TRIA03
      REAL(rkind) :: eFact
      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2
      NEG_P=ZERO
      ASPAR_DIAG=ZERO
      DO ICON = 1, CCON(IP)
        J=J+1
        IE     =  IE_CELL2(IP,ICON)
        IPOS   = POS_CELL2(IP,ICON)
        DO ID=1,NUMDIR
          DO idx=1,3
            K(:,idx)=K_CRFS_NUMSIG(idx,:,J)*COSTH(ID) + K_CRFS_NUMSIG(idx+3,:,J)*SINTH(ID) + K_CRFS_U(idx,J)
            CRFS(:,idx)=K_CRFS_NUMSIG(idx+6,:,J)*COSTH(ID) + K_CRFS_NUMSIG(idx+9,:,J)*SINTH(ID) + K_CRFS_U(idx+3,J)
          END DO
          KM = MIN(ZERO,K)
          KP = MAX(ZERO,K)
          DELTAL = CRFS - KP
          NM(:)=ONE/MIN(-THR,KM(:,1) + KM(:,2) + KM(:,3))
          TRIA03 = ONETHIRD * TRIA(IE)
          !
          IP1=INE(POS_TRICK(IPOS,1),IE)
          IP2=INE(POS_TRICK(IPOS,2),IE)
          DTK(:) =  KP(:,IPOS) * DT4A * IOBPD(ID,IP) * IOBWB(IP) * IOBDP(IP)
          TMP3(:)  =  DTK(:) * NM(:)
          ASPAR_DIAG(:,ID)=ASPAR_DIAG(:,ID) + TRIA03+DTK(:)- TMP3(:) * DELTAL(:,IPOS)
          eF(:) = -TMP3(:)*DELTAL(:,POS_TRICK(IPOS,1))
          NEG_P(:,ID)=NEG_P(:,ID)  + eF(:)*AC2(:,ID,IP1)
          eF(:) = -TMP3(:)*DELTAL(:,POS_TRICK(IPOS,2))
          NEG_P(:,ID)=NEG_P(:,ID)  + eF(:)*AC2(:,ID,IP2)
        END DO
      END DO
      IF (REFRACTION_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LTHBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPTHETA(IP,CAD)
          CP_THE = MAX(ZERO,CAD)
          CM_THE = MIN(ZERO,CAD)
          eFact=(DT4D/DDIR)*SI(IP)
          DO ID=1,NUMDIR
            ID1 = ID-1
            ID2 = ID+1
            IF (ID == 1) ID1 = NUMDIR
            IF (ID == NUMDIR) ID2 = 1
            NEG_P(:,ID)=NEG_P(:,ID) - eFact*CP_THE(:,ID1)*AC2(:,ID1,IP)
            NEG_P(:,ID)=NEG_P(:,ID) + eFact*CM_THE(:,ID2)*AC2(:,ID2,IP)
          END DO
          ASPAR_DIAG = ASPAR_DIAG + eFact * (CP_THE(:,:) - CM_THE(:,:))
        END IF
      END IF
      IF (FREQ_SHIFT_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LSIGBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPSIGMA(IP,CAS)
          eFact=DT4F*SI(IP)
          DO ID = 1, NUMDIR
            CASS(1:NUMSIG) = CAS(:,ID)
            CASS(0)     = 0.
            CASS(NUMSIG+1) = CASS(NUMSIG)
            CP_SIG = MAX(ZERO,CASS)
            CM_SIG = MIN(ZERO,CASS)
            DO IS=1,NUMSIG
              B_SIG(IS)=eFact*(CP_SIG(IS)/DS_INCR(IS-1) - CM_SIG(IS) /DS_INCR(IS))
            END DO
            DO IS=2,NUMSIG
              NEG_P(IS,ID)=NEG_P(IS,ID) - eFact*CP_SIG(IS-1)/DS_INCR(IS-1)*AC2(IS-1,ID,IP)
            END DO
            DO IS=1,NUMSIG-1
              NEG_P(IS,ID)=NEG_P(IS,ID) + eFact*CM_SIG(IS+1)/DS_INCR(IS)*AC2(IS+1,ID,IP)
            END DO
            B_SIG(NUMSIG) = B_SIG(NUMSIG) + eFact*CM_SIG(NUMSIG+1)/DS_INCR(NUMSIG) * TAIL_ARR(5)
            ASPAR_DIAG(:,ID)=ASPAR_DIAG(:,ID) + B_SIG
          END DO
        END IF
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE NEGATIVE_PART_E(J, IP, NEG_P, ASPAR_DIAG)
      IMPLICIT NONE
      INTEGER, intent(in) :: IP
      INTEGER, intent(inout) :: J
      REAL(rkind), intent(out) :: NEG_P(NUMSIG,NUMDIR)
      REAL(rkind), intent(out) :: ASPAR_DIAG(NUMSIG,NUMDIR)
      REAL(rkind) :: K(3), CRFS(3)
      REAL(rkind) :: DELTAL(3)
      REAL(rkind) :: KM(3), KP(3)
      REAL(rkind) :: NM, DWDH, WKDEP
      REAL(rkind) :: DTK, TMP3
      REAL(rkind) :: CAD(NUMDIR)
      REAL(rkind) :: CAS(NUMSIG,NUMDIR)
      REAL(rkind) :: CP_THE(NUMDIR), CM_THE(NUMDIR)
      REAL(rkind) :: CASS(0:NUMSIG+1), B_SIG(NUMSIG)
      REAL(rkind) :: CP_SIG(0:NUMSIG+1), CM_SIG(0:NUMSIG+1)
      INTEGER ICON, ID, IS, idx, IE, IPOS, IP1, IP2, TheVal
      INTEGER ID1, ID2
      INTEGER :: POS_TRICK(3,2)
      REAL(rkind) :: TRIA03
      REAL(rkind) :: eFact
      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2
      NEG_P=ZERO
      ASPAR_DIAG=ZERO
      DO ICON = 1, CCON(IP)
        J=J+1
        IE     =  IE_CELL2(IP,ICON)
        IPOS   = POS_CELL2(IP,ICON)
        DO ID=1,NUMDIR
          DO IS=1,NUMSIG
            DO idx=1,3
              K(idx)=K_CRFS_NUMSIG(idx,IS,J)*COSTH(ID) + K_CRFS_NUMSIG(idx+3,IS,J)*SINTH(ID) + K_CRFS_U(idx,J)
              CRFS(idx)=K_CRFS_NUMSIG(idx+6,IS,J)*COSTH(ID) + K_CRFS_NUMSIG(idx+9,IS,J)*SINTH(ID) + K_CRFS_U(idx+3,J)
            END DO
            KM = MIN(ZERO,K)
            KP = MAX(ZERO,K)
            DELTAL = CRFS - KP
            NM=ONE/MIN(-THR,SUM(KM))
            TRIA03 = ONETHIRD * TRIA(IE)
            !
            IP1=INE(POS_TRICK(IPOS,1),IE)
            IP2=INE(POS_TRICK(IPOS,2),IE)
            DTK =  KP(IPOS) * DT4A * IOBPD(ID,IP) * IOBWB(IP) * IOBDP(IP)
            TMP3  =  DTK * NM
            ASPAR_DIAG(IS,ID)=ASPAR_DIAG(IS,ID) + TRIA03+DTK - TMP3 * DELTAL(IPOS)
            NEG_P(IS,ID)=NEG_P(IS,ID) -TMP3*DELTAL(POS_TRICK(IPOS,1))*AC2(IS,ID,IP1)
            NEG_P(IS,ID)=NEG_P(IS,ID) - TMP3*DELTAL(POS_TRICK(IPOS,2))*AC2(IS,ID,IP2)
          END DO
        END DO
      END DO
      IF (REFRACTION_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LTHBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          eFact=(DT4D/DDIR)*SI(IP)
          DO IS = 1, NUMSIG
            WKDEP = WK(IS,IP) * DEP(IP)
            IF (WKDEP .LT. 13.) THEN
              DWDH = SPSIG(IS)/SINH(MIN(KDMAX,2.*WKDEP))
              DO ID = 1, NUMDIR
                CAD(ID) = DWDH * ( SINTH(ID)*DDEP(IP,1)-COSTH(ID)*DDEP(IP,2) )
              END DO
            ELSE
              CAD=ZERO
            ENDIF
            IF (LSTCU .OR. LSECU) THEN
              DO ID = 1, NUMDIR
                CAD(ID) = CAD(ID) + SIN2TH(ID)*DCUY(IP,1)-COS2TH(ID)*DCUX(IP,2)+SINCOSTH(ID)*( DCUX(IP,1)-DCUY(IP,2) )
              END DO
            END IF
            CP_THE = MAX(ZERO,CAD)
            CM_THE = MIN(ZERO,CAD)
            DO ID=1,NUMDIR
              ID1 = ID-1
              ID2 = ID+1
              IF (ID == 1) ID1 = NUMDIR
              IF (ID == NUMDIR) ID2 = 1
              NEG_P(IS,ID)=NEG_P(IS,ID) - eFact*CP_THE(ID1)*AC2(IS,ID1,IP)
              NEG_P(IS,ID)=NEG_P(IS,ID) + eFact*CM_THE(ID2)*AC2(IS,ID2,IP)
            END DO
            ASPAR_DIAG(IS,:) = ASPAR_DIAG(IS,:) + eFact*(CP_THE(:) - CM_THE(:))
          END DO
        END IF
      END IF
      IF (FREQ_SHIFT_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LSIGBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPSIGMA(IP,CAS)
          eFact=DT4F*SI(IP)
          DO ID = 1, NUMDIR
            CASS(1:NUMSIG) = CAS(:,ID)
            CASS(0)     = 0.
            CASS(NUMSIG+1) = CASS(NUMSIG)
            CP_SIG = MAX(ZERO,CASS)
            CM_SIG = MIN(ZERO,CASS)
            DO IS=1,NUMSIG
              B_SIG(IS)=eFact*(CP_SIG(IS)/DS_INCR(IS-1) - CM_SIG(IS) /DS_INCR(IS))
            END DO
            DO IS=2,NUMSIG
              NEG_P(IS,ID)=NEG_P(IS,ID) - eFact*CP_SIG(IS-1)/DS_INCR(IS-1)*AC2(IS-1,ID,IP)
            END DO
            DO IS=1,NUMSIG-1
              NEG_P(IS,ID)=NEG_P(IS,ID) + eFact*CM_SIG(IS+1)/DS_INCR(IS)*AC2(IS+1,ID,IP)
            END DO
            B_SIG(NUMSIG) = B_SIG(NUMSIG) + eFact*CM_SIG(NUMSIG+1)/DS_INCR(NUMSIG) * TAIL_ARR(5)
            ASPAR_DIAG(:,ID)=ASPAR_DIAG(:,ID) + B_SIG
          END DO
        END IF
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE NEGATIVE_PART_F(IP, NEG_P, ASPAR_DIAG)
      IMPLICIT NONE
      INTEGER, intent(in) :: IP
      REAL(rkind), intent(out) :: NEG_P(NUMSIG,NUMDIR)
      REAL(rkind), intent(out) :: ASPAR_DIAG(NUMSIG,NUMDIR)
      REAL(rkind) :: FL11_X, FL12_X, FL21_X, FL22_X, FL31_X, FL32_X
      REAL(rkind) :: FL11_Y, FL12_Y, FL21_Y, FL22_Y, FL31_Y, FL32_Y
      REAL(rkind) :: FL11_U, FL12_U, FL21_U, FL22_U, FL31_U, FL32_U
      REAL(rkind) :: CRFS(3), KM(3), K(3), TRIA03
      REAL(rkind) :: DELTAL(3)
      REAL(rkind) :: KP(3), NM, val1, val2
      REAL(rkind) :: K_X(3), K_Y(3), CRFS_X(3), CRFS_Y(3)
      REAL(rkind) :: CSX(3), CSY(3)
      REAL(rkind) :: CRFS_U(3), K_U(3), LAMBDA_UX, LAMBDA_UY
      REAL(rkind) :: UV_CUR(3,2)
      REAL(rkind) :: DTK, TMP3
      REAL(rkind) :: LAMBDA_X, LAMBDA_Y
      INTEGER     :: I1, I2, I3
      INTEGER     :: ID, IS, IE, IPOS
      INTEGER     :: I, ICON
      INTEGER     :: IPie, TheVal, IP1, IP2
      INTEGER     :: ID1, ID2
      REAL(rkind) :: CAD(NUMSIG,NUMDIR)
      REAL(rkind) :: CAS(NUMSIG,NUMDIR)
      REAL(rkind) :: CP_THE(NUMSIG,NUMDIR), CM_THE(NUMSIG,NUMDIR)
      REAL(rkind) :: CASS(0:NUMSIG+1), B_SIG(NUMSIG)
      REAL(rkind) :: CP_SIG(0:NUMSIG+1), CM_SIG(0:NUMSIG+1)
      REAL(rkind) :: eFact
      INTEGER :: POS_TRICK(3,2)
      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2

      NEG_P=ZERO
      ASPAR_DIAG=ZERO
      DO ICON = 1, CCON(IP)
        IE     =  IE_CELL2(IP,ICON)
        IPOS   = POS_CELL2(IP,ICON)
        I1 = INE(1,IE)
        I2 = INE(2,IE)
        I3 = INE(3,IE)
        IF (LSECU .OR. LSTCU) THEN
          IF (LSPHE) THEN
            DO I=1,3
              IPie=INE(I,IE)
              UV_CUR(I,:)=CURTXY(IPie,:)*INVSPHTRANS(IPie,:)
            END DO
          ELSE
            DO I=1,3
              IPie=INE(I,IE)
              UV_CUR(I,:)=CURTXY(IPie,:)
            END DO
          END IF
          LAMBDA_UX=ONESIXTH*(UV_CUR(1,1)+UV_CUR(2,1)+UV_CUR(3,1))
          LAMBDA_UY=ONESIXTH*(UV_CUR(1,2)+UV_CUR(2,2)+UV_CUR(3,2))
          K_U(1)  = LAMBDA_UX * IEN(1,IE) + LAMBDA_UY * IEN(2,IE)
          K_U(2)  = LAMBDA_UX * IEN(3,IE) + LAMBDA_UY * IEN(4,IE)
          K_U(3)  = LAMBDA_UX * IEN(5,IE) + LAMBDA_UY * IEN(6,IE)
          FL11_U = UV_CUR(2,1)*IEN(1,IE)+UV_CUR(2,2)*IEN(2,IE)
          FL12_U = UV_CUR(3,1)*IEN(1,IE)+UV_CUR(3,2)*IEN(2,IE)
          FL21_U = UV_CUR(3,1)*IEN(3,IE)+UV_CUR(3,2)*IEN(4,IE)
          FL22_U = UV_CUR(1,1)*IEN(3,IE)+UV_CUR(1,2)*IEN(4,IE)
          FL31_U = UV_CUR(1,1)*IEN(5,IE)+UV_CUR(1,2)*IEN(6,IE)
          FL32_U = UV_CUR(2,1)*IEN(5,IE)+UV_CUR(2,2)*IEN(6,IE)
          CRFS_U(1) = - ONESIXTH*(TWO *FL31_U + FL32_U + FL21_U + TWO * FL22_U)
          CRFS_U(2) = - ONESIXTH*(TWO *FL32_U + TWO * FL11_U + FL12_U + FL31_U)
          CRFS_U(3) = - ONESIXTH*(TWO *FL12_U + TWO * FL21_U + FL22_U + FL11_U)
        ELSE
          K_U=ZERO
          CRFS_U=ZERO
        END IF
        IP1=INE(POS_TRICK(IPOS,1),IE)
        IP2=INE(POS_TRICK(IPOS,2),IE)
        DO IS=1,NUMSIG
          IF (LSPHE) THEN
            DO I=1,3
              IPie=INE(I,IE)
              CSX(I)=CG(IS,IPie)*INVSPHTRANS(IPie,1)
              CSY(I)=CG(IS,IPie)*INVSPHTRANS(IPie,2)
            END DO
          ELSE
            DO I=1,3
              IPie=INE(I,IE)
              CSX(I)=CG(IS,IPie)
              CSY(I)=CG(IS,IPie)
            END DO
          END IF
          LAMBDA_X=ONESIXTH * (CSX(1) + CSX(2) + CSX(3))
          LAMBDA_Y=ONESIXTH * (CSY(1) + CSY(2) + CSY(3))
          K_X(1)  = LAMBDA_X * IEN(1,IE)
          K_X(2)  = LAMBDA_X * IEN(3,IE)
          K_X(3)  = LAMBDA_X * IEN(5,IE)
          K_Y(1)  = LAMBDA_Y * IEN(2,IE)
          K_Y(2)  = LAMBDA_Y * IEN(4,IE)
          K_Y(3)  = LAMBDA_Y * IEN(6,IE)

          FL11_X = CSX(2)*IEN(1,IE)
          FL12_X = CSX(3)*IEN(1,IE)
          FL21_X = CSX(3)*IEN(3,IE)
          FL22_X = CSX(1)*IEN(3,IE)
          FL31_X = CSX(1)*IEN(5,IE)
          FL32_X = CSX(2)*IEN(5,IE)
          FL11_Y = CSY(2)*IEN(2,IE)
          FL12_Y = CSY(3)*IEN(2,IE)
          FL21_Y = CSY(3)*IEN(4,IE)
          FL22_Y = CSY(1)*IEN(4,IE)
          FL31_Y = CSY(1)*IEN(6,IE)
          FL32_Y = CSY(2)*IEN(6,IE)

          CRFS_X(1)= - ONESIXTH*(TWO *FL31_X + FL32_X + FL21_X + TWO * FL22_X )
          CRFS_X(2)= - ONESIXTH*(TWO *FL32_X + TWO * FL11_X + FL12_X + FL31_X )
          CRFS_X(3)= - ONESIXTH*(TWO *FL12_X + TWO * FL21_X + FL22_X + FL11_X )
          CRFS_Y(1)= - ONESIXTH*(TWO *FL31_Y + FL32_Y + FL21_Y + TWO * FL22_Y )
          CRFS_Y(2)= - ONESIXTH*(TWO *FL32_Y + TWO * FL11_Y + FL12_Y + FL31_Y )
          CRFS_Y(3)= - ONESIXTH*(TWO *FL12_Y + TWO * FL21_Y + FL22_Y + FL11_Y )

          DO ID=1,NUMDIR
            DO I=1,3
              K(I)=K_X(I)*COSTH(ID) + K_Y(I)*SINTH(ID) + K_U(I)
              CRFS(I)=CRFS_X(I)*COSTH(ID) + CRFS_Y(I)*SINTH(ID) + CRFS_U(I)
            END DO
            KM = MIN(ZERO,K)
            KP = MAX(ZERO,K)
            DELTAL = CRFS - KP
            NM=ONE/MIN(-THR,SUM(KM))
            TRIA03 = ONETHIRD * TRIA(IE)
            !
            DTK =  KP(IPOS) * DT4A * IOBPD(ID,IP) * IOBWB(IP) * IOBDP(IP)
            TMP3  =  DTK * NM
            ASPAR_DIAG(IS,ID)=ASPAR_DIAG(IS,ID) + TRIA03+DTK- TMP3 * DELTAL(IPOS)
            val1=-TMP3*DELTAL(POS_TRICK(IPOS,1))
            val2=-TMP3*DELTAL(POS_TRICK(IPOS,2))
            NEG_P(IS,ID)=NEG_P(IS,ID) + val1*AC2(IS,ID,IP1)
            NEG_P(IS,ID)=NEG_P(IS,ID) + val2*AC2(IS,ID,IP2)
          END DO
        END DO
      END DO
      IF (REFRACTION_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LTHBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPTHETA(IP,CAD)
          CP_THE = MAX(ZERO,CAD)
          CM_THE = MIN(ZERO,CAD)
          eFact=(DT4D/DDIR)*SI(IP)
          DO ID=1,NUMDIR
            ID1 = ID-1
            ID2 = ID+1
            IF (ID == 1) ID1 = NUMDIR
            IF (ID == NUMDIR) ID2 = 1
            NEG_P(:,ID)=NEG_P(:,ID) - eFact*CP_THE(:,ID1)*AC2(:,ID1,IP)
            NEG_P(:,ID)=NEG_P(:,ID) + eFact*CM_THE(:,ID2)*AC2(:,ID2,IP)
          END DO
          ASPAR_DIAG = ASPAR_DIAG + eFact * (CP_THE(:,:) - CM_THE(:,:))
        END IF
      END IF
      IF (FREQ_SHIFT_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LSIGBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPSIGMA(IP,CAS)
          eFact=DT4F*SI(IP)
          DO ID = 1, NUMDIR
            CASS(1:NUMSIG) = CAS(:,ID)
            CASS(0)     = 0.
            CASS(NUMSIG+1) = CASS(NUMSIG)
            CP_SIG = MAX(ZERO,CASS)
            CM_SIG = MIN(ZERO,CASS)
            DO IS=1,NUMSIG
              B_SIG(IS)=eFact*(CP_SIG(IS)/DS_INCR(IS-1) - CM_SIG(IS) /DS_INCR(IS))
            END DO
            DO IS=2,NUMSIG
              NEG_P(IS,ID)=NEG_P(IS,ID) - eFact*CP_SIG(IS-1)/DS_INCR(IS-1)*AC2(IS-1,ID,IP)
            END DO
            DO IS=1,NUMSIG-1
              NEG_P(IS,ID)=NEG_P(IS,ID) + eFact*CM_SIG(IS+1)/DS_INCR(IS)*AC2(IS+1,ID,IP)
            END DO
            B_SIG(NUMSIG) = B_SIG(NUMSIG) + eFact*CM_SIG(NUMSIG+1)/DS_INCR(NUMSIG) * TAIL_ARR(5)
            ASPAR_DIAG(:,ID)=ASPAR_DIAG(:,ID) + B_SIG
          END DO
        END IF
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE NEGATIVE_PART_G(IP, NEG_P, ASPAR_DIAG)
      IMPLICIT NONE
      INTEGER, intent(in) :: IP
      REAL(rkind), intent(out) :: NEG_P(NUMSIG,NUMDIR)
      REAL(rkind), intent(out) :: ASPAR_DIAG(NUMSIG,NUMDIR)
      REAL(rkind) :: FL11_X, FL12_X, FL21_X, FL22_X, FL31_X, FL32_X
      REAL(rkind) :: FL11_Y, FL12_Y, FL21_Y, FL22_Y, FL31_Y, FL32_Y
      REAL(rkind) :: FL11_U, FL12_U, FL21_U, FL22_U, FL31_U, FL32_U
      REAL(rkind) :: CRFS(3), KM(3), K(3), TRIA03
      REAL(rkind) :: DELTAL(3)
      REAL(rkind) :: KP(3), NM, val1, val2
      REAL(rkind) :: K_X(3,NUMSIG), K_Y(3,NUMSIG), CRFS_X(3,NUMSIG), CRFS_Y(3,NUMSIG)
      REAL(rkind) :: CSX(3), CSY(3)
      REAL(rkind) :: CRFS_U(3), K_U(3), LAMBDA_UX, LAMBDA_UY
      REAL(rkind) :: UV_CUR(3,2)
      REAL(rkind) :: DTK, TMP3
      REAL(rkind) :: LAMBDA_X, LAMBDA_Y
      INTEGER     :: I1, I2, I3
      INTEGER     :: ID, IS, IE, IPOS
      INTEGER     :: I, ICON
      INTEGER     :: IPie, TheVal, IP1, IP2
      INTEGER     :: ID1, ID2
      REAL(rkind) :: CAD(NUMSIG,NUMDIR)
      REAL(rkind) :: CAS(NUMSIG,NUMDIR)
      REAL(rkind) :: CP_THE(NUMSIG,NUMDIR), CM_THE(NUMSIG,NUMDIR)
      REAL(rkind) :: CASS(0:NUMSIG+1), B_SIG(NUMSIG)
      REAL(rkind) :: CP_SIG(0:NUMSIG+1), CM_SIG(0:NUMSIG+1)
      REAL(rkind) :: eFact
      INTEGER :: POS_TRICK(3,2)
      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2

      NEG_P=ZERO
      ASPAR_DIAG=ZERO
      DO ICON = 1, CCON(IP)
        IE     =  IE_CELL2(IP,ICON)
        IPOS   = POS_CELL2(IP,ICON)
        I1 = INE(1,IE)
        I2 = INE(2,IE)
        I3 = INE(3,IE)
        IF (LSECU .OR. LSTCU) THEN
          IF (LSPHE) THEN
            DO I=1,3
              IPie=INE(I,IE)
              UV_CUR(I,:)=CURTXY(IPie,:)*INVSPHTRANS(IPie,:)
            END DO
          ELSE
            DO I=1,3
              IPie=INE(I,IE)
              UV_CUR(I,:)=CURTXY(IPie,:)
            END DO
          END IF
          LAMBDA_UX=ONESIXTH*(UV_CUR(1,1)+UV_CUR(2,1)+UV_CUR(3,1))
          LAMBDA_UY=ONESIXTH*(UV_CUR(1,2)+UV_CUR(2,2)+UV_CUR(3,2))
          K_U(1)  = LAMBDA_UX * IEN(1,IE) + LAMBDA_UY * IEN(2,IE)
          K_U(2)  = LAMBDA_UX * IEN(3,IE) + LAMBDA_UY * IEN(4,IE)
          K_U(3)  = LAMBDA_UX * IEN(5,IE) + LAMBDA_UY * IEN(6,IE)
          FL11_U = UV_CUR(2,1)*IEN(1,IE)+UV_CUR(2,2)*IEN(2,IE)
          FL12_U = UV_CUR(3,1)*IEN(1,IE)+UV_CUR(3,2)*IEN(2,IE)
          FL21_U = UV_CUR(3,1)*IEN(3,IE)+UV_CUR(3,2)*IEN(4,IE)
          FL22_U = UV_CUR(1,1)*IEN(3,IE)+UV_CUR(1,2)*IEN(4,IE)
          FL31_U = UV_CUR(1,1)*IEN(5,IE)+UV_CUR(1,2)*IEN(6,IE)
          FL32_U = UV_CUR(2,1)*IEN(5,IE)+UV_CUR(2,2)*IEN(6,IE)
          CRFS_U(1) = - ONESIXTH*(TWO *FL31_U + FL32_U + FL21_U + TWO * FL22_U)
          CRFS_U(2) = - ONESIXTH*(TWO *FL32_U + TWO * FL11_U + FL12_U + FL31_U)
          CRFS_U(3) = - ONESIXTH*(TWO *FL12_U + TWO * FL21_U + FL22_U + FL11_U)
        ELSE
          K_U=ZERO
          CRFS_U=ZERO
        END IF
        IP1=INE(POS_TRICK(IPOS,1),IE)
        IP2=INE(POS_TRICK(IPOS,2),IE)
        DO IS=1,NUMSIG
          IF (LSPHE) THEN
            DO I=1,3
              IPie=INE(I,IE)
              CSX(I)=CG(IS,IPie)*INVSPHTRANS(IPie,1)
              CSY(I)=CG(IS,IPie)*INVSPHTRANS(IPie,2)
            END DO
          ELSE
            DO I=1,3
              IPie=INE(I,IE)
              CSX(I)=CG(IS,IPie)
              CSY(I)=CG(IS,IPie)
            END DO
          END IF
          LAMBDA_X=ONESIXTH * (CSX(1) + CSX(2) + CSX(3))
          LAMBDA_Y=ONESIXTH * (CSY(1) + CSY(2) + CSY(3))
          K_X(1,IS)  = LAMBDA_X * IEN(1,IE)
          K_X(2,IS)  = LAMBDA_X * IEN(3,IE)
          K_X(3,IS)  = LAMBDA_X * IEN(5,IE)
          K_Y(1,IS)  = LAMBDA_Y * IEN(2,IE)
          K_Y(2,IS)  = LAMBDA_Y * IEN(4,IE)
          K_Y(3,IS)  = LAMBDA_Y * IEN(6,IE)

          FL11_X = CSX(2)*IEN(1,IE)
          FL12_X = CSX(3)*IEN(1,IE)
          FL21_X = CSX(3)*IEN(3,IE)
          FL22_X = CSX(1)*IEN(3,IE)
          FL31_X = CSX(1)*IEN(5,IE)
          FL32_X = CSX(2)*IEN(5,IE)
          FL11_Y = CSY(2)*IEN(2,IE)
          FL12_Y = CSY(3)*IEN(2,IE)
          FL21_Y = CSY(3)*IEN(4,IE)
          FL22_Y = CSY(1)*IEN(4,IE)
          FL31_Y = CSY(1)*IEN(6,IE)
          FL32_Y = CSY(2)*IEN(6,IE)

          CRFS_X(1,IS)= - ONESIXTH*(TWO*FL31_X + FL32_X + FL21_X + TWO * FL22_X)
          CRFS_X(2,IS)= - ONESIXTH*(TWO*FL32_X + TWO * FL11_X + FL12_X + FL31_X)
          CRFS_X(3,IS)= - ONESIXTH*(TWO*FL12_X + TWO * FL21_X + FL22_X + FL11_X)
          CRFS_Y(1,IS)= - ONESIXTH*(TWO*FL31_Y + FL32_Y + FL21_Y + TWO * FL22_Y)
          CRFS_Y(2,IS)= - ONESIXTH*(TWO*FL32_Y + TWO * FL11_Y + FL12_Y + FL31_Y)
          CRFS_Y(3,IS)= - ONESIXTH*(TWO*FL12_Y + TWO * FL21_Y + FL22_Y + FL11_Y)
        END DO
        DO ID=1,NUMDIR
          DO IS=1,NUMSIG
            DO I=1,3
              K(I)=K_X(I,IS)*COSTH(ID) + K_Y(I,IS)*SINTH(ID) + K_U(I)
              CRFS(I)=CRFS_X(I,IS)*COSTH(ID) + CRFS_Y(I,IS)*SINTH(ID) + CRFS_U(I)
            END DO
            KM = MIN(ZERO,K)
            KP = MAX(ZERO,K)
            DELTAL = CRFS - KP
            NM=ONE/MIN(-THR,SUM(KM))
            TRIA03 = ONETHIRD * TRIA(IE)
            !
            DTK =  KP(IPOS) * DT4A * IOBPD(ID,IP) * IOBWB(IP) * IOBDP(IP)
            TMP3  =  DTK * NM
            ASPAR_DIAG(IS,ID)=ASPAR_DIAG(IS,ID) + TRIA03+DTK- TMP3 * DELTAL(IPOS)
            val1=-TMP3*DELTAL(POS_TRICK(IPOS,1))
            val2=-TMP3*DELTAL(POS_TRICK(IPOS,2))
            NEG_P(IS,ID)=NEG_P(IS,ID) + val1*AC2(IS,ID,IP1)
            NEG_P(IS,ID)=NEG_P(IS,ID) + val2*AC2(IS,ID,IP2)
          END DO
        END DO
      END DO
      IF (REFRACTION_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LTHBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPTHETA(IP,CAD)
          CP_THE = MAX(ZERO,CAD)
          CM_THE = MIN(ZERO,CAD)
          eFact=(DT4D/DDIR)*SI(IP)
          DO ID=1,NUMDIR
            ID1 = ID-1
            ID2 = ID+1
            IF (ID == 1) ID1 = NUMDIR
            IF (ID == NUMDIR) ID2 = 1
            NEG_P(:,ID)=NEG_P(:,ID) - eFact*CP_THE(:,ID1)*AC2(:,ID1,IP)
            NEG_P(:,ID)=NEG_P(:,ID) + eFact*CM_THE(:,ID2)*AC2(:,ID2,IP)
          END DO
          ASPAR_DIAG = ASPAR_DIAG + eFact * (CP_THE(:,:) - CM_THE(:,:))
        END IF
      END IF
      IF (FREQ_SHIFT_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LSIGBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPSIGMA(IP,CAS)
          eFact=DT4F*SI(IP)
          DO ID = 1, NUMDIR
            CASS(1:NUMSIG) = CAS(:,ID)
            CASS(0)     = 0.
            CASS(NUMSIG+1) = CASS(NUMSIG)
            CP_SIG = MAX(ZERO,CASS)
            CM_SIG = MIN(ZERO,CASS)
            DO IS=1,NUMSIG
              B_SIG(IS)=eFact*(CP_SIG(IS)/DS_INCR(IS-1) - CM_SIG(IS) /DS_INCR(IS))
            END DO
            DO IS=2,NUMSIG
              NEG_P(IS,ID)=NEG_P(IS,ID) - eFact*CP_SIG(IS-1)/DS_INCR(IS-1)*AC2(IS-1,ID,IP)
            END DO
            DO IS=1,NUMSIG-1
              NEG_P(IS,ID)=NEG_P(IS,ID) + eFact*CM_SIG(IS+1)/DS_INCR(IS)*AC2(IS+1,ID,IP)
            END DO
            B_SIG(NUMSIG) = B_SIG(NUMSIG) + eFact*CM_SIG(NUMSIG+1)/DS_INCR(NUMSIG) * TAIL_ARR(5)
            ASPAR_DIAG(:,ID)=ASPAR_DIAG(:,ID) + B_SIG
          END DO
        END IF
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE NEGATIVE_PART_H(IP, NEG_P, ASPAR_DIAG)
      IMPLICIT NONE
      INTEGER, intent(in) :: IP
      REAL(rkind), intent(out) :: NEG_P(NUMSIG,NUMDIR)
      REAL(rkind), intent(out) :: ASPAR_DIAG(NUMSIG,NUMDIR)
      REAL(rkind) :: FL11_X, FL12_X, FL21_X, FL22_X, FL31_X, FL32_X
      REAL(rkind) :: FL11_Y, FL12_Y, FL21_Y, FL22_Y, FL31_Y, FL32_Y
      REAL(rkind) :: FL11_U, FL12_U, FL21_U, FL22_U, FL31_U, FL32_U
      REAL(rkind) :: CRFS(3), KM(3), K(3), TRIA03
      REAL(rkind) :: DELTAL(3)
      REAL(rkind) :: KP(3), NM, val1, val2
      REAL(rkind) :: K_X(3,NUMSIG), K_Y(3,NUMSIG), CRFS_X(3,NUMSIG), CRFS_Y(3,NUMSIG)
      REAL(rkind) :: CSX(3), CSY(3)
      REAL(rkind) :: CRFS_U(3), K_U(3), LAMBDA_UX, LAMBDA_UY
      REAL(rkind) :: UV_CUR(3,2)
      REAL(rkind) :: DTK, TMP3
      REAL(rkind) :: LAMBDA_X, LAMBDA_Y
      INTEGER     :: I1, I2, I3
      INTEGER     :: ID, IS, IE, IPOS
      INTEGER     :: I, ICON
      INTEGER     :: IPie, TheVal, IP1, IP2
      INTEGER     :: ID1, ID2
      REAL(rkind) :: CAS(NUMSIG,NUMDIR)
      REAL(rkind) :: CASS(0:NUMSIG+1), B_SIG(NUMSIG)
      REAL(rkind) :: CP_SIG(0:NUMSIG+1), CM_SIG(0:NUMSIG+1)
      REAL(rkind) :: eFact, eCAD, eCP_THE, eCM_THE
      REAL(rkind) :: CAD_U(NUMDIR), DWDH(NUMSIG), WKDEP

      INTEGER :: POS_TRICK(3,2)
      POS_TRICK(1,1) = 2
      POS_TRICK(1,2) = 3
      POS_TRICK(2,1) = 3
      POS_TRICK(2,2) = 1
      POS_TRICK(3,1) = 1
      POS_TRICK(3,2) = 2

      NEG_P=ZERO
      ASPAR_DIAG=ZERO
      DO ICON = 1, CCON(IP)
        IE     =  IE_CELL2(IP,ICON)
        IPOS   = POS_CELL2(IP,ICON)
        I1 = INE(1,IE)
        I2 = INE(2,IE)
        I3 = INE(3,IE)
        IF (LSECU .OR. LSTCU) THEN
          IF (LSPHE) THEN
            DO I=1,3
              IPie=INE(I,IE)
              UV_CUR(I,:)=CURTXY(IPie,:)*INVSPHTRANS(IPie,:)
            END DO
          ELSE
            DO I=1,3
              IPie=INE(I,IE)
              UV_CUR(I,:)=CURTXY(IPie,:)
            END DO
          END IF
          LAMBDA_UX=ONESIXTH*(UV_CUR(1,1)+UV_CUR(2,1)+UV_CUR(3,1))
          LAMBDA_UY=ONESIXTH*(UV_CUR(1,2)+UV_CUR(2,2)+UV_CUR(3,2))
          K_U(1)  = LAMBDA_UX * IEN(1,IE) + LAMBDA_UY * IEN(2,IE)
          K_U(2)  = LAMBDA_UX * IEN(3,IE) + LAMBDA_UY * IEN(4,IE)
          K_U(3)  = LAMBDA_UX * IEN(5,IE) + LAMBDA_UY * IEN(6,IE)
          FL11_U = UV_CUR(2,1)*IEN(1,IE)+UV_CUR(2,2)*IEN(2,IE)
          FL12_U = UV_CUR(3,1)*IEN(1,IE)+UV_CUR(3,2)*IEN(2,IE)
          FL21_U = UV_CUR(3,1)*IEN(3,IE)+UV_CUR(3,2)*IEN(4,IE)
          FL22_U = UV_CUR(1,1)*IEN(3,IE)+UV_CUR(1,2)*IEN(4,IE)
          FL31_U = UV_CUR(1,1)*IEN(5,IE)+UV_CUR(1,2)*IEN(6,IE)
          FL32_U = UV_CUR(2,1)*IEN(5,IE)+UV_CUR(2,2)*IEN(6,IE)
          CRFS_U(1) = - ONESIXTH*(TWO *FL31_U + FL32_U + FL21_U + TWO * FL22_U)
          CRFS_U(2) = - ONESIXTH*(TWO *FL32_U + TWO * FL11_U + FL12_U + FL31_U)
          CRFS_U(3) = - ONESIXTH*(TWO *FL12_U + TWO * FL21_U + FL22_U + FL11_U)
        ELSE
          K_U=ZERO
          CRFS_U=ZERO
        END IF
        IP1=INE(POS_TRICK(IPOS,1),IE)
        IP2=INE(POS_TRICK(IPOS,2),IE)
        DO IS=1,NUMSIG
          IF (LSPHE) THEN
            DO I=1,3
              IPie=INE(I,IE)
              CSX(I)=CG(IS,IPie)*INVSPHTRANS(IPie,1)
              CSY(I)=CG(IS,IPie)*INVSPHTRANS(IPie,2)
            END DO
          ELSE
            DO I=1,3
              IPie=INE(I,IE)
              CSX(I)=CG(IS,IPie)
              CSY(I)=CG(IS,IPie)
            END DO
          END IF
          LAMBDA_X=ONESIXTH * (CSX(1) + CSX(2) + CSX(3))
          LAMBDA_Y=ONESIXTH * (CSY(1) + CSY(2) + CSY(3))
          K_X(1,IS)  = LAMBDA_X * IEN(1,IE)
          K_X(2,IS)  = LAMBDA_X * IEN(3,IE)
          K_X(3,IS)  = LAMBDA_X * IEN(5,IE)
          K_Y(1,IS)  = LAMBDA_Y * IEN(2,IE)
          K_Y(2,IS)  = LAMBDA_Y * IEN(4,IE)
          K_Y(3,IS)  = LAMBDA_Y * IEN(6,IE)

          FL11_X = CSX(2)*IEN(1,IE)
          FL12_X = CSX(3)*IEN(1,IE)
          FL21_X = CSX(3)*IEN(3,IE)
          FL22_X = CSX(1)*IEN(3,IE)
          FL31_X = CSX(1)*IEN(5,IE)
          FL32_X = CSX(2)*IEN(5,IE)
          FL11_Y = CSY(2)*IEN(2,IE)
          FL12_Y = CSY(3)*IEN(2,IE)
          FL21_Y = CSY(3)*IEN(4,IE)
          FL22_Y = CSY(1)*IEN(4,IE)
          FL31_Y = CSY(1)*IEN(6,IE)
          FL32_Y = CSY(2)*IEN(6,IE)

          CRFS_X(1,IS)= - ONESIXTH*(TWO*FL31_X + FL32_X + FL21_X + TWO * FL22_X)
          CRFS_X(2,IS)= - ONESIXTH*(TWO*FL32_X + TWO * FL11_X + FL12_X + FL31_X)
          CRFS_X(3,IS)= - ONESIXTH*(TWO*FL12_X + TWO * FL21_X + FL22_X + FL11_X)
          CRFS_Y(1,IS)= - ONESIXTH*(TWO*FL31_Y + FL32_Y + FL21_Y + TWO * FL22_Y)
          CRFS_Y(2,IS)= - ONESIXTH*(TWO*FL32_Y + TWO * FL11_Y + FL12_Y + FL31_Y)
          CRFS_Y(3,IS)= - ONESIXTH*(TWO*FL12_Y + TWO * FL21_Y + FL22_Y + FL11_Y)
        END DO
        DO ID=1,NUMDIR
          DO IS=1,NUMSIG
            DO I=1,3
              K(I)=K_X(I,IS)*COSTH(ID) + K_Y(I,IS)*SINTH(ID) + K_U(I)
              CRFS(I)=CRFS_X(I,IS)*COSTH(ID) + CRFS_Y(I,IS)*SINTH(ID) + CRFS_U(I)
            END DO
            KM = MIN(ZERO,K)
            KP = MAX(ZERO,K)
            DELTAL = CRFS - KP
            NM=ONE/MIN(-THR,SUM(KM))
            TRIA03 = ONETHIRD * TRIA(IE)
            !
            DTK =  KP(IPOS) * DT4A * IOBPD(ID,IP) * IOBWB(IP) * IOBDP(IP)
            TMP3  =  DTK * NM
            ASPAR_DIAG(IS,ID)=ASPAR_DIAG(IS,ID) + TRIA03+DTK- TMP3 * DELTAL(IPOS)
            val1=-TMP3*DELTAL(POS_TRICK(IPOS,1))
            val2=-TMP3*DELTAL(POS_TRICK(IPOS,2))
            NEG_P(IS,ID)=NEG_P(IS,ID) + val1*AC2(IS,ID,IP1)
            NEG_P(IS,ID)=NEG_P(IS,ID) + val2*AC2(IS,ID,IP2)
          END DO
        END DO
      END DO
      IF (REFRACTION_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LTHBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          eFact=(DT4D/DDIR)*SI(IP)
          DO IS = 1, NUMSIG
            WKDEP = WK(IS,IP) * DEP(IP)
            IF (WKDEP .LT. 13.) THEN
              DWDH(IS) = SPSIG(IS)/SINH(MIN(KDMAX,2.*WKDEP))
              DO ID = 1, NUMDIR
              END DO
            ELSE
              DWDH(IS)=ZERO
            ENDIF
          END DO
          IF (LSTCU .OR. LSECU) THEN
            DO ID = 1, NUMDIR
              CAD_U(ID) = SIN2TH(ID)*DCUY(IP,1)-COS2TH(ID)*DCUX(IP,2)+SINCOSTH(ID)*( DCUX(IP,1)-DCUY(IP,2) )
            END DO
          ELSE
            CAD_U=ZERO
          END IF
          DO ID=1,NUMDIR
            ID1 = ID-1
            ID2 = ID+1
            IF (ID == 1) ID1 = NUMDIR
            IF (ID == NUMDIR) ID2 = 1
            DO IS=1,NUMSIG
              eCAD=DWDH(IS) * ( SINTH(ID)*DDEP(IP,1)-COSTH(ID)*DDEP(IP,2) ) + CAD_U(ID)
              eCP_THE=MAX(ZERO,eCAD)
              eCM_THE=MIN(ZERO,eCAD)
              ASPAR_DIAG(IS,ID) = ASPAR_DIAG(IS,ID) + eFact*(eCP_THE - eCM_THE)
              NEG_P(IS,ID2)=NEG_P(IS,ID2) - eFact*eCP_THE*AC2(IS,ID,IP)
              NEG_P(IS,ID1)=NEG_P(IS,ID1) + eFact*eCM_THE*AC2(IS,ID,IP)
            END DO
          END DO
        END IF
      END IF
      IF (FREQ_SHIFT_IMPL) THEN
        TheVal=1
        IF ((ABS(IOBP(IP)) .EQ. 1 .OR. ABS(IOBP(IP)) .EQ. 3) .AND. .NOT. LSIGBOUND) TheVal=0
        IF (DEP(IP) .LT. DMIN) TheVal=0
        IF (IOBP(IP) .EQ. 2) TheVal=0
        IF (TheVal .eq. 1) THEN
          CALL PROPSIGMA(IP,CAS)
          eFact=DT4F*SI(IP)
          DO ID = 1, NUMDIR
            CASS(1:NUMSIG) = CAS(:,ID)
            CASS(0)     = 0.
            CASS(NUMSIG+1) = CASS(NUMSIG)
            CP_SIG = MAX(ZERO,CASS)
            CM_SIG = MIN(ZERO,CASS)
            DO IS=1,NUMSIG
              B_SIG(IS)=eFact*(CP_SIG(IS)/DS_INCR(IS-1) - CM_SIG(IS) /DS_INCR(IS))
            END DO
            DO IS=2,NUMSIG
              NEG_P(IS,ID)=NEG_P(IS,ID) - eFact*CP_SIG(IS-1)/DS_INCR(IS-1)*AC2(IS-1,ID,IP)
            END DO
            DO IS=1,NUMSIG-1
              NEG_P(IS,ID)=NEG_P(IS,ID) + eFact*CM_SIG(IS+1)/DS_INCR(IS)*AC2(IS+1,ID,IP)
            END DO
            B_SIG(NUMSIG) = B_SIG(NUMSIG) + eFact*CM_SIG(NUMSIG+1)/DS_INCR(NUMSIG) * TAIL_ARR(5)
            ASPAR_DIAG(:,ID)=ASPAR_DIAG(:,ID) + B_SIG
          END DO
        END IF
      END IF
      END SUBROUTINE
      
      END SUBROUTINE
