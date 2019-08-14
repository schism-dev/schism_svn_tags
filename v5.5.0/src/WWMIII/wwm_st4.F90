#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ST4_PRE (IP, WALOC, PHI, DPHIDN, SSINE, DSSINE, SSDS, DSSDS, SSNL4, DSSNL4, SSINL)
      USE DATAPOOL
      USE W3SRC4MD
      IMPLICIT NONE

      INTEGER, INTENT(IN)        :: IP
      REAL(rkind), INTENT(IN)    :: WALOC(NUMSIG,NUMDIR)

      REAL(rkind), INTENT(OUT)   :: PHI(NUMSIG,NUMDIR), DPHIDN(NUMSIG,NUMDIR)
      REAL(rkind), INTENT(OUT)   :: SSINE(NUMSIG,NUMDIR), DSSINE(NUMSIG,NUMDIR) 
      REAL(rkind), INTENT(OUT)   :: SSDS(NUMSIG,NUMDIR),DSSDS(NUMSIG,NUMDIR)
      REAL(rkind), INTENT(OUT)   :: SSNL4(NUMSIG,NUMDIR),DSSNL4(NUMSIG,NUMDIR)
      REAL(rkind), INTENT(OUT)   :: SSINL(NUMSIG,NUMDIR)

      INTEGER      :: IS, ID, ITH, IK, IS0

      REAL(rkind)  :: AWW3(NSPEC)
      REAL(rkind)  :: VDDS(NSPEC), VSDS(NSPEC), BRLAMBDA(NSPEC)
      REAL(rkind)  :: WN2(NUMSIG*NUMDIR), WHITECAP(1:4)
      REAL(rkind)  :: VSIN(NSPEC), VDIN(NSPEC)

      REAL(rkind)  :: ETOT, FAVG, FMEAN1, WNMEAN, AS, SUMWALOC, FAVGWS
      REAL(rkind)  :: TAUWAX, TAUWAY, AMAX, FPM, WIND10, WINDTH
      REAL(rkind)  :: HS,SME01,SME10,KME01,KMWAM,KMWAM2

      DO IS = 1, NUMSIG
        DO ID = 1, NUMDIR
          AWW3(ID + (IS-1) * NUMDIR) = WALOC(IS,ID) * CG(IS,IP)
        END DO
      END DO

      DO IK=1, NK
        WN2(1+(IK-1)*NTH) = WK(IK,IP)
      END DO

      DO IK=1, NK
        IS0    = (IK-1)*NTH
        DO ITH=2, NTH
          WN2(ITH+IS0) = WN2(1+IS0)
        END DO
      END DO
!
! wind input
!
      TAUWX(IP)  = ZERO
      TAUWY(IP)  = ZERO               
      SSINL      = ZERO
      NUMSIG_HF(IP) = NUMSIG
      AS         = 0.
      BRLAMBDA   = ZERO

      IF (MESIN .GT. 0) THEN

        CALL SET_WIND( IP, WIND10, WINDTH )
        CALL SET_FRICTION( IP, WALOC, WIND10, WINDTH, FPM )
        LLWS=.TRUE.
#ifdef DEBUGSRC
        WRITE(740+myrank,*) '1: input value USTAR=', UFRIC(IP), ' USTDIR=', USTDIR(IP)
#endif
        CALL W3SPR4 ( AWW3, CG(:,IP), WK(:,IP), EMEAN(IP), FMEAN(IP), FMEAN1, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS(IP))
#ifdef DEBUGSRC
        WRITE(740+myrank,*) '1: out value USTAR=', UFRIC(IP), ' USTDIR=', USTDIR(IP)
        WRITE(740+myrank,*) '1: out value EMEAN=', EMEAN(IP), ' FMEAN=', FMEAN(IP)
        WRITE(740+myrank,*) '1: out value FMEAN1=', FMEAN1, ' WNMEAN=', WNMEAN
        WRITE(740+myrank,*) '1: out value CD=', CD(IP), ' Z0=', Z0(IP)
        WRITE(740+myrank,*) '1: out value ALPHA=', ALPHA_CH(IP), ' FMEANWS=', FMEANWS(IP)
#endif
        IF (EMEAN(IP) .LT. THR .AND. WIND10 .GT. THR) CALL SIN_LIN_CAV(IP,WINDTH,FPM,SSINL)
        CALL W3SIN4 ( IP, AWW3, CG(:,IP), WN2, WIND10, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, VSIN, VDIN, LLWS, BRLAMBDA)
#ifdef DEBUGSRC
        WRITE(740+myrank,*) '1: WINDTH=', WINDTH, ' Z0=', Z0(IP), ' CD=', CD(IP)
        WRITE(740+myrank,*) '1: UFRIC=', UFRIC(IP), 'WIND10=', WIND10, ' RHOAW=', RHOAW
        WRITE(740+myrank,*) '1: TAUWX=', TAUWX(IP), ' TAUWY=', TAUWY(IP)
        WRITE(740+myrank,*) '1: TAUWAX=', TAUWAX, ' TAUWAY=', TAUWAY
        WRITE(740+myrank,*) '1: W3SIN4min/max/sum(VSIN)=', minval(VSIN), maxval(VSIN), sum(VSIN)
        WRITE(740+myrank,*) '1: W3SIN4min/max/sum(VDIN)=', minval(VDIN), maxval(VDIN), sum(VDIN)
#endif
        CALL W3SPR4 ( AWW3, CG(:,IP), WK(:,IP), EMEAN(IP), FMEAN(IP), FMEAN1, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS(IP))  
        CALL W3SIN4 ( IP, AWW3, CG(:,IP), WN2, WIND10, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, VSIN, VDIN, LLWS, BRLAMBDA)
#ifdef DEBUGSRC
        WRITE(740+myrank,*) '2: W3SIN4min/max/sum(VSIN)=', minval(VSIN), maxval(VSIN), sum(VSIN)
        WRITE(740+myrank,*) '2: W3SIN4min/max/sum(VDIN)=', minval(VDIN), maxval(VDIN), sum(VDIN)
#endif
        CALL CONVERT_VS_VD_WWM(IP, VSIN, VDIN, SSINE, DSSINE)
      ENDIF

      IF (MESNL .GT. 0) THEN
         CALL MEAN_WAVE_PARAMETER(IP,WALOC,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2)
         CALL DIASNL4WW3(IP, KMWAM, WALOC, SSNL4, DSSNL4)
      END IF

      IF (MESDS .GT. 0) THEN
        CALL W3SDS4(AWW3,WK(:,IP),CG(:,IP),UFRIC(IP),USTDIR(IP),DEP(IP),VSDS,VDDS,BRLAMBDA,WHITECAP)

#ifdef DEBUGSRC
        WRITE(740+myrank,*) '2: W3SDS4min/max/sum(VSDS)=', minval(VSDS), maxval(VSDS), sum(VSDS)
        WRITE(740+myrank,*) '2: W3SDS4min/max/sum(VDDS)=', minval(VDDS), maxval(VDDS), sum(VDDS)
#endif
        CALL CONVERT_VS_VD_WWM(IP, VSDS, VDDS, SSDS, DSSDS)
      ENDIF
!
      PHI    = SSINL + SSINE  + SSNL4  + SSDS
      DPHIDN =         DSSINE + DSSNL4 + DSSDS
!
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE ST4_POST (IP, WALOC, SSINE, DSSINE, SSDS, DSSDS, SSINL)
        USE DATAPOOL
        USE W3SRC4MD
        IMPLICIT NONE
       
        INTEGER, INTENT(IN)        :: IP
        REAL(rkind), INTENT(IN)    :: WALOC(NUMSIG,NUMDIR)
       
        REAL(rkind), INTENT(OUT)   :: SSINE(NUMSIG,NUMDIR),DSSINE(NUMSIG,NUMDIR)
        REAL(rkind), INTENT(OUT)   :: SSDS(NUMSIG,NUMDIR),DSSDS(NUMSIG,NUMDIR)
        REAL(rkind), INTENT(OUT)   :: SSINL(NUMSIG,NUMDIR)

        INTEGER                    :: IS, ID, IK, ITH, ITH2, IS0

        REAL(rkind)                :: PHI(NUMSIG,NUMDIR), DPHIDN(NUMSIG,NUMDIR)
        REAL(rkind)                :: AWW3(NSPEC), WN2(NUMSIG*NUMDIR), BRLAMBDA(NSPEC)
        REAL(rkind)                :: DPHIDN1D(NSPEC), PHI1D(NSPEC), TMP_DS(NUMSIG)

        REAL(rkind)                :: ETOT, FAVG, FMEAN1, WNMEAN, AS, FAVGWS
        REAL(rkind)                :: TAUWAX, TAUWAY, AMAX, WIND10, WINDTH
        REAL(rkind)                :: WHITECAP(1:4), SUMWALOC, FPM

        DO IS = 1, NUMSIG
          DO ID = 1, NUMDIR
            AWW3(ID + (IS-1) * NUMDIR) = WALOC(IS,ID) * CG(IS,IP)
          END DO
        END DO

        DO IK=1, NK
          WN2(1+(IK-1)*NTH) = WK(IK,IP)
        END DO
        DO IK=1, NK
          IS0    = (IK-1)*NTH
          DO ITH=2, NTH
            WN2(ITH+IS0) = WN2(1+IS0)
          END DO
        END DO
!
! wind input
!
        AS      = 0.
        NUMSIG_HF(IP) = NUMSIG
        CALL SET_WIND( IP, WIND10, WINDTH )
        CALL SET_FRICTION( IP, WALOC, WIND10, WINDTH, FPM )
        CALL W3SPR4 ( AWW3, CG(:,IP), WK(:,IP), EMEAN(IP), FMEAN(IP), FMEAN1, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS(IP))
        IF (EMEAN(IP) .LT. THR .AND. WIND10 .GT. THR) CALL SIN_LIN_CAV(IP,WINDTH,FPM,SSINL)
        CALL W3SIN4 ( IP, AWW3, CG(:,IP), WN2,  WIND10, UFRIC(IP), RHOAW, AS, WINDTH, Z0(IP), CD(IP), TAUWX(IP), TAUWY(IP), TAUWAX, TAUWAY, PHI1D, DPHIDN1D, LLWS, BRLAMBDA)
        CALL W3SPR4 ( AWW3, CG(:,IP), WK(:,IP), EMEAN(IP), FMEAN(IP), FMEAN1, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS(IP))
        CALL CONVERT_VS_VD_WWM(IP, PHI1D, DPHIDN1D, SSINE, DSSINE)
!
! dissipation 
!
        CALL W3SDS4(AWW3,WK(:,IP),CG(:,IP),UFRIC(IP),USTDIR(IP),DEP(IP),PHI1D,DPHIDN1D,BRLAMBDA,WHITECAP)
        CALL CONVERT_VS_VD_WWM(IP,PHI1D,DPHIDN1D,SSDS,DSSDS)
!
! missing high freq. tail contribution -> 2do
!
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
