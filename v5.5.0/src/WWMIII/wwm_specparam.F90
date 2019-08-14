#include "wwm_functions.h"
!     Last change:  1     9 Jun 2004    1:44 am
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE STOKES_DRIFT_SURFACE_BAROTROPIC(IP,STOKESBOTTX,STOKESBOTTY,STOKESSURFX,STOKESSURFY,STOKESBAROX,STOKESBAROY)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP

         REAL(rkind), INTENT(OUT)   :: STOKESBOTTX, STOKESBOTTY
         REAL(rkind), INTENT(OUT)   :: STOKESSURFX, STOKESSURFY
         REAL(rkind), INTENT(OUT)   :: STOKESBAROX, STOKESBAROY
         INTEGER                    :: ID, IS
         REAL(rkind)                :: eQuot1, eQuot2, eProd1, eProd2
         REAL(rkind)                :: eMult, eWk, kD, eWkReal
         REAL(rkind)                :: eSinh2kd, eSinhkd, eSinhkd2
         REAL(rkind)                :: eSigma, eLoc, eUint, eVint
         REAL(rkind)                :: eDep
         REAL(rkind)                :: eProd0, eQuot0

         STOKESBOTTX=0
         STOKESBOTTY=0
         STOKESSURFX=0
         STOKESSURFY=0
         STOKESBAROX=0
         STOKESBAROY=0
         eDep=DEP(IP)
         DO IS=1,NUMSIG
           eMult=SPSIG(IS)*DDIR*DS_INCR(IS)
           eWk=WK(IS,IP)
           kD=MIN(KDMAX, eWk*eDep)
           eWkReal=kD/eDep
           eSinh2kd=MySINH(2*kD)
           eSinhkd=MySINH(kD)
           eSinhkd2=eSinhkd**2
           eSigma=SPSIG(IS)
           eUint=0
           eVint=0
           DO ID=1,NUMDIR
             eLoc=AC2(IS,ID,IP)*eMult
             eUint=eUint + eLoc*COSTH(ID)
             eVint=eVint + eLoc*SINTH(ID)
           END DO
           eQuot0=ONE/eSinhkd2
           eProd0=eSigma*eWkReal*eQuot0
           eQuot1=MyCOSH(2*kD)/eSinhkd2
           eProd1=eSigma*eWkReal*eQuot1
           eQuot2=(eSinh2kd/(2*kD))/eSinhkd2
           eProd2=eSigma*eWkReal*eQuot2
           STOKESBOTTX = STOKESBOTTX + eUint*eProd0
           STOKESBOTTY = STOKESBOTTY + eVint*eProd0
           STOKESSURFX = STOKESSURFX + eUint*eProd1
           STOKESSURFY = STOKESSURFY + eVint*eProd1
           STOKESBAROX = STOKESBAROX + eUint*eProd2
           STOKESBAROY = STOKESBAROY + eVint*eProd2
         END DO
      END SUBROUTINE
!**********************************************************************:
!*                                                                    *
!**********************************************************************
      SUBROUTINE STOKES_DRIFT_SURFACE_BAROTROPIC_LOC(WALOC,DEPLOC,WKLOC,  &
     &    STOKESBOTTX,STOKESBOTTY,STOKESSURFX,STOKESSURFY,STOKESBAROX,STOKESBAROY)

         USE DATAPOOL
         IMPLICIT NONE

         REAL(rkind), INTENT(IN)    :: WALOC(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(IN)    :: DEPLOC
         REAL(rkind), INTENT(IN)    :: WKLOC(NUMSIG)
         REAL(rkind), INTENT(OUT)   :: STOKESBOTTX, STOKESBOTTY
         REAL(rkind), INTENT(OUT)   :: STOKESSURFX, STOKESSURFY
         REAL(rkind), INTENT(OUT)   :: STOKESBAROX, STOKESBAROY
         INTEGER             :: ID, IS
         REAL(rkind)              :: eQuot1, eQuot2, eProd1, eProd2
         REAL(rkind)              :: eMult, eWk, kD, eWkReal
         REAL(rkind)              :: eSinh2kd, eSinhkd, eSinhkd2
         REAL(rkind)              :: eSigma, eLoc, eUint, eVint
         REAL(rkind)              :: eDep
         REAL(rkind)                :: eProd0, eQuot0

         STOKESBOTTX=0
         STOKESBOTTY=0
         STOKESSURFX=0
         STOKESSURFY=0
         STOKESBAROX=0
         STOKESBAROY=0
         eDep=DEPLOC
         DO IS=1,NUMSIG
           eMult=SPSIG(IS)*DDIR*DS_INCR(IS)
           eWk=WKLOC(IS)
           kD=MIN(KDMAX, eWk*eDep)
           eWkReal=kD/eDep
           eSinh2kd=MySINH(2*kD)
           eSinhkd=MySINH(kD)
           eSinhkd2=eSinhkd**2
           eSigma=SPSIG(IS)
           eUint=0
           eVint=0
           DO ID=1,NUMDIR
             eLoc=WALOC(IS,ID)*eMult
             eUint=eUint + eLoc*COSTH(ID)
             eVint=eVint + eLoc*SINTH(ID)
           END DO
           eQuot0=ONE/eSinhkd2
           eProd0=eSigma*eWkReal*eQuot0
           eQuot1=MyCOSH(2*kD)/eSinhkd2
           eProd1=eSigma*eWkReal*eQuot1
           eQuot2=(eSinh2kd/(2*kD))/eSinhkd2
           eProd2=eSigma*eWkReal*eQuot2
           STOKESBOTTX = STOKESBOTTX + eUint*eProd0
           STOKESBOTTY = STOKESBOTTY + eVint*eProd0
           STOKESSURFX = STOKESSURFX + eUint*eProd1
           STOKESSURFY = STOKESSURFY + eVint*eProd1
           STOKESBAROX = STOKESBAROX + eUint*eProd2
           STOKESBAROY = STOKESBAROY + eVint*eProd2
         END DO
      END SUBROUTINE
!**********************************************************************:
!*                                                                    *
!**********************************************************************
      SUBROUTINE MEAN_WAVE_PARAMETER(IP,WALOC,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2)
!AR: This must be replaced by ST4_PRE or the certain WAM routine we are not consistent here this routine needs urgen deletion 

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP

         REAL(rkind), INTENT(IN)    :: WALOC(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(OUT)   :: SME01, SME10, KME01
         REAL(rkind), INTENT(OUT)   :: KMWAM, KMWAM2
         REAL(rkind), INTENT(OUT)   :: HS

         INTEGER             :: ID, IS

         REAL(rkind)                :: ACTOT, ETOT
         REAL(rkind)                :: ETOT_SPSIG
         REAL(rkind)                :: ETOT_WK
         REAL(rkind)                :: ETOT_ISQ_WK
         REAL(rkind)                :: ETOT_SQ_WK

         REAL(rkind)                :: Y(NUMSIG), tmp(NUMSIG)
         REAL(rkind)                :: DS, ATAIL, ETAIL, ESIGTAIL
!         REAL(rkind)                :: dintspec, dintspec_y
!
! total energy ...
! 2do improve efficiency ... 
! 2do check integration style ... 
!
         !ETOT = DINTSPEC(IP,WALOC)
         ETOT = ZERO
         do id = 1, NUMDIR
           tmp(:) = WALOC(:,id) * spsig 
           ETOT = ETOT + tmp(1) * ONEHALF * ds_incr(1)*ddir
           do is = 2, NUMSIG
             ETOT = ETOT + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
           end do
           ETOT = ETOT + ONEHALF * tmp(NUMSIG) * ds_incr(NUMSIG)*ddir
         end do
!
! if etot too small skip ...
!
         if (etot .gt. thr) then
!
! integrals ... inlined ... for speed ...
!
           ACTOT = ZERO
           y = ONE/SPSIG
           do id = 1, NUMDIR
             tmp(:) = WALOC(:,id) * spsig * y
             ACTOT  = ACTOT + tmp(1) * ONEHALF * ds_incr(1)*ddir
             do is = 2, NUMSIG
               ACTOT = ACTOT + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
             end do
             ACTOT  = ACTOT + ONEHALF * tmp(NUMSIG) * ds_incr(NUMSIG)*ddir
           end do
           !tmp = ONE/SPSIG
           !ACTOT       = DINTSPEC_Y(IP,WALOC,tmp)
           ETOT_SPSIG = ZERO
           y = SIGPOW(:,1) 
           do id = 1, NUMDIR
             tmp(:) = WALOC(:,id) * spsig * y
             ETOT_SPSIG = ETOT_SPSIG + tmp(1) * ONEHALF * ds_incr(1)*ddir
             do is = 2, NUMSIG
               ETOT_SPSIG = ETOT_SPSIG + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
             end do
             ETOT_SPSIG = ETOT_SPSIG + ONEHALF * tmp(NUMSIG) * ds_incr(NUMSIG)*ddir
           end do
           !tmp = SIGPOW(:,1)
           !ETOT_SPSIG  = DINTSPEC_Y(IP,WALOC,tmp)
           ETOT_WK = ZERO
           y = WK(:,IP) 
           do id = 1, NUMDIR
             tmp(:) = WALOC(:,id) * spsig * y
             ETOT_WK = ETOT_WK + ONEHALF * tmp(1) * ds_incr(1)*ddir
             do is = 2, NUMSIG
               ETOT_WK = ETOT_WK + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
             end do
             ETOT_WK = ETOT_WK + ONEHALF * tmp(NUMSIG) * ds_incr(NUMSIG)*ddir
           end do
           !tmp = WK(:,IP)
           !ETOT_WK     = DINTSPEC_Y(IP,WALOC,tmp)
           ETOT_ISQ_WK = ZERO 
           y = ONE/SQRT(WK(:,IP))
           do id = 1, NUMDIR
             tmp(:) = WALOC(:,id) * spsig * y
             ETOT_ISQ_WK = ETOT_ISQ_WK + ONEHALF * tmp(1) * ds_incr(1)*ddir
             do is = 2, NUMSIG
               ETOT_ISQ_WK = ETOT_ISQ_WK+ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
             end do
             ETOT_ISQ_WK = ETOT_ISQ_WK+ONEHALF*tmp(NUMSIG) * ds_incr(NUMSIG)*ddir
           end do
           !tmp = ONE/SQRT(WK(:,IP))
           !ETOT_ISQ_WK = DINTSPEC_Y(IP,WALOC,tmp)
           ETOT_SQ_WK = ZERO
           y = SQRT(WK(:,IP)) 
           do id = 1, NUMDIR
             tmp(:) = WALOC(:,id) * spsig * y
             ETOT_SQ_WK = ETOT_SQ_WK + ONEHALF * tmp(1) * ds_incr(1)*ddir
             do is = 2, NUMSIG -1 
               ETOT_SQ_WK = ETOT_SQ_WK + ONEHALF*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
             end do
             ETOT_SQ_WK = ETOT_SQ_WK + ONEHALF*tmp(NUMSIG) * ds_incr(NUMSIG)*ddir
           end do
           !tmp = SQRT(WK(:,IP))
           !ETOT_SQ_WK  = DINTSPEC_Y(IP,WALOC,tmp) 
!
! tail factors ...
!
           DS          = SPSIG(NUMSIG) - SPSIG(NUMSIG-1)

           ATAIL       = SUM(WALOC(NUMSIG,:)) * SIGPOW(NUMSIG,1) * DDIR * DS
           ETAIL       = SUM(WALOC(NUMSIG,:)) * SIGPOW(NUMSIG,2) * DDIR * DS
           ESIGTAIL    = SUM(WALOC(NUMSIG,:)) * SIGPOW(NUMSIG,3) * DDIR * DS
!
! tail factors ...
!
           ACTOT       = ACTOT        + TAIL_ARR(5)  * ATAIL 
           ETOT        = ETOT         + TAIL_ARR(6)  * ETAIL 
           ETOT_SPSIG  = ETOT_SPSIG   + TAIL_ARR(7)  * ETAIL 
           ETOT_ISQ_WK = ETOT_ISQ_WK  + TAIL_ARR(5)  * ETAIL / (SQRT(WK(NUMSIG,IP)))
           ETOT_SQ_WK  = ETOT_SQ_WK   + TAIL_ARR(5)  * ETAIL * (SQRT(WK(NUMSIG,IP)))
           ETOT_WK     = ETOT_WK      + TAIL_ARR(8)  * ETAIL * WK(NUMSIG,IP)
!
! integral parameters ...
!
           HS          = MAX(ZERO,4.*SQRT(ETOT))
           IF (LMONO_OUT) HS = HS / SQRT(2.)
           SME01       = ETOT_SPSIG / ETOT
           KME01       = ETOT_WK / ETOT
           SME10       = ETOT / ACTOT
           KMWAM       = (ETOT/ETOT_ISQ_WK)**2
           KMWAM2      = (ETOT_SQ_WK/ETOT)**2

         else
!
! no or too less energy ...
!
           ETOT        = ZERO
           HS          = ZERO 
           SME01       = ZERO 
           KME01       = 10.0_rkind
           SME10       = ZERO 
           KMWAM       = 10.0_rkind
           KMWAM2      = 10.0_rkind

         end if 
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
!AR: This must be replaced by ST4_PRE or the certain WAM routine we are not consistent here this routine needs urgen deletion
      SUBROUTINE MEAN_PARAMETER_BDCONS(WALOC,HS,TM01,TM02)

      USE DATAPOOL
      IMPLICIT NONE

      REAL(rkind), INTENT(IN)    :: WALOC(NUMSIG,NUMDIR)
      REAL(rkind), INTENT(OUT)   :: HS,TM01,TM02

      INTEGER             :: ID, IS

      REAL(rkind)                :: Y(NUMSIG)
      REAL(rkind)                :: DS, ETAIL
      REAL(rkind)                :: OMEG2, EAD, ETOT
      REAL(rkind)                :: EFTAIL,PTAIL_ARR,EFTOT,ETAIL_ARR
      REAL(rkind)                :: EHFR,AHFR,ATAIL_ARR,EPTOT,APTOT
      REAL(rkind)                :: tmp(NUMSIG),actmp(NUMSIG)
!
! total energy ...
!
      Y = ZERO
      tmp = ZERO
      actmp = ZERO
      ETOT = ZERO
      do id = 1, NUMDIR
        tmp(:) = WALOC(:,id) * spsig
        ETOT = ETOT + tmp(1) * ONEHALF * ds_incr(1)*ddir
        do is = 2, NUMSIG
          ETOT = ETOT + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
        end do
        ETOT = ETOT + ONEHALF * tmp(NUMSIG) * ds_incr(NUMSIG)*ddir
      end do

      IF (ETOT .GT. THR) THEN
!
! tail ratios
!
         DS    = SPSIG(NUMSIG) - SPSIG(NUMSIG-1)
         ETAIL = SUM(WALOC(NUMSIG,:)) * SIGPOW(NUMSIG,2) * DDIR * DS
         ETOT  = ETOT + TAIL_ARR(6) * ETAIL

         HS = 4*SQRT(ETOT)

         APTOT = ZERO
         EPTOT = ZERO
         PTAIL_ARR = TAIL_ARR(1)
         ATAIL_ARR = ONE / (PTAIL_ARR * (ONE + PTAIL_ARR * (FRINTH-ONE)))
         PTAIL_ARR = TAIL_ARR(1) - ONE
         ETAIL_ARR = ONE / (PTAIL_ARR * (ONE + PTAIL_ARR * (FRINTH-ONE)))

         DO ID = 1, NUMDIR
           DO IS = 1, NUMSIG
             APTOT = APTOT + SPSIG(IS)    * WALOC(IS,ID)
             EPTOT = EPTOT + SIGPOW(IS,2) * WALOC(IS,ID)
           ENDDO
         ENDDO

         APTOT = APTOT * FRINTF
         EPTOT = EPTOT * FRINTF

         IF (NUMSIG .GT. 3  .AND. .NOT. LSIGMAX) THEN
           DO ID = 1, NUMDIR
           AHFR  = SPSIG(NUMSIG) * WALOC(NUMSIG,ID)
           APTOT = APTOT + ATAIL_ARR * AHFR
           EHFR  = SPSIG(NUMSIG) * AHFR
           EPTOT = EPTOT + ETAIL_ARR * EHFR
           ENDDO
         ENDIF

         IF (EPTOT .GT. ZERO) THEN
            TM01 = PI2 * APTOT / EPTOT
         ELSE
            TM01 = ZERO
         END IF

         ETOT  = ZERO
         EFTOT = ZERO
         PTAIL_ARR = TAIL_ARR(1) - ONE
         ETAIL  = ONE / (PTAIL_ARR * (ONE + PTAIL_ARR * (FRINTH-ONE)))
         PTAIL_ARR = TAIL_ARR(1) - 3.
         EFTAIL = ONE / (PTAIL_ARR * (ONE + PTAIL_ARR * (FRINTH-ONE)))
         DO ID=1, NUMDIR
            DO IS = 1, NUMSIG
              EAD  = SIGPOW(IS,2) * WALOC(IS,ID) * FRINTF
              OMEG2 = SIGPOW(IS,2)
              ETOT  = ETOT + EAD
              EFTOT = EFTOT + EAD * OMEG2
            ENDDO
            IF (NUMSIG .GT. 3  .AND. .NOT. LSIGMAX) THEN
              EAD  = SIGPOW(NUMSIG,2) * WALOC(NUMSIG,ID)
              ETOT  = ETOT  + ETAIL * EAD
              EFTOT = EFTOT + EFTAIL * OMEG2 * EAD
            ENDIF
         ENDDO
         IF (EFTOT .GT. sqrt(verysmall)) THEN
           TM02 = PI2 * SQRT(ETOT/EFTOT)
         ELSE
           TM02 = ZERO
         END IF

      ELSE
        HS = ZERO
        TM01 = ZERO
        TM02 = ZERO
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
!AR: This must be replaced by ST4_PRE or the certain WAM routine we are not consistent here this routine needs urgen deletion
      SUBROUTINE MEAN_PARAMETER(IP,WALOC,ISMAX,HS,TM01,TM02,TM10,KLM,WLM)

      USE DATAPOOL
      IMPLICIT NONE
!2do ... rewrite this integration ...
      INTEGER, INTENT(IN) :: IP,ISMAX

      REAL(rkind), INTENT(IN)    :: WALOC(NUMSIG,NUMDIR)
      REAL(rkind), INTENT(OUT)   :: HS,TM01,TM02,KLM,WLM,TM10

      INTEGER             :: ID, IS

      REAL(rkind)                :: Y(NUMSIG)
      REAL(rkind)                :: DS, ETAIL
      REAL(rkind)                :: OMEG2,OMEG,EAD,UXD, ETOT
      REAL(rkind)                :: EFTAIL,PTAIL_ARR,EFTOT,ETAIL_ARR
      REAL(rkind)                :: EHFR,AHFR,ATAIL_ARR,EPTOT,APTOT
      REAL(rkind)                :: CKTAIL, ETOT1, SIG22, EKTOT, CETAIL
      REAL(rkind)                :: tmp(NUMSIG),actmp(NUMSIG), SKK
!
! total energy ...
!
      Y = ZERO
      tmp = ZERO
      actmp = ZERO
      ETOT = ZERO
      do id = 1, NUMDIR
        tmp(:) = WALOC(:,id) * spsig
        ETOT = ETOT + tmp(1) * ONEHALF * ds_incr(1)*ddir
        do is = 2, NUMSIG
          ETOT = ETOT + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
        end do
        ETOT = ETOT + ONEHALF * tmp(NUMSIG) * ds_incr(NUMSIG)*ddir
      end do

      IF (ETOT .GT. THR) THEN
!
! tail ratios
!
         DS    = SPSIG(NUMSIG) - SPSIG(NUMSIG-1)
         ETAIL = SUM(WALOC(NUMSIG,:)) * SIGPOW(NUMSIG,2) * DDIR * DS
         ETOT  = ETOT + TAIL_ARR(6) * ETAIL

         HS = 4*SQRT(ETOT)

         APTOT = ZERO
         EPTOT = ZERO
         PTAIL_ARR = TAIL_ARR(1)
         ATAIL_ARR = ONE / (PTAIL_ARR * (ONE + PTAIL_ARR * (FRINTH-ONE)))
         PTAIL_ARR = TAIL_ARR(1) - ONE
         ETAIL_ARR = ONE / (PTAIL_ARR * (ONE + PTAIL_ARR * (FRINTH-ONE)))
         DO ID = 1, NUMDIR
           DO IS = 1, ISMAX
             APTOT = APTOT + SPSIG(IS)    * WALOC(IS,ID)
             EPTOT = EPTOT + SIGPOW(IS,2) * WALOC(IS,ID)
           ENDDO
         ENDDO

         APTOT = APTOT * FRINTF
         EPTOT = EPTOT * FRINTF

         IF (NUMSIG .GT. 3  .AND. .NOT. LSIGMAX) THEN
           DO ID = 1, NUMDIR
             AHFR  = SPSIG(NUMSIG) * WALOC(NUMSIG,ID)
             APTOT = APTOT + ATAIL_ARR * AHFR
             EHFR  = SPSIG(NUMSIG) * AHFR
             EPTOT = EPTOT + ETAIL_ARR * EHFR
           ENDDO
         ENDIF

         IF (EPTOT .GT. ZERO) THEN
            TM01 = PI2 * APTOT / EPTOT
         ELSE
            TM01 = ZERO
         END IF

         ETOT  = ZERO
         EFTOT = ZERO
         PTAIL_ARR = TAIL_ARR(1) - ONE
         ETAIL  = ONE / (PTAIL_ARR * (ONE + PTAIL_ARR * (FRINTH-ONE)))
         PTAIL_ARR = TAIL_ARR(1) - 3.
         EFTAIL = ONE / (PTAIL_ARR * (ONE + PTAIL_ARR * (FRINTH-ONE)))
         DO ID=1, NUMDIR
            IF (LSECU .OR. LSTCU) THEN
              UXD  = CURTXY(IP,1)*COSTH(ID) + CURTXY(IP,2)*SINTH(ID)
            ENDIF
            DO IS = 1, ISMAX
              EAD  = SIGPOW(IS,2) * WALOC(IS,ID) * FRINTF
              IF (LSECU .OR. LSTCU) THEN
                OMEG  = SPSIG(IS) + WK(IS,IP) * UXD
                OMEG2 = OMEG**2
              ELSE
                OMEG2 = SIGPOW(IS,2)
              ENDIF
              ETOT  = ETOT + EAD
              EFTOT = EFTOT + EAD * OMEG2
            ENDDO
            IF (NUMSIG .GT. 3  .AND. .NOT. LSIGMAX) THEN
              EAD  = SIGPOW(NUMSIG,2) * WALOC(NUMSIG,ID)
              ETOT  = ETOT  + ETAIL * EAD
              EFTOT = EFTOT + EFTAIL * OMEG2 * EAD
            ENDIF
         ENDDO
         IF (EFTOT .GT. sqrt(verysmall)) THEN
           TM02 = PI2 * SQRT(ETOT/EFTOT)
         ELSE
           TM02 = ZERO
         END IF

         ETOT1 = ZERO
         EKTOT = ZERO
!
! tail ratios same
!
         PTAIL_ARR = TAIL_ARR(1) - ONE
         CETAIL = ONE / (PTAIL_ARR * (ONE + PTAIL_ARR * (FRINTH-ONE)))
         PTAIL_ARR = TAIL_ARR(1) - ONE - 2*ONE
         CKTAIL = ONE / (PTAIL_ARR * (ONE + PTAIL_ARR * (FRINTH-ONE)))

         DO IS = 1, ISMAX
           SIG22 = SIGPOW(IS,2)
           SKK  = SIG22 * WK(IS,IP)
           DO ID = 1, NUMDIR
             ETOT1 = ETOT1 + SIG22 * WALOC(IS,ID)
             EKTOT = EKTOT + SKK * WALOC(IS,ID)
           ENDDO
         ENDDO
         ETOT1 = FRINTF * ETOT1
         EKTOT = FRINTF * EKTOT

         IF (NUMSIG .GT. 3) THEN
            DO ID=1,NUMDIR
              ETOT1 = ETOT1 + CETAIL * SIG22 * WALOC(NUMSIG,ID)
              EKTOT = EKTOT + CKTAIL * SKK * WALOC(NUMSIG,ID)
            ENDDO
         ENDIF

         IF (ETOT1.GT.VERYSMALL.AND.EKTOT.GT.VERYSMALL) THEN
            WLM = PI2 * (ETOT1/EKTOT)
            KLM = PI2/WLM
         ELSE
            KLM = 10.0_rkind
            WLM = ZERO
         ENDIF

         APTOT = 0.
         EPTOT = 0.
         DO ID=1, NUMDIR
           DO IS=1,ISMAX
             APTOT = APTOT + SPSIG(IS) * WALOC(IS,ID)
             EPTOT = EPTOT + SIGPOW(IS,2) * WALOC(IS,ID)
           ENDDO
         ENDDO
         APTOT = APTOT * FRINTF
         EPTOT = EPTOT * FRINTF
         IF (NUMSIG .GT. 3) THEN
           PTAIL_ARR = TAIL_ARR(1)
           ATAIL_ARR = 1. / (PTAIL_ARR * (1. + PTAIL_ARR * (FRINTH-1.)))
           PTAIL_ARR = TAIL_ARR(1) - 1.
           ETAIL_ARR = 1. / (PTAIL_ARR * (1. + PTAIL_ARR * (FRINTH-1.)))
           DO ID = 1, NUMDIR
             AHFR = SPSIG(NUMSIG) * WALOC(NUMSIG,ID)
             APTOT = APTOT + ATAIL_ARR * AHFR
             EHFR = SPSIG(NUMSIG) * AHFR
             EPTOT = EPTOT + ETAIL_ARR * EHFR
           ENDDO
         ENDIF
         TM10 = 2.*PI * APTOT / EPTOT

      ELSE

           HS = ZERO
           TM01 = ZERO
           TM02 = ZERO
           TM10 = ZERO
           KLM  = 10.0_rkind
           WLM  = ZERO

      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MEAN_PARAMETER_OUTPUT(IP,WALOC,HS,TM01,TM02,TM10,KLM,WLM)
      USE DATAPOOL
      USE W3SRC4MD
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: IP
      REAL(rkind), INTENT(IN)    :: WALOC(NUMSIG,NUMDIR)
      REAL(rkind), INTENT(OUT)   :: HS,TM01,TM02,KLM,WLM,TM10
      INTEGER ISMAX, ID, IS
      REAL(rkind) AWW3(NSPEC), FL3(NUMDIR,NUMSIG)
      REAL(rkind) WIND10, WINDTH, JAC
      REAL(rkind) FMEAN1, WNMEAN, AMAX
      REAL(rkind) F1MEAN, AKMEAN, XKMEAN
      ISMAX = NUMSIG
      CALL MEAN_PARAMETER(IP,WALOC,ISMAX,HS,TM01,TM02,TM10,KLM,WLM)
      IF (DEP(IP) .gt. 0) THEN
        IF (ISOURCE .eq. 1) THEN
          DO IS = 1, NUMSIG
            DO ID = 1, NUMDIR
              AWW3(ID + (IS-1) * NUMDIR) = WALOC(IS,ID) * CG(IS,IP)
            END DO
          END DO
          LLWS = .TRUE.
          CALL SET_WIND( IP, WIND10, WINDTH )
          CALL W3SPR4 ( AWW3, CG(:,IP), WK(:,IP), EMEAN(IP), FMEAN(IP), FMEAN1, WNMEAN, AMAX, WIND10, WINDTH, UFRIC(IP), USTDIR(IP), TAUWX(IP), TAUWY(IP), CD(IP), Z0(IP), ALPHA_CH(IP), LLWS, FMEANWS(IP))
          HS = 4.0_rkind * SQRT(EMEAN(IP))
        END IF
        IF (ISOURCE .eq. 2) THEN
          DO IS = 1, NUMSIG
            JAC = PI2 * SPSIG(IS)
            DO ID = 1, NUMDIR
              FL3(ID,IS) = WALOC(IS,ID) * JAC
            END DO
          END DO
          CALL FKMEAN_LOCAL(IP, FL3, EMEAN(IP), FMEAN(IP), F1MEAN, AKMEAN, XKMEAN)
          HS = 4.0_rkind * SQRT(EMEAN(IP))
        END IF
      ELSE
        HS = ZERO
      END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
!AR: This must be replaced by ST4_PRE or the certain WAM routine we are not consistent here this routine needs urgen deletion
      SUBROUTINE MEAN_WAVE_PARAMETER_LOC(WALOC,CURTXYLOC,DEPLOC,WKLOC,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2)

         USE DATAPOOL
         IMPLICIT NONE
 
         REAL(rkind), INTENT(IN)    :: WALOC(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(IN)    :: WKLOC(NUMSIG), DEPLOC
         REAL(rkind), INTENT(IN)    :: CURTXYLOC(2)
         REAL(rkind), INTENT(OUT)   :: SME01, SME10, KME01
         REAL(rkind), INTENT(OUT)   :: KMWAM, KMWAM2
         REAL(rkind), INTENT(OUT)   :: HS

         INTEGER             :: ID, IS

         REAL(rkind)                :: ACTOT, ETOT
         REAL(rkind)                :: ETOT_SPSIG
         REAL(rkind)                :: ETOT_WK
         REAL(rkind)                :: ETOT_ISQ_WK
         REAL(rkind)                :: ETOT_SQ_WK

         REAL(rkind)                :: Y(NUMSIG)
         REAL(rkind)                :: DS, ATAIL, ETAIL, ESIGTAIL
         REAL(rkind)                :: tmp(NUMSIG)
!
! total energy ...
!
         ETOT = ZERO
         do id = 1, NUMDIR
           tmp(:) = WALOC(:,id) * spsig
           ETOT = ETOT + tmp(1) * ONEHALF * ds_incr(1)*ddir
           do is = 2, NUMSIG
             ETOT = ETOT + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
           end do
           ETOT = ETOT + ONEHALF * tmp(NUMSIG) * ds_incr(NUMSIG)*ddir
         end do
!
! if etot too small skip ...
!
         if (etot .gt. verysmall) then

           ACTOT = ZERO
           y = ONE/SPSIG
           do id = 1, NUMDIR
             tmp(:) = WALOC(:,id) * spsig * y
             ACTOT  = ACTOT + tmp(1) * ONEHALF * ds_incr(1)*ddir
             do is = 2, NUMSIG
               ACTOT = ACTOT + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
             end do
             ACTOT  = ACTOT + ONEHALF * tmp(NUMSIG) * ds_incr(NUMSIG)*ddir
           end do

           ETOT_SPSIG = ZERO
           y = SIGPOW(:,1)
           do id = 1, NUMDIR
             tmp(:) = WALOC(:,id) * spsig * y
             ETOT_SPSIG = ETOT_SPSIG + tmp(1) * ONEHALF * ds_incr(1)*ddir
             do is = 2, NUMSIG
               ETOT_SPSIG = ETOT_SPSIG + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
             end do
             ETOT_SPSIG = ETOT_SPSIG + ONEHALF * tmp(NUMSIG) * ds_incr(NUMSIG)*ddir
           end do

           ETOT_WK = ZERO
           y = WKLOC(:)
           do id = 1, NUMDIR
             tmp(:) = WALOC(:,id) * spsig * y
             ETOT_WK = ETOT_WK + ONEHALF * tmp(1) * ds_incr(1)*ddir
             do is = 2, NUMSIG
               ETOT_WK = ETOT_WK + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
             end do
             ETOT_WK = ETOT_WK + ONEHALF * tmp(NUMSIG) * ds_incr(NUMSIG)*ddir
           end do

           ETOT_ISQ_WK = ZERO
           y = ONE/SQRT(WKLOC)
           do id = 1, NUMDIR
             tmp(:) = WALOC(:,id) * spsig * y
             ETOT_ISQ_WK = ETOT_ISQ_WK + ONEHALF * tmp(1) * ds_incr(1)*ddir
             do is = 2, NUMSIG -1 
               ETOT_ISQ_WK = ETOT_ISQ_WK+ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
             end do
             ETOT_ISQ_WK = ETOT_ISQ_WK+ONEHALF*tmp(NUMSIG) * ds_incr(NUMSIG)*ddir
           end do

           ETOT_SQ_WK = ZERO
           y = SQRT(WKLOC)
           do id = 1, NUMDIR
             tmp(:) = WALOC(:,id) * spsig * y
             ETOT_SQ_WK = ETOT_SQ_WK + ONEHALF * tmp(1) * ds_incr(1)*ddir
             do is = 2, NUMSIG -1 
               ETOT_SQ_WK = ETOT_SQ_WK + ONEHALF*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
             end do
             ETOT_SQ_WK = ETOT_SQ_WK + ONEHALF*tmp(NUMSIG) * ds_incr(NUMSIG)*ddir
           end do
!
           DS          = SPSIG(NUMSIG) - SPSIG(NUMSIG-1)

           ATAIL       = SUM(WALOC(NUMSIG,:)) * SIGPOW(NUMSIG,1) * DDIR * DS
           ETAIL       = SUM(WALOC(NUMSIG,:)) * SIGPOW(NUMSIG,2) * DDIR * DS
           ESIGTAIL    = SUM(WALOC(NUMSIG,:)) * SIGPOW(NUMSIG,3) * DDIR * DS
!
! tail factors ...
!
           ACTOT       = ACTOT        + TAIL_ARR(5)  * ATAIL 
           ETOT        = ETOT         + TAIL_ARR(6)  * ETAIL 
           ETOT_SPSIG  = ETOT_SPSIG   + TAIL_ARR(7)  * ETAIL 
           ETOT_ISQ_WK = ETOT_ISQ_WK  + TAIL_ARR(5)  * ETAIL / SQRT(WKLOC(NUMSIG))
           ETOT_SQ_WK  = ETOT_SQ_WK   + TAIL_ARR(5)  * ETAIL * SQRT(WKLOC(NUMSIG))
           ETOT_WK     = ETOT_WK      + TAIL_ARR(8)  * ETAIL * WKLOC(NUMSIG) 
!
! integral parameters ...
!
           HS          = MAX(ZERO,4.*SQRT(ETOT))
           SME01       = ETOT_SPSIG / ETOT
           KME01       = ETOT_WK / ETOT
           SME10       = ETOT / ACTOT
           KMWAM       = (ETOT/ETOT_ISQ_WK)**2
           KMWAM2      = (ETOT_SQ_WK/ETOT)**2

         else
!
! no or too less energy ...
!
           HS          = ZERO 
           SME01       = ZERO 
           KME01       = 10.0_rkind
           SME10       = ZERO 
           KMWAM       = 10.0_rkind
           KMWAM2      = 10.0_rkind

         end if 
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
!AR: This must be replaced by ST4_PRE or the certain WAM routine we are not consistent here this routine needs urgen deletion
      SUBROUTINE WAVE_CURRENT_PARAMETER(IP,WALOC,UBOT,ORBITAL,BOTEXPER,TMBOT,CALLFROM)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP
         REAL(rkind),    INTENT(IN) :: WALOC(NUMSIG,NUMDIR)
         CHARACTER(len=*), INTENT(IN) :: CALLFROM

         REAL(rkind), INTENT(OUT)   :: UBOT, ORBITAL, BOTEXPER, TMBOT

         INTEGER             :: ID, IS

         REAL(rkind)                :: ETOT_SKD
         REAL(rkind)                :: ETOT_SKDSIG, TMP(NUMSIG), Y(NUMSIG)

!
! integrals ...
!
         IF (DEP(IP) .LT. DMIN) RETURN
         
         ETOT_SKD    = ZERO
         ETOT_SKDSIG = ZERO

         y = ONE/SINH(MIN(KDMAX,WK(:,IP)*DEP(IP)))**2

         do id = 1, NUMDIR
           tmp(:) = WALOC(:,id) * spsig * y
           ETOT_SKD  = ETOT_SKD + tmp(1) * ONEHALF * ds_incr(1)*ddir
           do is = 2, NUMSIG -1 
             ETOT_SKD = ETOT_SKD + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
           end do
           ETOT_SKD = ETOT_SKD + tmp(NUMSIG) * ONEHALF * ds_incr(NUMSIG)*ddir
         end do
 
         y =  SIGPOW(:,2)*ONE/SINH(MIN(KDMAX,WK(:,IP)*DEP(IP)))**2

         do id = 1, NUMDIR
           tmp(:) = WALOC(:,id) * spsig * y
           ETOT_SKDSIG = ETOT_SKDSIG + tmp(1) * ONEHALF * ds_incr(1)*ddir
           do is = 2, NUMSIG -1 
             ETOT_SKDSIG = ETOT_SKDSIG + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
           end do 
           ETOT_SKDSIG = ETOT_SKDSIG + tmp(NUMSIG) * ONEHALF * ds_incr(NUMSIG)*ddir
         end do
 
         IF (ETOT_SKD .gt. verysmall) THEN 
!
! integral parameters ...
!
           UBOT        = SQRT(ETOT_SKD)
           ORBITAL     = SQRT(2*ETOT_SKD)
           BOTEXPER    = SQRT(2*ETOT_SKDSIG)
           TMBOT       = PI2*SQRT(ETOT_SKDSIG/ETOT_SKD)
           
         ELSE 
!
! no or too less energy ...
! 
           UBOT        = ZERO 
           ORBITAL     = ZERO 
           BOTEXPER    = ZERO 
           TMBOT       = ZERO 

         ENDIF 
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
!AR: This must be replaced by ST4_PRE or the certain WAM routine we are not consistent here this routine needs urgen deletion
      SUBROUTINE WAVE_CURRENT_PARAMETER_LOC(WALOC,CURTXYLOC,DEPLOC,WKLOC,UBOT,ORBITAL,BOTEXPER,TMBOT)
         USE DATAPOOL
         IMPLICIT NONE

         REAL(rkind), INTENT(IN)    :: WALOC(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(IN)    :: WKLOC(NUMSIG), DEPLOC
         REAL(rkind), INTENT(IN)    :: CURTXYLOC(2)

         REAL(rkind), INTENT(OUT)   :: UBOT, ORBITAL, BOTEXPER, TMBOT

         INTEGER             :: ID, IS

         REAL(rkind)                :: ETOT_SKD
         REAL(rkind)                :: ETOT_SKDSIG, TMP(NUMSIG), Y(NUMSIG)
!
! integrals ...
!
         IF (DEPLOC .LT. DMIN) RETURN

         ETOT_SKD    = ZERO
         ETOT_SKDSIG = ZERO

         y = ONE/SINH(MIN(KDMAX,WKLOC*DEPLOC))**2
         do id = 1, NUMDIR
           tmp(:) = WALOC(:,id) * spsig * y
           ETOT_SKD  = ETOT_SKD + tmp(1) * ONEHALF * ds_incr(1)*ddir
           do is = 2, NUMSIG -1
             ETOT_SKD = ETOT_SKD + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
           end do
           ETOT_SKD = ETOT_SKD + tmp(NUMSIG) * ONEHALF * ds_incr(NUMSIG)*ddir
         end do

         y =  SIGPOW(:,2)*ONE/SINH(MIN(KDMAX,WKLOC*DEPLOC))**2
         do id = 1, NUMDIR
           tmp(:) = WALOC(:,id) * spsig * y
           ETOT_SKDSIG = ETOT_SKDSIG + tmp(1) * ONEHALF * ds_incr(1)*ddir
           do is = 2, NUMSIG -1
             ETOT_SKDSIG = ETOT_SKDSIG + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
           end do
           ETOT_SKDSIG = tmp(NUMSIG) * ONEHALF * ds_incr(NUMSIG)*ddir
         end do

         IF (ETOT_SKD .gt. verysmall) THEN 
!
! integral parameters ...
!
           UBOT        = SQRT(ETOT_SKD)
           ORBITAL     = SQRT(2*ETOT_SKD)
           BOTEXPER    = SQRT(2*ETOT_SKDSIG)
           TMBOT       = PI2*SQRT(ETOT_SKDSIG/ETOT_SKD)
           
         ELSE 
!
! no or too less energy ...
! 
           UBOT        = ZERO 
           ORBITAL     = ZERO 
           BOTEXPER    = ZERO 
           TMBOT       = ZERO 

         ENDIF 

         !WRITE(*,'(9F15.4)') DEPLOC,SUM(WALOC),CURTXYLOC,SUM(WKLOC), ETOT_SKD, ETOT_SKDSIG 
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
!AR: This must be replaced by ST4_PRE or the certain WAM routine we are not consistent here this routine needs urgen deletion
      SUBROUTINE URSELL_NUMBER(HS,SME,DEPTH,URSELL)
         USE DATAPOOL, ONLY : G9, DMIN, verysmall, rkind, ONE, TWO, ZERO

         IMPLICIT NONE

         REAL(rkind), INTENT(IN)    :: HS, SME, DEPTH
         REAL(rkind), INTENT(OUT)   :: URSELL 

         IF (DEPTH .GT. DMIN .AND. SME .GT. verysmall .AND. HS .GT. verysmall) THEN
           URSELL = (G9 * HS)/(TWO*SQRT(TWO)*SME**2*DEPTH**2)
         ELSE
           URSELL = ZERO
         END IF 
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
!AR: This must be replaced by ST4_PRE or the certain WAM routine we are not consistent here this routine needs urgen deletion
      SUBROUTINE WINDSEASWELLSEP( IP, WALOC, TM_W, CGP_W, CP_W, TP_W, LP_W, HS_W, KP_W )
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)           :: IP
         REAL(rkind)   , INTENT(IN)    :: WALOC(NUMSIG,NUMDIR)
         REAL(rkind)   , INTENT(OUT)   :: TM_W, CGP_W, CP_W, TP_W, LP_W, HS_W, KP_W

         INTEGER                :: ID, IS

         REAL(rkind)                   :: ETOT
         REAL(rkind)                   :: VEC2RAD, EFTOT, OMEG, WINDTH
         REAL(rkind)                   :: EAD, DS, EHFR, EFTAIL, ETAIL, PTAIL_ARR, ACWIND(NUMSIG,NUMDIR)
         REAL(rkind)                   :: UXD, ETOTF3, ETOTF4, FP, WN_W, WVC, WKDEP_W

         WINDTH = VEC2RAD(WINDXY(IP,1),WINDXY(IP,2))

         ETOT = ZERO
         EFTAIL = ONE / (TAIL_ARR(1)-ONE)

         DO ID  = 1, NUMDIR            ! Calculate wind sea energy ... weak criterion
           DO IS = 1, NUMSIG
             WVC = MyREAL(SPSIG(IS)/WK(IS,IP))
             IF (  1.2*UFRIC(IP)*COS(SPDIR(ID)-WINDTH)*(28./WVC) .LT. ONE) THEN
               ACWIND(IS,ID) = ZERO   ! Swell
             ELSE
               ACWIND(IS,ID) = WALOC(IS,ID)  ! Wind Sea
             END IF
           END DO
           DO IS = 2, NUMSIG
              DS = SPSIG(IS) - SPSIG(IS-1)
              EAD = ONEHALF*(SPSIG(IS)*ACWIND(IS,ID)+SPSIG(IS-1)*ACWIND(IS-1,ID))*DS*DDIR
              ETOT = ETOT + EAD
           END DO
           IF (NUMSIG > 3) THEN
             EHFR = WALOC(NUMSIG,ID) * SPSIG(NUMSIG)
             ETOT = ETOT + DDIR * EHFR * SPSIG(NUMSIG) * EFTAIL
           ENDIF
         END DO

         IF (ETOT .LT. VERYSMALL) RETURN

         IF (ETOT > ZERO) THEN
           HS_W = 4.0*SQRT(ETOT)
         ELSE
           HS_W = ZERO
         END IF

         ETOTF3 = ZERO
         ETOTF4 = ZERO
         DO IS = 1, NUMSIG
           DO ID = 1, NUMDIR
             ETOTF3 = ETOTF3 + SPSIG(IS) * WALOC(IS,ID)**4 * DDIR * DS_BAND(IS)
             ETOTF4 = ETOTF4 +             WALOC(IS,ID)**4 * DDIR * DS_BAND(IS)
           END DO
         END DO
         IF(ETOTF4 .GT. VERYSMALL) THEN
            FP = ETOTF3/ETOTF4*PI2
            TP_W = ONE/FP/PI2
            CALL WAVEKCG(DEP(IP), FP, WN_W, CP_W, KP_W, CGP_W)
            !CALL ALL_FROM_TABLE(FP,DEP(IP),KP_W,CGP_W,WKDEP_W,WN_W,CP_W)
            LP_W  = PI2/KP_w
         ELSE
            FP    = ZERO 
            CP_W  = ZERO 
            KP_W  = 10.0_rkind
            CGP_W = ZERO 
            LP_W  = ZERO 
         END IF

         ETOT = ZERO
         EFTOT = ZERO
         PTAIL_ARR = TAIL_ARR(1) - ONE
         ETAIL  = ONE / (PTAIL_ARR * (ONE + PTAIL_ARR * (FRINTH-ONE)))
         PTAIL_ARR = TAIL_ARR(1) - 2.
         EFTAIL = ONE / (PTAIL_ARR * (ONE + PTAIL_ARR * (FRINTH-ONE)))
         DO ID = 1, NUMDIR
            UXD = CURTXY(IP,1)*COSTH(ID) + CURTXY(IP,2)*SINTH(ID)
            DO IS = 1, NUMSIG
              OMEG = SPSIG(IS) + WK(IS,IP) * UXD
              EAD = FRINTF * SIGPOW(IS,2) * ACWIND(IS,ID)
              ETOT = ETOT + EAD
              EFTOT = EFTOT + EAD * OMEG
            ENDDO
            IF (NUMSIG .GT. 3) THEN
              EAD = SIGPOW(NUMSIG,2) * ACWIND(NUMSIG,ID)
              ETOT = ETOT + ETAIL * EAD
             EFTOT = EFTOT + EFTAIL * OMEG * EAD
            ENDIF
         ENDDO
         IF (EFTOT.GT.ZERO) THEN
            TM_W = PI2 * ETOT / EFTOT
         ELSE
            TM_W = ZERO
         ENDIF
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
!AR: This can remain ...
      SUBROUTINE PEAK_PARAMETER(IP,WALOC,ISMAX,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)           :: IP, ISMAX
         REAL(rkind), INTENT(IN)       :: WALOC(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(OUT)      :: KPP, FPP, CPP, WNPP, CGPP, TPP, LPP, PEAKDSPR,PEAKDM,DPEAK
         REAL(rkind), INTENT(OUT)      :: TPPD,KPPD,CGPD,CPPD

         INTEGER                       :: IS, ID, IDIRM, ISIGMP
         REAL(rkind)                   :: HQUOT, HQUOTP, ETOTF3, ETOTF4, ETOTC4, ETOTS4, PEAKFF,WKDEPD,WNPD
         REAL(rkind)                   :: DEG, VEC2DEG, MAXAC, E1, E2, DS, EAD, ETOTT, CPWN
!
! Peak period continues version... Taken from Thesis Henriques Alves ... correct citation is given there ... :)
!
       MAXAC = MAXVAL(WALOC)

       IF (MAXAC .gt. VERYSMALL .AND.  DEP(IP) .GT. DMIN) THEN

         ETOTF3 = ZERO
         ETOTF4 = ZERO
         ETOTC4 = ZERO
         ETOTS4 = ZERO

         DO IS = 1, NUMSIG
           DO ID = 1, NUMDIR
             HQUOT  = WALOC(IS,ID)/MAXAC
             HQUOTP = HQUOT**4
             ETOTF3 = ETOTF3 + SPSIG(IS) * HQUOTP * DS_BAND(IS)
             ETOTF4 = ETOTF4 +             HQUOTP * DS_BAND(IS)
             ETOTC4 = ETOTC4 + COSTH(ID) * HQUOTP * DS_BAND(IS)
             ETOTS4 = ETOTS4 + SINTH(ID) * HQUOTP * DS_BAND(IS)
           END DO
         END DO

         IF(ETOTF4 .GT. VERYSMALL .AND. ETOTF4 .GT. VERYSMALL) THEN
           FPP    = ETOTF3/ETOTF4*PI2
           CALL WAVEKCG(DEP(IP), FPP, WNPP, CPP, KPP, CGPP)
           TPP    = PI2/FPP
           LPP    = PI2/KPP
           PEAKDM = VEC2DEG (ETOTC4, ETOTS4)
           CALL DEG2NAUT(PEAKDM,DEG,LNAUTOUT)
           PEAKDM = DEG
           PEAKFF = MIN(ONE,SQRT(ETOTC4*ETOTC4+ETOTS4*ETOTS4)/ETOTF4)
           PEAKDSPR = SQRT(MAX(ZERO,2.-2.*PEAKFF)) * 180./PI

         ELSE 

           FPP = ZERO
           KPP = 10.0_rkind
           CGPP = ZERO 
           WNPP = ZERO
           CPP = ZERO
           TPP = ZERO
           LPP = ZERO
           PEAKDM = ZERO
           PEAKDSPR = ZERO 

         END IF

         DPEAK = 1
         ETOTT = ZERO
         IDIRM = -1
         DO ID = 1, NUMDIR
            EAD = ZERO
            DO IS = 2, ISMAX
               DS = SPSIG(IS)-SPSIG(IS-1)
               E1 = SPSIG(IS-1)*WALOC(IS-1,ID)
               E2 = SPSIG(IS)*WALOC(IS,ID)
               EAD = EAD + DS*(E1+E2)
            END DO
            IF (EAD .GT. ETOTT) THEN
               ETOTT = EAD
               IDIRM = ID
            END IF
         END DO
         IF (IDIRM .GT. 0) THEN
           DPEAK    = SPDIR(IDIRM)
         ELSE
           DPEAK    = ZERO
         END IF

       ELSE

         FPP = ZERO
         KPP = 10.0_rkind
         CGPP = ZERO
         WNPP = ZERO
         CPP = ZERO
         TPP = ZERO
         LPP = ZERO
         PEAKDM = ZERO
         PEAKDSPR = ZERO
         DPEAK = ZERO

       END IF


       ETOTT = ZERO
       ISIGMP = -1
       DO IS = 1, NUMSIG
         EAD = ZERO
         DO ID = 1, NUMDIR
            EAD = EAD + SPSIG(IS)*WALOC(IS,ID)*DDIR
         ENDDO
         IF (EAD > ETOTT) THEN
           ETOTT = EAD
           ISIGMP = IS
          END IF
       END DO
       IF (ISIGMP > 0) THEN
          TPPD = ONE/(SPSIG(ISIGMP)/PI2)
          CALL WAVEKCG(DEP(IP), SPSIG(ISIGMP), CPWN, CPPD, KPPD, CGPD)
       ELSE
          TPPD = ZERO
          CPPD  = ZERO
          KPPD  = ZERO
          CGPD  = ZERO
       END IF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PEAK_PARAMETER_LOC(WALOC,DEPLOC,ISMAX,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)    :: ISMAX
         REAL(rkind), INTENT(IN)    :: WALOC(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(IN)    :: DEPLOC


         REAL(rkind)  , INTENT(OUT)    :: KPP,CPP,WNPP,CGPP,TPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD,FPP

         INTEGER                :: IS, ID, IDIRM, ISIGMP
         REAL(rkind)            :: HQUOT, HQUOTP, ETOTF3, ETOTF4, ETOTC4, ETOTS4, PEAKFF
         REAL(rkind)            :: DEG, VEC2DEG, MAXAC, E1, E2, DS, EAD, ETOTT,WKDEPD,WNPD
!
! Peak period continues version... Taken from Thesis Henriques Alves ... correct citation is given there ... :)
!
       MAXAC = MAXVAL(WALOC)
       IF (MAXAC .gt. VERYSMALL .AND. DEPLOC .GT. DMIN) THEN
         ETOTF3 = ZERO
         ETOTF4 = ZERO
         ETOTC4 = ZERO
         ETOTS4 = ZERO
         DO IS = 1, NUMSIG
           DO ID = 1, NUMDIR
             HQUOT  = WALOC(IS,ID)/MAXAC
             HQUOTP = HQUOT**4
             ETOTF3 = ETOTF3 + SPSIG(IS) * HQUOTP * DS_BAND(IS)
             ETOTF4 = ETOTF4 +             HQUOTP * DS_BAND(IS)
             ETOTC4 = ETOTC4 + COSTH(ID) * HQUOTP * DS_BAND(IS)
             ETOTS4 = ETOTS4 + SINTH(ID) * HQUOTP * DS_BAND(IS)
           END DO
         END DO
         IF(ETOTF4 .GT. VERYSMALL .AND. ETOTF4 .GT. VERYSMALL) THEN
           FPP    = ETOTF3/ETOTF4*PI2
           CALL WAVEKCG(DEPLOC, FPP, WNPP, CPP, KPP, CGPP)
           PEAKDM = VEC2DEG (ETOTC4, ETOTS4)
           TPP    = PI2/FPP
           LPP    = PI2/KPP
           CALL DEG2NAUT(PEAKDM,DEG,LNAUTOUT)
           PEAKDM = DEG
           PEAKFF = MIN(ONE,SQRT(MAX(ZERO,ETOTC4*ETOTC4+ETOTS4*ETOTS4))/ETOTF4)
           PEAKDSPR = SQRT(MAX(ZERO,2.-2.*PEAKFF)) * 180./PI
         ELSE
           FPP = ZERO
           KPP = 10.0_rkind
           CGPP = ZERO
           WNPP = ZERO
           CPP = ZERO
           TPP = ZERO
           LPP = ZERO
           PEAKDM = ZERO
           PEAKDSPR = ZERO
         END IF
         DPEAK = 1
         ETOTT = ZERO
         IDIRM = -1
         DO ID = 1, NUMDIR
            EAD = ZERO
            DO IS = 2, ISMAX
               DS = SPSIG(IS)-SPSIG(IS-1)
               E1 = SPSIG(IS-1)*WALOC(IS-1,ID)
               E2 = SPSIG(IS)*WALOC(IS,ID)
               EAD = EAD + DS*(E1+E2)
            END DO
            IF (EAD .GT. ETOTT) THEN
               ETOTT = EAD
               IDIRM = ID
            END IF
         END DO
         IF (IDIRM .GT. 0) THEN
           DPEAK    = SPDIR(IDIRM) * RADDEG
           CALL DEG2NAUT(DPEAK,DEG,LNAUTOUT)
           DPEAK = DEG
         ELSE
           DPEAK = ZERO
         END IF
       ELSE
         FPP = ZERO
         KPP = 10.0_rkind
         CGPP = ZERO
         WNPP = ZERO
         CPP = ZERO
         TPP = ZERO
         LPP = ZERO
         PEAKDM = ZERO
         PEAKDSPR = ZERO
         DPEAK = ZERO
       END IF

       ETOTT = ZERO
       ISIGMP = -1
       DO IS = 1, NUMSIG
         EAD = ZERO
         DO ID = 1, NUMDIR
            EAD = EAD + SPSIG(IS)*WALOC(IS,ID)*DDIR
         ENDDO
         IF (EAD > ETOTT) THEN
           ETOTT = EAD
           ISIGMP = IS
         END IF
       END DO
       IF (ISIGMP > 0) THEN
          TPPD = ONE/(SPSIG(ISIGMP)/PI2)
          !CALL WAVEKCG(DEPLOC, SPSIG(ISIGMP), CPPD, KPPD, CGPD)
          CALL ALL_FROM_TABLE(SPSIG(ISIGMP),DEPLOC,KPPD,CGPD,WKDEPD,WNPD,CPPD)
       ELSE
          TPPD = ZERO
          CPPD  = ZERO
          KPPD  = ZERO
          CGPD  = ZERO
       END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
!AR: This must be replaced by ST4_PRE or the certain WAM routine we are not consistent here this routine needs urgen deletion
      SUBROUTINE MEAN_DIRECTION_AND_SPREAD(IP,WALOC,ISMAX,ETOTS,ETOTC,DM,DSPR)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)    :: IP, ISMAX
         REAL(rkind),    INTENT(IN)    :: WALOC(NUMSIG,NUMDIR)
         REAL(rkind),    INTENT(OUT)   :: ETOTS, ETOTC, DM, DSPR

         INTEGER                       :: IS, ID
         REAL(rkind)                   :: DS, EDI, EAD, ETOT1, EHFR
         REAL(rkind)                   :: EFTAIL, VEC2DEG, DEG, FF

         ETOTC = ZERO
         ETOTS = ZERO
         ETOT1 = ZERO

         EFTAIL = ONE / (TAIL_ARR(1)-ONE)

         DO ID = 1, NUMDIR
           EAD = ZERO
             DO  IS = 2, ISMAX 
               DS  = SPSIG(IS)-SPSIG(IS-1)
               EDI = ONEHALF*(SPSIG(IS)*WALOC(IS,ID)+SPSIG(IS-1)*WALOC(IS-1,ID))*DS
               EAD = EAD + EDI
            END DO
            IF (NUMSIG .GT. 3) THEN
              EHFR = WALOC(NUMSIG,ID) * SPSIG(NUMSIG)
              EAD = EAD + EHFR * SPSIG(NUMSIG) * EFTAIL
            ENDIF
            EAD = EAD * DDIR
            ETOT1 = ETOT1 + EAD
            ETOTC  = ETOTC + EAD * COSTH(ID)
            ETOTS  = ETOTS + EAD * SINTH(ID)
         END DO

         IF (ETOT1 .GT. verysmall) THEN
           DM    = VEC2DEG (ETOTC, ETOTS)
           CALL DEG2NAUT(DM,DEG,LNAUTOUT)
           DM = DEG
           FF = MIN (ONE, SQRT(MAX(ZERO,ETOTC*ETOTC+ETOTS*ETOTS))/ETOT1)
           DSPR = SQRT(MAX(ZERO,2.-2.*FF)) * 180./PI
         ELSE
           FF   = ZERO
           DM   = ZERO
           DSPR = ZERO
         END IF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
!AR: This must be replaced by ST4_PRE or the certain WAM routine we are not consistent here this routine needs urgen deletion
      SUBROUTINE MEAN_DIRECTION_AND_SPREAD_LOC(WALOC,ISMAX,ETOTS,ETOTC,DM,DSPR)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)    :: ISMAX
         REAL(rkind),    INTENT(IN)    :: WALOC(NUMSIG,NUMDIR)
         REAL(rkind),    INTENT(OUT)   :: ETOTS, ETOTC, DM, DSPR

         INTEGER                :: IS, ID

         REAL(rkind)                   :: DS, EDI, EAD, ETOT1, EHFR
         REAL(rkind)                  :: EFTAIL, VEC2DEG, DEG, FF

         ETOTC = ZERO
         ETOTS = ZERO
         ETOT1  = ZERO

         EFTAIL = ONE / (TAIL_ARR(1)-ONE)

         DO ID = 1, NUMDIR
           EAD = ZERO
             DO  IS = 2, ISMAX 
               DS  = SPSIG(IS)-SPSIG(IS-1)
               EDI = ONEHALF*(SPSIG(IS)*WALOC(IS,ID)+SPSIG(IS-1)*WALOC(IS-1,ID))*DS
               EAD = EAD + EDI
            END DO
            IF (NUMSIG .GT. 3) THEN
              EHFR = WALOC(NUMSIG,ID) * SPSIG(NUMSIG)
              EAD = EAD + EHFR * SPSIG(NUMSIG) * EFTAIL
            ENDIF
            EAD = EAD * DDIR
            ETOT1 = ETOT1 + EAD
            ETOTC  = ETOTC + EAD * COSTH(ID)
            ETOTS  = ETOTS + EAD * SINTH(ID)
         END DO

         IF (ETOT1 .GT. verysmall ) THEN
           DM    = VEC2DEG (ETOTC, ETOTS)
           CALL DEG2NAUT(DM,DEG,LNAUTOUT)
           DM = DEG
           FF = MIN (ONE, SQRT(MAX(ZERO,ETOTC*ETOTC+ETOTS*ETOTS))/ETOT1)
           DSPR = SQRT(MAX(ZERO,2.-2.*FF)) * 180./PI
         ELSE
           FF   = ZERO
           DM   = ZERO
           DSPR = ZERO
         END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
!AR: This must be replaced by ST4_PRE or the certain WAM routine we are not consistent here this routine needs urgen deletion
      SUBROUTINE MEAN_PARAMETER_LOC(WALOC,CURTXYLOC,DEPLOC,WKLOC,ISMAX,HS,TM01,TM02,TM10,KLM,WLM)
      USE DATAPOOL
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: ISMAX
      REAL(rkind), INTENT(IN)    :: WALOC(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(IN)    :: WKLOC(NUMSIG), DEPLOC
         REAL(rkind), INTENT(IN)    :: CURTXYLOC(2)


         REAL(rkind), INTENT(OUT)   :: HS,TM01,TM02,KLM,WLM,TM10

         INTEGER             :: ID, IS

         REAL(rkind)                :: DS,ETAIL
         REAL(rkind)                :: OMEG2,OMEG,EAD,UXD, ETOT
         REAL(rkind)                :: EFTAIL,PTAIL_ARR,EFTOT,ETAIL_ARR
         REAL(rkind)                :: EHFR,AHFR,ATAIL_ARR,EPTOT,APTOT
         REAL(rkind)                :: SKK, CKTAIL, ETOT1, SIG22, EKTOT, CETAIL
         REAL(rkind)                :: tmp(NUMSIG)
!
! total energy ...
!
         ETOT = ZERO
         do id = 1, NUMDIR
           tmp(:) = WALOC(:,id) * spsig
           ETOT = ETOT + tmp(1) * ONEHALF * ds_incr(1)*ddir
           do is = 2, NUMSIG
             ETOT = ETOT + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
           end do
           ETOT = ETOT + ONEHALF * tmp(NUMSIG) * ds_incr(NUMSIG)*ddir
         end do

         IF (ETOT .GT. verysmall) THEN
!
! tail ratios same as in swan ...
!
         DS    = SPSIG(NUMSIG) - SPSIG(NUMSIG-1)
         ETAIL = SUM(WALOC(NUMSIG,:)) * SIGPOW(NUMSIG,2) * DDIR * DS
         ETOT  = ETOT + TAIL_ARR(6) * ETAIL

         HS = 4*SQRT(ETOT)

         APTOT = ZERO
         EPTOT = ZERO
!
! tail ratios same as in swan ...
!
         PTAIL_ARR = TAIL_ARR(1)
         ATAIL_ARR = ONE / (PTAIL_ARR * (ONE + PTAIL_ARR * (FRINTH-ONE)))
         PTAIL_ARR = TAIL_ARR(1) - ONE
         ETAIL_ARR = ONE / (PTAIL_ARR * (ONE + PTAIL_ARR * (FRINTH-ONE)))

         DO ID = 1, NUMDIR
           DO IS = 1, ISMAX
             APTOT = APTOT + SPSIG(IS)    * WALOC(IS,ID)
             EPTOT = EPTOT + SIGPOW(IS,2) * WALOC(IS,ID)
           ENDDO
         ENDDO

         APTOT = APTOT * FRINTF
         EPTOT = EPTOT * FRINTF

         IF (NUMSIG .GT. 3  .AND. .NOT. LSIGMAX) THEN
           DO ID = 1, NUMDIR
             AHFR  = SPSIG(NUMSIG) * WALOC(NUMSIG,ID)
             APTOT = APTOT + ATAIL_ARR * AHFR
             EHFR  = SPSIG(NUMSIG) * AHFR
             EPTOT = EPTOT + ETAIL_ARR * EHFR
           ENDDO
         ENDIF

         IF (EPTOT .GT. VERYSMALL) THEN
            TM01 = PI2 * APTOT / EPTOT
         ELSE
            TM01 = ZERO
         END IF

         ETOT  = ZERO
         EFTOT = ZERO

         PTAIL_ARR = TAIL_ARR(1) - ONE
         ETAIL  = ONE / (PTAIL_ARR * (ONE + PTAIL_ARR * (FRINTH-ONE)))
         PTAIL_ARR = TAIL_ARR(1) - 3.
         EFTAIL = ONE / (PTAIL_ARR * (ONE + PTAIL_ARR * (FRINTH-ONE)))
!
! tail ratios same as in swan ...
!
         DO ID=1, NUMDIR
            IF (LSECU .OR. LSTCU) THEN
              UXD  = CURTXYLOC(1)*COSTH(ID) + CURTXYLOC(2)*SINTH(ID)
            ENDIF
            DO IS = 1, ISMAX
              EAD  = SIGPOW(IS,2) * WALOC(IS,ID) * FRINTF
              IF (LSECU .OR. LSTCU) THEN
                OMEG  = SPSIG(IS) + WKLOC(IS) * UXD
                OMEG2 = OMEG**2
              ELSE
                OMEG2 = SIGPOW(IS,2)
              ENDIF
              ETOT  = ETOT + EAD
              EFTOT = EFTOT + EAD * OMEG2
            ENDDO
            IF (NUMSIG .GT. 3  .AND. .NOT. LSIGMAX) THEN
              EAD  = SIGPOW(NUMSIG,2) * WALOC(NUMSIG,ID)
              ETOT  = ETOT  + ETAIL * EAD
              EFTOT = EFTOT + EFTAIL * OMEG2 * EAD
            ENDIF
         ENDDO
         IF (EFTOT .GT. verysmall .AND. ETOT .GT. verysmall) THEN
           TM02 = PI2 * SQRT(ETOT/EFTOT)
         ELSE
           TM02 = ZERO
         END IF

         ETOT1 = ZERO
         EKTOT = ZERO
!
! tail ratios
!
         PTAIL_ARR = TAIL_ARR(1) - ONE
         CETAIL = ONE / (PTAIL_ARR * (ONE + PTAIL_ARR * (FRINTH-ONE)))
         PTAIL_ARR = TAIL_ARR(1) - ONE - 2.*ONE
         CKTAIL = ONE / (PTAIL_ARR * (ONE + PTAIL_ARR * (FRINTH-ONE)))

         DO IS = 1, ISMAX
           SIG22 = SIGPOW(IS,2)
           SKK  = SIG22 * WKLOC(IS)
           DO ID = 1, NUMDIR
             ETOT1 = ETOT1 + SIG22 * WALOC(IS,ID)
             EKTOT = EKTOT + SKK * WALOC(IS,ID)
           ENDDO
         ENDDO

         ETOT1 = FRINTF * ETOT1
         EKTOT = FRINTF * EKTOT

         IF (NUMSIG .GT. 3) THEN
            DO ID=1,NUMDIR
              ETOT1 = ETOT1 + CETAIL * SIG22 * WALOC(NUMSIG,ID)
              EKTOT = EKTOT + CKTAIL * SKK * WALOC(NUMSIG,ID)
            ENDDO
         ENDIF

         IF (ETOT1.GT.VERYSMALL.AND.EKTOT.GT.VERYSMALL) THEN
            WLM = PI2 * (ETOT1/EKTOT)
            KLM  = PI2/WLM
         ELSE
            KLM  = 10.0_rkind
            WLM = ZERO
         ENDIF

         APTOT = 0.
         EPTOT = 0.
         DO ID=1, NUMDIR
           DO IS=1,ISMAX
             APTOT = APTOT + SPSIG(IS) * WALOC(IS,ID)
             EPTOT = EPTOT + SIGPOW(IS,2) * WALOC(IS,ID)
           ENDDO
         ENDDO
         APTOT = APTOT * FRINTF
         EPTOT = EPTOT * FRINTF
         IF (NUMSIG .GT. 3) THEN
           PTAIL_ARR = TAIL_ARR(1)
           ATAIL_ARR = 1. / (PTAIL_ARR * (1. + PTAIL_ARR * (FRINTH-1.)))
           PTAIL_ARR = TAIL_ARR(1) - 1.
           ETAIL_ARR = 1. / (PTAIL_ARR * (1. + PTAIL_ARR * (FRINTH-1.)))
           DO ID = 1, NUMDIR
             AHFR = SPSIG(NUMSIG) * WALOC(NUMSIG,ID)
             APTOT = APTOT + ATAIL_ARR * AHFR
             EHFR = SPSIG(NUMSIG) * AHFR
             EPTOT = EPTOT + ETAIL_ARR * EHFR
           ENDDO
         ENDIF
         TM10 = PI2 * APTOT / EPTOT

         ELSE

           HS   = ZERO
           TM01 = ZERO
           TM02 = ZERO
           TM10 = ZERO
           KLM  = 10.0_rkind
           WLM  = ZERO

         END IF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
!AR: This is all computed in the soruce terms allready we need to get rid of all this shit!
      SUBROUTINE MEAN_FREQS(IP,WALOC,SME01,SME10,ETOTWS,LWINDSEA)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP

         REAL(rkind), INTENT(IN)    :: WALOC(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(OUT)   :: SME01, SME10, ETOTWS

         LOGICAL, INTENT(IN) :: LWINDSEA(NUMSIG,NUMDIR)

         INTEGER             :: ID, IS

         REAL(rkind)                :: ACTOT, ETOT
         REAL(rkind)                :: ETOT_SPSIG

         REAL(rkind)                :: Y(NUMSIG)
         REAL(rkind)                :: DS, ATAIL, ETAIL
         REAL(rkind)                :: tmp(NUMSIG)
!
! total energy ...
!
         ETOT = ZERO
         do id = 1, NUMDIR
           tmp(:) = WALOC(:,id) * spsig 
           ETOT = ETOT + ONEHALF * tmp(1) * ds_incr(1)*ddir
           do is = 2, NUMSIG
             IF (LWINDSEA(IS,ID)) ETOT = ETOT + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir 
           end do
           ETOT = ETOT + ONEHALF * tmp(NUMSIG) * ds_incr(NUMSIG)*ddir
         end do
!
! if etot too small skip ...
!
         if (etot .gt. verysmall) then
!
! integrals ... inlined ... for speed ...
!
           ACTOT = ZERO
           y = ONE/SPSIG
           do id = 1, NUMDIR
             tmp(:) = WALOC(:,id) * spsig * y
             ACTOT = ACTOT + ONEHALF * tmp(1) * ds_incr(1)*ddir
             do is = 2, NUMSIG
               IF (LWINDSEA(IS,ID)) ACTOT = ACTOT + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir 
             end do
             ACTOT = ACTOT + ONEHALF * tmp(NUMSIG) * ds_incr(NUMSIG)*ddir
           end do

           ETOT_SPSIG = ZERO
           y = SIGPOW(:,1) 
           do id = 1, NUMDIR
             tmp(:) = WALOC(:,id) * spsig * y
             ETOT_SPSIG = ETOT_SPSIG + ONEHALF * tmp(1) * ds_incr(1)*ddir
             do is = 2, NUMSIG
               IF (LWINDSEA(IS,ID)) ETOT_SPSIG = ETOT_SPSIG + ONEHALF*(tmp(is)+tmp(is-1))*ds_band(is)*ddir 
             end do
             ETOT_SPSIG = ETOT_SPSIG + ONEHALF * tmp(NUMSIG) * ds_incr(NUMSIG)*ddir
           end do
!
! tail factors ...
!
           DS          = SPSIG(NUMSIG) - SPSIG(NUMSIG-1)

           ATAIL       = SUM(WALOC(NUMSIG,:)) * SIGPOW(NUMSIG,1) * DDIR * DS
           ETAIL       = SUM(WALOC(NUMSIG,:)) * SIGPOW(NUMSIG,2) * DDIR * DS
!
! tail factors ...
!
           ACTOT       = ACTOT        + TAIL_ARR(5)  * ATAIL 
           ETOT        = ETOT         + TAIL_ARR(6)  * ETAIL 
           ETOT_SPSIG  = ETOT_SPSIG   + TAIL_ARR(7)  * ETAIL 
!
! integral parameters ...
!
           SME01       = ETOT_SPSIG / ETOT
           SME10       = ETOT / ACTOT
           ETOTWS      = ETOT

         else
!
! no or too less energy ...
!
           SME01       = ZERO 
           SME10       = ZERO 
           ETOTWS      = ZERO

         end if 
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MEAN_WAVEN(IP,WALOC,KME01,KMWAM,KMWAM2)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP

         REAL(rkind), INTENT(IN)    :: WALOC(NUMSIG,NUMDIR)
         REAL(rkind), INTENT(OUT)   :: KME01
         REAL(rkind), INTENT(OUT)   :: KMWAM, KMWAM2
         INTEGER             :: ID, IS


         REAL(rkind)                :: ACTOT, ETOT
         REAL(rkind)                :: ETOT_SPSIG
         REAL(rkind)                :: ETOT_WK
         REAL(rkind)                :: ETOT_ISQ_WK
         REAL(rkind)                :: ETOT_SQ_WK

         REAL(rkind)                :: Y(NUMSIG), tmp(NUMSIG)
         REAL(rkind)                :: DS, ATAIL, ETAIL, ESIGTAIL
!         REAL(rkind)                :: dintspec, dintspec_y
!
! total energy ...
!
         !ETOT = DINTSPEC(IP,WALOC)
         ETOT = ZERO
         do id = 1, NUMDIR
           tmp(:) = WALOC(:,id) * spsig
           do is = 2, NUMSIG
             ETOT = ETOT + ONEHALF*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
           end do
         end do
!
! if etot too small skip ...
!
         if (etot .gt. verysmall) then
!
! integrals ... inlined ... for speed ...
!
           ACTOT = ZERO
           y = ONE/SPSIG
           do id = 1, NUMDIR
             tmp(:) = WALOC(:,id) * spsig * y
             do is = 2, NUMSIG
               ACTOT = ACTOT + ONEHALF*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
             end do
           end do
           !tmp = ONE/SPSIG
           !ACTOT       = DINTSPEC_Y(IP,WALOC,tmp)
           ETOT_SPSIG = ZERO
           y = SIGPOW(:,1) 
           do id = 1, NUMDIR
             tmp(:) = WALOC(:,id) * spsig * y
             do is = 2, NUMSIG
               ETOT_SPSIG = ETOT_SPSIG + ONEHALF*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
             end do
           end do
           !tmp = SIGPOW(:,1)
           !ETOT_SPSIG  = DINTSPEC_Y(IP,WALOC,tmp)
           ETOT_WK = ZERO
           y = WK(:,IP) 
           do id = 1, NUMDIR
             tmp(:) = WALOC(:,id) * spsig * y
             do is = 2, NUMSIG
               ETOT_WK = ETOT_WK + ONEHALF*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
             end do
           end do
           !tmp = WK(:,IP)
           !ETOT_WK     = DINTSPEC_Y(IP,WALOC,tmp)
           ETOT_ISQ_WK = ZERO 
           y = ONE/SQRT(WK(:,IP))
           do id = 1, NUMDIR
             tmp(:) = WALOC(:,id) * spsig * y
             do is = 2, NUMSIG
               ETOT_ISQ_WK = ETOT_ISQ_WK + ONEHALF*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
             end do
           end do
           !tmp = ONE/SQRT(WK(:,IP))
           !ETOT_ISQ_WK = DINTSPEC_Y(IP,WALOC,tmp)
           ETOT_SQ_WK = ZERO
           y = SQRT(WK(:,IP)) 
           do id = 1, NUMDIR
             tmp(:) = WALOC(:,id) * spsig * y
             do is = 2, NUMSIG
               ETOT_SQ_WK = ETOT_SQ_WK + ONEHALF*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
             end do
           end do
           !tmp = SQRT(WK(:,IP))
           !ETOT_SQ_WK  = DINTSPEC_Y(IP,WALOC,tmp) 
!
! tail factors ...
!
           DS          = SPSIG(NUMSIG) - SPSIG(NUMSIG-1)

           ATAIL       = SUM(WALOC(NUMSIG,:)) * SIGPOW(NUMSIG,1) * DDIR * DS
           ETAIL       = SUM(WALOC(NUMSIG,:)) * SIGPOW(NUMSIG,2) * DDIR * DS
           ESIGTAIL    = SUM(WALOC(NUMSIG,:)) * SIGPOW(NUMSIG,3) * DDIR * DS
!
! tail factors ... borowed from SWAN
!
           ACTOT       = ACTOT        + TAIL_ARR(5)  * ATAIL 
           ETOT        = ETOT         + TAIL_ARR(6)  * ETAIL 
           ETOT_SPSIG  = ETOT_SPSIG   + TAIL_ARR(7)  * ETAIL 
           ETOT_ISQ_WK = ETOT_ISQ_WK  + TAIL_ARR(5)  * ETAIL / (SQRT(WK(NUMSIG,IP)))
           ETOT_SQ_WK  = ETOT_SQ_WK   + TAIL_ARR(5)  * ETAIL * (SQRT(WK(NUMSIG,IP)))
           ETOT_WK     = ETOT_WK      + TAIL_ARR(8)  * ETAIL * WK(NUMSIG,IP)
!
! integral parameters ...
!
           KME01       =  ETOT_WK / ETOT
           KMWAM       = (ETOT/ETOT_ISQ_WK)**2
           KMWAM2      = (ETOT_SQ_WK/ETOT)**2

         else
!
! no or too less energy ...
!
           KME01       = 10.0_rkind
           KMWAM       = 10.0_rkind
           KMWAM2      = 10.0_rkind

         end if 
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
