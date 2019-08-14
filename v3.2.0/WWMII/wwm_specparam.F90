!     Last change:  1     9 Jun 2004    1:44 am
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PARAMENG(IP,ACLOC,SME01,SME10,KME01,KMWAM,KMWAM2,WLM,URSELL,UBOT,ABRBOT,TMBOT,HS,ETOT,FP,TP,CP,KPP,LPP,DM,DSPR,PEAKDSPR,PEAKDM) 

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)    :: IP
         REAL   , INTENT(INOUT) :: ACLOC(MSC,MDC)

         REAL, INTENT(OUT)   :: SME01, SME10
         REAL, INTENT(OUT)   :: KME01, KMWAM, KMWAM2, WLM
         REAL, INTENT(OUT)   :: URSELL
         REAL, INTENT(OUT)   :: UBOT, ABRBOT, TMBOT
         REAL, INTENT(OUT)   :: HS, ETOT, KPP, FP, CP, DM, DSPR, TP, LPP
         REAL, INTENT(OUT)   :: PEAKDSPR, PEAKDM

         INTEGER             :: ID, IS

         REAL                :: SINHKD2(MSC), ACTOTDS(MSC), ETOTD0S(MSC)
         REAL                :: ETOTD1S(MSC)
         REAL                :: ETOTQKD(MSC), ETOTKSS(MSC), ETOTSKD(MSC)
         REAL                :: ETOTKDS(MSC), ETOTQKD2(MSC), ETOT_DSIG(MSC)

         REAL                :: ACTOTDSbis(MSC), ETOTD0Sbis(MSC)
         REAL                :: ETOTSKDbis(MSC), ETOTKDSbis(MSC)

         REAL                :: ACTOT
         REAL                :: DKTOT, EKTOT
         REAL                :: ETOTC4, ETOTS4, PEAKFF
         REAL                :: ETOT1, ESUMAC, HQUOTP, HQUOT
         REAL                :: EAD, DS, EHFR, EFTAIL, DKTOT2
         REAL                :: UB2, AB2, CGP, ETOTF3, ETOTF4
         REAL                :: ETOTS, ETOTC, UB2bis, AB2bis, WVN
         REAL                :: FF, DEG, EDI, CKTAIL, CETAIL
         REAL                :: VEC2DEG, PPTAIL, SKK, SIG2

         KMWAM   = 10.
         KMWAM2  = 10.
         KME01   = 10.
         SME01   = 10.
         SME10   = 10.
         HS      = 0.
         ABRBOT  = 0.001
         UBOT    = 0.
         TMBOT   = 0.

         ETOT = 0.0
         EFTAIL = 1.0 / (PTAIL(1)-1.0)
         ESUMAC = 0.
         DO IS = 1, MSC
           DO ID = 1, MDC
             ESUMAC = ESUMAC + ACLOC(IS,ID)
           END DO
         END DO

         IF (MSC .GE. 2) THEN
            DO ID = 1, MDC
              DO IS = 2, MSC
                 DS = SPSIG(IS) - SPSIG(IS-1)
                 EAD = 0.5*(SPSIG(IS)*ACLOC(IS,ID)+SPSIG(IS-1)*ACLOC(IS-1,ID))*DS*DDIR
                 ETOT = ETOT + EAD
              END DO
              IF (MSC > 3) THEN
                 EHFR = ACLOC(MSC,ID) * SPSIG(MSC)
                 ETOT = ETOT + DDIR * EHFR * SPSIG(MSC) * EFTAIL
              ENDIF
           END DO
         ELSE
           DS = SGHIGH - SGLOW
           DO ID = 1, MDC
              EAD = ACLOC(1,ID) * DS * DDIR
              ETOT = ETOT + EAD
           END DO
         END IF
!
! 2do ... check the influence of ETOT on the results ... this is the swan type integration which i do not like since it is only of 1st order ...
! I better like to use the trapezoid rule or rather simpson type integration, this is the next step for MSC .GT. 3
!
!         ETOT_DSIG(:)  = SUM(ACLOC(:,:),DIM=2) * SIGPOW(:,2) * FDIR 
!         ETOT1         = SUM(ETOT_DSIG)
!         ETOT1 = ETOT1 + ETOT_DSIG(MSC) * PTAIL(6) / FRINTF

!         ETOT = ETOT1
         HS = MAX(VERYSMALL, 4. * SQRT(ETOT))

         IF (ETOT .GT. THR) THEN

            SINHKD2(:) = SINH(MIN(30.,WK(IP,:)*DEP(IP)))**2
            ACTOTDS(:) = SUM(ACLOC(:,:),DIM=2) * SIGPOW(:,1) * FDIR
            ETOTD0S(:) = ACTOTDS(:) * SIGPOW(:,1)
            ETOTD1S(:) = ACTOTDS(:) * SIGPOW(:,2)
            ETOTQKD(:) = ETOTD0S(:) / SQRT(WK(IP,:))
            ETOTQKD2(:)= ETOTD0S(:) * SQRT(WK(IP,:))
            ETOTKSS(:) = ETOTD0S(:) * WK(IP,:)
            ETOTSKD(:) = ETOTD0S(:) / SINHKD2(:)
            ETOTKDS(:) = ETOTSKD(:) * SIGPOW(:,2)

            ACTOT = SUM(ACTOTDS)
            ETOT1 = SUM(ETOTD1S)
            DKTOT = SUM(ETOTQKD)
            DKTOT2= SUM(ETOTQKD2)
            EKTOT = SUM(ETOTKSS)
            UB2   = SUM(ETOTSKD)
            AB2   = SUM(ETOTKDS)

            ACTOTDSbis(:) = SUM(ACLOC(:,:)/ESUMAC,DIM=2) * SIGPOW(:,1)
            ETOTD0Sbis(:) = ACTOTDSbis(:) * SIGPOW(:,1)
            ETOTSKDbis(:) = ETOTD0Sbis(:) / SINHKD2(:)
            ETOTKDSbis(:) = ETOTSKDbis(:) * SIGPOW(:,2)
            UB2bis   = SUM(ETOTSKDbis)
            AB2bis   = SUM(ETOTKDSbis)

            ACTOT   = ACTOT  + PTAIL(5)  * ACTOTDS(MSC) / FRINTF
            ETOT1   = ETOT1  + PTAIL(7)  * ETOTD0S(MSC) * SIGPOW(MSC,1) / FRINTF
            DKTOT   = DKTOT  + PTAIL(5)  * ETOTD0S(MSC) / (SQRT(WK(IP,MSC)) * FRINTF)
            DKTOT2  = DKTOT2 + PTAIL(5)  * ETOTD0S(MSC) * (SQRT(WK(IP,MSC)) * FRINTF)
            EKTOT   = EKTOT  + PTAIL(8)  * ETOTD0S(MSC) * WK(IP,MSC) / FRINTF

            IF (ETOT > VERYSMALL) SME01  = ETOT1 / ETOT
            IF (ETOT > VERySMALL) KME01  = EKTOT / ETOT
            IF (ACTOT > VERySMALL) SME10  = ETOT / ACTOT
            IF (DKTOT > VERySMALL) KMWAM    = (ETOT/DKTOT)**2
            IF (DKTOT2 > VERySMALL) KMWAM2  = (DKTOT2/ETOT)**2
            IF (UB2   > VERYSMALL) UBOT   = SQRT(UB2)
            !IF (UB2   > SMALL) ORBITAL(IP)   = SQRT(UB2)
            IF (AB2   > VERYSMALL) ABRBOT = SQRT(2*AB2)
            IF (UB2bis .GT. THR) THEN
              IF (AB2bis/UB2bis > THR) THEN
                 TMBOT = PI2*SQRT(AB2bis/UB2bis)
              END IF
            END IF
            URSELL = (G9*HS) / (2*SQRT(2.)*SME01**2*DEP(IP)**2)
         ELSE

            HS           = THR
            ABRBOT       = THR
            UBOT         = THR
            TMBOT        = THR
            SME01        = 10.0
            SME10        = 10.0
            KME01        = 10.0
            KMWAM        = 10.0
            KMWAM2       = 10.0
            URSELL       = THR 

         END IF
!
! Peak period continues version... Taken from Thesis Alves ... correct citation is given there ... :)
!
!
! Peak period continues version... Taken from Thesis Alves ... correct citation is given there ... :)
!
         IF (ESUMAC.gt.VERYSMALL) THEN
            ETOTF3 = 0.
            ETOTF4 = 0.
            DO IS = 1, MSC
               DO ID = 1, MDC
                  HQUOT=ACLOC(IS,ID)/ESUMAC
                  HQUOTP=HQUOT**4
                  ETOTF3 = ETOTF3 + SPSIG(IS) * HQUOTP * DS_BAND(IS)
                  ETOTF4 = ETOTF4 +             HQUOTP * DS_BAND(IS)
               END DO
            END DO
            IF(ETOTF4 .GT. VERYSMALL) THEN
               FP = ETOTF3/ETOTF4
               CALL WAVEKCG(DEP(IP), FP, WVN, CP, KPP, CGP)
               TP = 1./(FP/PI2)
               LPP = 1./KPP*PI2
            ELSE
               TP  = 0. 
               FP  = 0.
               CP  = 0. 
               KPP = 10.
               CGP = 0. 
               LPP = 0. 
            END IF
            IF (ETOTF4 .gt. THR) THEN
!AR: Mathieu, I am not sure what is happening here but ETOTS4 is not estimated at all ... this cannot work 
               PEAKDM    = VEC2DEG (ETOTC4, ETOTS4)
               CALL DEG2NAUT(PEAKDM,DEG,LNAUTOUT)
               PEAKDM = DEG
               PEAKFF = MIN (1., SQRT(ETOTC4*ETOTC4+ETOTS4*ETOTS4)/ETOTF4)
               PEAKDSPR = SQRT(2.-2.*PEAKFF) * 180./PI
            ELSE
               FF = 0.
               PEAKDSPR = 0.
               PEAKDM = 0.
            END IF
         ELSE
            TP  = 0. 
            FP  = 0. 
            CP  = 0. 
            KPP = 10.
            CGP = 0. 
            LPP = 0. 
         END IF

!         WRITE(*,*) FP, ETOTF3, ETOTF4
         ETOTC = 0.
         ETOTS = 0.
         ETOT1  = 0.
         DO ID = 1, MDC
           EAD = 0.
           IF (MSC .GE. 2) THEN
             DO  IS = 2, MSC 
               DS  = SPSIG(IS)-SPSIG(IS-1)
               EDI = 0.5*(SPSIG(IS)*ACLOC(IS,ID)+SPSIG(IS-1)*ACLOC(IS-1,ID))*DS
               EAD = EAD + EDI
            END DO
             IF (MSC .GT. 3) THEN
               EHFR = ACLOC(MSC,ID) * SPSIG(MSC)
               EAD = EAD + EHFR * SPSIG(MSC) * EFTAIL
             ENDIF
             EAD = EAD * DDIR
             ETOT1 = ETOT1 + EAD
             ETOTC  = ETOTC + EAD * COSTH(ID)
             ETOTS  = ETOTS + EAD * SINTH(ID)
           ELSE
             DS = SGHIGH - SGLOW
             EAD = ACLOC(1,ID) * DS * DDIR
             EAD = EAD * DDIR
             ETOT1 = ETOT1 + EAD
             ETOTC  = ETOTC + EAD * COSTH(ID)
             ETOTS  = ETOTS + EAD * SINTH(ID)
           END IF
         END DO

         IF (ETOT > THR ) THEN
           DM    = VEC2DEG (ETOTC, ETOTS)
           CALL DEG2NAUT(DM,DEG,LNAUTOUT)
           DM = DEG
           FF = MIN (1., SQRT(ETOTC*ETOTC+ETOTS*ETOTS)/ETOT)
           DSPR = SQRT(2.-2*FF) * 180./PI
         ELSE
           FF = 0.
           DM = 0.
           DSPR = 0.
         END IF

         ETOT1  = 0.
         EKTOT = 0.
         DO IS=1, MSC
            SIG2 = (SPSIG(IS))**2
            SKK  = SIG2 * (WK(IP,IS))**1.!OUTPAR(3)
            DO ID=1,MDC
              ETOT1  = ETOT1 + SIG2 * ACLOC(IS,ID)
              EKTOT = EKTOT + SKK * ACLOC(IS,ID)
            ENDDO
         ENDDO
         ETOT1  = FRINTF * ETOT1
         EKTOT = FRINTF * EKTOT
         IF (MSC .GT. 3) THEN
            PPTAIL = PTAIL(1) - 1.
            CETAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
            PPTAIL = PTAIL(1) - 1. - 2*1.!OUTPAR(3)
            CKTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
            DO ID=1,MDC
              ETOT1   = ETOT1 + CETAIL * SIG2 * ACLOC(MSC,ID)
              EKTOT  = EKTOT + CKTAIL * SKK * ACLOC(MSC,ID)
            ENDDO
         ENDIF
         IF (ETOT.GT.0.) THEN
            WLM = PI2 * (ETOT1 / EKTOT) ** 1.!(1./OUTPAR(3))     
         ELSE
            WLM = 0.
         ENDIF

         RETURN
      END SUBROUTINE
!**********************************************************************:
!*                                                                    *
!**********************************************************************
      SUBROUTINE MEAN_WAVE_PARAMETER(IP,ACLOC,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP

         REAL, INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL, INTENT(OUT)   :: SME01, SME10, KME01
         REAL, INTENT(OUT)   :: KMWAM, KMWAM2
         REAL, INTENT(OUT)   :: HS

         INTEGER             :: ID, IS

         REAL                :: ACTOT, ETOT
         REAL                :: ETOT_SPSIG
         REAL                :: ETOT_WK
         REAL                :: ETOT_ISQ_WK
         REAL                :: ETOT_SKD
         REAL                :: ETOT_SQ_WK
         REAL                :: ETOT_SKDSIG

         REAL                :: Y(MSC)
         REAL                :: DS, ATAIL, ETAIL, ESIGTAIL
         REAL                :: dintspec, dintspec_y, tmp(msc)
!
! total energy ...
! 2do improve efficiency ... 
! 2do check integration style ... 
!
         !ETOT = DINTSPEC(IP,ACLOC)
         ETOT = 0.
         do id = 1, mdc
           tmp(:) = acloc(:,id) * spsig 
           ETOT = ETOT + tmp(1) * 0.5 * ds_incr(1)
           do is = 2, msc - 1
             ETOT = ETOT + 0.5*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
           end do
           ETOT = ETOT + 0.5 * tmp(msc) * ds_incr(msc) 
         end do
!
! if etot too small skip ...
!
         if (etot .gt. verysmall) then
!
! integrals ... inlined ... for speed ...
!
           ACTOT = 0.
           y = 1./SPSIG
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             ACTOT  = ACTOT + tmp(1) * 0.5 * ds_incr(1)
             do is = 2, msc - 1
               ACTOT = ACTOT + 0.5*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
             end do
             ACTOT  = ACTOT + 0.5 * tmp(msc) * ds_incr(msc)
           end do
           !tmp = 1./SPSIG
           !ACTOT       = DINTSPEC_Y(IP,ACLOC,tmp)
           ETOT_SPSIG = 0.
           y = SIGPOW(:,1) 
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             ETOT_SPSIG = ETOT_SPSIG + tmp(1) * 0.5 * ds_incr(1) 
             do is = 2, msc - 1
               ETOT_SPSIG = ETOT_SPSIG + 0.5*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
             end do
             ETOT_SPSIG = ETOT_SPSIG + 0.5 * tmp(msc) * ds_incr(msc)
           end do
           !tmp = SIGPOW(:,1)
           !ETOT_SPSIG  = DINTSPEC_Y(IP,ACLOC,tmp)
           ETOT_WK = 0.
           y = WK(IP,:) 
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             ETOT_WK = ETOT_WK + 0.5 * tmp(1) * ds_incr(1)
             do is = 2, msc - 1
               ETOT_WK = ETOT_WK + 0.5*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
             end do
             ETOT_WK = ETOT_WK + 0.5 * tmp(msc) * ds_incr(msc)
           end do
           !tmp = WK(IP,:)
           !ETOT_WK     = DINTSPEC_Y(IP,ACLOC,tmp)
           ETOT_ISQ_WK = 0. 
           y = 1./SQRT(WK(IP,:))
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             ETOT_ISQ_WK = ETOT_ISQ_WK + 0.5 * tmp(1) * ds_incr(1)
             do is = 2, msc - 1
               ETOT_ISQ_WK = ETOT_ISQ_WK+0.5*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
             end do
             ETOT_ISQ_WK = ETOT_ISQ_WK+0.5*tmp(msc) * ds_incr(msc)
           end do
           !tmp = 1./SQRT(WK(IP,:))
           !ETOT_ISQ_WK = DINTSPEC_Y(IP,ACLOC,tmp)
           ETOT_SQ_WK = 0.
           y = SQRT(WK(IP,:)) 
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             ETOT_SQ_WK = ETOT_SQ_WK + 0.5 * tmp(1) * ds_incr(1)
             do is = 2, msc -1 
               ETOT_SQ_WK = ETOT_SQ_WK + 0.5*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
             end do
             ETOT_SQ_WK = ETOT_SQ_WK + 0.5*tmp(msc) * ds_incr(msc)
           end do
           !tmp = SQRT(WK(IP,:))
           !ETOT_SQ_WK  = DINTSPEC_Y(IP,ACLOC,tmp) 
!
! tail factors ...
!
           DS          = SPSIG(MSC) - SPSIG(MSC-1)

           ATAIL       = SUM(ACLOC(MSC,:)) * SIGPOW(MSC,1) * DDIR * DS
           ETAIL       = SUM(ACLOC(MSC,:)) * SIGPOW(MSC,2) * DDIR * DS
           ESIGTAIL    = SUM(ACLOC(MSC,:)) * SIGPOW(MSC,3) * DDIR * DS
!
! tail factors ... borowed from SWAN
!
           ACTOT       = ACTOT        + PTAIL(5)  * ATAIL 
           ETOT        = ETOT         + PTAIL(6)  * ETAIL 
           ETOT_SPSIG  = ETOT_SPSIG   + PTAIL(7)  * ETAIL 
           ETOT_ISQ_WK = ETOT_ISQ_WK  + PTAIL(5)  * ETAIL / (SQRT(WK(IP,MSC)))
           ETOT_SQ_WK  = ETOT_SQ_WK   + PTAIL(5)  * ETAIL * (SQRT(WK(IP,MSC)))
           ETOT_WK     = ETOT_WK      + PTAIL(8)  * ETAIL * WK(IP,MSC)
!
! integral parameters ...
!
           HS          = MAX(0.,4.*SQRT(ETOT))
           SME01       = ETOT_SPSIG / ETOT
           KME01       = ETOT_WK / ETOT
           SME10       = ETOT / ACTOT
           KMWAM       = (ETOT/ETOT_ISQ_WK)**2
           KMWAM2      = (ETOT_SQ_WK/ETOT)**2

         else
!
! no or too less energy ...
!
           HS          = 0. 
           SME01       = 0. 
           KME01       = 10. 
           SME10       = 0. 
           KMWAM       = 10. 
           KMWAM2      = 10. 

         end if 

       RETURN
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MEAN_PARAMETER(IP,ACLOC,ISMAX,HS,TM01,TM02,KLM,WLM)

         USE DATAPOOL
         IMPLICIT NONE
!2do ... rewrite this integration ...
         INTEGER, INTENT(IN) :: IP,ISMAX

         REAL, INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL, INTENT(OUT)   :: HS,TM01,TM02,KLM,WLM

         INTEGER             :: ID, IS

         REAL                :: Y(MSC)
         REAL                :: DS,ATAIL,ETAIL, ESIGTAIL
         REAL                :: OMEG2,OMEG,EAD,UXD,UYD,ETOT
         REAL                :: EFTAIL,PPTAIL,EFTOT,EPTAIL
         REAL                :: EHFR,AHFR,APTAIL,EPTOT,APTOT
         REAL                :: SKK, CKTAIL, ETOT1, SIG22, EKTOT, CETAIL
         REAL                :: dintspec, dintspec_y, tmp(msc),actmp(msc)
!
! total energy ...
!
         Y = 0.
         tmp = 0.
         actmp = 0.

         ETOT = 0.
         do id = 1, mdc
           actmp(:) = acloc(:,id)
           tmp      = actmp * spsig
           do is = 2, msc
             ETOT = ETOT + 0.5*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
           end do
         end do

         IF (ETOT .GT. THR) THEN
!
! tail ratios same as in swan ...
!
         DS    = SPSIG(MSC) - SPSIG(MSC-1)
         ETAIL = SUM(ACLOC(MSC,:)) * SIGPOW(MSC,2) * DDIR * DS
         ETOT  = ETOT + PTAIL(6) * ETAIL

         HS = 4*SQRT(ETOT)

         APTOT = 0.
         EPTOT = 0.
!
! tail ratios same as in swan ...
!
         PPTAIL = PTAIL(1)
         APTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
         PPTAIL = PTAIL(1) - 1.
         EPTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))

         DO ID = 1, MDC
           DO IS = 1, ISMAX
             APTOT = APTOT + SPSIG(IS) *    ACLOC(IS,ID)
             EPTOT = EPTOT + SPSIG(IS)**2 * ACLOC(IS,ID)
           ENDDO
         ENDDO

         APTOT = APTOT * FRINTF
         EPTOT = EPTOT * FRINTF

         IF (MSC .GT. 3  .AND. .NOT. LSIGMAX) THEN
           DO ID = 1, MDC
           AHFR  = SPSIG(MSC) * ACLOC(MSC,ID)
           APTOT = APTOT + APTAIL * AHFR
           EHFR  = SPSIG(MSC) * AHFR
           EPTOT = EPTOT + EPTAIL * EHFR
           ENDDO
         ENDIF

         IF (EPTOT .GT. 0.) THEN
            TM01 = 2*PI * APTOT / EPTOT
         ELSE
            TM01 = 0.
         END IF


         ETOT  = 0.
         EFTOT = 0.
         PPTAIL = PTAIL(1) - 1.
         ETAIL  = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
         PPTAIL = PTAIL(1) - 3.
         EFTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
         DO ID=1, MDC
            IF (LSECU .OR. LSTCU) THEN
              UXD  = CURTXY(IP,1)*COSTH(ID) + CURTXY(IP,2)*SINTH(ID)
            ENDIF
            DO IS = 1, ISMAX
              EAD  = SPSIG(IS)**2 * ACLOC(IS,ID) * FRINTF
              IF (LSECU .OR. LSTCU) THEN
                OMEG  = SPSIG(IS) + WK(IP,IS) * UXD
                OMEG2 = OMEG**2
              ELSE
               OMEG2 = SPSIG(IS)**2
              ENDIF
              ETOT  = ETOT + EAD
              EFTOT = EFTOT + EAD * OMEG2
            ENDDO
            IF (MSC .GT. 3  .AND. .NOT. LSIGMAX) THEN
              EAD  = SPSIG(MSC)**2 * ACLOC(MSC,ID)
              ETOT  = ETOT  + ETAIL * EAD
              EFTOT = EFTOT + EFTAIL * OMEG2 * EAD
            ENDIF
         ENDDO
         IF (EFTOT .GT. sqrt(verysmall)) THEN
           TM02 = PI2 * SQRT(ETOT/EFTOT)
         ELSE
           TM02 = 0.
         END IF

         ETOT1 = 0.
         EKTOT = 0.
!
! tail ratios same as in swan ...
!
         PPTAIL = PTAIL(1) - 1.
         CETAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
         PPTAIL = PTAIL(1) - 1. - 2*1.
         CKTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))

         DO IS = 1, ISMAX
           SIG22 = (SPSIG(IS))**2
           SKK  = SIG22 * WK(IP,IS)
           DO ID = 1, MDC
             ETOT1 = ETOT1 + SIG22 * ACLOC(IS,ID)
             EKTOT = EKTOT + SKK * ACLOC(IS,ID)
           ENDDO
         ENDDO
         ETOT1 = FRINTF * ETOT1
         EKTOT = FRINTF * EKTOT

         IF (MSC .GT. 3) THEN
            DO ID=1,MDC
              ETOT1 = ETOT1 + CETAIL * SIG22 * ACLOC(MSC,ID)
              EKTOT = EKTOT + CKTAIL * SKK * ACLOC(MSC,ID)
            ENDDO
         ENDIF

         IF (ETOT1.GT.VERYSMALL.AND.EKTOT.GT.VERYSMALL) THEN
            WLM = PI2 * (ETOT1/EKTOT)
            KLM  = PI2/WLM
         ELSE
            KLM  = 10.
            WLM = 0.
         ENDIF

         ELSE

           HS = 0.
           TM01 = 0.
           TM02 = 0.
           KLM  = 10.
           WLM  = 0.

         END IF

       RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MEAN_WAVE_PARAMETER_LOC(ACLOC,CURTXYLOC,DEPLOC,WKLOC,HS,ETOT,SME01,SME10,KME01,KMWAM,KMWAM2)

         USE DATAPOOL
         IMPLICIT NONE
 
         REAL, INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL, INTENT(IN)    :: WKLOC(MSC), DEPLOC
         REAL, INTENT(IN)    :: CURTXYLOC(2)
         REAL, INTENT(OUT)   :: SME01, SME10, KME01
         REAL, INTENT(OUT)   :: KMWAM, KMWAM2
         REAL, INTENT(OUT)   :: HS

         INTEGER             :: ID, IS

         REAL                :: ACTOT, ETOT
         REAL                :: ETOT_SPSIG
         REAL                :: ETOT_WK
         REAL                :: ETOT_ISQ_WK
         REAL                :: ETOT_SKD
         REAL                :: ETOT_SQ_WK
         REAL                :: ETOT_SKDSIG

         REAL                :: Y(MSC)
         REAL                :: DS, ATAIL, ETAIL, ESIGTAIL
         REAL                :: dintspec, dintspec_y, tmp(msc)
!
! total energy ...
!
         ETOT = 0.
         do id = 1, mdc
           tmp(:) = acloc(:,id) * spsig
           ETOT = ETOT + tmp(1) * 0.5 * ds_incr(1)
           do is = 2, msc - 1
             ETOT = ETOT + 0.5*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
           end do
           ETOT = ETOT + 0.5 * tmp(msc) * ds_incr(msc)
         end do
!
! if etot too small skip ...
!
         if (etot .gt. verysmall) then

           ACTOT = 0.
           y = 1./SPSIG
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             ACTOT  = ACTOT + tmp(1) * 0.5 * ds_incr(1)
             do is = 2, msc - 1
               ACTOT = ACTOT + 0.5*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
             end do
             ACTOT  = ACTOT + 0.5 * tmp(msc) * ds_incr(msc)
           end do

           ETOT_SPSIG = 0.
           y = SIGPOW(:,1)
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             ETOT_SPSIG = ETOT_SPSIG + tmp(1) * 0.5 * ds_incr(1)
             do is = 2, msc - 1
               ETOT_SPSIG = ETOT_SPSIG + 0.5*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
             end do
             ETOT_SPSIG = ETOT_SPSIG + 0.5 * tmp(msc) * ds_incr(msc)
           end do

           ETOT_WK = 0.
           y = WKLOC(:)
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             ETOT_WK = ETOT_WK + 0.5 * tmp(1) * ds_incr(1)
             do is = 2, msc - 1
               ETOT_WK = ETOT_WK + 0.5*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
             end do
             ETOT_WK = ETOT_WK + 0.5 * tmp(msc) * ds_incr(msc)
           end do

           ETOT_ISQ_WK = 0.
           y = 1./SQRT(WKLOC)
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             ETOT_ISQ_WK = ETOT_ISQ_WK + 0.5 * tmp(1) * ds_incr(1)
             do is = 2, msc -1 
               ETOT_ISQ_WK = ETOT_ISQ_WK+0.5*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
             end do
             ETOT_ISQ_WK = ETOT_ISQ_WK+0.5*tmp(msc) * ds_incr(msc)
           end do

           ETOT_SQ_WK = 0.
           y = SQRT(WKLOC)
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             ETOT_SQ_WK = ETOT_SQ_WK + 0.5 * tmp(1) * ds_incr(1)
             do is = 2, msc -1 
               ETOT_SQ_WK = ETOT_SQ_WK + 0.5*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
             end do
             ETOT_SQ_WK = ETOT_SQ_WK + 0.5*tmp(msc) * ds_incr(msc)
           end do
!
           DS          = SPSIG(MSC) - SPSIG(MSC-1)

           ATAIL       = SUM(ACLOC(MSC,:)) * SIGPOW(MSC,1) * DDIR * DS
           ETAIL       = SUM(ACLOC(MSC,:)) * SIGPOW(MSC,2) * DDIR * DS
           ESIGTAIL    = SUM(ACLOC(MSC,:)) * SIGPOW(MSC,3) * DDIR * DS
!
! tail factors ... borowed from SWAN
!
           ACTOT       = ACTOT        + PTAIL(5)  * ATAIL 
           ETOT        = ETOT         + PTAIL(6)  * ETAIL 
           ETOT_SPSIG  = ETOT_SPSIG   + PTAIL(7)  * ETAIL 
           ETOT_ISQ_WK = ETOT_ISQ_WK  + PTAIL(5)  * ETAIL / SQRT(WKLOC(MSC))
           ETOT_SQ_WK  = ETOT_SQ_WK   + PTAIL(5)  * ETAIL * SQRT(WKLOC(MSC))
           ETOT_WK     = ETOT_WK      + PTAIL(8)  * ETAIL * WKLOC(MSC) 
!
! integral parameters ...
!
           HS          = MAX(0.,4.*SQRT(ETOT))
           SME01       = ETOT_SPSIG / ETOT
           KME01       = ETOT_WK / ETOT
           SME10       = ETOT / ACTOT
           KMWAM       = (ETOT/ETOT_ISQ_WK)**2
           KMWAM2      = (ETOT_SQ_WK/ETOT)**2

         else
!
! no or too less energy ...
!
           HS          = 0. 
           SME01       = 0. 
           KME01       = 10. 
           SME10       = 0. 
           KMWAM       = 10. 
           KMWAM2      = 10. 

         end if 

       RETURN
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAVE_CURRENT_PARAMETER(IP,ACLOC,UBOT,ORBITAL,BOTEXPER,TMBOT)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP
         REAL,    INTENT(IN) :: ACLOC(MSC,MDC)

         REAL, INTENT(OUT)   :: UBOT, ORBITAL, BOTEXPER, TMBOT

         INTEGER             :: ID, IS

         REAL                :: ETOT_SKD
         REAL                :: ETOT_SKDSIG, TMP(MSC), Y(MSC)

         REAL                :: dintspec, dintspec_y
!
! integrals ...
!
         IF (DEP(IP) .LT. DMIN) RETURN
         
         ETOT_SKD    = 0.
         ETOT_SKDSIG = 0.

         y = 1./SINH(MIN(20.,WK(IP,:)*DEP(IP)))**2

         do id = 1, mdc
           tmp(:) = acloc(:,id) * spsig * y
           ETOT_SKD  = ETOT_SKD + tmp(1) * 0.5 * ds_incr(1)
           do is = 2, msc -1 
             ETOT_SKD = ETOT_SKD + 0.5*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
           end do
           ETOT_SKD = ETOT_SKD + tmp(msc) * 0.5 * ds_incr(msc)
         end do

         y =  SIGPOW(:,2)*1./SINH(MIN(20.,WK(IP,:)*DEP(IP)))**2

         do id = 1, mdc
           tmp(:) = acloc(:,id) * spsig * y
           ETOT_SKDSIG = ETOT_SKDSIG + tmp(1) * 0.5 * ds_incr(1)
           do is = 2, msc -1 
             ETOT_SKDSIG = ETOT_SKDSIG + 0.5*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
           end do 
           ETOT_SKDSIG = tmp(msc) * 0.5 * ds_incr(msc)
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
           UBOT        = 0. 
           ORBITAL     = 0. 
           BOTEXPER    = 0. 
           TMBOT       = 0. 

         ENDIF 

       RETURN
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WAVE_CURRENT_PARAMETER_LOC(ACLOC,CURTXYLOC,DEPLOC,WKLOC,UBOT,ORBITAL,BOTEXPER,TMBOT)
         USE DATAPOOL
         IMPLICIT NONE

         REAL, INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL, INTENT(IN)    :: WKLOC(MSC), DEPLOC
         REAL, INTENT(IN)    :: CURTXYLOC(2)

         REAL, INTENT(OUT)   :: UBOT, ORBITAL, BOTEXPER, TMBOT

         INTEGER             :: ID, IS

         LOGICAL             :: ISINF

         REAL                :: ETOT_SKD
         REAL                :: ETOT_SKDSIG, TMP(MSC), Y(MSC)

         REAL                :: dintspec, dintspec_y
!
! integrals ...
!
         IF (DEPLOC .LT. DMIN) RETURN

         ETOT_SKD    = 0.
         ETOT_SKDSIG = 0.

         y = 1./SINH(MIN(20.,WKLOC*DEPLOC))**2
         do id = 1, mdc
           tmp(:) = acloc(:,id) * spsig * y
           ETOT_SKD  = ETOT_SKD + tmp(1) * 0.5 * ds_incr(1)
           do is = 2, msc -1
             ETOT_SKD = ETOT_SKD + 0.5*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
           end do
           ETOT_SKD = ETOT_SKD + tmp(msc) * 0.5 * ds_incr(msc)
         end do

         y =  SIGPOW(:,2)*1./SINH(MIN(20.,WKLOC*DEPLOC))**2
         do id = 1, mdc
           tmp(:) = acloc(:,id) * spsig * y
           ETOT_SKDSIG = ETOT_SKDSIG + tmp(1) * 0.5 * ds_incr(1)
           do is = 2, msc -1
             ETOT_SKDSIG = ETOT_SKDSIG + 0.5*(tmp(is)+tmp(is-1))*ds_band(is)*ddir
           end do
           ETOT_SKDSIG = tmp(msc) * 0.5 * ds_incr(msc)
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
           UBOT        = 0. 
           ORBITAL     = 0. 
           BOTEXPER    = 0. 
           TMBOT       = 0. 

         ENDIF 

         !WRITE(*,'(9F15.4)') DEPLOC,SUM(ACLOC),CURTXYLOC,SUM(WKLOC), ETOT_SKD, ETOT_SKDSIG 

       RETURN
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE URSELL_NUMBER(HS,SME,DEPTH,URSELL)
         USE DATAPOOL, ONLY : G9, DMIN, verysmall

         IMPLICIT NONE

         REAL, INTENT(IN)    :: HS, SME, DEPTH
         REAL, INTENT(OUT)   :: URSELL 

         IF (DEPTH .GT. DMIN .AND. SME .GT. verysmall .AND. HS .GT. verysmall) THEN
           URSELL = (G9 * HS) * 1./(2.*SQRT(2.)*SME**2*DEPTH**2)
         ELSE
           URSELL = 0.
         END IF 

       RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE WINDSEASWELLSEP( IP, ACLOC, TM_W, CGP_W, CP_W, TP_W, LP_W, HS_W, KP_W )

      ! swan routine

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)    :: IP
         REAL   , INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL   , INTENT(OUT)   :: TM_W, CGP_W, CP_W, TP_W, LP_W, HS_W, KP_W

         INTEGER                :: ID, IS

         REAL                   :: ETOT
         REAL                   :: VEC2RAD, EFTOT, OMEG, WINDTH
         REAL                   :: EAD, DS, EHFR, EFTAIL, ETAIL, PPTAIL, ACWIND(MSC,MDC)
         REAL                   :: UXD, ETOTF3, ETOTF4, FP, WN_W, WVC, WKDEP_W

         WINDTH = VEC2RAD(WINDXY(IP,1),WINDXY(IP,2))

         ETOT = 0.0
         EFTAIL = 1.0 / (PTAIL(1)-1.0)

         DO ID  = 1, MDC            ! Calculate wind sea energy ... weak criterion
           DO IS = 1, MSC
             WVC = REAL(SPSIG(IS)/WK(IP,IS))
             IF (  1.2*UFRIC(IP)*COS(SPDIR(ID)-WINDTH)*(28./WVC) .LT. 1.) THEN
               ACWIND(IS,ID) = 0.   ! Swell
             ELSE
               ACWIND(IS,ID) = ACLOC(IS,ID)  ! Wind Sea
             END IF
           END DO
           DO IS = 2, MSC
              DS = SPSIG(IS) - SPSIG(IS-1)
              EAD = 0.5*(SPSIG(IS)*ACWIND(IS,ID)+SPSIG(IS-1)*ACWIND(IS-1,ID))*DS*DDIR
              ETOT = ETOT + EAD
           END DO
           IF (MSC > 3) THEN
             EHFR = ACLOC(MSC,ID) * SPSIG(MSC)
             ETOT = ETOT + DDIR * EHFR * SPSIG(MSC) * EFTAIL
           ENDIF
         END DO

         IF (ETOT .LT. VERYSMALL) GOTO 101

         IF (ETOT > 0.0) THEN
           HS_W = 4.0*SQRT(ETOT)
         ELSE
           HS_W = 0.0
         END IF

         ETOTF3 = 0.
         ETOTF4 = 0.
         DO IS = 1, MSC
           DO ID = 1, MDC
             ETOTF3 = ETOTF3 + SPSIG(IS) * ACLOC(IS,ID)**4 * DDIR * DS_BAND(IS)
             ETOTF4 = ETOTF4 +             ACLOC(IS,ID)**4 * DDIR * DS_BAND(IS)
           END DO
         END DO
         IF(ETOTF4 .GT. VERYSMALL) THEN
            FP = ETOTF3/ETOTF4
            TP_W = 1./FP/PI2
            !CALL WAVEKCG(DEP(IP), FP, WN_W, CP_W, KP_W, CGP_W)
            CALL ALL_FROM_TABLE(FP,DEP(IP),KP_W,CGP_W,WKDEP_W,WN_W,CP_W)
            LP_W  = PI2/KP_w
         ELSE
            FP    = 0. 
            CP_W  = 0. 
            KP_W  = 10. 
            CGP_W = 0. 
            LP_W  = 0. 
         END IF

         ETOT = 0.
         EFTOT = 0.
         PPTAIL = PTAIL(1) - 1.
         ETAIL  = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
         PPTAIL = PTAIL(1) - 2.
         EFTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
         DO ID = 1, MDC
            UXD = CURTXY(IP,1)*COSTH(ID) + CURTXY(IP,2)*SINTH(ID)
            DO IS = 1, MSC
              OMEG = SPSIG(IS) + WK(IP,IS) * UXD
              EAD = FRINTF * SPSIG(IS)**2 * ACWIND(IS,ID)
              ETOT = ETOT + EAD
              EFTOT = EFTOT + EAD * OMEG
            ENDDO
            IF (MSC .GT. 3) THEN
              EAD = SPSIG(MSC)**2 * ACWIND(MSC,ID)
              ETOT = ETOT + ETAIL * EAD
             EFTOT = EFTOT + EFTAIL * OMEG * EAD
            ENDIF
         ENDDO
         IF (EFTOT.GT.0.) THEN
            TM_W = PI2 * ETOT / EFTOT
         ELSE
            TM_W = 0.
         ENDIF

 
101      RETURN
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PEAK_PARAMETER(IP,ACLOC,ISMAX,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)    :: IP, ISMAX
         REAL,    INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL  , INTENT(OUT)    :: KPP,CPP,FPP,WNPP,CGPP,TPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD

         INTEGER                :: IS, ID, IDIRM, ISIGMP
         REAL                   :: HQUOT, HQUOTP, ETOTF3, ETOTF4, ETOTC4, ETOTS4, PEAKFF,WKDEPD,WNPD
         REAL                   :: FF, DEG, VEC2DEG, MAXAC, E1, E2, DS, EAD, ETOTT, WKDEPP
!
! Peak period continues version... Taken from Thesis Henriques Alves ... correct citation is given there ... :)
!
       MAXAC = MAXVAL(ACLOC)

       IF (MAXAC .gt. VERYSMALL .AND.  DEP(IP) .GT. DMIN) THEN

         ETOTF3 = 0.
         ETOTF4 = 0.
         ETOTC4 = 0.
         ETOTS4 = 0.

         DO IS = 1, MSC
           DO ID = 1, MDC
             HQUOT  = ACLOC(IS,ID)/MAXAC
             HQUOTP = HQUOT**4
             ETOTF3 = ETOTF3 + SPSIG(IS) * HQUOTP * DS_BAND(IS)
             ETOTF4 = ETOTF4 +             HQUOTP * DS_BAND(IS)
             ETOTC4 = ETOTC4 + COSTH(ID) * HQUOTP * DS_BAND(IS)
             ETOTS4 = ETOTS4 + SINTH(ID) * HQUOTP * DS_BAND(IS)
           END DO
         END DO

         IF(ETOTF4 .GT. VERYSMALL .AND. ETOTF4 .GT. VERYSMALL) THEN

           FPP    = ETOTF3/ETOTF4
           !CALL WAVEKCG(DEP(IP), FPP, WNPP, CPP, KPP, CGPP)
           CALL ALL_FROM_TABLE(FPP,DEP(IP),KPP,CGPP,WKDEPP,WNPP,CPP)
           TPP    = PI2/FPP
           LPP    = PI2/KPP
           PEAKDM = VEC2DEG (ETOTC4, ETOTS4)
           CALL DEG2NAUT(PEAKDM,DEG,LNAUTOUT)
           PEAKDM = DEG
           PEAKFF = MIN(1.,SQRT(ETOTC4*ETOTC4+ETOTS4*ETOTS4)/ETOTF4)
           PEAKDSPR = SQRT(MAX(0.,2.-2.*PEAKFF)) * 180./PI

         ELSE 

           FPP = 0.
           KPP = 10.
           CGPP = 0. 
           WKDEPP = 0.
           WNPP = 0.
           CPP = 0.
           TPP = 0.
           LPP = 0.
           PEAKDM = 0.
           PEAKFF = 0.
           PEAKDSPR = 0. 

         END IF

         DPEAK = 1
         ETOTT = 0.0
         IDIRM = -1
         DO ID = 1, MDC
            EAD = 0.0
            DO IS = 2, ISMAX
               DS = SPSIG(IS)-SPSIG(IS-1)
               E1 = SPSIG(IS-1)*ACLOC(IS-1,ID)
               E2 = SPSIG(IS)*ACLOC(IS,ID)
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
           DPEAK    = 0.
         END IF

       ELSE

         FPP = 0.
         KPP = 10.
         CGPP = 0.
         WKDEPP = 0.
         WNPP = 0.
         CPP = 0.
         TPP = 0.
         LPP = 0.
         PEAKDM = 0.
         PEAKFF = 0.
         PEAKDSPR = 0.
         DPEAK = 0.

       END IF

       TPPD = 0.0
       CPPD = 0.0
       KPPD = 0.0
       CGPD = 0.0
       ETOTT = 0.0
       ISIGMP = -1
       DO IS = 1, MSC
         EAD = 0.0
         DO ID = 1, MDC
            EAD = EAD + SPSIG(IS)*ACLOC(IS,ID)*DDIR
         ENDDO
         IF (EAD > ETOTT) THEN
           ETOTT = EAD
           ISIGMP = IS
          END IF
        END DO
        IF (ISIGMP > 0) THEN
          TPPD = 1./(SPSIG(ISIGMP)/PI2)
          !CALL WAVEKCG(DEP(IP), SPSIG(ISIGMP), CPPD, KPPD, CGPD)
          CALL ALL_FROM_TABLE(SPSIG(ISIGMP),DEP(IP),KPPD,CGPD,WKDEPD,WNPD,CPPD)
        ELSE
          TPPD = 0.0
          CPPD  = 0.0
          KPPD  = 0.0
          CGPD  = 0.0
       END IF

       !WRITE(*,'(11F15.4)') FPP, KPP, CGPP, WKDEPP, WNPP, CPP, TPP, LPP, PEAKDM, PEAKFF, PEAKDSPR 

       RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PEAK_PARAMETER_LOC(ACLOC,DEPLOC,ISMAX,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)    :: ISMAX
         REAL, INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL, INTENT(IN)    :: DEPLOC


         REAL  , INTENT(OUT)    :: KPP,CPP,WNPP,CGPP,TPP,LPP,PEAKDSPR,PEAKDM,DPEAK,TPPD,KPPD,CGPD,CPPD,FPP

         INTEGER                :: IS, ID, IDIRM, ISIGMP
         REAL                   :: HQUOT, HQUOTP, ETOTF3, ETOTF4, ETOTC4, ETOTS4, PEAKFF, WKDEPP
         REAL                   :: FF, DEG, VEC2DEG, MAXAC, E1, E2, DS, EAD, ETOTT,WKDEPD,WNPD
!
! Peak period continues version... Taken from Thesis Henriques Alves ... correct citation is given there ... :)
!
       MAXAC = MAXVAL(ACLOC)

       IF (MAXAC .gt. VERYSMALL .AND. DEPLOC .GT. DMIN) THEN
         ETOTF3 = 0.
         ETOTF4 = 0.
         ETOTC4 = 0.
         ETOTS4 = 0.
         DO IS = 1, MSC
           DO ID = 1, MDC
             HQUOT  = ACLOC(IS,ID)/MAXAC
             HQUOTP = HQUOT**4
             ETOTF3 = ETOTF3 + SPSIG(IS) * HQUOTP * DS_BAND(IS)
             ETOTF4 = ETOTF4 +             HQUOTP * DS_BAND(IS)
             ETOTC4 = ETOTC4 + COSTH(ID) * HQUOTP * DS_BAND(IS)
             ETOTS4 = ETOTS4 + SINTH(ID) * HQUOTP * DS_BAND(IS)
           END DO
         END DO
         IF(ETOTF4 .GT. VERYSMALL .AND. ETOTF4 .GT. VERYSMALL) THEN
           FPP    = ETOTF3/ETOTF4
           !CALL WAVEKCG(DEPLOC, FPP, WNPP, CPP, KPP, CGPP)
           CALL ALL_FROM_TABLE(FPP,DEPLOC,KPP,CGPP,WKDEPP,WNPP,CPP) 
           PEAKDM = VEC2DEG (ETOTC4, ETOTS4)
           TPP    = PI2/FPP
           LPP    = PI2/KPP
           CALL DEG2NAUT(PEAKDM,DEG,LNAUTOUT)
           PEAKDM = DEG
           PEAKFF = MIN(1.,SQRT(MAX(0.,ETOTC4*ETOTC4+ETOTS4*ETOTS4))/ETOTF4)
           PEAKDSPR = SQRT(MAX(0.,2.-2.*PEAKFF)) * 180./PI
         ELSE
           FPP = 0.
           KPP = 10.
           CGPP = 0.
           WKDEPP = 0.
           WNPP = 0.
           CPP = 0.
           TPP = 0.
           LPP = 0.
           PEAKDM = 0.
           PEAKFF = 0.
           PEAKDSPR = 0.
         END IF
         DPEAK = 1
         ETOTT = 0.0
         IDIRM = -1
         DO ID = 1, MDC
            EAD = 0.0
            DO IS = 2, ISMAX
               DS = SPSIG(IS)-SPSIG(IS-1)
               E1 = SPSIG(IS-1)*ACLOC(IS-1,ID)
               E2 = SPSIG(IS)*ACLOC(IS,ID)
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
           DPEAK = 0.
         END IF
       ELSE
         FPP = 0.
         KPP = 10.
         CGPP = 0.
         WKDEPP = 0.
         WNPP = 0.
         CPP = 0.
         TPP = 0.
         LPP = 0.
         PEAKDM = 0.
         PEAKFF = 0.
         PEAKDSPR = 0.
         DPEAK = 0.
       END IF

       TPPD = 0.0
       CPPD = 0.0
       KPPD = 0.0
       CGPD = 0.0
       ETOTT = 0.0
       ISIGMP = -1
       DO IS = 1, MSC
         EAD = 0.0
         DO ID = 1, MDC
            EAD = EAD + SPSIG(IS)*ACLOC(IS,ID)*DDIR
         ENDDO
         IF (EAD > ETOTT) THEN
           ETOTT = EAD
           ISIGMP = IS
          END IF
        END DO
        IF (ISIGMP > 0) THEN
          TPPD = 1./(SPSIG(ISIGMP)/PI2)
          !CALL WAVEKCG(DEPLOC, SPSIG(ISIGMP), CPPD, KPPD, CGPD)
          CALL ALL_FROM_TABLE(SPSIG(ISIGMP),DEPLOC,KPPD,CGPD,WKDEPD,WNPD,CPPD)
        ELSE
          TPPD = 0.0
          CPPD  = 0.0
          KPPD  = 0.0
          CGPD  = 0.0
       END IF

       RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MEAN_DIRECTION_AND_SPREAD(IP,ACLOC,ISMAX,ETOTS,ETOTC,DM,DSPR)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)    :: IP, ISMAX
         REAL,    INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL,    INTENT(OUT)   :: ETOTS, ETOTC, DM, DSPR

         INTEGER                :: IS, ID

         REAL                   :: DS, EDI, EAD, ETOT1, EHFR
         REAL                   :: EFTAIL, VEC2DEG, DEG, FF

         ETOTC = 0.
         ETOTS = 0.
         ETOT1  = 0.

         EFTAIL = 1.0 / (PTAIL(1)-1.0)

         DO ID = 1, MDC
           EAD = 0.
             DO  IS = 2, ISMAX 
               DS  = SPSIG(IS)-SPSIG(IS-1)
               EDI = 0.5*(SPSIG(IS)*ACLOC(IS,ID)+SPSIG(IS-1)*ACLOC(IS-1,ID))*DS
               EAD = EAD + EDI
            END DO
            IF (MSC .GT. 3) THEN
              EHFR = ACLOC(MSC,ID) * SPSIG(MSC)
              EAD = EAD + EHFR * SPSIG(MSC) * EFTAIL
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
           FF = MIN (1., SQRT(MAX(0.,ETOTC*ETOTC+ETOTS*ETOTS))/ETOT1)
           DSPR = SQRT(MAX(0.,2.-2.*FF)) * 180./PI
         ELSE
           FF   = 0.
           DM   = 0.
           DSPR = 0.
         END IF

       RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MEAN_DIRECTION_AND_SPREAD_LOC(ACLOC,ISMAX,ETOTS,ETOTC,DM,DSPR)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)    :: ISMAX
         REAL,    INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL,    INTENT(OUT)   :: ETOTS, ETOTC, DM, DSPR

         INTEGER                :: IS, ID

         REAL                   :: DS, EDI, EAD, ETOT1, EHFR
         REAL                   :: EFTAIL, VEC2DEG, DEG, FF

         ETOTC = 0.
         ETOTS = 0.
         ETOT1  = 0.

         EFTAIL = 1.0 / (PTAIL(1)-1.0)

         DO ID = 1, MDC
           EAD = 0.
             DO  IS = 2, ISMAX 
               DS  = SPSIG(IS)-SPSIG(IS-1)
               EDI = 0.5*(SPSIG(IS)*ACLOC(IS,ID)+SPSIG(IS-1)*ACLOC(IS-1,ID))*DS
               EAD = EAD + EDI
            END DO
            IF (MSC .GT. 3) THEN
              EHFR = ACLOC(MSC,ID) * SPSIG(MSC)
              EAD = EAD + EHFR * SPSIG(MSC) * EFTAIL
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
           FF = MIN (1., SQRT(MAX(0.,ETOTC*ETOTC+ETOTS*ETOTS))/ETOT1)
           DSPR = SQRT(MAX(0.,2.-2.*FF)) * 180./PI
         ELSE
           FF   = 0.
           DM   = 0.
           DSPR = 0.
         END IF

       RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MEAN_PARAMETER_LOC(ACLOC,CURTXYLOC,DEPLOC,WKLOC,ISMAX,HS,TM01,TM02,KLM,WLM)
         USE DATAPOOL
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: ISMAX
         REAL, INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL, INTENT(IN)    :: WKLOC(MSC), DEPLOC
         REAL, INTENT(IN)    :: CURTXYLOC(2)


         REAL, INTENT(OUT)   :: HS,TM01,TM02,KLM,WLM

         INTEGER             :: ID, IS

         REAL                :: Y(MSC)
         REAL                :: DS,ATAIL,ETAIL, ESIGTAIL
         REAL                :: OMEG2,OMEG,EAD,UXD,UYD,ETOT
         REAL                :: EFTAIL,PPTAIL,EFTOT,EPTAIL
         REAL                :: EHFR,AHFR,APTAIL,EPTOT,APTOT
         REAL                :: SKK, CKTAIL, ETOT1, SIG22, EKTOT, CETAIL
         REAL                :: dintspec, dintspec_y, tmp(msc),actmp(msc)
!
! total energy ...
!
         ETOT = 0.
         do id = 1, mdc
           actmp(:) = acloc(:,id)
           tmp      = actmp * spsig
           do is = 2, msc
             ETOT = ETOT + 0.5*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
           end do
         end do

         IF (ETOT .GT. verysmall) THEN
!
! tail ratios same as in swan ...
!
         DS    = SPSIG(MSC) - SPSIG(MSC-1)
         ETAIL = SUM(ACLOC(MSC,:)) * SIGPOW(MSC,2) * DDIR * DS
         ETOT  = ETOT + PTAIL(6) * ETAIL

         HS = 4*SQRT(ETOT)

         APTOT = 0.
         EPTOT = 0.
!
! tail ratios same as in swan ...
!
         PPTAIL = PTAIL(1)
         APTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
         PPTAIL = PTAIL(1) - 1.
         EPTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))

         DO ID = 1, MDC
           DO IS = 1, ISMAX
             APTOT = APTOT + SPSIG(IS) *    ACLOC(IS,ID)
             EPTOT = EPTOT + SPSIG(IS)**2 * ACLOC(IS,ID)
           ENDDO
         ENDDO

         APTOT = APTOT * FRINTF
         EPTOT = EPTOT * FRINTF

         IF (MSC .GT. 3  .AND. .NOT. LSIGMAX) THEN
           DO ID = 1, MDC
           AHFR  = SPSIG(MSC) * ACLOC(MSC,ID)
           APTOT = APTOT + APTAIL * AHFR
           EHFR  = SPSIG(MSC) * AHFR
           EPTOT = EPTOT + EPTAIL * EHFR
           ENDDO
         ENDIF

         IF (EPTOT .GT. VERYSMALL) THEN
            TM01 = 2.*PI * APTOT / EPTOT
         ELSE
            TM01 = 0.
         END IF

         ETOT  = 0.
         EFTOT = 0.

         PPTAIL = PTAIL(1) - 1.
         ETAIL  = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
         PPTAIL = PTAIL(1) - 3.
         EFTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
!
! tail ratios same as in swan ...
!
         DO ID=1, MDC
            IF (LSECU .OR. LSTCU) THEN
              UXD  = CURTXYLOC(1)*COSTH(ID) + CURTXYLOC(2)*SINTH(ID)
            ENDIF
            DO IS = 1, ISMAX
              EAD  = SPSIG(IS)**2 * ACLOC(IS,ID) * FRINTF
              IF (LSECU .OR. LSTCU) THEN
                OMEG  = SPSIG(IS) + WKLOC(IS) * UXD
                OMEG2 = OMEG**2
              ELSE
               OMEG2 = SPSIG(IS)**2
              ENDIF
              ETOT  = ETOT + EAD
              EFTOT = EFTOT + EAD * OMEG2
            ENDDO
            IF (MSC .GT. 3  .AND. .NOT. LSIGMAX) THEN
              EAD  = SPSIG(MSC)**2 * ACLOC(MSC,ID)
              ETOT  = ETOT  + ETAIL * EAD
              EFTOT = EFTOT + EFTAIL * OMEG2 * EAD
            ENDIF
         ENDDO
         IF (EFTOT .GT. verysmall .AND. ETOT .GT. verysmall) THEN
           TM02 = PI2 * SQRT(ETOT/EFTOT)
         ELSE
           TM02 = 0.
         END IF

         ETOT1 = 0.
         EKTOT = 0.
!
! tail ratios same as in swan ...
!
         PPTAIL = PTAIL(1) - 1.
         CETAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
         PPTAIL = PTAIL(1) - 1. - 2.*1.
         CKTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))

         DO IS = 1, ISMAX
           SIG22 = (SPSIG(IS))**2
           SKK  = SIG22 * WKLOC(IS)
           DO ID = 1, MDC
             ETOT1 = ETOT1 + SIG22 * ACLOC(IS,ID)
             EKTOT = EKTOT + SKK * ACLOC(IS,ID)
           ENDDO
         ENDDO

         ETOT1 = FRINTF * ETOT1
         EKTOT = FRINTF * EKTOT

         IF (MSC .GT. 3) THEN
            DO ID=1,MDC
              ETOT1 = ETOT1 + CETAIL * SIG22 * ACLOC(MSC,ID)
              EKTOT = EKTOT + CKTAIL * SKK * ACLOC(MSC,ID)
            ENDDO
         ENDIF

         IF (ETOT1.GT.VERYSMALL.AND.EKTOT.GT.VERYSMALL) THEN
            WLM = PI2 * (ETOT1/EKTOT)
            KLM  = PI2/WLM
         ELSE
            KLM  = 10.
            WLM = 0.
         ENDIF

         ELSE

           HS   = 0.
           TM01 = 0.
           TM02 = 0.
           KLM  = 10.
           WLM  = 0.

         END IF

       RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MEAN_FREQS(IP,ACLOC,SME01,SME10,ETOTWS,LWINDSEA)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP

         REAL, INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL, INTENT(OUT)   :: SME01, SME10, ETOTWS

         LOGICAL, INTENT(IN) :: LWINDSEA(MSC,MDC)

         INTEGER             :: ID, IS

         REAL                :: ACTOT, ETOT
         REAL                :: ETOT_SPSIG
         REAL                :: ETOT_WK
         REAL                :: ETOT_ISQ_WK
         REAL                :: ETOT_SKD
         REAL                :: ETOT_SQ_WK
         REAL                :: ETOT_SKDSIG

         REAL                :: Y(MSC)
         REAL                :: DS, ATAIL, ETAIL, ESIGTAIL
         REAL                :: dintspec, dintspec_y, tmp(msc)
!
! total energy ...
!
         ETOT = 0.
         do id = 1, mdc
           tmp(:) = acloc(:,id) * spsig 
           ETOT = ETOT + 0.5 * tmp(1) * ds_incr(1)
           do is = 2, msc - 1
             IF (LWINDSEA(IS,ID)) ETOT = ETOT + 0.5*(tmp(is)+tmp(is-1))*ds_band(is)*ddir 
           end do
           ETOT = ETOT + 0.5 * tmp(msc) * ds_incr(msc)
         end do
!
! if etot too small skip ...
!
         if (etot .gt. verysmall) then
!
! integrals ... inlined ... for speed ...
!
           ACTOT = 0.
           y = 1./SPSIG
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             ACTOT = ACTOT + 0.5 * tmp(1) * ds_incr(1)
             do is = 2, msc - 1
               IF (LWINDSEA(IS,ID)) ACTOT = ACTOT + 0.5*(tmp(is)+tmp(is-1))*ds_band(is)*ddir 
             end do
             ACTOT = ACTOT + 0.5 * tmp(msc) * ds_incr(msc)
           end do

           ETOT_SPSIG = 0.
           y = SIGPOW(:,1) 
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             ETOT_SPSIG = ETOT_SPSIG + 0.5 * tmp(1) * ds_incr(1)
             do is = 2, msc
               IF (LWINDSEA(IS,ID)) ETOT_SPSIG = ETOT_SPSIG + 0.5*(tmp(is)+tmp(is-1))*ds_band(is)*ddir 
             end do
             ETOT_SPSIG = ETOT_SPSIG + 0.5 * tmp(msc) * ds_incr(msc)
           end do
!
! tail factors ...
!
           DS          = SPSIG(MSC) - SPSIG(MSC-1)

           ATAIL       = SUM(ACLOC(MSC,:)) * SIGPOW(MSC,1) * DDIR * DS
           ETAIL       = SUM(ACLOC(MSC,:)) * SIGPOW(MSC,2) * DDIR * DS
!
! tail factors ... borowed from SWAN
!
           ACTOT       = ACTOT        + PTAIL(5)  * ATAIL 
           ETOT        = ETOT         + PTAIL(6)  * ETAIL 
           ETOT_SPSIG  = ETOT_SPSIG   + PTAIL(7)  * ETAIL 
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
           SME01       = 0. 
           SME10       = 0. 
           ETOTWS      = 0.

         end if 

       RETURN
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MEAN_WAVEN(IP,ACLOC,KME01,KMWAM,KMWAM2)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: IP

         REAL, INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL, INTENT(OUT)   :: KME01
         REAL, INTENT(OUT)   :: KMWAM, KMWAM2

         INTEGER             :: ID, IS

         REAL                :: ACTOT, ETOT
         REAL                :: ETOT_SPSIG
         REAL                :: ETOT_WK
         REAL                :: ETOT_ISQ_WK
         REAL                :: ETOT_SKD
         REAL                :: ETOT_SQ_WK
         REAL                :: ETOT_SKDSIG

         REAL                :: Y(MSC)
         REAL                :: DS, ATAIL, ETAIL, ESIGTAIL
         REAL                :: dintspec, dintspec_y, tmp(msc)
!
! total energy ...
!
         !ETOT = DINTSPEC(IP,ACLOC)
         ETOT = 0.
         do id = 1, mdc
           tmp(:) = acloc(:,id) * spsig 
           do is = 2, msc
             ETOT = ETOT + 0.5*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
           end do
         end do
!
! if etot too small skip ...
!
         if (etot .gt. verysmall) then
!
! integrals ... inlined ... for speed ...
!
           ACTOT = 0.
           y = 1./SPSIG
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             do is = 2, msc
               ACTOT = ACTOT + 0.5*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
             end do
           end do
           !tmp = 1./SPSIG
           !ACTOT       = DINTSPEC_Y(IP,ACLOC,tmp)
           ETOT_SPSIG = 0.
           y = SIGPOW(:,1) 
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             do is = 2, msc
               ETOT_SPSIG = ETOT_SPSIG + 0.5*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
             end do
           end do
           !tmp = SIGPOW(:,1)
           !ETOT_SPSIG  = DINTSPEC_Y(IP,ACLOC,tmp)
           ETOT_WK = 0.
           y = WK(IP,:) 
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             do is = 2, msc
               ETOT_WK = ETOT_WK + 0.5*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
             end do
           end do
           !tmp = WK(IP,:)
           !ETOT_WK     = DINTSPEC_Y(IP,ACLOC,tmp)
           ETOT_ISQ_WK = 0. 
           y = 1./SQRT(WK(IP,:))
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             do is = 2, msc
               ETOT_ISQ_WK = ETOT_ISQ_WK + 0.5*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
             end do
           end do
           !tmp = 1./SQRT(WK(IP,:))
           !ETOT_ISQ_WK = DINTSPEC_Y(IP,ACLOC,tmp)
           ETOT_SQ_WK = 0.
           y = SQRT(WK(IP,:)) 
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             do is = 2, msc
               ETOT_SQ_WK = ETOT_SQ_WK + 0.5*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
             end do
           end do
           !tmp = SQRT(WK(IP,:))
           !ETOT_SQ_WK  = DINTSPEC_Y(IP,ACLOC,tmp) 
!
! tail factors ...
!
           DS          = SPSIG(MSC) - SPSIG(MSC-1)

           ATAIL       = SUM(ACLOC(MSC,:)) * SIGPOW(MSC,1) * DDIR * DS
           ETAIL       = SUM(ACLOC(MSC,:)) * SIGPOW(MSC,2) * DDIR * DS
           ESIGTAIL    = SUM(ACLOC(MSC,:)) * SIGPOW(MSC,3) * DDIR * DS
!
! tail factors ... borowed from SWAN
!
           ACTOT       = ACTOT        + PTAIL(5)  * ATAIL 
           ETOT        = ETOT         + PTAIL(6)  * ETAIL 
           ETOT_SPSIG  = ETOT_SPSIG   + PTAIL(7)  * ETAIL 
           ETOT_ISQ_WK = ETOT_ISQ_WK  + PTAIL(5)  * ETAIL / (SQRT(WK(IP,MSC)))
           ETOT_SQ_WK  = ETOT_SQ_WK   + PTAIL(5)  * ETAIL * (SQRT(WK(IP,MSC)))
           ETOT_WK     = ETOT_WK      + PTAIL(8)  * ETAIL * WK(IP,MSC)
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
           KME01       = 10. 
           KMWAM       = 10. 
           KMWAM2      = 10. 

         end if 

       RETURN
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
