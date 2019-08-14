!**********************************************************************
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
         if (etot .gt. small) then
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
           HS          = MAX(THR,4.*SQRT(ETOT))
           SME01       = ETOT_SPSIG / ETOT
           KME01       = ETOT_WK / ETOT
           SME10       = ETOT / ACTOT
           KMWAM       = (ETOT/ETOT_ISQ_WK)**2.0
           KMWAM2      = (ETOT_SQ_WK/ETOT)**2.0

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
           do is = 2, msc
             ETOT = ETOT + 0.5*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
           end do
         end do
!
! if etot too small skip ...
!
         if (etot .gt. small) then
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
           ETOT_SPSIG = 0.
           y = SIGPOW(:,1) 
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             do is = 2, msc
               ETOT_SPSIG = ETOT_SPSIG + 0.5*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
             end do
           end do
           !tmp = SIGPOW(:,1)
           ETOT_WK = 0.
           y = WKLOC(:) 
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             do is = 2, msc
               ETOT_WK = ETOT_WK + 0.5*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
             end do
           end do
           ETOT_ISQ_WK = 0. 
           y = 1./SQRT(WKLOC(:)) 
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             do is = 2, msc
               ETOT_ISQ_WK = ETOT_ISQ_WK + 0.5*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
             end do
           end do
           ETOT_SQ_WK = 0.
           y = SQRT(WKLOC(:)) 
           do id = 1, mdc
             tmp(:) = acloc(:,id) * spsig * y
             do is = 2, msc
               ETOT_SQ_WK = ETOT_SQ_WK + 0.5*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
             end do
           end do
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
           ETOT_ISQ_WK = ETOT_ISQ_WK  + PTAIL(5)  * ETAIL / SQRT(WKLOC(MSC))
           ETOT_SQ_WK  = ETOT_SQ_WK   + PTAIL(5)  * ETAIL * SQRT(WKLOC(MSC))
           ETOT_WK     = ETOT_WK      + PTAIL(8)  * ETAIL * WKLOC(MSC) 
!
! integral parameters ...
!
           HS          = MAX(THR,4.*SQRT(ETOT))
           SME01       = ETOT_SPSIG / ETOT
           KME01       = ETOT_WK / ETOT
           SME10       = ETOT / ACTOT
           KMWAM       = (ETOT/ETOT_ISQ_WK)**2.0
           KMWAM2      = (ETOT_SQ_WK/ETOT)**2.0

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
         ETOT_SKD    = 0.
         ETOT_SKDSIG = 0.

         y = 1./SINH(MIN(20.,WK(IP,:)*DEP(IP)))**2 
         do id = 1, mdc
           tmp(:) = acloc(:,id) * spsig * y
           do is = 2, msc
             ETOT_SKD = ETOT_SKD + 0.5*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
           end do
         end do
         !tmp         = 1./SINH(MIN(30.,WK(IP,:)*DEP(IP)))**2
         !ETOT_SKD    = DINTSPEC_Y(IP, ACLOC, tmp)
         y =  SIGPOW(:,2)/SINH(MIN(20.,WK(IP,:)*DEP(IP)))**2 
         do id = 1, mdc
           tmp(:) = acloc(:,id) * spsig * y
           do is = 2, msc
             ETOT_SKDSIG = ETOT_SKDSIG + 0.5*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
           end do
         end do
         !tmp         = SIGPOW(:,2)/SINH(MIN(30.,WK(IP,:)*DEP(IP)))**2
         !ETOT_SKDSIG = DINTSPEC_Y(IP, ACLOC, tmp)

         IF (ETOT_SKD .gt. SMALL) THEN 
!
! integral parameters ...
!
           UBOT        = SQRT(ETOT_SKD)
           ORBITAL     = SQRT(2.*ETOT_SKD)
           BOTEXPER    = SQRT(2.*ETOT_SKDSIG)
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
         ETOT_SKD    = 0.
         ETOT_SKDSIG = 0.

         y = 1./SINH(MIN(20.,WKLOC*DEPLOC))**2

         do id = 1, mdc
           tmp(:) = acloc(:,id) * spsig**4. * y 
           do is = 2, msc
             ETOT_SKD = ETOT_SKD + 0.5*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
           end do
         end do

         y =  SIGPOW(:,2)/SINH(MIN(20.,WKLOC(:)*DEPLOC))**2 

         do id = 1, mdc
           tmp(:) = acloc(:,id) * spsig * y
           do is = 2, msc
             ETOT_SKDSIG = ETOT_SKDSIG + 0.5*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
           end do
         end do

         IF (ETOT_SKD .gt. SMALL) THEN 
!
! integral parameters ...
!
           UBOT        = SQRT(ETOT_SKD)
           ORBITAL     = SQRT(2.*ETOT_SKD)
           BOTEXPER    = SQRT(2.*ETOT_SKDSIG)
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
      SUBROUTINE URSELL_NUMBER(HS,TP,DEPTH,URSELL)
         USE DATAPOOL, ONLY : G9, DMIN, SMALL

         IMPLICIT NONE

         REAL, INTENT(IN)    :: HS, TP, DEPTH
         REAL, INTENT(OUT)   :: URSELL 

         IF (DEPTH .GT. DMIN .AND. HS .GT. SMALL .AND. TP .GT. SMALL) THEN
           URSELL = (G9 * HS) / (2.*SQRT(2.)*TP**2*DEPTH**2) 
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

         IF (ETOT .LT. THR) GOTO 101

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
         IF(ETOTF4 .GT. THR) THEN
            FP = ETOTF3/ETOTF4
            TP_W = 1./FP/PI2
            !CALL WAVEKCG(DEP(IP), FP, WN_W, CP_W, KP_W, CGP_W)
            CALL ALL_FROM_TABLE(FP,DEP(IP),KP_W,CGP_W,WKDEP_W,WN_W,CP_W)
            LP_W  = PI2/KP_w
         ELSE
            FP  = SMALL
            CP_W  = SMALL
            KP_W  = SMALL
            CGP_W = SMALL
            LP_W = SMALL
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
      SUBROUTINE PEAK_PARAMETER(IP,ACLOC,ISMAX,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)    :: IP, ISMAX
         REAL,    INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL  , INTENT(OUT)    :: KPP, FPP, CPP, WNPP, CGPP, TPP, LPP, PEAKDSPR,PEAKDM,DPEAK

         INTEGER                :: IS, ID, IDIRM
         REAL                   :: HQUOT, HQUOTP, ETOTF3, ETOTF4, ETOTC4, ETOTS4, PEAKFF
         REAL                   :: FF, DEG, VEC2DEG, MAXAC, E1, E2, DS, EAD, ETOTT, WKDEPP
!
! Peak period continues version... Taken from Thesis Henriques Alves ... correct citation is given there ... :)
!
       MAXAC = MAXVAL(ACLOC)

       IF (MAXAC .gt. SMALL) THEN
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
         IF(ETOTF4 .GT. SMALL) THEN
           FPP    = ETOTF3/ETOTF4
           !CALL WAVEKCG(DEP(IP), FPP, WNPP, CPP, KPP, CGPP)
           CALL ALL_FROM_TABLE(FPP,DEP(IP),KPP,CGPP,WKDEPP,WNPP,CPP)
           TPP    = 1./(FPP/PI2)
           LPP    = 1./KPP*PI2
           PEAKDM = VEC2DEG (ETOTC4, ETOTS4)
           CALL DEG2NAUT(PEAKDM,DEG,LNAUTOUT)
           PEAKDM = DEG
           PEAKFF = MIN (1.,SQRT(ETOTC4*ETOTC4+ETOTS4*ETOTS4)/ETOTF4)
           PEAKDSPR = SQRT(2.-2.*PEAKFF) * 180./PI
         ELSE
           TPP  = 0.
           FPP  = 0.
           CPP  = 0.
           KPP  = 10.
           CGPP = 0.
           LPP  = 0.
           FF       = 0.
           PEAKDSPR = 0.
           PEAKDM   = 0.
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
         TPP  = 0. 
         FPP  = 0. 
         CPP  = 0. 
         KPP  = 10.
         CGPP = 0.
         LPP  = 0.
         FF       = 0.
         PEAKDSPR = 0.
         PEAKDM   = 0.
       END IF

       RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE PEAK_PARAMETER_LOC(ACLOC,DEPLOC,ISMAX,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK)

         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)    :: ISMAX
         REAL, INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL, INTENT(IN)    :: DEPLOC


         REAL  , INTENT(OUT)    :: KPP, FPP, CPP, WNPP, CGPP, TPP, LPP, PEAKDSPR,PEAKDM,DPEAK

         INTEGER                :: IS, ID, IDIRM
         REAL                   :: HQUOT, HQUOTP, ETOTF3, ETOTF4, ETOTC4, ETOTS4, PEAKFF, WKDEPP
         REAL                   :: FF, DEG, VEC2DEG, MAXAC, E1, E2, DS, EAD, ETOTT
!
! Peak period continues version... Taken from Thesis Henriques Alves ... correct citation is given there ... :)
!
       MAXAC = MAXVAL(ACLOC)

       IF (SUM(ACLOC) .gt. SMALL) THEN
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
         IF(ETOTF4 .GT. SMALL) THEN
           FPP    = ETOTF3/ETOTF4
           !CALL WAVEKCG(DEPLOC, FPP, WNPP, CPP, KPP, CGPP)
           CALL ALL_FROM_TABLE(FPP,DEPLOC,KPP,CGPP,WKDEPP,WNPP,CPP) 
           TPP    = 1./(FPP/PI2)
           LPP    = 1./KPP*PI2
           PEAKDM = VEC2DEG (ETOTC4, ETOTS4)
           CALL DEG2NAUT(PEAKDM,DEG,LNAUTOUT)
           PEAKDM = DEG
           PEAKFF = MIN(1.,SQRT(ETOTC4**2.+ETOTS4**2.)/ETOTF3)
           PEAKDSPR = -999.!SQRT(2.-2.*PEAKFF) * 180./PI
!AR: Needs some fix it gives always 0. back since PEAKFF is always 1 ... 
         ELSE
           TPP  = 0.
           FPP  = 0.
           CPP  = 0.
           KPP  = 10.
           CGPP = 0.
           LPP  = 0.
           FF       = 0.
           PEAKDSPR = 0.
           PEAKDM   = 0.
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
           DPEAK    = 0.
         END IF
       ELSE
         TPP  = 0. 
         FPP  = 0. 
         CPP  = 0. 
         KPP  = 10.
         CGPP = 0.
         LPP  = 0.
         FF       = 0.
         PEAKDSPR = 0.
         PEAKDM   = 0.
       END IF

       !WRITE(*,'(2F15.4,I10,10F8.2)') SUM(ACLOC),DEPLOC,ISMAX,FPP,TPP,CPP,WNPP,CGPP,KPP,LPP,PEAKDSPR,PEAKDM,DPEAK

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

         IF (ETOT1 .GT. SMALL ) THEN
           DM    = VEC2DEG (ETOTC, ETOTS)
           CALL DEG2NAUT(DM,DEG,LNAUTOUT)
           DM = DEG
           FF = MIN (1., SQRT(ETOTC*ETOTC+ETOTS*ETOTS)/ETOT1)
           DSPR = SQRT(2.-2.*FF) * 180./PI
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

         IF (ETOT1 .GT. SMALL ) THEN
           DM    = VEC2DEG (ETOTC, ETOTS)
           CALL DEG2NAUT(DM,DEG,LNAUTOUT)
           DM = DEG
           FF = MIN (1., SQRT(ETOTC*ETOTC+ETOTS*ETOTS)/ETOT1)
           DSPR = SQRT(2.-2.*FF) * 180./PI
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
      SUBROUTINE MEAN_WAVE_LENGTH(IP,ACLOC,ISMAX,WLM,KLM)
         USE DATAPOOL
         IMPLICIT NONE

         INTEGER, INTENT(IN)    :: IP, ISMAX
         REAL,    INTENT(IN)    :: ACLOC(MSC,MDC)
         REAL,    INTENT(OUT)   :: WLM, KLM

         INTEGER                :: IS, ID

         REAL                   :: SIG2, SKK, PPTAIL, CETAIL, CKTAIL
         REAL                   :: ETOT1, EKTOT, STEEPNESS

         ETOT1 = 0.
         EKTOT = 0.

         DO IS=1, MSC
           SIG2 = SIGPOW(IS,2)
           SKK  = SIG2 * WK(IP,IS)
           DO ID=1,MDC
             ETOT1  = ETOT1 + SIG2 * ACLOC(IS,ID)
             EKTOT  = EKTOT +   SKK * ACLOC(IS,ID)
           ENDDO
         ENDDO

         ETOT1  = FRINTF * ETOT1
         EKTOT  = FRINTF * EKTOT

         IF (MSC .GT. 3) THEN
            PPTAIL = PTAIL(1) - 1.
            CETAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
            PPTAIL = PTAIL(1) - 1. - 2.*1.
            CKTAIL = 1. / (PPTAIL * (1. + PPTAIL * (FRINTH-1.)))
            DO ID=1,MDC
              ETOT1 = ETOT1 + CETAIL * SIG2  * ACLOC(MSC,ID)
              EKTOT = EKTOT + CKTAIL * SKK   * ACLOC(MSC,ID)
            ENDDO
         ENDIF

         IF (EKTOT.GT.SMALL) THEN
            WLM = PI2 * (ETOT1/EKTOT)
            KLM = PI2 / WLM 
            STEEPNESS = 4.* SQRT(ETOT1*DDIR)/WLM
         ELSE
            WLM = 0.
            KLM = 0.
         ENDIF

       RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MEAN_PARAMETER(IP,ACLOC,ISMAX,HS,TM01,TM02,KLM,WLM)

         USE DATAPOOL
         IMPLICIT NONE

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
         ETOT = 0.
         do id = 1, mdc
           actmp(:) = acloc(:,id)
           tmp      = actmp * spsig
           do is = 2, msc
             ETOT = ETOT + 0.5*(tmp(is)+tmp(is-1))*ds_incr(is)*ddir
           end do
         end do
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

         IF (EPTOT .GT. THR) THEN
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
         IF (EFTOT .GT. SMALL .AND. ETOT .GT. SMALL) THEN
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

         IF (ETOT1.GT.THR.AND.EKTOT.GT.THR) THEN
            WLM = PI2 * (ETOT1/EKTOT)
            KLM  = PI2/WLM
         ELSE
            KLM  = 10.
            WLM = 0.
         ENDIF

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

         IF (EPTOT .GT. THR) THEN
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
         IF (EFTOT .GT. SMALL .AND. ETOT .GT. SMALL) THEN
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

         IF (ETOT1.GT.THR.AND.EKTOT.GT.THR) THEN
            WLM = PI2 * (ETOT1/EKTOT)
            KLM  = PI2/WLM
         ELSE
            KLM  = 10.
            WLM = 0.
         ENDIF

       RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
