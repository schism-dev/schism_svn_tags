!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DIFFRA_EXTENDED()
         USE DATAPOOL
#ifdef SELFE
        use elfe_msgp
#endif
         IMPLICIT NONE

         INTEGER :: IP, IS, ID

         REAL :: ETOT, EWKTOT, EWCTOT, ECGTOT, EAD
         REAL :: DFWAV
         REAL :: AUX
         REAL :: DELTA

         REAL :: EWK(MNP), EWC(MNP), ECG(MNP), ENG(MNP)
         REAL :: CCG(MNP)
         REAL :: DXENG(MNP), DYENG(MNP), DXXEN(MNP), DYYEN(MNP), DXYEN(MNP)
         REAL :: DXCCG(MNP), DYCCG(MNP)
         REAL :: DFCUR(MNP)
         REAL :: DFBOT(MNP)

         REAL :: DFCUT
         REAL :: ETOTC, ETOTS, DM
         REAL :: US
         REAL :: CAUX, CAUX2
         REAL :: NAUX, AUX1, AUX2

         REAL :: T1, T2, BETA
         
         DFBOT(:) = 0.0
         DFCUR(:) = 0.0
         DIFRM    = 0.0

         EWK = 0.
         EWC = 0.
         ECG = 0.
         ENG = 0.
         CCG = 0.
         
         DO IP = 1, MNP
           IF (DEP(IP) .LT. DMIN) THEN
             EWK(IP) = 10.
             EWC(IP) = 0.
             ECG(IP) = 0.
             ENG(IP) = 0.
             CCG(IP) = 0.
           ELSE
             ETOT = 0.0
             EWKTOT = 0.0
             EWCTOT = 0.0
             ECGTOT = 0.0
             DO IS = 1, MSC
               EAD = SUM(AC2(IP,IS,:))*DDIR*SIGPOW(IS,2)
               ETOT = ETOT + EAD
               EWKTOT = EWKTOT + WK(IP,IS) * EAD
               EWCTOT = EWCTOT + SPSIG(IS)/WK(IP,IS) * EAD 
               ECGTOT = ECGTOT + CG(IP,IS) * EAD
             END DO
             IF (ETOT .GT. SMALL) THEN 
               ETOT    = FRINTF * ETOT
               EWKTOT  = FRINTF * EWKTOT
               EWCTOT  = FRINTF * EWCTOT
               ECGTOT  = FRINTF * ECGTOT
               EWK(IP) = EWKTOT / ETOT
               EWC(IP) = EWCTOT / ETOT
               ECG(IP) = ECGTOT / ETOT
               ENG(IP) = SQRT(ETOT)
               CCG(IP) = EWC(IP) * ECG(IP)
             ELSE  
               EWK(IP) = 10. 
               EWC(IP) = 0. 
               ECG(IP) = 0. 
               ENG(IP) = 0. 
               CCG(IP) = 0. 
             END IF 
           END IF
         END DO

         CALL DIFFERENTIATE_XYDIR(   ENG, DXENG, DYENG )
         CALL DIFFERENTIATE_XYDIR( DXENG, DXXEN, DXYEN )
         CALL DIFFERENTIATE_XYDIR( DYENG, DXYEN, DYYEN )
         CALL DIFFERENTIATE_XYDIR(   CCG, DXCCG, DYCCG )

         CALL BOTEFCT( EWK, DFBOT ) 
        
         IF (LSTCU .OR. LSECU) CALL CUREFCT( MNP, DXENG, DYENG, CURTXY, DFCUR )

         DO IP = 1, MNP
           DIFRM(IP) = 1.0
           IF (ENG(IP) .GT. SMALL .AND. DEP(IP) .GT. DMIN .AND. IOBP(IP) .NE. 0) THEN 
             AUX = CCG(IP)*EWK(IP)*EWK(IP)
             BETA = CCG(IP)*EWK(IP)*EWK(IP)*ENG(IP)
             DFWAV = ( DXCCG(IP)*DXENG(IP)+DYCCG(IP)*DYENG(IP)+CCG(IP)*(DXXEN(IP)+DYYEN(IP)) ) / MAX(SMALL,ENG(IP))
             NAUX = ECG(IP) / MAX(SMALL,EWC(IP))
             IF (LSTCU .OR. LSECU) THEN
               ETOTC = 0.0
               ETOTS = 0.0
               DO ID = 1, MDC
                 EAD = SUM(AC2(IP,:,ID)*SIGPOW(:,2))*FRINTF*DDIR
                 ETOTC = ETOTC + EAD*COSTH(ID)
                 ETOTS = ETOTS + EAD*SINTH(ID)
               END DO
               DM = ATAN2(ETOTS,ETOTC)
               US = CURTXY(IP,1)*COS(DM)+CURTXY(IP,2)*SIN(DM)
               CAUX = US / MAX(SMALL,EWC(IP))
               DFCUT = (2.0/MAX(SMALL,NAUX)+NAUX*CAUX)*CAUX
             ELSE
               DFCUT = 0.0
               CAUX = 0.0
             END IF ! LSTCU .OR. LSECU
             CAUX2 = CAUX * CAUX
             DELTA = CAUX2*(1.0+CAUX)**2.0-NAUX*(CAUX2-NAUX)*(1.0+(DFWAV+DFBOT(IP)+DFCUR(IP))/MAX(SMALL,AUX)+DFCUT) 
             IF (BETA > SMALL) THEN
               IF (ENG(IP) > SMALL) THEN
                  AUX1 = DXCCG(IP)*DXENG(IP)+DYCCG(IP)*DYENG(IP)
                  AUX2 = CCG(IP)*(DXXEN(IP)+DYYEN(IP))
               ELSE
                  AUX1 = 0.
                  AUX2 = 0.
               END IF
               T1 = AUX1 + AUX2
               T2 = ENG(IP) * G9 * DFBOT(IP)
               DELTA = (T1/MAX(SMALL,BETA)) + (T2/MAX(SMALL,BETA))
               IF (DELTA < -1.) THEN
                 DIFRM(IP) = 1.
               ELSE
                 DIFRM(IP) = SQRT(1.+DELTA)
               END IF
             ELSE ! 
               DIFRM(IP) = 1.
             END IF
           ELSE
             DIFRM(IP) = 1.
           END IF ! ENG(IP) .GT. SMALL .AND. DEP(IP) .GT. DMIN .AND. IOBP(IP) .NE. 0
           IF (DIFRM(IP) .GT. 1.2) DIFRM(IP) = 1.2
           IF (DIFRM(IP) .LT. 0.8) DIFRM(IP) = 0.8
         END DO ! MNP

         CALL DIFFERENTIATE_XYDIR(DIFRM, DIFRX, DIFRY)       

         IF (.FALSE.) THEN
           OPEN(555, FILE  = 'ergdiffr.bin'  , FORM = 'UNFORMATTED')
           WRITE(555) RTIME
           WRITE(555)  (DIFRX(IP),DIFRY(IP),DIFRM(IP), IP = 1, MNP)
         END IF

         RETURN
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DIFFRA_SIMPLE()
         USE DATAPOOL
#ifdef SELFE
         use elfe_msgp
#endif
         IMPLICIT NONE

         INTEGER :: IP, IS
         REAL :: WVC, WVK, WVCG
         REAL :: ETOT, EWKTOT, ECGTOT, EAD
         REAL :: TMP, AUX, DELTA, TRANS_X(MNP), TRANS_Y(MNP)
         REAL :: EWK(MNP), ECG(MNP), ENG(MNP)
         REAL :: CGK(MNP)
         REAL :: DXENG(MNP), DYENG(MNP), DXXEN(MNP), DYYEN(MNP), DXYEN(MNP)
         REAL :: DXCGK(MNP), DYCGK(MNP)

         EWK = 0.
         ECG = 0.
         ENG = 0.
         CGK = 0.
         DO IP = 1, MNP
           IF (IOBP(IP) .EQ. 0 .AND. DEP(IP) .GT. DMIN) THEN 
             ETOT = 0.0
             EWKTOT = 0.0
             ECGTOT = 0.0
             DO IS = 1, MSC
               EAD = SUM(AC2(IP,IS,:))*DDIR*SIGPOW(IS,2)
               ETOT = ETOT + EAD
               EWKTOT = EWKTOT + WK(IP,IS) * EAD
               ECGTOT = ECGTOT + CG(IP,IS) * EAD
             END DO
             ETOT   = FRINTF * ETOT
             EWKTOT = FRINTF * EWKTOT
             ECGTOT = FRINTF * ECGTOT
             IF (ETOT .GT. SMALL) THEN
               EWK(IP) = EWKTOT / ETOT
               ECG(IP) = ECGTOT / ETOT
               ENG(IP) = SQRT(MAX(ETOT,SMALL))
             ELSE
               EWK(IP) = 0.
               ECG(IP) = 0. 
               ENG(IP) = 0.
             END IF 
             IF (EWK(IP) .GT. SMALL) THEN
               CGK(IP) = ECG(IP) / EWK(IP)
             ELSE
               CGK(IP) = 0. 
             END IF 
           ELSE
             EWK(IP) = 0.
             ECG(IP) = 0.
             ENG(IP) = 0.
             CGK(IP) = 0.
           END IF
         END DO

         CALL DIFFERENTIATE_XYDIR(ENG  , DXENG, DYENG)
         CALL DIFFERENTIATE_XYDIR(DXENG, DXXEN, DXYEN)
         CALL DIFFERENTIATE_XYDIR(DYENG, DXYEN, DYYEN)         
         CALL DIFFERENTIATE_XYDIR(CGK  , DXCGK, DYCGK)

         IF (LSPHE) THEN
           TRANS_X = 1.d0/(DEGRAD*REARTH*COS(YP*DEGRAD))
           TRANS_Y = 1.d0/(DEGRAD*REARTH) 
           DXENG = DXENG * TRANS_X 
           DYENG = DYENG * TRANS_Y 
           DXXEN = DXXEN * TRANS_X**2. 
           DYYEN = DYYEN * TRANS_Y**2.
           DXCGK = DXCGK * TRANS_X 
           DYCGK = DYCGK * TRANS_Y 
         END IF
         
         DO IP = 1, MNP
            IF (IOBP(IP) .EQ. 0 .AND. DEP(IP) .GT. DMIN .AND. ENG(IP) .GT. SMALL) THEN
              TMP = ECG(IP)*EWK(IP)*ENG(IP)
              IF (TMP > SMALL) THEN
                 AUX = DXCGK(IP)*DXENG(IP)+DYCGK(IP)*DYENG(IP)+CGK(IP)*(DXXEN(IP)+DYYEN(IP))
                 TMP = AUX/TMP
              ELSE
                 TMP = 0.0
              END IF
              IF (TMP < -1.0) THEN
                 DIFRM(IP) = 1.0
              ELSE
                 DIFRM(IP) = SQRT(1.0+TMP)
              END IF
            ELSE
              DIFRM(IP) = 1.0
            END IF
            IF (DIFRM(IP) .GT. 1.2) DIFRM(IP) = 1.2
            IF (DIFRM(IP) .LT. 0.8) DIFRM(IP) = 0.8
         END DO

         CALL DIFFERENTIATE_XYDIR(DIFRM, DIFRX, DIFRY)

         IF (LSPHE) THEN
           DIFRX = DIFRX * TRANS_X 
           DIFRY = DIFRY * TRANS_Y 
         END IF

         IF (.TRUE.) THEN
           OPEN(555, FILE  = 'ergdiffr.bin'  , FORM = 'UNFORMATTED')
           WRITE(555) RTIME
           WRITE(555)  (DIFRX(IP), DIFRY(IP),DIFRM(IP)-1., IP = 1, MNP)
         END IF

         !WRITE(WWMDBG%FHNDL,*) MAXVAL(DIFRM), MAXVAL(DIFRX), MAXVAL(DIFRY)
         !WRITE(WWMDBG%FHNDL,*) MINVAL(DIFRM), MINVAL(DIFRX), MINVAL(DIFRY)
         
         RETURN
      END SUBROUTINE 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE SMOOTH( BETA, MNP, XP, YP, VAR )
         IMPLICIT NONE

         INTEGER :: MNP
         REAL, INTENT(IN)    :: XP(MNP), YP(MNP)
         REAL, INTENT(INOUT) :: VAR(MNP)
         REAL, INTENT(IN)    :: BETA
         REAL :: VART(MNP)
         REAL :: SW, SWQ, DISX, DISY, DIST, DIS
         INTEGER :: I, J, IP, IC
         
         DO I = 1, MNP
            SW = 0.0
            SWQ = 0.0
            DO J = 1, MNP               
               DISX = (XP(I) - XP(J))**2.0
               DISY = (YP(I) - YP(J))**2.0
               DIST = DISX + DISY
               IF (DIST > TINY(1.)) THEN
                  DIS = SQRT(DIST)**BETA
               ELSE
                  DIS = 1.0
               END IF 
               SW = SW + DIS
               SWQ = SWQ + DIS*VAR(J)
            END DO
            VART(I) = SWQ / SW
         END DO
         VAR(:) = VART(:)

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
     SUBROUTINE BOTEFCT( EWK, DFBOT )
        USE DATAPOOL
#ifdef SELFE
        use elfe_glbl, only: errmsg
        use elfe_msgp
#endif
        IMPLICIT NONE  

        REAL, INTENT(IN)    :: EWK(MNP)
        REAL, INTENT(INOUT) :: DFBOT(MNP)

        REAL :: SLPH(MNP), CURH(MNP)
        REAL :: DXDEP(MNP) , DYDEP(MNP) 
        REAL :: DXXDEP(MNP), DXYDEP(MNP), DYYDEP(MNP)

        REAL :: KH, BOTFC, BOTFS

        INTEGER :: IP

        LOGICAL :: ISNAN, ISINF

!        CALL SMOOTH( -1.1, MNP, XP, YP, DEP )  

        CALL DIFFERENTIATE_XYDIR(DEP, DXDEP ,  DYDEP)
        CALL DIFFERENTIATE_XYDIR(DXDEP , DXXDEP, DXYDEP)
        CALL DIFFERENTIATE_XYDIR(DYDEP , DXYDEP, DYYDEP)

        SLPH = DXDEP**2.0 + DYDEP**2.0
        CURH = DXXDEP + DYYDEP
        DFBOT = 0. 

        DO IP = 1, MNP
          IF (EWK(IP) < SMALL) CYCLE  
          KH = EWK(IP)*DEP(IP)
          IF (KH > PI) CYCLE
          DFBOT(IP) = (BOTFC(KH)*CURH(IP)+BOTFS(KH)*EWK(IP)*SLPH(IP))*G9
          IF (ISINF(DFBOT(IP)) .OR. ISNAN(DFBOT(IP))) THEN
#ifdef SELFE
            WRITE(errmsg,*)'DFBOT is NaN', IP,KH, CURH(IP), BOTFS(KH), BOTFC(KH), EWK(IP), SLPH(IP)
            call parallel_abort(errmsg)
#else
            WRITE(*,*) 'DFBOT is NaN', IP
            WRITE(*,*) KH, CURH(IP), BOTFS(KH), BOTFC(KH), EWK(IP), SLPH(IP)
            STOP 'DFBOT'
#endif
          END IF
        END DO

        RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      REAL FUNCTION BOTFC(KH)
        USE DATAPOOL
#ifdef SELFE
        use elfe_glbl, only: errmsg
        use elfe_msgp
#endif
        IMPLICIT NONE

        REAL, INTENT(IN) :: KH
        REAL*8 :: DKH, AUX, AUX1
        REAL*8 :: COSHKH, COSH2KH, SINHKH, SINH2KH, SINH3KH

        DKH = DBLE(KH)

        SINHKH  = DSINH(MIN(300.d0,DKH))
        SINH2KH = DSINH(MIN(300.d0,2.d0*DKH))
        SINH3KH = DSINH(MIN(300.d0,3.d0*DKH))
        COSHKH  = DCOSH(MIN(300.d0,DKH))
        COSH2KH = DCOSH(MIN(300.d0,2.d0*DKH))        
        AUX = -4.*KH*COSHKH+SINH3KH+SINHKH+8.*(KH**2)*SINHKH
        AUX1 = 8.*COSHKH**3*(2.*KH+SINH2KH)
        BOTFC = REAL(AUX/MAX(DBLE(SMALL),AUX1) - KH*TANH(KH)/MAX(DBLE(SMALL),(2.*(COSHKH)**2)))

        IF (BOTFC .NE. BOTFC) THEN
#ifdef SELFE
          WRITE(errmsg,*)'BOTFC is NaN Aron', KH, AUX, AUX1, TANH(KH)
          call parallel_abort(errmsg)
#else
          WRITE(*,*) 'BOTFC is NaN', KH, AUX, AUX1, TANH(KH)
          WRITE(*,*) SINHKH, SINH2KH, SINH3KH, COSHKH, COSH2KH
          WRITE(*,*) AUX, AUX1, BOTFC, SINH(30.), COSH(30.)
          STOP 'BOTFC'
#endif
        END IF
        
        RETURN
      END FUNCTION
!**********************************************************************
!*                                                                    *
!**********************************************************************
      REAL FUNCTION BOTFS(KH)
        USE DATAPOOL, ONLY : THR8, SMALL
#ifdef SELFE
        use elfe_glbl, only: errmsg
        use elfe_msgp
#endif

        IMPLICIT NONE

        REAL, INTENT(IN)  :: KH
        REAL*8            :: SECH, SINH2KH, SINHKH, COSH2KH
        REAL*8            :: AUX, AUX1

        REAL*8            :: DKH
        
        DKH = DBLE(KH)

        IF (ABS(DCOSH(DKH)) > THR8) THEN
          SECH = 1. / DCOSH(DKH)
        ELSE
          SECH = 0. 
        END IF

        SINHKH  = DSINH(MIN(300.d0,DKH))
        SINH2KH = DSINH(MIN(300.d0,2.d0*DKH))
        COSH2KH = DCOSH(MIN(300.d0,2.*DKH))
        AUX     = SECH**2/MAX(THR8,(6.d0*(2.d0*DKH+SINH2KH)**3))
        AUX1    = 8.d0*(DKH**4.d0)+16.d0*(DKH**3.d0)*SINH2KH- &
              &   9.d0*(SINH2KH**2.d0*COSH2KH+12.d0*DKH*(1.d0+2.d0*(SINHKH)**4.d0)* &
              & (DKH+SINH2KH))
        BOTFS = REAL(AUX * AUX1)

        IF (BOTFS .NE.  BOTFS) THEN
#ifdef SELFE
          WRITE(errmsg,*)'BOTFS is NaN', BOTFS, AUX, AUX1
          call parallel_abort(errmsg)
#else
          WRITE(*,*) 'BOTFS is NaN', BOTFS, AUX, AUX1
          WRITE(*,*) KH
          WRITE(*,*) SECH, SINH2KH, SINH2KH**2.,COSH2KH,SINHKH
          WRITE(*,*) 9.d0*(SINH2KH**2.d0*COSH2KH+12.d0*KH*(1.d0+2.d0*(SINHKH)**4.d0)*(KH+SINH2KH))
          STOP 'BOTFS'
#endif
        ENDIF

        RETURN
      END FUNCTION      
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CUREFCT( MNP, DXENG, DYENG, CURT, DFCUR )     
         IMPLICIT NONE

         INTEGER, INTENT(IN) :: MNP
         REAL, INTENT(IN)    :: DXENG(MNP), DYENG(MNP), CURT(MNP,2)
         REAL, INTENT(INOUT) :: DFCUR(MNP)

         REAL                :: AUX(MNP), AUXX(MNP), AUXY(MNP)

         REAL                :: DXAUXX(MNP), DYAUXY(MNP)
         
         AUX(:) = CURT(:,1) * DXENG(:) + CURT(:,2) * DYENG(:)
         AUXX(:) = AUX(:) * CURT(:,1)
         AUXY(:) = AUX(:) * CURT(:,2)        
         
         CALL DIFFERENTIATE_XYDIR(AUXX, DXAUXX, AUX)
         CALL DIFFERENTIATE_XYDIR(AUXY, DYAUXY, AUX)
         
         DFCUR(:) = DXAUXX(:) + DYAUXY(:)
         
         RETURN
      END SUBROUTINE     
!**********************************************************************
!*                                                                    *
!**********************************************************************
      REAL FUNCTION BOTFC2(KH)
        USE DATAPOOL, ONLY : SMALL, LARGE
        IMPLICIT NONE
        REAL, INTENT(IN) :: KH
        REAL :: AUX, AUX1 
        REAL :: SINHKH, COSHKH, SINH2KH, SINH3KH

        IF (KH .GT. LARGE) THEN
          BOTFC2 = 0.
          RETURN
        END IF

        COSHKH  = COSH(MIN(30.,KH))
        SINHKH  = SINH(MIN(30.,KH))
        SINH2KH = SINH(MIN(30.,2.*KH))
        SINH3KH = SINH(MIN(30.,3.*KH))

        AUX = -4.0*KH*COSHKH + SINH3KH + SINHKH + 8.0*(KH**2.0)*SINHKH
        AUX1 = 8.0*COSHKH**3.0*(2.*KH + SINH2KH)
        BOTFC2 = AUX / MAX(SMALL,AUX1) - KH*TANH(KH) / (2.0*(COSHKH)**2.0)

        IF (BOTFC2 .NE. BOTFC2) THEN
           WRITE(*,*) 'BOTFC2'
           WRITE(*,*) SINHKH, COSHKH, SINH2KH, SINH3KH, KH
           WRITE(*,*) AUX, AUX1, KH*TANH(KH), (2.0*(COSHKH)**2.0)
           STOP  'BOTFC2'
        ENDIF

        RETURN
      END FUNCTION
!**********************************************************************
!*                                                                    *
!**********************************************************************
      REAL FUNCTION BOTFS2(KH)
        USE DATAPOOL, ONLY : SMALL, LARGE
        IMPLICIT NONE
        REAL, INTENT(IN) :: KH
        REAL :: SECH
        REAL :: AUX, AUX1, SINHKH, COSHKH, SINH2KH, SINH3KH, COSH2KH

        IF (KH .GT. LARGE) THEN
          BOTFS2 = 0.
          RETURN
        END IF 

        COSHKH  = COSH(MIN(30.,KH))
        COSH2KH = COSH(MIN(30.,2*KH))
        SINHKH  = SINH(MIN(30.,KH))
        SINH2KH = SINH(MIN(30.,2.*KH))
        SINH3KH = SINH(MIN(30.,3.*KH))

        SECH = 1.0 / MAX(SMALL,COSHKH)
        AUX = SECH**2.0 / MAX(SMALL, (6.0*(2.0*KH + SINH2KH)**3.0))
        IF (AUX .GT. SMALL) THEN
          AUX1 = 8.0*(KH**4.0) + 16.0*(KH**3.0)*SINH2KH &
              - 9.0*(SINH2KH)**2.0*COSH2KH &
              + 12.0*KH*(1.0 + 2*SINHKH**4.0)*(KH + SINH2KH)
        ELSE
          AUX1 = 0.
        END IF
        BOTFS2 = AUX * AUX1
        IF (BOTFS2 .NE. BOTFS2) THEN
          WRITE(*,*) 'BOTFS2'
          WRITE(*,*) COSHKH, COSH2KH, SINHKH, SINH2KH, SINH3KH
          WRITE(*,*) AUX, AUX1, KH, SECH, SECH**2., COSHKH
          STOP 'BOTFS2'
        END IF

        RETURN
      END FUNCTION
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CUREFCT2( MNP, DXENG, DYENG, CURT, DFCUR )
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: MNP
         REAL, INTENT(IN) :: DXENG(MNP), DYENG(MNP), CURT(MNP,2)
         REAL, INTENT(INOUT) :: DFCUR(MNP)
         REAL :: AUX(MNP), AUXX(MNP), AUXY(MNP)
         REAL :: DXAUXX(MNP), DYAUXY(MNP)

         AUX(:) = CURT(:,1) * DXENG(:) + CURT(:,2) * DYENG(:)
         AUXX(:) = AUX(:) * CURT(:,1)
         AUXY(:) = AUX(:) * CURT(:,2)

         CALL DIFFERENTIATE_XYDIR(AUXX, DXAUXX, AUX)
         CALL DIFFERENTIATE_XYDIR(AUXY, DYAUXY, AUX)

         DFCUR(:) = DXAUXX(:) + DYAUXY(:)

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE BOTEFCT2( MNP, EWK, DEP, DFBOT )
        USE DATAPOOL, ONLY : LARGE
        IMPLICIT NONE
        REAL, PARAMETER :: GRAV = 9.81
        INTEGER, INTENT(IN) :: MNP
        REAL, INTENT(IN) :: EWK(MNP)
        REAL, INTENT(IN) :: DEP(MNP)
        REAL, INTENT(INOUT) :: DFBOT(MNP)
        REAL :: SLPH(MNP), CURH(MNP)
        REAL :: DXDEP(MNP), DYDEP(MNP)
        REAL :: DXXDEP(MNP), DXYDEP(MNP), DYYDEP(MNP)
        REAL :: KH
        INTEGER :: IP
        REAL :: BOTFC2, BOTFS2

        CALL DIFFERENTIATE_XYDIR(DEP(1)  , DXDEP(1) ,  DYDEP(1))
        CALL DIFFERENTIATE_XYDIR(DXDEP(1), DXXDEP(1), DXYDEP(1))
        CALL DIFFERENTIATE_XYDIR(DYDEP(1), DXYDEP(1), DYYDEP(1))

        SLPH(:) = DXDEP(:)**2.0 + DYDEP(:)**2.0
        CURH(:) = DXXDEP(:) + DYYDEP(:)

        DO IP = 1, MNP
          KH = EWK(IP)*DEP(IP)
          DFBOT(IP) = (BOTFC2(KH)*CURH(IP)+BOTFS2(KH)*EWK(IP)*SLPH(IP))*GRAV
          IF (DFBOT(IP) .NE. DFBOT(IP)) THEN
            WRITE(*,*) 'DFBOT'
            WRITE(*,*) DFBOT(IP), BOTFC2(KH), BOTFS2(KH), CURH(IP), EWK(IP), SLPH(IP)
            STOP 'DFBOT'
          ENDIF
        END DO

        RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE DIFFRA_LIAU()
         USE DATAPOOL
!         USE TRIANGULATION_BOUNDARY
         IMPLICIT NONE
         INTEGER :: IP, IS, ID, IE
         REAL :: WVC, WVK, WVCG, WVKDEP, WVN
         REAL :: ETOT, EWKTOT, EWCTOT, ECGTOT, EAD
         REAL :: DFWAV
         REAL :: AUX
         REAL :: DELTA
         REAL :: EWK(MNP), EWC(MNP), ECG(MNP), ENG(MNP)
         REAL :: CCG(MNP)
         REAL :: DXENG(MNP), DYENG(MNP), DXXEN(MNP), DYYEN(MNP), DXYEN(MNP)
         REAL :: DXCCG(MNP), DYCCG(MNP)
         REAL :: DFCUR(MNP)
         REAL :: DFBOT(MNP)
         REAL :: DFCUT
         REAL :: ETOTC, ETOTS, DM
         REAL :: US
         REAL :: CAUX, CAUX2
         REAL :: NAUX

         DFBOT(:) = 0.0
         DFCUR(:) = 0.0

         DO IP = 1, MNP
            ETOT = 0.0
            EWKTOT = 0.0
            EWCTOT = 0.0
            ECGTOT = 0.0
            IF (DEP(IP) .GT. DMIN) THEN
              DO IS = 1, MSC
                CALL ALL_FROM_TABLE(SPSIG(IS),DEP(IP),WVK,WVCG,WVKDEP,WVN,WVC)
                EAD = SUM(AC2(IP,IS,:))*DDIR*SIGPOW(IS,2)
                ETOT = ETOT + EAD
                EWKTOT = EWKTOT + WVK *EAD
                EWCTOT = EWCTOT + WVC *EAD
                ECGTOT = ECGTOT + WVCG*EAD
              END DO
              ETOT   = FRINTF * ETOT
              EWKTOT = FRINTF * EWKTOT
              EWCTOT = FRINTF * EWCTOT
              ECGTOT = FRINTF * ECGTOT
              IF (ETOT .LT. SMALL) THEN
                EWK(IP) = 0.0
                EWC(IP) = 0.0
                ECG(IP) = 0.0
                ENG(IP) = 0.0
                CCG(IP) = 0.0
              ELSE
                IF (MSC > 3) THEN
                  ETOT   = ETOT   + PTAIL(2)*EAD
                  EWKTOT = EWKTOT + PTAIL(4)*EAD
                  EWCTOT = EWCTOT + PTAIL(4)*EAD
                  ECGTOT = ECGTOT + PTAIL(4)*EAD
                END IF
                EWK(IP) = EWKTOT / ETOT
                EWC(IP) = EWCTOT / ETOT
                ECG(IP) = ECGTOT / ETOT
                ENG(IP) = SQRT(MAX(ETOT,1.0E-8))
                CCG(IP) = EWC(IP) * ECG(IP)
              END IF
            ELSE
              EWK(IP) = 0.0
              EWC(IP) = 0.0
              ECG(IP) = 0.0
              ENG(IP) = 0.0
              CCG(IP) = 0.0
           END IF
         END DO

         CALL DIFFERENTIATE_XYDIR(ENG, DXENG, DYENG)
         CALL DIFFERENTIATE_XYDIR(DXENG, DXXEN, DXYEN)
         CALL DIFFERENTIATE_XYDIR(DYENG, DXYEN, DYYEN)
         CALL DIFFERENTIATE_XYDIR(CCG, DXCCG, DYCCG)

         IF (.TRUE.) CALL BOTEFCT2( MNP, EWK, DEP, DFBOT )

         IF (LSTCU .OR. LSECU) CALL CUREFCT2( MNP, DXENG, DYENG, CURTXY, DFCUR )

         DO IP = 1, MNP
            AUX = CCG(IP)*EWK(IP)*EWK(IP)
            IF ( AUX*ENG(IP) .GT. SMALL) THEN
               DFWAV = ( DXCCG(IP)*DXENG(IP)+DYCCG(IP)*DYENG(IP)+CCG(IP)*(DXXEN(IP)+DYYEN(IP)) ) / MAX(SMALL,ENG(IP))
               NAUX = ECG(IP) / MAX(SMALL,EWC(IP))
               IF (LSTCU .OR. LSECU) THEN
                 ETOTC = 0.0
                 ETOTS = 0.0
                 DO ID = 1, MDC
                   EAD = SUM(AC2(IP,:,ID)*SIGPOW(:,2))*FRINTF*DDIR
                   ETOTC = ETOTC + EAD*COS(SPDIR(ID))
                   ETOTS = ETOTS + EAD*SIN(SPDIR(ID))
                 END DO
                 DM = ATAN2(ETOTS,ETOTC)
                 US = CURTXY(IP,1)*COS(DM)+CURTXY(IP,2)*SIN(DM)
                 CAUX = US / EWC(IP)
                 DFCUT = (2.0/NAUX+NAUX*CAUX)*CAUX
               ELSE
                 DFCUT = 0.0
                 CAUX = 0.0
               END IF
               CAUX2 = CAUX * CAUX
               DELTA = CAUX2*(1.0+CAUX)**2.0-NAUX*(CAUX2-NAUX)*(1.0+(DFWAV+DFBOT(IP)+DFCUR(IP))/AUX+DFCUT)
               IF (DELTA <= 0.0) THEN
                 DIFRM(IP) = 1.0
               ELSE
                 DIFRM(IP) = 1.0/(CAUX2-NAUX)*(CAUX*(1.0+CAUX)-SQRT(DELTA))
               END IF
            ELSE
               DIFRM(IP) = 1.0
            END IF
            IF (DIFRM(IP) .GT. 1.2) DIFRM(IP) = 1.2
            IF (DIFRM(IP) .LT. 0.8) DIFRM(IP) = 0.8
         END DO
         CALL DIFFERENTIATE_XYDIR(DIFRM, DIFRX, DIFRY)

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************

