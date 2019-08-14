!**********************************************************************
!*                                                                    *
!**********************************************************************
#ifdef WWMONLY
       SUBROUTINE RADIATION_STRESS()
        USE DATAPOOL
        IMPLICIT NONE

        INTEGER :: IP,IL,IS,ID
        REAL*8  :: ACLOC(MSC,MDC)
        REAL*8  :: COSE2, SINE2, COSI2
        REAL*8  :: D_RSXX(MNP,2), D_RSXY(MNP,2), D_RSYY(MNP,2)
        REAL*8  :: ZSXX(MNP),ZSYY(MNP),ZSXY(MNP),EWK(MNP),EWS(MNP),EWN(MNP),ETOT(MNP),MDIR(MNP)
        REAL*8  :: m0, m0d, tmp, EWKTOT, EHFR, ELOC, EFTAIL, INCRZ, Z_MIN, ZZETA, DVEC2RAD, WN
        REAL*8  :: DS, D, KW, KD, SINH2KD, SINHKW, COSH2KW, COSHKW, COSHKD, ETOTS, ETOTC, EWSIG
        REAL    :: WNTMP,WKTMP,WCGTMP,WCTMP,WKDEPTMP
        REAL    :: WSTMP, DEPLOC

        SXX3D(:,:) = 0.
        SYY3D(:,:) = 0.
        SXY3D(:,:) = 0.

        EFTAIL = 1.0 / (PTAIL(1)-1.0)

        ETOT = 0.0
        MDIR = 0.0

        IF (LETOT) THEN
!AR: Estimate zeroth moment m0, mean wave direction, dominant wave number, dominant sigma ...
          DO IP = 1, MNP
            DEPLOC = MAX(DMIN,DEP(IP))
            ACLOC = DBLE(AC2(IP,:,:))
            m0    = 0.
            EWSIG  = 0.
            ETOTS  = 0.
            ETOTC  = 0.
            IF (MSC .GE. 2) THEN
              DO ID = 1, MDC
                m0d = 0.
                DO IS = 2, MSC
                  tmp = 0.5*(SPSIG(IS)*ACLOC(IS,ID)+SPSIG(IS-1)*ACLOC(IS-1,ID))*DS_INCR(IS)*DDIR
                  m0 = m0 + tmp
                  EWSIG  = EWSIG  + SPSIG(IS) * tmp
                  m0d = m0d + tmp
                END DO
                IF (MSC > 3) THEN
                  EHFR = ACLOC(MSC,ID) * SPSIG(MSC)
                  m0 = m0 + DDIR * EHFR * SPSIG(MSC) * EFTAIL
                endif
                ETOTC  = ETOTC + m0d * COS(SPDIR(ID))
                ETOTS  = ETOTS + m0d * SIN(SPDIR(ID))
              END DO
            ELSE
              DS = SGHIG - SGLOW
              DO ID = 1, MDC
                m0d = ACLOC(1,ID) * DS * DDIR
                m0 = m0 + m0d
              END DO
            END IF
            ETOT(IP) = m0
            IF (m0 .GT. small) then
              EWS(IP) = EWSIG/m0
            ELSE
              EWS(IP) = 0.d0
              EWN(IP) = 0.d0
              MDIR(IP) = 0.d0
              ETOT(IP) = 0.d0
              CYCLE
            ENDIF
            WSTMP = REAL(EWS(IP))
            CALL ALL_FROM_TABLE(WSTMP,DEPLOC,WKTMP,WCGTMP,WKDEPTMP,WNTMP,WCTMP) 
            EWN(IP) = DBLE(WNTMP)
            EWK(IP) = DBLE(WKTMP)
            MDIR(IP) = DVEC2RAD (ETOTC, ETOTS)
          END DO !IP
        END IF !LETOT

!AR: Here comes the whole story ... 
! Etot = 1/16 * Hs² = 1/8 * Hmono² => Hs² = 2 * Hmono² => Hs = sqrt(2) * Hmono => Hmono = Hs / SQRT(2) 
! Etot = 1/16 * Hs² = 1/16 * (4 * sqrt(m0))² = m0 
! Etot = 1/8 * Hmono² ... so the problem for the analytical solution evolved because we treat the Etot from Hs and Hmono there is a factor of 2 between this!
! Or in other words for the analytical solution we impose a Hs = X[m], we integrate m0 out of it and get Etot, since this Etot is a function of Hs and not Hmono
! it needs the factor of 2 between it! This should make now things clear forever. So the question is not how we calculate the total energy the question is 
! what is defined on the boundary that means we should always recalculate the boundary in terms of Hs =  SQRT(2) * Hmono !!!
! Or saying it again in other words our boundary conditions is wrong if we impose Hmono in wwminput.nml !!!

        IF (RADFLAG == 'LON') THEN
          RSXX = 0.
          RSXY = 0.
          RSYY = 0.
          DO IP = 1, MNP
            IF (.NOT. LETOT) THEN
              ACLOC = DBLE(AC2(IP,:,:))
              DO ID = 1, MDC
                DO IS = 2, MSC
                  ELOC  = (SPSIG(IS)*ACLOC(IS,ID)+SPSIG(IS-1)*ACLOC(IS-1,ID))*DS_INCR(IS)*DDIR
                  COSE2 = COS(SPDIR(ID))**2.
                  SINE2 = SIN(SPDIR(ID))**2.
                  COSI2 = COS(SPDIR(ID)) * SIN(SPDIR(ID))
                  WN    = DBLE(CG(IP,IS) / ( SPSIG(IS)/WK(IP,IS) ))
                  RSXX(IP) = RSXX(IP) + ( WN * COSE2 + WN - .5) * ELOC   ! Units = [ 1/s + 1/s - 1/s ] * m²s = m²
                  RSXY(IP) = RSXY(IP) + ( WN * COSI2                 ) * ELOC
                  RSYY(IP) = RSYY(IP) + ( WN * SINE2 + WN - .5) * ELOC
                ENDDO
              ENDDO
            ELSE
              RSXX(IP) =  ETOT(IP) * ( EWN(IP) - 0.5 + EWN(IP) * COS(MDIR(IP))**2.)
              RSXY(IP) =  ETOT(IP) *                   EWN(IP) * COS(MDIR(IP)) * SIN(MDIR(IP))
              RSYY(IP) =  ETOT(IP) * ( EWN(IP) - 0.5 + EWN(IP) * SIN(MDIR(IP))**2.)
            END IF
          END DO

          SXX3D = 0.
          SXY3D = 0.
          SYY3D = 0.
          DO IP = 1, MNP
            IF (DEP(IP) .GT. DMIN)  THEN
              SXX3D(:,IP) = DBLE(RSXX(IP)) * G9 !ccf
              SXY3D(:,IP) = DBLE(RSXY(IP)) * G9 !ccf
              SYY3D(:,IP) = DBLE(RSYY(IP)) * G9 !ccf 
            ELSE
              SXX3D(:,IP) = 0.d0
              SXY3D(:,IP) = 0.d0
              SYY3D(:,IP) = 0.d0
            END IF
          END DO
        ELSE IF (RADFLAG == 'XIA') THEN
          IF (LETOT) THEN
            DO IP = 1, MNP
              IF (ETOT(IP) .LT. SMALL .OR. EWK(IP) .LT. SMALL) CYCLE
              IF (DEP(IP) .LE. DMIN) CYCLE
              DO IL = 1, NLEV(IP) !NLVT ccf
                ZZETA        = SHYFZETA(IL,IP) + DEP(IP)
                IF (ZZETA .LT. 0) CYCLE
                KW           = EWK(IP) * ZZETA 
                KD           = EWK(IP) * DEP(IP)
                SINH2KD      = DSINH(MIN(300.d0,2.*KD))
                COSHKD       = DCOSH(MIN(300.d0,KD))
                SINHKW       = DSINH(MIN(300.d0,KW))
                COSH2KW      = DCOSH(MIN(300.d0,2.*KW))
                COSHKW       = DCOSH(MIN(300.d0,KW))
                SXX3D(IL,IP) = ETOT(IP) * EWK(IP) / SINH2KD * (COSH2KW + 1.) * COS(MDIR(IP))**2. - &
     &                         ETOT(IP) * EWK(IP) / SINH2KD * (COSH2KW - 1.)                      - &
     &                         ETOT(IP) * SHYFZETA(IL,IP) / DEP(IP)**2.                              + &
     &                         ETOT(IP) * KW * SINHKW / ( DEP(IP) * COSHKD ) - &
     &                         ETOT(IP) / DEP(IP)  *  (1. - COSHKW / COSHKD)
              END DO
            END DO
          ELSE
            SXX3D = 0.
            DO IP = 1, MNP
              IF (DEP(IP) .LT. DMIN) CYCLE
              ACLOC = DBLE(AC2(IP,:,:))
              DO IL = 1, NLVT
                ZZETA = SHYFZETA(IL,IP) + DEP(IP)
                IF (ZZETA .LT. 0) CYCLE
                DO IS = 1, MSC
                  KW           = WK(IP,IS) * ZZETA
                  KD           = WK(IP,IS) * DEP(IP)
                  SINH2KD      = DSINH(MIN(300.d0,2.*KD))
                  COSHKD       = DCOSH(MIN(300.d0,KD))
                  SINHKW       = DSINH(MIN(300.d0,KW))
                  COSH2KW      = DCOSH(MIN(300.d0,2.*KW))
                  COSHKW       = DCOSH(MIN(300.d0,KW))
                  DO ID = 1, MDC
                    !Dimension of ELOC = m^2
                    ELOC = AC2(IP,IS,ID) * SPSIG(IS)**2. * DDIR * FRINTF * 2. ! Here is the factor 2. same as mono
                    IF (ELOC .LT. SMALL) CYCLE
                      tmp          =-ELOC * WK(IP,IS) / SINH2KD * (COSH2KW - 1.) - &
     &                               ELOC * SHYFZETA(IL,IP) / DEP(IP)**2.        + &
     &                               ELOC * KW * SINHKW / ( DEP(IP) * COSHKD )   - &
     &                               ELOC / DEP(IP) *  (1. - COSHKW / COSHKD)
                      SXX3D(IL,IP) = SXX3D(IL,IP) + &
     &                               ELOC * WK(IP,IS) / SINH2KD * (COSH2KW + 1.) * COS(SPDIR(ID))**2.+tmp
                      SYY3D(IL,IP) = SYY3D(IL,IP) + &
     &                               ELOC * WK(IP,IS) / SINH2KD * (COSH2KW + 1.) * SIN(SPDIR(ID))**2.+tmp
                      SXY3D(IL,IP) = SXY3D(IL,IP) + &
     &                ELOC * WK(IP,IS) / SINH2KD * (COSH2KW + 1.) * SIN(SPDIR(ID))*COS(SPDIR(ID))
                  END DO
                END DO
              END DO
            END DO
          END IF
        END IF

        RSXX = 0.
        DO IP = 1, MNP
          IF (DEP(IP) .LE. DMIN) CYCLE
          DO IL = 2, NLVT
            RSXX(IP) = RSXX(IP) + 0.5*( SXX3D(IL,IP)+SXX3D(IL-1,IP) ) !* INCRZ/G9  ... put the right INCR in Z
          END DO
        END DO

        END SUBROUTINE
#elif SELFE
!**********************************************************************
!*                                                                    *
!**********************************************************************
       SUBROUTINE RADIATION_STRESS

        use elfe_glbl, only: iplg,errmsg
        USE elfe_msgp !, only : myrank,parallel_abort

        USE DATAPOOL
        IMPLICIT NONE

        INTEGER :: IP,IL,IS,ID
        REAL*8  :: ACLOC(MSC,MDC)
        REAL*8  :: COSE2, SINE2, COSI2
        REAL*8  :: D_RSXX(MNP,2), D_RSXY(MNP,2), D_RSYY(MNP,2)
        REAL*8  :: ZSXX(MNP),ZSYY(MNP),ZSXY(MNP),EWK(MNP),EWS(MNP),EWN(MNP),ETOT(MNP),MDIR(MNP)
        REAL*8  :: m0, m0d, tmp, EWKTOT, EHFR, ELOC, EFTAIL, INCRZ, Z_MIN, ZZETA, DVEC2RAD 
        REAL*8  :: DS, D, KW, KD, SINH2KD, SINHKW, COSH2KW, COSHKW, COSHKD, ETOTS, ETOTC, EWSIG, S11, S22
        REAL    :: WNTMP,WKTMP,WCGTMP,WCTMP,WN,WKDEPTMP
        REAL    :: WSTMP, DEPLOC

        INTEGER :: ND1,ND2
        REAL*8  :: SINHKD,FSS(NVRT,MNP),FCS(NVRT,MNP),FSC(NVRT,MNP),FCC(NVRT,MNP)
        REAL*8  :: dr_dxy(2,NVRT,nsa),HTOT,SXX3D0(NVRT,MNP),SYY3D0(NVRT,MNP),SXY3D0(NVRT,MNP), &
     &WILD1(NVRT,MNP),WILD2(NVRT,MNP),WILD3(2,NVRT,nsa),WILD4(3,NVRT,MNP),DSPX,DSPY, &
     &WILD5(10)


        SXX3D(:,:) = 0.
        SYY3D(:,:) = 0.
        SXY3D(:,:) = 0.

        EFTAIL = 1.0 / (PTAIL(1)-1.0)

        ETOT = 0.0
        MDIR = 0.0

        IF (LETOT) THEN
!AR: Estimate zeroth moment m0, mean wave direction, dominant wave number, dominant sigma ...
          DO IP = 1, MNP
            DEPLOC = MAX(DMIN,DEP(IP))
            ACLOC = DBLE(AC2(IP,:,:))
            m0    = 0.
            EWSIG  = 0.
            ETOTS  = 0.
            ETOTC  = 0.
            IF (MSC .GE. 2) THEN
              DO ID = 1, MDC
                m0d = 0.
                DO IS = 2, MSC
                  tmp = 0.5*(SPSIG(IS)*ACLOC(IS,ID)+SPSIG(IS-1)*ACLOC(IS-1,ID))*DS_INCR(IS)*DDIR
                  m0 = m0 + tmp
                  EWSIG  = EWSIG  + SPSIG(IS) * tmp
                  m0d = m0d + tmp
                END DO
                IF (MSC > 3) THEN
                  EHFR = ACLOC(MSC,ID) * SPSIG(MSC)
                  m0 = m0 + DDIR * EHFR * SPSIG(MSC) * EFTAIL
                endif
                ETOTC  = ETOTC + m0d * COS(SPDIR(ID))
                ETOTS  = ETOTS + m0d * SIN(SPDIR(ID))
              END DO
            ELSE
              DS = SGHIG - SGLOW
              DO ID = 1, MDC
                m0d = ACLOC(1,ID) * DS * DDIR
                m0 = m0 + m0d
              END DO
            END IF
            ETOT(IP) = m0
            IF (m0 .GT. small .and. .not. dep(ip) .lt. dmin) then
              EWS(IP) = EWSIG/m0
              WSTMP = EWSIG/m0
              CALL ALL_FROM_TABLE(WSTMP,DEPLOC,WKTMP,WCGTMP,WKDEPTMP,WNTMP,WCTMP)
              EWN(IP) = DBLE(WNTMP)
              EWK(IP) = DBLE(WKTMP)
              MDIR(IP) = DVEC2RAD (ETOTC, ETOTS)
            ELSE
              EWS(IP)  = 0. 
              EWN(IP)  = 0. 
              EWK(IP)  = 10. 
              MDIR(IP) = 0. 
            END IF 
          END DO !IP
        END IF !LETOT

!AR: Here comes the whole story ... 
! Etot = 1/16 * Hs² = 1/8 * Hmono² => Hs² = 2 * Hmono² => Hs = sqrt(2) * Hmono => Hmono = Hs / SQRT(2) 
! Etot = 1/16 * Hs² = 1/16 * (4 * sqrt(m0))² = m0 
! Etot = 1/8 * Hmono² ... so the problem for the analytical solution evolved because we treat the Etot from Hs and Hmono there is a factor of 2 between this!
! Or in other words for the analytical solution we impose a Hs = X[m], we integrate m0 out of it and get Etot, since this Etot is a function of Hs and not Hmono^X^O! it needs the factor of 2 between it! This should make now things clear forever. So the question is not how we calculate the total energy the question is 
! what is defined on the boundary that means we should always recalculate the boundary in terms of Hs =  SQRT(2) * Hmono !!!
! Or saying it again in other words our boundary conditions is wrong if we impose Hmono in wwminput.nml !!!

        SXX3D = 0.
        SXY3D = 0.
        SYY3D = 0.
        WWAVE_FORCE=0.
        IF (RADFLAG .EQ. 'LON') THEN
          RSXX = 0.
          RSXY = 0.
          RSYY = 0.
          DO IP = 1, MNP
            IF (.NOT. LETOT) THEN
              ACLOC = DBLE(AC2(IP,:,:)) 
              DO ID = 1, MDC
                DO IS = 2, MSC
                  ELOC  = 0.5 * (SPSIG(IS)*ACLOC(IS,ID)+SPSIG(IS-1)*ACLOC(IS-1,ID))*DS_INCR(IS)*DDIR
                  COSE2 = COS(SPDIR(ID))**2.
                  SINE2 = SIN(SPDIR(ID))**2.
                  COSI2 = COS(SPDIR(ID)) * SIN(SPDIR(ID))
                  WN    = DBLE(CG(IP,IS) / ( SPSIG(IS)/WK(IP,IS) )) 
                  RSXX(IP) = RSXX(IP) + ( WN * COSE2 + WN - .5) * ELOC   ! Units = [ 1/s + 1/s - 1/s ] * m²s = m²
                  RSXY(IP) = RSXY(IP) + ( WN * COSI2          ) * ELOC
                  RSYY(IP) = RSYY(IP) + ( WN * SINE2 + WN - .5) * ELOC
                ENDDO
              ENDDO
            ELSE IF (LETOT) THEN
              RSXX(IP) =  ETOT(IP) * (EWN(IP)*((EWK(IP)*SIN(MDIR(IP)))**2./EWK(IP)**2.+1.)-0.5)
              RSXY(IP) =  ETOT(IP) *  EWN(IP)* EWK(IP)*SIN(MDIR(IP))*EWK(IP)*COS(MDIR(IP))* 1./EWK(IP)
              RSYY(IP) =  ETOT(IP) * (EWN(IP)*((EWK(IP)*COS(MDIR(IP)))**2./EWK(IP)**2.+1.)-0.5)
            END IF 
          END DO

          DO IP = 1, MNP
            IF (DEP(IP) .GT. DMIN)  THEN
              SXX3D(:,IP) = DBLE(RSXX(IP) / DEP(IP)) * G9
              SXY3D(:,IP) = DBLE(RSXY(IP) / DEP(IP)) * G9
              SYY3D(:,IP) = DBLE(RSYY(IP) / DEP(IP)) * G9
            ELSE
              SXX3D(:,IP) = 0.d0
              SXY3D(:,IP) = 0.d0
              SYY3D(:,IP) = 0.d0
            END IF
          END DO
          !Store as double for force later
          SXX3D0 = SXX3D
          SXY3D0 = SXY3D
          SYY3D0 = SYY3D
        ELSE IF (RADFLAG .EQ. 'XIA') THEN
          IF (LETOT) THEN
            DO IP = 1, MNP
              IF (DEP(IP) .LT. DMIN .OR. IDRY(IP) .EQ. 1 .OR. EWK(IP) .LT. THR) CYCLE
              DO IL = KBP(IP), NVRT
                ZZETA = ZETA(IL,IP)-ZETA(KBP(IP),IP) !from bottom; 'z+D'
                KW           = EWK(IP) * ZZETA !k*(z+D)
                KD           = EWK(IP) * DEP(IP) !k*D
                SINH2KD      = DSINH(MIN(300.d0,2.*KD))
                COSHKD       = DCOSH(MIN(300.d0,KD))
                SINHKW       = DSINH(MIN(300.d0,KW))
                COSH2KW      = DCOSH(MIN(300.d0,2.*KW))
                COSHKW       = DCOSH(MIN(300.d0,KW))
                IF(ABS(SINH2KD) .LT. THR) THEN
                  write(errmsg,*)'R.S.: div by 0 (0);',iplg(IP),EWK(IP),DEP(IP)
                  call parallel_abort(errmsg)
                endif
                tmp         = -ETOT(IP) * EWK(IP) / SINH2KD * (COSH2KW - 1.)                      - &
     &                         ETOT(IP) * (ZETA(IL,IP)-ZETA(NVRT,IP)) / DEP(IP)**2                  + & 
     &                         ETOT(IP) * KW * SINHKW / ( DEP(IP) * COSHKD ) - &
     &                         ETOT(IP) / DEP(IP)  *  (1. - COSHKW / COSHKD)
                SXX3D(IL,IP) = (ETOT(IP) * EWK(IP) / SINH2KD * (COSH2KW + 1.) * COS(MDIR(IP))**2+tmp)* G9
                SYY3D(IL,IP) = (ETOT(IP) * EWK(IP) / SINH2KD * (COSH2KW + 1.) * SIN(MDIR(IP))**2+tmp)* G9
                SXY3D(IL,IP) = ETOT(IP) * EWK(IP) / SINH2KD * (COSH2KW + 1.) * COS(MDIR(IP))*SIN(MDIR(IP))* G9
              END DO !IL
            END DO !IP
          ELSE !not mono ... random waves version ... treat eveery wave packet like a single monochr. wave with the height of Hs 
            DO IP = 1, MNP
              IF (DEP(IP) .LT. DMIN .OR. IDRY(IP) .EQ. 1) CYCLE
              ACLOC = AC2(IP,:,:)
              DO IL = KBP(IP), NVRT
                ZZETA = ZETA(IL,IP)-ZETA(KBP(IP),IP) !from bottom
                DO IS = 1, MSC !freq
                  KW           = WK(IP,IS) * ZZETA !k*(z+D)
                  KD           = WK(IP,IS) * DEP(IP) !k*D
                  SINH2KD      = DSINH(MIN(300.d0,2.*KD))
                  COSHKD       = DCOSH(MIN(300.d0,KD))
                  SINHKW       = DSINH(MIN(300.d0,KW))
                  COSH2KW      = DCOSH(MIN(300.d0,2.*KW))
                  COSHKW       = DCOSH(MIN(300.d0,KW))
                  IF(ABS(SINH2KD) .LT. THR) call parallel_abort('R.S.: div 0 (0)')
                  DO ID = 1, MDC !direction
                    !Dimension of ELOC = m^2
                    ELOC = ACLOC(IS,ID) * SPSIG(IS)**2. * DDIR * FRINTF
                    IF (ELOC .LT. SMALL) CYCLE
                      tmp          =-ELOC * WK(IP,IS) / SINH2KD * (COSH2KW - 1.)            - &
     &                               ELOC * (ZETA(IL,IP)-ZETA(NVRT, IP)) / DEP(IP)**2.        + &
     &                               ELOC * KW * SINHKW / ( DEP(IP) * COSHKD ) - &
     &                               ELOC / DEP(IP) *  (1. - COSHKW / COSHKD)
                      tmp=tmp*G9

                      SXX3D(IL,IP) = SXX3D(IL,IP) + &
     &                               G9*ELOC * WK(IP,IS) / SINH2KD * (COSH2KW + 1.) * COS(SPDIR(ID))**2.+tmp
                      SYY3D(IL,IP) = SYY3D(IL,IP) + &
     &                               G9*ELOC * WK(IP,IS) / SINH2KD * (COSH2KW + 1.) * SIN(SPDIR(ID))**2.+tmp
                      SXY3D(IL,IP) = SXY3D(IL,IP)+G9*ELOC * WK(IP,IS) / SINH2KD * (COSH2KW + 1.) * SIN(SPDIR(ID))*COS(SPDIR(ID))
                  END DO
                END DO
              END DO !IL
            END DO !IP
          END IF !mono

          !Store as double for force later
          SXX3D0=SXX3D
          SXY3D0=SXY3D
          SYY3D0=SYY3D

!ZYL: Mellor 2003; random waves only; only good for traditional sigma coord.
!     did not include the last term in Warner et al. (2008) as it is not in Mellor 2005
        ELSE IF (RADFLAG .EQ. 'MEL') THEN
          IF(LETOT) call parallel_abort('R.S.: no mono for Mellor 88')
          IF(KZ/=1.or.THETA_F>1.e-4) call parallel_abort('R.S.: MEL must use sigma')
          !WILD1: \sum_{dir}{E} at nodes; WILD2: K*D at nodes
          WILD1=0; WILD2=0; WILD4=0 !for exceptions
          FSS=0; FCS=0; FSC=0; FCC=0
          DO IS = 1, MSC !freq
            DO IP = 1, MNP
              IF (DEP(IP) .LT. DMIN .OR. IDRY(IP) .EQ. 1) CYCLE
              KD           = WK(IP,IS) * DEP(IP) !k*D
              WILD2(:,IP)=KD
              SINHKD       = DSINH(MIN(300.d0,KD))
              COSHKD       = DCOSH(MIN(300.d0,KD))
              IF(ABS(SINHKD) .LT. THR) call parallel_abort('R.S.: div by 0 (1)')
              DO IL = KBP(IP), NVRT
                ZZETA = ZETA(IL,IP)-ZETA(KBP(IP),IP) !from bottom
                KW           = WK(IP,IS) * ZZETA !k*(z+D)
                SINHKW       = DSINH(MIN(300.d0,KW))
                COSHKW       = DCOSH(MIN(300.d0,KW))
                FSS(IL,IP)   = SINHKW/SINHKD
                FCS(IL,IP)   = COSHKW/SINHKD
                FSC(IL,IP)   = SINHKW/COSHKD
                FCC(IL,IP)   = COSHKW/COSHKD

                DO ID = 1, MDC !direction
                  !Dimension of ELOC = m^2
                  ELOC = AC2(IP,IS,ID) * SPSIG(IS)**2. * DDIR * FRINTF * 1.
                  IF (ELOC .LT. SMALL) CYCLE
                  WILD1(IL,IP)=WILD1(IL,IP)+ELOC
                  SXX3D(IL,IP) = SXX3D(IL,IP)+G9*ELOC*WK(IP,IS)*(FCS(IL,IP)*FCC(IL,IP)*(COS(SPDIR(ID))**2+1.)-FSS(IL,IP)*FCS(IL,IP)) 
                  SYY3D(IL,IP) = SYY3D(IL,IP)+G9*ELOC*WK(IP,IS)*(FCS(IL,IP)*FCC(IL,IP)*(SIN(SPDIR(ID))**2+1.)-FSS(IL,IP)*FCS(IL,IP))
                  SXY3D(IL,IP) = SXY3D(IL,IP)+G9*ELOC*WK(IP,IS)*FCS(IL,IP)*FCC(IL,IP)*SIN(SPDIR(ID))*COS(SPDIR(ID))
                END DO !ID
              END DO !IL
            END DO !IP

            !Compute horizontal derivatives of WILD[1,2], and temporarily store them in dr_dxy and WILD3
            !dr_dxy: d{sum{E}}/d{x,y} - not a function of sigma
            !WILD3: d{K*D}/d{x,y} - not a function of sigma
            !Can we sum up E or KD over freq. before differentiation to save exchange time?
            call hgrad_nodes(0,NVRT,MNP,nsa,WILD1,dr_dxy)
            call hgrad_nodes(0,NVRT,MNP,nsa,WILD2,WILD3)
            call exchange_s3d_2(dr_dxy)
            call exchange_s3d_2(WILD3)

            !Compute vertical derivatives of FCC etc. and store them in WILD4
            DO IP = 1, MNP
              IF (DEP(IP) .LT. DMIN .OR. IDRY(IP) .EQ. 1) CYCLE
              KD= WK(IP,IS) * DEP(IP) !k*D
              DO IL = KBP(IP), NVRT
                WILD4(1,IL,IP)=FSC(IL,IP)*KD !d{FCC}/d{sigma}
                WILD4(2,IL,IP)=FCS(IL,IP)*KD !d{FSS}/d{sigma}
                WILD4(3,IL,IP)=FSS(IL,IP)*KD !d{FCS}/d{sigma}
              END DO !IL
            ENDDO !IP

            !Vertical derivative parts of wave forces (more later for horizontal parts): D{SPX}/D{sigma}
            do IP=1,nsa !sides
              if(idry_s(IP)==1) CYCLE
              ND1=isidenode(IP,1) !1st node of the side
              ND2=isidenode(IP,2)
              HTOT=(eta2(ND1)+eta2(ND2)+DEP8(ND1)+DEP8(ND2))/2
              if(HTOT<=0) call parallel_abort('RADIATION_STRESS: (996)')
              DO IL=1,NVRT
                WILD5(1)=(WILD4(1,IL,ND1)+WILD4(1,IL,ND2))/2 !d{FCC}/d{sigma} at side
                WILD5(2)=(WILD4(2,IL,ND1)+WILD4(2,IL,ND2))/2 !d{FSS}/d{sigma} at side
                WILD5(3)=(WILD4(3,IL,ND1)+WILD4(3,IL,ND2))/2 !d{FCS}/d{sigma} at side
                WILD5(4)=(FCC(IL,ND1)+FCC(IL,ND2))/2 !FCC at side
                WILD5(5)=(FSS(IL,ND1)+FSS(IL,ND2))/2 !FSS at side
                WILD5(6)=(FCS(IL,ND1)+FCS(IL,ND2))/2 !FCS at side
                WILD5(7)=(WILD1(IL,ND1)+WILD1(IL,ND2))/2 !sum{E} at side
                WILD5(8)=(WILD2(IL,ND1)+WILD2(IL,ND2))/2 !K*D at side
                DSPX=dr_dxy(1,IL,IP)/2*((WILD5(1)-WILD5(2))*WILD5(5)+(WILD5(4)-WILD5(5))*WILD5(2))+&
     &WILD5(7)*WILD3(1,IL,IP)*((WILD5(1)-WILD5(2))*WILD5(6)*(1+SIGMACOR(IL))+ &
     &(WILD5(4)-WILD5(5))*WILD5(3)*(1+SIGMACOR(IL))+(WILD5(4)-WILD5(5))*WILD5(6))
                DSPY=dr_dxy(2,IL,IP)/2*((WILD5(1)-WILD5(2))*WILD5(5)+(WILD5(4)-WILD5(5))*WILD5(2))+&
     &WILD5(7)*WILD3(2,IL,IP)*((WILD5(1)-WILD5(2))*WILD5(6)*(1+SIGMACOR(IL))+ &
     &(WILD5(4)-WILD5(5))*WILD5(3)*(1+SIGMACOR(IL))+(WILD5(4)-WILD5(5))*WILD5(6))

                WWAVE_FORCE(IL,IP,1)=WWAVE_FORCE(IL,IP,1)+G9*DSPX/HTOT !m/s/s
                WWAVE_FORCE(IL,IP,2)=WWAVE_FORCE(IL,IP,2)+G9*DSPY/HTOT
              END DO !IL
            enddo !IP; sides
          END DO !IS

          !Store as double for force later
          SXX3D0=SXX3D
          SXY3D0=SXY3D
          SYY3D0=SYY3D

        ELSE
            call parallel_abort('R.S.: unknown R.S. model') 
        END IF !RADFLAG 

!       Integrate over depth for checking
        RSXX = 0.
        DO IP = 1, MNP
          IF (DEP(IP) .LE. DMIN .OR. IDRY(IP) .EQ. 1) CYCLE
          DO IL = KBP(IP)+1, NVRT 
            RSXX(IP) = RSXX(IP) + 0.5*( SXX3D(IL,IP)+SXX3D(IL-1,IP) ) * ABS((ZETA(IL,IP) - ZETA(IL-1,IP)))/G9
          END DO !IL
        END DO !IP

!       Computation in double precision here
!        IF (RADFLAG.EQ.'LON'.OR.RADFLAG.EQ.'XIA') THEN
!         SXX3D0() etc. should have dimension of m^2/s/s, defined at nodes and whole levels.
!         Use same arrays to temporarily store properly scaled Sxx etc
!         write(12,*)'Checking Sxx,Sxy,Syy:'
          do IP=1,MNP
            if(IDRY(IP)==1) then
              SXX3D0(:,IP)=0
              SYY3D0(:,IP)=0
              SXY3D0(:,IP)=0
              cycle
            endif

            do IL=KBP(IP),NVRT
!             D*(Sxx, Sxy, Syy)/rho in Xia et al. (2004) 
!             After this the dimension of sdbt should be m^3/s/s
              SXX3D0(IL,IP)=(DEP8(IP)+ETA2(IP))*SXX3D0(IL,IP) !D*Sxx/rho
              SXY3D0(IL,IP)=(DEP8(IP)+ETA2(IP))*SXY3D0(IL,IP) !D*Sxy/rho
              SYY3D0(IL,IP)=(DEP8(IP)+ETA2(IP))*SYY3D0(IL,IP) !D*Syy/rho
            enddo !k
          enddo !IP

!         Compute radiation stress force 
!         wwave_force(:,1:nsa,1:2) = Rsx, Rsy in my notes (the terms in momen. eq.)
!         and has a dimension of m/s/s
          call hgrad_nodes(0,NVRT,MNP,nsa,SXX3D0,dr_dxy)
          call exchange_s3d_2(dr_dxy)
          do IS=1,nsa
            if(idry_s(IS)==0) then
              HTOT=(eta2(isidenode(IS,1))+eta2(isidenode(IS,2))+DEP8(isidenode(IS,1))+DEP8(isidenode(IS,2)))/2
              if(HTOT<=0) call parallel_abort('RADIATION_STRESS: (999)')
              do il = 1, nvrt
                WWAVE_FORCE(il,IS,1)=WWAVE_FORCE(il,IS,1)-dr_dxy(1,il,IS)/HTOT
              end do
            endif
          enddo !IS

          call hgrad_nodes(0,NVRT,MNP,nsa,SYY3D0,dr_dxy)
          call exchange_s3d_2(dr_dxy)
          do IS=1,nsa
            if(idry_s(IS)==0) then
              HTOT=(eta2(isidenode(IS,1))+eta2(isidenode(IS,2))+DEP8(isidenode(IS,1))+DEP8(isidenode(IS,2)))/2
              if(HTOT<=0) call parallel_abort('RADIATION_STRESS: (998)')
              do il = 1, nvrt 
                WWAVE_FORCE(il,IS,2)=WWAVE_FORCE(il,IS,2)-dr_dxy(2,il,IS)/HTOT
              end do
            endif
          enddo !IS

          call hgrad_nodes(0,NVRT,MNP,nsa,SXY3D0,dr_dxy)
          call exchange_s3d_2(dr_dxy)

!Debug
!          write(12,*)'Checking R.S.'
          do IS=1,nsa
            if(idry_s(IS)==0) then
              HTOT=(eta2(isidenode(IS,1))+eta2(isidenode(IS,2))+DEP8(isidenode(IS,1))+DEP8(isidenode(IS,2)))/2
              if(HTOT<=0) call parallel_abort('RADIATION_STRESS: (997)')
              WWAVE_FORCE(:,IS,1)=WWAVE_FORCE(:,IS,1)-dr_dxy(2,:,IS)/HTOT
              WWAVE_FORCE(:,IS,2)=WWAVE_FORCE(:,IS,2)-dr_dxy(1,:,IS)/HTOT
            endif
          enddo !IS

        END SUBROUTINE
#endif 
!**********************************************************************
!*                                                                    *
!**********************************************************************
