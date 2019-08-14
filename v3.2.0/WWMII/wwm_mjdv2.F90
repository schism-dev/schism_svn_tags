!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CT2MJD(STIME,XMJD)
         IMPLICIT NONE
         CHARACTER(LEN=15), INTENT(IN) :: STIME
         DOUBLE PRECISION, INTENT(INOUT) :: XMJD
!
! ... FORMAT IS YYYYMMDD.HHMMSS , LENGTH IS 15
!
         INTEGER :: IY = 0
         INTEGER :: IM = 0
         INTEGER :: ID = 0
         INTEGER :: IH = 0
         INTEGER :: IMIN  = 0
         INTEGER :: ISEC  = 0
         INTEGER :: IFLAG = 1

         READ(STIME(1:4),  *,END=100,ERR=100) IY
         READ(STIME(5:6),  *,END=100,ERR=100) IM
         READ(STIME(7:8),  *,END=100,ERR=100) ID
         READ(STIME(10:11),*,END=100,ERR=100) IH
         READ(STIME(12:13),*,END=100,ERR=100) IMIN
         READ(STIME(14:15),*,END=100,ERR=100) ISEC

         CALL MJDYMD(XMJD  , IY    , IM    , ID    , IH    ,   &
     &               IMIN  , ISEC  , IFLAG                    )
         RETURN
         
100      XMJD = 0.0     
         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MJD2CT(XMJD,STIME)
         IMPLICIT NONE
         CHARACTER(LEN=15), INTENT(OUT) :: STIME
         DOUBLE PRECISION, INTENT(IN) :: XMJD
         DOUBLE PRECISION :: TMJD
         INTEGER :: IY, IM, ID, IH, IMIN, ISEC

         INTEGER :: IFLAG = 2

         TMJD = XMJD

         CALL MJDYMD( TMJD, IY, IM, ID, IH, IMIN, ISEC, IFLAG )

         WRITE(STIME(1:4),'(I4.4)') IY
         WRITE(STIME(5:6),'(I2.2)') IM
         WRITE(STIME(7:8),'(I2.2)') ID
         WRITE(STIME(9:9),'(A)') '.'
         WRITE(STIME(10:11),'(I2.2)') IH
         WRITE(STIME(12:13),'(I2.2)') IMIN
         WRITE(STIME(14:15),'(I2.2)') ISEC

         RETURN
      END SUBROUTINE

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE CU2SEC(UNITT, DT)
         USE DATAPOOL, ONLY : DBG
         IMPLICIT NONE
         CHARACTER(LEN=*), INTENT(IN) :: UNITT
         REAL*8, INTENT(INOUT) :: DT

         SELECT CASE (UNITT)
            CASE ('H', 'h', 'HR', 'hr')
               DT = DT * 3600.0
            CASE ('M', 'm', 'MIN', 'min')
               DT = DT * 60.0
            CASE ('S', 's', 'SEC', 'sec')
               DT = DT
            CASE DEFAULT
               WRITE(DBG%FHNDL,*) 'ERROR WRONG UNIT, UNIT = ', UNITT
               DT = 0.0
         END SELECT

         RETURN
      END SUBROUTINE

!**********************************************************************
!*                                                                    *
!**********************************************************************
      SUBROUTINE MJDYMD( XMJD, IY, IM, ID, IH, IMIN, ISEC, IFLAG )
        USE DATAPOOL, ONLY : SEC2DAY, DAY2SEC
!
! ... XMJD  : MODIFIED JULIAN DATE
! ... IY    : YEAR
! ... IM    : MONTH
! ... ID    : DAY
! ... IH    : HOUR
! ... IMIN  : MINUTE
! ... ISEC  : SECOND
! ... IFLAG : 1 -> YMDHMS TO MJD
!             2 -> MJD TO YMDHMS
! ... DATE MUST BE WITHIN THE YEARS MAR. 1, 1900 TO FEB. 28, 2100

#ifdef SELFE
         use elfe_msgp
#endif
         IMPLICIT DOUBLE PRECISION (A-H,O-Z)
!
         PARAMETER ( XJD0 = 2400000.5D0 )
         PARAMETER ( HALF =       0.5D0 )
         INTEGER IMONTH(12)
         DATA IMONTH /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/
!
! ... -----< YMDHMS TO MJD >-----
!
         IF (IFLAG.EQ.1) THEN

            Y = DFLOAT(IY - 1)

            IF (IM.GT.2) THEN
               M = IM
               Y = Y + 1
            ELSE
               M = IM + 12
            ENDIF

            XJD  = INT(365.25D0*Y) + INT(30.6001D0*(M+1)) - 15  &
     &           + 1720996.5D0     + ID
            XMJD = XJD - XJD0

            FSEC = DFLOAT(IH)*3600.D0 + DFLOAT(IMIN)*60.D0 + DFLOAT(ISEC)

            XMJD = XMJD + FSEC * SEC2DAY
!
! ... -----< MJD TO YMDHMS >-----
!
         ELSE IF (IFLAG.EQ.2) THEN

            MJD  = XMJD
            XJD  = DFLOAT(MJD) + XJD0
            C    = INT(XJD + HALF) + 1537
            ND   = INT((C - 122.1D0)/365.25D0 )
            E    = INT(365.25D0*ND)
            NF   = INT((C - E)/30.6001D0)

            IFR  = INT(XJD + HALF)
            FRC  = XJD + HALF - DFLOAT(IFR)
            ID   = C - E - INT(30.6001D0*NF) + FRC
            IM   = NF - 1 - 12*INT(NF/14)
            IY   = ND - 4715 - INT((7+IM)/10)

            SEC  = (XMJD-DFLOAT(MJD))*DAY2SEC
            ISEC = SEC
            IF ((SEC-ISEC).GT.0.5D0) ISEC = ISEC + 1
            IH   = ISEC/3600
            IMIN = (ISEC - IH*3600)/60
            ISEC = ISEC - IH*3600 - IMIN*60
!
! ... set 24:00 to 00:00
!
            IF (MOD(IY,4) == 0) IMONTH(2) = 29

            IF (IH == 24) THEN
               IH = 0
               ID = ID + 1
               IF (ID > IMONTH(IM)) THEN
                  ID = ID - IMONTH(IM)
                  IM = IM + 1
                  IF (IM > 12) THEN
                     IM = IM - 12
                     IY = IY + 1
                  END IF
               END IF
            END IF
         ELSE

#ifdef SELFE
            call parallel_abort('!!! ERROR IN <MJDYMD>. IFLAG SHOULD BE 1 OR 2.')
#else
            PRINT*,'!!! ERROR IN <MJDYMD>. IFLAG SHOULD BE 1 OR 2.'
            STOP 'MJDYMD'
#endif

         ENDIF

         RETURN
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
