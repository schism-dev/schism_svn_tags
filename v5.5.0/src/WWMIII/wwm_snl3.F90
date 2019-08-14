#include "wwm_functions.h"
!**********************************************************************
!*  Yasser Eldeberky, Nonlinear Transformation of Wave Spectra in     *
!*        the Nearshore Zone, PhD thesis, TU Delft                    *
!**********************************************************************
      SUBROUTINE TRIAD_ELDEBERKY(ip, hs, smespc, WALOC, ssnl3, dssnl3)
      use datapool
      implicit none
      integer, intent(in)        :: ip
      real(rkind), intent(in)    :: hs, smespc
      real(rkind), intent(inout) :: ssnl3(NUMSIG,NUMDIR), dssnl3(NUMSIG,NUMDIR)
      real(rkind), intent(in)    :: WALOC(NUMSIG,NUMDIR)
      integer id, is, ismax
      real(rkind)    aux1, aux2, biph, c0, cm, e0
      real(rkind)    em,ft, rint, sigpi, sinbph, stri
      real(rkind)    w0, wm, wn0, wnm, ursell, c1, c2, c3
      real(rkind)    eCont
      real(rkind) :: E(NUMSIG)
      real(rkind) :: SA(1:NUMSIG+TRI_ISP1,1:NUMDIR)

      IF (HS .LT. SMALL) RETURN
      CALL URSELL_NUMBER(HS,SMESPC,DEP(IP),URSELL) 
      IF ( URSELL .le. TRI_ARR(5) ) RETURN

      E  = ZERO
      SA = ZERO
      ISMAX = TRI_ISP1

      DO IS = TRI_ISP1, NUMSIG
       IF ( SPSIG(IS) .LT. ( TRI_ARR(2) * SMESPC) ) ISMAX = IS
      ENDDO

      c1     = G9*(DEP(IP)**2)
      c2     = (TWO/15._rkind)*G9*(DEP(IP)**4)
      c3     = (TWO/5._rkind)*(DEP(IP)**3)
      BIPH   = PIHALF*(MyTANH(TRI_ARR(4)/URSELL)-1.)
      SINBPH = ABS( SIN(BIPH) )

      DO ID = 1, NUMDIR
        E = WALOC(:,ID) * PI2 * SPSIG
        DO IS=TRI_ISBEGIN, ISMAX 
          W0  = SPSIG(IS)
          WN0 = WK(IS,IP)
          C0  = W0 / WN0
          EM  = TRI_WISM * E(IS+TRI_ISM1)      + TRI_WISM1 * E(IS+TRI_ISM)
          WM  = TRI_WISM * SPSIG(IS+TRI_ISM1)  + TRI_WISM1 * SPSIG(IS+TRI_ISM)
          WNM = TRI_WISM * WK(IS+TRI_ISM1,IP)  + TRI_WISM1 * WK(IS+TRI_ISM,IP)
          CM  = WM / WNM
          AUX1 = WNM**2 * ( G9 * DEP(IP) + TWO*CM**2 )
          AUX2 = WN0 * ( c1 + c2 * WN0**2 - c3 * W0**2) ! (m/s² * m + m/s² * m³*1/m² - 1/s² * m²)
          RINT = AUX1 / AUX2
          FT = TRI_ARR(1) * C0 * CG(IS,IP) * RINT**2 * SINBPH
          SA(IS,ID) = MAX(ZERO, FT*(EM*(EM - 2*E(IS)))) ! Here the max opeartor allows only shift to higher harmonics so it is non-conservative ?!?! Total crap ... must be deleted ...
        END DO
      END DO

      DO ID = 1, NUMDIR
        DO IS = 1, NUMSIG
          SIGPI = SPSIG(IS) * PI2
          STRI  = SA(IS,ID) - TWO * (TRI_WISP*SA(IS+TRI_ISP1,ID) + TRI_WISP1*SA(IS+TRI_ISP,ID))
          eCont = STRI / SIGPI
          IF (ABS(WALOC(IS,ID)) .GT. SMALL) THEN
            DSSNL3(IS,ID) = eCont
             SSNL3(IS,ID) = eCont / WALOC(IS,ID)
          ELSE
            DSSNL3(IS,ID) = ZERO
             SSNL3(IS,ID) = ZERO
          ENDIF
        END DO
      END DO

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
