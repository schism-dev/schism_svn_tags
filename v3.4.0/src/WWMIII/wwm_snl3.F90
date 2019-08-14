#include "wwm_functions.h"
!**********************************************************************
!*                                                                    *
!**********************************************************************
!2do add init function
      subroutine triadswan_new (ip, hs, smespc, acloc, imatra, imatda, ssnl3, dssnl3)
!
      use datapool
      implicit none

      integer, intent(in)        :: ip
      real(rkind), intent(in)    :: hs, smespc
      real(rkind), intent(out)   :: ssnl3(msc,mdc), dssnl3(msc,mdc)
      real(rkind), intent(in)    :: acloc(msc,mdc)
      real(rkind), intent(inout) :: imatra(msc,mdc), imatda(msc,mdc)
!
      integer i1, i2, id, is, ism, ism1, ismax, isp, isp1,  ij1, ij2, ires

      real(rkind)    aux1, aux2, biph, c0, cm, dep_2, dep_3, e0, bb
      real(rkind)    em,ft, rint, sigpi, sinbph, stri, wism, wism1 , fac1
      real(rkind)    wisp, wisp1,w0, wm, wn0, wnm,  xisln, ursell,facres,facscl,siglow
      
      real(rkind) :: E(MSC)
      real(rkind), allocatable :: sa(:,:)

      PTRIAD(1)  = 0.1
      PTRIAD(2)  = 2.2
      PTRIAD(3)  = 10.
      PTRIAD(4)  = 0.2
      PTRIAD(5)  = 0.01

      ssnl3 = zero 
      dssnl3 = zero

      IF (TRICO .GT. 0.)  PTRIAD(1) = TRICO
      IF (TRIRA .GT. 0.)  PTRIAD(2) = TRIRA
      IF (TRIURS .GT. 0.) PTRIAD(5) = TRIURS

      BB = 1./15.

      IF (HS .LT. SMALL) RETURN

      CALL URSELL_NUMBER(HS,SMESPC,DEP(IP),URSELL) 

      !if (ip == 1786) write(stat%fhndl,'(A20,I10,8F15.10)') 'URSELL',IP,DEP(IP),HS,SMESPC,URSELL,(G9 * HS),(TWO*SQRT(TWO)*SMESPC**2*DEP(IP)**2)

!      write(*,*) '---- calling snl3 -----', ip, iobp(ip)

      IJ2    = INT (FLOAT(MSC) / 2.)
      IJ1    = IJ2 - 1
      FAC1   = SPSIG(IJ2) / SPSIG(IJ1)
      IRES   = NINT ( LOG10( 2.) / LOG10( FAC1 ) )
      FACSCL = SPSIG(MSC)/SPSIG(MSC-IRES)

      IF (ABS(FACSCL-2.).GT.0.05) THEN
         FACRES = 10.**( LOG10(2.) / FLOAT(IRES) )
         SIGLOW   = SPSIG(MSC) / ( FACRES**(FLOAT(MSC-1) ) )
         !WRITE(DBG%FHNDL,*) 'CHECK RESOLUTION', IRES, FACSCL, FACRES, SIGLOW
      END IF

      DEP_2 = DEP(IP)**2
      DEP_3 = DEP(IP)**3
      I2     = INT (FLOAT(MSC) / 2.)
      I1     = I2 - 1
      XIS    = SPSIG(I2) / SPSIG(I1)
      XISLN  = LOG( XIS )
      ISP    = INT( LOG(2.) / XISLN )
      ISP1   = ISP + 1
      WISP   = (2. - XIS**ISP) / (XIS**ISP1 - XIS**ISP)
      WISP1  = 1. - WISP
      ISM    = INT( LOG(0.5) / XISLN )
      ISM1   = ISM - 1
      WISM   = (XIS**ISM -0.5) / (XIS**ISM - XIS**ISM1)
      WISM1  = 1. - WISM

      ALLOCATE (SA(1:MSC+ISP1,1:MDC))

      E  = 0.
      SA = 0.

      ISMAX = 1
      DO IS = 1, MSC
       IF ( SPSIG(IS) .LT. ( PTRIAD(2) * SMESPC) ) THEN
          ISMAX = IS
        ENDIF
      ENDDO
!      ISMAX = MIN( MSC, MAX ( ISMAX , ISP1 ) ) ! added fix the bug described below ...
      ISMAX = MAX ( ISMAX , ISP1 )
       
      IF ( URSELL .GT. PTRIAD(5) ) THEN

        BIPH   = (0.5*PI)*(MyTANH(PTRIAD(4)/URSELL)-1.)
        SINBPH = ABS( SIN(BIPH) )

        DO ID = 1, MDC
          E = ACLOC(:,ID) * PI2 * SPSIG
          DO IS = 1, ISMAX 
            E0  = E(IS)
            W0  = SPSIG(IS)
            WN0 = WK(IS,IP)
            C0  = W0 / WN0
            IF ( IS.GT.-ISM1 ) THEN
              EM  = WISM * E(IS+ISM1)      + WISM1 * E(IS+ISM)
              WM  = WISM * SPSIG(IS+ISM1)  + WISM1 * SPSIG(IS+ISM)
              WNM = WISM * WK(IS+ISM1,IP)  + WISM1 * WK(IS+ISM,IP)
              CM  = WM / WNM
            ELSE
               EM  = ZERO
               WM  = ZERO 
               WNM = ZERO 
               CM  = ZERO 
            END IF
            AUX1 = WNM**2 * ( G9 * DEP(IP) + 2.*CM**2 )
            AUX2 = WN0 * DEP(IP) * ( G9 * DEP(IP) + (2./15.) * G9 * DEP_3 * WN0**2 - (2./5.) * W0**2 * DEP_2 ) ! (m/s² * m + m/s² * m³*1/m² - 1/s² * m²)
            RINT = AUX1 / AUX2
            FT = PTRIAD(1) * C0 * CG(IS,IP) * RINT**2 * SINBPH
            SA(IS,ID) = MAX(ZERO, FT * ( EM * (EM - 2*E0 )))
            !IF (IP == 1786 .AND. SA(IS,ID) .GT. THR) WRITE(*,'(2I10,10F25.10)') IS, ID, DEP(IP), URSELL, PTRIAD(5), RINT**2, SINBPH, SA(IS,ID), FT, ( EM * (EM - 2*E0 ))
          END DO
        END DO

        DO IS = 1, MSC
          SIGPI = SPSIG(IS) * PI2
          DO ID = 1, MDC
            IF (ACLOC(IS,ID) .LT. THR) CYCLE
            STRI = SA(IS,ID) - 2.*(WISP  * SA(IS+ISP1,ID) + WISP1 * SA(IS+ISP,ID))
            IF (ABS(STRI) .LT. THR) CYCLE
            !IF (IP == 1786)  WRITE(*,'(2I10,4F15.10,I10)') IS, ID, STRI, SA(IS,ID), SA(IS+ISP1,ID) , SA(IS+ISP,ID), ISP+IS
            IF (ICOMP .GE. 2) THEN
              !SSNL3(IS,ID)  = STRI / SIGPI 
              IF (STRI .GT. 0.) THEN
                IMATRA(IS,ID) = IMATRA(IS,ID) + STRI / SIGPI
                SSNL3(IS,ID)  =  STRI / SIGPI 
              ELSE
                IMATDA(IS,ID) = IMATDA(IS,ID) - STRI / (ACLOC(IS,ID)*SIGPI)
                DSSNL3(IS,ID) =  -STRI/(ACLOC(IS,ID)*SIGPI)
              END IF
              !IF (IP == TESTNODE) write(*,'(3I10,2F20.10)') IP, IS, ID, SSNL3(IS,ID), DSSNL3(IS,ID)
            ELSE
              IMATRA(IS,ID) = IMATRA(IS,ID) + STRI / SIGPI
              IMATDA(IS,ID) = ZERO!IMATDA(IS,ID) + STRI / (ACLOC(IS,ID)*SIGPI)
              SSNL3(IS,ID)  = STRI / SIGPI
              DSSNL3(IS,ID) = STRI / (ACLOC(IS,ID)*SIGPI) 
            END IF
          END DO
        END DO
      END IF
 
      deallocate(sa)

      end subroutine 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      subroutine triad_dingemans (ip, acloc, imatra, imatda, ssnl3)
!
        use datapool
        implicit none

        integer, intent(in) :: ip
        real(rkind), intent(in)    :: acloc(msc,mdc)
        real(rkind), intent(inout) :: imatra(msc,mdc), imatda(msc,mdc)
        real(rkind), intent(out)   :: ssnl3(msc,mdc)
        integer             :: is, is2, id
        real(rkind)         :: ecloc(msc,mdc), e2(msc,mdc), d20
        real(rkind)         :: df, domega, omega, omega1, fac, z1a, z1b
        integer             :: j1, j2, jmin, j2abs

        do is = 1, msc
          do id = 1, mdc 
            ecloc(is,id) = acloc(is,id) * spsig(is) * ddir
          end do 
        end do

        ssnl3 = 0.

        do id = 1, mdc
          do is = 1,msc-1
            df = ((spsig(is+1) - spsig(is)))/pi2
            domega = pi2 * df
            omega = spsig(is) * pi2 
            jmin = nint(0.5*MyREAL(is))
            if (2*jmin .eq. is) then
               fac = 0.5 
            else
               fac = 1. 
            endif
            z1a = dep(ip) * (jmin*domega)**2 / g9
            z1b = dep(ip) * ((jmin-is-1)*domega)**2 / g9
            e2(is,id) = 0. 
            do is2 = jmin, msc
               j1    = is2 
               j2    = is2 - is 
               j2abs = iabs (j2)
!            zero frequencies are skipped
               if (j2 .eq. 0)  cycle
               omega1 = is2 * domega
               e2(is,id) = e2(is,id) + fac * d20(omega,omega1,dep(ip),z1a,z1b)* ecloc(j1,id)*ecloc(j2abs,id)
               fac    = 1. 
            end do
            e2(is,id) = e2(is,id) * df
          end do
        end do
      end subroutine
!**********************************************************************
!*                                                                    *
!**********************************************************************
      function d20(omega, omega1, h, z1a, z1b)
      use datapool, only : g9, rkind
      implicit none
!
!----------------------------------------------------------------------*
!     compute the second order spectral response function.             *
!--------------------------------- delft hydraulics -- gkl -- 880419 --*
!
      real(rkind), intent(in) :: omega, omega1, h
      real(rkind), intent(inout) :: z1a, z1b
      real(rkind) :: k, k1, k2
      real(rkind) :: omega2
      real(rkind) :: aome2
      real(rkind) :: cthkh, cthk1h, cthk2h 
      real(rkind) :: d2, d20
!
      omega2 = omega - omega1
      aome2  = abs (omega2)
      z1a = abs(z1a)
      z1b = abs(z1b)
      call dispu2 (omega1, h, z1a, k1)
      call dispu2 (aome2,  h, z1b, k2)
      if (omega2 .lt. 0.) k2 = - k2
      k = k1 + k2
      cthk1h = 1. / MyTANH (k1*h)
      cthk2h = 1. / MyTANH (k2*h)
      cthkh  = 1. / MyTANH (k*h)
      d2 = 0.5 *  (omega1**2 + omega2**2 + omega1*omega2 -              &
     &             omega1 * omega2 * cthk1h * cthk2h - omega *          & 
     &             omega1 * cthkh  * cthk1h - omega  * omega2 *         &
     &             cthkh  * cthk2h)
      d20 = d2 / (g9 * (1. - omega**2 * cthkh / (g9*k)))
      d20 = 2. * d20**2
!
      return
      end function 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      subroutine dispu2(omega,d,z1,k)
      use datapool, only : rkind
      implicit none

      real(rkind), intent(in)    :: omega, d
      real(rkind), intent(inout) :: z1
      real(rkind), intent(out)   :: k
      real(rkind), parameter     :: eps = 0.0001
      real(rkind), parameter     :: g = 9.81

      real(rkind) :: z0, z2, fak1, fak2, sig 

      z0 = d*omega*omega/g
   10    sig = MyTANH(z1)
         fak1 = z1*sig
         fak2 = z1 + sig*(1.-fak1)
         z2 = z1 + (z0-fak1)/fak2
         if (abs((z2-z1)/z2).gt.eps) goto 40
         goto 60
   40    z1 = z2
         goto 10
   60 k = z2/d
      z1 = z2
      return
      end subroutine
!**********************************************************************
!*                                                                    *
!**********************************************************************
      function delta(ip, is, is1, is2) result(res)
      use datapool, only : rkind, wk
      implicit none

      integer, intent(in)        :: ip, is, is1, is2
      real(rkind)                :: res !return

        res = wk(ip,is) - wk(ip,is1) - wk(ip,is2)

!        write(*,*) 'delta', res 
        
      end function delta
!**********************************************************************
!*                                                                    *
!**********************************************************************
      function ddelta_dx(ip, is, is1, is2, id) result(res)
      use datapool, only : rkind, sinth, costh, dwkdx, dwkdy
      implicit none

      integer, intent(in)        :: ip, is, is1, is2, id
      real(rkind)                :: res !return

      real(rkind)                :: ka_1abl, kb_1abl, kc_1abl

        ka_1abl = costh(id) * dwkdx(ip,is) + sinth(id) * dwkdy(ip,is)
        kb_1abl = costh(id) * dwkdx(ip,is1) + sinth(id) * dwkdy(ip,is1)
        kc_1abl = costh(id) * dwkdx(ip,is2) + sinth(id) * dwkdy(ip,is2)
      
        res = ka_1abl - kb_1abl - kc_1abl

      if (ip == 157) write(*,*) 'ddelta_dx', res
        
      end function ddelta_dx
!**********************************************************************
!*                                                                    *
!**********************************************************************
      function w(ip, is, is1, is2, n1, emf) result(res)
      use datapool, only : rkind, wk, cg, spsig, g9, one
      implicit none
     
      integer, intent(in)        :: ip, is, is1, is2, n1, emf
      real(rkind)                :: res !return

!!! mapping var.
      real(rkind)                :: omegaa  
      real(rkind)                :: omegab 
      real(rkind)                :: omegac 
      real(rkind)                :: cga  
      real(rkind)                :: cgb  
      real(rkind)                :: cgc  
      real(rkind)                :: ka  
      real(rkind)                :: kb  
      real(rkind)                :: kc 
!!! function declaration
      real(rkind)                :: tau
      
      real(rkind)                :: a,b,c,d, e
      
!       n1=0 means +
!       n1=1 means -
!       em=0 without eldeberky & madsen
!       em=1 with eldeberky & madsen

        omegaa = spsig(is)
        omegab = spsig(is)
        omegac = spsig(is)

        cga    = cg(ip,is)
        cgc    = cg(ip,is1)
        cgb    = cg(ip,is2)

        ka      = wk(ip,is)
        kb      = wk(ip,is1)
        kc      = wk(ip,is2)

!AR: reduce the +- stuff ...
!
        res =  one / (8 * sqrt(cga * cgb *cgc))  *  &        
     &         ( (-1)**n1 * (2 - emf * tau(ip, is, is1, is2, n1)) * kb * kc +  & 
     &         ( 1 - emf * tau(ip, is, is1, is2, n1)) * ((omegab * omegac)**2 / g9**2) + & 
     &         ( kb**2 * omegac / omegaa) + ( (-1)**n1 * kc**2 * omegab / omegaa) + ( (-1)**(n1+1) * & 
     &         ( (-1)**(n1+1) * ( 1- emf * tau(ip, is, is1, is2,n1)) * ((omegaa**2 * omegab * omegac) / g9**2) ) ) )

        a  = one / (8 * sqrt(cga * cgb *cgc))
        b  = (-1)**n1 * (2 - emf * tau(ip, is, is1, is2, n1)) * kb * kc
        c  = ( 1 - emf * tau(ip, is, is1, is2, n1)) * ((omegab * omegac)**2 / g9**2)
        d  = ( kb**2 * omegac / omegaa) + ( (-1)**n1 * kc**2 * omegab / omegaa)  
        e  = ( (-1)**(n1+1) * ( (-1)**(n1+1) * ( 1- emf * tau(ip, is, is1, is2,n1)) * ((omegaa**2 * omegab * omegac) / g9**2) )) 

        if (ip == 157) write(*,'(A20,7F15.10)') 'w    ', res, tau(ip, is, is1, is2,n1), a, b, c, d, e
        
      end function w
!**********************************************************************
!*                                                                    *
!**********************************************************************
      function tau(ip, is, is1, is2, n1) result(res)
      use datapool, only : rkind, wk, cg, spsig
      implicit none

      integer, intent(in)        :: ip, is, is1, is2, n1
      real(rkind)                :: res !return

      real(rkind)                :: cga 
      real(rkind)                :: ka  
      real(rkind)                :: kb  
      real(rkind)                :: kc  
      real(rkind)                :: omega_a  

        omega_a = spsig(is)

        cga    = cg(ip,is)

        ka      = wk(ip,is)
        kb      = wk(ip,is1)
        kc      = wk(ip,is2)
      
        res = 2 * cga * ( ka + (-1)**(n1+1) * kb - kc ) / omega_a

      end function tau
!**********************************************************************
!*                                                                    *
!**********************************************************************
      function dwdx(ip, is, is1, is2, id, n1, switch) result(res)
      use datapool, only : rkind, wk, cg, spsig, g9, one, costh, sinth, dcgdx, dcgdy, dwkdx, dwkdy, one, two, three
      implicit none
     
      integer, intent(in)        :: ip, is, is1, is2, id, n1
      real(rkind)                :: res !return

      integer                    :: switch

      real(rkind)                :: omegaa 
      real(rkind)                :: omegab 
      real(rkind)                :: omegac 
      real(rkind)                :: cga 
      real(rkind)                :: cgb 
      real(rkind)                :: cgc 
      real(rkind)                :: cga_
      real(rkind)                :: cgb_
      real(rkind)                :: cgc_
      real(rkind)                :: ka 
      real(rkind)                :: kb 
      real(rkind)                :: kc 
      real(rkind)                :: ka_ 
      real(rkind)                :: kb_ 
      real(rkind)                :: kc_ 

        omegaa = spsig(is)
        omegab = spsig(is)
        omegac = spsig(is)

        cga    = cg(ip,is)
        cgc    = cg(ip,is1)
        cgb    = cg(ip,is2)

        ka      = wk(ip,is)
        kb      = wk(ip,is1)
        kc      = wk(ip,is2)

        cga_ = costh(id) * dcgdx(ip,is) + sinth(id) * dcgdy(ip,is) 
        cgb_ = costh(id) * dcgdx(ip,is1) + sinth(id) * dcgdy(ip,is1)
        cgc_ = costh(id) * dcgdx(ip,is2) + sinth(id) * dcgdy(ip,is2)

        ka_ = costh(id) * dwkdx(ip,is) + sinth(id) * dwkdy(ip,is)
        kb_ = costh(id) * dwkdx(ip,is1) + sinth(id) * dwkdy(ip,is1)
        kc_ = costh(id) * dwkdx(ip,is2) + sinth(id) * dwkdy(ip,is2)

      
        res = (((4*Cga*Cgb*Cgc*((-1)**n1*((kb_*kc + kb*kc_)*omegaa +kc*kc_*omegab) +kb*kb_*omegac))/&
               omegaa -(Cga*Cgb_*Cgc +Cgb*(Cga_*Cgc + Cga*Cgc_))*(2*(-1)**n1*kb*kc -((-1)**n1*omegaa**2*omegab*omegac)/g9**2+&
              (omegab**2*omegac**2)/g9**2 +((-1)**n1*kc**2*omegab +kb**2*omegac)/omegaa)))/(16.*(Cga*Cgb*Cgc)**1.5)

        !if (ip == 157) write(*,*) 'dwdx', res 

      end function dwdx        
!**********************************************************************
!*                                                                    *
!**********************************************************************
      function k(ip, is, is1, is2, id, is3, is4, is5, n1, emf) result(res)
      use datapool, only : rkind, g9, one, stat
      implicit none

      integer, intent(in)        :: ip, is, is1, is2, is3, is4, is5, id, n1, emf
      real(rkind)                :: res ! return
      real(rkind)                :: delta, dwdx, w, ddelta_dx
      
       res = one / ((delta(ip, is3, is4, is5))**2) * dwdx(ip,is,is1,is2,id,n1,emf) - (w(ip,is,is1,is2,n1,emf) / ((delta(ip, is, is1, is2))**3)) * ddelta_dx(ip, is3, is4, is5, id)

      if (ip == 157) write(stat%fhndl,*) 'k', res

      end function k
!**********************************************************************
!*                                                                    *
!**********************************************************************
      function j(ip, is, is1, is2, id) result(res)
      use datapool, only : rkind, one, stat
      implicit none
     
      integer, intent(in)        :: ip, is, is1, is2, id 
      real(rkind)                :: res ! return 
      !!! function declaration 
      real(rkind)                :: delta, ddelta_dx
           
        res = - ( one / (delta(ip, is, is1, is2)**3)  ) * ddelta_dx(ip, is, is1, is2, id)

       if (ip == 157) write(stat%fhndl,*) 'j', res 
        
      end function j
!**********************************************************************
!*                                                                    *
!**********************************************************************
     subroutine snl3ta(ip,snl3,dsnl3)
     use datapool, only : rkind, msc, mdc, ac2, ZERO, spsig, cg, frintf, ddir, fr, stat
     implicit none

     real(rkind), intent(out) :: snl3(msc,mdc), dsnl3(msc,mdc)
     integer, intent(in)      :: ip
     integer :: is, is1, is2, id
     real(rkind) :: SUPER, SUB, f, f1, f2, SUPERD, SUBD, k, w, JAC

     do id = 1, mdc
       snl3(:,id) = zero
       dsnl3(:,id) = zero
#ifdef SNL3_DOUBLE
       do is = 1, msc
         SUPER = ZERO; SUPERD = ZERO
         SUB   = ZERO; SUBD = ZERO
         f = (ac2(ip,is,id) * spsig(is)**2. * frintf * ddir)* cg(ip,is)
         do is1 = 1, msc
           f1 = (ac2(ip,is1,id) * spsig(is1)**2. * frintf * ddir) * cg(ip,is1)
           do is2 = 1, msc
             f2 = (ac2(ip,is2,id) * spsig(is2)**2. * frintf * ddir) * cg(ip,is2)
             SUPER =   SUPER +  ( k(ip, is, is1, is2, id, is, is1, is2, 0, 0) * f1 * f2 &
     &                       +    k(ip, is1, is2, is, id, is, is1, is2, 1, 0) * f1 * f &
     &                       +    k(ip, is2, is1, is, id, is, is1, is2, 1, 0) * f2 * f ) * w(ip, is, is1, is2, 0, 0) * kron_delta(is,is1+is2)
             SUPERD = SUPERD +  ( k(ip, is1, is2, is, id, is, is1, is2, 1, 0) * f1 + &
     &                            k(ip, is2, is1, is, id, is, is1, is2, 1, 0) * f2 ) * w(ip, is, is1, is2, 0, 0) * kron_delta(is,is1+is2) 
           enddo
         enddo
         SUPER = 4 * SUPER
         do is1 = 1, msc
           f1 = ac2(ip,is1,id)
           do is2 = 1, msc
             f1 = ac2(ip,is2,id)
             SUB = SUB +   ( k(ip, is, is1, is2, id, is2, is, is1, 1, 0) * f1 * f2  +  k(ip, is1, is, is2, id, is2, is, is1, 1, 0) * f2 * f + &
      &                      k(ip, is2, is1, is, id, is2, is, is1, 0, 0) * f1 * f ) *  w(ip, is, is1, is2, 1, 0) * kron_delta(is2,is+is1)
             SUBD = SUBD + ( k(ip, is, is1, is2, id, is2, is, is1, 1, 0) * f1 * f2  +  k(ip, is1, is, is2, id, is2, is, is1, 1, 0) * f2 * f + &
      &                      k(ip, is2, is1, is, id, is2, is, is1, 0, 0) * f1 * f ) *  w(ip, is, is1, is2, 1, 0) * kron_delta(is2,is+is1)
           enddo
         enddo
#else
       do is = 2, msc 
         SUPER = ZERO; SUPERD = ZERO
         f = (ac2(ip,is,id) * spsig(is)**2. * frintf * ddir) * cg(ip,is)
         do is1 = 1, is - 1
           is2 = is - is1
           f1 = (ac2(ip,is1,id) * spsig(is1)**2. * frintf * ddir) * cg(ip,is1)
           f2 = (ac2(ip,is2,id) * spsig(is2)**2. * frintf * ddir) * cg(ip,is2)

           SUPER =   SUPER    + ( k(ip, is , is1, is2, id, is, is1, is2, 0, 0) * f1 * f2 &
     &                        +   k(ip, is1, is2, is , id, is, is1, is2, 1, 0) * f1 * f &
     &                        +   k(ip, is2, is1, is , id, is, is1, is2, 1, 0) * f2 * f ) * w(ip, is, is1, is2, 0, 0) 

           SUPERD = SUPERD    + ( k(ip, is1, is2, is , id, is, is1, is2, 1, 0) * f1  &
     &                        +   k(ip, is2, is1, is , id, is, is1, is2, 1, 0) * f2 ) *     w(ip, is, is1, is2, 0, 0)
         enddo
         SUPER = 4 * SUPER
         SUPERD = 4 * SUPERD
         JAC = 1./(SPSIG(IS)*DDIR*FRINTF)
         snl3(is,id) = SUPER * JAC
         dsnl3(is,id) = SUPERD
         if (ip == 157) write(*,'(2I10,A10,3F20.15)') is, id, 'super',fr(is),super,superd
       end do ! is 
       do is = 1, msc-1
         SUB   = ZERO; SUBD = ZERO
         f = (ac2(ip,is,id) * spsig(is)**2. * frintf * ddir) * cg(ip,is)
         do is1 = 1, msc - is
           is2 = is + is1
           f1 = (ac2(ip,is1,id) * spsig(is1)**2. * frintf * ddir) * cg(ip,is1)
           f2 = (ac2(ip,is2,id) * spsig(is2)**2. * frintf * ddir) * cg(ip,is2)

             SUB =    SUB + ( k(ip, is , is1, is2, id, is2, is, is1, 1, 0) * f1 * f2 &
      &                   +   k(ip, is1, is , is2, id, is2, is, is1, 1, 0) * f2 * f + &
      &                       k(ip, is2, is1, is , id, is2, is, is1, 0, 0) * f1 * f ) * w(ip, is, is1, is2, 1, 0) 

             SUBD = SUBD + (  k(ip, is , is1, is2, id, is2, is, is1, 1, 0) * f1 * f2  &
      &                  +    k(ip, is1, is , is2, id, is2, is, is1, 1, 0) * f2 * f  &
      &                  +    k(ip, is2, is1, is , id, is2, is, is1, 0, 0) * f1 * f ) * w(ip, is, is1, is2, 1, 0) 
         enddo ! is1 
         SUB = 8 * SUB
         SUBD = 8 * SUBD 
         if (ip == 157) write(*,'(2I10,A10,3F20.15)') is,id,'sub',fr(is),sub,subd
         JAC = 1./(SPSIG(IS)*DDIR*FRINTF)
         snl3(is,id) = snl3(is,id) + SUB * JAC
         dsnl3(is,id) = dsnl3(is,id) + SUBD
#endif
       enddo ! is 
     enddo ! id 

     if (ip == 157) then
       do is = 1, msc
         write(stat%fhndl,'(2i10,3f15.10)') ip, is, fr(is), sum(snl3(is,:)) * DDIR
       enddo
     endif

     end subroutine snl3ta
!**********************************************************************
!*                                                                    *
!**********************************************************************
      function kron_delta(i, j) result(res)
      use datapool, only : stat
      implicit none

      integer, intent(in)        :: i,j
      integer                    :: res ! return 

      res = int((float(i+j)-abs(i-j)))/(float((i+j)+abs(i-j))) 

      write(stat%fhndl,*) 'kron-delta', i, j, res

      end function kron_delta 
!**********************************************************************
!*                                                                    *
!**********************************************************************
      subroutine snl31 (ip, hs, smespc, acloc, imatra, imatda, ssnl3)
!
      use datapool
      implicit none

      integer, intent(in)        :: ip
      real(rkind), intent(in)    :: hs, smespc
      real(rkind), intent(out)   :: ssnl3(msc,mdc)
      real(rkind), intent(in)    :: acloc(msc,mdc)
      real(rkind), intent(inout) :: imatra(msc,mdc), imatda(msc,mdc)

      INTEGER :: IS, ID, I1, I2
      REAL(rkind)    :: BIPH, ASINB, XISTRI
      REAL(rkind)    :: ED(MSC)
      INTEGER :: IS1, IS2, ISMAX
      INTEGER :: IRES
      REAL(rkind)    :: AUX1, AUX2, FAC
      REAL(rkind)    :: JAC, DNL3IS1, DNL3IS2
      REAL(rkind)    :: NL3IS1, NL3IS2
      REAL(rkind)    :: cgl(msc),cl(msc),wkl(msc), AA
      REAL(rkind)    :: URSELL, BB, DEP1, DEP2, DEP3
      REAL(rkind)    :: SF3P(MSC,MDC), SF3M(MSC,MDC)

      PTRIAD(1)  = 1. 
      PTRIAD(2)  = 2.5
      PTRIAD(3)  = 10.
      PTRIAD(4)  = 0.2
      PTRIAD(5)  = 0.01

      BB   = ONE/15._rkind
      AA   = TWO/FIVE
      DEP1 = DEP(IP)
      DEP2 = DEP1**2
      DEP3 = DEP2*DEP1

      wkl = wk(ip,:)
      cgl = cg(ip,:)
      cl  = wc(ip,:)

      IF (TRICO .GT. 0.)  PTRIAD(1) = TRICO
      IF (TRIRA .GT. 0.)  PTRIAD(2) = TRIRA
      IF (TRIURS .GT. 0.) PTRIAD(5) = TRIURS

      CALL URSELL_NUMBER(HS,SMESPC,DEP(IP),URSELL)

      IF (URSELL > PTRIAD(5)) THEN
        I2     = INT(DBLE(MSC)/TWO)
        I1     = I2-1
        XISTRI = SPSIG(I2) / SPSIG(I1)
        IRES   = NINT(LOG(TWO)/LOG(XISTRI))
        ISMAX  = 1
        DO IS = 1, MSC
          IF (SPSIG(IS) < (PTRIAD(2)*SMESPC)) ISMAX = IS
        END DO
        ISMAX = MAX(ISMAX,IRES+1)
        BIPH  = PI/TWO*(MyTANH(0.2/URSELL)-ONE)
        ASINB = ABS(SIN(BIPH))
        DO ID = 1, MDC
          ED = ACLOC(:,ID)*SPSIG
          DO IS = 1, ISMAX-IRES
            IS1     = IS+IRES
            IS2     = IS
            AUX1    = WKL(IS2)**2*(G9*DEP1+TWO*CL(IS2)**2)
            AUX2    = WKL(IS1)*DEP1*(G9*DEP1+TWO*BB*G9*DEP3*WKL(IS1)**2-AA*SIGPOW(IS,2)*DEP2)
            JAC     = AUX1/AUX2
            FAC     = PTRIAD(1)*PI2*CL(IS1)*CGL(IS1)*JAC**2*ASINB
            NL3IS1  = MAX(ZERO,FAC*(ED(IS2)*ED(IS2)-TWO*ED(IS2)*ED(IS1)))
            DNL3IS1 = MAX(ZERO,FAC*(ED(IS2)-TWO*ED(IS1)))
            NL3IS2  = TWO*NL3IS1
            DNL3IS2 = TWO*DNL3IS1
            IF (ICOMP .LT. 2) THEN
              IMATRA(IS1,ID) = IMATRA(IS1,ID) + NL3IS1 / SPSIG(IS1)
              IMATRA(IS2,ID) = IMATRA(IS2,ID) - NL3IS2 / SPSIG(IS2)
              IMATDA(IS1,ID) = IMATDA(IS1,ID) + DNL3IS1 
              IMATDA(IS2,ID) = IMATDA(IS2,ID) - DNL3IS2 
            ELSE
              IMATRA(IS1,ID) = IMATRA(IS1,ID) + NL3IS1 / SPSIG(IS1)
              IMATRA(IS2,ID) = IMATRA(IS2,ID) - NL3IS2 / SPSIG(IS2)
              IMATDA(IS1,ID) = IMATDA(IS1,ID) - DNL3IS1 
              IMATDA(IS2,ID) = IMATDA(IS2,ID) + DNL3IS2 
            ENDIF
            SF3P(IS1,ID) = SF3P(IS1,ID) + NL3IS1
            SF3M(IS2,ID) = SF3M(IS2,ID) - NL3IS2
          END DO
        END DO
        SSNL3 = SF3P+SF3M
      END IF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      subroutine snl32 (ip, hs, smespc, acloc, imatra, imatda, ssnl3)
!
      use datapool
      implicit none

      integer, intent(in)        :: ip
      real(rkind), intent(in)    :: hs, smespc
      real(rkind), intent(out)   :: ssnl3(msc,mdc)
      real(rkind), intent(in)    :: acloc(msc,mdc)
      real(rkind), intent(inout) :: imatra(msc,mdc), imatda(msc,mdc)
      
      INTEGER           :: I, J, I1, I2
      INTEGER           :: IRES, ISMAX
      INTEGER           :: IS, ID
      REAL(rkind)              :: PTTRIAD(5)
      REAL(rkind)              :: DEP_2, DEP_3, BB, BIPH
      REAL(rkind)              :: TMN, URS
      REAL(rkind)              :: WiT,WjT,WNi,WNj,CGi,CGj,JACi,JACj,XISTRI
      REAL(rkind)              :: ALPH, BETA, SINBPH, FT, RHV_i, RHV_j, DIA_i, DIA_j
      REAL(rkind)              :: E(MSC), URSELL

      LOGICAL  :: TRIEXP

      PTTRIAD(1)  = 0.1 
      PTTRIAD(2)  = 2.2 
      PTTRIAD(3)  = 10. 
      PTTRIAD(4)  = 0.2
      PTTRIAD(5)  = 0.01

      IF (TRICO .GT. 0.)  PTTRIAD(1) = TRICO
      IF (TRIRA .GT. 0.)  PTTRIAD(2) = TRIRA
      IF (TRIURS .GT. 0.) PTTRIAD(5) = TRIURS

      CALL URSELL_NUMBER(HS,SMESPC,DEP(IP),URSELL)

      DEP_2 = DEP(IP)**2
      DEP_3 = DEP(IP)**3
      BB     = 1. / 15.
      
      TRIEXP = .FALSE.

      I2     = INT (FLOAT(MSC) / 2.)
      I1     = I2 - 1
      XISTRI = SPSIG(I2) / SPSIG(I1)
      IRES   = NINT ( LOG( 2.) / LOG ( XISTRI ) )
!
      ISMAX = 1
      DO IS = 1, MSC
        IF ( SPSIG(IS) .LT. ( PTTRIAD(2) * SMESPC) ) THEN
          ISMAX = IS
        ENDIF
      ENDDO 
!
      ISMAX = MAX ( ISMAX , IRES + 1 )
!      
      TMN = PI2 / SMESPC
      URS = MIN ( URSELL , TEN )
      IF (URS .LT. PTTRIAD(5)) THEN
        BIPH = 0.0
      ELSE
        BIPH = - PI2 / 8. * ( LOG10(URS) + 1.)
      ENDIF

      SINBPH = ABS( SIN(BIPH) )
      IF ( URS .GE. PTTRIAD(5) ) THEN
        DO ID = 1, MDC
          DO IS = 1, MSC
            E(IS)  = ACLOC(IS,ID) * PI2 * SPSIG(IS)
          ENDDO
          DO I = 1, ISMAX-IRES
            J     = I + IRES
            WiT   = SPSIG(I)
            WjT   = SPSIG(J)
            WNi   = WK(I,IP)
            WNj   = WK(J,IP)
            CGi   = CG(I,IP)
            CGj   = CG(J,IP)
            JACi  = PI2 * WiT
            JACj  = PI2 * WjT
            ALPH = 4. * WNi**2 * ( 0.5 + ( WiT**2 / ( WNi**2 * G9 * DEP(IP) ) ) )
            BETA  = -2. * WNj * ( G9 * DEP(IP) + 2. * BB * G9 * DEP_3 * WNj**2 - ( BB + 1./3. ) * WjT**2 * DEP_2   )
            FT    = PTTRIAD(1) * CGj * SINBPH * ( G9 * ALPH / BETA )**2
            RHV_i = 0.
            RHV_j = 0.
            DIA_i = 0.
            DIA_j = 0.
            IF ( TRIEXP ) THEN
              RHV_j = FT * ( (WjT/WNj) * E(i) * E(i) - 2. * (WiT/WNi) * E(j) * E(i) )
              RHV_i = RHV_j
              IF ( RHV_j .LE. 0. ) THEN
                RHV_j = 0.
                RHV_i = 0.
              ENDIF
              IMATRA(i,ID) = IMATRA(i,ID) - RHV_i / JACi
              IMATRA(j,ID) = IMATRA(j,ID) + 0.5 * RHV_j / JACj
              ssnl3(i,ID) = 0.5 * RHV_j / JACj
            ELSE
              RHV_j = FT * ((WjT/WNj) * E(i) * E(i) - 2. * (WiT/WNi) * E(j) * E(i))
              DIA_i = FT * ((WjT/WNj) * E(i)        - 2. * (WiT/WNi) * E(j))
              IF ( RHV_j .LE. 0. ) THEN
                RHV_j = 0.
                RHV_i = 0.
                DIA_j = 0.
                DIA_i = 0.
              ENDIF
              IMATRA(i,ID) = IMATRA(i,ID) - 0.5 * RHV_j / JACj
              IMATDA(j,ID) = IMATDA(j,ID) - DIA_i
              ssnl3(i,ID) = 0.5 * RHV_j / JACj
            ENDIF
          ENDDO
        ENDDO
      ENDIF
      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      subroutine snl33 (ip, hs, smespc, acloc, imatra, imatda, ssnl3)
!
      use datapool
      implicit none

      integer, intent(in)        :: ip
      real(rkind), intent(in)    :: hs, smespc
      real(rkind), intent(out)   :: ssnl3(msc,mdc)
      real(rkind), intent(in)    :: acloc(msc,mdc)
      real(rkind), intent(inout) :: imatra(msc,mdc), imatda(msc,mdc)

      INTEGER           :: I, J, I1, I2
      INTEGER           :: IRES, ISMAX
      INTEGER           :: IS, ID
      REAL(rkind)              :: PTTRIAD(5)
      REAL(rkind)              :: DEP_2, DEP_3, BIPH, RINT
      REAL(rkind)              :: AUX1, AUX2
      REAL(rkind)              :: WiT,WjT,WNi,WNj,CGi,CGj,JACi,JACj,XISTRI,Ci,Cj
      REAL(rkind)              :: SINBPH, FT, RHV_i, RHV_j, DIA_i, DIA_j
      REAL(rkind)              :: E(MSC), URSELL

      LOGICAL  :: TRIEXP

      PTTRIAD(1)  = 0.25
      PTTRIAD(2)  = 2.5
      PTTRIAD(3)  = 10.
      PTTRIAD(4)  = 0.2
      PTTRIAD(5)  = 0.01

      DEP_2 = DEP(IP)
      DEP_3 = DEP(IP)

      TRIEXP = .FALSE.

      IF (TRICO .GT. 0.)  PTTRIAD(1) = TRICO
      IF (TRIRA .GT. 0.)  PTTRIAD(2) = TRIRA
      IF (TRIURS .GT. 0.) PTTRIAD(5) = TRIURS

      CALL URSELL_NUMBER(HS,SMESPC,DEP(IP),URSELL)

      E(:) = 0.
      I2     = INT (FLOAT(MSC) / 2.)
      I1     = I2 - 1
      XISTRI = SPSIG(I2) / SPSIG(I1)
      IRES   = NINT ( LOG( 2.) / LOG ( XISTRI ) )
      ISMAX = 1

      DO IS = 1, MSC
       IF ( SPSIG(IS) .LT. ( PTTRIAD(2) * SMESPC) ) THEN
          ISMAX = IS
        ENDIF
      ENDDO

      ISMAX = MAX ( ISMAX , IRES + 1 )

      IF ( URSELL .GE. PTTRIAD(5) ) THEN

        BIPH   = (0.5*PI)*(MyTANH(PTTRIAD(4)/URSELL)-1)
        SINBPH = ABS( SIN(BIPH) )
        DO ID = 1, MDC
          DO IS = 1, MSC
            E(IS)  = ACLOC(IS,ID) * PI2 * SPSIG(IS)
          ENDDO
          DO I = 1, ISMAX-IRES
            J   = I + IRES
            WiT  = SPSIG(I)
            WjT  = SPSIG(J)
            WNi = WK(I,IP)
            WNj = WK(J,IP)
            Ci  = WiT / WNi
            Cj  = WjT / WNj
            CGi = CG(i,IP)
            CGj = CG(j,IP)
            JACi = PI2 * WiT
            JACj = PI2 * WjT
            AUX1 = WNi**2 * ( G9 * DEP(IP) + 2. * Ci**2 )
            AUX2 = WNj * DEP(IP) * ( G9 * DEP(IP) + (2./15.) * G9 * DEP_3 * WNj**2 - (2./ 5.) * WjT**2 * DEP_2 )
            RINT = AUX1 / AUX2
            FT = PTTRIAD(1) * Cj * CGj * RINT**2 * SINBPH
            RHV_i = 0.
            RHV_j = 0.
            DIA_i = 0.
            DIA_j = 0.
            IF ( TRIEXP ) THEN
              RHV_j = FT * (  E(i) * E(i)  - 2. * E(j) * E(i) )
              RHV_i = RHV_j
              IF ( RHV_j .LE. 0. ) THEN
                RHV_j = 0.
                RHV_i = 0.
              ENDIF
              IMATRA(i,ID) = IMATRA(i,ID) - 2. * RHV_i / JACi
              IMATRA(j,ID) = IMATRA(j,ID) + RHV_j / JACj
              ssnl3(i,id) = - 2. * RHV_i / JACi
            ELSE
              RHV_j = FT * ( E(i) * E(i) - 2. * E(j) * E(i) )
              DIA_i = FT * ( E(i) - 2. * E(j) )
              IF ( RHV_j .LE. 0. ) THEN
                RHV_j = 0.
                RHV_i = 0.
                DIA_j = 0.
                DIA_i = 0.
              ENDIF
              IMATRA(j,ID) = IMATRA(j,ID) + RHV_j / JACj
              IMATDA(i,ID) = IMATDA(i,ID) - 2. * DIA_i
              ssnl3(j,id) = RHV_j / JACj
            ENDIF
          ENDDO
        ENDDO
      ENDIF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      subroutine triadswan (ip, hs, smespc, acloc, imatra, imatda, ssnl3)
!
      use datapool
      implicit none

      integer, intent(in)        :: ip
      real(rkind), intent(in)    :: hs, smespc
      real(rkind), intent(out)   :: ssnl3(msc,mdc)
      real(rkind), intent(in)    :: acloc(msc,mdc)
      real(rkind), intent(inout) :: imatra(msc,mdc), imatda(msc,mdc)

      INTEGER I1, I2, ID, IS, ISM, ISM1, ISMAX, ISP, ISP1, IJ1, IJ2, IRES
      REAL(rkind)    AUX1, AUX2, BIPH, C0, CM, DEP_2, DEP_3, E0, ACTMP(MDC,MSC), URSELL
      REAL(rkind)    EM,FT, RINT, SIGPI, SINBPH, STRI, WISM, WISM1, FAC1, FACSCL, FACRES, SIGLOW, IMATDATMP(MDC,MSC)
      REAL(rkind)    WISP, WISP1,W0, WM, WN0, WNM,  XISLN, IMATRATMP(MDC,MSC)

      REAL(rkind) :: E(MSC)
      REAL(rkind), ALLOCATABLE :: SA(:,:)

      PTRIAD(1)  = 0.25
      PTRIAD(2)  = 2.5
      PTRIAD(3)  = 10.
      PTRIAD(4)  = 0.2
      PTRIAD(5)  = 0.01

      IF (TRICO .GT. 0.) PTRIAD(1) = TRICO
      IF (TRIRA .GT. 0.) PTRIAD(2) = TRIRA
      IF (TRIURS .GT. 0.) PTRIAD(5) = TRIURS

      CALL URSELL_NUMBER(HS,SMESPC,DEP(IP),URSELL)

      IMATRATMP = 0.
      IMATDATMP = 0.

      ACTMP = 0.
      DO IS = 1, MSC
        DO ID = 1, MDC
          ACTMP(ID,IS) = ACLOC(IS,ID)
        END DO
      END DO

      IJ2    = INT (FLOAT(MSC) / 2.)
      IJ1    = IJ2 - 1
      FAC1   = SPSIG(IJ2) / SPSIG(IJ1)
      IRES   = NINT ( LOG10( 2.) / LOG10( FAC1 ) )
      FACSCL = SPSIG(MSC)/SPSIG(MSC-IRES)

      IF (ABS(FACSCL-2.).GT.0.05) THEN
         FACRES = 10.**( LOG10(2.) / FLOAT(IRES) )
         SIGLOW   = SPSIG(MSC) / ( FACRES**(FLOAT(MSC-1) ) )
!         WRITE (*,*) 'CHECK RESOLUTION', IRES, FACSCL, FACRES, SIGLOW
      END IF
      DEP_2 = DEP(IP)**2
      DEP_3 = DEP(IP)**3
      I2     = INT (FLOAT(MSC) / 2.)
      I1     = I2 - 1
      XIS    = SPSIG(I2) / SPSIG(I1)
      XISLN  = LOG( XIS )
      ISP    = INT( LOG(2.) / XISLN )
      ISP1   = ISP + 1
      WISP   = (2. - XIS**ISP) / (XIS**ISP1 - XIS**ISP)
      WISP1  = 1. - WISP
      ISM    = INT( LOG(0.5) / XISLN )
      ISM1   = ISM - 1
      WISM   = (XIS**ISM -0.5) / (XIS**ISM - XIS**ISM1)
      WISM1  = 1. - WISM

!      IF (IP == 1786) THEN
!        WRITE(*,*) DEP_2, DEP_3, I2, I1, XIS, XISLN, ISP, ISP1, WISP, WISP1, ISM, ISM1, WISM, WISM1
!        WRITE(*,*) SUM(IMATRA), SUM(IMATDA), SUM(SSNL3)
!      ENDIF

      ALLOCATE (SA(1:MDC,1:MSC+ISP1))
      E  = 0.
      SA = 0.
      ISMAX = 1
      DO IS = 1, MSC
       IF ( SPSIG(IS) .LT. ( PTRIAD(2) * SMESPC) ) THEN
          ISMAX = IS
        ENDIF
      ENDDO

      ISMAX = MAX ( ISMAX , ISP1 )
      ISMAX = MIN( MSC, MAX ( ISMAX , ISP1 ) ) ! added fix the bug described below ...

      IF ( URSELL .GT. PTRIAD(5) ) THEN
        BIPH   = (0.5*PI)*(MyTANH(PTRIAD(4)/URSELL)-1.)
        SINBPH = ABS( SIN(BIPH) )
        DO ID = 1, MDC
           DO IS = 1, MSC
              E(IS) = ACTMP(ID,IS) * 2. * PI * SPSIG(IS)
           END DO
           DO IS = 1, ISMAX
              E0  = E(IS)
              W0  = SPSIG(IS)
              WN0 = WK(IS,IP)
              C0  = W0 / WN0
              IF ( IS.GT.-ISM1 ) THEN
                 EM  = WISM * E(IS+ISM1)      + WISM1 * E(IS+ISM)
                 WM  = WISM * SPSIG(IS+ISM1)  + WISM1 * SPSIG(IS+ISM)
                 WNM = WISM * WK(IS+ISM1,IP)  + WISM1 * WK(IS+ISM,IP)
                 CM  = WM / WNM
              ELSE
                 EM  = 0.
                 WM  = 0.
                 WNM = 0.
                 CM  = 0.
              END IF
              AUX1 = WNM**2 * ( G9 * DEP(IP) + 2.*CM**2 )
              AUX2 = WN0 * DEP(IP) * ( G9 * DEP(IP) + (2./15.) * G9 * DEP_3 * WN0**2 -(2./ 5.) * W0**2 * DEP_2 )
              RINT = AUX1 / AUX2
              FT = PTRIAD(1) * C0 * CG(IS,IP) * RINT**2 * SINBPH
              SA(ID,IS) = MAX(ZERO, FT * ( EM * EM - 2. * EM * E0 ))
           END DO
        END DO
        IF (ICOMP .LT. 2) THEN
          DO IS = 1, MSC
             SIGPI = SPSIG(IS) * 2. * PI
             DO ID = 1, MDC
                STRI = SA(ID,IS) - 2.*(WISP  * SA(ID,IS+ISP1) + WISP1 * SA(ID,IS+ISP ))
                ssnl3(is,id) = STRI
                IMATRA(IS,ID) =  IMATRA(IS,ID) + STRI / SIGPI
                IMATDA(IS,ID) =  IMATDA(IS,ID) + STRI / MAX(1.E-18_rkind, ACLOC(IS,ID)*SIGPI)
             END DO
          END DO
        ELSE
          DO IS = 1, MSC
             SIGPI = SPSIG(IS) * 2. * PI
             DO ID = 1, MDC
                STRI = SA(ID,IS) - 2.*(WISP  * SA(ID,IS+ISP1) + WISP1 * SA(ID,IS+ISP ))
                ssnl3(is,id) = STRI
                IF (STRI .GT. ZERO) THEN
                  IMATRA(IS,ID) =  IMATRA(IS,ID) + STRI / SIGPI
                ELSE
                  IMATDA(IS,ID) =  IMATDA(IS,ID) - STRI / MAX(1.E-18_rkind, ACLOC(IS,ID)*SIGPI)
                ENDIF
             END DO
          END DO
        ENDIF
      END IF

      DEALLOCATE (SA)

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
      subroutine triad_polnikov (ip, hs, smespc, acloc, imatra, imatda, ssnl3)
!
      use datapool
      implicit none

      integer, intent(in)        :: ip
      real(rkind), intent(in)    :: hs, smespc
      real(rkind), intent(out)   :: ssnl3(msc,mdc)
      real(rkind), intent(in)    :: acloc(msc,mdc)
      real(rkind), intent(inout) :: imatra(msc,mdc), imatda(msc,mdc)
!
      INTEGER               :: J, IS, ID, IT, ITER
      REAL*8                :: KI, KJ , ETOT
      REAL*8                :: PROPFAK, DBETA, DSIGMA
      REAL*8                :: BETA_0(MSC), BETA_1(MSC), SNL3(MSC,MDC)
      REAL*8                :: DELTAK(MSC), DK, H, TMP1, TMP2, TMP3
      REAL*8                :: PART1(MSC), PART2(MSC), EPS(MSC)
      REAL*8                :: SUMAC, KJKIKJ, KJKJKI

      ITER = 10

      ETOT = hs**2/FOUR

      DELTAK(1) = WK(IP,1)
      DO IS = 2, MSC
        DELTAK(IS) = WK(IS,IP) - WK(IS-1,IP)
      END DO

      H  = DEP(IP)

      IF (ETOT .GT. THR) THEN

        DO ID = 1, MDC
          DO IS = 1, MSC

            BETA_0(IS) = 0.1 * SPSIG(IS)/PI2
            KI = WK(IS,IP)
            DK = DELTAK(IS)
            PROPFAK = 9. * KI * SQRT(G9) * DK / (32. * PI * SQRT(H))

            DO IT = 1, ITER

              BETA_1(IS) = 0.
              PART1(IS)  = 0.
              PART2(IS)  = 0.

              DO J = 1, IS - 1
                SUMAC   = DBLE(( ACLOC(J,ID) * CG(J,IP) - ACLOC(IS-J,ID) * CG(IS-J,IP) ))
                KJ = WK(J,IP)
                KJKIKJ = KJ*(KI-KJ)
                KJKJKI = KJ*(KJ-KI)
                DSIGMA  = 0.5 * SQRT(G9*H**5.) * ABS(KI*KJKIKJ)
                DBETA   = BETA_0(IS) / ( PI * DSIGMA**2. + BETA_0(IS)**2.  )
                PART1(IS) = PART1(IS) +  KJKIKJ * DBETA * SUMAC
                IF (PART1(IS) .NE. PART1(IS)) WRITE (*,*) 'PART1', PART1(IS), KJKIKJ, DBETA, DSIGMA, BETA_0(IS)
              END DO
              DO J = IS+1, MSC
                SUMAC   = DBLE(( ACLOC(J-IS,ID) * CG(J-IS,IP) - ACLOC(J,ID) * CG(J,IP) ))
                KJ = WK(J,IP)
                KJKIKJ = KJ*(KI-KJ)
                KJKJKI = KJ*(KJ-KI)
                DSIGMA  = 0.5 * SQRT(G9*H**5.) * ABS(KI*KJKIKJ)
                DBETA   = BETA_0(IS) / ( PI * DSIGMA**2. + BETA_0(IS)**2.  )
                PART2(IS) = PART2(IS) +  KJKJKI * DBETA * SUMAC
                IF (PART2(IS) .NE. PART2(IS)) WRITE (*,*)'PART2', PART2(IS), KJKIKJ, DBETA, DSIGMA, BETA_0(IS)
              END DO
              BETA_1(IS) = PROPFAK * ( PART1(IS) -  2. * PART2(IS) )
              IF (ABS(BETA_1(IS)) .GT. THR) THEN
                EPS(IS) =  ( ABS(   BETA_0(IS) - BETA_1(IS)  ) ) / ABS(BETA_0(IS))
              ELSE
                EPS(IS) = 0.
              END IF
              IF (ABS(EPS(IS)) .LT.  10E-3) THEN
                EXIT
              END IF
              BETA_0(IS) = BETA_1(IS)
            END DO ! ... IT

            PART1(IS)  = 0.
            PART2(IS)  = 0.
            DO J = 1, IS - 1
              IF (BETA_0(IS) .LT. THR) CYCLE
              KJ = WK(J,IP)
              KJKIKJ = KJ*(KI-KJ)
              KJKJKI = KJ*(KJ-KI)
              DSIGMA  = 0.5 * SQRT(G9*H**5.) * ABS(KI*KJKIKJ)
              DBETA   = BETA_0(IS) / ( PI * DSIGMA**2. + BETA_0(IS)**2.  )
              TMP1    = DBLE(ACLOC(J,ID)*CG(J,IP)*ACLOC(IS-J,ID)*CG(IS-J,IP))
              TMP2    = DBLE(ACLOC(J,ID)*CG(J,IP)+ACLOC(IS-J,ID)*CG(IS-J,IP))
              TMP3    = (TMP1 - ACLOC(IS,ID)*CG(IS,IP)*TMP2)
              PART1(IS) = PART1(IS) +  KJKIKJ * DBETA * TMP3
              IF (PART1(IS) .NE. PART1(IS)) WRITE (*,*) 'PART1', PART1(IS), KJKIKJ, DBETA, DSIGMA, BETA_0(IS)
            END DO
            DO J = IS+1, MSC
              IF (BETA_0(IS) .LT. THR) CYCLE
              KJ = WK(J,IP)
              KJKIKJ = KJ*(KI-KJ)
              KJKJKI = KJ*(KJ-KI)
              DSIGMA  = 0.5 * SQRT(G9*H**5.) * ABS(KI*KJKIKJ)
              DBETA   = BETA_0(IS) / ( PI * DSIGMA**2. + BETA_0(IS)**2.  )
              TMP1    = DBLE(ACLOC(IS,ID)*CG(IS,IP)*ACLOC(J-IS,ID)*CG(J-IS,IP))
              TMP2    = DBLE(ACLOC(IS,ID)*CG(IS,IP)+ACLOC(J-IS,ID)*CG(J-IS,IP))
              TMP3    = (TMP1 - ACLOC(J,ID)*CG(J,IP)*TMP2)
              PART2(IS) = PART2(IS) +  KJKIKJ * DBETA * TMP3
              IF (PART2(IS) .NE. PART2(IS)) WRITE (*,*)'PART2', PART2(IS), KJKIKJ, DBETA, DSIGMA, BETA_0(IS)
            END DO
            SNL3(IS,ID) = PROPFAK * ( PART1(IS) - 2. *  PART2(IS) )
          END DO ! ... IS
        END DO ! ... ID

        DO IS = 1, MSC
          DO ID = 1, MDC
            IMATRA(IS,ID) = IMATRA(IS,ID) + SNL3(IS,ID)
            ssnl3(is,id) = SNL3(IS,ID)
            IF (ACLOC(IS,ID) .GT. 0.) THEN
              IMATDA(IS,ID) = IMATDA(IS,ID) + SNL3(IS,ID)/ACLOC(IS,ID)
            ELSE
              IMATDA(IS,ID) = 0.
            END IF
          END DO
        END DO

      END IF

      END SUBROUTINE
!**********************************************************************
!*                                                                    *
!**********************************************************************
