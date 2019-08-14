!===============================================================================
!===============================================================================
! MORSELFE MISCELLANEOUS SUBROUTINES
!
! subroutine sed_settleveloc
! subroutine sed_taucrit
! subroutine sed_write_debug
!
!===============================================================================
!===============================================================================

      SUBROUTINE sed_settleveloc(ised)
!--------------------------------------------------------------------!
! This subroutine computes settling velocity from Soulsby (1997)     !
!                                                                    !
! Author: florian ganthy (fganthy@lnec.pt ; florian.ganthy@gmail.com)!
! Date:   21/12/2012                                                 !
!                                                                    !
! History:                                                           !
!                                                                    !
!--------------------------------------------------------------------!

      USE elfe_glbl, ONLY : rkind,errmsg
      USE elfe_msgp, ONLY: parallel_abort
      USE sed_mod,   ONLY : rhom,g,Srho,Sd50,Wsed

      IMPLICIT NONE
      SAVE

!- Local variables --------------------------------------------------!
      
      INTEGER,INTENT(IN)    :: ised

      REAL(rkind),PARAMETER :: nu=1.36d-6 ! Cinematic viscosity
      REAL(rkind)           :: ratio,dstar,tmp

!- Start Statement --------------------------------------------------!

! - Preliminary parameters
      ratio    = Srho(ised)/rhom
      dstar    = Sd50(ised) * (g*(ratio-1.d0)/nu**2.d0)**(1.d0/3d0)

! - Settling velocity (in m.s-1)
      tmp=10.36d0**2.d0+1.049d0*dstar**3.d0
      if(tmp<0) then
        write(errmsg,*)'SED3D; sed_misc_subs: tmp<0;',tmp
        call parallel_abort(errmsg)
      endif
      !Wsed(ised) = (nu/Sd50(ised)) * ((10.36d0**2.d0 + 1.049d0*      &
      !&                                dstar**3.d0)**0.5d0 - 10.36d0)
      Wsed(ised) = (nu/Sd50(ised)) * ( sqrt(tmp) - 10.36d0)
      if(Wsed(ised)<=0) then
        write(errmsg,*)'SED3D; sed_misc_subs: Wsed<=0; ',Wsed(ised)
        call parallel_abort(errmsg)
      endif

!--------------------------------------------------------------------!
      END SUBROUTINE sed_settleveloc    

!===============================================================================
!===============================================================================

      SUBROUTINE sed_taucrit(ised)
!--------------------------------------------------------------------!
! This subroutine computes critical shear stres for erosion from     !
! Soulsby (1997)                                                     !
!                                                                    !
! Author: florian ganthy (fganthy@lnec.pt ; florian.ganthy@gmail.com)!
! Date:   21/12/2012                                                 !
!                                                                    !
! History:                                                           !
!                                                                    !
!--------------------------------------------------------------------!

      USE elfe_glbl, ONLY : rkind,errmsg
      use elfe_msgp, only: parallel_abort
      USE sed_mod,   ONLY : rhom,g,Srho,Sd50,tau_ce

      IMPLICIT NONE
      SAVE

!- Local variables --------------------------------------------------!
      
      INTEGER,INTENT(IN)    :: ised

      REAL(rkind),PARAMETER :: nu=1.36d-6 ! Cinematic viscosity
      REAL(rkind)           :: ratio,dstar,theta_cr

!- Start Statement --------------------------------------------------!

! - Preliminary parameters
      ratio    = Srho(ised)/rhom
      dstar    = Sd50(ised) * (g*(ratio-1.d0)/nu**2.d0)**(1.d0/3d0)
      if(1.d0+1.2d0*dstar==0) call parallel_abort('SED3D; sed_misc_subs: dev. by 0')
      theta_cr = 0.3d0/(1.d0+1.2d0*dstar) + 0.055d0 *                &
      &         (1.d0-EXP(-0.02d0*dstar))

! - Critical shear stress (N.m-2)
      tau_ce(ised) = theta_cr*g*(Srho(ised)-rhom)*Sd50(ised)

      if(tau_ce(ised)<=0) then
        write(errmsg,*)'SED3D; sed_misc_subs: tau_ce(ised)<=0; ',tau_ce(ised)
        call parallel_abort(errmsg)      
      endif

!--------------------------------------------------------------------!
      END SUBROUTINE sed_taucrit

!===============================================================================
!===============================================================================

      SUBROUTINE sed_write_debug(it)
!--------------------------------------------------------------------!
! This subroutine writes debug or additional information to a file   !
! e.g. mirror.out                                                    !
!                                                                    !
! Author: Ligia Pinto                                                !
! Date: xx/xx/xxxx                                                   !
!                                                                    !
! History: 2012/12 - F.Ganthy : form homogenisation of sediments     !
!          routines                                                  !
!                                                                    !
!--------------------------------------------------------------------!

      USE sed_mod
      USE elfe_glbl, ONLY : nea,ntracers,dt
      USE elfe_msgp, ONLY : myrank

      IMPLICIT NONE

      INTEGER,INTENT(IN)     :: it
!- Local variables --------------------------------------------------!

      INTEGER :: debug_file,i,j,k

!- Start Statement --------------------------------------------------!

!     In outputs/nonfatal_0000
      debug_file = 12

      IF(myrank==0) THEN
        WRITE(debug_file,*)'SED: sed_write_debug'
        WRITE(debug_file,*)'it: ',it
        WRITE(debug_file,*)'dt: ',dt
        WRITE(debug_file,*)'nstp: ',nstp
        WRITE(debug_file,*)'nnew: ',nnew

        DO i=1,nea
            WRITE(debug_file,*)'nea:',i,' bustr:',real(bustr(i)),      &
     &            ' bvstr:',real(bvstr(i)),' tau_w:',real(tau_w(i))
        ENDDO !nea

        DO i=1,Nbed
          DO j=1,nea
            DO k=1,ntracers
              WRITE(debug_file,*)'Nbed:',i,      &
     &                ' nea:',j,' ntracers:',k,' bed_frac:',real(bed_frac(i,j,k)),   &
     &                ' bed_mass(1):',real(bed_mass(i,j,1,k)),                       &
     &                ' bed_mass(2):',real(bed_mass(i,j,2,k))
            ENDDO !k
          ENDDO !j
        ENDDO !i

        DO i=1,Nbed
          DO j=1,nea
            WRITE(debug_file,*) 'Nbed:',i,' nea:',j,                                       &
     &              ' bed_thick:',real(bed_thick(j)),' bed(ithck):',real(bed(i,j,ithck)),  &
     &              ' bed(iaged):',real(bed(i,j,iaged)),' bed(iporo):',real(bed(i,j,iporo))
          ENDDO !j
        ENDDO !i

        DO i=1,nea
          WRITE(debug_file,*)'nea:',i,' bottom(itauc):',real(bottom(i,itauc)),                  &
     &            ' bottom(isd50):',real(bottom(i,isd50)),' bottom(iwsed):',real(bottom(i,iwsed)),&
     &            ' bottom(idens):',real(bottom(i,idens)),' bottom(iactv):',real(bottom(i,iactv))
        ENDDO !i

!!           DO i=1,nea
!!             WRITE(debug_file,*)'nea:',i,                    &
!!     &              ' bottom(nea,itauc):',bottom(i,itauc),                     &
!!     &              ' bottom(nea,isd50):',bottom(i,isd50),                     &
!!     &              ' bottom(nea,iwsed):',bottom(i,iwsed),                     &
!!     &              ' bottom(nea,idens):',bottom(i,idens)
!!           ENDDO !nea

      ENDIF !myrank

!--------------------------------------------------------------------!
      END SUBROUTINE sed_write_debug

