      SUBROUTINE set_vbc(bustr,bvstr)
!!======================================================================
!! August, 2007                                                        !
!!========================================================  Ligia Pinto=
!!                                                                     !
!! This subroutine is from ROMS                                        !
!!                                                                     !
!================================================== Hernan G. Arango ===
!  Copyright (c) 2002-2007 The ROMS/TOMS Group                         !
!    Licensed under a MIT/X style license                              !
!    See License_ROMS.txt                                              !
!=======================================================================
!                                                                      !
!  This module sets vertical boundary conditons for tracers.           !
!  tracers.                                                            !
!                                                                      !
!=======================================================================
!
      USE sed_param
      USE ocean_mod
      USE elfe_glbl, only: nvrt,nea,dfv,idry_e,kbe,nm,uu2,vv2,ze
      use elfe_msgp, only: myrank,parallel_abort

      IMPLICIT NONE
      SAVE
!
!***********************************************************************
!
#ifndef BBL_MODEL
      real(r8) :: bustr(nea)
      real(r8) :: bvstr(nea)
#endif

!
!  Local variable declarations.
!
      integer :: i,j,n1,n2,n3
      real(r8) :: cff1, cff2, cff3, cff4, cff5,tmp

#if defined UV_LOGDRAG
      real(r8) :: wrk
#endif

!     cdb_max and  cdb_min defined in sed_param

#ifndef BBL_MODEL
!
!-----------------------------------------------------------------------
!  Set kinematic bottom momentum flux (m2/s2).
!-----------------------------------------------------------------------

#ifdef UV_LOGDRAG
!
!  Set logarithmic bottom stress.
!
      DO i=1,nea
        if(idry_e(i)==1) cycle

!        tmp=z_w(kbe(i),i)-z_w(kbe(i)-1,i)
        tmp=ze(kbe(i)+1,i)-ze(kbe(i),i)                
        if(tmp>Zob(i)) then
           if(tmp<=0.or.Zob(i)<=0) call parallel_abort('SEDIMENT: Cd failed')
!'
           cff1=1.0_r8/dlog(tmp/Zob(i))
           cff2=vonKar*vonKar*cff1*cff1
           wrk=MIN(Cdb_max,MAX(Cdb_min,cff2))
        endif

        n1=nm(i,1)
        n2=nm(i,2)
        n3=nm(i,3)
!        if((z_w(kbe(i),i)-z_w(kbe(i)-1,i))>Zob(i)) then   !LLP
        if((ze(kbe(i)+1,i)-ze(kbe(i),i))>Zob(i)) then      !LLP
          cff3=(vv2(kbe(i)+1,n1)+vv2(kbe(i)+1,n2)+vv2(kbe(i)+1,n3))/3
          cff4=(uu2(kbe(i)+1,n1)+uu2(kbe(i)+1,n2)+uu2(kbe(i)+1,n3))/3
          cff5=SQRT(cff4*cff4+cff3*cff3)
          bustr(i)=wrk*cff4*cff5 !bottom stress
          bvstr(i)=wrk*cff3*cff5
        else
          cff3=(vv2(kbe(i)+2,n1)+vv2(kbe(i)+2,n2)+vv2(kbe(i)+2,n3))/3- &
     &         (vv2(kbe(i)+1,n1)+vv2(kbe(i)+1,n2)+vv2(kbe(i)+1,n3))/3
          cff4=(uu2(kbe(i)+2,n1)+uu2(kbe(i)+2,n2)+uu2(kbe(i)+2,n3))/3- &
     &         (uu2(kbe(i)+1,n1)+uu2(kbe(i)+1,n2)+uu2(kbe(i)+1,n3))/3

          cff5=(dfv(n1,kbe(i)+1)+dfv(n2,kbe(i)+1)+dfv(n3,kbe(i)+1)+ &
     & dfv(n1,kbe(i)+2)+dfv(n2,kbe(i)+2)+dfv(n3,kbe(i)+2))/6
!          tmp=z_w(kbe(i)+1,i)-z_w(kbe(i),i)    !LLP
          tmp=ze(kbe(i)+2,i)-ze(kbe(i)+1,i)    !LLP
          if(tmp<=0) call parallel_abort('SEDIMENT: div. by 0')
          bustr(i)=cff5*cff4/tmp !bottom stress
          bvstr(i)=cff5*cff3/tmp
        endif
      END DO !i

#elif defined UV_QDRAG
!
!  Set quadratic bottom stress.
!
      DO i=1,nea
        if(idry_e(i)==1) cycle
        n1=nm(i,1)
        n2=nm(i,2)
        n3=nm(i,3)
        cff1=(vv2(kbe(i)+1,n1)+vv2(kbe(i)+1,n2)+vv2(kbe(i)+1,n3))/3
        cff2=(uu2(kbe(i)+1,n1)+uu2(kbe(i)+1,n2)+uu2(kbe(i)+1,n3))/3

        cff3=SQRT(cff2*cff2+cff1*cff1)
!Error: rdrg2 not defined
!LLP user defined, if defined valor must be read from ...
        bustr(i)=rdrg2*cff2*cff3
        bvstr(i)=rdrg2*cff1*cff3
      END DO

#elif defined UV_LDRAG
!
!  Set linear bottom stress.
!
      DO i=1,nea
        if(idry_e(i)==1) cycle
        n1=nm(i,1)
        n2=nm(i,2)
        n3=nm(i,3)
        cff1=(vv2(kbe(i)+1,n1)+vv2(kbe(i)+1,n2)+vv2(kbe(i)+1,n3)/3
        cff2=(uu2(kbe(i)+1,n1)+uu2(kbe(i)+1,n2)+uu2(kbe(i)+1,n3)/3
!Error: rdrg not defined
!LLP user defined 
        bustr(i)=rdrg*cff2
        bvstr(i)=rdrg*cff1
      END DO
#endif

!ifndef BBL_MODEL
#endif
      END SUBROUTINE set_vbc
