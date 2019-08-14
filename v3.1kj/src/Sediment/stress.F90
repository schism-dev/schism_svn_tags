!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine stress(cosv,sinv,dhx,dhy,alphas,sed_angle,tauc0,tauc,i)

!c-----------------------------------------------------------------------
!c      		    ,
!c      AUTHOR:		Andre Fortunato - 96-10-22
!c
!c      PURPOSE:	
!c
!c-----------------------------------------------------------------------

!ZYL
      use elfe_glbl, only: errmsg
      use elfe_msgp, only: parallel_abort
      implicit none
      real*8 dhx,dhy
      real*8 cosv,sinv,alphas,tauc
      real*8 cosas,sinas,cosat,tan_alphat,aux
      real*8 tauc0,sed_angle
      integer i

!c-----------------------------------------------------------------------

!      alphas = datan(-(dhx*cosv+dhy*sinv)) !alpha(s)
      if (abs(alphas)>datan(sed_angle))alphas=datan(sed_angle)*SIGN(1.0,alphas)
      cosas = dcos(alphas)
      sinas = dsin(alphas)
      tan_alphat = cosas*(dhx*sinv-dhy*cosv)
      cosat = dcos(datan(tan_alphat))
      aux = sed_angle*sed_angle-tan_alphat*tan_alphat
      if (aux > 0.d0) then
         tauc = tauc0*(cosas*cosat*dsqrt(1-((tan_alphat*tan_alphat)/(sed_angle*sed_angle)))+(sinas/sed_angle))
      else
         tauc=tauc0*(sinas/sed_angle)
      endif

!ZYL
       if (isnan(tauc)==-1) then
            write(errmsg,*)'tauc is NaN "i"',tauc,i,tauc0,sed_angle,cosas,cosat,tan_alphat,sinas
            call parallel_abort(errmsg)
       endif


      return
      end subroutine stress
