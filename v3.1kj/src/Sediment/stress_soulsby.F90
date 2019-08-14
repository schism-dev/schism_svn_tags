!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	subroutine	stress_soulsby(cosv,sinv,dhx,dhy,alphas,sed_angle,  &
     &			       tauc0,tauc,i)

!-----------------------------------------------------------------------
!			    ,
!	AUTHOR:		Andre Fortunato - 96-10-22 modified by
!                       Ligia Pinto - 09-09-15
!	PURPOSE:	
!
!-----------------------------------------------------------------------
!
!ZYL
        use elfe_glbl, only: errmsg
        use elfe_msgp, only: parallel_abort        

	real*8		dhx,dhy
	real*8		cosv,sinv,alphas,tauc
        real*8		cosas,sinas,cosat,tan_alphat,aux,aux1
	real*8		tauc0,sed_angle,sinat,aux2
        integer         i

!c-----------------------------------------------------------------------

!	alphas = datan(-(dhx*cosv+dhy*sinv)) !angulo alpha(s)
! Limitar o inclinacao a 0.9*(sed_angle) ou seja o angulo a atan(sed_angle*0.9)
!        aux=min(abs(alphas),datan(sed_angle*0.9))
!        alphas = sign(aux,alphas)
	if (abs(alphas)>datan(sed_angle))alphas=datan(sed_angle)*SIGN(1.0,alphas)
        cosas = dcos(alphas)
	sinas = dsin(alphas)
	tan_alphat = cosas*(dhx*sinv-dhy*cosv)
	cosat = dcos(datan(tan_alphat))
        sinat = dsin(datan(tan_alphat))
        tauc = tauc0*((sinas*cosat+dsqrt((cosas*cosas*sed_angle*sed_angle)-(sinat*sinat*sinas*sinas)))/(sed_angle))
        if (isnan(tauc)==-1) then
          write(errmsg,*)'tauc is NaN "i"',tauc
          call parallel_abort(errmsg)
        endif

	return
	end
