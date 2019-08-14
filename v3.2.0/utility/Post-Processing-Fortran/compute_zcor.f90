!     Routine to compute z coordinates for SELFE's SZ vertical system
!
!     ifort -Bstatic -O3 -o compute_zcor compute_zcor.f90

!     Driver
!      implicit real*8(a-h,o-z)
!      parameter(nvrt=11)
!      dimension ztot(nvrt),sigma(nvrt),zcor(nvrt)
!
!      h0=0.01; h_s=100; h_c=5; theta_b=1; theta_f=1.e-5
!      kz=3
!      ztot(1:kz-1)=(/-1000.,-500./)
!      ztot(kz)=-h_s
!      print*, 'Input eta, dp'
!      read*, eta, dp
!
!      nsig=nvrt-kz+1
!      do k=kz,nvrt
!        kin=k-kz+1
!        sigma(kin)=real(kin-1)/(nsig-1)-1
!      enddo !k
!      print*, 'sigma=',sigma(1:nsig)
!
!      call zcor_SZ(dp,eta,h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot,sigma,zcor,idry,kbp)
!      print*, 'idry, kbp =',idry,kbp
!      print*, 'zcor =',zcor(kbp:nvrt)
!
!      stop
!      end

!     Routine to compute z coordinates for SELFE's SZ vertical system
!     Inputs:
!             dp: depth;
!             eta: elevation;
!             h0: min. depth;
!             h_s: transition depth between S and Z layers;
!             h_c: transition depth between S and sigma
!             theta_b, theta_f: spacing const. in S coordinate system;
!             nvrt: total # of vertical levels (S+Z);
!             kz: # of Z levels (1 if pure S);
!             ztot(1:kz):  z coordinates for Z levels; note that ztot(kz)=-h_s;
!             sigma(1:nvrt-kz+1): sigma coordinates for S (or sigma) levels; note that sigma(1)=-1, sigma(nvrt-kz+1)=0;
!     Outputs:
!             idry: wet (0) or dry (1) flag;
!             kbp: bottom index (0 if dry);
!             zcor(kbp:nvrt): z coordinates (undefined if dry);    
      subroutine zcor_SZ(dp,eta,h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot,sigma,zcor,idry,kbp)
      implicit real*8(a-h,o-z)
      integer, intent(in) :: kz,nvrt
      real*8, intent(in) :: dp,eta,h0,h_s,h_c,theta_b,theta_f,ztot(nvrt),sigma(nvrt)
      integer, intent(out) :: idry,kbp
      real*8, intent(out) :: zcor(nvrt)

      real*8 :: cs(nvrt)

!     Sanity check
      if(nvrt<3) stop 'nvrt too small'
      if(kz<1.or.kz>nvrt-2) stop 'kz wrong'
      if(h_c<5.or.h_c>=h_s) then !large h_c to avoid 2nd type abnormaty
        print*, 'h_c needs to be larger:',h_c; stop
      endif
      if(theta_b<0.or.theta_b>1) then
        print*, 'Wrong theta_b:',theta_b; stop
      endif
      if(theta_f<=0) then
        print*, 'Wrong theta_f:',theta_f; stop
      endif

!     Pre-compute constants
      s_con1=sinh(theta_f)
      nsig=nvrt-kz+1 !# of S levels 
      do k=1,nsig
        cs(k)=(1-theta_b)*sinh(theta_f*sigma(k))/sinh(theta_f)+ &
     &theta_b*(tanh(theta_f*(sigma(k)+0.5))-tanh(theta_f*0.5))/2/tanh(theta_f*0.5)
      enddo !k

      if(eta<=h0-h_s) then
        stop 'Deep depth dry'
      else if(eta+dp<=h0) then
        idry=1; kbp=0
      else !wet
!       S-levels
        idry=0
        hmod=min(h_s,dp)
        do k=kz,nvrt
          kin=k-kz+1
          if(hmod<=h_c) then
            zcor(k)=sigma(kin)*(hmod+eta)+eta
          else if(eta<=-h_c-(hmod-h_c)*theta_f/s_con1) then !hmod(i)>h_c>=0
            print*, 'Pls choose a larger h_c (2):',eta,h_c
            stop
          else
            zcor(k)=eta*(1+sigma(kin))+h_c*sigma(kin)+(hmod-h_c)*cs(kin)
          endif
        enddo !k=kz,nvrt

!         z-levels
        if(dp<=h_s) then
          kbp=kz
        else !bottom index 
          kbp=0 !flag
          do k=1,kz-1
            if(-dp>=ztot(k).and.-dp<ztot(k+1)) then
              kbp=k
              exit
            endif
          enddo !k
          if(kbp==0) then
            print*, 'Cannot find a bottom level for node (3):',i
            stop
          endif

          if(kbp>=kz.or.kbp<1) then
            print*, 'Impossible 92:',kbp,kz
            stop
          endif
          zcor(kbp)=-dp
          do k=kbp+1,kz-1
            zcor(k)=ztot(k)
          enddo !k
        endif

        do k=kbp+1,nvrt
          if(zcor(k)-zcor(k-1)<=0) then
            write(*,*)'Inverted z-levels at:',i,k,zcor(k)-zcor(k-1),eta,hmod
            stop
          endif
        enddo !k
      endif !wet ot dry

      end subroutine zcor_SZ
