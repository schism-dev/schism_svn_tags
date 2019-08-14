!   Copyright 2014 College of William and Mary
!
!   Licensed under the Apache License, Version 2.0 (the "License");
!   you may not use this file except in compliance with the License.
!   You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
!   Unless required by applicable law or agreed to in writing, software
!   distributed under the License is distributed on an "AS IS" BASIS,
!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!   See the License for the specific language governing permissions and
!   limitations under the License.

!     Interpolate *3D.th onto new vgrid.in (no change in dt) for 3D variables only
!     Not up to date with ivcor=1.
!     Input: hgrid.gr3;
!            th.old (binary, time starting from dt); 
!            vint.in: 1st line: total # of days, h0 (min depth in param.in), ivs (1 for scalar; 2 for vector)
!                     2nd line: nond: # of open bnd nodes (lump all); 
!                     3rd and onwards: iond(1:nond) (list of # of open bnd nodes)
!            vgrid.in.0: old vgrid (if 2D, 1 line with '2')
!            vgrid.in: new vgrid (if 2D, 1 line with '2')
!     Output: th.new (binary, time starting from dt)

!     ifort -Bstatic -assume byterecl -O2 -o vint_3Dth vint_3Dth.f90
      program riverforcing
      parameter(nbyte=4)
      allocatable :: th0(:,:,:),th(:,:,:),iond(:)
      allocatable :: ztot(:),sigma(:),zcor(:,:),zcor0(:,:),kbp0(:)
      allocatable :: dp(:),kindx(:,:),zrat(:,:)
      dimension wild(10)

!      h0=0.05

      open(21,file='vint.in',status='old')
      read(21,*)rndays,h0,ivs
      if(ivs/=1.and.ivs/=2) stop 'wrong input ivs'
      read(21,*)nond0
      allocate(iond(nond0),stat=istat)
      do i=1,nond0
        read(21,*)iond(i)
      enddo !i
      close(21)

!     hgrid for depths
      open(14,file='hgrid.gr3',status='old')
      read(14,*); read(14,*)ne,np
      allocate(dp(np))
      do i=1,np 
        read(14,*)j,xtmp,ytmp,dp(i)
      enddo !i
      close(14)

!     Old vgrid
      open(19,file='vgrid.in.0',status='old')
      read(19,*) nvrt0
      allocate(kbp0(nond0),stat=istat)

      if(nvrt0==2) then
        kbp0=1
      else !3D
        rewind(19)
        read(19,*) nvrt0,kz,h_s
        allocate(ztot(nvrt0),sigma(nvrt0),zcor0(nvrt0,nond0),stat=istat)

        read(19,*) !for adding comment "Z levels"
        do k=1,kz-1
          read(19,*)j,ztot(k)
          if(k>1.and.ztot(k)<=ztot(k-1).or.ztot(k)>=-h_s) then
            write(*,*)'z-level inverted:',k
            stop
          endif
        enddo !k
        read(19,*) !level kz
!       In case kz=1, there is only 1 ztot(1)=-h_s
        ztot(kz)=-h_s

        nsig=nvrt0-kz+1 !# of S levels (including "bottom" & f.s.)
        read(19,*) !for adding comment "S levels"
        read(19,*)h_c,theta_b,theta_f

        sigma(1)=-1 !bottom
        sigma(nsig)=0 !surface
        read(19,*) !level kz
        do k=kz+1,nvrt0-1
          kin=k-kz+1
          read(19,*) j,sigma(kin)
          if(sigma(kin)<=sigma(kin-1).or.sigma(kin)>=0) then
            write(*,*)'Check sigma levels at:',k,sigma(kin),sigma(kin-1)
           stop
          endif
        enddo
        read(19,*) !level nvrt0
        close(19)

!       Compute zcor0
        do i=1,nond0
          nd=iond(i)
          call zcor_SZ(dp(nd),0.,h0,h_s,h_c,theta_b,theta_f,kz, &
     &nvrt0,ztot,sigma,zcor0(:,i),idry,kbp0(i))
          if(idry==1) stop 'Dry...'
        enddo !i
      endif !2D or 3D

!     New vgrid
      open(19,file='vgrid.in',status='old')
      read(19,*) nvrt

      if(nvrt>2.and.nvrt0>2) then !both 3D 
        rewind(19)
        read(19,*) nvrt,kz,h_s
        deallocate(ztot,sigma)
        allocate(ztot(nvrt),sigma(nvrt),zcor(nvrt,nond0), &
     &kindx(nvrt,nond0),zrat(nvrt,nond0),stat=istat)

        read(19,*) !for adding comment "Z levels"
        do k=1,kz-1
          read(19,*)j,ztot(k)
          if(k>1.and.ztot(k)<=ztot(k-1).or.ztot(k)>=-h_s) then
            write(*,*)'z-level inverted:',k
            stop
          endif
        enddo !k
        read(19,*) !level kz
!       In case kz=1, there is only 1 ztot(1)=-h_s
        ztot(kz)=-h_s

        nsig=nvrt-kz+1 !# of S levels (including "bottom" & f.s.)
        read(19,*) !for adding comment "S levels"
        read(19,*)h_c,theta_b,theta_f

        sigma(1)=-1 !bottom
        sigma(nsig)=0 !surface
        read(19,*) !level kz
        do k=kz+1,nvrt-1
          kin=k-kz+1
          read(19,*) j,sigma(kin)
          if(sigma(kin)<=sigma(kin-1).or.sigma(kin)>=0) then
            write(*,*)'Check sigma levels at:',k,sigma(kin),sigma(kin-1)
            stop
          endif
        enddo
        read(19,*) !level nvrt
        close(19)

!       Compute zcor
        do i=1,nond0
          nd=iond(i)
          call zcor_SZ(dp(nd),0.,h0,h_s,h_c,theta_b,theta_f,kz, &
     &nvrt,ztot,sigma,zcor(:,i),idry,kbp)
          if(idry==1) stop 'Dry...'

         !Compute kindx and zrat (distance from kindx)
         do k=kbp,nvrt
           if(zcor(k,i)<=zcor0(kbp0(i),i)) then
             kindx(k,i)=kbp0(i); zrat(k,i)=0
           else if(zcor(k,i)>=zcor0(nvrt0,i)) then
             kindx(k,i)=nvrt0-1; zrat(k,i)=1
           else
             kindx(k,i)=0
             do kk=kbp0(i),nvrt0-1
               if(zcor(k,i)>=zcor0(kk,i).and.zcor(k,i)<=zcor0(kk+1,i)) then
                 kindx(k,i)=kk
                 tmp=(zcor(k,i)-zcor0(kk,i))/(zcor0(kk+1,i)-zcor0(kk,i))
                 if(tmp<0.or.tmp>1) then
                   write(*,*)'Outbound:',i,k,zcor(k,i),tmp
                   stop
                 endif
                 zrat(k,i)=tmp
                 exit
               endif
             enddo !kk  
             if(kindx(k,i)==0) then
               write(*,*)'Failed to find a level:',i,nd,k,dp(nd),zcor(k,i), &
     &zcor0(kbp0(i):nvrt,i)
               stop
             endif
           endif !zcor
         enddo !k
       
         kindx(1:kbp-1,i)=kindx(kbp,i)
         zrat(1:kbp-1,i)=zrat(kbp,i)
        enddo !i=1,nond0
      endif !both 3D

      allocate(th0(ivs,nvrt0,nond0),th(ivs,nvrt,nond0),stat=istat)
      if(istat/=0) stop 'Failed too alloc.'

      irecl0=nbyte*(1+nond0*nvrt0*ivs)
      irecl=nbyte*(1+nond0*nvrt*ivs)
      open(17,file='th.old',access='direct',recl=irecl0,status='old')
      open(18,file='th.new',access='direct',recl=irecl,status='replace')
      read(17,rec=1)dt0,th0(:,:,:)
      nt0=rndays*86400/dt0+1.e-5

      do it=1,nt0
        read(17,rec=it)time,th0(:,:,:)
        print*, 'Time(days)=',time/86400

        if(nvrt0==2) then
          do i=1,nond0
            do m=1,ivs
              th(m,:,i)=th0(m,1,i)
            enddo !m
          enddo !i
        else if(nvrt==2) then !nvrt0/=2
          do i=1,nond0
            nd=iond(i)
            wild=0
            do k=kbp0(i),nvrt0-1
              wild(1:ivs)=wild(1:ivs)+(th0(1:ivs,k,i)+th0(1:ivs,k+1,i))/2*(zcor0(k+1,i)-zcor0(k,i))
            enddo !k
            if(dp(nd)<=h0) stop 'Dry open node'
            do m=1,ivs
              th(m,:,i)=wild(m)/dp(nd)
            enddo !m

            !Debug
            if(it==nt0) then
              write(99,*)'Node#, time:',i,time/86400
              write(99,*)'Old:',th0(:,:,i)
              write(99,*)'New:',th(:,:,i)
            endif
          enddo !i
        else !both 3D
          do i=1,nond0
            do k=1,nvrt
              kin=kindx(k,i)
              th(1:ivs,k,i)=th0(1:ivs,kin,i)*(1-zrat(k,i))+th0(1:ivs,kin+1,i)*zrat(k,i)
            enddo !k
          enddo !i
        endif !2D/3D

        write(18,rec=it)time,th(:,:,:)
      enddo !it

      stop
      end 

!!     Routine to compute z coordinates for SELFE's SZ vertical system
!!     Inputs:
!!             dp: depth;
!!             eta: elevation;
!!             h0: min. depth;
!!             h_s: transition depth between S and Z layers;
!!             h_c: transition depth between S and sigma
!!             theta_b, theta_f: spacing const. in S coordinate system;
!!             nvrt: total # of vertical levels (S+Z);
!!             kz: # of Z levels (1 if pure S);
!!             ztot(1:kz):  z coordinates for Z levels; note that ztot(kz)=-h_s;
!!             sigma(1:nvrt-kz+1): sigma coordinates for S (or sigma) levels; note that sigma(1)=-1, sigma(nvrt-kz+1)=0;
!!     Outputs:
!!             idry: wet (0) or dry (1) flag;
!!             kbp: bottom index (0 if dry);
!!             zcor(kbp:nvrt): z coordinates (undefined if dry);    
!      subroutine zcor_SZ(dp,eta,h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot,sigma,zcor,idry,kbp)
!      implicit real(a-h,o-z)
!      integer, intent(in) :: kz,nvrt
!      real, intent(in) :: dp,eta,h0,h_s,h_c,theta_b,theta_f,ztot(nvrt),sigma(nvrt)
!      integer, intent(out) :: idry,kbp
!      real, intent(out) :: zcor(nvrt)
!
!      real :: cs(nvrt)
!
!!     Sanity check
!      if(nvrt<3) stop 'nvrt too small'
!      if(kz<1.or.kz>nvrt-2) stop 'kz wrong'
!      if(h_c<5.or.h_c>=h_s) then !large h_c to avoid 2nd type abnormaty
!        print*, 'h_c needs to be larger:',h_c; stop
!      endif
!      if(theta_b<0.or.theta_b>1) then
!        print*, 'Wrong theta_b:',theta_b; stop
!      endif
!      if(theta_f<=0) then
!        print*, 'Wrong theta_f:',theta_f; stop
!      endif
!
!!     Pre-compute constants
!      s_con1=sinh(theta_f)
!      nsig=nvrt-kz+1 !# of S levels 
!      do k=1,nsig
!        cs(k)=(1-theta_b)*sinh(theta_f*sigma(k))/sinh(theta_f)+ &
!     &theta_b*(tanh(theta_f*(sigma(k)+0.5))-tanh(theta_f*0.5))/2/tanh(theta_f*0.5)
!      enddo !k
!
!      if(eta<=h0-h_s) then
!        stop 'Deep depth dry'
!      else if(eta+dp<=h0) then
!        idry=1; kbp=0
!      else !wet
!!       S-levels
!        idry=0
!        hmod=min(h_s,dp)
!        do k=kz,nvrt
!          kin=k-kz+1
!          if(hmod<=h_c) then
!            zcor(k)=sigma(kin)*(hmod+eta)+eta
!          else if(eta<=-h_c-(hmod-h_c)*theta_f/s_con1) then !hmod(i)>h_c>=0
!            print*, 'Pls choose a larger h_c (2):',eta,h_c
!            stop
!          else
!            zcor(k)=eta*(1+sigma(kin))+h_c*sigma(kin)+(hmod-h_c)*cs(kin)
!          endif
!        enddo !k=kz,nvrt
!
!!         z-levels
!        if(dp<=h_s) then
!          kbp=kz
!        else !bottom index 
!          kbp=0 !flag
!          do k=1,kz-1
!            if(-dp>=ztot(k).and.-dp<ztot(k+1)) then
!              kbp=k
!              exit
!            endif
!          enddo !k
!          if(kbp==0) then
!            print*, 'Cannot find a bottom level for node (3):',i
!            stop
!          endif
!
!          if(kbp>=kz.or.kbp<1) then
!            print*, 'Impossible 92:',kbp,kz
!            stop
!          endif
!          zcor(kbp)=-dp
!          do k=kbp+1,kz-1
!            zcor(k)=ztot(k)
!          enddo !k
!        endif
!
!        do k=kbp+1,nvrt
!          if(zcor(k)-zcor(k-1)<=0) then
!            write(*,*)'Inverted z-levels at:',i,k,zcor(k)-zcor(k-1),eta,hmod
!            stop
!          endif
!        enddo !k
!      endif !wet ot dry
!
!      end subroutine zcor_SZ
