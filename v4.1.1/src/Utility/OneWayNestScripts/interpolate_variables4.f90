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

!
!********************************************************************************
!										
!	Interpolate 2 or 3D variables from SELFE results (format v5.0) to get *3D.th (binary).		
!       Use generic method to find parents.
!       Changes from interpolate_variables2: abandoned bucket sort; vgrid.fg can be SZ now.
!       Changes from interpolate_variables3: *3D.th now binary; added ivcor=1.
!										
!       Inputs: 
!          (1) bg.gr3 (for more precise x,y): hgrid.gr3 from large-domain run;
! 	   (2) fg.gr3: hgrid.gr3 from small-domain run; boundary segments 
!                     (specified in interpolate_variables.in) that need *3D.th 
!                      may be lumped;
!          (3) vgrid.bg: background vgrid.in - not needed if the large-domain run is 2D;
!	   (4) vgrid.fg: vgrid.in from the small-domain run; if 2D run, put 2 lines in this file:
!                        2\n 2\n;
!          (5) interpolate_variables.in (see sample in this dir):
!              1st line: ifile rndays - ifile=1: generate elev2D.th; =2: salt3D.th 
!                        and temp3D.th; =3: uv3D.th); rndays is the # of days needed;
!              2nd line: total # of open bnd seg. (that need *3D.th), list of segments
!          (6) binary outputs from large-domain run (*_elev.61 or *_temp.63 and *_salt.63, or *_hvel.64)
!       Outputs: 
!          (1) *3D.th for SELFE, depending on the choice in 
!                  interpolate_variables.in 
!          (2) fort.11: fatal errors; fort.12: non-fatal errors.
!
!       If you need to change the time step in these outputs, use timeint_3Dth2.f90 
!									
!       ifort -Bstatic -O3 -assume byterecl -CB -o interpolate_variables4 interpolate_variables4.f90 ~/SELFE/svn/trunk/src/Utility/Post-Processing-Fortran/compute_zcor.f90
!       pgf90 -O2 -mcmodel=medium -Mbounds -Bstatic -o interpolate_variables4 interpolate_variables4.f90 ~/SELFE/svn/trunk/src/Utility/Post-Processing-Fortran/compute_zcor.f90
!********************************************************************************
!
      program interpolate
!      implicit real*8(a-h,o-z)
      parameter(nbyte=4)
      character*30 file63
      character*12 it_char
      character*48 start_time,version,variable_nm,variable_dim
      character*48 data_format
      !bg arrays
      allocatable :: x(:),y(:),dp(:),kbp(:),kbp00(:),area(:),eta2(:),tnd(:,:),snd(:,:), &
     &ztot(:),sigma2(:),sigma_lcl(:,:),z(:,:),icum(:,:,:),ze(:)
      integer, allocatable :: elnode(:,:)
      !fg arrays: whole grid
      allocatable :: x_fg(:),y_fg(:),dp_fg(:),ztot_fg(:),sigma_fg(:),sigma_lcl_fg(:,:), &
     &kbp_fg(:),cs_fg(:)
      !fg arrays: open bnd part only
      allocatable :: xfg(:),yfg(:),dpfg(:),zfg(:,:),etafg(:),ratmin(:),kbpfg(:),imap(:)
      allocatable :: tfg(:,:),sfg(:,:),iob(:),nond(:),iond(:,:),iparen(:),iparen_b(:), &
     &arco(:,:)

      open(21,file='interpolate_variables.in',status='old')
      read(21,*)ifile,rndays
      read(21,*)nob
      allocate(iob(nob))
      rewind(21)
      read(21,*)
      read(21,*)nob,iob(:) !points to fg grid
      close(21)

!..   Read bg.gr3 for more precise x,y
      open(14,file='bg.gr3',status='old')
      read(14,*)
      read(14,*)ne,np
      allocate(x(np),y(np),dp(np),elnode(3,ne),kbp(np),kbp00(np),area(ne), &
     &eta2(np),stat=istat)
      if(istat/=0) stop 'Failed to alloc (1)'
      do i=1,np
        read(14,*)j,x(i),y(i),dp(i)
      enddo !i
      close(14)

!...  Read header from binary files
      if(ifile==1) then
        file63='elev.61'
      else if(ifile==2) then
        file63='temp.63'
      else !=3
        file63='hvel.64'
      endif
!      step_nu=900

      open(63,file='1_'//file63,status='old',access='direct',recl=nbyte)
      if(ifile==2) open(65,file='1_salt.63',status='old',access='direct',recl=nbyte)
      irec=0
      do m=1,48/nbyte
        read(63,rec=irec+m) data_format(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) version(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) start_time(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) variable_nm(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) variable_dim(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte

!      write(*,'(a48)')data_format
!      write(*,'(a48)')version
!      write(*,'(a48)')start_time
!      write(*,'(a48)')variable_nm
!      write(*,'(a48)')variable_dim

      read(63,rec=irec+1) nrec
      read(63,rec=irec+2) dtout
!      it_nu=step_nu/dtout
!      if(it_nu<1) then
!        write(*,*)'Nudging time too small:',it_nu
!        stop
!      endif
      read(63,rec=irec+3) nspool
      read(63,rec=irec+4) ivs
      read(63,rec=irec+5) i23d
!      read(63,rec=irec+6) float1 !vpos
      irec=irec+5

!     Vertical grid
      read(63,rec=irec+1) nvrt
!      read(63,rec=irec+2) kz
      read(63,rec=irec+3) h0
!      read(63,rec=irec+4) h_s
!      read(63,rec=irec+5) h_c2
!      read(63,rec=irec+6) theta_b2
!      read(63,rec=irec+7) theta_f2
!      irec=irec+7
!      if(nvrt.gt.mnv) then
!        write(*,*)'Too many vertical levels',nvrt
!        stop
!      endif
!      do k=1,kz-1
!        read(63,rec=irec+k) ztot(k)
!      enddo
!      do k=kz,nvrt
!        kin=k-kz+1
!        read(63,rec=irec+k) sigma2(kin)
!      enddo
!      irec=irec+nvrt
      irec=irec+nvrt+7

!     Horizontal grid
      read(63,rec=irec+1) np1
      read(63,rec=irec+2) ne1
      if(np1/=np.or.ne1/=ne) then
        print*, 'Mismatch in dimension between binary and hgrid.gr3'
        stop
      endif
      irec=irec+2
      do m=1,np
        read(63,rec=irec+1)float1
        read(63,rec=irec+2)float1
        read(63,rec=irec+3)float1
        read(63,rec=irec+4)kbp00(m)
        irec=irec+4
      enddo !m=1,np
      do m=1,ne
        read(63,rec=irec+1)i34
        if(i34/=3) then
          write(*,*)'Quad'
          stop
        endif
        irec=irec+1
        do mm=1,3 !i34(m)
          read(63,rec=irec+1)elnode(mm,m)
          irec=irec+1
        enddo !mm

        n1=elnode(1,m)
        n2=elnode(2,m)
        n3=elnode(3,m)
        area(m)=signa(x(n1),x(n2),x(n3),y(n1),y(n2),y(n3))
        if(area(m)<=0) then
          write(11,*)'Negative area at',m
          stop
        endif
      enddo !m
      irec0=irec

!     Read in vgrid.in from bg grid
      allocate(ztot(nvrt),sigma2(nvrt),sigma_lcl(nvrt,np),z(np,nvrt), &
     &tnd(np,nvrt),snd(np,nvrt),icum(np,nvrt,2),ze(nvrt),stat=istat)
      if(istat/=0) stop 'Failed to alloc (2)'

      if(nvrt==2) then !2D
        ivcor=2; kz=1; nsig=2; h_s=1.e6; h_c=h_s 
        theta_b=0; theta_f=1.e-4;
        do i=1,np
          if(dp(i)<=h0) then
            kbp(i)=0
          else
            kbp(i)=1
            z(i,1)=-dp(i)
            z(i,nvrt)=0
          endif
        enddo !i
      else !3D
        call get_vgrid('vgrid.bg',np,nvrt,ivcor,kz,h_s,h_c,theta_b,theta_f,ztot,sigma2,sigma_lcl,kbp)
        do i=1,np
          if(ivcor==2) then
            call zcor_SZ(dp(i),0.,h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot,sigma2,z(i,:),idry,kbp(i))
          else if(ivcor==1) then
            if(dp(i)<=h0) then
              kbp(i)=0
            else
              z(i,kbp(i):nvrt)=dp(i)*sigma_lcl(kbp(i):nvrt,i)
            endif
          else
            write(11,*)'Unknow ivcor:',ivcor
            stop
          endif
        enddo !i
      endif !2/3D

!     Read in fg.gr3
      open(13,file='fg.gr3',status='old')
      read(13,*)
      read(13,*)ne_fg,np_fg
      allocate(x_fg(np_fg),y_fg(np_fg),dp_fg(np_fg))
      do i=1,np_fg
        read(13,*)j,x_fg(i),y_fg(i),dp_fg(i)
!        if(dp_fg(i)<=0) then
!          write(*,*)'Fg pt has negative depth:',i
!          stop
!        endif
      enddo !i
      do i=1,ne_fg; read(13,*); enddo
      read(13,*)nope
      read(13,*)neta
      allocate(nond(nope))
      do i=1,nope
        read(13,*)nond(i)
        do j=1,nond(i)
          read(13,*) !iond
        enddo !j
      enddo !i
      rewind(13)
      mnope=maxval(nond)
      allocate(iond(nope,mnope))
      do i=1,np_fg+ne_fg+2+2; read(13,*); enddo
      do i=1,nope
        read(13,*) !nond(i)
        do j=1,nond(i)
          read(13,*) iond(i,j)
        enddo !j
      enddo !i
      close(13)

      npfg=sum(nond(iob(1:nob))) !# of b.c. nodes
      print*, '# of relevant open bnd nodes=',npfg

      allocate(imap(npfg),xfg(npfg),yfg(npfg),dpfg(npfg),etafg(npfg),ratmin(npfg), &
     &iparen(npfg),iparen_b(npfg),arco(npfg,3))

      npfg=0
      do i=1,nob
        ibnd=iob(i)
        do j=1,nond(ibnd)
          npfg=npfg+1
          nd=iond(ibnd,j)
          imap(npfg)=nd
          xfg(npfg)=x_fg(nd)
          yfg(npfg)=y_fg(nd)
          dpfg(npfg)=dp_fg(nd)
        enddo !j 
      enddo !i

!     Read in vgrid.fg 
      open(19,file='vgrid.fg',status='old')
      read(19,*); read(19,*)nvrt_fg
      close(19)

      allocate(ztot_fg(nvrt_fg),sigma_fg(nvrt_fg),sigma_lcl_fg(nvrt_fg,np_fg),kbp_fg(np_fg), &
     &zfg(npfg,nvrt_fg),kbpfg(npfg),tfg(npfg,nvrt_fg),sfg(npfg,nvrt_fg),cs_fg(nvrt_fg))

      if(nvrt_fg==2) then !2D
!        ivcor=2; kz_fg=1; nsig_fg=2; h_s_fg=1.e6; h_c_fg=h_s
!        theta_b_fg=0; theta_f_fg=1.e-4;
        do i=1,npfg
          if(dpfg(i)<=h0) then
            kbpfg(i)=0
          else
            kbpfg(i)=1
            zfg(i,1)=-dpfg(i)
            zfg(i,nvrt_fg)=0
          endif
        enddo !i
      else !3D
        call get_vgrid('vgrid.fg',np_fg,nvrt_fg,ivcor,kz_fg,h_s_fg,h_c_fg, &
     &theta_b_fg,theta_f_fg,ztot_fg,sigma_fg,sigma_lcl_fg,kbp_fg)

        !zcor
        do i=1,npfg
          if(ivcor==2) then !SZ
            call zcor_SZ(dpfg(i),0.,h0,h_s_fg,h_c_fg,theta_b_fg,theta_f_fg,kz_fg, &
     &nvrt_fg,ztot_fg,sigma_fg,zfg(i,:),idry2,kbpfg(i))
          else !=1
            if(dpfg(i)<=h0) then
              kbpfg(i)=0
            else
              nd=imap(i)
              kbpfg(i)=kbp_fg(nd)
              zfg(i,kbpfg(i):nvrt_fg)=dpfg(i)*sigma_lcl_fg(kbpfg(i):nvrt_fg,nd)
           endif
          endif !ivcor
          !Extend (not for dry nodes)
          zfg(i,1:kbpfg(i)-1)=zfg(i,kbpfg(i))
        enddo !i=1,np_fg
      endif !nvrt_fg

      if(1==2) then
!-----------------------------------
!     Read in vgrid.fg (SELFE)
      open(19,file='vgrid.fg',status='old')
      read(19,*) nvrt_fg,kz_fg,h_s_fg !kz>=1
      if(nvrt_fg>mnv.or.nvrt_fg<3) then
        write(*,*)'nvrt_fg > mnv or nvrt_fg<3'
        stop
      endif
      if(kz_fg<1.or.kz_fg>nvrt_fg-2) then
        write(*,*)'Wrong kz_fg:',kz_fg
        stop
      endif
      if(h_s_fg<10) then
        write(*,*)'h_s_fg needs to be larger:',h_s_fg
        stop
      endif

!     # of z-levels excluding "bottom" at h_s
      read(19,*) !for adding comment "Z levels"
      do k=1,kz_fg-1
        read(19,*)j,ztot_fg(k)
        if(k>1.and.ztot_fg(k)<=ztot_fg(k-1).or.ztot_fg(k)>=-h_s_fg) then
         write(*,*)'z-level inverted:',k
         stop
        endif
      enddo !k
      read(19,*) !level kz       
!     In case kz=1, there is only 1 ztot(1)=-h_s
      ztot_fg(kz_fg)=-h_s_fg

      nsig_fg=nvrt_fg-kz_fg+1 !# of S levels (including "bottom" & f.s.)
      read(19,*) !for adding comment "S levels"
      read(19,*)h_c_fg,theta_b_fg,theta_f_fg
      if(h_c_fg<5) then !large h_c to avoid 2nd type abnormaty
        write(*,*)'h_c_fg needs to be larger:',h_c_fg
        stop
      endif
      if(theta_b_fg<0.or.theta_b_fg>1) then
        write(*,*)'Wrong theta_b_fg:',theta_b_fg
        stop
      endif
      if(theta_f_fg<=0) then 
        write(*,*)'Wrong theta_f_fg:',theta_f_fg
        stop
      endif

      sigma_fg(1)=-1 !bottom
      sigma_fg(nsig_fg)=0 !surface
      read(19,*) !level kz_fg
      do k=kz_fg+1,nvrt_fg-1
        kin=k-kz_fg+1
        read(19,*) j,sigma_fg(kin)
        if(sigma_fg(kin)<=sigma_fg(kin-1).or.sigma_fg(kin)>=0) then
          write(*,*)'Check sigma levels at:',k,sigma_fg(kin),sigma_fg(kin-1)
          stop
        endif
      enddo
      read(19,*) !level nvrt
      close(19)

!     Compute C(s) and C'(s)
      do k=1,nsig_fg
        cs_fg(k)=(1-theta_b_fg)*sinh(theta_f_fg*sigma_fg(k))/sinh(theta_f_fg)+ &
     &theta_b_fg*(tanh(theta_f_fg*(sigma_fg(k)+0.5))-tanh(theta_f_fg*0.5))/2/tanh(theta_f_fg*0.5)
      enddo !k=1,nvrt
      close(19)

!...  Compute z-cor assuming eta=0
      do i=1,npfg
        do k=kz_fg,nvrt_fg
          kin=k-kz_fg+1
          hmod2=min(dpfg(i),h_s_fg) !dpfg(i)>0
          if(hmod2<=h_c_fg) then 
            zfg(i,k)=sigma_fg(kin)*dpfg(i)
          else
            zfg(i,k)=h_c_fg*sigma_fg(kin)+(hmod2-h_c_fg)*cs_fg(kin)
          endif
        enddo !k
        do k=1,kz_fg-1
          zfg(i,k)=ztot_fg(k)
        enddo !k
      enddo !i
!-----------------------------------
      endif !1==2

!...  Find parent elements
      arco=0
      iparen=0 !flags
      iparen_b=0 !flags; back-up parent based on min ratio
      ratmin=1.e25 !min. area ratio for fg pt
      loop1: do ie=1,ne
        iexit=1 !flag
        do i=1,npfg
          if(iparen(i)/=0) cycle

          iexit=0
          n1=elnode(1,ie); n2=elnode(2,ie); n3=elnode(3,ie)
          ar1=signa(xfg(i),x(n2),x(n3),yfg(i),y(n2),y(n3))  
          ar2=signa(x(n1),xfg(i),x(n3),y(n1),yfg(i),y(n3))
          ar3=signa(x(n1),x(n2),xfg(i),y(n1),y(n2),yfg(i))
          bb=abs(ar1)+abs(ar2)+abs(ar3)
          aa=abs(signa(x(n1),x(n2),x(n3),y(n1),y(n2),y(n3)))
          rat=abs(bb-aa)/aa
          if(rat<ratmin(i)) then
            ratmin(i)=rat
            iparen_b(i)=ie
          endif
          if(rat<1.e-4) then
            iparen(i)=ie
            arco(i,1)=max(0.,min(1.,ar1/aa))
            arco(i,2)=max(0.,min(1.,ar2/aa))
            if(arco(i,1)+arco(i,2)>1) then
              arco(i,3)=0
              arco(i,2)=1-arco(i,1)
            else
              arco(i,3)=1-arco(i,1)-arco(i,2)
            endif
          endif
        enddo !i=1,npfg
        if(iexit==1) exit
      end do loop1 !ie
      
!     Back-up parents
      icount=0
      write(12,*)'Following pts have no parents:'
      do i=1,npfg
        if(iparen(i)==0) then
          icount=icount+1
          if(iparen_b(i)==0) then
            write(11,*) 'No back-up:',i
            stop
          endif
          write(12,*)i,iparen_b(i)
          iparen(i)=iparen_b(i)
          arco(i,1)=1./3; arco(i,2)=1./3; arco(i,3)=1./3
        endif
      enddo !i
!      print*, icount,' pts have no immediate parent elements...:'

!...  Compute relative record # for a node and level for 3D outputs
!...
      icount=0
      do i=1,np
        do k=max0(1,kbp00(i)),nvrt
          do m=1,ivs
            icount=icount+1
            icum(i,k,m)=icount
          enddo !m
        enddo !k
      enddo !i=1,np

!...  Time iteration
!...
!     Open outputs
      if(ifile==1) then
        nrecl=nbyte*(1+npfg)
        open(16,file='elev2D.th',access='direct',recl=nrecl,status='replace')
      else if(ifile==2) then
        nrecl=nbyte*(1+npfg*nvrt_fg)
        open(17,file='temp3D.th',access='direct',recl=nrecl,status='replace')
        open(18,file='salt3D.th',access='direct',recl=nrecl,status='replace')
      else !uv
        nrecl=nbyte*(1+npfg*nvrt_fg*2)
        open(17,file='uv3D.th',access='direct',recl=nrecl,status='replace')
      endif
!'
      irec_out=0
      it_tot=0
      do iday=1,max(1,int(rndays*86400/dtout/nrec+0.1)) !1+int(rndays*86400/dtout/nrec+0.1)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      write(it_char,'(i12)')iday
      open(63,file=adjustl(it_char//'_'//file63),status='old',access='direct',recl=nbyte)
      if(ifile==2) open(65,file=it_char//'_salt.63',status='old',access='direct',recl=nbyte)

      irec=irec0
      do it1=1,nrec
        it_tot=it_tot+1
        read(63,rec=irec+1) float1
        time=it_tot*dtout !float1
        read(63,rec=irec+2) it
        irec=irec+2

        print*, 'processing time (days) = ',time/86400,' from stack:',iday

!        do i=1,np
!          read(63,rec=irec+i) float1
!        enddo !i
        irec=irec+np

        if(i23d==2) then
          do i=1,npfg
            ie=iparen(i)
            do j=1,3
              nd=elnode(j,ie)
              read(63,rec=irec+nd) eta2(nd)
            enddo !j
          enddo !i
          irec=irec+np

!         Do interpolation
          do i=1,npfg
            ie=iparen(i)
            iflag=0
            etafg(i)=0
            do j=1,3
              nd=elnode(j,ie)

              if(eta2(nd)+dp(nd)<h0) iflag=1
              etafg(i)=etafg(i)+eta2(nd)*arco(i,j)
            enddo !j
!           Bg dry nodes
            if(iflag==1) write(12,*)'Warning: parent element has dry node:',i
!'
          enddo !i

!         Output
!          write(16,*)real(time)
!          do i=1,npfg
!            write(16,*)i,real(etafg(i))
!          enddo !i
          write(16,rec=irec_out+1)time,(etafg(i),i=1,npfg)
          irec_out=irec_out+1
        else !i23d=3 
          do i=1,npfg
            ie=iparen(i)
            do j=1,3
              nd=elnode(j,ie)
              do k=max0(1,kbp00(nd)),nvrt
                if(ifile==2) then !ST
                  read(63,rec=irec+icum(nd,k,1)) tnd(nd,k)
                  read(65,rec=irec+icum(nd,k,1)) snd(nd,k)
                else !uv
                  read(63,rec=irec+icum(nd,k,1)) tnd(nd,k) !use for u
                  read(63,rec=irec+icum(nd,k,2)) snd(nd,k)
                endif
              enddo !k
            enddo !j
          enddo !i
          irec=irec+icum(np,nvrt,ivs)

!         Do interpolation (no true vertical interpolation)
!          print*, 'Day ',iday
          do i=1,npfg
            ie=iparen(i)
            kbe=maxval(kbp(elnode(:,ie)))
            if(kbe==0) then
              write(11,*)'All dry'
              stop
            endif
            do k=kbe,nvrt
              ze(k)=-1.e6
              do j=1,3
                nd=elnode(j,ie)
                if(kbp(nd)/=0.and.z(nd,k)>ze(k)) ze(k)=z(nd,k)
              enddo !j
            enddo !k

            !Debug
            !write(97,*)'Time=',time/86400,i,ie,kbe,ze(kbe:nvrt)

            do k=1,nvrt_fg
              if(zfg(i,k)<=ze(kbe)) then
                klev=kbe
              else if(zfg(i,k)>=ze(nvrt)) then
                klev=nvrt
              else
                klev=0
                do kk=kbe,nvrt-1
                  if(zfg(i,k)>=ze(kk).and.zfg(i,k)<=ze(kk+1)) then
                    klev=kk
                    exit
                  endif
                enddo !kk
                if(klev==0) then
                  write(11,*)'Cannot find a level:',i,k,zfg(i,k)
                  stop
                endif
              endif

              tfg(i,k)=0 !for u if ifile=3
              sfg(i,k)=0 !for v if ifile=3
              iflag=0
              tanc=-99; sanc=-99 !back-up values
              do j=1,3
                nd=elnode(j,ie)
                 
                if(ifile==2.and.tnd(nd,klev)<0) then
                  iflag=1
                else
                  tanc=tnd(nd,klev)
                  sanc=snd(nd,klev)
                endif
                tfg(i,k)=tfg(i,k)+tnd(nd,klev)*arco(i,j)
                sfg(i,k)=sfg(i,k)+snd(nd,klev)*arco(i,j)
              enddo !j
!             Bg dry nodes
              if(iflag==1) then !ifile=2
                if(tanc<0) then 
                  write(11,*)'Parent element dry:',i
                  stop
                endif
                tfg(i,k)=tanc
                sfg(i,k)=sanc
              endif

!             Debug
!              write(97,*)'bg:',k,klev,tnd(elnode(:,ie),klev),snd(elnode(:,ie),klev)
!              write(97,*)'fg:',arco(i,:),tfg(i,k),sfg(i,k)

            enddo !k=1,nvrt_fg
!            if(i==1) then
!              write(98,*)kbe,(ze(k),k=kbe,nvrt)
!              write(98,*)(zfg(i,k),k=1,nvrt_fg)
!            endif
          enddo !i=1,npfg

!         Output
          if(ifile==2) then
            write(17,rec=irec_out+1)time,((tfg(i,k),k=1,nvrt_fg),i=1,npfg)
            write(18,rec=irec_out+1)time,((sfg(i,k),k=1,nvrt_fg),i=1,npfg)
          else !uv
            write(17,rec=irec_out+1)time,((tfg(i,k),sfg(i,k),k=1,nvrt_fg),i=1,npfg)
            !Debug
            !do i=1,npfg
            !  write(99,*)'it1=',it1,xfg(i),yfg(i),dpfg(i),kbpfg(i),time,tfg(i,:),sfg(i,:)
            !enddo !i
          endif
          irec_out=irec_out+1
        endif !i23d
        print*, 'finished writing time (days) = ',time/86400
      enddo !it1=1,nrec

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      enddo !iday=1,

      stop
      end

!      function signa(x1,x2,x3,y1,y2,y3)
!!...  Compute signed area formed by pts 1,2,3
!!      implicit real*8(a-h,o-z)
!
!      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2
!      
!      return
!      end

