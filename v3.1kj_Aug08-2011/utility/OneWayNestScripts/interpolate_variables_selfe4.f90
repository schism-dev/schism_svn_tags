!
!********************************************************************************
!										
!	Interpolate 2 or 3D variables from SELFE results (format v5.0) to get *3D.th (binary).		
!       Use generic method to find parents.
!       Changes from interpolate_variables_selfe2: abandoned bucket sort; vgrid.fg can be SZ now.
!       Changes from interpolate_variables_selfe3: *3D.th now binary.
!										
!       Inputs: 
!              (1) bg.gr3 (for more precise x,y): hgrid.gr3 from large-domain run;
! 	       (2) fg.bp: boundary points (in build point format) that need *3D.th 
!                         for small-domain run; if multiple boundary segments
!                         have same type of b.c., you may lump all boundary points into
!                         this file (in proper order as bctides.in). 
!                         You may use gen_fg.f90 to generate this input;
!	       (3) vgrid.fg: vgrid.in from the small-domain run;
!              (4) interpolate_variables_selfe.in (see sample in this dir):
!                  ifile ndays - ifile=1: generate elev3D.th; =2: salt3D.th 
!                  and temp3D.th; =3: uv3D.th); ndays is the # of days needed;
!              (5) binary outputs from large-domain run (*.6?)
!       Outputs: 
!              (1) *3D.th for SELFE, depending on the choice in 
!                  interpolate_variables_selfe.in 
!              (2) fort.11: fatal errors; fort.12: non-fatal errors.
!
!       If you need to change the time step in these outputs, use timeint_3Dth2.f90 
!									
!       ifort -Bstatic -O3 -assume byterecl -o interpolate_variables_selfe4 interpolate_variables_selfe4.f90
!********************************************************************************
!
      program interpolate
      implicit real*8(a-h,o-z)
      parameter(mnp=80000) !max. # of nodes in fg/bg grids
      parameter(mne=160000) !max. # of elements fg/bg grids
      parameter(mns=200000) !max. # of sides fg/bg grids
      parameter(mnv=100) !max. # of vertical levels in fg/bg grids
      parameter(mnei=20) !max. # of neighboring elements for each node
!      parameter(mnxy=51)
!      parameter(mnelem=3000)
      parameter(nbyte=4)
      character*30 file63
      character*12 it_char
      character*48 start_time,version,variable_nm,variable_dim
      character*48 data_format
      dimension dp(mnp),kbp(mnp)
      dimension x(mnp),y(mnp),xfg(mnp),yfg(mnp),dpfg(mnp)
      dimension sigma_fg(mnv),cs(mnv),cs_fg(mnv),zfg(mnp,mnv),ratmin(mnp),ztot_fg(mnv)
      dimension nm(mne,4),i34(mne),area(mne),eta2(mnp),tnd(mnp,mnv),snd(mnp,mnv),tfg(mnp,mnv),sfg(mnp,mnv)
      !dimension ielem(mnxy,mnxy,mnelem),nelem(mnxy,mnxy),iparen(mnp),iparen_b(mnp),arco(mnp,4)
      dimension iparen(mnp),iparen_b(mnp),arco(mnp,4)
      dimension nxx(4,4,3),nm_fg(mne,4),i34_fg(mne),nne_fg(mnp),ine_fg(mnp,mnei),ic3_fg(mne,4)
      dimension isidenode(mns,2),js(mne,4),is(mns,2),xcj(mns),ycj(mns)
      dimension tsd(mns,mnv),ssd(mns,mnv),iglobal(mnp),etafg(mnp)
      dimension ze(mnv),z(mnp,mnv),icum(mnp,mnv,2)
      real h_0,h_s,h_c2,theta_b2,theta_f2,ztot(mnv),sigma2(mnv)
      real float1,float2
      
!..   Read hgrid.gr3 for more precise x,y
      open(14,file='bg.gr3',status='old')
      read(14,*)
      read(14,*)ne,np
      if(np.gt.mnp.or.ne.gt.mne) then
        write(*,*)'Too many nodes/elements',np,ne
        stop
      endif
      do i=1,np
        read(14,*)j,x(i),y(i),dp(i)
      enddo !i
      close(14)

!...  Read header from binary files
!...
      open(21,file='interpolate_variables_selfe.in',status='old')
      read(21,*)ifile,ndays
      close(21)
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
      read(63,rec=irec+2) float1 
      dtout=float1
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
      read(63,rec=irec+2) kz
      read(63,rec=irec+3) h_0
      read(63,rec=irec+4) h_s
      read(63,rec=irec+5) h_c2
      read(63,rec=irec+6) theta_b2
      read(63,rec=irec+7) theta_f2
      irec=irec+7
      if(nvrt.gt.mnv) then
        write(*,*)'Too many vertical levels',nvrt
        stop
      endif
      do k=1,kz-1
        read(63,rec=irec+k) ztot(k)
      enddo
      do k=kz,nvrt
        kin=k-kz+1
        read(63,rec=irec+k) sigma2(kin)
      enddo
      irec=irec+nvrt

!     Compute z-cor (not defined for initially dry nodes)
      do i=1,np
        if(dp(i)<=h_0) then
          kbp(i)=0 !dry
        else if(dp(i)<=h_s) then
          kbp(i)=kz
        else
          kbp(i)=0 !flag
          do k=1,kz-1
            if(-dp(i)>=ztot(k).and.-dp(i)<ztot(k+1)) then
              kbp(i)=k
              exit
            endif
          enddo !k
          if(kbp(i)==0) then
            write(*,*)'Cannot find a bottom level for node:',i
            stop
          endif
        endif
      enddo !i

      nsig=nvrt-kz+1
      do k=1,nsig
        cs(k)=(1-theta_b2)*sinh(theta_f2*sigma2(k))/sinh(theta_f2)+ &
     &theta_b2*(tanh(theta_f2*(sigma2(k)+0.5))-tanh(theta_f2*0.5))/2/tanh(theta_f2*0.5)
      enddo !k=1,nsig

      do i=1,np
        if(kbp(i)==0) cycle

!       Wet node
        do k=kz,nvrt
          kin=k-kz+1
          hmod2=min(dp(i),h_s)
          if(hmod2<=h_c2) then
            z(i,k)=sigma2(kin)*hmod2
          else
            z(i,k)=h_c*sigma2(kin)+(hmod2-h_c)*cs(kin)
          endif
        enddo !k

        z(i,kbp(i))=-dp(i)
        do k=kbp(i)+1,kz-1
          z(i,k)=ztot(k)
        enddo !k
      enddo !i


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
        read(63,rec=irec+4)kbp(m)
        irec=irec+4
      enddo !m=1,np
      do m=1,ne
        read(63,rec=irec+1)i34(m)
        if(i34(m)/=3) then
          write(*,*)'Quad'
          stop
        endif
        irec=irec+1
        do mm=1,i34(m)
          read(63,rec=irec+1)nm(m,mm)
          irec=irec+1
        enddo !mm

        n1=nm(m,1)
        n2=nm(m,2)
        n3=nm(m,3)
        area(m)=signa(x(n1),x(n2),x(n3),y(n1),y(n2),y(n3))
        if(area(m)<=0) then
          write(11,*)'Negative area at',m
          stop
        endif

      enddo !m
      irec0=irec

!     Read in fg.bp
      open(13,file='fg.bp',status='old')
      read(13,*)
      read(13,*)npfg
      if(npfg>mnp) then
        write(*,*)'Too many nodes in fg grid',npfg
        stop
      endif
      do i=1,npfg
        read(13,*)iglobal(i),xfg(i),yfg(i),dpfg(i)
        if(dpfg(i)<=0) then
          write(*,*)'Fg pt has negative depth:',i
          stop
        endif
      enddo !i
      close(13)

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
      do k=1,nsig
        cs_fg(k)=(1-theta_b_fg)*dsinh(theta_f_fg*sigma_fg(k))/dsinh(theta_f_fg)+ &
     &theta_b_fg*(dtanh(theta_f_fg*(sigma_fg(k)+0.5))-dtanh(theta_f_fg*0.5))/2/dtanh(theta_f_fg*0.5)
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
          n1=nm(ie,1); n2=nm(ie,2); n3=nm(ie,3)
          ar1=signa(xfg(i),x(n2),x(n3),yfg(i),y(n2),y(n3))  
          ar2=signa(x(n1),xfg(i),x(n3),y(n1),yfg(i),y(n3))
          ar3=signa(x(n1),x(n2),xfg(i),y(n1),y(n2),yfg(i))
          bb=dabs(ar1)+dabs(ar2)+dabs(ar3)
          aa=dabs(signa(x(n1),x(n2),x(n3),y(n1),y(n2),y(n3)))
          rat=dabs(bb-aa)/aa
          if(rat<ratmin(i)) then
            ratmin(i)=rat
            iparen_b(i)=ie
          endif
          if(rat<1.e-4) then
            iparen(i)=ie
            arco(i,1)=dmax1(0.d0,dmin1(1.d0,ar1/aa))
            arco(i,2)=dmax1(0.d0,dmin1(1.d0,ar2/aa))
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
        do k=max0(1,kbp(i)),nvrt
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
        open(16,file='elev3D.th',access='direct',recl=nrecl,status='replace')
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
      do iday=1,max0(1,ndays)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      write(it_char,'(i12)')iday
      open(63,file=it_char//'_'//file63,status='old',access='direct',recl=nbyte)
      if(ifile==2) open(65,file=it_char//'_salt.63',status='old',access='direct',recl=nbyte)

      irec=irec0
      do it1=1,nrec
        it_tot=it_tot+1
        read(63,rec=irec+1) float1
        time=it_tot*dtout !float1
        read(63,rec=irec+2) it
        irec=irec+2

        print*, 'time in days = ',time/86400

!        do i=1,np
!          read(63,rec=irec+i) float1
!        enddo !i
        irec=irec+np

        if(i23d==2) then
          do i=1,npfg
            ie=iparen(i)
            do j=1,3
              nd=nm(ie,j)
              read(63,rec=irec+nd) float1
              eta2(nd)=float1
            enddo !j
          enddo !i
          irec=irec+np

!         Do interpolation
          do i=1,npfg
            ie=iparen(i)
            iflag=0
            etafg(i)=0
            do j=1,i34(ie)
              nd=nm(ie,j)

              if(eta2(nd)+dp(nd)<h_0) iflag=1
              etafg(i)=etafg(i)+eta2(nd)*arco(i,j)
            enddo !j
!           Bg dry nodes
            if(iflag==1) write(12,*)'Warning: parent element has dry node:',i
!'
          enddo !i

!         Output
!          write(16,*)real(time)
!          do i=1,npfg
!            write(16,*)iglobal(i),real(etafg(i))
!          enddo !i
          write(16,rec=irec_out+1)real(time),(real(etafg(i)),i=1,npfg)
          irec_out=irec_out+1
        else !i23d=3 
          do i=1,npfg
            ie=iparen(i)
            do j=1,3
              nd=nm(ie,j)
              do k=max0(1,kbp(nd)),nvrt
                if(ifile==2) then !ST
                  read(63,rec=irec+icum(nd,k,1)) float1
                  read(65,rec=irec+icum(nd,k,1)) float2
                  tnd(nd,k)=float1
                  snd(nd,k)=float2
                else !uv
                  read(63,rec=irec+icum(nd,k,1)) float1
                  read(63,rec=irec+icum(nd,k,2)) float2
                  tnd(nd,k)=float1 !use for u
                  snd(nd,k)=float2 !for v
                endif
              enddo !k
            enddo !j
          enddo !i
          irec=irec+icum(np,nvrt,ivs)

!         Do interpolation (no true vertical interpolation)
!          print*, 'Day ',iday
          do i=1,npfg
            ie=iparen(i)
            kbe=max0(kbp(nm(ie,1)),kbp(nm(ie,2)),kbp(nm(ie,3)))
            if(kbe==0) then
              write(11,*)'All dry'
              stop
            endif
            do k=kbe,nvrt
              ze(k)=-1.e6
              do j=1,3
                nd=nm(ie,j)
                if(kbp(nd)/=0.and.z(nd,k)>ze(k)) ze(k)=z(nd,k)
              enddo !j
            enddo !k

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
              do j=1,i34(ie)
                nd=nm(ie,j)
                 
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
!              if(i==1) write(97,*)k,klev

            enddo !k=1,nvrt_fg
!            if(i==1) then
!              write(98,*)kbe,(ze(k),k=kbe,nvrt)
!              write(98,*)(zfg(i,k),k=1,nvrt_fg)
!            endif
          enddo !i=1,npfg

!         Output
!          write(17,*)real(time)
!          if(ifile==2) write(18,*)real(time)
!          do i=1,npfg
          if(ifile==2) then
            write(17,rec=irec_out+1)real(time),((real(tfg(i,k)),k=1,nvrt_fg),i=1,npfg)
            write(18,rec=irec_out+1)real(time),((real(sfg(i,k)),k=1,nvrt_fg),i=1,npfg)
          else !uv
            write(17,rec=irec_out+1)real(time),((real(tfg(i,k)),real(sfg(i,k)),k=1,nvrt_fg),i=1,npfg)
          endif
          irec_out=irec_out+1
!          enddo !i
        endif !i23d
      enddo !it1=1,nrec

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      enddo !iday=1,ndays

      stop
      end

      function signa(x1,x2,x3,y1,y2,y3)
!...  Compute signed area formed by pts 1,2,3
      implicit real*8(a-h,o-z)

      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2
      
      return
      end

