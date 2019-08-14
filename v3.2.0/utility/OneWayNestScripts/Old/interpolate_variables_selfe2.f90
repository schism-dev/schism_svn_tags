!
!********************************************************************************
!										
!	Interpolate 2 or 3D variables from SELFE results (format v5.0).		
!       Use bucket sort to find parents.					
!       Changes from interpolate_variables_selfe: add hvel.
!										
!       Input: binary files; 							
!              bg.gr3 (for more precise x,y);					
! 	       fg.bp: grid points that need interpolation from bg.gr3 		
!	       vgrid.fg: vgrid for fg.bp (SELFE);				
!              interpolate_variables_elcirc.in: ifiletype (1: elev; 2: S,T; 3: uv);
!                                               ndays;				
!       Output: ifile=1, fort.16 (elev); ifile=2, fort.17 and fort.18 (T,S); 
!               ifile=3, fort.17 (uv).
!									
!       Coompilation: ifort -Bstatic -O3 -assume byterecl -o interpolate_variables_selfe2 interpolate_variables_selfe2.f90
!********************************************************************************
!
      program interpolate
      implicit real*8(a-h,o-z)
      parameter(mnp=80000)
      parameter(mne=160000)
      parameter(mns=200000)
      parameter(mnv=100)
      parameter(mnei=20)
      parameter(mnxy=51)
      parameter(mnelem=3000)
      parameter(nbyte=4)
      character*30 file63
      character*12 it_char
      character*48 start_time,version,variable_nm,variable_dim
      character*48 data_format
      dimension dp(mnp),kbp(mnp) !,kfp(mnp)
      dimension x(mnp),y(mnp),xfg(mnp),yfg(mnp),dpfg(mnp),sigma(mnv),cs(mnv),dcs(mnv),zfg(mnp,mnv)
      dimension nm(mne,4),i34(mne),area(mne),eta2(mnp),tnd(mnp,mnv),snd(mnp,mnv),tfg(mnp,mnv),sfg(mnp,mnv)
      dimension ielem(mnxy,mnxy,mnelem),nelem(mnxy,mnxy),iparen(mnp),iparen_b(mnp),arco(mnp,4)
      dimension ixm(2),iym(2),xm(2),ym(2),xi(mnxy),yi(mnxy),k0(mnv)
      dimension nxx(4,4,3),nm_fg(mne,4),i34_fg(mne),nne_fg(mnp),ine_fg(mnp,mnei),ic3_fg(mne,4)
      dimension isidenode(mns,2),js(mne,4),is(mns,2),xcj(mns),ycj(mns)
      dimension tsd(mns,mnv),ssd(mns,mnv),iglobal(mnp),etafg(mnp)
      dimension ze(mnv),z(mnp,mnv)
      real h_0,h_s,h_c2,theta_b2,theta_f2,ztot(mnv),sigma2(mnv)
      real float1,float2
      
!...  Define bins for sorting elements
      nx=51; ny=31
      if(nx<=1.or.ny<=1) then
        write(*,*)'Increase nx, ny'
        stop
      endif
      if(nx>mnxy.or.ny>mnxy) then
        write(*,*)'Increase mnxy'
        stop
      endif
      xi(1)=-99041; xi(6)=250847; xi(21)=333820; xi(46)=372411; xi(51)=628640
      do i=2,5
        xi(i)=xi(1)+(i-1)*(xi(6)-xi(1))/5
      enddo !i
      do i=7,20
        xi(i)=xi(6)+(i-6)*(xi(21)-xi(6))/15
      enddo !i
      do i=22,45
        xi(i)=xi(21)+(i-21)*(xi(46)-xi(21))/25
      enddo !i
      do i=47,50
        xi(i)=xi(46)+(i-46)*(xi(51)-xi(46))/5
      enddo !i
      yi(1)=-987619; yi(6)=175366; yi(26)=365380; yi(31)=757954
      do i=2,5
        yi(i)=yi(1)+(i-1)*(yi(6)-yi(1))/5
      enddo !i
      do i=7,25
        yi(i)=yi(6)+(i-6)*(yi(26)-yi(6))/20
      enddo !i
      do i=27,30
        yi(i)=yi(26)+(i-26)*(yi(31)-yi(26))/5
      enddo !i

!...  End inputting bin lines

!     Check
      do i=1,nx-1
        if(xi(i)>=xi(i+1)) then
          write(*,*)'Inverted x lines:',i
          stop
        endif
      enddo !i
      do i=1,ny-1
        if(yi(i)>=yi(i+1)) then
          write(*,*)'Inverted y lines:',i
          stop
        endif
      enddo !i

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
      open(21,file='interpolate_variables_elcirc.in',status='old')
      read(21,*)ifiletype,ndays
      if(ifiletype==1) then
        file63='elev.61'
      else if(ifiletype==2) then
        file63='temp.63'
      else !=3
        file63='hvel.64'
      endif
!      step_nu=900

      open(63,file='1_'//file63,status='old',access='direct',recl=nbyte)
      open(65,file='1_salt.63',status='old',access='direct',recl=nbyte)
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

!...  Find parent element using buckets
!...  Sort elements into bins
      nelem=0
      isuc=1 !success flag
      do i=1,ne
        xm(1)=dmin1(x(nm(i,1)),x(nm(i,2)),x(nm(i,3)))
        xm(2)=dmax1(x(nm(i,1)),x(nm(i,2)),x(nm(i,3)))
        ym(1)=dmin1(y(nm(i,1)),y(nm(i,2)),y(nm(i,3)))
        ym(2)=dmax1(y(nm(i,1)),y(nm(i,2)),y(nm(i,3)))
        if(i34(i)==4) then
          xm(1)=dmin1(xm(1),x(nm(i,4)))
          xm(2)=dmax1(xm(2),x(nm(i,4)))
          ym(1)=dmin1(ym(1),y(nm(i,4)))
          ym(2)=dmax1(ym(2),y(nm(i,4)))
        endif
          
        do l=1,2
          ixm(l)=0
          do ix=1,nx-1
            if(xm(l)>=xi(ix).and.xm(l)<=xi(ix+1)) then
              ixm(l)=ix
              exit
            endif
          enddo !ix
          if(ixm(l)==0) then
            write(*,*)'Cannot find in x-direction:',i
            stop
          endif

          iym(l)=0
          do iy=1,ny-1
            if(ym(l)>=yi(iy).and.ym(l)<=yi(iy+1)) then
              iym(l)=iy
              exit
            endif
          enddo !iy
          if(iym(l)==0) then
            write(*,*)'Cannot find in y-direction:',i
            stop
          endif
        enddo !l=1,2
        if(ixm(1)>ixm(2).or.iym(1)>iym(2)) then
          write(*,*)'Null box:',ixm(1),ixm(2),iym(1),iym(2)
          stop
        endif

        do ix=ixm(1),ixm(2)
          do iy=iym(1),iym(2)
            nelem(ix,iy)=nelem(ix,iy)+1
            if(nelem(ix,iy)>mnelem) isuc=0
            ielem(ix,iy,min0(mnelem,nelem(ix,iy)))=i
          enddo !iy
        enddo !ix
      enddo !i=1,ne

      open(10,file='box_summary.out')
      write(10,*)'# of elements in each bucket'
      write(10,*)'ith bin is between line i and i+1'
      do iy=1,ny-1
        maxn=0
        do ix=1,nx-1
          if(nelem(ix,iy)>maxn) maxn=nelem(ix,iy)
        enddo !ix
        write(10,*)'Row ',iy,' ,max # of elements=',maxn
        write(10,*)(nelem(ix,iy),ix=1,nx-1)
      enddo !iy
      write(10,*)'First element in buck (1,1):',ielem(1,1,1)
      close(10)
      if(isuc==0) then
        print*, 'Too many elements'
        stop
      else
        print*, 'done bucket sorting...'
      endif

!...  Find parent elements
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
      read(19,*) ivcor
      if(ivcor==1) then !traditional sigma
        h_c=0
      else if(ivcor==2) then !s1
        read(19,*)h_c,theta_b,theta_f
        if(h_c<0) then
          write(*,*)'Negative h_c:',h_c
          stop
        endif
        if(theta_b<0.or.theta_b>1) then
          write(*,*)'Wrong theta_b:',theta_b
          stop
        endif
        if(theta_f<=0.or.theta_f>20) then
          write(*,*)'Wrong theta_f:',theta_f 
          stop
        endif
!       Pre-compute constants
        s_con1=dsinh(theta_f)
      else if(ivcor==3) then !s2
        read(19,*)h_c,itheta_1,itheta_2
        if(itheta_1<=0.or.itheta_2<=0) then
          write(*,*)'Wrong itheta_1,2:',itheta_1,itheta_2
          stop
        endif
      else
        write(*,*)'Unknown sigma:',ivcor
        stop
      endif

      read(19,*) nvrt_fg
      if(nvrt_fg>mnv.or.nvrt_fg<3) then
        write(*,*)'nvrt > mnv or nvrt<3'
        stop
      endif

      sigma(1)=-1 !bottom
      sigma(nvrt_fg)=0 !surface
      do k=2,nvrt_fg-1
        read(19,*) j,sigma(k)
        if(sigma(k)<=sigma(k-1).or.sigma(k)>=0) then
          write(*,*)'Check sigma levels at:',k
          stop
        endif
      enddo
      close(19)

!     Compute C(s) and C'(s)
      do k=1,nvrt_fg
        if(ivcor==1) then
          cs(k)=sigma(k); dcs(k)=1
        else if(ivcor==2) then
          cs(k)=(1-theta_b)*dsinh(theta_f*sigma(k))/dsinh(theta_f)+ &
     &theta_b*(dtanh(theta_f*(sigma(k)+0.5))-dtanh(theta_f*0.5))/2/dtanh(theta_f*0.5)
          dcs(k)=(1-theta_b)*theta_f*dcosh(theta_f*sigma(k))/dsinh(theta_f)+ &
     &theta_b*theta_f/2/dtanh(theta_f*0.5)/dcosh(theta_f*(sigma(k)+0.5))**2
        else if(ivcor==3) then
          cs(k)=-(1-(1+sigma(k))**(2*itheta_1))**itheta_2
          if(k==nvrt_fg.and.itheta_2==1) then
            dcs(k)=2*itheta_1
          else
            dcs(k)=2*itheta_1*itheta_2*(1+sigma(k))**(2*itheta_1-1)*(1-(1+sigma(k))**(2*itheta_1))**(itheta_2-1)
          endif
        endif
      enddo !k=1,nvrt_fg

!...  Compute z-cor assuming eta=0
      do i=1,npfg
        do k=1,nvrt_fg
          if(dpfg(i)<=h_c) then !dpfg(i)>0
!            small3=1.e-2 !small number
!            hhat2=h_c+small3 !>0
!            z_til=h_c*sigma(k)+small3*cs(k)
!            sig_til=z_til/hhat2
!            if(sig_til<-1.or.sig_til>0) then
!              write(*,*)'Back-up failed (0):',i,k,sig_til,z_til,hhat2,cs(k),sigma(k)
!              stop
!            endif
            zfg(i,k)=sigma(k)*dpfg(i) !sig_til*dpfg(i)
          else
            zfg(i,k)=h_c*sigma(k)+(dpfg(i)-h_c)*cs(k)
          endif
        enddo !k
      enddo !i

      arco=0
      iparen=0 !flags
      iparen_b=0 !flags
      loop1: do i=1,npfg
        ratmin=1.e25 !min. area ratio for pt i
        ix1=0; ix2=0
        do ix=1,nx-1
          if(xfg(i)>=xi(ix).and.xfg(i)<=xi(ix+1)) then
            if(xfg(i)==xi(ix)) then
              ix1=max0(ix-1,1); ix2=ix
            else if(xfg(i)==xi(ix+1)) then
              ix1=ix; ix2=max0(ix+1,nx-1)
            else
              ix1=ix; ix2=ix
            endif
            exit
          endif
        enddo !ix
        if(ix1==0.or.ix2==0) then
          write(*,*)'Point out of x-bound:',i
          stop
        endif
        iy1=0; iy2=0
        do iy=1,ny-1
          if(yfg(i)>=yi(iy).and.yfg(i)<=yi(iy+1)) then
            if(yfg(i)==yi(iy)) then
              iy1=max0(iy-1,1); iy2=iy
            else if(yfg(i)==yi(iy+1)) then
              iy1=iy; iy2=max0(iy+1,ny-1)
            else
              iy1=iy; iy2=iy
            endif
            exit
          endif
        enddo !iy
        if(iy1==0.or.iy2==0) then
          write(*,*)'Point out of y-bound:',i
          stop
        endif

!       Check if empty boxes
        icount=0
        do ix=ix1,ix2
          do iy=iy1,iy2
            icount=icount+nelem(ix,iy)
          enddo !iy
        enddo !ix
        if(icount==0) then
!          print*, i,' has empty box'
          ix1=1; ix2=1; iy1=1; iy2=1
        endif

        do ix=ix1,ix2
          do iy=iy1,iy2
            if(icount==0) then
              k2=ne
            else
              k2=nelem(ix,iy)
            endif
            do k=1,k2
              if(icount==0) then
                ie=k
              else
                ie=ielem(ix,iy,k)
              endif
              n1=nm(ie,1); n2=nm(ie,2); n3=nm(ie,3)
              ar1=signa(xfg(i),x(n2),x(n3),yfg(i),y(n2),y(n3))  
              ar2=signa(x(n1),xfg(i),x(n3),y(n1),yfg(i),y(n3))
              ar3=signa(x(n1),x(n2),xfg(i),y(n1),y(n2),yfg(i))
              bb=dabs(ar1)+dabs(ar2)+dabs(ar3)
              aa=dabs(signa(x(n1),x(n2),x(n3),y(n1),y(n2),y(n3)))
              rat=dabs(bb-aa)/aa
              if(rat<ratmin) then
                ratmin=rat
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
                cycle loop1
              endif
              if(i34(ie)==4) then
                n4=nm(ie,4)
                ar1=signa(xfg(i),x(n3),x(n4),yfg(i),y(n3),y(n4))  
                ar2=signa(x(n1),xfg(i),x(n4),y(n1),yfg(i),y(n4))
                ar3=signa(x(n1),x(n3),xfg(i),y(n1),y(n3),yfg(i))
                bb=dabs(ar1)+dabs(ar2)+dabs(ar3)
                aa=dabs(signa(x(n1),x(n3),x(n4),y(n1),y(n3),y(n4)))
                rat=dabs(bb-aa)/aa
                if(rat<ratmin) then
                  ratmin=rat
                  iparen_b(i)=ie
                endif
                if(rat<1.e-4) then
                  iparen(i)=ie
!                 Put in the form \sum_{k=1}^{4} {s(k)*arco(k)}
                  arco(i,1)=dmax1(0.d0,dmin1(1.d0,ar1/aa))
                  arco(i,2)=0
                  arco(i,3)=dmax1(0.d0,dmin1(1.d0,ar2/aa))
                  if(arco(i,1)+arco(i,3)>1) then
                    arco(i,4)=0
                    arco(i,3)=1-arco(i,1)
                  else
                    arco(i,4)=1-arco(i,1)-arco(i,3)
                  endif
                  cycle loop1
                endif
              endif !quad
            enddo !k
          enddo !iy
        enddo !ix
      end do loop1 !i
      
!     Back-up parents
      icount=0
      do i=1,npfg
        if(iparen(i)==0) then
          icount=icount+1
          if(iparen_b(i)==0) then
            print*, 'No back-up:',i
            stop
          endif
!          write(88,*)i,iparen_b(i)
          iparen(i)=iparen_b(i)
          if(i34(iparen(i))==3) then
            arco(i,1)=1./3; arco(i,2)=1./3; arco(i,3)=1./3
          else
            arco(i,1)=0.25; arco(i,2)=0.25; arco(i,3)=0.25; arco(i,4)=0.25
          endif
        endif
      enddo !i
!      print*, icount,' pts have no immediate parent elements...:'

!...  Time iteration
!...
      do iday=1,max0(1,ndays)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      write(it_char,'(i12)')iday
      open(63,file=it_char//'_'//file63,status='old',access='direct',recl=nbyte)
      open(65,file=it_char//'_salt.63',status='old',access='direct',recl=nbyte)

      irec=irec0
      do it1=1,nrec
        read(63,rec=irec+1) float1
        time=float1
        read(63,rec=irec+2) it
        irec=irec+2

!        print*, 'time in days = ',time/86400

        do i=1,np
          read(63,rec=irec+i) float1
        enddo !i
        irec=irec+np

        if(i23d==2) then
          do i=1,np
            read(63,rec=irec+1) float1
            eta2(i)=float1
            irec=irec+1
          enddo !i

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
            if(iflag==1) print*, 'Warning: parent element has dry node:',i
          enddo !i

!         Output
          write(16,*)real(time)
          do i=1,npfg
            write(16,*)iglobal(i),real(etafg(i))
          enddo !i

        else !i23d=3 
          do i=1,np
            do k=max0(1,kbp(i)),nvrt
              if(ifile==2) then !ST
                read(63,rec=irec+1) float1
                read(65,rec=irec+1) float2
                irec=irec+1
                tnd(i,k)=float1
                snd(i,k)=float2
              else !uv
                read(63,rec=irec+1) float1
                read(63,rec=irec+2) float2
                irec=irec+2
                tnd(i,k)=float1 !use for u
                snd(i,k)=float2 !for v
              endif
            enddo !k
          enddo !i

!         Do interpolation (no vertical interpolation)
!          print*, 'Day ',iday
          do i=1,npfg
            ie=iparen(i)
            kbe=max0(kbp(nm(ie,1)),kbp(nm(ie,2)),kbp(nm(ie,3)))
            if(kbe==0) then
              write(*,*)'All dry'
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
                  write(*,*)'Cannot find a level:',i,k,zfg(i,k)
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
                  print*, 'Parent element dry:',i
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
          write(17,*)real(time)
          if(ifile==2) write(18,*)real(time)
          do i=1,npfg
            if(ifile==2) then
              write(17,*)iglobal(i),(real(tfg(i,k)),k=1,nvrt_fg)
              write(18,*)iglobal(i),(real(sfg(i,k)),k=1,nvrt_fg)
            else !uv
              write(17,*)iglobal(i),(real(tfg(i,k)),real(sfg(i,k)),k=1,nvrt_fg)
            endif
          enddo !i
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

