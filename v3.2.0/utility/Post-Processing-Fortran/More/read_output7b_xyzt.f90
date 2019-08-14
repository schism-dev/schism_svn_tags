!
!****************************************************************************************
!											*
!	read_output7b with (x,y,z,time) read in from station.xyzt (z>=0 is distance from F.S.; time in sec) 
!         for 3D variables (surface values for 2D variables). Interpolation in time.
!       Inputs: (1) binary files;
!               (2) station.xyzt: make sure all times are after 1st record (to ensure interpolation in time); 
!                                 pad extra days before and after if necessary.
!               (3) screen inputs: file63; invalid value (for out of domain, dry etc); start corie day for output
!       Outputs: fort.1[89]; fort.11 (fatal errors); fort.12: nonfatal errors.
!											*
!       ifort -Bstatic -assume byterecl -O3 -o read_output7b_xyzt read_output7b_xyzt.f90
!****************************************************************************************
!
      program read_out
      parameter(nbyte=4)
!      parameter(mnp=80000)
!      parameter(mnp2=80000) !# of pts in station.xyzt
!      parameter(mne=160000)
!      parameter(mnv=100)
      character*30 file63
      character*12 it_char
      character*48 start_time,version,variable_nm,variable_dim
      character*48 data_format
!      dimension sigma(mnv),cs(mnv),ztot(mnv),x(mnp),y(mnp),dp(mnp),kbp00(mnp),kfp(mnp)
!      dimension nm(mne,4),out(3,mnv,2),out2(2,mnv,2),icum(mnp,mnv,2),eta2(mnp),node3(mnp2,3),arco(mnp2,3)
!      dimension ztmp(mnv),x00(mnp2),y00(mnp2),iep(mnp2),z00(mnp2),t00(mnp2)
!      dimension iday(2,mnp2),irecord(2,mnp2),times(2,mnp2)
      allocatable ::  sigma(:),cs(:),ztot(:),x(:),y(:),dp(:),kbp00(:),kfp(:)
      allocatable ::  nm(:,:),out(:,:,:),out2(:,:,:),icum(:,:,:),eta2(:),node3(:,:),arco(:,:)
      allocatable ::  ztmp(:),x00(:),y00(:),iep(:),z00(:),t00(:)
      allocatable ::  iday(:,:),irecord(:,:),times(:,:)
      dimension swild(3),out3(2,2),out4(2)
      
      print*, 'Input file to read from (without *_):'
      read(*,'(a30)')file63
      
!     Invliad number used for 3D variables: below bottom; dry spot; no parents
      print*, 'Input values to be used for invalid place:'
      read(*,*)rjunk

      print*, 'Input starting corie day:'
      read(*,*)start_corie
      
      open(10,file='station.xyzt',status='old')
      read(10,*) 
      read(10,*) nxy
      allocate(x00(nxy),y00(nxy),z00(nxy),t00(nxy),stat=istat)
      if(istat/=0) stop 'Falied to allocate (1)'

      do i=1,nxy
        read(10,*)j,x00(i),y00(i),z00(i),t00(i) !z00>=0 from F.S.
        if(z00(i)<0) then
          write(*,*)'Invalid z value:',i; stop
        endif
      enddo !i
      close(10)

!...  Header
!...
      open(63,file='1_'//file63,status='old',access='direct',recl=nbyte)
      irec=0
      do m=1,48/nbyte
        read(63,rec=irec+m) data_format(nbyte*(m-1)+1:nbyte*m)
      enddo
      if(data_format.ne.'DataFormat v5.0') then
        print*, 'This code reads only v5.0:  ',data_format
        stop
      endif
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
      read(63,rec=irec+3) nspool
      read(63,rec=irec+4) ivs
      read(63,rec=irec+5) i23d
      irec=irec+5

!      print*, 'i23d=',i23d,' nrec= ',nrec

!     Vertical grid
      read(63,rec=irec+1) nvrt
      read(63,rec=irec+2) kz
      read(63,rec=irec+3) h0
      read(63,rec=irec+4) h_s
      read(63,rec=irec+5) h_c
      read(63,rec=irec+6) theta_b
      read(63,rec=irec+7) theta_f
      irec=irec+7
      allocate(ztot(nvrt),sigma(nvrt),cs(nvrt),ztmp(nvrt),stat=istat)
      if(istat/=0) stop 'Falied to allocate (2)'

      do k=1,kz-1
        read(63,rec=irec+k) ztot(k)
      enddo
      do k=kz,nvrt
        kin=k-kz+1
        read(63,rec=irec+k) sigma(kin)
        cs(kin)=(1-theta_b)*sinh(theta_f*sigma(kin))/sinh(theta_f)+ &
     &theta_b*(tanh(theta_f*(sigma(kin)+0.5))-tanh(theta_f*0.5))/2/tanh(theta_f*0.5)
      enddo
      irec=irec+nvrt

!     Horizontal grid
      read(63,rec=irec+1) np
      read(63,rec=irec+2) ne
      irec=irec+2
      allocate(x(np),y(np),dp(np),kbp00(np),kfp(np), &
     &nm(ne,4),out(3,nvrt,2),out2(2,nvrt,2),icum(np,nvrt,2),eta2(np),node3(nxy,3), &
     &arco(nxy,3),iep(nxy),iday(2,nxy),irecord(2,nxy),times(2,nxy),stat=istat)
      if(istat/=0) stop 'Falied to allocate (3)'

      do m=1,np
        read(63,rec=irec+1)x(m)
        read(63,rec=irec+2)y(m)
        read(63,rec=irec+3)dp(m)
        read(63,rec=irec+4)kbp00(m)
        irec=irec+4
      enddo !m=1,np
      do m=1,ne
        read(63,rec=irec+1)i34
        irec=irec+1
        do mm=1,i34
          read(63,rec=irec+1)nm(m,mm)
          irec=irec+1
        enddo !mm
      enddo !m
      irec0=irec

!      print*, 'last element',(nm(ne,j),j=1,3)

!...  Find parent element for (x00,y00)
      iep=0
      arco=1./3 !initialize for pts without parents
      do l=1,nxy
        node3(l,1:3)=nm(1,1:3) !initialize for pts without parents
      enddo !l

      do i=1,ne
        do l=1,nxy
          aa=0
          ar=0 !area
          do j=1,3
            j1=j+1
            j2=j+2
            if(j1>3) j1=j1-3
            if(j2>3) j2=j2-3
            n0=nm(i,j)
            n1=nm(i,j1)
            n2=nm(i,j2)
            swild(j)=signa(x(n1),x(n2),x00(l),y(n1),y(n2),y00(l)) !temporary storage
            aa=aa+abs(swild(j))
            if(j==1) ar=signa(x(n1),x(n2),x(n0),y(n1),y(n2),y(n0))
          enddo !j
          if(ar<=0) then
            print*, 'Negative area:',ar
            stop
          endif
          ae=abs(aa-ar)/ar
          if(ae<=1.e-5) then
            iep(l)=i
            node3(l,1:3)=nm(i,1:3)
            arco(l,1:3)=swild(1:3)/ar
            arco(l,1)=max(0.,min(1.,arco(l,1)))
            arco(l,2)=max(0.,min(1.,arco(l,2)))
            if(arco(l,1)+arco(l,2)>1) then 
              arco(l,3)=0
              arco(l,2)=1-arco(l,1)
            else
              arco(l,3)=1-arco(l,1)-arco(l,2)
            endif
            cycle
          endif
        enddo !l; build pts

        ifl=0 !flag
        do l=1,nxy
          if(iep(l)==0) then
            ifl=1
            exit
          endif
        enddo !l
        if(ifl==0) exit
      enddo !i=1,ne

      do j=1,nxy
        if(iep(j)==0) then
          write(12,*)'Cannot find a parent for pt:',j,x00(j),y00(j)
!          stop
        endif
      enddo !j

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

!...  Compute stack and record # for each pt
      do i=1,nxy
!       Check if time is before first record
        if(t00(i)<dtout) then
          write(11,*)'Time before first record; try to padd extra day:',i,t00(i)
          stop
        endif

!       Lower and upper bound days and record #s for t00
        iday(1,i)=(t00(i)-dtout)/nrec/dtout+1
        if(iday(1,i)<1) then
          write(11,*)'Impossible'; stop
        else
          irecord(1,i)=(t00(i)-(iday(1,i)-1)*nrec*dtout)/dtout
          times(1,i)=((iday(1,i)-1)*nrec+irecord(1,i))*dtout
          iday(2,i)=t00(i)/nrec/dtout+1
          irecord(2,i)=(t00(i)-(iday(2,i)-1)*nrec*dtout)/dtout+1
          times(2,i)=((iday(2,i)-1)*nrec+irecord(2,i))*dtout
        endif

        if(irecord(1,i)>nrec.or.irecord(2,i)>nrec) then
          write(11,*)'Record # overflow: ',i,irecord(:,i)
          stop
        endif
        if(t00(i)<times(1,i).or.t00(i)>times(2,i)) then
          write(11,*)'Wrong time bounds:',i,t00(i),times(:,i),iday(:,i),irecord(:,i)
          stop
        endif
      enddo !i=1,nxy

!...  Time iteration
!...
      do i=1,nxy
        loop1: do l=1,2 !2 times
          write(it_char,'(i12)')iday(l,i)
          open(63,file=it_char//'_'//file63,status='old',access='direct',recl=nbyte)

          if(i23d==2) then
            irec=irec0+(irecord(l,i)-1)*(2+np+np*ivs)
          else
            irec=irec0+(irecord(l,i)-1)*(2+np+icum(np,nvrt,ivs))
          endif

!          read(63,rec=irec+1) time
!          read(63,rec=irec+2) it
          irec=irec+2

!         Only read in 3 nodes to speed up
          do j=1,3
            nd=node3(i,j)
            read(63,rec=irec+nd) eta2(nd)
          enddo !j
!          do ii=1,np
!            read(63,rec=irec+ii) eta2(ii)
!          enddo !ii
          irec=irec+np

          out2(l,:,:)=0
          out3(l,:)=0
          if(i23d==2) then
            do j=1,3 !nodes
              nd=node3(i,j)
              do m=1,ivs
                read(63,rec=irec+(nd-1)*ivs+m) tmp
                out2(l,1,m)=out2(l,1,m)+arco(i,j)*tmp
              enddo !m
            enddo !j
!           irec=irec+np*ivs
          else !i23d=3 
            do j=1,3 !nodes
              nd=node3(i,j)
              do k=max0(1,kbp00(nd)),nvrt
                do m=1,ivs
                  read(63,rec=irec+icum(nd,k,m)) out(j,k,m)
                enddo !m
              enddo !k
            enddo !j
!            irec=irec+icum(np,nvrt,ivs)

!           Do interpolation
            etal=0; dep=0; idry=0
            do j=1,3
              nd=node3(i,j)
              if(eta2(nd)+dp(nd)<h0) idry=1
              etal=etal+arco(i,j)*eta2(nd)
              dep=dep+arco(i,j)*dp(nd)
!             Debug
!              write(11,*)i,j,nd,dp(nd),arco(i,j)
            enddo !j
            if(idry==1) then
              out3(:,:)=rjunk
              exit loop1
            else !element wet
!             Compute z-coordinates
              do k=kz,nvrt
                kin=k-kz+1
                hmod2=min(dep,h_s)
                if(hmod2<=h_c) then
                  ztmp(k)=sigma(kin)*(hmod2+etal)+etal
                else if(etal<=-h_c-(hmod2-h_c)*theta_f/sinh(theta_f)) then
                  write(11,*)'Pls choose a larger h_c (2):',etal,h_c
                  stop
                else
                  ztmp(k)=etal*(1+sigma(kin))+h_c*sigma(kin)+(hmod2-h_c)*cs(kin)
                endif

!               Following to prevent underflow
                if(k==kz) ztmp(k)=-hmod2
                if(k==nvrt) ztmp(k)=etal
              enddo !k

              if(dep<=h_s) then
                kbpl=kz
              else !z levels
!               Find bottom index
                kbpl=0
                do k=1,kz-1
                  if(-dep>=ztot(k).and.-dep<ztot(k+1)) then
                    kbpl=k
                    exit
                  endif
                enddo !k
                if(kbpl==0) then
                  write(11,*)'Cannot find a bottom level:',dep,i
                  stop
                endif
                ztmp(kbpl)=-dep
                do k=kbpl+1,kz-1
                  ztmp(k)=ztot(k)
                enddo !k
              endif

              do k=kbpl+1,nvrt
                if(ztmp(k)-ztmp(k-1)<=0) then
                  write(11,*)'Inverted z-level:',etal,dep,ztmp(k)-ztmp(k-1)
                  stop
                endif
              enddo !k
       
              do k=kbpl,nvrt
                do m=1,ivs
                  do j=1,3
                    nd=node3(i,j)
                    kin=max(k,kbp00(nd))
                    out2(l,k,m)=out2(l,k,m)+arco(i,j)*out(j,kin,m)
                  enddo !j
                enddo !m
              enddo !k

!             Interplate in vertical
              k0=0
              do k=kbpl,nvrt-1
                if(ztmp(nvrt)-z00(i)>=ztmp(k).and.ztmp(nvrt)-z00(i)<=ztmp(k+1)) then
                  k0=k
                  rat=(ztmp(nvrt)-z00(i)-ztmp(k))/(ztmp(k+1)-ztmp(k))
                  exit
                endif
              enddo !k
              if(k0==0) then
                out3(:,:)=rjunk
                exit loop1
!               write(12,*)'Warning: failed to find a vertical level:',it,i
              else
                do m=1,ivs
                  out3(l,m)=out2(l,k0,m)*(1-rat)+out2(l,k0+1,m)*rat
                enddo !m
              endif
            endif !dry/wet
          endif !i23d
        enddo loop1 !l=1,2; 2 times

!       Interpolate in time
        trat=(t00(i)-times(1,i))/(times(2,i)-times(1,i)) !must be [0,1]
        if(i23d==2) then
          if(iep(i)==0) then !no parents
            out4(1:ivs)=rjunk
          else
            out4(1:ivs)=out2(1,1,1:ivs)*(1-trat)+out2(2,1,1:ivs)*trat
          endif
          write(18,'(e16.8,8(1x,f15.3))')t00(i)/86400+start_corie,out4(1:ivs),x00(i),y00(i)
        else !3D
          if(iep(i)==0) then !no parents
            out4(1:ivs)=rjunk
          else
            out4(1:ivs)=out3(1,1:ivs)*(1-trat)+out3(2,1:ivs)*trat
          endif
          write(18,'(e16.8,4(1x,f15.3))')t00(i)/86400+start_corie,out4(1),x00(i),y00(i),z00(i)
          if(ivs==2) write(19,'(e16.8,1x,f15.3)')t00(i)/86400+start_corie,out4(2)
        endif
      enddo !i=1,nxy

      stop
      end

      function signa(x1,x2,x3,y1,y2,y3)
!...  Compute signed area formed by pts 1,2,3

      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2

      return
      end

