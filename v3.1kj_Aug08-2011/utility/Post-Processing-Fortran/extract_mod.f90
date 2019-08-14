!     Basic routines for reading SELFE binary outputs
!     Authors: Charles Seaton, Joseph Zhang, Paul Turner
!     Date: April 2011

!     Global data pool
      module extract_mod
        parameter(nbyte=4)
        character(len=12), save :: it_char
        character(len=48), save :: start_time,version,variable_nm,variable_dim,data_format
        integer, save :: itransect,ics,ifs,nxy,irec0,ivs,nrec,i23d,nspool,nvrt,kz,np,ne, &
     &ixy_or_xyz,ialong_S,klev0,zout,it_tot,debug
        real, save :: theta_f,theta_b,h_c,h_s,h0,dtout,time,fill_in

        real, save, allocatable :: sigma(:),cs(:),ztot(:),x(:),y(:),dp(:),arco(:,:), &
                                   ztmp(:),x00(:),y00(:),z00(:),t00(:),rl2min(:), &
                                   varout(:,:,:,:),outtime(:),eta2(:),varout2(:,:,:),varout3(:,:,:)
        integer, save, allocatable :: kbp00(:),kfp(:),nm(:,:),icum(:,:,:),node3(:,:),iep(:)
      end module extract_mod

!================================================================
!================================================================
      function signa(x1,x2,x3,y1,y2,y3)
!...  Compute signed area formed by pts 1,2,3

      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2

      end function signa

!================================================================
!================================================================
      subroutine readheader(fname)
      use extract_mod
      
      character(len=500),intent(in) :: fname

      open(63,file=fname,status='old',access='direct',recl=nbyte)
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

      read(63,rec=irec+1) nrec
      read(63,rec=irec+2) dtout
      read(63,rec=irec+3) nspool
      read(63,rec=irec+4) ivs
      read(63,rec=irec+5) i23d
      irec=irec+5

!      write(*,'(a48)')data_format
!      write(*,'(a48)')version
!      write(*,'(a48)')start_time
!      write(*,'(a48)')variable_nm
!      write(*,'(a48)')variable_dim
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

!      print*, 'nvrt= ',nvrt
!      print*, 'kz= ',kz
!      print*, 'h0= ',h0
!      print*, 'h_s= ',h_s
!      print*, 'h_c= ',h_c
!      print*, 'theta_b= ',theta_b
!      print*, 'theta_f= ',theta_f

      if(allocated(ztot)) deallocate(ztot)
      if(allocated(sigma)) deallocate(sigma)
      if(allocated(cs)) deallocate(cs)
      if(allocated(ztmp)) deallocate(ztmp)
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

      if(allocated(x)) deallocate(x)
      if(allocated(y)) deallocate(y)
      if(allocated(dp)) deallocate(dp)
      if(allocated(kbp00)) deallocate(kbp00)
      if(allocated(kfp)) deallocate(kfp)
      if(allocated(eta2)) deallocate(eta2)
      if(allocated(nm)) deallocate(nm)
      if(allocated(icum)) deallocate(icum)
      allocate(x(np),y(np),dp(np),kbp00(np),kfp(np),nm(ne,4), &
     &icum(np,nvrt,2),eta2(np),stat=istat)
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

!      print*, 'last node',x(np),y(np),dp(np),kbp00(np)
!      print*, 'last element',(nm(ne,j),j=1,3)

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

      end subroutine readheader

!================================================================
!================================================================
!     Find parent elements for points
      subroutine find_parents
      use extract_mod

      real :: swild(3)

      if(allocated(rl2min)) deallocate(rl2min)
      if(allocated(node3)) deallocate(node3)
      if(allocated(arco)) deallocate(arco)
      if(allocated(iep)) deallocate(iep)
      allocate(rl2min(nxy),node3(nxy,3),arco(nxy,3),iep(nxy),stat=istat)
      if(istat/=0) stop 'Falied to allocate (10)'

!...  Find parent element for (x00,y00)
      iep=0
      if(ics==1) then !Cartesian
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
              if(j==1) then
                ar=signa(x(n1),x(n2),x(n0),y(n1),y(n2),y(n0))
                if(ar<=0) then
                  print*, 'Negative area:',ar,i,n1,n2,n0,x(n1),y(n1),x(n2),y(n2),x(n0),y(n0)
                  stop
                endif
              endif
            enddo !j
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
      else !lat/lon
        rl2min=1.e25 !min distance^2
        do i=1,np
          do l=1,nxy
            rl2=(x(i)-x00(l))**2+(y(i)-y00(l))**2
            if(rl2<rl2min(l)) then
              rl2min(l)=rl2
              iep(l)=1 !actual elem. # not used
              node3(l,1:3)=i
              arco(l,1:3)=1./3
            endif
          enddo !l=1,nxy
        enddo !i=1,np
      endif !ics

      do j=1,nxy
        if(iep(j)==0) then
          print*, 'Cannot find a parent for pt:',j,x00(j),y00(j)
          stop
        endif
      enddo !j

      end subroutine find_parents

!================================================================
!================================================================
!     Read binary files at a given time step
      subroutine readdata(fname)
      use extract_mod

      character(len=*), intent(in) :: fname

      real, allocatable :: out(:,:,:,:),out2(:,:,:)
 
      allocate(out(nxy,3,nvrt,2),out2(nxy,nvrt,2),stat=istat)
      if(istat/=0) stop 'Falied to allocate (7)'

!     Check transect
      if(itransect==1.and.i23D==2) then
        print*, '2D variables do not have transect'
        stop
      endif

!     Init. for abnormal cases
      varout=fill_in

      open(63,file=fname,status='old',access='direct',recl=nbyte)

      irec=irec0
      do it1=1,nrec
        read(63,rec=irec+1) time
        read(63,rec=irec+2) it
        irec=irec+2
!        it_tot=it_tot+1
!        time=it_tot*dtout
        outtime(it1)=time

!        print*, 'time=',time/86400

        do i=1,nxy
          do j=1,3
            nd=node3(i,j)
            read(63,rec=irec+nd) eta2(nd)
          enddo !j
        enddo !i
        irec=irec+np

        out2=0
        if(i23d==2) then !2D variable
          do i=1,nxy
            do j=1,3 !nodes
              nd=node3(i,j)
              do m=1,ivs
                read(63,rec=irec+(nd-1)*ivs+m) tmp
                out2(i,1,m)=out2(i,1,m)+arco(i,j)*tmp
              enddo !m
            enddo !j
          enddo !i
          irec=irec+np*ivs
          do i=1,nxy
            varout(1:ivs,1,i,it1)=out2(i,1,1:ivs)
          enddo !i
        else !i23d=3 
          do i=1,nxy
            do j=1,3 !nodes
              nd=node3(i,j)
              do k=max0(1,kbp00(nd)),nvrt
                do m=1,ivs
                  read(63,rec=irec+icum(nd,k,m)) out(i,j,k,m)
                enddo !m
              enddo !k
            enddo !j
          enddo !i
          irec=irec+icum(np,nvrt,ivs)

!         Do interpolation
          do i=1,nxy
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
            else !element wet
!             Compute z-coordinates
              do k=kz,nvrt
                kin=k-kz+1
                hmod2=min(dep,h_s)
                if(hmod2<=h_c) then
                  ztmp(k)=sigma(kin)*(hmod2+etal)+etal
                else if(etal<=-h_c-(hmod2-h_c)*theta_f/sinh(theta_f)) then
                  write(*,*)'Pls choose a larger h_c (2):',etal,h_c
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
                  write(*,*)'Cannot find a bottom level:',dep,i
                  stop
                endif
                ztmp(kbpl)=-dep
                do k=kbpl+1,kz-1
                  ztmp(k)=ztot(k)
                enddo !k
              endif

              do k=kbpl+1,nvrt
                if(ztmp(k)-ztmp(k-1)<=0) then
                  write(*,*)'Inverted z-level:',etal,dep,ztmp(k)-ztmp(k-1)
                  stop
                endif
              enddo !k
       
              do k=kbpl,nvrt
                do m=1,ivs
                  do j=1,3
                    nd=node3(i,j)
                    kin=max(k,kbp00(nd))
                    out2(i,k,m)=out2(i,k,m)+arco(i,j)*out(i,j,kin,m)
                  enddo !j
                enddo !m
!                write(65,*)i,k,ztmp(k),(out2(i,k,m),m=1,ivs)     
              enddo !k

              if(itransect==1) then !output transect
                do k=kbpl,nvrt
                  do m=1,ivs
                    varout(m,k,i,it1)=out2(i,k,m)
                  enddo !m
                  varout(3,k,i,it1)=ztmp(k) !save depths
                enddo !k
              else !not transect
!               Interplate in vertical
                k0=0
                do k=kbpl,nvrt-1
                  if(ifs==0) then !z coord
                    ztmp2=z00(i)
                  else !relative to F.S.
                    ztmp2=ztmp(nvrt)-z00(i)
                  endif !ifs

                  if(ztmp2>=ztmp(k).and.ztmp2<=ztmp(k+1)) then
                    k0=k
                    rat=(ztmp2-ztmp(k))/(ztmp(k+1)-ztmp(k))
                    exit
                  endif
                enddo !k
                if(k0==0) then
!                  write(12,*)'Warning: failed to find a vertical level:',it,i
                else
                  do m=1,ivs
                    !out3(i,m)=out2(i,k0,m)*(1-rat)+out2(i,k0+1,m)*rat
                    varout(m,1,i,it1)=out2(i,k0,m)*(1-rat)+out2(i,k0+1,m)*rat
                  enddo !m
                endif
              endif !itransect
            endif !dry/wet
          enddo !i=1,nxy
!          write(18,'(e16.8,6000(1x,f12.3))')time,(out3(i,1),i=1,nxy)
!          if(ivs==2) write(19,'(e16.8,6000(1x,f12.3))')time,(out3(i,2),i=1,nxy)

        endif !i23d
      enddo !it1=1,nrec

      deallocate(out,out2)

      end subroutine readdata

!================================================================
!================================================================
!     Read binary files and do interpolate at (x,y,z,t)
      subroutine find_xyzt(basedir,binary_file)
      use extract_mod

      character(len=*), intent(in) :: basedir,binary_file

      integer, allocatable :: iday(:,:),irecord(:,:)
      real, allocatable :: times(:,:),out(:,:,:),out2(:,:,:),out3(:,:,:)

      allocate(iday(2,nxy),irecord(2,nxy),times(2,nxy),out(3,nvrt,2), &
     &out2(2,nvrt,2),out3(2,2,nvrt),stat=istat)
      if(istat/=0) stop 'Falied to allocate (6)'

!     Init. for abnormal cases
      out3=fill_in

!...  Compute stack and record # for each pt
      do i=1,nxy
!       Check if time is before first record
        if(t00(i)<dtout) then
          write(*,*)'Time before first record; try to padd extra day:',i,t00(i)
          stop
        endif

!       Lower and upper bound days and record #s for t00
        iday(1,i)=(t00(i)-dtout)/nrec/dtout+1
        if(iday(1,i)<1) then
          write(*,*)'Impossible'; stop
        else
          irecord(1,i)=(t00(i)-(iday(1,i)-1)*nrec*dtout)/dtout
          times(1,i)=((iday(1,i)-1)*nrec+irecord(1,i))*dtout
          iday(2,i)=t00(i)/nrec/dtout+1
          irecord(2,i)=(t00(i)-(iday(2,i)-1)*nrec*dtout)/dtout+1
          times(2,i)=((iday(2,i)-1)*nrec+irecord(2,i))*dtout
        endif

        if(irecord(1,i)>nrec.or.irecord(2,i)>nrec) then
          write(*,*)'Record # overflow: ',i,irecord(:,i)
          stop
        endif
        if(t00(i)<times(1,i).or.t00(i)>times(2,i)) then
          write(*,*)'Wrong time bounds:',i,t00(i),times(:,i),iday(:,i),irecord(:,i)
          stop
        endif
      enddo !i=1,nxy

      do i=1,nxy
        loop1: do l=1,2 !2 times
          write(it_char,'(i12)')iday(l,i)
          it_char=adjustl(it_char)
          leng=len_trim(it_char)
          open(63,file=trim(basedir)//it_char(1:leng)//'_'//binary_file, &
     &status='old',access='direct',recl=nbyte)

          if(i23d==2) then
            irec=irec0+(irecord(l,i)-1)*(2+np+np*ivs)
          else
            irec=irec0+(irecord(l,i)-1)*(2+np+icum(np,nvrt,ivs))
          endif
          irec=irec+2 !skip 2 records

!         Only read in 3 nodes to speed up
          do j=1,3
            nd=node3(i,j)
            read(63,rec=irec+nd) eta2(nd)
          enddo !j
          irec=irec+np

          out2(l,:,:)=0
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
              !out3=rjunk
              exit loop1
            else !element wet
!             Compute z-coordinates
              do k=kz,nvrt
                kin=k-kz+1
                hmod2=min(dep,h_s)
                if(hmod2<=h_c) then
                  ztmp(k)=sigma(kin)*(hmod2+etal)+etal
                else if(etal<=-h_c-(hmod2-h_c)*theta_f/sinh(theta_f)) then
                  write(*,*)'Pls choose a larger h_c (2):',etal,h_c
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
                  write(*,*)'Inverted z-level:',etal,dep,ztmp(k)-ztmp(k-1)
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

              if(ixy_or_xyz==1) then !xy format
                do m=1,ivs
                  do k=kbpl,nvrt
                    out3(l,m,k)=out2(l,k,m)
                    if(m==1) varout2(3,k,i)=ztmp(k) !save depths
                  enddo !k
                enddo !m
              else !xyz format
!               Interplate in vertical
                k0=0
                do k=kbpl,nvrt-1
                  if(ifs==0) then !relative MSL
                    ztmp2=z00(i)
                  else !relative to F.S.
                    ztmp2=ztmp(nvrt)-z00(i)
                  endif !ifs

                  if(ztmp2>=ztmp(k).and.ztmp2<=ztmp(k+1)) then
                    k0=k
                    rat=(ztmp2-ztmp(k))/(ztmp(k+1)-ztmp(k))
                    exit
                  endif
                enddo !k
                if(k0==0) then
                  !out3=rjunk
                  exit loop1
!                 write(12,*)'Warning: failed to find a vertical level:',it,i
                else
                  do m=1,ivs
                    out3(l,m,1)=out2(l,k0,m)*(1-rat)+out2(l,k0+1,m)*rat
                  enddo !m
                endif
              endif !ixy_or_xyz
            endif !dry/wet
          endif !i23d
        enddo loop1 !l=1,2; 2 times

!       Interpolate in time
        trat=(t00(i)-times(1,i))/(times(2,i)-times(1,i)) !must be [0,1]
        if(i23d==2) then
          varout2(1:ivs,1,i)=out2(1,1,1:ivs)*(1-trat)+out2(2,1,1:ivs)*trat
          !write(18,'(e16.8,8(1x,f15.3))')t00(i)/86400+start_corie,out4(1:ivs),x00(i),y00(i)
        else !3D
          if(ixy_or_xyz==1) then !xy format
            varout2(1:ivs,:,i)=out3(1,1:ivs,:)*(1-trat)+out3(2,1:ivs,:)*trat
          else !xyz format
            varout2(1:ivs,1,i)=out3(1,1:ivs,1)*(1-trat)+out3(2,1:ivs,1)*trat
          endif
          !write(18,'(e16.8,4(1x,f15.3))')t00(i)/86400+start_corie,out4(1),x00(i),y00(i),z00(i)
          !if(ivs==2) write(19,'(e16.8,1x,f15.3)')t00(i)/86400+start_corie,out4(2)
        endif
      enddo !i=1,nxy

      close(63)

      deallocate(iday,irecord,times,out,out2,out3)

      end subroutine find_xyzt

!================================================================
!================================================================
!     Read binary files at a given time step at all nodes
      subroutine readslab(fname)
      use extract_mod

      character(len=*), intent(in) :: fname

      real, allocatable :: out2(:,:,:)
 
      allocate(out2(ivs,nvrt,np),stat=istat)
      if(istat/=0) stop 'Falied to allocate (11)'

!     Init. for abnormal cases
      varout3=fill_in

      open(63,file=fname,status='old',access='direct',recl=nbyte)

      irec=irec0
      do it1=1,nrec
        read(63,rec=irec+1) time
        read(63,rec=irec+2) it
        irec=irec+2
        outtime(it1)=time

        do i=1,np
          read(63,rec=irec+i) eta2(i)
        enddo !i
        irec=irec+np

        if(i23d==2) then !2D variable
          do i=1,np
            do m=1,ivs
              read(63,rec=irec+1) varout3(m,i,it1)
              irec=irec+1
            enddo !m
          enddo !i
        else !i23d=3 
          do i=1,np
            do k=max0(1,kbp00(i)),nvrt
              do m=1,ivs
                read(63,rec=irec+1) out2(m,k,i)
                irec=irec+1
              enddo !m
            enddo !k
          enddo !i

          if(ialong_S==1) then !along a S level
            varout3(:,:,it1)=out2(:,klev0,:)
          else 
!           Not along a S level; do vertical interpolation
            do i=1,np
              if(eta2(i)+dp(i)<h0) then
                idry=1
              else
                idry=0
              endif
              etal=eta2(i)
              dep=dp(i)

              if(idry==1) then
!                print*, 'Dry:',i,eta2(i),dp(i),h0
              else !wet
!               Compute z-coordinates
                do k=kz,nvrt
                  kin=k-kz+1
                  hmod2=min(dep,h_s)
                  if(hmod2<=h_c) then
                    ztmp(k)=sigma(kin)*(hmod2+etal)+etal
                  else if(etal<=-h_c-(hmod2-h_c)*theta_f/sinh(theta_f)) then
                    write(*,*)'Pls choose a larger h_c (4):',etal,h_c
                    stop
                  else
                    ztmp(k)=etal*(1+sigma(kin))+h_c*sigma(kin)+(hmod2-h_c)*cs(kin)
                  endif

!                 Following to prevent underflow
                  if(k==kz) ztmp(k)=-hmod2
                  if(k==nvrt) ztmp(k)=etal
                enddo !k

                if(dep<=h_s) then
                  kbpl=kz
                else !z levels
!                 Find bottom index
                  kbpl=0
                  do k=1,kz-1
                    if(-dep>=ztot(k).and.-dep<ztot(k+1)) then
                      kbpl=k
                      exit
                    endif
                  enddo !k
                  if(kbpl==0) then
                    write(*,*)'Cannot find a bottom level (1):',dep,i
                    stop
                  endif
                  ztmp(kbpl)=-dep
                  do k=kbpl+1,kz-1
                    ztmp(k)=ztot(k)
                  enddo !k
                endif

                do k=kbpl+1,nvrt
                  if(ztmp(k)-ztmp(k-1)<=0) then
                    write(*,*)'Inverted z-level (1):',etal,dep,ztmp(k)-ztmp(k-1)
                    stop
                  endif
                enddo !k
       
!               Interplate in vertical
                k0=0
                do k=kbpl,nvrt-1
                  if(ifs==0) then !z coord
                    ztmp2=zout
                  else !relative to F.S.
                    ztmp2=ztmp(nvrt)-zout
                  endif !ifs

                  if(ztmp2>=ztmp(k).and.ztmp2<=ztmp(k+1)) then
                    k0=k
                    rat=(ztmp2-ztmp(k))/(ztmp(k+1)-ztmp(k))
                    exit
                  endif
                enddo !k
                if(k0==0) then
!                  write(12,*)'Warning: failed to find a vertical level:',it,i
                else
                  do m=1,ivs
                    varout3(m,i,it1)=out2(m,k0,i)*(1-rat)+out2(m,k0+1,i)*rat
                  enddo !m
                endif
              endif !dry/wet
            enddo !i=1,np
          endif !ialong_S
        endif !i23d
      enddo !it1=1,nrec

      deallocate(out2)

      end subroutine readslab

!================================================================
!================================================================
      subroutine read_station(fname)
!     Read in station build point file
      use extract_mod

      character(len=*), intent(in) :: fname

      open(10,file=fname,status='old')
!     itransect=0: time series; 1: transect
!     ics - 1 (Cartesian coordinates) or 2 (lat/lon) used in the binary outputs
!     ifs - 0 (input z are z coordinates) or 1 (input z are relative to free surface;
!           and in this case, make sure all z>=0). This flag has no effects on
!           2D variables or transect
      read(10,*)itransect,ics,ifs
      if(itransect/=0.and.itransect/=1.or.ics/=1.and.ics/=2.or. &
     &ifs/=0.and.ifs/=1) stop 'Wrong flags in bp_file'
      read(10,*) nxy

      if(allocated(x00)) deallocate(x00)
      if(allocated(y00)) deallocate(y00)
      if(allocated(z00)) deallocate(z00)
      allocate(x00(nxy),y00(nxy),z00(nxy),stat=istat)
      if(istat/=0) stop 'Falied to allocate (1)'

      do i=1,nxy
        read(10,*)j,x00(i),y00(i),z00(i)
        if(ifs==1.and.z00(i)<0) then
          write(*,*)'Invalid z value:',i; stop
        endif
      enddo !i
      close(10)

      end subroutine read_station

!================================================================
!================================================================
      subroutine read_xyzt(fname)
!     Read in station point file
      use extract_mod
      character(len=*), intent(in) :: fname

      open(10,file=fname,status='old')
!     ixy_or_xyz - 1: xyt format; 2: xyzt format
!     ics - 1 (Cartesian coordinates) or 2 (lat/lon) used in the binary outputs
!     ifs - 0 (input z are z coordinates) or 1 (input z are relative to free surface;
!           and in this case, make sure all z>=0). This flag has no effects on
!           2D variables or transect
      read(10,*)ixy_or_xyz,ics,ifs
      if(ixy_or_xyz/=1.and.ixy_or_xyz/=2.or.ics/=1.and.ics/=2.or. &
     &ifs/=0.and.ifs/=1) then
        print*, 'Wrong header for input ',xyzt_file,ixy_or_xyz,ics,ifs
        stop
      endif

      read(10,*) nxy

      if(allocated(x00)) deallocate(x00)
      if(allocated(y00)) deallocate(y00)
      if(allocated(z00)) deallocate(z00)
      if(allocated(t00)) deallocate(t00)
      allocate(x00(nxy),y00(nxy),z00(nxy),t00(nxy),stat=istat)
      if(istat/=0) stop 'Falied to allocate (5)'

      do i=1,nxy
        if(ixy_or_xyz==1) then !xyt format
          read(10,*)j,x00(i),y00(i),t00(i)
        else !xyzt format
          read(10,*)j,x00(i),y00(i),z00(i),t00(i) !z00>=0 from F.S.
          if(ifs==1.and.z00(i)<0) then
            write(*,*)'Invalid z value:',i; stop
          endif
        endif !ixy_or_xyz
      enddo !i
      close(10)

      end subroutine read_xyzt

!================================================================
!================================================================
