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
!										*
!	Read binary format v5.0 (SZ models) and compute average field 
!       at a particular Z-level (use above surface/below bottom to get surface/bottom).
!       Skip dry times for 3D variables.
!       Works for mixed quad/tri.
!       Input: *.6?; vgrid.in; iday1,iday2, irec_start,irec_end; z-coor.
!       Output: average.out (gredit for scalar or xmgr5 format for vector)
!										*
!       ifort -O2 -mcmodel=medium -assume byterecl -o compute_average2 compute_average2.f90 ~/SCHISM/svn/trunk/src/Utility/UtilLib/compute_zcor.f90
!       pgf90 -O2 -mcmodel=medium -Mbounds -o compute_average2 compute_average2.f90 ~/SCHISM/svn/trunk/src/Utility/UtilLib/compute_zcor.f90
!       gfortran -ffree-line-length-none -O2 -o compute_average2 compute_average2.f90 ~/SCHISM/svn/trunk/src/Utility/UtilLib/compute_zcor.f90
!********************************************************************************
!
      program read_out
      parameter(nbyte=4)
      character*30 file63
      character*12 it_char
      character*48 start_time,version,variable_nm,variable_dim
      character*48 data_format
  
      integer,allocatable :: i34(:),elnode(:,:)
      allocatable :: sigma(:),cs(:),ztot(:),x(:),y(:),dp(:),kbp00(:),kfp(:)
      allocatable :: out(:,:,:),icum(:,:,:),eta2(:),ztmp(:,:)
      allocatable :: xcj(:),ycj(:),dps(:),kbs(:),xctr(:),yctr(:),dpe(:),kbe(:)
      allocatable :: idry(:),outs(:,:,:),residual(:,:),icounter(:)
      allocatable :: ic3(:,:),isdel(:,:),isidenode(:,:),kbp(:),sigma_lcl(:,:)
      !long int for large files
      integer(kind=8) :: irec
      
      print*, 'Input file to read from (without *_; e.g. salt.63):'
      read(*,'(a30)')file63
      
      print*, 'Input start and end stack #s to read:'
      read(*,*) iday1,iday2

      print*, 'Input start and end record #s in start|end stack respectively:'
!'
      read(*,*) irec_start,irec_end

      print*, 'Input z-coord. (<=0 below MSL):'
      read(*,*) z00

      file63=adjustl(file63)
      len_file63=len_trim(file63)

!...  Header
!...
      write(it_char,'(i12)')iday1
      open(63,file=adjustl(it_char//'_'//file63(1:len_file63)),status='old',access='direct',recl=nbyte)
!'
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

      write(*,'(a48)')data_format
      write(*,'(a48)')version
      write(*,'(a48)')start_time
      write(*,'(a48)')variable_nm
      write(*,'(a48)')variable_dim

      read(63,rec=irec+1) nrec
      read(63,rec=irec+2) dtout
      read(63,rec=irec+3) nspool
      read(63,rec=irec+4) ivs
      read(63,rec=irec+5) i23d
      irec=irec+5

      print*, 'i23d=',i23d,' nrec= ',nrec

!     Vertical grid (obsolete)
      read(63,rec=irec+1) nvrt
!      read(63,rec=irec+2) kz
      read(63,rec=irec+3) h0
!      read(63,rec=irec+4) h_s
!      read(63,rec=irec+5) h_c
!      read(63,rec=irec+6) theta_b
!      read(63,rec=irec+7) theta_f
!      irec=irec+7
!      allocate(ztot(nvrt),sigma(nvrt),cs(nvrt),stat=istat)
!      if(istat/=0) stop 'Falied to allocate (1)'
!
!      do k=1,kz-1
!        read(63,rec=irec+k) ztot(k)
!      enddo
!      do k=kz,nvrt
!        kin=k-kz+1
!        read(63,rec=irec+k) sigma(kin)
!        cs(kin)=(1-theta_b)*sinh(theta_f*sigma(kin))/sinh(theta_f)+ &
!     &theta_b*(tanh(theta_f*(sigma(kin)+0.5))-tanh(theta_f*0.5))/2/tanh(theta_f*0.5)
!      enddo

!     Horizontal grid
      irec=irec+7+nvrt
      read(63,rec=irec+1) np !could be ns,ne also
      read(63,rec=irec+2) ne
      irec=irec+2
      allocate(x(np),y(np),dp(np),kbp00(np),kfp(np),i34(ne),elnode(4,ne),out(np,nvrt,2),idry(np), &
     &eta2(np),xctr(ne),yctr(ne),dpe(ne),kbe(ne),ztmp(nvrt,np),residual(np,2),icounter(np),stat=istat)
      if(istat/=0) stop 'Falied to allocate (2)'

      do m=1,np
        read(63,rec=irec+1)x(m)
        read(63,rec=irec+2)y(m)
        read(63,rec=irec+3)dp(m)
        read(63,rec=irec+4)kbp00(m)
        irec=irec+4
      enddo !m=1,np

      do m=1,ne
        read(63,rec=irec+1)i34(m)
        irec=irec+1
        do mm=1,i34(m)
          read(63,rec=irec+1)elnode(mm,m)
          irec=irec+1
        enddo !mm
      enddo !m
      irec0=irec

!     Read vgrid.in
      allocate(ztot(nvrt),sigma(nvrt),sigma_lcl(nvrt,np),kbp(np))
      call get_vgrid('vgrid.in',np,nvrt,ivcor,kz,h_s,h_c,theta_b,theta_f,ztot,sigma,sigma_lcl,kbp)

      print*, 'last element',elnode(1:i34(ne),ne)

!...  Time iteration
!...
      out=-99 !init.
      ztmp=-99
      residual=0
!      ndays=iday2-iday1+1
      icounter=0 !counter for each node (wet/dry)
      do iday=iday1,iday2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      write(it_char,'(i12)')iday
      open(63,file=adjustl(it_char//'_'//file63(1:len_file63)),status='old',access='direct',recl=nbyte)
!'

      irec=irec0
      do it1=1,nrec
        read(63,rec=irec+1) time
        read(63,rec=irec+2) it
        irec=irec+2

        !print*, 'time in days =',time/86400

        do i=1,np
          read(63,rec=irec+i) eta2(i)
        enddo !i
        irec=irec+np

        if(i23d==2) then !2D
          do i=1,np
            do m=1,ivs
              read(63,rec=irec+1) out(i,1,m)
              irec=irec+1
            enddo !m
          enddo !i

          if(.not.(it1==iday1.and.it1<irec_start.or.it1==iday2.and.it1>irec_end)) then
            do i=1,np
              icounter(i)=icounter(i)+1
              residual(i,1:ivs)=residual(i,1:ivs)+out(i,1,1:ivs) !/ndays/nrec
            enddo !i
          endif
        else !if(i23d==3) then !3D 
          do i=1,np
            do k=max0(1,kbp00(i)),nvrt
              do m=1,ivs
                read(63,rec=irec+1) out(i,k,m)
                irec=irec+1
              enddo !m
            enddo !k
          enddo !i

          !Compute z coordinates
          do i=1,np
            if(ivcor==1) then !localized
              if(dp(i)+eta2(i)<=h0) then
                idry(i)=1
              else !wet
                idry(i)=0
                do k=kbp(i),nvrt
                  ztmp(k,i)=(eta2(i)+dp(i))*sigma_lcl(k,i)+eta2(i)
                enddo !k
              endif !dp
            else if(ivcor==2) then !SZ
              call zcor_SZ(dp(i),eta2(i),h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot, &
     &sigma,ztmp(:,i),idry(i),kbpl)
              kbp(i)=kbpl
            endif !ivcor

            if(idry(i)==0) then !wet
              !Interplate in vertical
              if(z00>=ztmp(nvrt,i)) then !above F.S.
                k0=nvrt-1; rat=1
              else if(z00<=ztmp(kbp(i),i)) then !below bottom; extrapolate
                k0=kbp(i); rat=0
              else !above bottom; cannot be above F.S.
                k0=0
                do k=kbp(i),nvrt-1
                  if(z00>=ztmp(k,i).and.z00<=ztmp(k+1,i)) then
                    k0=k
                    rat=(z00-ztmp(k,i))/(ztmp(k+1,i)-ztmp(k,i))
                    exit
                  endif
                enddo !k
              endif !ztmp

              if(k0==0) then
                write(*,*)'read_output7b_xyz: failed to find a vertical level:',it1,i,ifs,z2,ztmp(:,i)
!'
                stop
              endif

              if(.not.(it1==iday1.and.it1<irec_start.or.it1==iday2.and.it1>irec_end)) then
                icounter(i)=icounter(i)+1
                do m=1,ivs
                  tmp=out(i,k0,m)*(1-rat)+out(i,k0+1,m)*rat
                  residual(i,m)=residual(i,m)+tmp
                enddo !m
              endif !not
            endif !idry
          enddo !i=1,np
        endif !i23d
      enddo !it1=1,nrec

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      enddo !iday

!     Output
      do i=1,np
        if(icounter(i)==0) then
          residual(i,:)=0
        else
          residual(i,:)=residual(i,:)/icounter(i)
        endif !icounter
      enddo !i

      open(65,file='average.out')

      !For wave dir
!      file63=adjustl(file63)
!      leng=len_trim(file63)
      if(file63(1:len_file63).eq.'wwm_16.61') then
        !Read in sub-samples
        open(20,file='nodeflags.bp',status='old')
        read(20,*); read(20,*)

        eta2=residual(:,1) !temp save
        do i=1,np
          read(20,*)j,tmp,tmp,tmp2
          residual(i,1)=-sin(eta2(i)/180*pi)
          residual(i,2)=-cos(eta2(i)/180*pi)
          if(nint(tmp2)==0) residual(i,1:2)=0
        enddo !i
        ivs=2 !xyuv format
      endif !file63

      if(ivs==1) then
        write(65,*)
        write(65,*)ne,np
        do i=1,np
          write(65,*)i,x(i),y(i),residual(i,1)
        enddo !i
        do i=1,ne
          write(65,*)i,i34(i),elnode(1:i34(i),i)
        enddo !i
      else !vectors
        do i=1,np
          write(65,*)x(i),y(i),residual(i,1:2)
        enddo !i
      endif
      close(65)

      print*, 'Finished!'

      stop
      end
