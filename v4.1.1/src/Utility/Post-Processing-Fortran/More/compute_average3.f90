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
!	Read binary format v5.0 (SZ models) and compute Z0 and average surface vel.
!       Input: screen [elev.61 or hvel.64; iday1,iday2 etc]
!       Output: residual.out (*.61 or 62 type file)
!										*
!********************************************************************************
!     pgf90 -O2 -mcmodel=medium  -Bstatic -o compute_average3 compute_average3.f90

      program read_out
      parameter(nbyte=4)
!      parameter(mnp=150000)
!      parameter(mne=300000)
!      parameter(mnv=60)
      character*30 file63
      character*12 it_char
      character*48 start_time,version,variable_nm,variable_dim
      character*48 data_format
      allocatable :: ztot(:),sigma(:),x(:),y(:),dp(:),kbp00(:)
      allocatable :: out(:,:),icum(:,:,:),eta2(:),residual(:,:,:)
      integer,allocatable :: elnode(:,:)
      print*, 'Input file to read from (elev.61 or hvel.64):'
      read(*,'(a30)')file63
      
!      print*, 'Input file pre-fix index (start and end) to read:'
!      read(*,*) iday1,iday2
!      if(iday2<iday1) then
!        print*, 'end day < start day'
!        stop
!      endif

      print*, 'Input start and end output time in days:'
      read(*,*)start_day,end_day

      print*, 'Input output frequency in days:'
      read(*,*)stride_day

      print*, 'Input averaging window in days:'
      read(*,*)window
      window_sec=window*86400

      open(65,file='residual.out',access='direct',recl=nbyte)
      
!...  Header
!...
      open(63,file='1_'//file63,status='old',access='direct',recl=nbyte)
      irec=0
      do m=1,48/nbyte
        read(63,rec=irec+m) data_format(nbyte*(m-1)+1:nbyte*m)
        write(65,rec=irec+m) data_format(nbyte*(m-1)+1:nbyte*m)
      enddo
      if(data_format.ne.'DataFormat v5.0') then
        print*, 'This code reads only v5.0:  ',data_format
        stop
      endif
      irec=irec+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) version(nbyte*(m-1)+1:nbyte*m)
        write(65,rec=irec+m) version(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) start_time(nbyte*(m-1)+1:nbyte*m)
        write(65,rec=irec+m) start_time(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) variable_nm(nbyte*(m-1)+1:nbyte*m)
        write(65,rec=irec+m) variable_nm(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte
      do m=1,48/nbyte
        read(63,rec=irec+m) variable_dim(nbyte*(m-1)+1:nbyte*m)
        write(65,rec=irec+m) variable_dim(nbyte*(m-1)+1:nbyte*m)
      enddo
      irec=irec+48/nbyte
      irec2=irec !for output only

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

      print*, 'i23d=',i23d,' ; dtout=',dtout

!     Compute averaging window and output times
      iwindow=window_sec/dtout
      iwindow_half=iwindow/2 !# of records on either side
      print*, 'iwindow_half=',iwindow_half

!     Make output times at interger of dtout
      itime1=start_day*86400/dtout
      itime1=max(itime1,iwindow_half+1) !1st record is not t=0
      istride=stride_day*86400/dtout
      print*, 'Start output record # and stride=',itime1,istride

!     Compute 1st and last stack
      iday1=(itime1-iwindow_half-1)/nrec+1
      nout=(end_day*86400/dtout-iwindow_half-itime1)/istride+1 
      if(nout>200) stop 'nout>200'
      iday2=(itime1+(nout-1)*istride-1)/nrec+1
      print*, 'Start and last stack #=',iday1,iday2
      print*, 'nout=',nout
      print*, '1st and last output time in days=',itime1*dtout/86400,(itime1+(nout-1)*istride)*dtout/86400
      print*, 'Output times (days)=',((itime1+(i-1)*istride)*dtout/86400,i=1,nout)

      write(65,rec=irec2+1) nout
      write(65,rec=irec2+2) istride*dtout
      write(65,rec=irec2+3) 1
      write(65,rec=irec2+4) ivs
      write(65,rec=irec2+5) 2
      irec2=irec2+5

!     Vertical grid
      read(63,rec=irec+1) nvrt
      allocate(ztot(nvrt),sigma(nvrt))

      read(63,rec=irec+2) kz
      read(63,rec=irec+3) h_0
      read(63,rec=irec+4) h_s
      read(63,rec=irec+5) h_c
      read(63,rec=irec+6) theta_b
      read(63,rec=irec+7) theta_f
      irec=irec+7

      write(65,rec=irec2+1) nvrt
      write(65,rec=irec2+2) kz
      write(65,rec=irec2+3) h_0
      write(65,rec=irec2+4) h_s
      write(65,rec=irec2+5) h_c
      write(65,rec=irec2+6) theta_b
      write(65,rec=irec2+7) theta_f
      irec2=irec2+7

      do k=1,kz-1
        read(63,rec=irec+k) ztot(k)
        write(65,rec=irec2+k) ztot(k)
      enddo
      do k=kz,nvrt
        kin=k-kz+1
        read(63,rec=irec+k) sigma(kin)
        write(65,rec=irec2+k) sigma(kin)
      enddo
      irec=irec+nvrt
      irec2=irec2+nvrt

!     Horizontal grid
      read(63,rec=irec+1) np
      write(65,rec=irec2+1) np
      read(63,rec=irec+2) ne
      write(65,rec=irec2+2) ne
      irec=irec+2
      irec2=irec2+2

      !Allocate
      allocate(x(np),y(np),dp(np),kbp00(np),elnode(4,ne),out(nvrt,2),icum(np,nvrt,2),eta2(np),residual(2,nout,np),stat=istat)
      if(istat/=0) stop 'Failed to alloc (1)'

      do m=1,np
        read(63,rec=irec+1)x(m)
        read(63,rec=irec+2)y(m)
        read(63,rec=irec+3)dp(m)
        read(63,rec=irec+4)kbp00(m)
        write(65,rec=irec2+1)x(m)
        write(65,rec=irec2+2)y(m)
        write(65,rec=irec2+3)dp(m)
        write(65,rec=irec2+4)kbp00(m)
        irec=irec+4
        irec2=irec2+4
      enddo !m=1,np
      do m=1,ne
        read(63,rec=irec+1)i34
        write(65,rec=irec2+1)i34
        irec=irec+1
        irec2=irec2+1
        do mm=1,i34
          read(63,rec=irec+1)elnode(mm,m)
          write(65,rec=irec2+1)elnode(mm,m)
          irec=irec+1
          irec2=irec2+1
        enddo !mm
      enddo !m
      irec0=irec
      irec02=irec2

      print*, 'last element',(elnode(j,ne),j=1,3)

!...  Compute relative record # for a node and level
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
      residual=0
      it_tot=(iday1-1)*nrec !time record #
      do iday=iday1,iday2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      write(it_char,'(i12)')iday
      open(63,file=it_char//'_'//file63,status='old',access='direct',recl=nbyte)

      irec=irec0
      do it1=1,nrec
        read(63,rec=irec+1) time
        read(63,rec=irec+2) it
        irec=irec+2
        it_tot=it_tot+1

        print*, 'time (days),it_tot=',time/86400,it_tot

        do i=1,np
          read(63,rec=irec+i) eta2(i)
        enddo !i
        irec=irec+np

        if(i23d.eq.2) then
!          do m=1,ivs
!            read(63,rec=irec+(node-1)*ivs+m) out(1,m)
!          enddo !m
          irec=irec+np*ivs
!         Z0
          do i=1,np
            !Search for all output times
            do j=1,nout
              ileft=itime1+(j-1)*istride-iwindow_half !left end pt for avering
              iright=itime1+(j-1)*istride+iwindow_half
              if(ileft>it_tot) exit
            
              if(it_tot>=ileft.and.it_tot<=iright) residual(1,j,i)=residual(1,j,i)+eta2(i)/(2*iwindow_half+1)
            enddo !j
          enddo !i
        else !i23d=3 
          do i=1,np
            do k=max0(1,kbp00(i)),nvrt
              do m=1,ivs
                read(63,rec=irec+1) out(k,m)
                irec=irec+1
              enddo !m
            enddo !k

            !Average
            do j=1,nout
              ileft=itime1+(j-1)*istride-iwindow_half !left end pt for avering
              iright=itime1+(j-1)*istride+iwindow_half
              if(ileft>it_tot) exit

              if(it_tot>=ileft.and.it_tot<=iright) then
                residual(1:2,j,i)=residual(1:2,j,i)+out(nvrt,1:2)/(2*iwindow_half+1)
              endif
            enddo !j
          enddo !i
        endif !i23d
      enddo !it1=1,nrec

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      enddo !iday

!     Output
      do it=1,nout
        write(65,rec=irec2+1) (itime1+(it-1)*istride)*dtout
        write(65,rec=irec2+2) it
        irec2=irec2+2
        do i=1,np
          write(65,rec=irec2+i) 0. 
        enddo !i
        irec2=irec2+np
        if(i23d.eq.2) then
          do i=1,np
            write(65,rec=irec2+1) residual(1,it,i)
            irec2=irec2+1
          enddo !i
        else
          do i=1,np
            write(65,rec=irec2+1) residual(1,it,i)
            write(65,rec=irec2+2) residual(2,it,i)
            irec2=irec2+2
          enddo !i
        endif
      enddo !it

      stop
      end
