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
!       at a particular level.
!       Input: *.6?; iday1,iday2; level index (lev0)
!       Output: average.out (gredit for scalar or xmgr5 format for vector)
!										*
!       pgf90 -O2 -mcmodel=medium  -Bstatic -o compute_average2 compute_average2.f90
!********************************************************************************
!
      program read_out
      parameter(nbyte=4)
      parameter(mnp=200000)
      parameter(mne=400000)
      parameter(mnv=50)
      character*30 file63
      character*12 it_char
      character*48 start_time,version,variable_nm,variable_dim
      character*48 data_format
      inteter :: elnode(4,mne)
      dimension ztot(mnv),sigma(mnv),x(mnp),y(mnp),dp(mnp),kbp00(mnp) !,kfp(mnp)
      dimension out(mnv,2),icum(mnp,mnv,2),eta2(mnp)
      dimension residual(mnp,2)
      
      print*, 'Input file to read from (without *_):'
      read(*,'(a30)')file63
      
      print*, 'Input file pre-fix indices (start and end) to read:'
      read(*,*) iday1,iday2
      if(iday2<iday1) then
        print*, 'end day < start day'
        stop
      endif

      print*, 'Input level # for 3D variables:'
      read(*,*) lev0

      
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

      print*, 'i23d=',i23d,' vpos=',vpos

!     Vertical grid
      read(63,rec=irec+1) nvrt
      read(63,rec=irec+2) kz
      read(63,rec=irec+3) h_0
      read(63,rec=irec+4) h_s
      read(63,rec=irec+5) h_c
      read(63,rec=irec+6) theta_b
      read(63,rec=irec+7) theta_f
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
        read(63,rec=irec+k) sigma(kin)
      enddo
      irec=irec+nvrt

!     Horizontal grid
      read(63,rec=irec+1) np
      read(63,rec=irec+2) ne
      irec=irec+2
      if(np.gt.mnp.or.ne.gt.mne) then
        write(*,*)'Too many nodes/elements',np,ne
        stop
      endif
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
          read(63,rec=irec+1)elnode(mm,m)
          irec=irec+1
        enddo !mm
      enddo !m
      irec0=irec

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
      ndays=iday2-iday1+1
      do iday=iday1,iday2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      write(it_char,'(i12)')iday
      open(63,file=it_char//'_'//file63,status='old',access='direct',recl=nbyte)

      irec=irec0
      do it1=1,nrec
        read(63,rec=irec+1) time
        read(63,rec=irec+2) it
        irec=irec+2

        print*, 'time=',time/86400

        do i=1,np
          read(63,rec=irec+i) eta2(i)
        enddo !i
        irec=irec+np

        if(i23d.eq.2) then
          do i=1,np
            do m=1,ivs
              read(63,rec=irec+1) out(1,m)
              irec=irec+1
            enddo !m
            residual(i,1:ivs)=residual(i,1:ivs)+out(1,1:ivs)/ndays/nrec
          enddo !i
        else !i23d=3 
          do i=1,np
            do k=max0(1,kbp00(i)),nvrt
              do m=1,ivs
                read(63,rec=irec+1) out(k,m)
                irec=irec+1
              enddo !m
            enddo !k
            residual(i,1:ivs)=residual(i,1:ivs)+out(lev0,1:ivs)/ndays/nrec
          enddo !i
        endif !i23d
      enddo !it1=1,nrec

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      enddo !iday

!     Output
      open(65,file='average.out')
      if(ivs==1) then
        write(65,*)
        write(65,*)ne,np
        do i=1,np
          write(65,*)i,x(i),y(i),residual(i,1)
        enddo !i
        do i=1,ne
          write(65,*)i,3,(elnode(1,i),l=1,3)
        enddo !i
      else !vectors
        do i=1,np
          write(65,*)x(i),y(i),residual(i,1:2)
        enddo !i
      endif
      close(65)

      stop
      end
