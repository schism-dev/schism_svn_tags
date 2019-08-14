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
!****************************************************************************************

!	Read binary format v5.0 (hybrid S-Z models) for multiple files at all nodes	
!       and output at a time and a fixed z-level (3D variables).
!       Works for mixed tri/quad outputs.
!       Inputs: screen.out; binary file; vgrid.in (in this dir or ../)
!       Outputs: extract_slab.gr3 (magni. for vectors)
!       History: (1) added non-standard outputs (April 2012) - transparent to most scripts
!               as format is same; (2) added ivcor=1 (Dec 2013); (3)
!               added quads (Nov. 2014)
!****************************************************************************************
!     pgf90 -O2 -mcmodel=medium -Mbounds -Bstatic -o read_output7_allnodes_slab read_output7_allnodes_slab.f90 ~/SELFE/svn/trunk/src/Utility/UtilLib/compute_zcor.f90
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
      allocatable :: zs(:,:),ze(:,:),idry(:),outs(:,:,:),out3(:,:)
      allocatable :: ic3(:,:),isdel(:,:),isidenode(:,:),kbp(:),sigma_lcl(:,:)
      
      print*, 'Input file to read from (without *_; e.g. salt.63):'
      read(*,'(a30)')file63
      
      print*, 'Input time out in sec (MUST be one of the output step!):'
      read(*,*)timeout 

      print*, 'Input z coord. for output (<0 for below MSL):'
      read(*,*)zout 

      open(65,file='extract_slab.gr3',status='replace')
      
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

      print*, 'i23d=',i23d,' nrec= ',nrec
     
!     Calc stacks
      iday2=int((timeout-dtout)/nrec/dtout+1.e-6)+1
      if(iday2<0) then
        print*, 'timeout is before 1st record:',timeout
        stop
      endif
      print*, 'Found stack #:',iday2

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
     !&icum(np,nvrt,2),eta2(np),stat=istat)
     &eta2(np),xctr(ne),yctr(ne),dpe(ne),kbe(ne),ztmp(nvrt,np),ze(nvrt,ne), &
     &out3(2,np),stat=istat)
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
      do iday=iday2,iday2
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      write(it_char,'(i12)')iday
      open(63,file=it_char//'_'//file63,status='old',access='direct',recl=nbyte)

      irec=irec0
      do it1=1,nrec
        read(63,rec=irec+1) time
        read(63,rec=irec+2) it
        irec=irec+2

        print*, 'time in days =',time/86400

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

          !Output: time, 2D variable at all nodes
          if(abs(timeout-time)<1.e-5) then
            write(65,*)'Day: ',time/86400
            write(65,*)ne,np
            do i=1,np
              if(ivs==1) then
                write(65,*)i,x(i),y(i),out(i,1,1)
              else
                write(65,*)i,x(i),y(i),sqrt(out(i,1,1)**2+out(i,1,2)**2)
              endif
            enddo !i
            do i=1,ne
              write(65,*)i,i34(i),elnode(1:i34(i),i)
            enddo !i
            stop
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
              endif !wet/dry
            else if(ivcor==2) then !SZ
              call zcor_SZ(dp(i),eta2(i),h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot, &
     &sigma,ztmp(:,i),idry(i),kbpl)
            endif
          enddo !i=1,np

          !Output
          if(abs(timeout-time)<1.e-5) then
            print*, 'Final output time (days):',time/86400
            rjunk=-9999. !fill in value 
            do i=1,np
              if(idry(i)==1) then
                out3(1:ivs,i)=rjunk
              else !wet
                !Vertical interp
                kbpl=kbp(i)
                if(zout>=ztmp(nvrt,i)) then !above F.S.
                  !k0=nvrt-1; rat=1
                  out3(1:ivs,i)=out(i,nvrt,1:ivs)
                else if(zout<=ztmp(kbpl,i)) then !below bottom
                  !k0=kbpl; rat=0
                  out3(1:ivs,i)=rjunk
                else !above bottom; cannot be above F.S.
                  k0=0
                  do k=kbpl,nvrt-1
                    if(zout>=ztmp(k,i).and.zout<=ztmp(k+1,i)) then
                      k0=k
                      rat=(zout-ztmp(k,i))/(ztmp(k+1,i)-ztmp(k,i))
                      exit
                    endif
                  enddo !k

                  if(k0==0) then
                    write(*,*)'read_output7b_xyz: failed to find a vertical level:',it1,i,ifs,zout,ztmp(:,i)
!'
                    stop
                  endif

                  do m=1,ivs
                    out3(m,i)=out(i,k0,m)*(1-rat)+out(i,k0+1,m)*rat
                  enddo !m
                endif !ztmp
              endif !idry
            enddo !i=1,np

            write(65,*)'Day: ',time/86400
            write(65,*)ne,np
            do i=1,np
              if(ivs==1) then
                write(65,*)i,x(i),y(i),out3(1,i)
              else
                write(65,*)i,x(i),y(i),sqrt(out3(1,i)**2+out3(2,i)**2)
              endif
            enddo !i=1,np

            do i=1,ne
              write(65,*)i,i34(i),elnode(1:i34(i),i)
            enddo !i
            stop

          endif !abs
        endif !i23d
      enddo !it1=1,nrec

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      enddo !iday=

      print*, 'Finished!'

      stop
      end
