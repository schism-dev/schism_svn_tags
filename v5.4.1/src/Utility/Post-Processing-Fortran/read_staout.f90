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


!     Read 3D station outputs (hvel, S,T etc) from SELFE and output at a level at a station
!     This simple script shows the format for 3D station outputs.
!     Will extrapolate below bottom/above surface.
!     Inputs: outputs/staout_*; read_staout.in
!     Output: fort.18

!     pgf90 -O2 -mcmodel=medium  -Bstatic -o read_staout read_staout.f90

      character(len=40) :: fname
      allocatable :: sta_out_gb(:),sta_out3d_gb(:,:),zta_out3d_gb(:,:)

!     read_staout.in
      open(9,file='read_staout.in',status='old')
!      print*, 'Input file name (staout_[5-9]):'
      read(9,'(a40)') fname
!      print*, 'Input runtime in days:'
      read(9,*) rndays
!      print*, 'Input output time step in sec:'
      read(9,*) dtout
!      print*, 'Input total # of stations:'
      read(9,*) nsta 
!      print*, 'Input station index to extract:'
      read(9,*) ista
!      print*, 'Input total # of levels (nvrt):'
      read(9,*) nvrt
!      print*, 'Is the z input relative to F.S. (1) or MSL (0)?'
      read(9,*) ifs
!      print*, 'Input z:'
      read(9,*) z00 !If ifs=0, z00 generally <0
      if(ifs==1.and.z00<0) stop 'z00<0'

      allocate(sta_out_gb(nsta),sta_out3d_gb(nvrt,nsta),zta_out3d_gb(nvrt,nsta))

      open(10,file='outputs/'//fname,status='old')
      nstep=rndays*86400/dtout+0.1
      do it=1,nstep
        read(10,*)time,sta_out_gb(:)
        read(10,*)time,sta_out3d_gb(:,:),zta_out3d_gb(:,:)

!       Do vertical interpolation
        if(ifs==1) then
          zfs=z00
        else
          zfs=zta_out3d_gb(nvrt,ista)-z00
        endif

        if(zfs>=zta_out3d_gb(nvrt,ista)-zta_out3d_gb(1,ista)) then
          k0=1; zrat=0
        else if(zfs<=0) then
          k0=nvrt-1; zrat=1
        else 
          k0=0
          do k=1,nvrt-1
            zz=zta_out3d_gb(nvrt,ista)-zfs
            if(zz>=zta_out3d_gb(k,ista).and.zz<=zta_out3d_gb(k+1,ista)) then
              k0=k
              zrat=(zz-zta_out3d_gb(k,ista))/(zta_out3d_gb(k+1,ista)-zta_out3d_gb(k,ista))
              exit
            endif
          enddo !k
          if(k0==0) then
            print*, 'Cannot find a level:',it,zfs,zta_out3d_gb(:,ista)
          endif
        endif !zfs

        varout=sta_out3d_gb(k0,ista)*(1-zrat)+sta_out3d_gb(k0+1,ista)*zrat
        write(18,*)time/86400,varout
      enddo !it
      close(10)
       
      stop
      end
