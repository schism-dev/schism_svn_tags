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
!											*
!	read_output7b with (x,y,time) read in from station.xyt (time in sec); e.g., casts; 
!       for 3D variables (surface values for 2D variables). Interpolation in time, and
!       add extra times before and after to examine phase errors.
!       Inputs: (1) binary files;
!               (2) station.xyt: make sure all times (in sec) are after 1st record (to ensure interpolation in time); 
!                                 pad extra days before and after if necessary.
!               (3) read_output7b_xyt.in: 
!                   file63\n 
!                   invalid value (for out of domain, dry etc)\n 
!                   start corie day (0 usually) \n
!                   window (hours) b4 and after the cast, stride (hrs) - used to 
!                   examine the phase error. If window=0, each cast is repeated twice 
!                   there are no other extra casts.
!                (4) vgrid.in (in this dir or ../)
!       Outputs: fort.1[89]; fort.11 (fatal errors); fort.12: nonfatal errors.
!                The total # of 'virtual' casts for each actual cast is 2*window/stride+2
!									
!       pgf90 -O2 -mcmodel=medium -Mbounds -Bstatic -o read_output7b_xyt read_output7b_xyt.f90 ~/SELFE/svn/trunk/src/Utility/UtilLib/compute_zcor.f90
!****************************************************************************************
      program read_output_xyt
      character*30  :: file    !< base file name (e.g. elev.61) of files to be read 
      character*80  :: stationfile  !< station file, e.g. station.xyzt
      character*80  :: outfile  !< output file name
      character*80  :: vgrid   !< vertical grid file 
      real          :: fill    !< fill value for na values
      real          :: window
      real          :: stride

      call read_output7b_xyt_input(file,vgrid,stationfile,outfile,fill,window,stride)
      call read_output7b_xyt(file,vgrid,stationfile,outfile,fill,window,stride)
     
      end program


      subroutine read_output7b_xyt_input(file,vgrid,stationfile,outfile,fill,window,wstride)
      use argparse
      implicit none
      character(LEN=*)  :: file    !< base file name (e.g. elev.61) of files to be read 
      character*80      :: infile  !< input file containing options for this routine
      character(LEN=*)  :: stationfile  !< input file, e.g. station.xyzt
      character(LEN=*)  :: vgrid 
      character(LEN=*)  :: outfile
      real(kind=4)  :: fill    !< fill value for na values
      real          :: window
      real          :: wstride
      integer       :: comcount      
 
      cmd_name = "read_output_xyt"
      call cla_init(cmd_name,"Read vertical profiles of model output at particular times and locations as with a CTD cast. " //&
                             "Purpose of the window is to allow for small phase errors by grabbing a few time points " // &
                             "before and after the cast.The total # of 'virtual' casts is 2*window/stride+2." // &
                             "Note that for window = 0 this produces two (redundant) casts exactly at the requested time.")
      call cla_register('-f','--file', 'base input file (e.g. salt.63) without block number prefex', cla_char,'')
      call cla_register('-v','--vgrid', 'path to vertical grid file (typically vgrid.in or ../vgrid.in', cla_char,'vgrid.in')
      call cla_register('-o','--out', 'path to output file', cla_char,'extract.out')
      call cla_register('-i','--in', 'input file containing options (possibly overriden at command line)',cla_char,'')
      call cla_register('-s','--stations','input file listing stations (first line is #, second is totoal number of stations,then id x y t)',cla_char,'station.xyt')
      call cla_register('-x','--fill', 'fill value for invalid data', cla_float,'-999999.')
      call cla_register('-w','--window','window around selected times over which to select values', cla_float ,'0')
      call cla_register('-r','--stride','stride window over which to select values due to possible phase differences', cla_float  ,'1')
      call cla_validate

      comcount = command_argument_count()
      if (comcount == 0) then
          infile = "read_output7b_xyt.in"
      else
          call cla_get("--in",infile)
      end if
      
      if (len_trim(infile)>1) then
      open(10,file=trim(infile),status='old')
     
!      print*, 'Input file to read from (without *_):'
      read(10,*)file
!     Fill number used for invalid 3D variables: below bottom; dry spot; no parents
!      print*, 'Input values to be used for invalid place:'
      read(10,*)fill
      read(10,*)window,wstride !in hours
      close(10)
      ! override with (explicitly provided) command line values)
      if (cla_key_present("--file"))   call cla_get("--file",file)
      if (cla_key_present("--fill"))   call cla_get("--fill",fill)
      if (cla_key_present("--window")) call cla_get("--window",window)
      if (cla_key_present("--stride")) call cla_get("--stride",wstride)
      else
      call cla_get("--file",file)
      call cla_get("--fill",fill)
      call cla_get("--window",window)
      call cla_get("--stride",wstride)
      end if
      call cla_get("--stations",stationfile)
      call cla_get("--vgrid",vgrid)
      call cla_get("--out",outfile)


      end subroutine



      subroutine read_output7b_xyt(file63,vgrid,stationfile,outfile,fill,window,wstride)
      parameter(nbyte=4) 
      character(LEN=*) :: vgrid
      character(LEN=*) :: stationfile
      character(LEN=*) :: outfile
      real             :: window
      real             :: wstride
      real             :: start_corie
      character(LEN=*) :: file63
      character*12 it_char
      character*48 start_time,version,variable_nm,variable_dim
      character*48 data_format
      integer,allocatable :: elnode(:,:)
      dimension swild(3),out3(2,2)
      allocatable :: sigma(:),cs(:),ztot(:),x(:),y(:),dp(:),kbp00(:),kfp(:)
      allocatable :: out(:,:,:),out2(:,:,:),icum(:,:,:),eta2(:),node3(:,:),arco(:,:)
      allocatable :: ztmp(:),x00(:),y00(:),iep(:),out4(:,:),t00(:)
      allocatable :: iday(:,:),irecord(:,:),times(:,:)
      allocatable :: sigma_lcl(:,:),kbp(:),ztmp2(:,:)
      rjunk=fill
      start_corie = 0.0
      if(wstride==0) stop 'wstride=0'
      nextra=2*window/wstride+1 !extra casts in addition to each cast (for phase error)
      
      if ( window ==0 ) then
          nextra= 0
      endif
     
      
      
      print*,file63,vgrid,stationfile,outfile,fill,window,wstride
      
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
      
!     Read in station.xyt
      open(10,file=trim(stationfile),status='old')
      read(10,*) 
      read(10,*) nxy
      nxy2=nxy*(1+nextra)
      allocate(x00(nxy2),y00(nxy2),t00(nxy2),stat=istat)
      if(istat/=0) stop 'Falied to allocate (1)'

      do i=1,nxy
        read(10,*)k,xtmp,ytmp,ttmp
!       Check if time is before first record
        if(ttmp<dtout) then
          write(11,*)'Time before first record; try to pad extra day:',i,ttmp
          stop
        endif

        indx=(i-1)*(1+nextra)+1
        x00(indx)=xtmp
        y00(indx)=ytmp
        t00(indx)=ttmp !time in sec
        !Add extra casts
        do j=1,nextra
          indx=indx+1
          x00(indx)=xtmp
          y00(indx)=ytmp
          t00(indx)=max(dtout,ttmp-window*3600+(j-1)*wstride*3600) !also need to ensure it's <last time
        enddo !j
       
      enddo !i
      close(10)
      nxy=nxy2
      
!      print*, 'i23d=',i23d,' nrec= ',nrec

!     Vertical grid (obsolete)
      read(63,rec=irec+1) nvrt
!      read(63,rec=irec+2) kz
      read(63,rec=irec+3) h0
!      read(63,rec=irec+4) h_s
!      read(63,rec=irec+5) h_c
!      read(63,rec=irec+6) theta_b
!      read(63,rec=irec+7) theta_f
!      irec=irec+7
!      allocate(ztot(nvrt),sigma(nvrt),cs(nvrt),ztmp(nvrt),out(3,nvrt,2), &
!     &out2(2,nvrt,2),out4(nvrt,2),stat=istat)
!      if(istat/=0) stop 'Falied to allocate (2)'
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
      irec=irec+nvrt+7
      read(63,rec=irec+1) np
      read(63,rec=irec+2) ne
      irec=irec+2
      allocate(x(np),y(np),dp(np),kbp00(np),kfp(np),elnode(4,ne), &
     &icum(np,nvrt,2),eta2(np),node3(nxy,3),arco(nxy,3),iep(nxy),iday(2,nxy), &
     &irecord(2,nxy),times(2,nxy),stat=istat)
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
          read(63,rec=irec+1)elnode(mm,m)
          irec=irec+1
        enddo !mm
      enddo !m
      irec0=irec

!      print*, 'last element',(elnode(j,ne),j=1,3)

!     Read in vgrid.in
      allocate(ztot(nvrt),sigma(nvrt),sigma_lcl(nvrt,np),kbp(np))
      if(nvrt>2) call get_vgrid(trim(vgrid),np,nvrt,ivcor,kz,h_s,h_c,theta_b,theta_f,ztot,sigma,sigma_lcl,kbp)

      allocate(ztmp(nvrt),ztmp2(nvrt,3),out(3,nvrt,2),out2(2,nvrt,2),out4(nvrt,2),stat=istat)
      if(istat/=0) stop 'Falied to allocate (2)'

!...  Find parent element for (x00,y00)
      iep=0
      arco=1./3 !initialize for pts without parents
      do l=1,nxy
        node3(l,1:3)=elnode(1:3,1) !initialize for pts without parents
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
            n0=elnode(j,i)
            n1=elnode(j1,i)
            n2=elnode(j2,i)
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
            node3(l,1:3)=elnode(1:3,i)
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

      if (ivs == 1)then
          open(18,file=trim(outfile),status='unknown')
      else if (ivs==2) then
          open(18,file="1_"//trim(outfile),status='unknown')
          open(19,file="2_"//trim(outfile),status='unknown')
      else
          stop "Problem reading dimensionality of variable (ivs) from file"
      end if

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
!        if(t00(i)<dtout) then
!          write(11,*)'Time before first record; try to padd extra day (0):',i,t00(i)
!!'
!          stop
!        endif

!       Lower and upper bound days and record #s for t00(i)
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

!         Only read in 3 nodes for efficiency
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
              !Compute z-coordinates
              if(ivcor==1) then !PWS
                do j=1,3
                  nd=node3(i,j)
                  do k=kbp(nd)+1,nvrt-1
                    ztmp2(k,j)=(eta2(nd)+dp(nd))*sigma_lcl(k,nd)+eta2(nd)
                  enddo !k
                  ztmp2(kbp(nd),j)=-dp(nd) !to avoid underflow
                  ztmp2(nvrt,j)=eta2(nd) !to avoid underflow
                enddo !j

                ztmp=0
                kbpl=minval(kbp(node3(i,1:3)))
                do k=kbpl,nvrt
                  do j=1,3
                    nd=node3(i,j)
                    ztmp(k)=ztmp(k)+arco(i,j)*ztmp2(max(k,kbp(nd)),j)
                  enddo !j
                enddo !k

              else if(ivcor==2) then !SZ
                call zcor_SZ(dep,etal,h0,h_s,h_c,theta_b,theta_f,kz,nvrt,ztot, &
     &sigma,ztmp(:),idry2,kbpl)
              endif

              if(1==2) then
!----------------------------------------
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
!------------------------------------------
              endif !1==2
       
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
!              k0=0
!              do k=kbpl,nvrt-1
!                if(ztmp(nvrt)-z00(i)>=ztmp(k).and.ztmp(nvrt)-z00(i)<=ztmp(k+1)) then
!                  k0=k
!                  rat=(ztmp(nvrt)-z00(i)-ztmp(k))/(ztmp(k+1)-ztmp(k))
!                  exit
!                endif
!              enddo !k
!              if(k0==0) then
!                out3(:,:)=rjunk
!                exit loop1
!               write(12,*)'Warning: failed to find a vertical level:',it,i
!              else
!                do m=1,ivs
!                  out3(l,m)=out2(l,k0,m)*(1-rat)+out2(l,k0+1,m)*rat
!                enddo !m
!              endif
            endif !dry/wet
          endif !i23d
        enddo loop1 !l=1,2; 2 times

!       Interpolate in time
        trat=(t00(i)-times(1,i))/(times(2,i)-times(1,i)) !must be [0,1]
        if(i23d==2) then
          if(iep(i)==0) then !no parents
            out4(1,1:ivs)=rjunk
          else
            out4(1,1:ivs)=out2(1,1,1:ivs)*(1-trat)+out2(2,1,1:ivs)*trat
          endif
          write(18,'(e16.8,2(1x,f12.3))')t00(i)/86400.,out4(1,1:ivs)
        else !3D
          if(iep(i)==0) then !no parents
            out4(:,1:ivs)=rjunk
          else
            out4(kbpl:nvrt,1:ivs)=out2(1,kbpl:nvrt,1:ivs)*(1-trat)+out2(2,kbpl:nvrt,1:ivs)*trat
          endif

          do k=nvrt,kbpl,-1
            !First of each cast suite is at the actual cast time (followed by b4 and after)
            write(18,'(i6,4(1x,f12.3))')i,out4(k,1),ztmp(k)-ztmp(nvrt),ztmp(k),t00(i)/86400.
            if(ivs==2) write(19,'(i6,4(1x,f12.3))')i,out4(k,2),ztmp(k)-ztmp(nvrt),ztmp(k),t00(i)/86400.
          enddo !k
        endif
      enddo !i=1,nxy
      
      
      deallocate(x00,y00,t00,stat=istat)
      deallocate(x,y,dp,kbp00,kfp,elnode,icum,eta2,node3,arco,iep,iday,irecord,times,stat=istat)
      deallocate(ztot,sigma,sigma_lcl,kbp,stat=istat)
      deallocate(ztmp,ztmp2,out,out2,out4,stat=istat)
      close(18)
      if (i23d/=2) then
          close(19)
      endif
      stop
      end

      function signa(x1,x2,x3,y1,y2,y3)
!...  Compute signed area formed by pts 1,2,3

      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2

      return
      end

