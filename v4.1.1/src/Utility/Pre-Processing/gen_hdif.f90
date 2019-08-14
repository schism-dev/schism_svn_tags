!     Generate hdif.gr3 (or hvis.gr3) with hdif=area * gamma/dt (where 
!     gamma is a dimensionless constant between [0,0.5))
!     Inputs: screen (gamma1,dt; hdif_max); hgrid.gr3
!     Output: hdif.gr3

!     ifort -Bstatic -O3 -o gen_hdif gen_hdif.f90
!     pgf90 -O2 -mcmodel=medium -O2 -o gen_hdif gen_hdif.f90

      implicit real*8(a-h,o-z)
      parameter(mnp=200000)
      parameter(mne=400000)
      parameter(mnei=30)
      integer :: elnode(3,mne)
      dimension x(mnp),y(mnp),dp(mnp),area(mne)
      dimension nne(mnp),indel(mnei,mnp),hdif(mnp),hdif_e(mne)

      print*, 'Input gamma \in [0,0.5) and time step:' !stricter than necessary as dt' in transport is smaller
      read*, gamma1,dt
      if(gamma1>=0.5.or.gamma1<0) then
        print*, 'Choose a smaller gamma: ',gamma1
        stop
      endif
      print*, 'Input max. hdif:'
      read*, hdif_max

      open(14,file='hgrid.gr3')
      open(13,file='hdif.gr3')
      read(14,*)
      read(14,*) ne,np
      if(ne>mne.or.np>mnp) then
        write(*,*)'Increase mne/mnp',mne,mnp,ne,np
        stop
      endif

      do i=1,np
        read(14,*) j,x(i),y(i) !,dp(i)
      enddo !i=1,np
      do i=1,ne
        read(14,*) j,l,(elnode(k,i),k=1,l)
        if(l/=3) then
          write(*,*)'SELFE cannot handle quads:',i
          stop
        endif
        n1=elnode(1,i)
        n2=elnode(2,i)
        n3=elnode(3,i)
        area(i)=signa(x(n1),x(n2),x(n3),y(n1),y(n2),y(n3))
        if(area(i)<=0) then
          write(*,*)'Negative area at',i
          stop
        endif
!        write(13,*)i,dmin1(hdif_max,gamma1*area(i))
        hdif_e(i)=gamma1*area(i)/dt
      enddo !i=1,ne      

!     Neighborhood
      nne=0
      do i=1,ne
        do j=1,3
          nd=elnode(j,i)
          nne(nd)=nne(nd)+1
          if(nne(nd)>mnei) then
            write(*,*)'Too many neighbors',nd
            stop
          endif
          indel(nne(nd),nd)=i
        enddo
      enddo

      hdif=0
      do i=1,np
        do j=1,nne(i)
          ie=indel(j,i)
          hdif(i)=hdif(i)+hdif_e(ie)/nne(i)
        enddo !j
          hdif(i)=dmin1(hdif(i),hdif_max)
      enddo !i

      write(13,*)'gamma1= ',real(gamma1), ' ,dt= ',real(dt),' ,hdif_max= ',real(hdif_max)
!'
      write(13,*)ne,np
      do i=1,np
        write(13,*)i,real(x(i)),real(y(i)),real(hdif(i))
      enddo !i
      do i=1,ne
        write(13,*) i,3,elnode(1:3,i)
      enddo !i

      print*, 'Min. diffusivity/viscosity= ',minval(hdif(1:np))

      stop
      end

      function signa(x1,x2,x3,y1,y2,y3)
!...  Compute signed area formed by pts 1,2,3
      implicit real*8(a-h,o-z)

      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2

      return
      end

