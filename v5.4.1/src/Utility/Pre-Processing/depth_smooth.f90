!      subroutine depth_smooth_hannah2
!
!	this program smooths depths according to the Hannah-Wright ratio
!	of (delta h)/h_min for all nodes in a *.gr3 file
!	This ratio is calculated over each element (delta h= h_max - h_min)
!	and compared with a prescribed
!	upper limit.  (h is taken as min  depth in the elem) If the ratio exceeds the limit, the largest and
!	smallest depths are made smaller and larger repectively by
!	2% of the larger value. This adjustment (approx) conserves the volume
!	of each triangle
!	This procedure is iterated until all elements satisfy the ratio test
!	or the max. change in depths among all nodes during two consecutive iterations 
!       is less than a prescribed value (diff_h_max).
!       The final *.ngh file is output
!	This smoothing is only done inside a prescribed polygon

!zyl
!     Input: hgrid.old (mixed tri/quads); region.gr3; screen inputs: upper,dpmin,diff_h_max 
!     Output: hgrid.new (node order & connectivity unchanged); depth_diff.gr3 (difference in depths)

!     ifort -Bstatic -O2 -o depth_smooth depth_smooth.f90
!     pgf90 -O2 -mcmodel=medium -Mbounds -Bstatic -o depth_smooth depth_smooth.f90
      program smooth
      implicit real*8(a-h,o-z)
      parameter(mnp=2000000)
      parameter(mne=4000000)
      dimension x(mnp),y(mnp),dp(mnp),dp0(mnp),dp00(mnp)
      dimension i34(mne),nm(mne,4),iregion(mnp)

!      upper=0.1 !upper limit of the ratio
!      dpmin=200. !cut-off depth (shallower places will not be changed)
!      diff_h_max=1 !the max. change in depths among all nodes during two consecutive iterations
      print*, 'Input upper,dpmin, diff_h_max:'
      read*, upper,dpmin,diff_h_max
      if(upper<0.or.dpmin<0.or.diff_h_max<=0) then
        write(*,*)'Wrong input for upper,dpmin,diff_h_max'
        stop
      endif

      open(14,file='hgrid.old',status='old')
      open(13,file='region.gr3',status='old')
!     Read in grid file
      read(14,*)
      read(14,*)ne,np
      if(np.gt.mnp.or.ne.gt.mne) then
        write(*,*)'Increase mnp/mne'
        stop
      endif
      do i=1,np
        read(14,*)j,x(i),y(i),dp(i)
        dp00(i)=dp(i)
        dp0(i)=dp(i)
      enddo !i
      do i=1,ne
        read(14,*)j,i34(i),nm(i,1:i34(i))
      enddo !i
      close(14)
      
!     read in the region where the smoothing will occur
      read(13,*)
      read(13,*)
      do i=1,np
        read(13,*)j,xtmp,ytmp,tmp
        iregion(i)=tmp
        if(iregion(i).ne.0.and.iregion(i).ne.1) then
          write(*,*)'Wrong region flag:',i,iregion(i)
          stop
        endif
      enddo !i
      close(13)

!     Iterate calculating (delta h)/h for each triangle
      iter=0
16    continue
      iter=iter+1
      if(iter.gt.10000) then
        write(*,*)'too many iterations'
        stop
      endif
      rmax=0
      do 13 i=1,ne
        hmax=-100
        hmin=1.e25
        do j=1,i34(i) !nodes
          nd=nm(i,j)
          if(iregion(nd).eq.0) go to 13
          if(dp(nd).gt.hmax) then
            nmax=nd
            hmax=dp(nd)
          endif
          if(dp(nd).lt.hmin) then
            nmin=nd
            hmin=dp(nd)
          endif
        enddo !j
        if(hmin.lt.dpmin) go to 13

        if(hmax.lt.hmin) then
          write(*,*)'Weird depths:',hmax,hmin
          stop
        endif
        ratio=(hmax-hmin)/hmin
        if(ratio.gt.rmax) rmax=ratio
        if(ratio.gt.upper) then
          change=0.02*hmax
          dp(nmax)=dp(nmax)-change 
          dp(nmin)=dp(nmin)+change 
        endif                
13     continue !i=1,ne  
       write(*,*)'iter= ',iter,' max. ratio= ',rmax

       do i=1,np
         if(dabs(dp(i)-dp0(i)).ge.diff_h_max) then
           do j=1,np
             dp0(j)=dp(j)
           enddo !j
           go to 16
         endif        
       enddo !i

!     Output
      write(*,*)'Iteration converged in ',iter
      open(16,file='hgrid.new')
      write(16,*)upper,dpmin,diff_h_max
      write(16,*)ne,np
      do i=1,np
        write(16,'(i8,2(1x,e20.12),1x,f9.3)')i,x(i),y(i),dp(i)
      enddo !i
      do i=1,ne
        write(16,*)i,i34(i),nm(i,1:i34(i))
      enddo !i

      open(15,file='depth_diff.gr3')
      write(15,*)'Depth=new-old'
      write(15,*)ne,np
      do i=1,np
        write(15,'(i8,2(1x,e20.12),1x,f9.3)')i,x(i),y(i),dp(i)-dp00(i)
      enddo !i
      do i=1,ne
        write(15,*)i,i34(i),nm(i,1:i34(i))
      enddo !i

      stop
      end
