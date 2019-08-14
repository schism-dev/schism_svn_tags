!     Extract time series at all nodes along a horizontal slab for 2D or 3D variables
!     Authors: Charles Seaton, Joseph Zhang, Paul Turner
!     Date: April 2011

!     Fatal errors may result in some abnormal cases (e.g.,
!     cannot find a containing triangle for a build point etc) - search the source code
!     for more details.

      program extract_slab
      use extract_mod
      
      character(len=500) :: basedir,fname
      character(len=50) :: binary_file

!     Fill in value used for abnormal case (e.g., below bottom, dry etc)
      fill_in=-9999.

!     basedir - directory where binary and build point files are stored
      basedir='./RUN05g/outputs/'

!     binary_file - 'elev.61', 'hvel.64' etc; it should be inside basedir
      binary_file='salt.63'

!     istack[1,2] - start and end stack # for binary outputs;
      istack1=2; istack2=2

!     Flag to indicate whether extracting along a S level (0: no; 1 yes)
      ialong_S=0
 
!     If along a S level, specify level index (no effects on 2D variables)
      if(ialong_S==1) then
        klev0=1
      endif

!     If not along a S level, specify z value and a flag to indicate whether input z value
!     is relative to MSL (ifs=0) or free surface (ifs=1). These flags have 
!     no effects on 2D variables
      if(ialong_S==0) then
        ifs=0
        zout=-5.
        if(ifs==1.and.zout<0) stop 'ifs==1.and.zout<0'
      endif

!     Read binary header
      write(it_char,'(i12)')istack1
      it_char=adjustl(it_char)
      leng=len_trim(it_char)
      fname = trim(basedir)//'/'//it_char(1:leng)//'_'//binary_file
      !print *, 'Binary name: ',fname
      call readheader(fname)

!     Allocate outputs arrays
      allocate(varout3(ivs,np,nrec),outtime(nrec),stat=istat)
      if(istat/=0) stop 'Failed to allocate (8)'

!...  Time iteration
      do istack=istack1,istack2
!       Name of the appropriate binary file
        write(it_char,'(i12)')istack
        it_char=adjustl(it_char)
        leng=len_trim(it_char)
        fname = trim(basedir)//it_char(1:leng)//'_'//binary_file
        call readslab(fname)
!       At this point, varout3(1:ivs,1:np,1:nrec) is the output with
!       times given by outtime(1:nrec) (in sec), where nrec is # of time steps within 
!       each stack, np is the # of nodes, ivs=1 indicates
!       scalar (1) output; ivs=2 indicates vector outputs (e.g. 1=u; 2=v). 

        if(istack==istack2) then
          do i=1,np
            write(18,*)i,x(i),y(i),varout3(1:ivs,i,nrec)
          enddo !i
        endif
      enddo !istack

      end program extract_slab

