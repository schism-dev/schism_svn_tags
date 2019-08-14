!     Convert fluxflag.gr3 or tvd.gr3 to .prop that can be displayed in gredit
!     Inputs: <file name>
!     Output: .prop 

!     ifort -Bstatic -o gr3_to_prop gr3_to_prop.f90
!     pgf90 -O2 -mcmodel=medium  -Bstatic -o gr3_to_prop gr3_to_prop.f90

      allocatable :: xnd(:),ynd(:),dp(:)
      integer,allocatable :: elnode(:,:)
      character(len=100) :: fname

      print*, 'Input file name to convert:'
      read*, fname
      print*, 'Take the min (0; tvd.gr3) or max(1; fluxflag.gr3) of the 3 nodes? '
      read*, imm

      fname=adjustl(fname)
      lfname=len_trim(fname)

      open(14,file=fname,status='old')
      open(13,file=fname(1:lfname-4)//'.prop',status='replace')
      read(14,*); read(14,*)ne,np
      allocate(xnd(np),ynd(np),dp(np),elnode(4,ne))
      do i=1,np
        read(14,*)j,xnd(i),ynd(i),dp(i)
      enddo !i
      do i=1,ne
        read(14,*)j,k,elnode(1:k,i)
      enddo !i
      close(14)

      do i=1,ne
        if(imm==0) then !min
          dpe=minval(dp(elnode(:,i)))
        else
          dpe=maxval(dp(elnode(:,i)))
        endif
        write(13,*)i,dpe
      enddo !i

      stop
      end
