!     Basic routines for reading SCHISM netcdf outputs
!     Author: Joseph Zhang
!     Date: Sept 2017

!     Routines in this module
!     function signa
!     subroutine readheader

    module extract_mod
    use netcdf
    implicit none
    public

    character(len=48), save :: start_time 
    integer, save :: itime_id,ielev_id,nrec,nvrt,kz,np,ne,ns,ncid2
    real, save :: h0,fill_in,dtout
    real, save, allocatable :: x(:),y(:),dp(:)
    integer, save, allocatable :: kbp00(:),i34(:),elnode(:,:)
    real*8,allocatable :: timeout2(:)

      contains
!================================================================
!================================================================
      function signa(x1,x2,x3,y1,y2,y3)
!...  Compute signed area formed by pts 1,2,3
      real :: signa,x1,x2,x3,y1,y2,y3

      signa=((x1-x3)*(y2-y3)-(x2-x3)*(y1-y3))/2

      end function signa

!================================================================
!     Read in static info from nc
!     Returned vars: ne,np,ns,nrec,start_time,[x y dp](np),elnode,i34,nvrt,
!                    itime_id,ielev_id,h0,dtout
!================================================================
      subroutine readheader(fname)
      character(len=*),intent(in) :: fname !schout_*.nc
      integer :: i,varid1,varid2,dimids(3),istat,nvtx,iret

      iret=nf90_open(trim(adjustl(fname)),OR(NF90_NETCDF4,NF90_NOWRITE),ncid2)
      iret=nf90_inq_dimid(ncid2,'nSCHISM_hgrid_edge',i)
      iret=nf90_Inquire_Dimension(ncid2,i,len=ns)
      iret=nf90_inq_dimid(ncid2,'nSCHISM_vgrid_layers',i)
      iret=nf90_Inquire_Dimension(ncid2,i,len=nvrt)
      iret=nf90_inq_varid(ncid2,'minimum_depth',varid1)
      iret=nf90_get_var(ncid2,varid1,h0)
      iret=nf90_inq_varid(ncid2,'SCHISM_hgrid_face_nodes',varid1)
      iret=nf90_Inquire_Variable(ncid2,varid1,dimids=dimids(1:2))
      iret=nf90_Inquire_Dimension(ncid2,dimids(1),len=nvtx)
      iret=nf90_Inquire_Dimension(ncid2,dimids(2),len=ne)
      if(nvtx/=4) stop 'readheader: vtx/=4'
      iret=nf90_inq_varid(ncid2,'SCHISM_hgrid_node_x',varid2)
      iret=nf90_Inquire_Variable(ncid2,varid2,dimids=dimids)
      iret=nf90_Inquire_Dimension(ncid2,dimids(1),len=np)
      iret=nf90_inq_varid(ncid2,'time',itime_id)
      iret=nf90_Inquire_Variable(ncid2,itime_id,dimids=dimids)
      iret=nf90_Inquire_Dimension(ncid2,dimids(1),len=nrec)
      iret=nf90_get_att(ncid2,itime_id,'base_date',start_time)
      iret=nf90_inq_varid(ncid2,'elev',ielev_id)
      if(iret.ne.NF90_NOERR) then
        print*, nf90_strerror(iret); stop 'readheader: error reading header'
      endif

      if(allocated(x)) deallocate(x)
      if(allocated(y)) deallocate(y)
      if(allocated(dp)) deallocate(dp)
      if(allocated(kbp00)) deallocate(kbp00)
      if(allocated(elnode)) deallocate(elnode)
      if(allocated(i34)) deallocate(i34)
      if(allocated(timeout2)) deallocate(timeout2)
      allocate(x(np),y(np),dp(np),kbp00(np),i34(ne),elnode(4,ne),timeout2(nrec),stat=istat)
      if(istat/=0) stop 'readheader: failed to allocate (3)'
      iret=nf90_get_var(ncid2,varid1,elnode)
      iret=nf90_get_var(ncid2,varid2,x)
      iret=nf90_get_var(ncid2,itime_id,timeout2,(/1/),(/nrec/))
      dtout=timeout2(2)-timeout2(1)

      iret=nf90_inq_varid(ncid2,'SCHISM_hgrid_node_y',varid1)
      iret=nf90_get_var(ncid2,varid1,y)
      iret=nf90_inq_varid(ncid2,'depth',varid1)
      iret=nf90_get_var(ncid2,varid1,dp)
      !iret=nf90_inq_varid(ncid2,'node_bottom_index',varid1)
      !iret=nf90_get_var(ncid2,varid1,kbp00)

      iret=nf90_close(ncid2)

      !print*, 'nc dim:',nvrt,np,ne,nrec,start_time

      !Calc i34
      i34=4 !init
      do i=1,ne
        if(elnode(4,i)<0) i34(i)=3
      enddo !i

      end subroutine readheader
!================================================================
!================================================================
    end module extract_mod
