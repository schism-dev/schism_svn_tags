   MODULE WWM_MPI
   implicit none
#if defined SELFE || defined MPI
   include 'mpif.h'
#endif

#ifdef SELFE
   use elfe_msgp
#endif

    CONTAINS
    subroutine wwm_abort(string,error)
    implicit none
    character(*),optional,intent(in) :: string !string to print
    integer,optional,intent(in) :: error       !mpi errorcode
    integer :: ierror,i
#ifdef SELFE
    call parallel_abort(string, error)
#elif defined MPI
!    call mpi_abort(comm,0,ierror)
#else
    write(*,*) string
    STOP
#endif
   end subroutine wwm_abort

   END MODULE WWM_MPI
