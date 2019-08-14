include ${PETSC_DIR}/conf/variables
PROG = $(HOME)/bin/wwmadv
#PROG = wwmadv

FFLAGS = -DST42 

# Put MPI = on for compiling with mpi, nothing if no MPI is wished.
MPI = 
# Put PETSC = on for compiling with petsc, nothing if no petsc is wished.
PETSC = 
# Put NETCDF = on for compiling with netcdf, nothing if no netcdf is wished.
NETCDF = on
# Put OPENMP = on for compiling with openmp, nothing if no openmp is wished.
OPENMP = 
# Put OPENMP = on for compiling with openmp, nothing if no openmp is wished.
DARKO =
# Put OPENMP = on for compiling with openmp, nothing if no openmp is wished.
SNL4_TSA = 

#F90 = gfortran
F90 = mpif90
#F90 = mpif90-vt
#F90 = ifort
#F90 = $(FLINKER)
#F90 = pgf90
#F90 = pathf95
#F90 = f77



#INTEL
#save opti
#F90OPTI = -fpp -traceback -g
#F90OPTI = -O1 -axSSE4.2 -traceback -g  -check uninit -check bounds -check pointers -warn interfaces,nouncalled -gen-interface
#F90OPTI = -fpp -O2 -stand f03 -assume realloc_lhs -fstack-protector -assume protect_parens
#F90OPTI = -warn interfaces,nouncalled -fpp -gen-interface -g -traceback -check uninit -check bounds -check pointers
#F90OPTI = -warn interfaces,nouncalled -fpp -gen-interface -g -traceback -check all
#F90OPTI = -fpp -O2 -stand f03 -assume realloc_lhs -check all -traceback -fstack-protector -assume protect_parens -implicitnone
#extreme aggressive opti ...
#F90OPTI = -g -traceback -O1 -axSSE4.2 -unroll-aggressive -assume byterecl #-check all -warn all
#F90OPTI = -g -traceback -O1 -axSSE4.2 -unroll-aggressive -assume byterecl 
F90OPTI = -g -traceback -O5 -axSSE4.2 -unroll-aggressive -assume byterecl 
#aggressive opti
#F90OPTI = -O1 -fp-model fast -opt-multi-version-aggressive -axSSE4.2 -stack_temps -unroll-aggressive -vec-guard-write -traceback -g
#aggressive opti no openmp
#F90OPTI = -O1 -axSSE4.2 -fpp
#debug 1
#F90OPTI = -O1 -mp1 -fpp -traceback -g -traceback -check pointers -check bound -check uninit 
#debug 2
#F90OPTI =  -fpp -g -traceback -check all
#aggresive debug openmp
#F90OPTI = -fpp -g -traceback -check all -openmp
#debug no openmp
#F90OPTI = -fpp -g -traceback -fp-model precise -check all 
#PGI
#F90OPTI = -g -O1 -Mprof=dwarf -Mprof=time -Minfo=ccff -mp
#F90OPTI = -fastsse -Munroll -Minline=reshape -Minline=levels:10 -Mvect -Mipa=fast Mipa=levels:10 -Mlarge_arrays -mp
#F90OPTI  = -g -Mbounds -Mchkfpstk -Mchkptr -Mchkstk
#GFORTRAN
#F90OPTI = -g -O3 -ffree-line-length-0 -fbacktrace -fopenmp -Wall

#LDFLAGS = -Mprof=dwarf -Mprof=lines -Mprof=func -Mprof=time Minfo=ccff
#LDFLAGS = -lpthread

#ITT_BASE = /home/aron/intel/itt

#LIBS = #libgoto_penrynp-r1.26.a #sparskit_pgi.a libgoto.a #sparskit.a

#ARCH  = $(shell uname -m | sed -e s/x86_64/32e/ -e s/i.86/32/ -e s/ia// )
#LPATH = $(ITT_BASE)/tcheck/bin/$(ARCH)



OBJSELFE = elfe_glbl.o elfe_msgp.o solver_subs.o io_subs.o grid_subs.o misc_subs.o
OBJWWM = wwm_datapl.o wwm_petscpool.o  wwm_m_constants.o wwm_m_fileio.o wwm_serv_xnl4v5.o wwm_mod_xnl4v5.o wwm_petsc_seriell.o wwm_petsc_parallel.o \
       wwm_petsc_block.o wwm_petsc_controller.o wwm_aux.o wwm_mjdv2.o wwm_blasaux.o wwm_sparskit.o \
       wwm_ardhuin_old.o wwm_wind.o wwm_ardhuin_new.o wwm_dislin.o wwm_diclass.o \
       wwm_gridcf.o wwm_hotfile.o wwm_parall_solver.o wwm_initio.o wwm_input.o wwm_bdcons.o \
       wwm_specparam.o wwm_windinput.o wwm_dissipation.o wwm_snl3.o wwm_snl4.o wwm_snl4_wrt.o wwm_snl42.o wwm_babanin.o wwm_sourceterms.o \
       wwm_specint.o wwm_nums1d.o wwm_numtheta.o wwm_numsigma.o wwm_netcdf.o wwm_fluctsplit.o \
       wwm_ecmwf.o wwm_output.o wwm_compute.o wwm_diffrac.o wwm_coupl_roms.o wwm_coupl_shyfem.o wwm_coupl_selfe.o wwm_coupl_timor.o wwm_main.o
#       dislin.o diclass.o \

#LIBS = -L/home/aron/opt/dislin/lib/dislin.a

#LIBS = /home/aron/opt/dislin/lib/dislin-10.2.so

ifdef NETCDF
  all: $(NETC) $(PROG)
  LIBS += -L$(NETCDF_LIBDIR) -lnetcdf
  FFLAGS +=  -DNCDF -I$(NETCDF_INCDIR)
else
  all: $(PROG)
endif

ifdef MPI
  MPI_INCLUDE = -I/opt/mpich2_1.4_intel/include/
  LIBS += $(METIS_PATH)/libparmetis.a $(METIS_PATH)/libmetis.a
  MPIFLAG  = -DWWM_MPI -D'MPIVERSION=2' -DUSE_WWM -DWWM_SOLVER -DDEBUG
  OBJS= $(OBJSELFE) $(OBJWWM)
else
  MPIFLAG = 
  OBJS= $(OBJWWM)
endif

ifdef OPENMP 
  OPENMPFLAG  = -openmp 
endif

ifdef PETSC
  MPIFLAG += -DPETSC ${PETSC_FC_INCLUDES} 
  LIBS  += ${PETSC_LIB}  
endif

ifdef DARKO 
  DARKOFLAG = -DDARKO 
endif

ifdef SNL4_TSA 
  SNL4_TSAFLAG = -DSNL4_TSA
endif

MPI_INCLUDE = -I/opt/mpich2_1.4_intel/include/  


F90FLAGS = ${MPIFLAG} ${FFLAGS} ${F90OPTI} ${OPENMPFLAG} ${DARKOFLAG} ${SNL4_TSAFLAG} ${MPI_INCLUDE}

NETC = NETC

$(NETC):
	cp -f $(NETCDF_INCDIR)/netcdf.mod .
	cp -f $(NETCDF_INCDIR)/typesizes.mod .

$(PROG): $(OBJS)
	$(F90) -o $(PROG) $(F90FLAGS) $(OBJS) $(LIBS) 

clean:
	rm -f *.o *.oo *.obj *.ipo *.mod *.map *__genmod.f90 *.ilk *.pdb $(PROG)

cleanall:
	rm -f *.*~ *.spg *.smb *.o *.oo *.obj *.ipo *.mod *.map *.ilk *.pdb *genmod* $(PROG)

.SUFFIXES: $(SUFFIXES) .F90 .f .ftn

.F90.o:
	$(F90) $(F90FLAGS) $(TCFLAG)   -c $<

.ftn.o:
	$(F90) $(F90FLAGS) -c $<

