include ${PETSC_DIR}/conf/variables
PROG = $(HOME)/bin/wwmadv

OBJS = wwm_datapl.o wwm_petscpool.o wwm_aux.o wwm_mjdv2.o wwm_blasaux.o wwm_sparskit.o \
       wwm_ardhuin_old.o wwm_ardhuin_new.o wwm_dislin.o wwm_diclass.o \
       wwm_gridcf.o wwm_initio.o wwm_input.o wwm_bdcons.o \
       wwm_specparam.o wwm_windinput.o wwm_dissipation.o wwm_snl3.o wwm_snl4.o wwm_snl42.o wwm_babanin.o wwm_sourceterms.o \
       wwm_specint.o wwm_nums1d.o wwm_numtheta.o wwm_numsigma.o wwm_fluctsplit.o \
       wwm_output.o wwm_compute.o wwm_radiation.o wwm_diffrac.o wwm_main.o 
#       dislin.o diclass.o \

#LIBS = -L$(NETCDF_LIBDIR) -lnetcdf -L/home/aron/opt/dislin/lib/dislin.a

LIBS = /home/aron/opt/dislin/lib/dislin-10.2.so  -L$(NETCDF_LIBDIR) -lnetcdf ${PETSC_DIR}/lib/libpetsc.a

F90 = ifort
#F90 = $(FLINKER)
#F90 = gfortran
#F90 = pgf90
#F90 = pathf95
#F90 = f77

#INTEL
#save opti
#F90FLAGS = -fpp -DWWMONLY -DNCDF -DST42 -traceback -g  -check pointers -check uninit -check bounds
#extreme aggressive opti ...
#F90FLAGS = -g -traceback -fpp -DST42 -DWWMONLY -DNCDF -DPETSC -O5 -axSSE4.2 -unroll-aggressive -vec-guard-write -openmp #-check all -warn all
#aggressive opti
F90FLAGS = -O1 -axSSE4.2 -fpp -DWWMONLY -DNCDF -DST42 -DPETSC -traceback -g -openmp
#aggressive opti no openmp
#F90FLAGS = -O1 -axSSE4.2 -fpp -DWWMONLY -DST42 -DNCDF
#debug 1
#F90FLAGS = -O1 -mp1 -fpp -DWWMONLY -DNCDF -DPETSC -DST42 -traceback -g -traceback -check pointers -check bound -check uninit 
#debug 2
#F90FLAGS =  -fpp -DWWMONLY -DST4 -DNCDF -g -traceback -check all
#aggresive debug openmp
#F90FLAGS = -fpp -DWWMONLY -DNCDF -DST4  -g -traceback -check all -openmp
#debug no openmp
#F90FLAGS = -fpp -DWWMONLY -DNCDF -g -traceback -fp-model precise -check all 
#PGI
#F90FLAGS = -g -O1 -Mprof=dwarf -Mprof=time -Minfo=ccff -mp -DWWMONLY
#F90FLAGS = -fastsse -DWWMONLY  -Munroll -Minline=reshape -Minline=levels:10 -Mvect -Mipa=fast Mipa=levels:10 -Mlarge_arrays -mp
#F90FLAGS  = -g -Mbounds -Mchkfpstk -Mchkptr -Mchkstk -DWWMONLY
#GFORTRAN
#F90FLAGS = -g -O3 -DWWMONLY -ffree-line-length-132 -fbacktrace -fopenmp

#LDFLAGS = -Mprof=dwarf -Mprof=lines -Mprof=func -Mprof=time Minfo=ccff
#TCFLAG  = -tcheck
#LDFLAGS = -lpthread

#ITT_BASE = /home/aron/intel/itt

#LIBS = #libgoto_penrynp-r1.26.a #sparskit_pgi.a libgoto.a #sparskit.a

#ARCH  = $(shell uname -m | sed -e s/x86_64/32e/ -e s/i.86/32/ -e s/ia// )
#LPATH = $(ITT_BASE)/tcheck/bin/$(ARCH)

PETSC_INCLUDE = -I${PETSC_DIR}/include 
MPI_INCLUDE = -I/home/aron/opt/mpich2_1.3_intel/include/  

F90FLAGS += ${PETSC_INCLUDE}
F90FLAGS += ${MPI_INCLUDE}

NETC = NETC

all: $(NETC) $(PROG)

$(NETC):
#	cp -f $(NETCDF_INCDIR)/netcdf.mod .
#	cp -f $(NETCDF_INCDIR)/typesizes.mod .

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) -o $(PROG) $(F90FLAGS) $(TCFLAG) $(OBJS) $(LIBS) $(LDFLAGS) ${PETSC_LIB}

clean:
	rm *.o *.oo *.obj *.ipo *.mod *.map *.ilk *.pdb $(PROG)

cleanall:
	rm  *.*~ *.spg *.smb *.o *.oo *.obj *.ipo *.mod *.map *.ilk *.pdb $(PROG)

.SUFFIXES: $(SUFFIXES) .F90 .f .ftn

.F90.o:
	$(F90) $(F90FLAGS) $(TCFLAG)   -c $<

.ftn.o:
	$(F90) $(F90FLAGS) -c $<
