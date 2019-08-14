PROG = $(HOME)/bin/stat3

OBJS = stat_pool.f90 stat3.f90 \
#       dislin.o diclass.o \

#LIBS = /opt/dislin/lib/dislin-9.5.so
#LIBS = /opt/dislin/libdislin.so

F90 = ifort
#F90 = $(FLINKER)
#F90 = gfortran
#F90 = pgf90
#F90 = pathf95
#F90 = f77

#INTEL
#save opti
F90FLAGS = -O1 -axSSE4.2 -fpp -DWWMONLY -DNCDF -DST4 -traceback -g #-check all
#extreme aggressive opti ...
#F90FLAGS = -O5 -fp-model fast -opt-multi-version-aggressive -axSSE4.2 -no-heap-arrays -fpp -DST4 -DWWMONLY -DNCDF -stack_temps -unroll-aggressive -vec-guard-write #-check all -warn all
#aggressive opti
#F90FLAGS = -O1 -fp-model fast -opt-multi-version-aggressive -axSSE4.2 -no-heap-arrays -fpp -DWWMONLY -DNCDF -DST4 -stack_temps -unroll-aggressive -vec-guard-write -traceback -g -openmp
#aggressive opti no openmp
#F90FLAGS = -fast -axSSE4.2 -no-heap-arrays -fpp -DWWMONLY -stack_temps -unroll-aggressive -vec-guard-write
#debug 1
#F90FLAGS = -O5 -mp1 -fpp -DWWMONLY -traceback -g -axSSE4.2 -openmp -g -traceback -check all
#debug 2
#F90FLAGS =  -fpp -DWWMONLY -DST4 -DNCDF -g -traceback -check all
#aggresive debug openmp
#F90FLAGS = -fpp -DWWMONLY -DNCDF -DST4  -g -traceback -check all -openmp
#debug no openmp
#F90FLAGS = -g -traceback -check all 
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

F90FLAGS += ${PETSC_INCLUDE}

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
