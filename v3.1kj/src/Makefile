################################################################################
# Parallel SELFE Makefile
#
# User makes environment settings for particular OS / PLATFORM / COMPILER / MPI
# below. These settings are based on the following options.
#
# Compiler name:
#   FCS: Serial compiler (for utilities)
#   FCP: Parallel compiler
#   FLD: Linker (in general same as parallel compiler)
#
# Compilation flags
#   FCSFLAGS: Flags for serial compilation
#   FCPFLAGS: Flags for parallel compilation (including all pre-processing flags)
#   FLDFLAGS: Flags for linker (e.g., -O2)
#
# Preprocessor flags:
#   DEBUG: Enable debugging code
#   ORDERED_SUM: Enable globally ordered sums & dot-products for bit reproducibility
#     of state quantities independent of number of processors (note: this can
#     significantly degrade performance);
#   INCLUDE_TIMING: Enable wallclock timing of code (note: this can have slight
#     effect on performance);
#   MPI_VERSION = 1 or 2: Version of MPI (try 2 first, if compile fails due to mpi
#     related errors then switch to version 1;
#
# Libraries (needed for parallel code)
#   MTSLIBS: Flags for linking ParMeTiS/MeTiS libaries
#   ALTLIBS: Flags for linking alternate solver libraries (LAPACK or ITPACK,
#            these are just for testing)
################################################################################


################################################################################
# Comment out the following line if you don't want the build dependent
# upon the Makefile, otherwise use the name of the Makefile
################################################################################
#MAKEFILE = Makefile

#ENV         = PLEIADES
#ENV         = LONI
#ENV         = RANGER
#ENV         = VIMS
#ENV         = XEON
#ENV         = GCMPICH2 #GFortran
#ENV         = CANOPUS
ENV         = SIRIUS

################################################################################
# Executable
################################################################################
EXEC   = pelfe_3.1kj_sirius_WWM

################################################################################
# Define pre-processing flags.
# Comment out unwanted flags.
################################################################################
# -DSELFE always on
  PPFLAGS = -DSELFE 

# Wave models need a coupler
# SWAN (not active)
#  USE_SWAN = yes
#  PPFLAGS := $(PPFLAGS) -DUSE_SWAN
  #MPI version
#  SWAN_DIR = /disk/ambcs03/0/yinglong/SWAN/swan4072_MPI/

# Precip/evaporation model
#  PPFLAGS := $(PPFLAGS) -DPREC_EVAP 

# MM5 in heat exchange model
#  PPFLAGS := $(PPFLAGS) -DMM5

# GOTM turbulence closure model
#  USE_GOTM = yes
#  PPFLAGS := $(PPFLAGS) -DUSE_GOTM 

# Comment out if timer report is not needed.
#  PPFLAGS := $(PPFLAGS) -DINCLUDE_TIMING

# Debug mode (more time consuming)
#  PPFLAGS := $(PPFLAGS) -DDEBUG 

# Wind wave model WWM
  USE_WWM = yes
  PPFLAGS := $(PPFLAGS) -DUSE_WWM -DNCDF -DST4

# Harmonic analysis tool
#  USE_HA = yes
#  PPFLAGS := $(PPFLAGS) -DUSE_HA

# Select only _one_ model from below

# Ecological model - NAPZD Spitz (internal use only)
#  USE_NAPZD = yes
#  PPFLAGS := $(PPFLAGS) -DUSE_NAPZD

# Or:
# Ecological model (EcoSim)
#  USE_ECO = yes
#  PPFLAGS := $(PPFLAGS) -DUSE_ECO

# Or:
# CE-QUAL-ICM
# USE_ICM = yes
# PPFLAGS := $(PPFLAGS) -DUSE_ICM

# Or:
# Sediment model 
#  USE_SED = yes
#  PPFLAGS := $(PPFLAGS) -DUSE_SED

  ##Bedload 
#  PPFLAGS := $(PPFLAGS) -DBEDLOAD

  ##Bedload - MPM model
#  PPFLAGS := $(PPFLAGS) -DBEDLOAD_MPM

  ##slope formulation
#  PPFLAGS := $(PPFLAGS) -DDAMGAARD
#  PPFLAGS := $(PPFLAGS) -DDELFT
#  PPFLAGS := $(PPFLAGS) -DCARMO

  ##Bedload - VR model
#  PPFLAGS := $(PPFLAGS) -DBEDLOAD_VR

  ##Suspended load
#  PPFLAGS := $(PPFLAGS) -DSUSPLOAD

  ##boundary conditions for WENO
  ## default strictly monotonic
#  PPFLAGS:= $(PPFLAGS) -DLINEAR_CONTINUATION
#  PPFLAGS:= $(PPFLAGS) -DNEUMANN

  ## Morphology
#  PPFLAGS := $(PPFLAGS) -DSED_MORPH

  ## Choose one drag formulation from the following 3 choices (only 1st one is functional now)
#  PPFLAGS := $(PPFLAGS) -DUV_LOGDRAG
  #PPFLAGS := $(PPFLAGS) -DUV_QDRAG
  #PPFLAGS := $(PPFLAGS) -DUV_LDRAG

  ##sediment density in eqstate
#   PPFLAGS:= $(PPFLAGS) -DDENSED

# Or:
# Oil spill model (not active)
#  USE_OIL = yes
#  PPFLAGS := $(PPFLAGS) -DUSE_OIL

# Add -DNO_TR_15581 like below for allocatable array problem in sflux_subs.F90
# PPFLAGS := $(PPFLAGS) -DNO_TR_15581

# For openMPI compiler, search for "USE_OPEN64" below for compiler flags

# Obsolete flags: use USE_WRAP flag to avoid problems in ParMetis lib (calling C from FORTRAN)
# PPFLAGS := $(PPFLAGS) -DUSE_WRAP 

################################################################################
# Environment for NASA/PLEIADES (SGI ALTIX)
################################################################################
ifeq ($(ENV),PLEIADES)
  FCP = ifort
  FLD = $(FCP)
  PPFLAGS := $(PPFLAGS) -DMPIVERSION=1
  FCPFLAGS = $(PPFLAGS) -O2 -mcmodel=medium  -Bstatic -assume byterecl
  FLDFLAGS = -O2 -mcmodel=medium -shared-intel -lmpi #for final linking of object files
  #Libraries
  MTSLIBS = -L/home1/yjzhang/ParMetis-3.1-Sep2010/ -lparmetis -lmetis
  CDFLIBS = -L/nasa/netcdf/3.6.0/intel/lib -lnetcdf
  CDFMOD = -I/nasa/netcdf/3.6.0/intel/include/ # modules for netcdf
  ifdef USE_GOTM
    GTMMOD =  -I/home3/yjzhang/GOTM_Intel64/modules/IFORT/ #modules
    GTMLIBS = -L/home3/yjzhang/GOTM_Intel64/lib/IFORT/ -lturbulence_prod -lutil_prod
  else
    GTMMOD =
    GTMLIBS =
  endif
endif

################################################################################
# Environment for TACC/Ranger (need to remove USE_WRAP)
################################################################################
ifeq ($(ENV),RANGER)
#  USE_WRAP = yes
  FCP = mpif90 #MPI compiler
  FLD = $(FCP)
# MPI vserion (1 or 2)
  PPFLAGS := $(PPFLAGS) -DMPIVERSION=1 #-DUSE_WRAP
  FCPFLAGS = $(PPFLAGS) -O2 -mcmodel=medium  -Bstatic #MPI code
  FLDFLAGS = -O2 -mcmodel=medium #for final linking of object files
  #####Libraries
  MTSLIBS = -L/share/home/01621/zhangy/Joseph/ParMetis-3.1-Sep2010/ -lparmetis -lmetis
  CDFLIBS = -L/opt/apps/pgi7_2/netcdf/3.6.2/lib/ -lnetcdf
  CDFMOD = -I/opt/apps/pgi7_2/netcdf/3.6.2/include/ # modules for netcdf
  ifdef USE_GOTM
    GTMMOD =  -I/work/01555/harryw/gotm-4.0.0/modules/PGF90/ #modules
    GTMLIBS = -L/work/01555/harryw/gotm-4.0.0/lib/PGF90/ -lturbulence_prod -lutil_prod
  else
    GTMMOD =
    GTMLIBS =
  endif
endif

################################################################################
# Environment for LONI / 64 bit /  Intel Compiler / MPICH2 (LONI)
################################################################################
ifeq ($(ENV),LONI)
  FCP = mpif90 #MPI compiler
  FLD = $(FCP)
# MPI vserion (1 or 2)
  PPFLAGS := $(PPFLAGS) -DMPIVERSION=2
  FCPFLAGS = $(PPFLAGS) -O2 -Bstatic -assume byterecl  #MPI code
  FLDFLAGS = -O2   #for final linking of object files
  #####Libraries
  MTSLIBS = -L/home/hvwang/ParMetis-3.1-64bit -lparmetis -lmetis
  CDFLIBS = -L/usr/local/packages/netcdf/4.0/intel-11.1/lib -lnetcdf
  CDFMOD = -I/usr/local/packages/netcdf/4.0/intel-11.1/include # modules for netcdf
  ifdef USE_GOTM
    GTMMOD =  -I/home/users/yinglong/GOTM/gotm-3.2.5/32bit/gotm-3.2.5/modules/IFORT/ #modules
    GTMLIBS = -L/home/users/yinglong/GOTM/gotm-3.2.5/32bit/gotm-3.2.5/lib/IFORT/ -lturbulence_prod  -lutil_prod
  else
    GTMMOD =
    GTMLIBS =
  endif
endif

################################################################################
# Environment for TYPHOON cluster of VIMS-College of William and Mary
################################################################################
ifeq ($(ENV),VIMS)
  FCP = mpif90 #MPI compiler
  FLD = $(FCP)
# MPI vserion (1 or 2)
  PPFLAGS := $(PPFLAGS) -DMPIVERSION=2
  FCPFLAGS = $(PPFLAGS) -O2 -C -Bstatic -assume byterecl  #MPI code
  FLDFLAGS = -O2 -C #for final linking of object files
  #####Libraries
  MTSLIBS = -L/home/hvwang/ParMetis-3.1-64bit -lparmetis -lmetis
  CDFLIBS = -L/usr/local/packages/netcdf/4.0/intel-11.1/lib -lnetcdf
  CDFMOD = -I/usr/local/packages/netcdf/4.0/intel-11.1/include # modules for netcdf
  ifdef USE_GOTM
    GTMMOD =
    GTMLIBS =
  else
    GTMMOD =
    GTMLIBS =
  endif
endif

################################################################################
# Environment for Linux / 32 bit /  Intel Compiler / MPICH
################################################################################
ifeq ($(ENV),XEON)
  FCP = mpif90 #MPI compiler
  FLD = $(FCP)
# MPI vserion (1 or 2) 
  PPFLAGS := $(PPFLAGS) -DMPIVERSION=1 
  FCPFLAGS = $(PPFLAGS) -O2 -Bstatic -assume byterecl #MPI code
  FLDFLAGS = -O2 #for final linking of object files
  #####Libraries
  MTSLIBS = -L/home/users/yinglong/SELFE/ParMetis-3.1-32bit -lparmetis -lmetis
  CDFLIBS = -L/usr/local/netcdf/lib -lnetcdf 
  CDFMOD = -I/usr/local/netcdf/include # modules for netcdf
  ifdef USE_GOTM
    GTMMOD =  -I/home/users/yinglong/GOTM/gotm-3.2.5/32bit/gotm-3.2.5/modules/IFORT/ #modules
    GTMLIBS = -L/home/users/yinglong/GOTM/gotm-3.2.5/32bit/gotm-3.2.5/lib/IFORT/ -lturbulence_prod  -lutil_prod 
  else
    GTMMOD =
    GTMLIBS =
  endif
endif

################################################################################
# Environment for GCMPICH2 / 64 bit /  GNU Compiler / MPICH2 (GCMPICH2)
# See gfortran_notes for other details
################################################################################
ifeq ($(ENV),GCMPICH2)
  FCP = mpif90 -f90=/usr/bin/gfortran #MPI compiler
  FLD = $(FCP)
# MPI vserion (1 or 2)
  PPFLAGS := $(PPFLAGS) -DMPIVERSION=2
  FCPFLAGS = $(PPFLAGS) -O2 -Bstatic # -assume byterecl  #MPI code
  FLDFLAGS = -O2   #for final linking of object files
  #####Libraries
  MTSLIBS = -L/home/kyle/selfe/modified/V3.1d/ParMetis-3.1-Sep2010.gnucmp -lparmetis -lmetis
  CDFLIBS = -L/home/kyle/selfe/modified/V3.1d/netcdf/lib -lnetcdf
  CDFMOD = -I/home/kyle/selfe/modified/V3.1d/netcdf/include # modules for netcdf
  ifdef USE_GOTM
    GTMMOD =
    GTMLIBS =
  else
    GTMMOD =
    GTMLIBS =
  endif
endif

################################################################################
# Environment for Linux / 64 bit /  Intel Compiler / MPICH (canopus)
################################################################################
ifeq ($(ENV),CANOPUS)
  #FCP = /home/yinglong/mpich2/bin/mpif90
  FCP = mpif90
  FLD = $(FCP)
# MPI vserion (1 or 2) 
  PPFLAGS := $(PPFLAGS) -DMPIVERSION=2 
# USE_OPEN64 = yes
  ifdef USE_OPEN64
    FCPFLAGS = $(PPFLAGS) -DUSE_OPEN64 -O2 -OPT:Olimit=0 #MPI code
  else
    FCPFLAGS = $(PPFLAGS) -O2 -Bstatic -assume byterecl #MPI code
  endif
  FLDFLAGS = -O2 #for final linking of object files
  #####Libraries
  MTSLIBS = -L/home/users/yinglong/SELFE/ParMetis-3.1-64bit/Intel -lparmetis -lmetis
  #MTSLIBS = -L/home/users/yinglong/SELFE/ParMetis-3.1-Sep2010/ -lparmetis -lmetis
  CDFLIBS = -L/usr/local/lib -lnetcdf 
  CDFMOD = -I/usr/local/include # modules for netcdf
  ifdef USE_GOTM
    GTMMOD =  -I/home/users/yinglong/GOTM/gotm-3.2.5/Intel64/modules/IFORT/ #modules
    GTMLIBS = -L/home/users/yinglong/GOTM/gotm-3.2.5/Intel64/lib/IFORT/ -lturbulence_prod  -lutil_prod 
  else
    GTMMOD =
    GTMLIBS =
  endif
endif

################################################################################
# Environment for SIRIUS cluster
################################################################################
ifeq ($(ENV),SIRIUS)
#  USE_WRAP = yes
  FCP = mpif90
  FLD = $(FCP)
# MPI vserion (1 or 2) 
  PPFLAGS := $(PPFLAGS) -DMPIVERSION=2 #-DUSE_WRAP
  #-CB is much slower to compile
  #FCPFLAGS = $(PPFLAGS) -O2 -CB -Bstatic -assume byterecl #MPI code; check bound
  FCPFLAGS = $(PPFLAGS) -O2 -Bstatic -assume byterecl #MPI code
  FLDFLAGS = -O2 #for final linking of object files
  #####Libraries
#  MTSLIBS = -L/share/apps/ParMetis/ -lparmetis -lmetis
  MTSLIBS = -L/home/users/yinglong/SELFE/ParMetis-3.1-Sep2010/ -lparmetis -lmetis
  CDFLIBS = -L/share/apps/netcdf/lib/ -lnetcdf -L/opt/intel/fce/10.1.015/lib/ -lirc
  CDFMOD = -I/share/apps/netcdf/include/ # modules for netcdf
  ifdef USE_GOTM
    GTMMOD =  -I/home/users/yinglong/GOTM/gotm-3.2.5/Intel64/modules/IFORT/ #modules
    GTMLIBS = -L/home/users/yinglong/GOTM/gotm-3.2.5/Intel64/lib/IFORT/ -lturbulence_prod  -lutil_prod
  else
    GTMMOD =
    GTMLIBS =
  endif
endif

################################################################################
# Object files; modules compile first; modules that depend on other modules compile last.
################################################################################
#Search path
VPATH = Hydro/ EcoSim/ ICM/ Sediment/ WWMII/ \
        /disk/ambcs03/0/yinglong/SELFE/from_Durski_Sept2009/V3.0c/ 

OBJS = elfe_glbl.o elfe_msgp.o 

ifdef USE_HA
  OBJS := $(OBJS) harm.o
endif

ifdef USE_ECO
  OBJS := $(OBJS) bio_param.o eclight.o biology.o bio_init.o ecosim.o read_ecoin.o specir.o
endif

ifdef USE_ICM
  OBJS := $(OBJS) icm_param.o icm_biology.o  icm.o  read_icm_input.o
endif

ifdef USE_NAPZD
  OBJS := $(OBJS) biology_napzd.o napzd_spitz.o read_napzd_input.o
endif

ifdef USE_SED
  OBJS := $(OBJS) sed_param.o sed_mod.o ocean_mod.o init_sed.o read_sed_input.o sed_init.o \
          set_vbc.o stress.o sediment_v8.o stress_soulsby.o
endif

ifdef USE_WWM
OBJS := $(OBJS) wwm_datapl.o wwm_petscpool.o wwm_aux.o wwm_mjdv2.o wwm_sparskit.o wwm_blasaux.o \
       wwm_ardhuin.o \
       wwm_gridcf.o wwm_input.o wwm_initio.o wwm_bdcons.o \
       wwm_specparam.o wwm_windinput.o wwm_dissipation.o wwm_snl3.o wwm_snl4.o wwm_babanin.o wwm_sourceterms.o \
       wwm_sourceterms.o wwm_specint.o wwm_nums1d.o wwm_numtheta.o wwm_numsigma.o wwm_fluctsplit.o \
       wwm_output.o wwm_compute.o wwm_radiation.o wwm_diffrac.o wwm_main.o
endif

OBJS := $(OBJS) grid_subs.o elfe_main.o bktrk_subs.o solver_subs.o \
        io_subs.o misc_subs.o transport.o sflux_9c.o lap.o 

################################################################################
# Library list
################################################################################
LIBS = $(MTSLIBS) $(CDFLIBS) $(GTMLIBS)

################################################################################
# Build commands for SELFE
################################################################################
#ifdef USE_SWAN
#  FLDFLAGS := $(FLDFLAGS) -ipo -xW -mp1 -W0 -assume byterecl -traceback
#endif

$(EXEC): $(OBJS) $(MAKEFILE)
	$(FLD) $(FLDFLAGS) $(OBJS) $(LIBS) -o $(EXEC)
	rm -f *.o *.mod 

%.o: %.F90 $(MAKEFILE)
	$(FCP) -c $(FCPFLAGS) $(GTMMOD) $(CDFMOD) $<

################################################################################
# Build commands for utilities
################################################################################


################################################################################
# Clean up or Clobber
################################################################################
# Clean up objects and modules
clean:
	rm -f *.F90 *.c *.o *.mod $(EXEC) 
