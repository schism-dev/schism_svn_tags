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

PPFLAGS = -DSELFE 

########### Ensure that the user's Make.defs.local exists and include it

make_defs_local_exists := $(wildcard ../mk/Make.defs.local)
ifeq ($(strip $(make_defs_local_exists)),)
   $(error Make.defs.local does not exist in the /mk directory (parallel to this /src directory). See readme)
endif


-include ../mk/Make.defs.local
include ../mk/Make.defs.config


ifeq ($(findstring GNU Fortran,$(shell $(FCP) --version)),GNU Fortran)
  $(print Here)
  MODULEOUTDIR = -J
else
  MODULEOUTDIR = -module 
endif


############
# Derived compiler flags
###########

ifneq ($(USE_OPEN64),)
    ifeq ($(USE_64), yes)
      AVOID_ADV_WRITE := $(USE_OPEN64)
      PPFLAGS := $(PPFLAGS) -DAVOID_ADV_WRITE -DUSE_OPEN64
    endif
endif




################################################################################
# Object files; modules compile first; modules that depend on other modules compile last.
################################################################################

# rules to build objects for executable programs and libraries

#Search path
VPATH = Hydro/ EcoSim/ ICM/ Sediment/ WWMII/ \
        /disk/ambcs03/0/yinglong/SELFE/from_Durski_Sept2009/V3.0c/ 

OBJS = elfe_glbl.o elfe_msgp.o 

ifdef USE_DEBUG
  PPFLAGS := $(PPFLAGS) -DDEBUG 
endif

ifdef USE_TIMER
  PPFLAGS := $(PPFLAGS) -DINCLUDE_TIMING
endif

ifdef USE_TIMOR
    PPFLAGS := $(PPFLAGS) -DUSE_TIMOR
endif

ifdef USE_GOTM
     PPFLAGS := $(PPFLAGS) -DUSE_GOTM
endif

ifdef USE_HA
  PPFLAGS := $(PPFLAGS) -DUSE_HA
  OBJS := $(OBJS) harm.o
endif

ifdef USE_ECO
  PPFLAGS := $(PPFLAGS) -DUSE_ECO
  OBJS := $(OBJS) bio_param.o eclight.o biology.o bio_init.o ecosim.o read_ecoin.o specir.o
endif

ifdef USE_ICM
  PPFLAGS := $(PPFLAGS) -DUSE_ICM
  OBJS := $(OBJS) icm_param.o icm_biology.o  icm.o  read_icm_input.o
endif

ifdef USE_NAPZD
  PPFLAGS := $(PPFLAGS) -DUSE_NAPZD
  OBJS := $(OBJS) biology_napzd.o napzd_spitz.o read_napzd_input.o
endif

ifdef USE_SED
  PPFLAGS := $(PPFLAGS) -DUSE_SED
  OBJS := $(OBJS) sed_param.o sed_mod.o ocean_mod.o init_sed.o read_sed_input.o sed_init.o \
          set_vbc.o stress.o sediment_v8.o stress_soulsby.o
endif

ifdef USE_OIL
    PPFLAGS := $(PPFLAGS) -DUSE_OIL
endif

ifdef USE_WWM
# please take care we have now ST41 (ardhuin_old) and ST42 (ardhuin_new)
  PPFLAGS := $(PPFLAGS) -DUSE_WWM -DNCDF -DST42
  OBJS := $(OBJS) wwm_datapl.o wwm_petscpool.o wwm_aux.o wwm_mjdv2.o wwm_sparskit.o wwm_blasaux.o \
       wwm_ardhuin.o \
       wwm_gridcf.o wwm_input.o wwm_initio.o wwm_bdcons.o \
       wwm_specparam.o wwm_windinput.o wwm_dissipation.o wwm_snl3.o wwm_snl4.o wwm_babanin.o wwm_sourceterms.o \
       wwm_specint.o wwm_nums1d.o wwm_numtheta.o wwm_numsigma.o wwm_fluctsplit.o \
       wwm_output.o wwm_compute.o wwm_radiation.o wwm_diffrac.o wwm_main.o
endif

ifdef USE_TIMOR
#Add TIMOR routines here
endif

OBJS := $(OBJS) grid_subs.o elfe_main.o bktrk_subs.o solver_subs.o \
        io_subs.o misc_subs.o transport.o sflux_9c.o lap.o hydraulic.o


# Names of the object files with configuration intermediate file (.o,.mod) output directory prepended
CONFIG_OBJS :=  $(addprefix o/$(config)/, $(OBJS))

# Names of source files with relative path as found on $VPATH
FULLSRC := $(foreach f,$(subst .o,.F90,$(OBJS)),$(firstword $(wildcard $(VPATH:%=%/$f)) not-found:$f))

################################################################################
# Library list
################################################################################
LIBS = $(MTSLIBS) $(CDFLIBS) $(GTMLIBS)

################################################################################
# Build commands for SELFE
################################################################################

.PRECIOUS: $(SELFE_HOME)/mk/.check.$(config) o/$(config)/%.o o/$(config)/.dummy


# Allow for the override of executable name by the user
ifeq ($(EXEC),)
    EXEC := pelfe.$(config).ex
endif


.PHONY: pelfe

pelfe: $(EXEC) o/$(config)/.dummy


.PHONY: all
all: pelfe utility


PARMETISDIR := ./ParMetis-3.1-Sep2010

# This is the main executable target
$(EXEC): $(CONFIG_OBJS) $(PARMETISDIR)/parmetis.a $(PARMETISDIR)/metis.a
	$(FLD) $(FLDFLAGS) $(CONFIG_OBJS) $(LIBS) ${MODULEOUTDIR}. -I. -o $@


$(PARMETISDIR)/%.a:
	$(MAKE) --directory=$(PARMETISDIR)		



#ifdef USE_SWAN
#  FLDFLAGS := $(FLDFLAGS) -ipo -xW -mp1 -W0 -assume byterecl -traceback
#endif

.PHONY: depend
depend: o/$(config)/.dummy o/$(config)/.depend ../bin/.dummy
	@echo Dependencies generated


o/$(config)/.depend: o/$(config)/.dummy $(FULLSRC)
	@echo Building depend list
	@ ../mk/sfmakedepend.pl -d o/$(config) -I $(VPATH) -f o/$(config)/.depend $(FULLSRC) > depend.out

ifneq "$(MAKECMDGOALS)" "clean"
  -include o/$(config)/.depend
endif


o/$(config)/%.o: %.F90 o/$(config)/.dummy o/$(config)/.depend
	$(QUIET)$(FCP) -c $(FCPFLAGS)  ${MODULEOUTDIR}o/$(config) -Io/$(config) $(GTTMOD) $(CDFMOD) $< -o $@


# These needed to keep the compiler from using an archaic default rule for building .o and .depend files from .mod
%.o : %.mod

%.f90 : %.mod

%.depend : %.mod
	$(error You need to do a 'make depend')


%.mod:
	@touch $@


################# Utility
.PHONY: utility

util_dirs := Combining_Scripts OneWayNestScripts Particle_Tracking Pre-processing Sflux_nc
util_dirs := $(addprefix utility/,$(util_dirs))
VPATH := $(VPATH) $(util_dirs)

src_extension := .f90

# list of source files in all the other source directories
ifneq ($(util_dirs),)
  util_srcs := $(notdir $(wildcard $(foreach ext,$(src_extension),$(addsuffix /*$(ext),$(util_dirs)))))
endif

# remove a couple of items that don't compile as executables... gotta find out what to do about these
#util_srcs := $(filter-out extract_mod.f90 compute_zcor.f90,$(util_srcs))



util_execs := $(addprefix ../bin/,$(subst $(src_extension),.ex,$(util_srcs)))


util_objs := $(subst $(src_extension),.o,$(util_srcs))
util_config_objs := $(addprefix o/$(config)/, $(util_objs))

utility: $(util_execs)  o/$(config)/.dummy

#	@echo building $(util_execs) $(util_srcs)

o/$(config)/%.o : %.f90 o/$(config)/.dummy  ../bin/.dummy
	$(FCS) -c $(FCPFLAGS)  ${MODULEOUTDIR}o/$(config) -Io/$(config) $(GTTMOD) $(CDFMOD) $< -o $@

../bin/%.ex: o/$(config)/%.o  ../bin/.dummy
	$(FLD) $(FLDFLAGS) $< $(LIBS) ${MODULEOUTDIR}. -I. -o $@


#################################################################################

# This file exists solely to force creation of the subdirectories that
# hold the dependency, object and preprocessor files.
# NOTE: we cant use a dependency for this because .dummy will always
#       be out-of-date with respect to the subdirs because they are
#       modified after .dummy is created.  But since there are no
#       dependencies, this rule only runs when the file doesn't exist.
outdir_exists := $(wildcard o/$config)
make_defs_local_exists := $(wildcard ../mk/Make.defs.local)

o/$(config)/.dummy:
	@echo $(make_defs_local_exists)
        ifeq ( $(strip $(make_defs_local_exists)), )
             $(error Make.defs.local does not exist in the /mk directory. Please create one using another mk/Make.defs.xxx as a template.  \
             If you put your Make.defs.xxx file in svn, you can make a symbolic link Make.defs.local pointing to it.[1])
        endif
	$(ECHO)mkdir -p o/$(config)
	@touch o/$(config)/.dummy


../bin/.dummy:
	$(ECHO)mkdir -p ../bin
	@touch ../bin/.dummy

../mk/Make.defs.local :
	$(error Make.defs.local does not exist in the /mk directory. Please create one using another mk/Make.defs.xxx as a template. \n \
             If you put your Make.defs.xxx file in svn, you can make a symbolic link Make.defs.local pointing to it.[1])



################################################################################
# Clean up or Clobber
################################################################################
# Clean up objects and modules. Clean kills all configurations.

.PHONY: clean
clean:
	$(QUIET)rm -f $(EXEC) -f *.F90 *.c *.o *.mod $(EXEC)
	$(QUIET)rm -rf o/*
	$(QUIET)rm -rf ../bin/*



