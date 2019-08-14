FCP = ifort
F2PY= /usr/bin/f2py
FLAGS=-fPIC -assume byterecl -check all -debug all -diag-enable -C -warn
FCPFLAGS = -assume byterecl 
PROG1 = extract_timeseries
PROG2 = extract_xyzt
PROG3 = extract_slab
PROG4 = extract_mod
OBJ1 = extract_mod.o $(PROG1).o
OBJ2 = extract_mod.o $(PROG2).o
OBJ3 = extract_mod.o $(PROG3).o
OBJ4 = extract_mod.f90

all: $(PROG1) $(PROG2) $(PROG3)

$(PROG1): $(OBJ1) $(MAKEFILE)
	$(FCP) $(FCPFLAGS) $(OBJ1) -o $@

$(PROG2): $(OBJ2) $(MAKEFILE)
	$(FCP) $(FCPFLAGS) $(OBJ2) -o $@

$(PROG3): $(OBJ3) $(MAKEFILE)
	$(FCP) $(FCPFLAGS) $(OBJ3) -o $@
#	rm -f *.o *.mod

$(PROG4): $(OBJ4) $(MAKEFILE)
	#$(FCP) $(FCPFLAGS) $(OBJ4) -o $@ 
	python2.6 $(F2PY) -m $@ -h $(PROG4).pyf $< --overwrite-signature
	python2.6 $(F2PY) --fcompiler=intelem --debug-capi -DF2PY_REPORT_ATEXIT -DF2PY_REPORT_ON_ARRAY_COPY=1 -DDEBUG_COPY_ND_ARRAY -c  $(PROG4).pyf $< --f90flags="$(FLAGS)"
 



%.o: %.f90 $(MAKEFILE)
	$(FCP) -c $(FCPFLAGS) $<

clean:
	rm -f *.o *.mod $(PROG1) $(PROG2) $(PROG3)
