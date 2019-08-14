FCP = ifort
FCPFLAGS = -O2 -Bstatic -assume byterecl
PROG1 = extract_timeseries
PROG2 = extract_xyzt
PROG3 = extract_slab
OBJ1 = extract_mod.o $(PROG1).o
OBJ2 = extract_mod.o $(PROG2).o
OBJ3 = extract_mod.o $(PROG3).o

all: $(PROG1) $(PROG2) $(PROG3)

$(PROG1): $(OBJ1) $(MAKEFILE)
	$(FCP) $(FCPFLAGS) $(OBJ1) -o $@

$(PROG2): $(OBJ2) $(MAKEFILE)
	$(FCP) $(FCPFLAGS) $(OBJ2) -o $@

$(PROG3): $(OBJ3) $(MAKEFILE)
	$(FCP) $(FCPFLAGS) $(OBJ3) -o $@
	rm -f *.o *.mod

%.o: %.f90 $(MAKEFILE)
	$(FCP) -c $(FCPFLAGS) $<

clean:
	rm -f *.o *.mod $(PROG1) $(PROG2) $(PROG3)
