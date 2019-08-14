#
# Generic Makefile for single file Fortran codes
#
# Example: make readncom8b

%:
	ifort -g -O3 -Bstatic -o $* $*.f90 -Vaxlib -I/usr/local/include -L/usr/local/lib -lnetcdf
