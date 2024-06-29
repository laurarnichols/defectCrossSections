# Makefile for Modules

include ../make.sys

#location of needed modules
MODFLAGS = -I.

MODULES = constants.o miscUtilities.o readInputFiles.o generalComputations.o errorsAndMPI.o cell.o base.o hungarian.o
DEBUGFLAGS = -check all -g -warn all -traceback -gen-interfaces 

all : commonmodules.a

commonmodules.a : $(MODULES)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

%.o : %.f90
	$(mpif90) -O3 -c $<
	#$(mpif90) $(DEBUGFLAGS) -O0 -c -assume byterecl -fpp $<

clean : 
	- /bin/rm -f *.o *.a *.mod
