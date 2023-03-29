# Makefile for Modules

include ../make.sys

#location of needed modules
MODFLAGS = -I.

MODULES = constants.o miscUtilities.o readInputFiles.o generalComputations.o errorsAndMPI.o cell.o

all : commonmodules.a

commonmodules.a : $(MODULES)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

%.o : %.f90
	$(mpif90) -O3 -c $<

clean : 
	- /bin/rm -f *.o *.a *.mod
