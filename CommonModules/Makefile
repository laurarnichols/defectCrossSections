# Makefile for Modules

include ../make.sys

#location of needed modules
MODFLAGS = -I.

MODULES = miscUtilities.o readInputFiles.o

all : commonmodules.a

commonmodules.a : $(MODULES)
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@

%.o : %.f90
	$(mpif90) -O3 -c $<

clean : 
	- /bin/rm -f *.o *.a *.f90 *.mod
