
-include ../make.sys

all : kpointBinningByEnergy
	
kpointBinningByEnergy : kpointBinningByEnergy.o
	$(f90) -o kpointBinningByEnergy.x kpointBinningByEnergy.o

	@echo "" ;\
	echo "" ;\
	echo "kpointBinningByEnergy compiled successfully ! " ;\
	echo "" ;\
	echo "" ;

%.o : %.f90
	#$(f90) -c $<
	$(f90) -Wall -Wextra -fbacktrace -Wconversion -pedantic -fbounds-check -ffpe-trap=zero,overflow,underflow,denormal,invalid -O3 -c -g $<

clean :

	@/bin/rm -f *.x *.o *.mod .DS_Store ;

