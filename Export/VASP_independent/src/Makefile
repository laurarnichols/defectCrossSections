
-include ../../../make.sys

# location of needed modules
MODFLAGS= -I$(CommonModules_srcPath) -I.

EXPORTOBJS = Export_VASP_module.o Export_VASP_main.o
COMMONMODS = $(CommonModules_srcPath)/commonmodules.a

all : Export_VASPInd.x 

Export_VASPInd.x : mods $(COMMONMODS) $(EXPORTOBJS)
	$(mpif90) $(MODFLAGS) -o $@	$(EXPORTOBJS) $(COMMONMODS)
	- ( cd $(Home_Path)/bin ; ln -fs $(Export_VASPInd_srcPath)/$@ . )
	
	@echo "" ; \
	echo "" ; \
	echo "Module 'Export_QE-VASP' compiled successfully ! " ; \
	echo "" ; \
	echo "" ;

mods :
	cd $(CommonModules_srcPath) ; \
	make

%.o : %.f90
	$(mpif90) $(MODFLAGS) -O3 -c -fpp -D__MPI $<

clean :

	@echo "" ; \
	echo "Cleaning module 'Export_VASPInd' ... " ; \
	/bin/rm -f *.x *.o *.mod *~ .DS_Store ; \
	echo "Cleaning module 'Export_VASPInd' done !" ; \
	echo "" ;
