
-include ../../../make.sys

# location of needed modules
MODFLAGS= -I$(CommonModules_srcPath) -I.

EXPORTOBJS = Export_VASP_module.o Export_VASP_main.o
COMMONMODS = $(CommonModules_srcPath)/commonmodules.a

DEBUGFLAGS = -check all -g -warn all -traceback -gen-interfaces 

all : Export_VASP.x 

Export_VASP.x : mods $(COMMONMODS) $(EXPORTOBJS)
	$(mpif90) $(MODFLAGS) -o $@	$(EXPORTOBJS) $(COMMONMODS)
	- ( cd $(Home_Path)/bin ; ln -fs $(Export_VASP_srcPath)/$@ . )
	
	@echo "" ; \
	echo "" ; \
	echo "Module 'Export_VASP' compiled successfully ! " ; \
	echo "" ; \
	echo "" ;

mods :
	cd $(CommonModules_srcPath) ; \
	make

%.o : %.f90
	$(mpif90) $(MODFLAGS) -O3 -c -assume byterecl -fpp $<
	#$(mpif90) $(MODFLAGS) $(DEBUGFLAGS) -O0 -c -assume byterecl -fpp $<

clean :

	@echo "" ; \
	echo "Cleaning module 'Export_VASP' ... " ; \
	/bin/rm -f *.x *.o *.mod *~ .DS_Store ; \
	echo "Cleaning module 'Export_VASP' done !" ; \
	echo "" ;
