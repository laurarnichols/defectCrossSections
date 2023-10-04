
-include ../../make.sys

# location of needed modules
MODFLAGS= -I$(CommonModules_srcPath) -I$(EnergyTabulator_srcPath) -I.

TME_OBJS = TME_module.o TME_main.o
COMMONMODS = $(CommonModules_srcPath)/commonmodules.a $(EnergyTabulator_srcPath)/EnergyTabulator_module.o

DEBUGFLAGS = -check all -g -warn all -traceback -gen-interfaces 

all : TME.x 
	
TME.x : mods $(COMMONMODS) $(TME_OBJS)
	$(mpif90) $(MODFLAGS) -o $@ $(TME_OBJS) $(COMMONMODS)
	- ( cd $(Home_Path)/bin ; ln -fs $(TME_srcPath)/$@ . )

	@echo "" ; \
	echo "" ; \
	echo "Module 'TME' compiled successfully !" ; \
	echo "" ; \
	echo "" ; 

mods :
	cd $(CommonModules_srcPath) ; \
	make ; \
  cd $(EnergyTabulator_srcPath) ; \
  make

%.o : %.f90
	#$(mpif90) $(MODFLAGS) -O3 -c -assume byterecl -fpp $<
	$(mpif90) $(MODFLAGS) $(DEBUGFLAGS) -O0 -c -assume byterecl -fpp $<

clean :

	@echo "" ; \
	echo "Cleaning module 'TME' ... " ; \
	/bin/rm -f *.x *.o *.mod .DS_Store ; \
	echo "Cleaning module 'TME' done !" ; \
	echo "" ;

