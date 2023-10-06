
-include ../../make.sys

# location of needed modules
MODFLAGS= -I$(CommonModules_srcPath) -I$(EnergyTabulator_srcPath) -I$(TME_srcPath) -I.

DEBUGFLAGS = -check all -g -warn all -traceback -gen-interfaces 

LSF_OBJS = LSF_module.o LSF_main.o
COMMONMODS = $(CommonModules_srcPath)/commonmodules.a $(EnergyTabulator_srcPath)/EnergyTabulator_module.o $(TME_srcPath)/TME_module.o

all : LSF.x 
	
LSF.x : mods $(COMMONMODS) $(LSF_OBJS)
	$(mpif90) $(MODFLAGS) -o $@ $(LSF_OBJS) $(COMMONMODS)
	- ( cd $(Home_Path)/bin ; ln -fs $(LSF_srcPath)/$@ . )

	@echo "" ; \
	echo "" ; \
	echo "Module 'LSF' compiled successfully !" ; \
	echo "" ; \
	echo "" ;

mods :
	cd $(CommonModules_srcPath) ; \
	make ; \
  cd $(TME_srcPath) ; \
	make ; \
  cd $(EnergyTabulator_srcPath) ; \
  make

%.o : %.f90
	$(mpif90) $(MODFLAGS) -O3 -c -assume byterecl -fpp $<
	#$(mpif90) $(MODFLAGS) $(DEBUGFLAGS) -O0 -c -assume byterecl -fpp $<

clean :

	@echo "" ; \
	echo "Cleaning module 'LSF' ... " ; \
	/bin/rm -f *.x *.o *.mod .DS_Store ; \
	echo "Cleaning module 'LSF' done !" ; \
	echo "" ; \

