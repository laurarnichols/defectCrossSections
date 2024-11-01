
-include ../../make.sys

# location of needed modules
MODFLAGS= -I$(CommonModules_srcPath) -I$(TME_srcPath) -I.

DEBUGFLAGS = -check all -g -warn all -traceback -gen-interfaces 

ENERGYTAB_OBJS = EnergyTabulator_module.o EnergyTabulator_main.o
COMMONMODS = $(TME_srcPath)/TME_module.o $(CommonModules_srcPath)/commonmodules.a

all : EnergyTabulator.x 
	
EnergyTabulator.x : mods $(COMMONMODS) $(ENERGYTAB_OBJS)
	$(mpif90) $(MODFLAGS) -o EnergyTabulator.x $(ENERGYTAB_OBJS) $(COMMONMODS)
	- ( cd $(Home_Path)/bin ; ln -fs $(EnergyTabulator_srcPath)/$@ . )

	@echo "" ; \
	echo "" ; \
	echo "Program 'EnergyTabulator' compiled successfully !" ; \
	echo "" ; \
	echo "" ; 

mods :
	cd $(CommonModules_srcPath) ; \
	make ; \
  cd $(TME_srcPath) ; \
  make TMENoMods

%.o : %.f90
	$(mpif90) $(MODFLAGS) -O3 -c -assume byterecl -fpp $<
	#$(mpif90) $(MODFLAGS) $(DEBUGFLAGS) -O0 -c -assume byterecl -fpp $<

clean :

	@echo "" ; \
	echo "Cleaning folder 'EnergyTabulator' ... " ; \
	/bin/rm -f *.x *.o *.mod .DS_Store ; \
	echo "Cleaning folder 'EnergyTabulator' done !" ; \
	echo "" ;

