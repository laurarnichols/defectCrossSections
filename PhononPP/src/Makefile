
-include ../../make.sys

# location of needed modules
MODFLAGS= -I$(CommonModules_srcPath) -I$(EnergyTabulator_srcPath) -I.

DEBUGFLAGS = -check all -g -warn all -traceback -gen-interfaces 

PHONONPP_OBJS = PhononPP_module.o PhononPP_main.o
COMMONMODS = $(EnergyTabulator_srcPath)/EnergyTabulator_module.o $(CommonModules_srcPath)/commonmodules.a

all : PhononPP.x 
	
PhononPP.x : mods $(COMMONMODS) $(PHONONPP_OBJS)
	$(mpif90) $(MODFLAGS) -o $@ $(PHONONPP_OBJS) $(COMMONMODS)
	- ( cd $(Home_Path)/bin ; ln -fs $(PhononPP_srcPath)/$@ . )

	@echo "" ; \
	echo "" ; \
	echo "Program 'PhononPP' compiled successfully !" ; \
	echo "" ; \
	echo "" ; 

mods :
	cd $(CommonModules_srcPath) ; \
	make ; \
  cd $(EnergyTabulator_srcPath) ; \
  make

%.o : %.f90
	$(mpif90) $(MODFLAGS) -O3 -c -assume byterecl -fpp $<
	#$(mpif90) $(MODFLAGS) $(DEBUGFLAGS) -O0 -c -assume byterecl -fpp $<

clean :

	@echo "" ; \
	echo "Cleaning folder 'PhononPP' ... " ; \
	/bin/rm -f *.x *.o *.mod .DS_Store ; \
	echo "Cleaning folder 'PhononPP' done !" ; \
	echo "" ;

