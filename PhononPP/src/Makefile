
-include ../../make.sys

# location of needed modules
MODFLAGS= -I$(CommonModules_srcPath) -I.

SHIFTER_OBJS = Shifter_module.o Shifter_main.o
COMMONMODS = $(CommonModules_srcPath)/commonmodules.a

all : Shifter.x 
	
Shifter.x : mods $(COMMONMODS) $(SHIFTER_OBJS)
	$(mpif90) $(MODFLAGS) -o Shifter.x $(SHIFTER_OBJS) $(COMMONMODS)
	- ( cd $(Home_Path)/bin ; ln -fs $(PhononPP_srcPath)/$@ . )

	@echo "" ; \
	echo "" ; \
	echo "Program 'Shifter' compiled successfully !" ; \
	echo "" ; \
	echo "" ; 

mods :
	cd $(CommonModules_srcPath) ; \
	make

%.o : %.f90
	$(mpif90) $(MODFLAGS) -O3 -c -assume byterecl -fpp $<

clean :

	@echo "" ; \
	echo "Cleaning folder 'PhononPP' ... " ; \
	/bin/rm -f *.x *.o *.mod .DS_Store ; \
	echo "Cleaning folder 'PhononPP' done !" ; \
	echo "" ;

