
-include ../../make.sys

# location of needed modules
MODFLAGS= -I$(CommonModules_srcPath) -I.

GVEL_OBJS = GVel_module.o GVel_main.o
COMMONMODS = $(CommonModules_srcPath)/commonmodules.a

all : GVel.x 
	
GVel.x : mods $(COMMONMODS) $(GVEL_OBJS)
	$(mpif90) $(MODFLAGS) -o GVel.x $(GVEL_OBJS) $(COMMONMODS)
	- ( cd $(Home_Path)/bin ; ln -fs $(GVel_srcPath)/$@ . )

	@echo "" ; \
	echo "" ; \
	echo "Program 'GVel' compiled successfully !" ; \
	echo "" ; \
	echo "" ; 

mods :
	cd $(CommonModules_srcPath) ; \
	make

%.o : %.f90
	$(mpif90) $(MODFLAGS) -O3 -c -assume byterecl -fpp $<

clean :

	@echo "" ; \
	echo "Cleaning folder 'GVel' ... " ; \
	/bin/rm -f *.x *.o *.mod .DS_Store ; \
	echo "Cleaning folder 'GVel' done !" ; \
	echo "" ;

