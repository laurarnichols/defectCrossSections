
-include ../../make.sys

# location of needed modules
MODFLAGS = -I../../CommonModules -I.

Mj_OBJS = Mj_Module_v1.o Mj_Main.o

CommonMods = $(CommonModules_srcPath)/commonmodules.a

all : Mj.x
	
Mj.x : mods $(CommonMods) $(Mj_OBJS)
	$(mpif90) $(LDFLAGS) $(MODFLAGS) -o $@ $(Mj_OBJS) $(CommonMods)
	- ( cd $(Home_Path)/bin ; ln -fs $(Mj_srcPath)/$@ . )

	@echo "" ; \
	echo "" ; \
	echo "Module 'Mj' compiled successfully !" ; \
	echo "" ; \
	echo "" ;

mods :
	cd ../../CommonModules ; \
	make

%.o : %.f90
	$(mpif90) $(MODFLAGS) -O3 -c $<

clean :

	@echo "" ; \
	echo "Cleaning module 'Mj' ... " ; \
	/bin/rm -f *.x *.o *.mod .DS_Store ; \
	echo "Cleaning module 'Mj' done !" ; \
	echo "" ; \

include make.depend
