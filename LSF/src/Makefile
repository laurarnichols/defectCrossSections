
-include ../../make.sys

# location of needed modules
MODFLAGS = -I$(CommonModules_srcPath) -I.

LSF_OBJS = LSF_module.o LSF_main.o

CommonMods = $(CommonModules_srcPath)/commonmodules.a

all : LSF.x 
	
LSF.x : mods $(CommonMods) $(LSF_OBJS)
	$(mpif90) $(LDFLAGS) $(MODFLAGS) -o $@ $(LSF_OBJS) $(CommonMods)
	- ( cd $(Home_Path)/bin ; ln -fs $(LSF_srcPath)/$@ . )

	@echo "" ; \
	echo "" ; \
	echo "Module 'LSF' compiled successfully !" ; \
	echo "" ; \
	echo "" ;

mods :
	cd $(CommonModules_srcPath) ; \
	make

%.o : %.f90
	$(mpif90) $(MODFLAGS) -c -g $<

clean :

	@echo "" ; \
	echo "Cleaning module 'LSF' ... " ; \
	/bin/rm -f *.x *.o *.mod .DS_Store ; \
	echo "Cleaning module 'LSF' done !" ; \
	echo "" ; \

