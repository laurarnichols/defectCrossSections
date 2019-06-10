
-include ../../../make.sys

LSF0_OBJS = LSF_zeroth_Module_v35.o LSF_zeroth_Main.o

all : LSF0.x
	
LSF0.x : $(LSF0_OBJS)
	$(mpif90) -o LSF0.x $(LSF0_OBJS)
	- ( cd $(Home_Path)/bin ; ln -fs $(LSF0_srcPath)/$@ . )

	@echo "" ; \
	echo "" ; \
	echo "Module 'LSF0' compiled successfully !" ; \
	echo "" ; \
	echo "" ;

%.o : %.f90
	$(mpif90) -c $<

clean :

	@echo "" ; \
	echo "Cleaning module 'LSF0' ... " ; \
	/bin/rm -f *.x *.o *.mod .DS_Store ; \
	echo "Cleaning module 'LSF0' done !" ; \
	echo "" ; \

