
-include ../../../make.sys

LSF1_OBJS = LSF_linear_Module_v1.o LSF_linear_Main.o

all : LSF1.x
	
LSF1.x : $(LSF1_OBJS)
	$(mpif90) -o LSF1.x $(LSF1_OBJS)
	- ( cd $(Home_Path)/bin ; ln -fs $(LSF1_srcPath)/$@ . )

	@echo "" ; \
	echo "" ; \
	echo "Module 'LSF1' compiled successfully !" ; \
	echo "" ; \
	echo "" ;

%.o : %.f90
	$(mpif90) -c $<

clean :

	@echo "" ; \
	echo "Cleaning module 'LSF1' ... " ; \
	/bin/rm -f *.x *.o *.mod .DS_Store ; \
	echo "Cleaning module 'LSF1' done !" ; \
	echo "" ; \

