
-include ../../make.sys

Mj_OBJS = Mj_Module_v1.o Mj_Main.o

all : Mj.x
	
Mj.x : $(Mj_OBJS)
	$(mpif90) -o Mj.x $(Mj_OBJS)
	- ( cd $(Home_Path)/bin ; ln -fs $(Mj_srcPath)/$@ . )

	@echo "" ; \
	echo "" ; \
	echo "Module 'Mj' compiled successfully !" ; \
	echo "" ; \
	echo "" ;

%.o : %.f90
	$(mpif90) -O3 -c $<

clean :

	@echo "" ; \
	echo "Cleaning module 'Mj' ... " ; \
	/bin/rm -f *.x *.o *.mod .DS_Store ; \
	echo "Cleaning module 'Mj' done !" ; \
	echo "" ; \

