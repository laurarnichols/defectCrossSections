
-include ../../make.sys

TME_OBJS = TME_Module_v28.o TME_Main_v9.o

all : TME.x 
	
TME.x : $(TME_OBJS)
	$(mpif90) -o TME.x $(TME_OBJS)
	- ( cd $(Home_Path)/bin ; ln -fs $(TME_srcPath)/$@ . )

	@echo "" ; \
	echo "" ; \
	echo "Module 'TME' compiled successfully !" ; \
	echo "" ; \
	echo "" ; 

%.o : %.f90
	$(mpif90) -O3 -c $<

clean :

	@echo "" ; \
	echo "Cleaning module 'TME' ... " ; \
	/bin/rm -f *.x *.o *.mod .DS_Store ; \
	echo "Cleaning module 'TME' done !" ; \
	echo "" ;

