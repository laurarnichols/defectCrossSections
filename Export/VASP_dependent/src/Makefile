
-include ../../../make.sys

all : 
	cp ExportFromVASP.f90 $(VASP_Path)/src ; \
	cd $(VASP_Path)/src ; \
	mv main.F origMain ; \
	mv ExportFromVASP.f90 main.F ; \
	cd ../ ; \
	make std ; \
	cp bin/vasp_std $(Export_VASPDep_srcPath)/ExportFromVASP.x  ; \
	cd $(Home_Path)/bin ; \
	ln -fs $(Export_VASPDep_srcPath)/ExportFromVASP.x . ; \
	cd $(VASP_Path)/src ; \
	mv origMain main.F ; \
	cd ../ ; \
	make std ; \
	echo "" ; \
	echo "" ; \
	echo "Module 'ExportFromVASP' compiled successfully ! " ; \
	echo "" ; \
	echo "" ;

clean :

	@echo "" ; \
	echo "Cleaning module 'ExportFromVASP' ... " ; \
	/bin/rm -f *.x *.o *.mod *~ .DS_Store ; \
	echo "Cleaning module 'ExportFromVASP' done !" ; \
	echo "" ;
