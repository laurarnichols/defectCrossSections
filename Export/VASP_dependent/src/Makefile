
-include ../../../make.sys

all : 
	mv $(VASP_Path)/src/main.F $(VASP_Path)/src/origMain ; \
	cp $(Export_VASPDep_srcPath)/ExportFromVASP.f90 $(VASP_Path)/src/main.F ; \
	cd $(VASP_Path) ; \
	make std || (mv $(VASP_Path)/src/origMain $(VASP_Path)/src/main.F; exit 1); \

	@mv $(VASP_Path)/src/origMain $(VASP_Path)/src/main.F ; \
	cp $(VASP_Path)/bin/vasp_std $(Export_VASPDep_srcPath)/ExportFromVASP.x  ; \
	cd $(Home_Path)/bin ; \
	ln -fs $(Export_VASPDep_srcPath)/ExportFromVASP.x . ; \
	cd $(VASP_Path) ; \
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
