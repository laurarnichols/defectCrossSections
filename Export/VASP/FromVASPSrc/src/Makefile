
-include ../../../../make.sys

all : 
	cp $(ExportFromVASPSrc_srcPath)/Export_FromVASPSrc.f90 $(VASPCompilation_Path)/src/main.F ; \
	cp $(ExportFromVASPSrc_srcPath)/poscar.F $(VASPCompilation_Path)/src/poscar.F ; \
	cp $(ExportFromVASPSrc_srcPath)/base.F $(VASPCompilation_Path)/src/base.F ; \
	cp $(ExportFromVASPSrc_srcPath)/fileio.F $(VASPCompilation_Path)/src/fileio.F ; \
	cp $(ExportFromVASPSrc_srcPath)/pseudo.F $(VASPCompilation_Path)/src/pseudo.F ; \
	cp $(ExportFromVASPSrc_srcPath)/mkpoints.F $(VASPCompilation_Path)/src/mkpoints.F ; \
	cd $(VASPCompilation_Path) ; \
	make std || exit 1; \

	@cp $(VASPCompilation_Path)/bin/vasp_std $(ExportFromVASPSrc_srcPath)/Export_FromVASPSrc.x  ; \
	cd $(Home_Path)/bin ; \
	ln -fs $(ExportFromVASPSrc_srcPath)/Export_FromVASPSrc.x . ; \
	echo "" ; \
	echo "" ; \
	echo "Module 'ExportFromVASPSrc' compiled successfully ! " ; \
	echo "" ; \
	echo "" ;

clean :

	@echo "" ; \
	echo "Cleaning module 'ExportFromVASPSrc' ... " ; \
	/bin/rm -f *.x *.o *.mod *~ .DS_Store ; \
	echo "Cleaning module 'ExportFromVASPSrc' done !" ; \
	echo "" ;
