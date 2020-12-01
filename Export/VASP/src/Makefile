
-include ../../../make.sys
-include $(QE-5.3.0_Path)/make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)$(QE-5.3.0_Path)/iotk/src $(MOD_FLAG)$(QE-5.3.0_Path)/Modules \
          $(MOD_FLAG)$(QE-5.3.0_Path)/PW/src $(MOD_FLAG)$(CommonModules_srcPath) $(MOD_FLAG).

PPOBJS = \
$(QE-5.3.0_Path)/PP/src/openfil_pp.o\
$(QE-5.3.0_Path)/PP/src/stop_pp.o

EXPORTOBJS = Export_VASP_module.o Export_VASP_main.o
PWOBJS = $(QE-5.3.0_Path)/PW/src/libpw.a
QEOBJS = $(QE-5.3.0_Path)/Modules/libqemod.a $(QE-5.3.0_Path)/FFTXlib/libqefft.a
CommonMods = $(CommonModules_srcPath)/commonmodules.a

MODULES = $(PWOBJS) $(QEOBJS) $(CommonMods)

TLDEPS= bindir mods libs pw

LIBOBJS = $(QE-5.3.0_Path)/flib/ptools.a $(QE-5.3.0_Path)/flib/flib.a $(QE-5.3.0_Path)/clib/clib.a $(QE-5.3.0_Path)/iotk/src/libiotk.a

DEBUGFLAGS = -fcheck=all -Wall -Wextra -fbacktrace -Wconversion -pedantic -fbounds-check -ffpe-trap=zero,overflow,underflow,denormal,invalid -Og -g

all : tldeps Export_VASP.x 

Export_VASP.x : mods $(EXPORTOBJS) $(PPOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) $(DEBUGFLAGS) -o $@ \
		$(EXPORTOBJS) $(PPOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd $(Home_Path)/bin ; ln -fs $(Export_VASP_srcPath)/$@ . )
	
	@echo "" ; \
	echo "" ; \
	echo "Module 'Export_QE-VASP' compiled successfully ! " ; \
	echo "" ; \
	echo "" ;

mods :
	cd $(CommonModules_srcPath) ; \
	make

tldeps:
	test -n "$(TLDEPS)" && ( cd $(QE-5.3.0_Path) ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :

	@echo "" ; \
	echo "Cleaning module 'Export_VASP' ... " ; \
	/bin/rm -f *.x *.o *.mod *~ .DS_Store ; \
	echo "Cleaning module 'Export_VASP' done !" ; \
	echo "" ;
