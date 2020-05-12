
-include ../../../make.sys
-include $(QE-5.3.0_Path)/make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)$(QE-5.3.0_Path)/iotk/src $(MOD_FLAG)$(QE-5.3.0_Path)/Modules \
          $(MOD_FLAG)$(QE-5.3.0_Path)/PW/src $(MOD_FLAG).

PPOBJS = \
$(QE-5.3.0_Path)/PP/src/openfil_pp.o\
$(QE-5.3.0_Path)/PP/src/stop_pp.o



PWOBJS = $(QE-5.3.0_Path)/PW/src/libpw.a
QEOBJS = $(QE-5.3.0_Path)/Modules/libqemod.a $(QE-5.3.0_Path)/FFTXlib/libqefft.a

MODULES = $(PWOBJS) $(QEOBJS)

TLDEPS= bindir mods libs pw

LIBOBJS = $(QE-5.3.0_Path)/flib/ptools.a $(QE-5.3.0_Path)/flib/flib.a $(QE-5.3.0_Path)/clib/clib.a $(QE-5.3.0_Path)/iotk/src/libiotk.a

all : tldeps Export_QE-5.3.0.x 

Export_QE-5.3.0.x : Export_QE-5.3.0_v3.o $(PPOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		Export_QE-5.3.0_v3.o $(PPOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd $(Home_Path)/bin ; ln -fs $(Export_QE-5.3.0_srcPath)/$@ . )
	
	@echo "" ; \
	echo "" ; \
	echo "Module 'Export_QE-5.3.0' compiled successfully ! " ; \
	echo "" ; \
	echo "" ;

tldeps:
	test -n "$(TLDEPS)" && ( cd $(QE-5.3.0_Path) ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :

	@echo "" ; \
	echo "Cleaning module 'Export_QE-5.3.0' ... " ; \
	/bin/rm -f *.x *.o *.mod *~ .DS_Store ; \
	echo "Cleaning module 'Export_QE-5.3.0' done !" ; \
	echo "" ;
