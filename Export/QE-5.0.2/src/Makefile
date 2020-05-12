
-include ../../../../make.sys
-include $(QE-5.0.2_Path)/make.sys

# location of needed modules
MODFLAGS= $(MOD_FLAG)$(QE-5.0.2_Path)/iotk/src $(MOD_FLAG)$(QE-5.0.2_Path)/Modules \
          $(MOD_FLAG)$(QE-5.0.2_Path)/PW/src $(MOD_FLAG).

PPOBJS = \
$(QE-5.0.2_Path)/PP/src/openfil_pp.o\
$(QE-5.0.2_Path)/PP/src/stop_pp.o

PWOBJS = $(QE-5.0.2_Path)/PW/src/libpw.a
QEOBJS = $(QE-5.0.2_Path)/Modules/libqemod.a

MODULES = $(PWOBJS) $(QEOBJS)

TLDEPS= bindir mods libs pw

LIBOBJS = $(QE-5.0.2_Path)/flib/ptools.a $(QE-5.0.2_Path)/flib/flib.a $(QE-5.0.2_Path)/clib/clib.a $(QE-5.0.2_Path)/iotk/src/libiotk.a

all : tldeps Export_QE-5.0.2.x 

Export_QE-5.0.2.x : Export_QE-5.0.2_v3.o $(PPOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		Export_QE-5.0.2_v3.o $(PPOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd $(Home_Path)/bin ; ln -fs $(Home_Path)/$(Export_QE-5.0.2_srcPath)/$@ . )
	
	@echo "" ; \
	echo "" ; \
	echo "Module 'Export_QE-5.0.2' compiled successfully ! " ; \
	echo "" ; \
	echo "" ;

tldeps:
	test -n "$(TLDEPS)" && ( cd $(QE-5.0.2_Path) ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :

	@echo "" ; \
	echo "Cleaning module 'Export_QE-5.0.2' ... " ; \
	/bin/rm -f *.x *.o *.mod *~ .DS_Store ; \
	echo "Cleaning module 'Export_QE-5.0.2' done !" ; \
	echo "" ;
