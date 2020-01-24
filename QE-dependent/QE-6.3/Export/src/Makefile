
-include ../../../../make.sys
-include $(QE-6.3_Path)/make.inc

# location of needed modules
#MODFLAGS= $(MOD_FLAG)$(QE-6.3_Path)/iotk/src $(MOD_FLAG)$(QE-6.3_Path)/Modules \
          $(MOD_FLAG)$(QE-6.3_Path)/PW/src $(MOD_FLAG)$(QE-6.3_Path)/UtilXlib \
          $(MOD_FLAG)$(QE-6.3_Path)/KS_Solvers/Davidson \
          $(MOD_FLAG)$(QE-6.3_Path)/KS_Solvers/CG $(MOD_FLAG)$(QE-6.3_Path)/dft-d3/ \
          $(MOD_FLAG).
MODFLAGS= $(BASEMOD_FLAGS) \
          $(MOD_FLAG)$(QE-6.3_Path)/PW/src \
          $(MOD_FLAG)$(QE-6.3_Path)/UtilXlib \
          $(MOD_FLAG)$(QE-6.3_Path)/KS_Solvers/Davidson \
          $(MOD_FLAG)$(QE-6.3_Path)/KS_Solvers/CG \
          $(MOD_FLAG)$(QE-6.3_Path)/dft-d3/ \
          $(MOD_FLAG) .


PPOBJS = \
$(QE-6.3_Path)/PP/src/openfil_pp.o\
$(QE-6.3_Path)/PP/src/stop_pp.o



PWOBJS = $(QE-6.3_Path)/PW/src/libpw.a
#QEOBJS = $(QE-6.3_Path)/Modules/libqemod.a $(QE-6.3_Path)/FFTXlib/libqefft.a $(QE-6.3_Path)/UtilXlib/libutil.a \
         $(QE-6.3_Path)/LAXlib/libqela.a $(QE-6.3_Path)/KS_Solvers/Davidson/libdavid.a \
         $(QE-6.3_Path)/KS_Solvers/CG/libcg.a

QEOBJS = $(QE-6.3_Path)/Modules/libqemod.a $(QE-6.3_Path)/KS_Solvers/Davidson/libdavid.a $(QE-6.3_Path)/KS_Solvers/CG/libcg.a \
         $(QE-6.3_Path)/FFTXlib/libqefft.a $(QE-6.3_Path)/LAXlib/libqela.a $(QE-6.3_Path)/UtilXlib/libutil.a

MODULES = $(PWOBJS) $(QEOBJS)

TLDEPS= bindir mods libs pw

# LIBOBJS = $(QE-6.3_Path)/flib/ptools.a $(QE-6.3_Path)/flib/flib.a $(QE-6.3_Path)/clib/clib.a $(QE-6.3_Path)/iotk/src/libiotk.a

all : tldeps Export_QE-6.3.x 

Export_QE-6.3.x : Export_QE-6.3.o $(PPOBJS) $(MODULES) $(LIBOBJS)
	$(LD) $(LDFLAGS) -o $@ \
		Export_QE-6.3.o $(PPOBJS) $(MODULES) $(LIBOBJS) $(LIBS)
	- ( cd $(Home_Path)/bin ; ln -fs $(Export_QE-6.3_srcPath)/$@ . )
	
	@echo "" ; \
	echo "" ; \
	echo "Module 'Export_QE-6.3' compiled successfully ! " ; \
	echo "" ; \
	echo "" ;

tldeps:
	test -n "$(TLDEPS)" && ( cd $(QE-6.3_Path) ; $(MAKE) $(MFLAGS) $(TLDEPS) || exit 1) || :

clean :

	@echo "" ; \
	echo "Cleaning module 'Export_QE-6.3' ... " ; \
	/bin/rm -f *.x *.o *.mod *~ .DS_Store ; \
	echo "Cleaning module 'Export_QE-6.3' done !" ; \
	echo "" ;
