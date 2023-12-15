

########################################################################
#                                                                      #
# Please change the following declarations to meet yours system setup. #
#                                                                      #
########################################################################

# Machine-dependent compiler wrappers
f90 = ifort
mpif90 = ftn


########################################################################
#                                                                      #
#       Please do NOT make any change in the following lines!          #
#                                                                      #
########################################################################


all : initialize Export_VASP EnergyTabulator TME PhononPP LSF

help :

	@echo ""
	@echo ""
	@echo " Edit this Makefile with the intel compiler wrappers for your machine. "
	@echo ""
	@echo ""
	@echo " Please use one of the following commands :"
	@echo ""
	@echo "    make Export                   to build the VASP Export program."
	@echo "    make EnergyTabulator          to build the energy tabulation program."
	@echo "    make TME                      to build the Transition Matrix Elements (TME) program."
	@echo "    make PhononPP                 to build the phonon post-processing program."
	@echo "    make LSF                      to build the Line Shape Function (LSF) program."
	@echo ""
	@echo ""
	@echo "    make clean                    to clean all the modules of the package."
	@echo "    make cleanExport              to clean the VASP Export program."
	@echo "    make cleanEnergyTabulator     to clean the energy tabulation program."
	@echo "    make cleanTME                 to clean the Transition Matrix Elements (TME) program."
	@echo "    make cleanPhononPP            to clean the phonon post-processing program."
	@echo "    make cleanLSF                 to clean the Line Shape Function (LSF) program."
	@echo ""
	@echo ""


Export_VASP_srcPath     = Export/VASP/src
EnergyTabulator_srcPath = EnergyTabulator/src
TME_srcPath             = TME/src
LSF_srcPath             = LSF/src
PhononPP_srcPath        = PhononPP/src

bin = './bin'

initialize :

	@echo "" > make.sys ; \
	echo "Home_Path               = " $(PWD) >> make.sys ; \
	echo "Export_VASP_srcPath     = " $(PWD)/$(Export_VASP_srcPath) >> make.sys ; \
	echo "EnergyTabulator_srcPath = " $(PWD)/$(EnergyTabulator_srcPath) >> make.sys ; \
	echo "TME_srcPath             = " $(PWD)/$(TME_srcPath) >> make.sys ; \
	echo "PhononPP_srcPath        = " $(PWD)/$(PhononPP_srcPath) >> make.sys ; \
	echo "LSF_srcPath             = " $(PWD)/$(LSF_srcPath) >> make.sys ; \
	echo "CommonModules_srcPath   = " $(PWD)/$(CommonModules_srcPath) >> make.sys ; \
	echo "" >> make.sys ; \
	echo "f90    = "$(f90) >> make.sys ; \
	echo "mpif90 = "$(mpif90) >> make.sys ; \
	echo "" >> make.sys ; \
  	echo "# ar command and flags - for most architectures: AR = ar, ARFLAGS = ruv" >> make.sys ; \
  	echo "" >> make.sys ; \
  	echo "AR = ar" >> make.sys ; \
  	echo "ARFLAGS = ruv" >> make.sys ; \
  	echo "" >> make.sys ; \
  	echo "# ranlib command" >> make.sys ; \
  	echo "" >> make.sys ; \
  	echo "RANLIB = ranlib" >> make.sys
#
	@if test ! -d $(bin) ; then \
		mkdir $(bin) ; \
	fi


# Define the targets
TARGETS := Export_VASP EnergyTabulator TME PhononPP LSF

# Define a pattern rule for all targets
$(TARGETS) : initialize

  @cd $($@_srcPath) ; \
    make all

clean : cleanInitialization cleanExport_VASP cleanEnergyTabulator cleanGVel cleanTME cleanPhononPP cleanMj cleanLSF cleanSigma

cleanInitialization :

	@if test -d bin ; then \
		/bin/rm -rf bin ; \
	fi
	@if test -e make.sys ; then \
		/bin/rm -f make.sys ; \
	fi

# Define a pattern rule for cleaning targets
clean% : initialize

  @cd $($@_srcPath) ; \
    make clean
