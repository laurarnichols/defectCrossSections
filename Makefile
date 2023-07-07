

########################################################################
#                                                                      #
# Please change the following declarations to meet yours system setup. #
#                                                                      #
########################################################################

QE-5.0.2_Path = ${HOME}/q-e-qe-5.0.2

QE-5.3.0_Path = ${HOME}/q-e-qe-5.3

QE-6.3_Path = ${HOME}/qe-6.3

f90 = ifort

mpif90 = ftn


########################################################################
#                                                                      #
#       Please do NOT make any change in the following lines!          #
#                                                                      #
########################################################################


all : initialize QE-5.0.2_dependent QE-5.3.0_dependent QE-6.3_dependent Export_VASP EnergyTabulator TME LSF Mj Sigma

help :

	@echo ""
	@echo ""
	@echo " If you want to built Quantum Espresso dependent modules, "
	@echo " make sure that PW and PP packages are already built."
	@echo ""
	@echo " Then edit this Makefile to meet yours system setup. "
	@echo " You have to define the path to the Quantum Espresso home directory,"
	@echo " and the fortran compilers you want to be used."
	@echo ""
	@echo ""
	@echo " Please use one of the following commands :"
	@echo ""
	@echo "    make all_QE-5.0.2             to built all the modules of the package using Quantum Espresso 5.0.2."
	@echo "    make all_QE-5.3.0             to built all the modules of the package using Quantum Espresso 5.3.0."
	@echo "    make all_QE-6.3               to built all the modules of the package using Quantum Espresso 6.3."
	@echo "    make QE-5.0.2_dependent       to built all the Quantum Espresso 5.0.2 dependent modules."
	@echo "    make QE-5.3.0_dependent       to built all the Quantum Espresso 5.3.0 dependent modules."
	@echo "    make QE-6.3_dependent         to built all the Quantum Espresso 6.3 dependent modules."
	@echo "    make Export_QE-5.0.2          to built the Quantum Espresso 5.0.2 dependent Export module."
	@echo "    make Export_QE-5.3.0          to built the Quantum Espresso 5.3.0 dependent Export module."
	@echo "    make Export_QE-6.3            to built the Quantum Espresso 6.3 dependent Export module."
	@echo "    make Export_VASP              to built the VASP Export code"
	@echo "    make EnergyTabulator          to built the energy tabulation code."
	@echo "    make GVel                     to built the group velocity code."
	@echo "    make TME                      to built the Transition Matrix Elements (TME) module."
	@echo "    make PhononPP                 to built the phonon post-processing module."
	@echo "    make Mj                       to built the Mj module."
	@echo "    make LSF                      to built the Line Shape Function (LSF) module."
	@echo "    make Sigma                    to built the Cross Section (Sigma) module."
	@echo ""
	@echo ""
	@echo "    make clean_all_QE-5.3.0            to clean all the modules of the package using Quantum Espresso 5.3.0."
	@echo "    make clean_all_QE-6.3              to clean all the modules of the package using Quantum Espresso 6.3."
	@echo "    make clean_QE-5.3.0_dependent      to clean all the Quantum Espresso 5.3.0 dependent modules."
	@echo "    make clean_QE-6.3_dependent        to clean all the Quantum Espresso 6.3 dependent modules."
	@echo "    make cleanExport_QE-5.3.0          to clean the Quantum Espresso 5.3.0 dependent Export module."
	@echo "    make cleanExport_QE-6.3            to clean the Quantum Espresso 6.3 dependent Export module."
	@echo "    make cleanPhononPP                 to clean the phonon post-processing module."
	@echo "    make cleanMj                       to clean the Mj module."
	@echo "    make cleanEnergyTabulator          to clean the energy tabulation code."
	@echo "    make cleanGVel                     to clean the group velocity code."
	@echo "    make cleanTME                      to clean the Transition Matrix Elements (TME) module."
	@echo "    make cleanLSF                      to clean the Line Shape Function (LSF) module."
	@echo "    make cleanSigma                    to clean the Cross Section (Sigma) module."
	@echo ""
	@echo ""


Export_QE-5.0.2_srcPath = Export/QE/5.0.2/src
Export_QE-5.3.0_srcPath = Export/QE/5.3/src
Export_QE-6.3_srcPath   = Export/QE/6.3/src
Export_VASP_srcPath     = Export/VASP/src
EnergyTabulator_srcPath = EnergyTabulator/src
GVel_srcPath            = GVel/src
TME_srcPath             = TME/src
Mj_srcPath              = Mj/src
LSF_srcPath             = LSF/src
Sigma_srcPath           = Sigma/src
CommonModules_srcPath   = CommonModules
PhononPP_srcPath        = PhononPP/src

bin = './bin'

all_QE-5.0.2 : initialize QE-5.0.2_dependent base

all_QE-5.3.0 : initialize QE-5.3.0_dependent base

all_QE-6.3 : initialize QE-6.3_dependent base

all_VASP : initialize Export_VASP base

QE-5.0.2_dependent : initialize Export_QE-5.0.2

QE-5.3.0_dependent : initialize Export_QE-5.3.0

QE-6.3_dependent : initialize Export_QE-6.3

base : EnergyTabulator GVel TME LSF Mj PhononPP Sigma

initialize :

	@echo "" > make.sys ; \
	echo "Home_Path               = " $(PWD) >> make.sys ; \
	echo "QE-5.0.2_Path           = " $(QE-5.0.2_Path) >> make.sys ; \
	echo "QE-5.3.0_Path           = " $(QE-5.3.0_Path) >> make.sys ; \
	echo "QE-6.3_Path             = " $(QE-6.3_Path) >> make.sys ; \
	echo "Export_QE-5.0.2_srcPath = " $(PWD)/$(Export_QE-5.0.2_srcPath) >> make.sys ; \
	echo "Export_QE-5.3.0_srcPath = " $(PWD)/$(Export_QE-5.3.0_srcPath) >> make.sys ; \
	echo "Export_QE-6.3_srcPath   = " $(PWD)/$(Export_QE-6.3_srcPath) >> make.sys ; \
	echo "Export_VASP_srcPath     = " $(PWD)/$(Export_VASP_srcPath) >> make.sys ; \
	echo "EnergyTabulator_srcPath = " $(PWD)/$(EnergyTabulator_srcPath) >> make.sys ; \
	echo "GVel_srcPath            = " $(PWD)/$(GVel_srcPath) >> make.sys ; \
	echo "TME_srcPath             = " $(PWD)/$(TME_srcPath) >> make.sys ; \
	echo "PhononPP_srcPath        = " $(PWD)/$(PhononPP_srcPath) >> make.sys ; \
	echo "Mj_srcPath              = " $(PWD)/$(Mj_srcPath) >> make.sys ; \
	echo "LSF_srcPath             = " $(PWD)/$(LSF_srcPath) >> make.sys ; \
	echo "Sigma_srcPath           = " $(PWD)/$(Sigma_srcPath) >> make.sys ; \
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


Export_QE-5.0.2 : initialize

	@cd $(Export_QE-5.0.2_srcPath) ; \
		make all

Export_QE-5.3.0 : initialize

	@cd $(Export_QE-5.3.0_srcPath) ; \
		make all

Export_QE-6.3 : initialize

	@cd $(Export_QE-6.3_srcPath) ; \
		make all

Export_VASP : initialize

	@cd $(Export_VASP_srcPath) ; \
		make all

EnergyTabulator : initialize

	@cd $(EnergyTabulator_srcPath) ; \
        	make all

GVel : initialize

	@cd $(GVel_srcPath) ; \
        	make all

TME : initialize

	@cd $(TME_srcPath) ; \
        	make all

PhononPP : initialize

	@cd $(PhononPP_srcPath) ; \
        	make all

Mj : initialize

	@cd $(Mj_srcPath) ; \
        	make all

LSF : initialize

	@cd $(LSF_srcPath) ; \
        	make all

Sigma : initialize

	@cd $(Sigma_srcPath) ; \
        	make all

clean : cleanExport_VASP cleanExport_QE-5.0.2 cleanExport_QE-5.3.0 cleanExport_QE-6.3 cleanBase

clean_all_QE-5.0.2 : cleanExport_QE-5.0.2 cleanBase

clean_all_QE-5.3.0 : cleanExport_QE-5.3.0 cleanBase

clean_all_QE-6.3 : cleanExport_QE-6.3 cleanBase

cleanBase : cleanInitialization cleanEnergyTabulator cleanGVel cleanTME cleanPhononPP cleanMj cleanLSF cleanSigma

cleanInitialization :

	@if test -d bin ; then \
		/bin/rm -rf bin ; \
	fi
	@if test -e make.sys ; then \
		/bin/rm -f make.sys ; \
	fi

cleanExport_QE-5.0.2 :

	@cd $(Export_QE-5.0.2_srcPath) ; \
        	make clean

cleanExport_QE-5.3.0 :

	@cd $(Export_QE-5.3.0_srcPath) ; \
        	make clean

cleanExport_QE-6.3 :

	@cd $(Export_QE-6.3_srcPath) ; \
        	make clean

cleanExport_VASP :

	@cd $(Export_VASP_srcPath) ; \
        	make clean

cleanEnergyTabulator :

	@cd $(EnergyTabulator_srcPath) ; \
		make clean

cleanGVel :

	@cd $(GVel_srcPath) ; \
		make clean

cleanTME :

	@cd $(TME_srcPath) ; \
		make clean

cleanPhononPP :

	@cd $(PhononPP_srcPath) ; \
		make clean

cleanMj :

	@cd $(Mj_srcPath) ; \
		make clean


cleanLSF :

	@cd $(LSF_srcPath) ; \
		make clean

cleanSigma :

	@cd $(Sigma_srcPath) ; \
		make clean

