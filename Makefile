

########################################################################
#                                                                      #
# Please change the following declarations to meet yours system setup. #
#                                                                      #
########################################################################

QE-5.3.0_Path = ${HOME}/q-e-qe-5.3

f90 = ftn

mpif90 = ftn


########################################################################
#                                                                      #
#       Please do NOT make any change in the following lines!          #
#                                                                      #
########################################################################


default :

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
	@echo "    make all_QE-5.3.0             to built all the modules of the package using Quantum Espresso 5.3.0."
	@echo "    make QE-5.3.0_dependent       to built all the Quantum Espresso 5.3.0 dependent modules."
	@echo "    make Export_QE-5.3.0          to built the Quantum Espresso 5.3.0 dependent Export module."
	@echo "    make TME                      to built the Transition Matrix Elements (TME) module."
	@echo "    make Mj                       to built the Mj module."
	@echo "    make LSF                      to built the Line Shape Function (LSF) module."
	@echo "    make Sigma                    to built the Cross Section (Sigma) module."
	@echo ""
	@echo ""
	@echo "    make clean_all_QE-5.3.0            to clean all the modules of the package using Quantum Espresso 5.3.0."
	@echo "    make clean_QE-5.3.0_dependent      to clean all the Quantum Espresso 5.3.0 dependent modules."
	@echo "    make cleanExport_QE-5.3.0          to clean the Quantum Espresso 5.3.0 dependent Export module."
	@echo "    make cleanMj                       to clean the Mj module."
	@echo "    make cleanTME                      to clean the Transition Matrix Elements (TME) module."
	@echo "    make cleanLSF                      to clean the Line Shape Function (LSF) module."
	@echo "    make cleanSigma                    to clean the Cross Section (Sigma) module."
	@echo ""
	@echo ""


Export_QE-5.3.0_srcPath = QE-dependent/QE-5.3.0/Export/src
TME_srcPath    = TME/src
Mj_srcPath     = Mj/src
LSF_srcPath    = LSF/src
LSF0_srcPath   = $(LSF_srcPath)/zerothOrder
LSF1_srcPath   = $(LSF_srcPath)/linearTerms
Sigma_srcPath  = Sigma/src
CommonModules_srcPath  = CommonModules

bin = './bin'

all_QE-5.3.0 : initialize QE-5.3.0_dependent TME LSF0 Mj LSF1 Sigma

QE-5.3.0_dependent : initialize Export_QE-5.3.0

initialize :

	@echo "" > make.sys ; \
	echo "Home_Path      = " $(PWD) >> make.sys ; \
	echo "QE-5.3.0_Path           = " $(QE-5.3.0_Path) >> make.sys ; \
	echo "Export_QE-5.3.0_srcPath = " $(PWD)/$(Export_QE-5.3.0_srcPath) >> make.sys ; \
	echo "TME_srcPath    = " $(PWD)/$(TME_srcPath) >> make.sys ; \
	echo "Mj_srcPath    = " $(PWD)/$(Mj_srcPath) >> make.sys ; \
	echo "LSF_srcPath = " $(PWD)/$(LSF_srcPath) >> make.sys ; \
	echo "LSF0_srcPath = " $(PWD)/$(LSF0_srcPath) >> make.sys ; \
	echo "LSF1_srcPath = " $(PWD)/$(LSF1_srcPath) >> make.sys ; \
	echo "Sigma_srcPath = " $(PWD)/$(Sigma_srcPath) >> make.sys ; \
	echo "CommonModules_srcPath = " $(PWD)/$(CommonModules_srcPath) >> make.sys ; \
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


Export_QE-5.3.0 : initialize

	@cd $(Export_QE-5.3.0_srcPath) ; \
		make all

TME : initialize

	@cd $(TME_srcPath) ; \
        	make all

Mj : initialize

	@cd $(Mj_srcPath) ; \
        	make all

LSF : LSF0 LSF1


LSF0 : initialize

	@cd $(LSF0_srcPath) ; \
        	make all

LSF1 : initialize

	@cd $(LSF1_srcPath) ; \
        	make all

Sigma : initialize

	@cd $(Sigma_srcPath) ; \
        	make all


clean_all_QE-5.3.0 : clean_QE-5.3.0_dependent cleanTME cleanInitialization cleanMj cleanLSF0 cleanLSF1 cleanSigma

clean_QE-5.3.0_dependent : cleanExport_QE-5.3.0

cleanInitialization :

	@if test -d bin ; then \
		/bin/rm -rf bin ; \
	fi
	@if test -e make.sys ; then \
		/bin/rm -f make.sys ; \
	fi

cleanExport_QE-5.3.0 :

	@cd $(Export_QE-5.3.0_srcPath) ; \
        	make clean

cleanTME :

	@cd $(TME_srcPath) ; \
		make clean

cleanMj :

	@cd $(Mj_srcPath) ; \
		make clean


cleanLSF : cleanLSF0 cleanLSF1

cleanLSF0 :

	@cd $(LSF0_srcPath) ; \
		make clean

cleanLSF1 :

	@cd $(LSF1_srcPath) ; \
		make clean

cleanSigma :

	@cd $(Sigma_srcPath) ; \
		make clean

