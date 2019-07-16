title: Todo
author: Laura Nichols
date: 07/16/2019

* Binary file ../TME/src/.TME_Module_v28.f90.swp matches
* `../TME/src/TME_Main_v9.f90:2:` Finish documentation for main program
* `../TME/src/TME_Main_v9.f90:46:` Figure out if need to allocate space for arrays so soon
* `../TME/src/TME_Main_v9.f90:58:` Figure out if SD and PC `numOfGvecs` should be the same
* `../TME/src/TME_Main_v9.f90:324:` Figure out if should be solid defect volume or pristine
* `../TME/src/TME_Main_v9.f90:325:` Are pristine and solid defect volume the same?
* `../TME/src/TME_Module_v28.f90:34:` Change I/O from file to console so that usage matches that of QE
* `../TME/src/TME_Module_v28.f90:170:` Consider changing `atom` type to `element` since it holds more than one atom
* `../TME/src/TME_Module_v28.f90:409:` Change `checkInitialization()` to have arguments to make clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:551:` Remove everything with `ki` and `kf` because never used
* `../TME/src/TME_Module_v28.f90:836:` Look more into how AE and PS wavefunctions are combined to further understand this
* `../TME/src/TME_Module_v28.f90:837:` Move this behavior to another subroutine for clarity
* `../TME/src/TME_Module_v28.f90:846:` Figure out if differences in PC and SD `F1` calculations are intentional
* `../TME/src/TME_Module_v28.f90:847:` Figure out if should be `(wps_i wae_j - wae_i wps_j)r_{ab}`
* `../TME/src/TME_Module_v28.f90:886:` Figure out if intentional to only use `JMAX` from SD input
* `../TME/src/TME_Module_v28.f90:1027:` Change if statement to use `int2str` subroutine
* `../TME/src/TME_Module_v28.f90:1050:` Document `calculatePWsOverlap()`
* `../TME/src/TME_Module_v28.f90:1079:` Document `readWfcPC()`
* `../TME/src/TME_Module_v28.f90:1135:` Document `readWfcSD()`
* `../TME/src/TME_Module_v28.f90:1136:` Figure out difference between PC and SD `readWfc` and possibly merge
* `../TME/src/TME_Module_v28.f90:1192:` Document `readProjectionsPC()`
* `../TME/src/TME_Module_v28.f90:1225:` Document `readProjectionsSD()`
* `../TME/src/TME_Module_v28.f90:1226:` Figure out the difference between PC and SD `readProjections` and possibly merge
* `../TME/src/TME_Module_v28.f90:1259:` Document `projectBetaPCwfcSD()`
* `../TME/src/TME_Module_v28.f90:1321:` Document `projectBetaSDwfcPC()`
* `../TME/src/TME_Module_v28.f90:1322:` Figure out the difference between PC and SD `projectBeta_wfc_` and possibly merge
* `../TME/src/TME_Module_v28.f90:1384:` Document `pawCorrectionPsiPC()`
* `../TME/src/TME_Module_v28.f90:1450:` Document `pawCorrectionSDPhi()`
* `../TME/src/TME_Module_v28.f90:1451:` Figure out the difference between PC and SD `pawCorrectionPsi` and possibly merge
* `../TME/src/TME_Module_v28.f90:1515:` Document `pawCorrectionKPC()`
* `../TME/src/TME_Module_v28.f90:1614:` Document `pawCorrectionSDK()`
* `../TME/src/TME_Module_v28.f90:1615:` Figure out the difference between PC and SD `pawCorrection_K` and possibly merge
* `../TME/src/TME_Module_v28.f90:1712:` Document `pawCorrection()`
* `../TME/src/TME_Module_v28.f90:1713:` Figure out difference between `pawCorrection()` and the PC `pawCorrection` functions
* `../TME/src/TME_Module_v28.f90:2062:` Document `bessel_j()`
* `../TME/src/TME_Module_v28.f90:2092:` Document `writeResults()` @endto
* `../TME/src/TME_Module_v28.f90:2154:` Document `readUfis()`
* `../TME/src/TME_Module_v28.f90:2208:` Document `calculateVFiElements()`
* `../TME/src/TME_Module_v28.f90:2330:` Document `readEigenvalues()`
