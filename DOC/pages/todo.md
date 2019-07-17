title: Todo
author: Laura Nichols
date: 07/17/2019

* Binary file ../TME/src/.TME_Module_v28.f90.swo matches
* `../TME/src/TME_Main_v9.f90:2:` Finish documentation for main program
* `../TME/src/TME_Main_v9.f90:46:` Figure out if need to allocate space for arrays so soon
* `../TME/src/TME_Main_v9.f90:58:` Figure out if SD and PC `numOfGvecs` should be the same
* `../TME/src/TME_Main_v9.f90:325:` Figure out if should be solid defect volume or pristine
* `../TME/src/TME_Main_v9.f90:326:` Are pristine and solid defect volume the same?
* `../TME/src/TME_Module_v28.f90:34:` Change I/O from file to console so that usage matches that of QE
* `../TME/src/TME_Module_v28.f90:168:` Consider changing `atom` type to `element` since it holds more than one atom
* `../TME/src/TME_Module_v28.f90:408:` Change `checkInitialization()` to have arguments to make clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:550:` Remove everything with `ki` and `kf` because never used
* `../TME/src/TME_Module_v28.f90:835:` Look more into how AE and PS wavefunctions are combined to further understand this
* `../TME/src/TME_Module_v28.f90:836:` Move this behavior to another subroutine for clarity
* `../TME/src/TME_Module_v28.f90:845:` Figure out if differences in PC and SD `F1` calculations are intentional
* `../TME/src/TME_Module_v28.f90:846:` Figure out if should be `(wps_i wae_j - wae_i wps_j)r_{ab}`
* `../TME/src/TME_Module_v28.f90:885:` Figure out if intentional to only use `JMAX` from SD input
* `../TME/src/TME_Module_v28.f90:1041:` Document `calculatePWsOverlap()`
* `../TME/src/TME_Module_v28.f90:1160:` Document `readWfcSD()`
* `../TME/src/TME_Module_v28.f90:1161:` Figure out difference between PC and SD `readWfc` and possibly merge
* `../TME/src/TME_Module_v28.f90:1217:` Document `readProjectionsPC()`
* `../TME/src/TME_Module_v28.f90:1250:` Document `readProjectionsSD()`
* `../TME/src/TME_Module_v28.f90:1251:` Figure out the difference between PC and SD `readProjections` and possibly merge
* `../TME/src/TME_Module_v28.f90:1284:` Document `projectBetaPCwfcSD()`
* `../TME/src/TME_Module_v28.f90:1346:` Document `projectBetaSDwfcPC()`
* `../TME/src/TME_Module_v28.f90:1347:` Figure out the difference between PC and SD `projectBeta_wfc_` and possibly merge
* `../TME/src/TME_Module_v28.f90:1409:` Document `pawCorrectionPsiPC()`
* `../TME/src/TME_Module_v28.f90:1475:` Document `pawCorrectionSDPhi()`
* `../TME/src/TME_Module_v28.f90:1476:` Figure out the difference between PC and SD `pawCorrectionPsi` and possibly merge
* `../TME/src/TME_Module_v28.f90:1540:` Document `pawCorrectionKPC()`
* `../TME/src/TME_Module_v28.f90:1639:` Document `pawCorrectionSDK()`
* `../TME/src/TME_Module_v28.f90:1640:` Figure out the difference between PC and SD `pawCorrection_K` and possibly merge
* `../TME/src/TME_Module_v28.f90:1737:` Document `pawCorrection()`
* `../TME/src/TME_Module_v28.f90:1738:` Figure out difference between `pawCorrection()` and the PC `pawCorrection` functions
* `../TME/src/TME_Module_v28.f90:2087:` Document `bessel_j()`
* `../TME/src/TME_Module_v28.f90:2117:` Document `writeResults()` @endto
* `../TME/src/TME_Module_v28.f90:2179:` Document `readUfis()`
* `../TME/src/TME_Module_v28.f90:2233:` Document `calculateVFiElements()`
* `../TME/src/TME_Module_v28.f90:2355:` Document `readEigenvalues()`
