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
* `../TME/src/TME_Module_v28.f90:1042:` Document `calculatePWsOverlap()`
* `../TME/src/TME_Module_v28.f90:1071:` Document `readWfcPC()`
* `../TME/src/TME_Module_v28.f90:1127:` Document `readWfcSD()`
* `../TME/src/TME_Module_v28.f90:1128:` Figure out difference between PC and SD `readWfc` and possibly merge
* `../TME/src/TME_Module_v28.f90:1184:` Document `readProjectionsPC()`
* `../TME/src/TME_Module_v28.f90:1217:` Document `readProjectionsSD()`
* `../TME/src/TME_Module_v28.f90:1218:` Figure out the difference between PC and SD `readProjections` and possibly merge
* `../TME/src/TME_Module_v28.f90:1251:` Document `projectBetaPCwfcSD()`
* `../TME/src/TME_Module_v28.f90:1313:` Document `projectBetaSDwfcPC()`
* `../TME/src/TME_Module_v28.f90:1314:` Figure out the difference between PC and SD `projectBeta_wfc_` and possibly merge
* `../TME/src/TME_Module_v28.f90:1376:` Document `pawCorrectionPsiPC()`
* `../TME/src/TME_Module_v28.f90:1442:` Document `pawCorrectionSDPhi()`
* `../TME/src/TME_Module_v28.f90:1443:` Figure out the difference between PC and SD `pawCorrectionPsi` and possibly merge
* `../TME/src/TME_Module_v28.f90:1507:` Document `pawCorrectionKPC()`
* `../TME/src/TME_Module_v28.f90:1606:` Document `pawCorrectionSDK()`
* `../TME/src/TME_Module_v28.f90:1607:` Figure out the difference between PC and SD `pawCorrection_K` and possibly merge
* `../TME/src/TME_Module_v28.f90:1704:` Document `pawCorrection()`
* `../TME/src/TME_Module_v28.f90:1705:` Figure out difference between `pawCorrection()` and the PC `pawCorrection` functions
* `../TME/src/TME_Module_v28.f90:2054:` Document `bessel_j()`
* `../TME/src/TME_Module_v28.f90:2084:` Document `writeResults()` @endto
* `../TME/src/TME_Module_v28.f90:2146:` Document `readUfis()`
* `../TME/src/TME_Module_v28.f90:2200:` Document `calculateVFiElements()`
* `../TME/src/TME_Module_v28.f90:2322:` Document `readEigenvalues()`
