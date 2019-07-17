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
* `../TME/src/TME_Module_v28.f90:166:` Consider changing `atom` type to `element` since it holds more than one atom
* `../TME/src/TME_Module_v28.f90:406:` Change `checkInitialization()` to have arguments to make clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:548:` Remove everything with `ki` and `kf` because never used
* `../TME/src/TME_Module_v28.f90:833:` Look more into how AE and PS wavefunctions are combined to further understand this
* `../TME/src/TME_Module_v28.f90:834:` Move this behavior to another subroutine for clarity
* `../TME/src/TME_Module_v28.f90:843:` Figure out if differences in PC and SD `F1` calculations are intentional
* `../TME/src/TME_Module_v28.f90:844:` Figure out if should be `(wps_i wae_j - wae_i wps_j)r_{ab}`
* `../TME/src/TME_Module_v28.f90:883:` Figure out if intentional to only use `JMAX` from SD input
* `../TME/src/TME_Module_v28.f90:1039:` Document `calculatePWsOverlap()`
* `../TME/src/TME_Module_v28.f90:1161:` Document `readWfcSD()`
* `../TME/src/TME_Module_v28.f90:1162:` Figure out difference between PC and SD `readWfc` and possibly merge
* `../TME/src/TME_Module_v28.f90:1219:` Document `readProjectionsPC()`
* `../TME/src/TME_Module_v28.f90:1252:` Document `readProjectionsSD()`
* `../TME/src/TME_Module_v28.f90:1253:` Figure out the difference between PC and SD `readProjections` and possibly merge
* `../TME/src/TME_Module_v28.f90:1286:` Document `projectBetaPCwfcSD()`
* `../TME/src/TME_Module_v28.f90:1349:` Document `projectBetaSDwfcPC()`
* `../TME/src/TME_Module_v28.f90:1350:` Figure out the difference between PC and SD `projectBeta_wfc_` and possibly merge
* `../TME/src/TME_Module_v28.f90:1413:` Document `pawCorrectionPsiPC()`
* `../TME/src/TME_Module_v28.f90:1479:` Document `pawCorrectionSDPhi()`
* `../TME/src/TME_Module_v28.f90:1480:` Figure out the difference between PC and SD `pawCorrectionPsi` and possibly merge
* `../TME/src/TME_Module_v28.f90:1544:` Document `pawCorrectionKPC()`
* `../TME/src/TME_Module_v28.f90:1643:` Document `pawCorrectionSDK()`
* `../TME/src/TME_Module_v28.f90:1644:` Figure out the difference between PC and SD `pawCorrection_K` and possibly merge
* `../TME/src/TME_Module_v28.f90:1741:` Document `pawCorrection()`
* `../TME/src/TME_Module_v28.f90:1742:` Figure out difference between `pawCorrection()` and the PC `pawCorrection` functions
* `../TME/src/TME_Module_v28.f90:2091:` Document `bessel_j()`
* `../TME/src/TME_Module_v28.f90:2121:` Document `writeResults()` @endto
* `../TME/src/TME_Module_v28.f90:2183:` Document `readUfis()`
* `../TME/src/TME_Module_v28.f90:2237:` Document `calculateVFiElements()`
* `../TME/src/TME_Module_v28.f90:2359:` Document `readEigenvalues()`
