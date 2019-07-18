title: Todo
author: Laura Nichols
date: 07/18/2019

* `../TME/src/TME_Main_v9.f90:2:` Finish documentation for main program
* `../TME/src/TME_Main_v9.f90:46:` Figure out if need to allocate space for arrays so soon
* `../TME/src/TME_Main_v9.f90:58:` Figure out if SD and PC `numOfGvecs` should be the same
* `../TME/src/TME_Main_v9.f90:325:` Figure out if should be solid defect volume or pristine
* `../TME/src/TME_Main_v9.f90:326:` Are pristine and solid defect volume the same?
* `../TME/src/TME_Module_v28.f90:34:` Change I/O from file to console so that usage matches that of QE
* `../TME/src/TME_Module_v28.f90:164:` Consider changing `atom` type to `element` since it holds more than one atom
* `../TME/src/TME_Module_v28.f90:405:` Change `checkInitialization()` to have arguments to make clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:547:` Remove everything with `ki` and `kf` because never used
* `../TME/src/TME_Module_v28.f90:832:` Look more into how AE and PS wavefunctions are combined to further understand this
* `../TME/src/TME_Module_v28.f90:833:` Move this behavior to another subroutine for clarity
* `../TME/src/TME_Module_v28.f90:842:` Figure out if differences in PC and SD `F1` calculations are intentional
* `../TME/src/TME_Module_v28.f90:843:` Figure out if should be `(wps_i wae_j - wae_i wps_j)r_{ab}`
* `../TME/src/TME_Module_v28.f90:882:` Figure out if intentional to only use `JMAX` from SD input
* `../TME/src/TME_Module_v28.f90:1038:` Document `calculatePWsOverlap()`
* `../TME/src/TME_Module_v28.f90:1065:` Figure out what `Ufi` is supposed to be
* `../TME/src/TME_Module_v28.f90:1182:` Document `readProjectionsPC()`
* `../TME/src/TME_Module_v28.f90:1215:` Document `readProjectionsSD()`
* `../TME/src/TME_Module_v28.f90:1216:` Figure out the difference between PC and SD `readProjections` and possibly merge
* `../TME/src/TME_Module_v28.f90:1371:` Document `pawCorrectionPsiPC()`
* `../TME/src/TME_Module_v28.f90:1437:` Document `pawCorrectionSDPhi()`
* `../TME/src/TME_Module_v28.f90:1438:` Figure out the difference between PC and SD `pawCorrectionPsi` and possibly merge
* `../TME/src/TME_Module_v28.f90:1502:` Document `pawCorrectionKPC()`
* `../TME/src/TME_Module_v28.f90:1601:` Document `pawCorrectionSDK()`
* `../TME/src/TME_Module_v28.f90:1602:` Figure out the difference between PC and SD `pawCorrection_K` and possibly merge
* `../TME/src/TME_Module_v28.f90:1699:` Document `pawCorrection()`
* `../TME/src/TME_Module_v28.f90:1700:` Figure out difference between `pawCorrection()` and the PC `pawCorrection` functions
* `../TME/src/TME_Module_v28.f90:2049:` Document `bessel_j()`
* `../TME/src/TME_Module_v28.f90:2079:` Document `writeResults()` @endto
* `../TME/src/TME_Module_v28.f90:2141:` Document `readUfis()`
* `../TME/src/TME_Module_v28.f90:2195:` Document `calculateVFiElements()`
* `../TME/src/TME_Module_v28.f90:2317:` Document `readEigenvalues()`
