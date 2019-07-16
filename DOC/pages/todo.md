title: Todo
author: Laura Nichols
date: 07/16/2019

* `../TME/src/TME_Main_v9.f90:2:` Finish documentation for main program
* `../TME/src/TME_Main_v9.f90:319:` Figure out if should be solid defect volume or pristine
* `../TME/src/TME_Main_v9.f90:320:` Are pristine and solid defect volume the same?
* `../TME/src/TME_Module_v28.f90:34:` Change I/O from file to console so that usage matches that of QE
* `../TME/src/TME_Module_v28.f90:170:` Consider changing `atom` type to `element` since it holds more than one atom
* `../TME/src/TME_Module_v28.f90:409:` Change `checkInitialization()` to have arguments to make clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:551:` Remove everything with `ki` and `kf` because never used
* `../TME/src/TME_Module_v28.f90:763: Change `(system%posIon(j,ni) , j = 1,3)` to `system%posIon(1:`3,ni)` in `readQEExport()` for clarity
* `../TME/src/TME_Module_v28.f90:837:` Look more into how AE and PS wavefunctions are combined to further understand this
* `../TME/src/TME_Module_v28.f90:838:` Move this behavior to another subroutine for clarity
* `../TME/src/TME_Module_v28.f90:847:` Figure out if differences in PC and SD `F1` calculations are intentional
* `../TME/src/TME_Module_v28.f90:848:` Figure out if should be `(wps_i wae_j - wae_i wps_j)r_{ab}`
* `../TME/src/TME_Module_v28.f90:887:` Figure out if intentional to only use `JMAX` from SD input
* `../TME/src/TME_Module_v28.f90:922:` Document `readPWsSet()`
* `../TME/src/TME_Module_v28.f90:952:` Document `distributePWstoProcs()`
* `../TME/src/TME_Module_v28.f90:982:` Document `checkIfCalculated()`
* `../TME/src/TME_Module_v28.f90:991:` Change if statement to use `int2str` subroutine
* `../TME/src/TME_Module_v28.f90:1012:` Document `calculatePWsOverlap()`
* `../TME/src/TME_Module_v28.f90:1041:` Document `readWfcPC()`
* `../TME/src/TME_Module_v28.f90:1097:` Document `readWfcSD()`
* `../TME/src/TME_Module_v28.f90:1098:` Figure out difference between PC and SD `readWfc` and possibly merge
* `../TME/src/TME_Module_v28.f90:1154:` Document `readProjectionsPC()`
* `../TME/src/TME_Module_v28.f90:1187:` Document `readProjectionsSD()`
* `../TME/src/TME_Module_v28.f90:1188:` Figure out the difference between PC and SD `readProjections` and possibly merge
* `../TME/src/TME_Module_v28.f90:1221:` Document `projectBetaPCwfcSD()`
* `../TME/src/TME_Module_v28.f90:1283:` Document `projectBetaSDwfcPC()`
* `../TME/src/TME_Module_v28.f90:1284:` Figure out the difference between PC and SD `projectBeta_wfc_` and possibly merge
* `../TME/src/TME_Module_v28.f90:1346:` Document `pawCorrectionPsiPC()`
* `../TME/src/TME_Module_v28.f90:1412:` Document `pawCorrectionSDPhi()`
* `../TME/src/TME_Module_v28.f90:1413:` Figure out the difference between PC and SD `pawCorrectionPsi` and possibly merge
* `../TME/src/TME_Module_v28.f90:1477:` Document `pawCorrectionKPC()`
* `../TME/src/TME_Module_v28.f90:1576:` Document `pawCorrectionSDK()`
* `../TME/src/TME_Module_v28.f90:1577:` Figure out the difference between PC and SD `pawCorrection_K` and possibly merge
* `../TME/src/TME_Module_v28.f90:1674:` Document `pawCorrection()`
* `../TME/src/TME_Module_v28.f90:1675:` Figure out difference between `pawCorrection()` and the PC `pawCorrection` functions
* `../TME/src/TME_Module_v28.f90:2024:` Document `bessel_j()`
* `../TME/src/TME_Module_v28.f90:2054:` Document `writeResults()` @endto
* `../TME/src/TME_Module_v28.f90:2116:` Document `readUfis()`
* `../TME/src/TME_Module_v28.f90:2170:` Document `calculateVFiElements()`
* `../TME/src/TME_Module_v28.f90:2292:` Document `readEigenvalues()`
* `../TME/src/TME_Module_v28.f90:2429:` Document `int2str()`
