title: Todo
author: Laura Nichols
date: 07/15/2019

* `../TME/src/TME_Main_v9.f90:2:` Finish documentation for main program
* `../TME/src/TME_Main_v9.f90:319:` Figure out if should be solid defect volume or pristine
* `../TME/src/TME_Main_v9.f90:320:` Are pristine and solid defect volume the same?
* `../TME/src/TME_Module_v28.f90:34:` Change I/O from file to console so that usage matches that of QE
* `../TME/src/TME_Module_v28.f90:170:` Consider changing `atom` type to `element` since it holds more than one atom
* `../TME/src/TME_Module_v28.f90:409:` Change `checkInitialization()` to have arguments to make clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:551:` Remove everything with `ki` and `kf` because never used
* `../TME/src/TME_Module_v28.f90:833: Change `(system%posIon(j,ni) , j = 1,3)` to `system%posIon(1:`3,ni)` in `readQEExport()` for clarity
* `../TME/src/TME_Module_v28.f90:907:` Look more into how AE and PS wavefunctions are combined to further understand this
* `../TME/src/TME_Module_v28.f90:908:` Move this behavior to another subroutine for clarity
* `../TME/src/TME_Module_v28.f90:917:` Figure out if differences in PC and SD `F1` calculations are intentional
* `../TME/src/TME_Module_v28.f90:918:` Figure out if should be `(wps_i wae_j - wae_i wps_j)r_{ab}`
* `../TME/src/TME_Module_v28.f90:957:` Figure out if intentional to only use `JMAX` from SD input
* `../TME/src/TME_Module_v28.f90:992:` Document `distributePWstoProcs()`
* `../TME/src/TME_Module_v28.f90:1022:` Document `int2str()`
* `../TME/src/TME_Module_v28.f90:1063:` Document `readPWsSet()`
* `../TME/src/TME_Module_v28.f90:1093:` Document `readWfcPC()`
* `../TME/src/TME_Module_v28.f90:1149:` Document `projectBetaPCwfcSD()`
* `../TME/src/TME_Module_v28.f90:1211:` Document `readWfcSD()`
* `../TME/src/TME_Module_v28.f90:1212:` Figure out difference between PC and SD `readWfc` and possibly merge
* `../TME/src/TME_Module_v28.f90:1268:` Document `calculatePWsOverlap()`
* `../TME/src/TME_Module_v28.f90:1297:` Document `readProjectionsPC()`
* `../TME/src/TME_Module_v28.f90:1330:` Document `readProjectionsSD()`
* `../TME/src/TME_Module_v28.f90:1331:` Figure out the difference between PC and SD `readProjections` and possibly merge
* `../TME/src/TME_Module_v28.f90:1364:` Document `projectBetaSDwfcPC()`
* `../TME/src/TME_Module_v28.f90:1365:` Figure out the difference between PC and SD `projectBeta_wfc_` and possibly merge
* `../TME/src/TME_Module_v28.f90:1427:` Document `pawCorrection()`
* `../TME/src/TME_Module_v28.f90:1526:` Document `pawCorrectionSDK()`
* `../TME/src/TME_Module_v28.f90:1527:` Figure out the difference between PC and SD `pawCorrection_K` and possibly merge
* `../TME/src/TME_Module_v28.f90:1624:` Document `pawCorrectionPsiPC()`
* `../TME/src/TME_Module_v28.f90:1690:` Document `pawCorrectionSDPhi()`
* `../TME/src/TME_Module_v28.f90:1691:` Figure out the difference between PC and SD `pawCorrectionPsi` and possibly merge
* `../TME/src/TME_Module_v28.f90:1755:` Document `pawCorrection()`
* `../TME/src/TME_Module_v28.f90:1756:` Figure out difference between `pawCorrection()` and the PC `pawCorrection` functions
* `../TME/src/TME_Module_v28.f90:1829:` Document `readEigenvalues()`
* `../TME/src/TME_Module_v28.f90:1876:` Document `calculateVFiElements()`
* `../TME/src/TME_Module_v28.f90:1998:` Document `checkIfCalculated()`
* `../TME/src/TME_Module_v28.f90:2007:` Change if statement to use `int2str` subroutine
* `../TME/src/TME_Module_v28.f90:2028:` Document `readUfis()`
* `../TME/src/TME_Module_v28.f90:2082:` Document `writeResults()` @endto
* `../TME/src/TME_Module_v28.f90:2144:` Document `bessel_j()`
