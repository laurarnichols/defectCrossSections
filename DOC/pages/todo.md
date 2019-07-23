title: Todo
author: Laura Nichols
date: 07/23/2019

* Binary file ../TME/src/.TME_Module_v28.f90.swp matches
* `../TME/src/TME_Main_v9.f90:2:` Finish documentation for main program
* `../TME/src/TME_Main_v9.f90:46:` Figure out if need to allocate space for arrays so soon
* `../TME/src/TME_Main_v9.f90:58:` Figure out if SD and PC `numOfGvecs` should be the same
* `../TME/src/TME_Main_v9.f90:329:` Figure out if should be solid defect volume or pristine
* `../TME/src/TME_Main_v9.f90:330:` Are pristine and solid defect volume the same?
* `../TME/src/TME_Module_v28.f90:34:` Change I/O from file to console so that usage matches that of QE
* `../TME/src/TME_Module_v28.f90:158:` Consider changing `atom` type to `element` since it holds more than one atom
* `../TME/src/TME_Module_v28.f90:403:` Change `checkInitialization()` to have arguments to make clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:545:` Remove everything with `ki` and `kf` because never used
* `../TME/src/TME_Module_v28.f90:832:` Look more into how AE and PS wavefunctions are combined to further understand this
* `../TME/src/TME_Module_v28.f90:833:` Move this behavior to another subroutine for clarity
* `../TME/src/TME_Module_v28.f90:842:` Figure out if differences in PC and SD `F1` calculations are intentional
* `../TME/src/TME_Module_v28.f90:843:` Figure out if should be `(wae_i wae_j - wps_i wps_j)r_{ab}`
* `../TME/src/TME_Module_v28.f90:844:` Figure out if first term in each should be conjugated for inner product form
* `../TME/src/TME_Module_v28.f90:845:` Figure out if `rab` plays role of \(dr\) within augmentation sphere
* `../TME/src/TME_Module_v28.f90:884:` Figure out if intentional to only use `JMAX` from SD input
* `../TME/src/TME_Module_v28.f90:1040:` Document `calculatePWsOverlap()`
* `../TME/src/TME_Module_v28.f90:1067:` Figure out what `Ufi` is supposed to be
* `../TME/src/TME_Module_v28.f90:1218:` Get actual perfect crystal and solid defect output to test
* `../TME/src/TME_Module_v28.f90:1219:` Figure out if loop should be over `solidDefect` or `perfectCrystal`
* `../TME/src/TME_Module_v28.f90:1220:` Look into `nSpins` to figure out if it is needed
* `../TME/src/TME_Module_v28.f90:1412:` Figure out the significance of \(l = l^{\prime}\) and \(m = m^{\prime}\)
* `../TME/src/TME_Module_v28.f90:1478:` Document `pawCorrectionKPC()`
* `../TME/src/TME_Module_v28.f90:1526:` Figure out if this output slows things down significantly
* `../TME/src/TME_Module_v28.f90:1527:` Figure out if formula gives accurate representation of time left
* `../TME/src/TME_Module_v28.f90:1614:` Document `pawCorrectionSDK()`
* `../TME/src/TME_Module_v28.f90:1615:` Figure out the difference between PC and SD `pawCorrection_K` and possibly merge
* `../TME/src/TME_Module_v28.f90:1712:` Document `pawCorrection()`
* `../TME/src/TME_Module_v28.f90:1713:` Figure out difference between `pawCorrection()` and the PC `pawCorrection` functions
* `../TME/src/TME_Module_v28.f90:2095:` Document `writeResults()` @endto
* `../TME/src/TME_Module_v28.f90:2157:` Document `readUfis()`
* `../TME/src/TME_Module_v28.f90:2211:` Document `calculateVFiElements()`
* `../TME/src/TME_Module_v28.f90:2333:` Document `readEigenvalues()`
