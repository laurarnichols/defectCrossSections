title: Todo
author: Laura Nichols
date: 07/22/2019

* `../TME/src/TME_Main_v9.f90:2:` Finish documentation for main program
* `../TME/src/TME_Main_v9.f90:46:` Figure out if need to allocate space for arrays so soon
* `../TME/src/TME_Main_v9.f90:58:` Figure out if SD and PC `numOfGvecs` should be the same
* `../TME/src/TME_Main_v9.f90:328:` Figure out if should be solid defect volume or pristine
* `../TME/src/TME_Main_v9.f90:329:` Are pristine and solid defect volume the same?
* `../TME/src/TME_Module_v28.f90:34:` Change I/O from file to console so that usage matches that of QE
* `../TME/src/TME_Module_v28.f90:158:` Consider changing `atom` type to `element` since it holds more than one atom
* `../TME/src/TME_Module_v28.f90:401:` Change `checkInitialization()` to have arguments to make clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:543:` Remove everything with `ki` and `kf` because never used
* `../TME/src/TME_Module_v28.f90:830:` Look more into how AE and PS wavefunctions are combined to further understand this
* `../TME/src/TME_Module_v28.f90:831:` Move this behavior to another subroutine for clarity
* `../TME/src/TME_Module_v28.f90:840:` Figure out if differences in PC and SD `F1` calculations are intentional
* `../TME/src/TME_Module_v28.f90:841:` Figure out if should be `(wae_i wae_j - wps_i wps_j)r_{ab}`
* `../TME/src/TME_Module_v28.f90:842:` Figure out if first term in each should be conjugated for inner product form
* `../TME/src/TME_Module_v28.f90:843:` Figure out if `rab` plays role of \(dr\) within augmentation sphere
* `../TME/src/TME_Module_v28.f90:882:` Figure out if intentional to only use `JMAX` from SD input
* `../TME/src/TME_Module_v28.f90:1038:` Document `calculatePWsOverlap()`
* `../TME/src/TME_Module_v28.f90:1065:` Figure out what `Ufi` is supposed to be
* `../TME/src/TME_Module_v28.f90:1216:` Get actual perfect crystal and solid defect output to test
* `../TME/src/TME_Module_v28.f90:1217:` Figure out if loop should be over `solidDefect` or `perfectCrystal`
* `../TME/src/TME_Module_v28.f90:1218:` Look into `nSpins` to figure out if it is needed
* `../TME/src/TME_Module_v28.f90:1410:` Figure out the significance of \(l = l^{\prime}\) and \(m = m^{\prime}\)
* `../TME/src/TME_Module_v28.f90:1476:` Document `pawCorrectionKPC()`
* `../TME/src/TME_Module_v28.f90:1575:` Document `pawCorrectionSDK()`
* `../TME/src/TME_Module_v28.f90:1576:` Figure out the difference between PC and SD `pawCorrection_K` and possibly merge
* `../TME/src/TME_Module_v28.f90:1673:` Document `pawCorrection()`
* `../TME/src/TME_Module_v28.f90:1674:` Figure out difference between `pawCorrection()` and the PC `pawCorrection` functions
* `../TME/src/TME_Module_v28.f90:2070:` Document `writeResults()` @endto
* `../TME/src/TME_Module_v28.f90:2132:` Document `readUfis()`
* `../TME/src/TME_Module_v28.f90:2186:` Document `calculateVFiElements()`
* `../TME/src/TME_Module_v28.f90:2308:` Document `readEigenvalues()`
