title: Todo
author: Laura Nichols
date: 07/23/2019

* `../TME/src/PC.txt:2:` Document `pawCorrectionKPC()`
* `../TME/src/SD.txt:2:` Document `pawCorrectionSDK()`
* `../TME/src/SD.txt:3:` Figure out the difference between PC and SD `pawCorrection_K` and possibly merge
* `../TME/src/TME_Main_v9.f90:2:` Finish documentation for main program
* `../TME/src/TME_Main_v9.f90:46:` Figure out if need to allocate space for arrays so soon
* `../TME/src/TME_Main_v9.f90:58:` Figure out if SD and PC `numOfGvecs` should be the same
* `../TME/src/TME_Main_v9.f90:329:` Figure out if should be solid defect volume or pristine
* `../TME/src/TME_Main_v9.f90:330:` Are pristine and solid defect volume the same?
* `../TME/src/TME_Module_v28.f90:34:` Change I/O from file to console so that usage matches that of QE
* `../TME/src/TME_Module_v28.f90:156:` Consider changing `atom` type to `element` since it holds more than one atom
* `../TME/src/TME_Module_v28.f90:402:` Change `checkInitialization()` to have arguments to make clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:544:` Remove everything with `ki` and `kf` because never used
* `../TME/src/TME_Module_v28.f90:831:` Look more into how AE and PS wavefunctions are combined to further understand this
* `../TME/src/TME_Module_v28.f90:832:` Move this behavior to another subroutine for clarity
* `../TME/src/TME_Module_v28.f90:841:` Figure out if differences in PC and SD `F1` calculations are intentional
* `../TME/src/TME_Module_v28.f90:842:` Figure out if should be `(wae_i wae_j - wps_i wps_j)r_{ab}`
* `../TME/src/TME_Module_v28.f90:843:` Figure out if first term in each should be conjugated for inner product form
* `../TME/src/TME_Module_v28.f90:844:` Figure out if `rab` plays role of \(dr\) within augmentation sphere
* `../TME/src/TME_Module_v28.f90:883:` Figure out if intentional to only use `JMAX` from SD input
* `../TME/src/TME_Module_v28.f90:1039:` Document `calculatePWsOverlap()`
* `../TME/src/TME_Module_v28.f90:1066:` Figure out what `Ufi` is supposed to be
* `../TME/src/TME_Module_v28.f90:1217:` Get actual perfect crystal and solid defect output to test
* `../TME/src/TME_Module_v28.f90:1218:` Figure out if loop should be over `solidDefect` or `perfectCrystal`
* `../TME/src/TME_Module_v28.f90:1219:` Look into `nSpins` to figure out if it is needed
* `../TME/src/TME_Module_v28.f90:1237:` Figure out what this subroutine really does
* `../TME/src/TME_Module_v28.f90:1359:` Figure out what this subroutine really does
* `../TME/src/TME_Module_v28.f90:1411:` Figure out the significance of \(l = l^{\prime}\) and \(m = m^{\prime}\)
* `../TME/src/TME_Module_v28.f90:1477:` Figure out what this subroutine really does
* `../TME/src/TME_Module_v28.f90:1545:` Figure out if this output slows things down significantly
* `../TME/src/TME_Module_v28.f90:1546:` Figure out if formula gives accurate representation of time left
* `../TME/src/TME_Module_v28.f90:1590:` Figure out if this should be perfect crystal
* `../TME/src/TME_Module_v28.f90:1591:` Figure out significance of "qr" point
* `../TME/src/TME_Module_v28.f90:1594:` Test if can just directly store in each atom type's `bes_J_qr`
* `../TME/src/TME_Module_v28.f90:1614:` Figure out if this should be `gDotR`
* `../TME/src/TME_Module_v28.f90:1617:` Figure out why this is called `ATOMIC_CENTER`
* `../TME/src/TME_Module_v28.f90:1664:` Document `pawCorrectionSDK()`
* `../TME/src/TME_Module_v28.f90:1665:` Figure out the difference between PC and SD `pawCorrection_K` and possibly merge
* `../TME/src/TME_Module_v28.f90:1763:` Document `pawCorrection()`
* `../TME/src/TME_Module_v28.f90:1764:` Figure out difference between `pawCorrection()` and the PC `pawCorrection` functions
* `../TME/src/TME_Module_v28.f90:2127:` Document `writeResults()` @endto
* `../TME/src/TME_Module_v28.f90:2189:` Document `readUfis()`
* `../TME/src/TME_Module_v28.f90:2243:` Document `calculateVFiElements()`
* `../TME/src/TME_Module_v28.f90:2365:` Document `readEigenvalues()`
