title: Todo
author: Laura Nichols
date: 07/09/2019

* `../TME/src/TME_Main_v9.f90:317:` Figure out if should be solid defect volume or pristine
* `../TME/src/TME_Main_v9.f90:318:` Are pristine and solid defect volume the same?
* `../TME/src/TME_Module_v28.f90:34:` Change I/O from file to console so that usage matches that of QE
* `../TME/src/TME_Module_v28.f90:190:` Consider changing `atom` type to `element` since it holds more than one atom
* `../TME/src/TME_Module_v28.f90:230:` Use crystal type instead of all of the explicit variables
* `../TME/src/TME_Module_v28.f90:375:` Change `readInput()` to have arguments to make clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:387:` Figure out what the difference in PC and SD is
* `../TME/src/TME_Module_v28.f90:410:` Change `checkInitialization()` to have arguments to make clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:418:` Change `file_exists` to `fileExists` in `checkInitialization()`
* `../TME/src/TME_Module_v28.f90:604:` Check if there is any kind of check on `ki` and `kf`. Why was this commented out?
* `../TME/src/TME_Module_v28.f90:674:` Change `readInputPC()` to have arguments so that it is clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:757:` Add information about these variables to top
* `../TME/src/TME_Module_v28.f90:863: Change `(system%posIon(j,ni) , j = 1,3)` to `system%posIon(1:`3,ni)` in `readQEExport()` for clarity
* `../TME/src/TME_Module_v28.f90:953:` Look more into how AE and PS wavefunctions are combined to further understand this
* `../TME/src/TME_Module_v28.f90:954:` Move this behavior to another subroutine for clarity
* `../TME/src/TME_Module_v28.f90:963:` Figure out if differences in PC and SD `F1` calculations are intentional
* `../TME/src/TME_Module_v28.f90:964:` Figure out if should be `(wps_i wae_j - wae_i wps_j)r_{ab}`
* `../TME/src/TME_Module_v28.f90:1003:` Figure out if intentional to only use `JMAX` from SD input
