title: Todo
author: Laura Nichols
date: 07/08/2019

* `../TME/src/TME_Module_v28.f90:191:` Consider changing `atom` type to `element` since it holds more than one atom
* `../TME/src/TME_Module_v28.f90:231:` Use crystal type instead of all of the explicit variables
* `../TME/src/TME_Module_v28.f90:376:` Change `readInput()` to have arguments to make clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:388:` Figure out what the difference in PC and SD is
* `../TME/src/TME_Module_v28.f90:411:` Change `checkInitialization()` to have arguments to make clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:419:` Change `file_exists` to `fileExists` in `checkInitialization()`
* `../TME/src/TME_Module_v28.f90:605:` Check if there is any kind of check on `ki` and `kf`. Why was this commented out?
* `../TME/src/TME_Module_v28.f90:675:` Change `readInputPC()` to have arguments so that it is clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:758:` Add information about these variables to top
* `../TME/src/TME_Module_v28.f90:885: Change `(system%posIon(j,ni) , j = 1,3)` to `system%posIon(1:`3,ni)` in `readQEExport()` for clarity
* `../TME/src/TME_Module_v28.f90:975:` Look more into how AE and PS wavefunctions are combined to further understand this
* `../TME/src/TME_Module_v28.f90:976:` Move this behavior to another subroutine for clarity
* `../TME/src/TME_Module_v28.f90:985:` Figure out if differences in PC and SD `F1` calculations are intentional
* `../TME/src/TME_Module_v28.f90:986:` Figure out if should be `(wps_i wae_j - wae_i wps_j)r_{ab}`
* `../TME/src/TME_Module_v28.f90:1023:` Figure out if intentional to only use `JMAX` from SD input
