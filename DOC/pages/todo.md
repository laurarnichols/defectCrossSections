title: Todo
author: Laura Nichols
date: 07/09/2019

* `../TME/src/TME_Module_v28.f90:34:` Change I/O from file to console so that usage matches that of QE
* `../TME/src/TME_Module_v28.f90:192:` Consider changing `atom` type to `element` since it holds more than one atom
* `../TME/src/TME_Module_v28.f90:232:` Use crystal type instead of all of the explicit variables
* `../TME/src/TME_Module_v28.f90:377:` Change `readInput()` to have arguments to make clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:389:` Figure out what the difference in PC and SD is
* `../TME/src/TME_Module_v28.f90:412:` Change `checkInitialization()` to have arguments to make clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:420:` Change `file_exists` to `fileExists` in `checkInitialization()`
* `../TME/src/TME_Module_v28.f90:606:` Check if there is any kind of check on `ki` and `kf`. Why was this commented out?
* `../TME/src/TME_Module_v28.f90:676:` Change `readInputPC()` to have arguments so that it is clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:759:` Add information about these variables to top
* `../TME/src/TME_Module_v28.f90:886: Change `(system%posIon(j,ni) , j = 1,3)` to `system%posIon(1:`3,ni)` in `readQEExport()` for clarity
* `../TME/src/TME_Module_v28.f90:976:` Look more into how AE and PS wavefunctions are combined to further understand this
* `../TME/src/TME_Module_v28.f90:977:` Move this behavior to another subroutine for clarity
* `../TME/src/TME_Module_v28.f90:986:` Figure out if differences in PC and SD `F1` calculations are intentional
* `../TME/src/TME_Module_v28.f90:987:` Figure out if should be `(wps_i wae_j - wae_i wps_j)r_{ab}`
* `../TME/src/TME_Module_v28.f90:1024:` Figure out if intentional to only use `JMAX` from SD input
