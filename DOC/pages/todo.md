title: Todo
author: Laura Nichols
date: 07/11/2019

* `../TME/src/TME_Main_v9.f90:318:` Figure out if should be solid defect volume or pristine
* `../TME/src/TME_Main_v9.f90:319:` Are pristine and solid defect volume the same?
* `../TME/src/TME_Module_v28.f90:34:` Change I/O from file to console so that usage matches that of QE
* `../TME/src/TME_Module_v28.f90:170:` Consider changing `atom` type to `element` since it holds more than one atom
* `../TME/src/TME_Module_v28.f90:409:` Change `checkInitialization()` to have arguments to make clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:577:` Remove everything with `ki` and `kf` because never used
* `../TME/src/TME_Module_v28.f90:714:` Change `readInputPC()` to have arguments so that it is clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:797:` Add information about these variables to top
* `../TME/src/TME_Module_v28.f90:865: Change `(system%posIon(j,ni) , j = 1,3)` to `system%posIon(1:`3,ni)` in `readQEExport()` for clarity
* `../TME/src/TME_Module_v28.f90:939:` Look more into how AE and PS wavefunctions are combined to further understand this
* `../TME/src/TME_Module_v28.f90:940:` Move this behavior to another subroutine for clarity
* `../TME/src/TME_Module_v28.f90:949:` Figure out if differences in PC and SD `F1` calculations are intentional
* `../TME/src/TME_Module_v28.f90:950:` Figure out if should be `(wps_i wae_j - wae_i wps_j)r_{ab}`
* `../TME/src/TME_Module_v28.f90:989:` Figure out if intentional to only use `JMAX` from SD input
