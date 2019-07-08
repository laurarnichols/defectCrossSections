title: Todo
author: Laura Nichols
date: 07/08/2019

* `../TME/src/TME_Module_v28.f90:217:` Use crystal type instead of all of the explicit variables
* `../TME/src/TME_Module_v28.f90:245:` Consider changing `atom` type to `element` since it holds more than one atom
* `../TME/src/TME_Module_v28.f90:374:` Change `readInput()` to have arguments to make clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:384:` Figure out what the difference in PC and SD is
* `../TME/src/TME_Module_v28.f90:409:` Change `checkInitialization()` to have arguments to make clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:417:` Change `file_exists` to `fileExists` in `checkInitialization()`
* `../TME/src/TME_Module_v28.f90:603:` Check if there is any kind of check on `ki` and `kf`. Why was this commented out?
* `../TME/src/TME_Module_v28.f90:674:` Change `readInputPC()` to have arguments so that it is clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:781:` Add information about these variables to top
* `../TME/src/TME_Module_v28.f90:908: Change `(posIon(j,ni) , j = 1,3)` to `posIon(1:`3,ni)` in `readQEExport()` for clarity
* `../TME/src/TME_Module_v28.f90:998:` Look more into how AE and PS wavefunctions are combined to further understand this
* `../TME/src/TME_Module_v28.f90:999:` Move this behavior to another subroutine for clarity
* `../TME/src/TME_Module_v28.f90:1334:` Combine `readInputSD()` and `readInputPC()`
