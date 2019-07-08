title: Todo
author: Laura Nichols
date: 07/08/2019

* Binary file ../TME/src/.TME_Module_v28.f90.swp matches
* `../TME/src/TME_Module_v28.f90:244:` Consider changing `atom` type to `element` since it holds more than one atom
* `../TME/src/TME_Module_v28.f90:374:` Change `readInput()` to have arguments to make clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:384:` Figure out what the difference in PC and SD is
* `../TME/src/TME_Module_v28.f90:406:` Change `checkInitialization()` to have arguments to make clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:414:` Change `file_exists` to `fileExists` in `checkInitialization()`
* `../TME/src/TME_Module_v28.f90:600:` Check if there is any kind of check on `ki` and `kf`. Why was this commented out?
* `../TME/src/TME_Module_v28.f90:670:` Change `readInputPC()` to have arguments so that it is clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:759:` Add information about these variables to top
* `../TME/src/TME_Module_v28.f90:812: Change `(posIonPC(j,ni) , j = 1,3)` to `posIonPC(1:`3,ni)` in `readInputPC()` for clarity
* `../TME/src/TME_Module_v28.f90:884:` Look more into how AE and PS wavefunctions are combined to further understand this
* `../TME/src/TME_Module_v28.f90:885:` Move this behavior to another subroutine for clarity
* `../TME/src/TME_Module_v28.f90:1220:` Combine `readInputSD()` and `readInputPC()`
