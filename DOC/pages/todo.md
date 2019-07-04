title: Todo
author: Laura Nichols
date: 07/04/2019

* `../TME/src/TME_Module_v28.f90:247:` Consider changing `atom` type to `element` since it holds more than one atom
* `../TME/src/TME_Module_v28.f90:317:` Change `initializeCalculation()` to have arguments to make clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:366:` Change `readInput()` to have arguments to make clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:376:` Figure out what the difference in PC and SD is
* `../TME/src/TME_Module_v28.f90:398:` Change `checkInitialization()` to have arguments to make clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:406:` Change `file_exists` to `fileExists` in `checkInitialization()`
* `../TME/src/TME_Module_v28.f90:592:` Check if there is any kind of check on `ki` and `kf`. Why was this commented out?
* `../TME/src/TME_Module_v28.f90:662:` Change `readInputPC()` to have arguments so that it is clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:678:` Change `file_exists` to `fileExists` in `readInputPC()`
* `../TME/src/TME_Module_v28.f90:708:` Add information about these variables to top
* `../TME/src/TME_Module_v28.f90:761: Change `(posIonPC(j,ni) , j = 1,3)` to `posIonPC(1:`3,ni)` in `readInputPC()` for clarity
* `../TME/src/TME_Module_v28.f90:833:` Look more into how AE and PS wavefunctions are combined to further understand this
* `../TME/src/TME_Module_v28.f90:834:` Move this behavior to another subroutine for clarity
* `../TME/src/TME_Module_v28.f90:1169:` Combine `readInputSD()` and `readInputPC()`
