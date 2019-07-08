title: Todo
author: Laura Nichols
date: 07/08/2019

* `../TME/src/TME_Module_v28.f90:244:` Consider changing `atom` type to `element` since it holds more than one atom
* `../TME/src/TME_Module_v28.f90:373:` Change `readInput()` to have arguments to make clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:383:` Figure out what the difference in PC and SD is
* `../TME/src/TME_Module_v28.f90:405:` Change `checkInitialization()` to have arguments to make clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:413:` Change `file_exists` to `fileExists` in `checkInitialization()`
* `../TME/src/TME_Module_v28.f90:599:` Check if there is any kind of check on `ki` and `kf`. Why was this commented out?
* `../TME/src/TME_Module_v28.f90:669:` Change `readInputPC()` to have arguments so that it is clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:760:` Add information about these variables to top
* `../TME/src/TME_Module_v28.f90:826: Change `(posIonPC(j,ni) , j = 1,3)` to `posIonPC(1:`3,ni)` in `readInputPC()` for clarity
* `../TME/src/TME_Module_v28.f90:898:` Look more into how AE and PS wavefunctions are combined to further understand this
* `../TME/src/TME_Module_v28.f90:899:` Move this behavior to another subroutine for clarity
* `../TME/src/TME_Module_v28.f90:1234:` Combine `readInputSD()` and `readInputPC()`
