title: Todo
author: Laura Nichols
date: 07/03/2019

* `../TME/src/TME_Module_v28.f90:245:` Consider changing `atom` type to `element` since it holds more than one atom
* `../TME/src/TME_Module_v28.f90:307:` Change `readInput()` to have arguments to make clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:313:` Change `file_exists` to `fileExists` in `readInput()`
* `../TME/src/TME_Module_v28.f90:338:` Figure out what the difference in PC and SD is
* `../TME/src/TME_Module_v28.f90:357:` Change `initialize()` to have arguments to make clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:392:` Change `checkInitialization()` to have arguments to make clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:400:` Change `file_exists` to `fileExists` in `checkInitialization()`
* `../TME/src/TME_Module_v28.f90:586:` Check if there is any kind of check on `ki` and `kf`. Why was this commented out?
* `../TME/src/TME_Module_v28.f90:656:` Change `readInputPC()` to have arguments so that it is clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:671:` Change `file_exists` to `fileExists` in `readInputPC()`
* `../TME/src/TME_Module_v28.f90:701:` Add information about these variables to top
* `../TME/src/TME_Module_v28.f90:754: Change `(posIonPC(j,ni) , j = 1,3)` to `posIonPC(1:`3,ni)` in `readInputPC()` for clarity
* `../TME/src/TME_Module_v28.f90:826:` Look more into how AE and PS wavefunctions are combined to further understand this
* `../TME/src/TME_Module_v28.f90:827:` Move this behavior to another subroutine for clarity
