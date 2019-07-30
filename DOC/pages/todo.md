title: Todo
author: Laura Nichols
date: 07/30/2019

* `../Mj/src/Mj_Main.f90:13:` Make sure that there is an end timer
* `../Mj/src/Mj_Module_v1.f90:29:` Make sure default value is set for `qPoint` 
* `../Mj/src/Mj_Module_v1.f90:426:` Figure out why there is no \(N\) in this equation in the code 
* `../Mj/src/Mj_Module_v1.f90:436:` Figure out if this needs to be another variable 
* `../Mj/src/Mj_Module_v1.f90:459:` Frequency should actually be in eV/\(\hbar\). Check that that's the case. 
* `../Mj/src/Mj_Module_v1.f90:475:` Figure out why this is assigned to another variable 
* `../Mj/src/Mj_Module_v1.f90:480:` Figure out if still need `di`, `bk`, and `dk` 
* `../Mj/src/Mj_Module_v1.f90:487:` Figure out if need to have this in loop. Why change `nm` in `iknb`? 
* `../Mj/src/Mj_Module_v1.f90:490:` Figure out if should send `nm` as it is immediately modified and not used here 
* `../Mj/src/Mj_Module_v1.f90:494:` Possibly change `besOrderNofModeM` to `modBesOrderNofModeM` 
* `../Mj/src/Mj_Module_v1.f90:495: Figure out why this loop is here. Can not just pass `besOrderNofModeM(:`,j)` 
* `../Mj/src/Mj_Module_v1.f90:514:` Change this to a more efficient algorithm 
* `../Mj/src/Mj_Module_v1.f90:578:` Figure out if expect `modeI` and `modeF` to represent index of magnitude of argument `x` 
* `../Mj/src/Mj_Module_v1.f90:684:` Make this loop more clear 
* `../Sigma/src/Sigma_Module_v4.f90:176:` Merge these dummy characters
* `../Sigma/src/Sigma_Module_v4.f90:224:` Merge dummy variables
* `../Sigma/src/Sigma_Module_v4.f90:288:` Figure out where this file is opened
* `../Sigma/src/Sigma_Module_v4.f90:299:` Figure out where this file is opened
* `../Sigma/src/Sigma_Module_v4.f90:325:` Figure out where this file is opened
* `../Sigma/src/Sigma_Module_v4.f90:329:` Figure out what `abCM` is
* `../TME/src/TME_Main_v9.f90:2:` Add detailed math derivation and summary for main program
* `../TME/src/TME_Main_v9.f90:48:` Figure out if need to allocate space for arrays so soon
* `../TME/src/TME_Main_v9.f90:60:` Figure out if SD and PC `numOfGvecs` should be the same
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
* `../TME/src/TME_Module_v28.f90:1213:` Get actual perfect crystal and solid defect output to test
* `../TME/src/TME_Module_v28.f90:1214:` Figure out if loop should be over `solidDefect` or `perfectCrystal`
* `../TME/src/TME_Module_v28.f90:1215:` Look into `nSpins` to figure out if it is needed
* `../TME/src/TME_Module_v28.f90:1233:` Figure out what this subroutine really does
* `../TME/src/TME_Module_v28.f90:1355:` Figure out what this subroutine really does
* `../TME/src/TME_Module_v28.f90:1407:` Figure out the significance of \(l = l^{\prime}\) and \(m = m^{\prime}\)
* `../TME/src/TME_Module_v28.f90:1439:` Figure out why the difference between SD and PC
* `../TME/src/TME_Module_v28.f90:1452:` Figure out why the difference between SD and PC
* `../TME/src/TME_Module_v28.f90:1489:` Figure out what this subroutine really does
* `../TME/src/TME_Module_v28.f90:1561:` Figure out if this output slows things down significantly
* `../TME/src/TME_Module_v28.f90:1562:` Figure out if formula gives accurate representation of time left
* `../TME/src/TME_Module_v28.f90:1606:` Figure out if this should be `system`
* `../TME/src/TME_Module_v28.f90:1607:` Figure out significance of "qr" point
* `../TME/src/TME_Module_v28.f90:1610:` Test if can just directly store in each atom type's `bes_J_qr`
* `../TME/src/TME_Module_v28.f90:1630:` Figure out if this should be `gDotR`
* `../TME/src/TME_Module_v28.f90:1632:` Figure out why this is called `ATOMIC_CENTER`
* `../TME/src/TME_Module_v28.f90:1633:` Figure out why the difference between SD and PC
* `../TME/src/TME_Module_v28.f90:1662:` Figure out why the difference between SD and PC
* `../TME/src/TME_Module_v28.f90:1677:` Figure out why the difference between SD and PC
* `../TME/src/TME_Module_v28.f90:2172:` Figure out what the purpose of this function is. For plotting?
* `../TME/src/TME_Module_v28.f90:2259:` Figure out why `DHifMin` is needed
* `../TME/src/TME_Module_v28.f90:2260:` Figure out why used difference rather than `==
* `../TME/src/TME_Module_v28.f90:2265:` Figure out why this test is here. All of these should be positive, right?
* `../TME/src/TME_Module_v28.f90:2319:` Figure out where unit 12 file is opened and what it is
* `../TME/src/TME_Module_v28.f90:2354:` Figure out why use `eMin + iE*eBin` rather than `DE`
