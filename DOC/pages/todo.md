title: Todo
author: Laura Nichols
date: 12/01/2020

* `../LSF/src/zerothOrder/LSF_zeroth_Main.f90:107:` Redo the loop for less than 5 phonons to be more clear and streamlined
* `../LSF/src/zerothOrder/LSF_zeroth_Main.f90:166:` Change this to have `size(iEbinsByBands)`
* `../LSF/src/zerothOrder/LSF_zeroth_Main.f90:171:` Add `else iEbinsByPhonons = iEbinsByBands` to remove if below
* `../LSF/src/zerothOrder/LSF_zeroth_Main.f90:172:` Maybe change variable names to be clearer
* `../LSF/src/zerothOrder/LSF_zeroth_Main.f90:185:` Figure out how getting `de` is "calculating DOS" and if not where DOS is
* `../LSF/src/zerothOrder/LSF_zeroth_Main.f90:186:` Figure out why DOS isn't in sum as in formula
* `../LSF/src/zerothOrder/LSF_zeroth_Main.f90:302:` Move the behavior of splitting up Monte Carlo steps to a subroutine
* `../LSF/src/zerothOrder/LSF_zeroth_Main.f90:331:` Move this to a subroutine
* `../LSF/src/zerothOrder/LSF_zeroth_Main.f90:337:` Replace this with `binomialCoefficient(kPhonons-1, kPhonons-nBands)`
* `../LSF/src/zerothOrder/LSF_zeroth_Module_v35.f90:192:` Figure out why increase `minimumNumberOfPhonons` by 1
* `../LSF/src/zerothOrder/LSF_zeroth_Module_v35.f90:227:` Change `nOfqPoints` to `numberOfQPoints` 
* `../LSF/src/zerothOrder/LSF_zeroth_Module_v35.f90:532:` Remove all of these comments 
* `../LSF/src/zerothOrder/LSF_zeroth_Module_v35.f90:596:` Figure out if `iModeFs(myid)` has a max of `nModes-3` or `nModes-nBands+1`
* `../LSF/src/zerothOrder/LSF_zeroth_Module_v35.f90:621:` Figure out what the purpose of `ic` is 
* `../LSF/src/zerothOrder/LSF_zeroth_Module_v35.f90:685:` Change this to merge if statements 
* `../LSF/src/zerothOrder/LSF_zeroth_Module_v35.f90:699:` Figure out why don't just exit here because will be multiplying by 0 
* `../LSF/src/zerothOrder/LSF_zeroth_Module_v35.f90:719:` Redo `besRatio` if statement to be more clear that it is if/else 
* `../LSF/src/zerothOrder/LSF_zeroth_Module_v35.f90:765:` Send slice instead of using `other` 
* `../LSF/src/zerothOrder/LSF_zeroth_Module_v35.f90:778:` Figure out if there is a better way in general to do this 
* `../LSF/src/zerothOrder/LSF_zeroth_Module_v35.f90:779:` Write a recursive function to replace explicit loops 
* `../LSF/src/zerothOrder/LSF_zeroth_Module_v35.f90:780:` Fix typo in `distrubute` 
* `../LSF/src/zerothOrder/LSF_zeroth_Module_v35.f90:781:` Change `l` to `nBands` and `m` to `kPhonons` or something similar 
* `../LSF/src/zerothOrder/LSF_zeroth_Module_v35.f90:822:` Add a condition to exit inner loop if `i1 + i2 + i3 > m` 
* `../LSF/src/zerothOrder/LSF_zeroth_Module_v35.f90:831:` Add a condition to exit inner loop if `i > size of pj0s` 
* `../LSF/src/zerothOrder/LSF_zeroth_Module_v35.f90:1187:` Make sure that Monte Carlo makes sense here 
* `../LSF/src/zerothOrder/LSF_zeroth_Module_v35.f90:1188:` Figure out if there are any methods that would be better/faster 
* `../LSF/src/zerothOrder/LSF_zeroth_Module_v35.f90:1224:` Figure out a better way to do this as it is crazy inefficient 
* `../LSF/src/zerothOrder/LSF_zeroth_Module_v35.f90:1229:` Remove as not needed 
* `../LSF/src/zerothOrder/LSF_zeroth_Module_v35.f90:1250:` Fix the possible bug here 
* `../LSF/src/zerothOrder/LSF_zeroth_Module_v35.f90:1261:` Move some of this to another subroutine 
* `../LSF/src/zerothOrder/LSF_zeroth_Module_v35.f90:1343:` Figure out why array is reversed 
* `../LSF/src/zerothOrder/LSF_zeroth_Module_v35.f90:1383:` Figure out if there is a clearer or faster way to do this 
* `../LSF/src/zerothOrder/LSF_zeroth_Module_v35.f90:1384:` Maybe move this to a subroutine 
* `../LSF/src/zerothOrder/LSF_zeroth_Module_v35.f90:1409:` Switch the `.true.` and `.false.` assignments to make more sense 
* `../LSF/src/zerothOrder/LSF_zeroth_Module_v35.f90:1560:` Replace this with `binomialCoefficient(kPhonons-1, kPhonons-nBands)` 
* `../LSF/src/zerothOrder/LSF_zeroth_Module_v35.f90:1627:` Change the name of this subroutine to just `writeLSF` 
* `../LSF/src/zerothOrder/LSF_zeroth_Module_v35.f90:1628:` Remove all of the extra stuff from this subroutine 
* `../LSF/src/zerothOrder/LSF_zeroth_Module_v35.f90:1727:` Merge this with `parallelIsFsBy4` 
* `../LSF/src/zerothOrder/LSF_zeroth_Module_v35.f90:1737:` Write a binomialCoefficient function 
* `../LSF/src/zerothOrder/LSF_zeroth_Module_v35.f90:1764:` Change this to use available states instead of totalStates 
* `../Mj/src/Mj_Main.f90:17:` Make sure that there is an end timer
* `../Mj/src/Mj_Module_v1.f90:25:` Make sure default value is set for `qPoint` 
* `../Mj/src/Mj_Module_v1.f90:269:` Figure out if expect `modeI` and `modeF` to represent index of magnitude of argument `x` 
* `../Mj/src/Mj_Module_v1.f90:381:` Make this loop more clear 
* `../Sigma/src/Sigma_Module_v4.f90:176:` Merge these dummy characters
* `../Sigma/src/Sigma_Module_v4.f90:224:` Merge dummy variables
* `../Sigma/src/Sigma_Module_v4.f90:290:` Figure out where this file is opened
* `../Sigma/src/Sigma_Module_v4.f90:301:` Figure out where this file is opened
* `../Sigma/src/Sigma_Module_v4.f90:327:` Figure out where this file is opened
* `../Sigma/src/Sigma_Module_v4.f90:331:` Figure out what `abCM` is
* `../TME/src/TME_Main_v9.f90:2:` Add detailed math derivation and summary for main program
* `../TME/src/TME_Main_v9.f90:47:` Figure out if need to allocate space for arrays so soon
* `../TME/src/TME_Main_v9.f90:64:` Figure out if SD and PC `numOfGvecs` should be the same
* `../TME/src/TME_Module_v28.f90:34:` Change I/O from file to console so that usage matches that of QE
* `../TME/src/TME_Module_v28.f90:152:` Consider changing `atom` type to `element` since it holds more than one atom
* `../TME/src/TME_Module_v28.f90:410:` Change `checkInitialization()` to have arguments to make clear that these variables are getting changed
* `../TME/src/TME_Module_v28.f90:552:` Remove everything with `ki` and `kf` because never used
* `../TME/src/TME_Module_v28.f90:838:` Look more into how AE and PS wavefunctions are combined to further understand this
* `../TME/src/TME_Module_v28.f90:839:` Move this behavior to another subroutine for clarity
* `../TME/src/TME_Module_v28.f90:848:` Figure out if differences in PC and SD `F1` calculations are intentional
* `../TME/src/TME_Module_v28.f90:849:` Figure out if should be `(wae_i wae_j - wps_i wps_j)r_{ab}`
* `../TME/src/TME_Module_v28.f90:850:` Figure out if first term in each should be conjugated for inner product form
* `../TME/src/TME_Module_v28.f90:851:` Figure out if `rab` plays role of \(dr\) within augmentation sphere
* `../TME/src/TME_Module_v28.f90:890:` Figure out if intentional to only use `JMAX` from SD input
* `../TME/src/TME_Module_v28.f90:1220:` Get actual perfect crystal and solid defect output to test
* `../TME/src/TME_Module_v28.f90:1221:` Figure out if loop should be over `solidDefect` or `perfectCrystal`
* `../TME/src/TME_Module_v28.f90:1222:` Look into `nSpins` to figure out if it is needed
* `../TME/src/TME_Module_v28.f90:1240:` Figure out what this subroutine really does
* `../TME/src/TME_Module_v28.f90:1344:` Figure out what this subroutine really does
* `../TME/src/TME_Module_v28.f90:1396:` Figure out the significance of \(l = l^{\prime}\) and \(m = m^{\prime}\)
* `../TME/src/TME_Module_v28.f90:1428:` Figure out why the difference between SD and PC
* `../TME/src/TME_Module_v28.f90:1441:` Figure out why the difference between SD and PC
* `../TME/src/TME_Module_v28.f90:1478:` Figure out what this subroutine really does
* `../TME/src/TME_Module_v28.f90:1550:` Figure out if this output slows things down significantly
* `../TME/src/TME_Module_v28.f90:1551:` Figure out if formula gives accurate representation of time left
* `../TME/src/TME_Module_v28.f90:1595:` Figure out if this should be `system`
* `../TME/src/TME_Module_v28.f90:1596:` Figure out significance of "qr" point
* `../TME/src/TME_Module_v28.f90:1599:` Test if can just directly store in each atom type's `bes_J_qr`
* `../TME/src/TME_Module_v28.f90:1619:` Figure out if this should be `gDotR`
* `../TME/src/TME_Module_v28.f90:1621:` Figure out why this is called `ATOMIC_CENTER`
* `../TME/src/TME_Module_v28.f90:1622:` Figure out why the difference between SD and PC
* `../TME/src/TME_Module_v28.f90:1651:` Figure out why the difference between SD and PC
* `../TME/src/TME_Module_v28.f90:1666:` Figure out why the difference between SD and PC
* `../TME/src/TME_Module_v28.f90:2163:` Figure out what the purpose of this function is. For plotting?
* `../TME/src/TME_Module_v28.f90:2252:` Figure out why `DHifMin` is needed
* `../TME/src/TME_Module_v28.f90:2253:` Figure out why used difference rather than `==
* `../TME/src/TME_Module_v28.f90:2258:` Figure out why this test is here. All of these should be positive, right?
* `../TME/src/TME_Module_v28.f90:2312:` Figure out where unit 12 file is opened and what it is
* `../TME/src/TME_Module_v28.f90:2347:` Figure out why use `eMin + iE*eBin` rather than `DE`
