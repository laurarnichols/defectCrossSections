# readInputs()

* Check if file output exists. If it does, delete it
* Open a new output file
* [`call initialize()`](initialize.md)
* Read the input variables (`MjInput` namelist; see [Mj_module_v1.md](Mj_module_v1.md))
* [`call checkAndUpdateInput()`](checkAndUpdateInput.md)
* [`call readPhonons()`](readPhonons.md)
* [`call readAtomicPositions()`](readAtomicPositions.md)
