# checkInitialization()

* Define booleans for if a file exists and whether or not to abort the program
* Set the default value of `abortExecution` to `.false.`
* Output line to output file with `Inputs: ` to say we are going through the inputs
* If the SD export directory (`exportDirSD`) is blank
  * Output an error and set `abortExecution` to `.true.`
* Else if the SD export directory doesn't exist
  * Output an error and set `abortExecution` to `.true.`
* Output what was given as the SD export directory
* Do the same thing for the PC export directory (`exportDirPC`)
* If the elements path (`elementsPath`) is blank
  * Output a warning message and set the default value
* If the elements path directory doesn't already exist
  * Create the directory
* Output the value of elements path
* Check to see if `iBandIinit`, `iBandIfinal`, `iBandFinit`, or `iBandFfinal` are still less than zero
  * If any of them still have the default value, output an error message and set `abortExecution` to `.true.`
* Output the values of `iBandIinit`, `iBandIfinal`, `iBandFinit`, and `iBandFfinal`
* If `calculateVfis` is `true` and `iBandFinit` and `iBandFfinal` are not equal
  * Output an error message and set `abortExecution` to `.true.`
* Output the value of `calculateVfis`
* If the `VfisOutput` file name is blank
  * Output a warning message and set the default value
* Output the value of `VfisOutput`
* If `eBin` is still less than zero
  * Output a warning message and set the default value
* Output the value of `eBin`
* Convert `eBin` from eV to Hartree
