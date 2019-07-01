# readInput()

* Define a boolean for whether or not a file exists
* Start timer using `t0`
* Check if output file already exists and delete it if it does
* Open a new output file
* [`call initialize()`](initialize.md) to set default values for input values
* Read input variables from command line (or from file if use `< TME_Input.in`)
* [`call checkInitialization()`](checkInitialization.md) to check that all required variables were input and have values that make sense
* [`call readInputPC()`](readInputPC.md) to ????????
* [`call readInputSD()`](readInputSD.md) to ????????
* Calculate the number of plane waves?? (`numOfPWs`) as the `max` of `numOfPWsPC` and `numOfPWsSD`
