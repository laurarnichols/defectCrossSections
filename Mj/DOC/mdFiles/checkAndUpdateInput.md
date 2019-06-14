# checkAndUpdateInput()

* Set the default value of `abortExecution` to false
* Write out a blank line to the output file
* Check the if the certain variables have the default values set in [`initialize()`](initialize.md)
* If a variable still have the default value, output an error message and set `abortExecution` to true
* If the final mode (`modeF`) is less than the initial mode (`modeI`), output an error message and set `abortExecution` to true
* If any of the variables have the default value still (`abortExecution` is true), abort the program
