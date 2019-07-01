# Transition Matrix Elements (TME)

* Pull in `mpi` and `declarations` modules
* Declare real scalar variables
  * `t1`
  * `t2`
* Initialize mpi environment
* Get the id of the current process
* Get the number of processes
* Allocate space for `nPWsI` and `nPWsF`
* If root process
  * Start a timer
  * [`call readInput()`](readInput.md)
  * [`call readPWsSet()`](readPWsSet.md)
  * Allocate space for variables
  * Initialize entire `Ufi` matrix to complex, double zero
  * [`call distributePWsToProcs(numOfGvecs, numprocs)`](distributePWsToProcs.md)
  * Initialize number of initial and final plane waves to zero for each process
  * For each process, calculate the number of initial and final plane waves
* Broadcast variables from root process to all other processes
* For each `numOfTypesPC`
  * Broadcast variables in `atomsPC` structure array (see [this post](https://stackoverflow.com/questions/8100131/what-does-mean-do-in-fortran) for an explanation of the use of `%`)
