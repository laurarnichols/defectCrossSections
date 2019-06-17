# Calculate Carrier Cross Sections
## How to Clone
Make sure you have `git` installed and type
```git
git clone https://github.com/laurarnichols/carrierCrossSections.git
```

## How to Compile
### Install Dependencies
The commands for Ubuntu are

* `sudo apt-get update`
* `sudo apt install gfortran`
* `sudo apt install libmpich-dev`
* `sudo apt install libopenmpi-dev`

### Set up Quantum Espresso
* Download Quantum Espresso 5.3 from the [GitHub repo](https://github.com/QEF/q-e/releases?after=qe-6.2.0) (click on `tar.gz` under `qe-5.3`)
* Change to whatever directory you saved the tar file in on the CLI
* Decompress the tar file using `tar -xvzf q-e-qe-5.3.tar.gz`
* Change into the `q-e-qe-5.3` directory and run `./configure`
* Confirm that the last line of output is `configure: success`
* Run `make pw pp` to build the PW and PP packages

### Make Your Target
* Change into the `carrierCrossSections` directory
* Open the `Makefile` and edit the path to Quantum Espresso to match your system 

_Note: Make sure that your path does not have a `/` at the end or there will be an error_
* You should now be able to make the target you want (e.g. `make all_QE-5.3.0`)
* For a list of some possible targets, read through the `Makefile` or type `make`

## How to Run
I don't know how to run yet, so this section contains notes to try to figure that out.

* TME stands for transition matrix element
* LSF stands for line shape function
* Sigma is likely the cross section calculation
* Mj is likely calculating equation 46 (for linear phonon TMEs)
* If cross section (sigma) is goal (zeroth order):
  * TME and LSF --> Sigma ?
  * QE --> TME ?
  
## Input and Output

* QE:
  * Input:
    * `scf.in`
  * Output: 
    * `scf.out`
    * output directory (e.g. `tmp`)
* Export:
  * Input:
    * output directory
  * Output:
    * export directory
* TME:
  * Input:
    * export directory
  * Output: 
    * VfisVsE
* LSF:
  * Input:
    * phonons_disp.dat --> where does this come from?
  * Output:
    * lsfVsEwithUpTo<%i>phonons
* Sigma
  * Input:
    * VfisVsE
    * lsfVsEwithUpTo<%i>phonons
   
  
