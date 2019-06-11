# Calculate Carrier Cross Sections
## How to Clone
Make sure you have `git` installed and type
```markdown
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
* Open the `Makefile` and edit the path to Quantum Espresso to match your system 

_Note: Make sure that your path does not have a `/` at the end or there will be an error_
* Generate the `make.sys` file by typing `make initialize` 
* Copy over the `make.sys` file generated in the `carrierCrossSections` directory to the `q-e-qe-5.3` folder
* Change into the `q-e-qe-5.3` directory and run `./configure`
* Confirm that the last line of output is `configure: success`
* Run `make pw pp` to build the PW and PP packages

### Make Your Target
* Change back into the base directory
* You should now be able to make the target you want (e.g. `make all_QE-5.3.0`)
* For a list of some possible targets, read through the `Makefile` or type `make`
