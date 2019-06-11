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
