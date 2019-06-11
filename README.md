# Calculate Carrier Cross Sections
## How to Clone
Make sure you have `git` installed and type
```markdown
git clone https://github.com/laurarnichols/carrierCrossSections.git
```

## How to Compile
First, make sure that you have all of the dependencies installed. The commands for Ubuntu are

* `sudo apt-get update`
* `sudo apt install gfortran`
* `sudo apt install libmpich-dev`
* `sudo apt install libopenmpi-dev`

Generate the `make.sys` file by typing `make all_QE-5.3.0`. The first round will fail, but you can copy over the 
`make.sys` file generated in the `carrierCrossSections` folder to the `q-e-qe-5.3` folder.
