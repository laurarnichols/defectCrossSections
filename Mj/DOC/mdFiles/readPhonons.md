# readPhonons()

## Define variables
* Define scalar integers for loop indices:
	* `iAtom`
	* `iMode`
	* `iQ`
* Define scalar reals (double precision):
	* `dummyD`
	* `freqInTHz`
* Define scalar character:
	* `dummyC`

## Read from file
* Open the file given as `phononsInput`
* Read the number of q points (`nOfqPoints`) and the number of atoms (`nAtoms`) from the file
* Write those variables to the output file
* Calculate the number of modes as `nModes = 3*nAtoms - 3`
* Read a blank line from the file
* Allocate space for the ?? array (`atomD`) and the ?? array (`atomM`)
* Set all of the values in those arrays to 0.0
* For each atom, read in a line where the first 3 values go in a column of the `atomD` array and the last value goes in the `atomM` array
* Read a blank line from the file
* Allocate space for the ?? (`phonQ`), phonon frequency (`phonF`), and phonon displacement (`phonD`) arrays
* Set all of the values in those arrays to 0.0
* For each integer up to the number of q points:
	* Read in a line like `q = ( 0.0 0.0 0.0 )` and put only the numbers in a column in the `phonQ` array
	* For each integer up to the number of modes:
		* Read a blank line from the file
		* Read in the frequency in THz from the file
		* Convert the frequency to Hartree
		* Read a trash line (column headers)
		* For each atom:
			* Read a line where you throw away the first 3 values and store the last 3 values (phonon displacements) in the `phonD` array
* Close the input file
