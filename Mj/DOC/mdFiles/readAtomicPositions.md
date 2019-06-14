# readAtomicPositions()

* Define a scaler integer, `iAtom`, to be used in a loop
* Open the file designated as `equilibriumAtomicPositions`
* Allocate space for the `elements` and `atomPosition` arrays
* Set all values in the `atomPosition` array to 0.0
* For each atom:
	* Read in a line where the first field is the element name, then the next 3 are the coordinate position
* Close the file
