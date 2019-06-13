# MjModule
 
## Define global variables
Integer parameters:
* `dp    = selected_real_kind(15, 307)`
* `int32 = selected_int_kind(5)`
* `int64 = selected_int_kind(15)`
* `iostd = 16, un = 3`

Real (double precision) parameters:
* `pi = 3.1415926535897932_dp`
* `twopi = 2.0_dp*pi`
* `abCM = 0.529177219217e-8_dp`
* `THzToHartree = 1.0_dp/6579.683920729_dp`
* `HartreeToEv  = 27.21138386_dp`
* `eVToHartree  = 1.0_dp/27.21138386_dp`

Character parameters:
* output = `status`

Scalar integers:
* `nAtoms`
* `nOfqPoints`
* `nModes`
* `ios`
* `modeI`
* `modeF`
* `qPoint`

Scalar reals (double precision):
* `ti`
* `tf`
* `t1`
* `t2`
* `temperature`
* `kT`
* `maxDisplacement`

Scalar characters:
* `phononsInput`
* `equilibriumAtomicPositions`
* `newAtomicPositions`
* `QEInput`

Scalar logicals:
* `file_exists`
* `readQEInput`

Allocatable (vector/array) of integers:
* `s2L`

Allocatable (vector/array) of reals (double precision):
* `atomD`
* `atomM`
* `phonQ`
* `phonF`
* `genCoord`
* `atomPosition`
* `newAtomicPosition`
* `wby2kT`
* `phonD`
* `x`
* `Sj`
* `coth`
* `besOrderNofModeM`

Allocatable (vector/array) of characters:
* `elements`

Define [namelist](https://docs.oracle.com/cd/E19957-01/805-4939/6j4m0vnc6/index.html) (group of variables) called `MjInput`, containing
* `QEInput`
* `phononsInput`
* `temperature`
* `equilibriumAtomicPositions`
* `modeI`
* `modeF`
* `qPoint`
* `maxDisplacement`

## Define Subroutines
* [`readInputs()`](readInputs.md)
* [`initialize()`](initialize.md)
* [`checkAndUpdateInput()`](checkAndUpdateInput.md)
* [`readPhonons()`](readPhonons.md)
* [`readAtomicPositions()`](readAtomicPositions.md)
* [`computeGeneralizedDisplacements()`](computeGeneralizedDisplacements.md)
* [`computeVariables()`](computeVariables.md)
* [`arrangeLargerToSmaller()`](arrangeLargerToSmaller.md)
* [`displaceAtoms()`](displaceAtoms.md)
* [`writeNewAtomicPositions()`](writeNewAtomicPositions.md)
* [`exportQEInput()`](exportQEInput.md)
* [`readMjs()`](readMjs.md)
* [`iknb ( n, x, nm, bi) !, di, bk, dk )`](iknb.md)
* [`iknb2(n,x,nm,bi,di,bk,dk)`](iknb2.md)
* [`msta1(x, mp)`](msta1.md) 
* [`msta2(x, n, mp)`](msta2.md) 
* [`envj(n, x)`](envj.md) 
