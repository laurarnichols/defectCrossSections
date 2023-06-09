# Transition matrix element (`TME`)

The `TME` program outputs matrix elements for a range of initial and final band states using all-electron overlaps from the PAW method and the appropriate energy difference to get the correct matrix element. The `Export` code must be run on all VASP calculations before `TME` can be run. The code assumes that the WZP method is used, where carriers are excited to achieve desired charge states rather than adding/removing electrons using the jellium method. 

The code also assumes that all atom types in the `PC` system are also in the `SD` system in the same order and that any atoms in the SD system not in the PC system are at the end.

Output files are 
* `allElecOverlap.isp.ik` -- for each band in range and a given spin and k-point: initial and final band, all-electron overlaps, and matrix elements
* `output` -- information about input variables
* `VfisVsEofKpt`/`VfisVsE` -- optional k-binned matrix elements
* stdout -- timing information and status updates

The same TME code is used for both the zeroth-order and the first-order matrix elements. For both orders, the overlap part comes from from the TME code and the energy difference comes from the [`EnergyTabulator`](../EnergyTabulator) code.

## Zeroth-order

The zeroth-order term has a single matrix element that is shown in the GaN paper to be $$M_{\text{e}}^{\text{BO}} = \langle \phi_f | \psi_i^0 \rangle (E_f - E_i),$$ where $E_f - E_i$ is the total electronic energy difference of the defect crystal before and after capture and $|\phi_f\rangle$ and $|\psi_i^0\rangle$ are the single-quasiparticle orbitals of the final defect state and the perfect-crystal initial state, respectively. 

The `TME` input file for the zeroth-order should look like
```f90
&TME_Input
  ! Which order of matrix element to calculate
  order = 0
  
  ! Systems used for overlaps
  exportDirSD = 'path-to-final-charge-state-initial-positions-export'
  exportDirPC = 'path-to-perfect-crystal-export'
  
  ! Band range for overlaps
  iBandIinit = integer						! lowest initial-state band
  iBandIfinal = integer						! highest initial-state band
  iBandFinit = integer						! lowest final-state band
  iBandFfinal = integer						! highest final-state band
    ! Note: code logic only currently tested for single final band state (iBandFinit = iBandFfinal)
  
  ! Output info
  elementsPath = 'path-to-store-overlap-files' 			! default './TMEs'
  calculateVfis = .true. or .false.				! if k-binned overlaps should be calculated and output; default .false.
  eBin = real							! size of energy bin if used
  VfisOutput = 'path-to-write-k-binned-overlaps' 		! default './VfisVsE'
/
```
_Note: Do not alter the `&TME_Input` or `/` lines at the beginning and end of the file. They represent a namelist and fortran will not recognize the group of variables without this specific format_

## First-order

The first-order term has a matrix element for each mode: $$M_j = \frac{\varepsilon_i - \varepsilon_f}{\delta q_j} \langle \phi_f^j | \phi_i\rangle,$$ where $\varepsilon_i - \varepsilon_f$ is the eigenvalue energy difference between the initial and final band states and $|\phi_f^j\rangle$ is the final-state defect orbital after displacing the atoms along the phonon directions by a small $\delta q_j$. The formalism assumes that the phonon modes are the same in the initial and final electronic states, so it doesn't matter which electronic state is used as long as the same charge state is used for both input wave functions and the atoms are in the initial positions. In practice, the ground state is usually simplest and fastest.

The `TME` input file for the first-order term should look like
```f90
&TME_Input
  ! Which order of matrix element to calculate
  order = 1
  
  ! Systems used for overlaps
  exportDirSD = 'path-to-displaced-defect-export'
  exportDirPC = 'path-to-undisplaced-defect-export'
  
  ! Band range for overlaps
  iBandIinit = integer						! lowest initial-state band
  iBandIfinal = integer						! highest initial-state band
  iBandFinit = integer						! lowest final-state band
  iBandFfinal = integer						! highest final-state band
    ! Note: code logic only currently tested for single final band state (iBandFinit = iBandFfinal)

  ! Parameters to get dq_j
  dqFName = 'path-to-dq-file-from-shifter'
  phononModeJ = integer           ! phonon-mode index
  
  ! Output info
  elementsPath = 'path-to-store-overlap-files' 			! default './TMEs'
  calculateVfis = .true. or .false.				! if k-binned overlaps should be calculated and output; default .false.
  eBin = real							! size of energy bin if used
  VfisOutput = 'path-to-write-k-binned-overlaps' 		! default './VfisVsE'
/
```
_Note: Do not alter the `&TME_Input` or `/` lines at the beginning and end of the file. They represent a namelist and fortran will not recognize the group of variables without this specific format_

## Running and parallelization

To run the code, use something like 
```bash
aprun -n num-procs path-to-package/bin/TME.x -nk num-k-pools < TME.in > TME.out
```

All processors get divided into k-point pools, then plane waves get divided over the processes in each pool. There must be the same number of processes per pool, so `num-procs` must be evenly divisible by `num-k-pools`. The k-points, $G$ vectors, and $G+k$ vectors are split up sequentially (i.e., k-points get assigned to pool 1 first then pool 2, etc.). 
