# Transition matrix element (`TME`)

The `TME` program outputs matrix elements based on the all-electron overlap `<braWfc|ketWfc>` for a range of initial and final band states using all-electron overlaps from the PAW method and the appropriate energy difference to get the correct matrix element. The `Export` code must be run on all VASP calculations and the energies must be tabulated with the `EnergyTabulator` code before `TME` can be run. The code assumes that the WZP method is used, where carriers are excited to achieve desired charge states rather than adding/removing electrons using the jellium method. 

The code assumes that the order of the atom types in the two systems matches, with the possible addition of atom types at the end of one of the files.

Output files are 
* `allElecOverlap.isp.ik` -- for each band in range and a given spin and k-point: initial and final band, all-electron overlaps, and matrix elements
* stdout -- timing information and status updates

The same TME code is used for both the zeroth-order and the first-order matrix elements. For both orders, the overlap part comes from from the `TME` code and the energy difference comes from the [`EnergyTabulator`](../EnergyTabulator) code. The band parameters are also read directly from the energy table.

## Zeroth-order

The zeroth-order term has a single matrix element that is shown in the GaN paper to be $$M_{\text{e}}^{\text{BO}} = \langle \phi_f | \psi_i^0 \rangle (E_f - E_i),$$ where $E_f - E_i$ is the total electronic energy difference of the defect crystal before and after capture and $|\phi_f\rangle$ and $|\psi_i^0\rangle$ are the single-quasiparticle orbitals of the final defect state and the perfect-crystal initial state, respectively. 

The `TME` input file for the zeroth-order should look like
```f90
&TME_Input
  ! Which order of matrix element to calculate
  order = 0

  ! Spin channel selection; default unused
  ispSelect = integer  
  
  ! Systems used for overlaps
  braExportDir    = 'path-to-final-charge-state-initial-positions-export'
  ketExportDir    = 'path-to-perfect-crystal-export'
  energyTableDir = 'path-to-energy-tables'
  
  ! Output info
  outputDir = 'path-to-store-overlap-files' 			! default './TMEs'
/
```
_Note: Do not alter the `&TME_Input` or `/` lines at the beginning and end of the file. They represent a namelist and fortran will not recognize the group of variables without this specific format_

## First-order

The first-order term has a matrix element for each mode: $$M_j = \frac{\varepsilon_i - \varepsilon_f}{\delta q_j} \langle \phi_f^j | \phi_i\rangle,$$ where $\varepsilon_i - \varepsilon_f$ is the eigenvalue energy difference between the initial and final band states and $|\phi_f^j\rangle$ is the final-state defect orbital after displacing the atoms along the phonon directions by a small $\delta q_j$. The formalism assumes that the phonon modes are the same in the initial and final electronic states, so it doesn't matter which electronic state is used as long as the same charge state is used for both input wave functions and the atoms are in the initial positions. In practice, the ground state is usually simplest and fastest.

The equation above for the matrix element assumes that the wave functions are orthogonal, but the numerical overlaps of the non-displaced wave functions $\langle \phi_{l'}|\phi_l\rangle$ are not actually $\delta_{ll'}$. To increase the numerical accuracy of the matrix elements, it is best to subtract the non-displaced wave function overlaps: $$M_j = \frac{\varepsilon_i - \varepsilon_f}{\delta q_j} \left( \langle \phi_f^j | \phi_i\rangle - \langle \phi_f|\phi_i\rangle \right).$$
These overlaps can be found by inputting the non-displaced system for the PC and SD inputs. These baseline overlaps can then be subtracted from the overlaps from all of the mode-shifted overlaps.

The `TME` input file for the first-order term should look like
```f90
&TME_Input
  ! Which order of matrix element to calculate
  order = 1

  ! Spin channel selection; default unused
  ispSelect = integer                                   
  
  ! Systems used for overlaps
  braExportDir    = 'path-to-displaced-defect-export'
  ketExportDir    = 'path-to-undisplaced-defect-export'
  energyTableDir = 'path-to-energy-tables'

  ! Baseline info:
  subtractBaseline = .true. or .false. ! If "orthogonal" overlap should be subtracted
  baselineDir = 'path-to-non-displaced-overlaps'

  ! Parameters to get dq_j
  dqFName = 'path-to-dq-file-from-shifter'
  phononModeJ = integer           ! phonon-mode index
  
  ! Output info
  outputDir = 'path-to-store-overlap-files' 			! default './TMEs'
/
```
_Note: Do not alter the `&TME_Input` or `/` lines at the beginning and end of the file. They represent a namelist and fortran will not recognize the group of variables without this specific format_

## Running and parallelization

To run the code, use something like 
```bash
aprun -n num-procs path-to-package/bin/TME.x -nk num-k-pools < TME.in > TME.out
```

All processors get divided into k-point pools, then plane waves get divided over the processes in each pool. There must be the same number of processes per pool, so `num-procs` must be evenly divisible by `num-k-pools`. The k-points, $G$ vectors, and $G+k$ vectors are split up sequentially (i.e., k-points get assigned to pool 1 first then pool 2, etc.). 
