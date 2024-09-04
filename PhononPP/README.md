# Phonon post-processing

**Assumes phonons only at Gamma point!** 

This code does post-processing of the phonons and relaxed structures. The current options available are (set by logical variables)
* `calcSj` -- whether or not to calculate the Huang-Rhys factor $S_j$
* `diffOmega` -- if $S_j$ should be tabulated for two different $\omega_j$'s (only implemented for single displacement/capture)
* `calcDq` -- if the `dq.txt` file representing $\delta q_j$ for the first-order-derivative displacements should be output
* `generateShiftedPOSCARs` -- if POSCAR files shifted from given equilibium positions along the phonon modes should be output. These are used to get the first-order matrix element.
* `calcMaxDisp` -- if the maximum displacement across modes between a pair of atoms given a `shift` should be calculated
* `calcDeltaNj` -- if change of occupation numbers for carrier approach, transition, and departure should be calculated (only for scattering)

The code is structured so that each of the options above can be chosen independently. Each of the options and the inputs needed are given in more detail below. 

The inputs should be put in a namelist like
```f90
&inputParams
  ...
/
```
within a file (e.g., `PhononPP.in`) that is passed into the code like `mpirun PhononPP.x < PhononPP.in > PhononPP.out`. Parallelization is implemented, but is usually not needed. Although you may be able to get away with running without submitting to the queue because the run time is in seconds, the array sizes needed are not insignificant, so it is best practice to submit it to the queue. 

## Different calculations and inputs

All calculations require reading the `mesh.yaml` Phonopy output file. The path to the file (`phononFName`) and a threshold for the absolute value of the frequencies (`freqThresh`) are required. By default, `freqThresh = 0.5`. This threshold is needed to determine which modes to skip (translational modes) and which modes to consider. We use the absolute value in case there are negative-frequency modes that would indicate that the structure used for the calculation of phonons is unstable and at a saddle point relative to displacement along the eigenvector of those negative-frequency modes. 

The code also always calculates the thermal occupation numbers $n_j$ based on an input `temperature`. The default temperature used is 300 K.

For all calculations that require projection onto the phonon eigenvectors, when `diffOmega = .true.` the user has the option to set `dqEigvecsFinal` to `.true.` or `.false.` to determine if the final (prime) mode eigenvectors will be used for the projections (`dqEigvecsFinal = .true.`) or the initial-mode eigenvectors (`.false.`). The default is `.true.`. Note that the `diffOmega` option has not been implemented for the scattering/multiple displacements case because the theory is not clear on what the different frequencies would be and we only consider a single set of frequencies for scattering. 

### `calcSj`

The Huang-Rhys factor, $S_j$, is needed for the line-shape-function integration. It is related to the probability of energy getting transferred into the phonon mode $j$. For scattering, we use $S_j$ directly as a weight to distribute the total energy transfer across the modes. $S_j$ is output in descending order so that the modes that participate the most are at the top of the output file. 

To calculate $S_j$, set `calcSj = .true.` and set `singleDisp` (default `.true.`) for if there is a single displacement (i.e., only two relaxed positions) or relaxations for each of the possible states. 

For a single displacement (`singleDisp = .true.`), the two relaxed geometries must be provided using `initPOSCARFName` and `finalPOSCARFName`. $S_j$ is output in the `Sj.out` file.

For multiple displacements (`singleDisp = .false.`), the path to the energy table output by `EnergyTabulator` (`energyTableDir`) must be provided, as well as the selected spin channel (`ispSelect`) corresponding to that energy table. `EnergyTabulator` should be run before `PhononPP` in this case so that only the transitions with the required inputs and with allowed energies are considered. The path to the base directory (`allStatesBaseDir_relaxed`) where there are subfolders for each of the relaxed states must also be given. It is assumed that the subfolder naming convention is `k<ik>_b<ib>` and lines up with those used for `EnergyTabulator`. For multiple displacements, there are multiple $S_j$ output files named `Sj.k<iki>_b<ibi>.k<ikf>_b<ibf>.out` based on the initial/final bands and k-points and are stored in a subfolder called `transitions`.

When calculating $S_j$, there is an option to calculate using a single frequency from a single phonon file or two frequencies using two phonon files. Both $S_j$ and $S_j'$ are calculated using the same $\Delta q_j$ projection of the displacement onto the initial-phonon eigenvectors, and both are output in the same file. Currently, the only option for this output file is for the modes to be sorted in order of descending $S_j$. This mode is controlled by the logical `diffOmega`, which is by default `.false.`.

There is also a file `Sj.analysis.out` created. For capture, this will just be one row of data indicating the max $S_j/S_j'$ and the number above a threshold `SjThresh` given by the user. The default value of `SjThresh` is `0.1`. For scattering, there are independent $S_j$'s for each transition, so there will be a row with the relevant data points for each transition, with the indices of the involved states given in the first column. 

## `calcDeltaNj`

**Important:** To use this option, you must set `singleDisp = .false.` and `calcSj = .true.`.

For the scattering problem (`singleDisp = .false.`), we need to know how the occupation numbers are changed as the carrier approaches, transitions, and departs. The change in occupation numbers for each process is given by 
```math
\Delta n_j = \frac{\Delta E}{\hbar \omega_j} \frac{\omega_j S_j}{\sum_{j'} \omega_{j'} S_{j'}},
```
where $\Delta E$ is the energy difference resulting from each process (tabulated in `EnergyTabulator`) and each process has its own $S_j$ corresponding to a different change in equilibrium positions $\Delta q_j$. For scattering, the displacement also depends on the electronic state because we relax each state separately. 

To run this option, set `calcDeltaNj = .true.`, `singleDisp = .false.`, `calcSj = .true.`, and pass in the path to the total-energy calculations for each electronic state in the start/ground-state positions using `allStatesBaseDir_startPos`. It is assumed that the subfolder naming convention is `k<ik>_b<ib>` and lines up with those used for `EnergyTabulator`. The output files are named `deltaNj.k<iki>_b<ibi>.k<ikf>_b<ibf>.out` based on the initial/final bands and k-points and are stored in a subfolder called `deltaNjs`.

### First-order items

The first-order matrix element is calculated by shifting the atomic positions along the phonon eigenvectors to calculate a derivative of the electronic wave functions with respect to the nuclear positions. The magnitude of the shift across all atoms is set by `shift` (default 0.01 A). This process is only performed if one of the following options are true: `calcDq`, `generateShiftedPOSCARs`, `calcMaxDisp`. For all of the options, the starting-positions file (`basePOSCARFName`) is required. 

With the `calcDq = .true.` option, the shift magnitude is converted into terms of the generalized coordinate $\delta q_j$ and output in the `dqFName` file (default `./dq.txt`). 

For the `generateShiftedPOSCARs = .true.` option, the POSCARS are output with names `<prefix>-<j>`, where `prefix` is an input option (default `ph_POSCAR`).

The `calcMaxDisp` option calculates the maximum displacement between two atoms indexed in the phonon file and `basePOSCARFName` file by `disp2AtomInd = int1, int2`. 
