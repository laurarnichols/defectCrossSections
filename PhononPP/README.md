# Phonon post-processing

**Assumes phonons only at Gamma point!** 

This code does post-processing of the phonons and relaxed structures. The current options available are (set by logical variables)
* `calcSj` -- whether or not to calculate the Huang-Rhys factor $S_j$
* `diffOmega` -- if $S_j$ should be tabulated for two different $\omega_j$'s
* `calcDq` -- if the `dq.txt` file representing $\delta q_j$ for the first-order-derivative displacements should be output
* `generateShiftedPOSCARs` -- if POSCAR files shifted from given equilibium positions along the phonon modes should be output. These are used to get the first-order matrix element.
* `calcMaxDisp` -- if the maximum displacement across modes between a pair of atoms given a `shift` should be calculated

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

### `calcSj`

The Huang-Rhys factor, $S_j$, is needed for the line-shape-function integration. It is related to the probability of energy getting transferred into the phonon mode $j$. For scattering, we use $S_j$ directly as a weight to distribute the total energy transfer across the modes. $S_j$ is output in descending order so that the modes that participate the most are at the top of the output file. 

To calculate $S_j$, set `calcSj = .true.` and set `singleDisp` (default `.true.`) for if there is a single displacement (i.e., only two relaxed positions) or relaxations for each of the possible states. 

For a single displacement (`singleDisp = .true.`), the two relaxed geometries must be provided using `initPOSCARFName` and `finalPOSCARFName`. $S_j$ is output in the `Sj.out` file.

For multiple displacements (`singleDisp = .false.`), the path to the energy table output by `EnergyTabulator` (`energyTableDir`) must be provided. `EnergyTabulator` should be run before `PhononPP` in this case so that only the transitions with the required inputs and with allowed energies are considered. The path to the base directory (`CONTCARsBaseDir`) where there are subfolders for each of the relaxed states must also be given. It is assumed that the subfolder naming convention is `k<ik>_b<ib>` and lines up with those used for `EnergyTabulator`. For multiple displacements, there are multiple $S_j$ output files named `Sj.k<iki>_b<ibi>.k<ikf>_b<ibf>.out` based on the initial/final bands and k-points.

_Note: It is also currently assumed that the allowed states are not dependent on the spin channel and that `energyTable.1` exisits._

When calculating $S_j$, there is an option to calculate using a single frequency from a single phonon file or two frequencies using two phonon files. Both $S_j$ and $S_j'$ are calculated using the same $\Delta q_j$ projection of the displacement onto the initial-phonon eigenvectors, and both are output in the same file. Currently, the only option for this output file is for the modes to be sorted in order of descending $S_j$. This mode is controlled by the logical `diffOmega`, which is by default `.false.`.

### First-order items

The first-order matrix element is calculated by shifting the atomic positions along the phonon eigenvectors to calculate a derivative of the electronic wave functions with respect to the nuclear positions. The magnitude of the shift across all atoms is set by `shift` (default 0.01 A). This process is only performed if one of the following options are true: `calcDq`, `generateShiftedPOSCARs`, `calcMaxDisp`. For all of the options, the starting-positions file (`basePOSCARFName`) is required. 

With the `calcDq = .true.` option, the shift magnitude is converted into terms of the generalized coordinate $\delta q_j$ and output in the `dqFName` file (default `./dq.txt`). 

For the `generateShiftedPOSCARs = .true.` option, the POSCARS are output with names `<prefix>-<j>`, where `prefix` is an input option (default `ph_POSCAR`).

The `calcMaxDisp` option calculates the maximum displacement between two atoms indexed in the phonon file and `basePOSCARFName` file by `disp2AtomInd = int1, int2`. 
