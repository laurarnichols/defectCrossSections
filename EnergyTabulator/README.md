# Energy Tabulator

There are 3 different energy differences that need to be tabulated to calculate the transition rates:
* Delta function: total potential energy difference (electronic + ion-ion)
* Zeroth-order matrix element: total (many-body) electronic-only energy difference
* First-order matrix element: eigenvalue difference

It is also convenient to go ahead and calculate various energies for each state that could be used for plotting. Capture into a state within the gap and scattering between various band states require the energies to be calculated differently, so the energies are tabulated independently for each case.

Instead of calculating these energies in all of the different programs in different ways, this program tabulates all of the required energies so that they can easily be read by other programs. Following the format of the other programs, this code will output an energy table file `energyTable.isp.ik` for each spin (`isp`) (or the single selected spin) and k-point (`ik`) for capture. The spins and k-points come from the system that is used for the eigenvalues. For the scattering case where inter-k transitions are allowed, a single file is output for each spin channel.

## Capture

For capture, the energies needed for the transition rate are calculated as below.
* Delta function: total energy difference between relaxed charge states, plus additional carrier energy above WZP reference-carrier energy
* Zeroth-order matrix element: same as delta function but using total energy difference between charge states in the initial positions
* First-order matrix element: eigenvalue difference between initial and final state

The code extracts the total energies from the three different systems (`exportDirInitInit`, `exportDirFinalInit`, and `exportDirFinalFinal`) and the eigenvalues from the specified eigenvalue system (`exportDirEigs`). It is best to use HSE calclations to get the total energies, even if you can't do HSE for all of the matrix elements. If you can't use HSE for the energy differences and/or need to include some additional energy correction in the total energy differences, that can be done through `eCorrectTot`. The eigenvalues are best from the ground-state HSE as well, but if that cannot be used (e.g., if many k-points are used and that isn't feasible), PBE-level eigenvalues can be used. Ground-state eigenvalues should always be used because those are most accurate. Within the conduction band/valence band, PBE eigenvalue differences are okay, but you may need to include corrections to the reference carrier (`eCorrectEigRef`) (e.g., if your reference carrier is on the other side of the gap and the band gap from the PBE level needs to be corrected).

The first-order term uses an eigenvalue difference between the initial and final states. For capture, this includes the distance between the band states and the defect level, but this distance changes in the different charge states. Additionally, in finite supercell sizes there is a ficticious dispersion of the defect level. To address both of these issues, the distance between the lowest conduction band and the defect level at $\Gamma$ in the *initial-state HSE calculation* (`exportDirInitInit`) is used, and the lowest conduction band at $\Gamma$ is used as a reference point for the rest of the bands across all k-points. The same reference is used for the delta-function and zeroth-order eigenvalue differences. The order for eigenvalue differences is switched for electron carriers and hole carriers to represent the energy change of the actual electron making the transition.

*Note:* The code currently uses the first k-point as a reference and does not test that it is the $\Gamma$ point. There is also no current way to correct the energy between the defect and the lowest conduction band for the first-order term.

The code also outputs a `dEPlot.isp.ik` file with the index of each state and the positive electronic energy difference in eV. This was the energy that we used for plotting in the paper, but other values can easily be added here based on the energy you are interested in.

## Scattering
For scattering, the energies needed for the transition rate are calculated as below.
* Delta function: total energy difference between each of the different relaxations corresponding to a carrier being in each of the different k-points/bands allowed
* Zeroth-order matrix element: eigenvalue difference between final and initial state
* First-order matrix element: eigenvalue difference between initial and final state

This setup assumes that the defect system has been relaxed with a carrier in each of the possible band states (half in each channel so as to not consider spin coupling to the defect levels). Each of the calculations should be stored in subdirectories with the pattern `k<ik>_b<ib>`, where `ik` is the k-point index and `ib` is the band index. The k-point indices should match those of the system that the eigenvalues comes from.

For scattering, the eigenvalues should come from the perfect crystal system. The bands for the perfect crystal and defect systems will not necessarily line up. It is assumed that the band bounds and labels correspond to the defect system, while the eigenvalues are read and indexed with an optional integer shift `ibShift_eig`. 

For capture, all band states are calculated and output because it is assumed that all of the bands are present in the calculations. However, each of the different k-point/band states requires an additional relaxation for scattering, so the code just skips the states where the required export files do not exist. Beyond that, the energies are only output for the states where there is a non-zero electronic eigenvalue difference (i.e., the matrix elements are nonzero) and a non-zero energy transfer (positive or negative). This behavior is chosen so that, instead of looping over all possible k-points and bands in the given range, subsequent programs can read in the energy table and only consider the states that are physically meaningful.

The code also ouputs `dEPlot.isp` with the eigenvalue difference from the reference band (`refBand`) at the first k-point and the total energy relative to the state with the lowest total energy. The minimum energy and state where it occurs are also output in the header. 

## Inputs
The inputs should be put in a namelist like
```f90
&inputParams
  ...
/
```
within a file (e.g., `EnergyTabulator.in`) that is passed into the code like `mpirun EnergyTabulator.x < EnergyTabulator.in > EnergyTabulator.out`.

### Variables for capture and scattering
* `energyTableDir`
  * Default `./`
  * Where the energy table(s) should be stored
* `captured`
  * `.true.` (default) or `.false.`
  * If the carrier is captured (i.e., if there is a single final state in the gap)
* `elecCarrier`
  * `.true.` (default) or `.false.`
  * If the carrier is an electron vs a hole (i.e., if the carrier is in the conduction band or the valence band)
* `exportDirEigs`
  * Path to the export directory for the system where the band-state eigenvalues come from
  * Number of spins and k-points is read from this system
  * For capture, the eigenvalues should come from the ground state defect system.
  * For scattering, they should come from the ground state perfect crystal.
* `eCorrectEigRef`
  * Default `0.0` eV
  * Optional correction to eigenvalue difference from reference value
* `ispSelect` 
  * `1` or `2` (default is unused)
  * Gives the option to choose a single spin channel rather than looping
  * Must match up with the number of spin channels in the `exportDirEigs` system
* `iBandIinit`, `iBandIfinal`, `iBandFinit`, `iBandFfinal`
  * Loop boundaries for bands (initial states go from `iBandIinit` to `iBandIfinal` and final states go from `iBandFinit` to `iBandFfinal`)
  * For capture, it is currently assumed that there is a single final state so that `iBandFinit = iBandFfinal`.
  * It is assumed that these bands corresponds to the defect bands (see scattering inputs below for more info).
* `refBand`
  * Reference band to measure eigenvalues from
  * Should be the band where the WZP carrier is for capture and the CBM/VBM for scattering

### Variables for capture only
* `exportDirInitInit`, `exportDirFinalInit`, and `exportDirFinalFinal`
  * Path to export directories for systems with initial charge state/initial positions, final charge state/initial positions, and final charge state/final positions, respectively.
  * The initial and final positions come from relaxations of the initial and final charge states.
* `eCorrectTot`
  * Default `0.0` eV
  * Optional correction on total-energy differences
  * There are no bounds on this variable, so make sure to set a reasonable value.

### Variables for scattering only
* `allStatesBaseDir`
  * Path to the directory that holds the various subdirectories for each of the different states
  * It is assumed that the subdirectories have the form `k<ik>_b<ib>` and that there is a separate relaxation/total energy and export for each possible state
* `singleStateExportDir`
  * Name of the export directory within each subfolder (e.g., `'export'`)
* `ikIinit`, `ikIfinal`, `ikFinit`, and `ikFfinal`
  * Bounds for k-points considered for states
  * Indexed to match the subdirectory names within the `allStatesBaseDir` folder and the k-points as indexed in the export of the system used for eigenvalues (should be the perfect crystal for scattering)
* `ibShift_eig`
  * Default `0`
  * Optional integer shift of bands to line up bands from defect system (used to label each of the subdirectories) and those from the perfect-crystal system (where the eigenvalues are read)
  * It is assumed that the band bounds `iBandIinit`, `iBandIfinal`, `iBandFinit`, and `iBandFfinal` correspond to the defect system.
  * The eigenvalues are then read and indexed by `ib+ibShift_eig`.



