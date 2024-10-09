# Line-shape function (`LSF`)

The `LSF` program takes in matrix elements, an energy table, the Huang-Rhys factor, and the phonon frequencies and occupation numbers and performs the time-domain integration to get the final transition rate. Before this code can run, you need to run the `TME` program for the proper-order matrix elements, the `EnergyTabulator` program to get the band bounds and table of energies for the delta function, and the `PhononPP` program that post-processes the phonons to get the Huang-Rhys factors and phonon frequencies and occupations. 

The output files are 
  * `transitionRate.isp.ik`/`transitionRate.isp`/`transitionRate.isp.iRt`, where `iRt` is the real-time-integration index
     * initial band index and transition rate (1/s) (k-point index not present for scattering)
     * stored in `transRateOutDir`, which defaults to `'./transitionRates'`
  * `nj.iRt.out`, where `iRt` is the real-time-integration index
     * new occupations and occupation rates of change for scattering with the `generateNewOccupations` option
     * *Note that the first step is the initial occupations and the rate of change at the first step.*
     * stored in `njNewOutDir`, which defaults to `'./njNew'`
  * stdout -- timing information and status updates

## Inputs

The `LSF` input file should look like
```f90
&inputParams
 ...
/
```
with appropriate input variables specified. _Note: Do not alter the `&inputParams` or `/` lines at the beginning and end of the file. They represent a namelist and fortran will not recognize the group of variables without this specific format_

### Base

Every calculation must have the following inputs:
* `captured` -- default `.true.`; if the carrier is captured or not; affects how inputs are read and states are treated
* `iSpin` -- default `1`; spin index (`1` or `2`)
* `order` -- no default; what order to calculate (`0` or `1`)
* `energyTableDir` -- path to `EnergyTabulator` output directory
   *  $\delta$-function energy from `energyTable.isp.ik`/`energyTable.isp` is always used.
   *  Matrix-element energy is used if `newEnergyTable = .true.`
   *  Energy to band edge from `dEPlot.isp` used for scattering (`captured = .false.`) and `generateNewOccupations = .true.`
* `njBaseInput` -- path to `nj` file to be used for the starting occupations for each mode
* `SjBaseDir` -- path to directory holding Huang-Rhys files `Sj.out` (or indexed by transition for scattering) output by `PhononPP`
* `dtau` -- default `1e-4` **(Hartree atomic units)**; time step size for time-domain integration of phonon overlap integral
* `hbarGamma` -- default `0.0`; determines smearing of integrand based on $\exp(-\gamma \tau)$ to achieve convergence
* `smearingExpTolerance` -- default `0.0`; max time for integration is determined by this threshold set on the smearing exponential
* `diffOmega` -- default `.false.`; if there is a different frequency before and after transition; only currently tested for capture; should be consistent with `PhononPP` output
* `transRateOutDir` -- default `'./transitionRates'`; where to output transition-rate files

*There is an option to input `SjThresh`, but that input is not fully tested, so don't use it!* There is also an `oldFormat` option that can be useful if you are updating the code and needing to run `LSF` with two different formats. This option will most likely not be used, though.

### Matrix elements

To specify the path to the zeroth-order (`order = 0`) matrix elements, simply use `matrixElementDir = 'path-to-matrix-element-files'`. The specify the path to the first-order (`order = 1`) matrix elements, set `MjBaseDir = 'path-to-main-dir-for-first-order-matrix-elements'`. The code assumes that the matrix elements for each mode are in separate sub-folders with names `prefixXXX` with the `prefix` variable set on input and the `XXX` representing the mode index with the length set by the maximum width of the mode indices. The suffix length is set automatically. *Within each folder*, the matrix elements should be located at the path `matrixElementDir` (default `'./'`).

Running the `TME` code to get all of the first-order matrix elements is very expensive, but the expensive part is calculating all of the wave function overlaps needed. However, it can be valuable to use the same overlaps but get results for different energies or $\Delta q_j$'s. You can do this using the `newEnergyTable` and `rereadDq` options. If `newEnergyTable = .true.`, only the overlaps will be read from the `allElecOverlap.isp.ik`/`allElecOverlap.isp` files, and the energies for the matrix elements will be read from the file(s) in the `energyTableDir`. For the first-order term, the default behavior for `newEnergyTable` is to read $\Delta q_j$ from the header of the matrix element file(s), but you can also set `rereadDq = .true.` to read the $\Delta q_j$'s from the `dqInput` file. 

With a change in the frequencies after capture (`diffOmega = .true.`), the modes are indexed by the initial-state mode indices. The modes in the final-state phonons have to be matched up to those in the initial state because the vibrations are sorted on frequency in the `mesh.yaml` file and will not be in the same order. This lining-up of the modes means that if you want to use the overlaps to compare using one set of phonons vs two, the indexing of the matrix elements might not match the single-mode outputs from `PhononPP`. For this case, you can use `reSortMEs = .true.` to shuffle the indices for the matrix elements and `dqInput` if used. The mapping between the mode indices is given in the `optimalPairs.out` file from `PhononPP`, so with `reSortMEs = .true.`, you must provide the path to that file in `optimalPairsInput`.

### Scattering

For scattering (`captured = .false.`), there are several additional options:
* `addDeltaNj` -- default `.false.`; use the occupations adjusted after the initial carrier approach for calculating the transition rate
* `deltaNjBaseDir` -- path to the changes in occupations for carrier approach, transition, and departure (output by `PhononPP`)
* `generateNewOccupations` -- if average over transitions should be performed and new occupations should be calculated
* `dt` -- default `1e-4` **(s)**; real-time step to be multiplied by average rate of change of occupations to get $\Delta \bar{n}_j$
* `nRealTimeSteps` -- default `1`; number of real time steps to take if `generateNewOccupations = .true.`
* `carrierDensityInput` -- path to carrier density file used to perform integration over initial states
* `energyAvgWindow` -- default `1e-2` eV; energy window over which to average noisy carrier-density input
* `njNewOutDir` -- default `'./njNew'`; where to store new occupation files

Note that all of our code is in Hartree atomic units until the calculation of the new occupations. The energies from the band edge given in `dEPlot.isp` are in eV and the carrier density is in $\mathrm{cm}^{-3}$ $\mathrm{eV}^{-1}$, so the `energyAvgWindow` for evaluating the carrier density at the energies in `dEPlot.isp` is in eV. After integration over the initial states multiplied by the carrier density as a function of energy, the resulting rate of change of the occupations is in $\mathrm{s}^{-1}$, so the time step `dt` for the real-time integration of the occupations is in s. 

## Equations

### Zeroth-order

The zeroth-order transition rate is 
```math
\Gamma\_i^{(0)} = \frac{2}{\hbar^2} \text{Re} \int_0^{+\infty} |M_{\text{e}}^{\text{BO}}|^2 G^{(0)}(\tau) e^{i W_{if} \tau/\hbar - \gamma \tau} d\tau,
``` 
where 
```math
G^{(0)}(\tau) = e^{i \sum_j D_j^{(0)}(\tau)}.
```
 In general, 
 ```math
D_j^{(0)}(\tau) = -2 \left[ iS_j^{-1} A_j(\tau) - {S_j'}^{-1} \cot(\omega_j'\tau/2) \right]^{-1},
```
with 
```math
S_j = \frac{\omega_j}{2\hbar} (\Delta q_j)^2,
```
```math
S_j' = \frac{\omega_j'}{2\hbar} (\Delta q_j)^2,
```
and 
```math
A_j(\tau) = \frac{e^{i\omega_j\tau}(\bar{n}_j + 1) + \bar{n}_j}{e^{i\omega_j\tau}(\bar{n}_j + 1) - \bar{n}_j}.
```
When $\omega_j = \omega_j'$, $D_j^{(0)}(\tau)$ simplifies to 
```math
D_j^{(0)}(\tau) = S_j/i  \left[ (\bar{n}_j + 1)e^{i\omega_j \tau} + \bar{n}_j e^{-i\omega_j \tau} - (2\bar{n}_j + 1) \right].
```

### First-order

The first-order transition rate is 
```math
\Gamma_i^{(1)} = \frac{2}{\hbar^2} \text{Re} \int_0^{+\infty} \left[ \sum_j |M_j|^2 D_j^{(1)}(\tau) \right] G^{(0)}(\tau) e^{iW_{if} \tau/\hbar - \gamma \tau} d\tau,
```
with 
```math
D_j^{(1)} = -\left[ \frac{\Delta q_j}{2 S_j} \right]^2 D_j^{(0)}(\tau) \left( D_j^{(0)}(\tau) [A_j(\tau)]^2 - \frac{1}{2} \frac{\omega_j}{\omega_j'} \left[ A_j(\tau) \cot\frac{\omega_j' \tau}{2}  - \frac{\omega_j \cot(\omega_j'\tau/2) - i\omega_j' A_j(\tau)}{\omega_j A_j(\tau) + i\omega_j' \cot(\omega_j'\tau/2)}\right] \right).
```
When $\omega_j = \omega_j'$, $D_j^{(1)}(\tau)$ simplifies to 
```math
D_j^{(1)}(\tau) = \frac{1}{2} \frac{\hbar}{\omega_j} \left[ (\bar{n}_j + 1)e^{i\omega_j \tau} + \bar{n}_j e^{-i\omega_j \tau} + S_j \left( 1 - (\bar{n}_j + 1)e^{i\omega_j \tau} + \bar{n}_j e^{-i\omega_j \tau} \right)^2 \right].
```

