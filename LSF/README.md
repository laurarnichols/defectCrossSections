# Line-shape function (`LSF`)

The `LSF` program takes in matrix elements, an energy table, the Huang-Rhys factor, and the phonon frequencies and occupation numbers and performs the time-domain integration to get the final transition rate. Before this code can run, you need to run the `TME` program for the proper-order matrix elements, the `EnergyTabulator` program to get the band bounds and table of energies for the delta function, and the `PhononPP` program that post-processes the phonons to get the Huang-Rhys factors and phonon frequencies and occupations. 

The output files are 
  * `transitionRate.isp.ik`/`transitionRate.isp` -- initial band index and transition rate (1/s) (k-point index not present for scattering)
  * stdout -- timing information and status updates

Every calculation must have the following inputs (don't have to be explicitly set if want default behavior):
* `iSpin` -- default `1`; spin index (`1` or `2`)
* `order` -- no default; what order to calculate (`0` or `1`)
* `dt` -- default `1e-4`; time step size for integration
* `hbarGamma` -- default `0.0`; determines smearing of integrand based on $\exp(-\gamma \tau)$ to achieve convergence
* `smearingExpTolerance` -- default `0.0`; max time for integration is determined by this threshold set on the smearing exponential
* `captured` -- default `.true.`; if the carrier is captured or not; affects how inputs are read and states are treated
* `diffOmega` -- default `.false.`; if there is a different frequency before and after transition; only currently tested for capture; should be consistent with `PhononPP` output
* `newEnergyTable` -- default `.false.`; if an energy table other than the one input to `TME` should be used for matrix elements; only currently implemented for capture
* `oldFormat` -- default `.false.`; if old format of capture matrix elements should be used; this option should not be needed in most cases
* `energyTableDir` -- path to `EnergyTabulator` output directory; needed for $\delta$-function energy and if using new energy table for matrix elements
* `njInput` -- path to single `nj` file to be used for mode occupations
* `outputDir` -- default `'./'`; where to output transition-rate files
* `PhononPPDir` -- path to `PhononPP` output directory with `Sj.out`

The `LSF` input file should look like
```f90
&inputParams
 ...
/
```
with appropriate input variables specified. _Note: Do not alter the `&inputParams` or `/` lines at the beginning and end of the file. They represent a namelist and fortran will not recognize the group of variables without this specific format_

There is also an option in the code to set a threshold for inclusion of modes in the line-shape function based on a threshold for $`S_j`$ using `SjThresh`. I tried to test this but got weird results, so I would not suggest using it in its current form. 

The way the transition rate is calculated and the way the paths are given to the matrix elements differs for the zeroth- and first-order terms, as described below.

## Zeroth-order

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

To specify the path to the zeroth-order matrix elements, simply use `matrixElementDir = 'path-to-matrix-element-files'`.

## First-order

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

The specify the path to the first-order matrix elements, set `MjBaseDir = 'path-to-main-dir-for-first-order-matrix-elements'`. The code assumes that the matrix elements for each mode are in separate sub-folders with names `prefixXXX` with the `prefix` variable set on input and the `XXX` representing the mode index with the length set by the maximum width of the mode indices. The suffix length is set automatically. *Within each folder*, the matrix elements should be located at the path `matrixElementDir` (default `'./'`).

