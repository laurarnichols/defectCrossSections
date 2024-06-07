# Line-shape function (`LSF`)

The `LSF` program takes in matrix elements, an energy table, the Huang-Rhys factor, and the phonon frequencies and performs the time-domain integration to get the final transition rate. Before this code can run, you need to run the `TME` program for the proper-order matrix elements, the `EnergyTabulator` program to get the band bounds and table of energies for the delta function, and the `Sj` program that gives the Huang-Rhys factors and the phonon frequencies. 

The output files are 
  * `transitionRate.isp.ik` -- initial band index and transition rate (1/s)
  * stdout -- timing information and status updates

The same code is used for both the zeroth-order and first-order transition rates. The `diffOmega` option allows for calculation of the line-shape function with different initial and final frequencies. With `diffOmega = .true.`, the PhononPP output file `Sj.out` should be tabulated using the `diffOmega = .true.` option as well.

## Zeroth-order

The zeroth-order transition rate is $$\Gamma\_i^{(0)} = \frac{2}{\hbar^2} \text{Re} \int\_0^{+\infty} |M_{\text{e}}^{\text{BO}}|^2 G^{(0)}(\tau) e^{i W_{if} \tau/\hbar - \gamma \tau} d\tau,$$ where $$G^{(0)}(\tau) = e^{iD_j^{(0)}(\tau)}.$$ In general, $$D_j^{(0)}(\tau) = -2 \left[ iS_j^{-1} A_j(\tau) - {S_j'}^{-1} \cot(\omega_j'\tau/2) \right]^{-1},$$ with $$S\_j = \frac{\omega\_j}{2\hbar} (\Delta q\_j)^2,$$ $$S\_j' = \frac{\omega\_j'}{2\hbar} (\Delta q\_j)^2,$$ and $$A_j(\tau) = \frac{e^{i\omega_j\tau}(\bar{n}_j + 1) + \bar{n}_j}{e^{i\omega_j\tau}(\bar{n}_j + 1) - \bar{n}_j}.$$ When $\omega_j = \omega_j'$, $D_j^{(0)}(\tau)$ simplifies to $$D_j^{(0)}(\tau) = S\_j/i  \left[ (\bar{n}\_j + 1)e^{i\omega\_j \tau} + \bar{n}\_j e^{-i\omega\_j \tau} - (2\bar{n}\_j + 1) \right]$$



The `LSF` input file for the zeroth-order should look like
```f90
&inputParams

 order = 0

 iSpin = integer ! Spin channel of interest
 
 energyTableDir   = 'path-to-energy-tables' 
 matrixElementDir = 'path-to-matrix-element-files' 
 SjInput          = 'path-to-Sj-output/Sj.out' 

 outputDir = 'path-to-store-transition-rates'

 temperature          = real    ! (K)
 dt                   = real    ! time step to use (ps) 
 hbarGamma            = real    ! (meV)
 smearingExpTolerance = real    ! tolerance for e^{\gamma t} to determine max time
/
```

_Note: Do not alter the `&inputParams` or `/` lines at the beginning and end of the file. They represent a namelist and fortran will not recognize the group of variables without this specific format_

## First-order

The first-order transition rate is $$\Gamma\_i^{(1)} = \frac{2}{\hbar^2} \text{Re} \int\_0^{+\infty} \left[ \sum\_j |M\_j|^2 D\_j^{(1)}(\tau) \right] G^{(0)}(\tau) e^{iW_{if} \tau/\hbar - \gamma \tau} d\tau,$$ with $$D_j^{(1)} = -\left[ \frac{\Delta q_j}{2 S_j} \right]^2 D_j^{(0)}(\tau) \left( D_j^{(0)}(\tau) [A_j(\tau)]^2 - \frac{1}{2} \frac{\omega_j}{\omega_j'} \left[ A_j(\tau) \cot\frac{\omega_j' \tau}{2}  - \frac{\omega_j \cot(\omega_j'\tau/2) - i\omega_j' A_j(\tau)}{\omega_j A_j(\tau) + i\omega_j' \cot(\omega_j'\tau/2)}\right] \right).$$ When $\omega_j = \omega_j'$, $D_j^{(1)}(\tau)$ simplifies to $$D_j^{(1)}(\tau) = \frac{1}{2} \frac{\hbar}{\omega_j} \left[ (\bar{n}\_j + 1)e^{i\omega\_j \tau} + \bar{n}\_j e^{-i\omega\_j \tau} + S_j \left( 1 - (\bar{n}\_j + 1)e^{i\omega\_j \tau} + \bar{n}\_j e^{-i\omega\_j \tau} \right)^2 \right].$$

The `LSF` input file for the first-order should look like
```f90
&inputParams

 order = 1

 iSpin = integer ! Spin channel of interest
 
 energyTableDir   = 'path-to-energy-tables' 
 matrixElementDir = 'path-to-matrix-element-files' 
 MjBaseDir        = 'path-to-main-dir-for-first-order-matrix-elements'
 prefix           = 'prefix-on-each-subdir-for-matrix-element-calculation'
 SjInput          = 'path-to-Sj-output/Sj.out' 

 outputDir = 'path-to-store-transition-rates'

 temperature          = real    ! (K)
 dt                   = real    ! time step to use (ps) 
 hbarGamma            = real    ! (meV)
 smearingExpTolerance = real    ! tolerance for e^{\gamma t} to determine max time

/
```

_Note: Do not alter the `&inputParams` or `/` lines at the beginning and end of the file. They represent a namelist and fortran will not recognize the group of variables without this specific format_

