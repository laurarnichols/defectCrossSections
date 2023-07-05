# Line-shape function (`LSF`)

The `LSF` program takes in matrix elements, an energy table, the Huang-Rhys factor, and the phonon frequencies and performs the time-domain integration to get the final transition rate. Before this code can run, you need to run the `TME` program for the proper-order matrix elements, the `EnergyTabulator` program to get the table of energies for plotting and for the delta function, and the `Sj` program that gives the Huang-Rhys factors and the phonon frequencies. 

The output files are 
  * `transitionRate.txt` -- initial band index, energy to plot (eV), and transition rate (1/s)
  * stdout -- timing information and status updates

The same code is used for both the zeroth-order and first-order transition rates. 

## Zeroth-order

The zeroth-order transition rate is $$\Gamma\_i^{(0)} = \frac{2}{\hbar^2} \text{Re} \int\_0^{+\infty} |M_{\text{e}}^{\text{BO}}|^2 G^{(0)}(t) e^{iE_{if}^{\text{DFT},if}t/\hbar - \gamma t} dt,$$ where $$G^{(0)}(t) = \exp  \left\( \sum\_j S\_j  \left[ (\bar{n}\_j + 1)e^{i\omega\_j t} + \bar{n}\_j e^{-i\omega\_j t} - (2\bar{n}\_j + 1) \right] \right\) $$ and $S\_j$ is the Huang-Rhys factor $$S\_j = \frac{\omega\_j}{2\hbar} (\Delta q\_j)^2.$$

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

 iBandIinit   = integer         ! initial state initial band
 iBandIfinal  = integer         ! initial state final band
 iBandFinit   = integer         ! final state initial band
 iBandFfinal  = integer         ! final state final band

/
```

_Note: Do not alter the `&inputParams` or `/` lines at the beginning and end of the file. They represent a namelist and fortran will not recognize the group of variables without this specific format_

## First-order

The first-order transition rate is $$\Gamma\_i^{(1)} = \frac{1}{2\hbar^2} \text{Re} \int\_0^{+\infty} \left[ \sum\_j |M\_j|^2 A\_j(t) \right] G^{(0)}(t) e^{iE_{if}^{\text{DFT},if}t/\hbar - \gamma t} dt.$$

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

 iBandIinit   = integer         ! initial state initial band
 iBandIfinal  = integer         ! initial state final band
 iBandFinit   = integer         ! final state initial band
 iBandFfinal  = integer         ! final state final band

/
```

_Note: Do not alter the `&inputParams` or `/` lines at the beginning and end of the file. They represent a namelist and fortran will not recognize the group of variables without this specific format_

