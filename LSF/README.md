# Line-shape function (`LSF`)

The `LSF` program takes in matrix elements, an energy table, the Huang-Rhys factor, and the phonon frequencies and performs the time-domain integration to get the final transition rate. Before this code can run, you need to run the `TME` program for the proper-order matrix elements, the `EnergyTabulator` program to get the band bounds and table of energies for the delta function, and the `Sj` program that gives the Huang-Rhys factors and the phonon frequencies. 

The output files are 
  * `transitionRate.isp.ik` -- initial band index and transition rate (1/s)
  * stdout -- timing information and status updates

The same code is used for both the zeroth-order and first-order transition rates. 

The inputs should be put in a nameslist like 
```f90
&inputParams
   ...
/
```
within a file (e.g., `LSF.in`) that is passed into the code like `mpirun LSF.x < LSF.in > LSF.out`. The independent time-domain integration steps are divided over the processes. If there are not enough processes, there will be integer overflow when calculating the number of steps per process. 

For both orders, the following inputs are requried:
```f90
&inputParams
 iSpin = integer ! Spin channel of interest
 
 energyTableDir   = 'path-to-energy-tables' 
 matrixElementDir = 'path-to-matrix-element-files' 
 SjInput          = 'path-to-Sj-output/Sj.out'

 temperature          = real    ! (K)
 dt                   = real    ! time step to use (ps) 
 hbarGamma            = real    ! (meV)
 smearingExpTolerance = real    ! tolerance for e^{\gamma t} to determine max time
/
```


## Zeroth-order

The zeroth-order transition rate is $$\Gamma\_i^{(0)} = \frac{2}{\hbar^2} \text{Re} \int\_0^{+\infty} |M_{\text{e}}^{\text{BO}}|^2 G^{(0)}(\tau) e^{i W_{if} \tau/\hbar - \gamma \tau} d\tau,$$ where $$G^{(0)}(\tau) = \exp  \left\( \sum\_j S\_j  \left[ (\bar{n}\_j + 1)e^{i\omega\_j \tau} + \bar{n}\_j e^{-i\omega\_j \tau} - (2\bar{n}\_j + 1) \right] \right\), $$ and $S\_j$ is the Huang-Rhys factor $$S\_j = \frac{\omega\_j}{2\hbar} (\Delta q\_j)^2.$$ 

This option is selected by using `order=0` in the input file.

## First-order

The first-order transition rate is $$\Gamma\_i^{(1)} = \frac{1}{2\hbar^2} \text{Re} \int\_0^{+\infty} \left[ \sum\_j |M\_j|^2 A\_j(\tau) \right] G^{(0)}(\tau) e^{iW_{if} \tau/\hbar - \gamma \tau} d\tau.$$

In addition to the required inputs above, the following inputs must also be included in order to point to all of the different first-order matrix elements:
```f90
&inputParams
 order = 1
 MjBaseDir        = 'path-to-main-dir-for-first-order-matrix-elements'
 prefix           = 'prefix-on-each-subdir-for-matrix-element-calculation'
 matrixElementDir = 'path-to-matrix-element-files-within-each-subdir' 
/
```
It is assumed that the first-order matrix elements are all in a single main folder `MjBaseDir` that holds subfolders `<prefix>xxxx`, where `xxxx` represents a four-digit number indicating the phonon-mode index. Within each of the subfolders, the path to the matrix element (possibly `./`) is given by `matrixElementDir`.

