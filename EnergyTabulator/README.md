# Energy Tabulator

There are 4 different energy differences that need to be tabulated:
* Delta function: total energy difference between individually-relaxed charge states, plus additional carrier energy above WZP reference-carrier energy
* Zeroth-order matrix element: same as delta function but using total energy difference between charge states in the relaxed initial positions
* First-order matrix element: eigenvalue difference between initial and final state
* Plotting: eigenvalue difference between initial state and CBM (for electrons) or VBM (for holes)

Instead of calculating these energies in all of the different programs in different ways, this program tabulates all of the required energies so that they can easily be read by other programs. Following the format of the other programs, this code will output an energy table file `energyTable.isp.ik` for each spin (`isp`) and k-point (`ik`). The spins and k-points come from the system that is used for the eigenvalues. 

The inputs should look like
```f90
&inputParams
  ! System used for eigenvalues
  exportDirEigs = 'path-to-system-export-for-eigenvalues'

  ! Systems used for total energy differences
  exportDirInitInit = 'path-to-relaxed-initial-charge-state-export'
  exportDirFinalInit = 'path-to-final-charge-state-initial-positions-export'
  exportDirFinalFinal = 'path-to-relaxed-final-charge-state-export'
  
  ! Parameters for change in energy
  eCorrectTot = real				  ! size of total-energy correction in eV; default 0.0
  eCorrectEigF = real				  ! size of correction to eig diff to final state in eV; default 0.0
  eCorrectEigRef = real			  ! size of correction to eig diff to ref carrier in eV; default 0.0
  refBand = integer						! band location of WZP reference carrier
  CBMorVBMBand = integer      ! band for CBM (for electrons) or VBM (for holes)
  
  ! Band range for overlaps
  iBandIinit = integer						! lowest initial-state band
  iBandIfinal = integer						! highest initial-state band
  iBandFinit = integer						! lowest final-state band
  iBandFfinal = integer						! highest final-state band
  
  ! Output info
  outputDir = 'path-to-store-energy-tables' 			! default './'
/
```
_Note: Do not alter the `&inputParams` or `/` lines at the beginning and end of the file. They represent a namelist and fortran will not recognize the group of variables without this specific format_

The code extracts the total energies from the three different systems and the eigenvalues from the specified eigenvalue system. It is best to use HSE calclations to get the energies, even if you can't do HSE for all of the matrix elements. If you can't use HSE for the energy differences and/or need to include some additional energy correction in the total energy differences, that can be done through `eCorrectTot`. The eigenvalues are best from the ground-state HSE as well, but if that cannot be used (e.g., if many k-points are used and that isn't feasible), PBE-level eigenvalues can be used. Within the conduction band/valence band, PBE eigenvalue differences are okay, but you may need to include corrections to the reference carrier (`eCorrectEigRef`) and final state (`eCorrectEigF`).

The energies are calculated including potential energy corrections as
```f90
dETotWRelax = eTotFinalFinal - eTotInitInit + eCorrectTot
dETotElecOnly = eTotFinalInit - eTotInitInit + eCorrectTot

do ibf = iBandFinit, iBandFfinal
  do ibi = iBandIinit, iBandIfinal
    dEDelta = dETotWRelax - abs(eigvI(ibi) - refEig + eCorrectEigRef)
    dEZeroth = dETotElecOnly - abs(eigvI(ibi) - refEig + eCorrectEigRef)
    dEFirst = abs(eigvI(ibi) - eigvF(ibf) + eCorrectEigF)
    dEPlot = abs(eigvI(ibi) - eigCBMorVBM)
  enddo
enddo
```
where `eigvI` are initial-state eigenvalues indexed by `ibi` and `eigvF` are final-state eigenvalues indexed by `ibf`.
