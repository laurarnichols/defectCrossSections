# Energy Tabulator

There are 4 different energy differences that need to be tabulated:
* Delta function: total energy difference between relaxed charge states, plus additional carrier energy above WZP reference-carrier energy
* Zeroth-order matrix element: same as delta function but using total energy difference between charge states in the initial positions
* First-order matrix element: eigenvalue difference between initial and final state
* Plotting: same as zeroth-order but positive and in eV

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

  ! Physical problem details
  captured = .true. or .false.      ! If the carrier is captured
  elecCarrier = .true. or .false.   ! If the carrier is an electron (vs a hole)

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
  energyTableDir = 'path-to-store-energy-tables' 			! default './'
/
```
_Note: Do not alter the `&inputParams` or `/` lines at the beginning and end of the file. They represent a namelist and fortran will not recognize the group of variables without this specific format_

The code extracts the total energies from the three different systems (`exportDirInitInit`, `exportDirFinalInit`, and `exportDirFinalFinal`) and the eigenvalues from the specified eigenvalue system (`exportDirEigs`). It is best to use HSE calclations to get the total energies, even if you can't do HSE for all of the matrix elements. If you can't use HSE for the energy differences and/or need to include some additional energy correction in the total energy differences, that can be done through `eCorrectTot`. The eigenvalues are best from the ground-state HSE as well, but if that cannot be used (e.g., if many k-points are used and that isn't feasible), PBE-level eigenvalues can be used. Ground-state eigenvalues should always be used because those are most accurate. Within the conduction band/valence band, PBE eigenvalue differences are okay, but you may need to include corrections to the reference carrier (`eCorrectEigRef`) (e.g., if you reference carrier is on the other side of the gap and the band gap from the PBE level needs to be corrected.

The first-order term uses an eigenvalue difference between the initial and final states. For capture, this includes the distance between the band states and the defect level, but this distance changes in the different charge states. Additionally, in finite supercell sizes there is a ficticious dispersion of the defect level. To address both of these issues, the distance between the lowest conduction band and the defect level at $\Gamma$ in the *initial-state HSE calculation* (`exportDirInitInit`) is used, and the lowest conduction band at $\Gamma$ is used as a reference point for the rest of the bands across all k-points. The same reference is used for the delta-function and zeroth-order eigenvalue differences. The order for eigenvalue differences is switched for electron carriers and hole carriers to represent the energy change of the actual electron making the transition.

*Note:* The code currently uses the first k-point as a reference and does not test that it is the $\Gamma$ point. There is also no current way to correct the energy between the defect and the lowest conduction band for the first-order term.


The energies are calculated including potential energy corrections as
```f90
dETotWRelax = eTotFinalFinal - eTotInitInit + eCorrectTot
dETotElecOnly = eTotFinalInit - eTotInitInit + eCorrectTot

do ibf = iBandFinit, iBandFfinal
  do ibi = iBandIinit, iBandIfinal
    if(elecCarrier) then
      dEEigRef = eigvI(ibi) - refEig + eCorrectEigRef
    else
      dEEigRef = refEig - eigvI(ibi) + eCorrectEigRef
    endif

    dEDelta = dETotWRelax - dEEigRef
    dEZeroth = dETotElecOnly - dEEigRef

    if(captured) then
      dEFirst = dEEigRef + dEEigRefDefect
    else
      dEFirst = abs(eigvI(ibi) - eigvF(ibf))
    endif

    dEPlot = dEPlot = abs(dEZeroth)/eVToHartree
  enddo
enddo
```
where `eigvI` are initial-state eigenvalues indexed by `ibi`, `eigvF` are final-state eigenvalues indexed by `ibf`, and `dEEigRefDefect` is the distance between the defect and the lowest conduction band at $\Gamma$ from the initial-state HSE calculation (used for capture only).
