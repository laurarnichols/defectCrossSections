# H Release

The hydrogen release problem builds on top of the capture formalism and code. It treats every small step out of the H as a capture problem with a defined "relaxation" and steps the H out until it is released. The final time before release gives the average lifetime ($\tau$) of a H in the given defect, which can be directly related to the number of activated defects over time $$N_d(t) = N_d^{\infty} (1 - e^{-t/\tau}).$$

## Algorithm Flow

<div align="center">
  <img width="50%" src="https://github.com/laurarnichols/defectCrossSections/assets/32521892/fcd9dadd-4359-471c-b118-3ec3bca9bf25" alt="H release flow chart">
</div>
<br/>

First, create a new input geometry that only includes an infinitesimal displacement of the H in the direction of release. The capture code projects that displacement onto each one of the phonon eigenvectors and defines the Huang-Rhys factor $$S_j = \frac{\omega_j}{2\hbar} (\Delta q_j)^2,$$ which determines how strongly each mode contributes to the multiphonon energy dissipation. Get the average energy transfer rate for a given initial state from the capture code, then integrate over initial states to get the total energy transfer rate: $$\frac{dE}{dt} = \int n(E_i) P(E_i) dE.$$ Use $S_j$ as a weight to split up the total energy transfer rate into transfer rates into each of the phonon modes, $dE_j/dt$. Using $\hbar \omega_j$ for each mode, convert the energy transfer rates to a rate of change of the occupation number for each mode, $d\bar{n}\_j/dt$. Using maps of the total energy as a function of each phonon generalized coordinate $q_j$, convert $d\bar{n}\_j/dt$ to $d(\Delta q_j)/dt$. Combine all of the components to get a Cartesian-space displacement rate, then determine the time step needed to achieve the initial H displacement. The new Cartesian-space displacement will now include movement of some of the surrounding atoms. Iterate this process until the displacement of those atoms does not significantly change from one step to another. Once self-consistency is reached, repeat the whole process with a new infinitesimal displacement of H until it is released.

## Todos

- [ ] Add sum over final electronic states to get total scattering rate for carrier in initial state (`Sigma`)
- [ ] Add option to include eigenvalue energy difference in sum over final electronic states to get average energy transfer rate (`Sigma`)
- [ ] Figure out what code needs to be repeated from the capture process and what can be a tabulated result
