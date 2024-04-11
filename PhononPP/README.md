# Phonon post-processing

Shifts atom positions by a magnitude of `shift` along each phonon-mode eigenvector direction, and calculate the Huang-Rhys factor. 

**Assumes phonons only at Gamma point!** 

There are two cases that `PhononPP` can be run with: a single displacement (e.g., capture) or a set of displacements based on bands/k-points (e.g., scattering). For a single displacement, the inputs should look like
```
&inputParams
  singleDisp = .true.             ! If there's a single displacement to consider, default .true.

  ! Required to read phonons:
  phononFName      = 'path-to-phonopy-output-yaml-file'        ! Default: './mesh.yaml'
  freqThresh = frequency-threshold       ! Threshold for keeping frequencies, default 0.5 (same units as Phononpy)

  ! Calculating Sj
  calcSj           = logical                                   ! Default: .true.; if Sj should be calculated
  initPOSCARFName  = 'path-to-init-poscar'                     ! Default: './POSCAR_init'
  finalPOSCARFName = 'path-to-final-poscar'                    ! Default: './POSCAR_final'

  ! Calculating dq:
  calcDq  = logical                                   ! Default: .true.; if dq output file should be generated
  dqFName = 'file-to-output-generalized-coord-norms'  ! Default: './dq.txt'

  ! Generating shifted POSCARs:
  generateShiftedPOSCARs = logical                         ! Logical, default .true.
  prefix                 = 'prefix-for-shifted-poscars'    ! Default: './ph_POSCAR'
  
  ! Calculating max displacement:
  calcMaxDisp = logical               ! Default .false.; if max displacement between two atoms across modes should be calculated
  disp2AtomInd     = int1, int2            ! Two integers for indices of atoms to get relative displacement for each mode             

  ! Needed to get dq and shifted POSCARs and calculating max displacement:
  basePOSCARFName = 'path-to-base-poscar'           ! Default: './POSCAR_init'; base to shift from
  shift           = shift-magnitude                 ! Real (Angstrom), default 0.01 A
/
```


For a range of displacements, the inputs should look like
```
&inputParams
  singleDisp = .false.             ! If there's a single displacement to consider, default .true.

  ! Required to read phonons:
  phononFName      = 'path-to-phonopy-output-yaml-file'        ! Default: './mesh.yaml'
  freqThresh = frequency-threshold       ! Threshold for keeping frequencies, default 0.5 (same units as Phononpy)

  ! Calculating Sj
  calcSj           = logical                                   ! Default: .true.; if Sj should be calculated
  CONTCARsBaseDir  = 'path-to-base-dir-for-relaxed-positions'   
  iBandIinit       = integer						! lowest initial-state band
  iBandIfinal      = integer						! highest initial-state band
  iBandFinit       = integer						! lowest final-state band
  iBandFfinal      = integer						! highest final-state band

  ikIinit          = integer                        ! lowest initial-state k-point
  ikIfinal         = integer                        ! highest initial-state k-point
  ikFinit          = integer                        ! lowest final-state k-point
  ikFfinal         = integer                        ! highest final-state k-point

  ! Calculating dq:
  calcDq  = logical                                   ! Default: .true.; if dq output file should be generated
  dqFName = 'file-to-output-generalized-coord-norms'  ! Default: './dq.txt'

  ! Generating shifted POSCARs:
  generateShiftedPOSCARs = logical                         ! Logical, default .true.
  prefix                 = 'prefix-for-shifted-poscars'    ! Default: './ph_POSCAR'

  ! Calculating max displacement:
  calcMaxDisp = logical               ! Default .false.; if max displacement between two atoms across modes should be calculated
  disp2AtomInd     = int1, int2            ! Two integers for indices of atoms to get relative displacement for each mode             

  ! Needed to get dq and shifted POSCARs and calculating max displacement:
  basePOSCARFName = 'path-to-base-poscar'           ! Default: './POSCAR_init'; base to shift from
  shift           = shift-magnitude                 ! Real (Angstrom), default 0.01 A
/
```


To run, use
```
aprun -n num-procs path-to-package/bin/Shifter.x < shifter.in
```
in a PBS/SLURM script. Parallelization is implemented, but is usually not needed. Although you may be able to get away with running without submitting to the queue because the run time is in seconds, the array sizes needed are not insignificant, so it is best practice to submit it to the queue. 

Output:
* `3*nAtoms-3` shifted POSCAR files with name `prefix_j`, where j is the mode number padded by leading zeros
* `Sj.k<iki>_b<ibi>.k<ikf>_b<ibf>.out`, which contains the Huang-Rhys factor and the phonon-mode frequencies
* `dq.txt`, which has the $\delta q\_j$ used in the first-order matrix elements

_Note: The first 3 modes in the yaml file are ignored because those are the acoustic modes and we only want the optical modes._
