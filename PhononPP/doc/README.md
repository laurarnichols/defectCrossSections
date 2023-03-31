# Phonon post-processing

## Shifter

Shifts atom positions by a magnitude of `shift` along each phonon-mode eigenvector direction. **Assumes phonons only at Gamma point!** 

Input should look like
```
&inputParams
  poscarFName = 'path-to-input-poscar'               ! Default: './POSCAR'
  phononFName = 'path-to-phonopy-output-yaml-file'   ! Default: './mesh.yaml'
  prefix = 'prefix-for-shifted-poscars'              ! Default: './ph_POSCAR'
  nAtoms = n-atoms                                   ! Integer
  shift = shift-magnitude                            ! Real (Angstrom), default 0.01 A
/
```

Run time is seconds, so you should be able to use
```
path-to-package/bin/Shifter.x < shifter.in
```
without submitting to the queue, but parallelization is also implemented in case there is a large number of modes, so you can also use
```
aprun -n num-procs path-to-package/bin/Shifter.x < shifter.in
```
in a PBS/SLURM script.

Output is `3*nAtoms-3` shifted POSCAR files with name `prefix_j`, where j is the mode number padded by leading zeros. 

_Note: The first 3 modes in the yaml file are ignored because those are the acoustic modes and we only want the optical modes._
