title: `input` file (exported from [[pw_export_for_tme(program)]])
author: Laura Nichols
date: 7/3/2019

This file has the following format:
```
# Cell volume (a.u.)^3. Format: '(ES24.15E3)'
  2.653020000000000E+002
# Number of K-points. Format: '(i10)'
         2
# ik, groundState, ngk_g(ik), wk(ik), xk(1:3,ik). Format: '(3i10,4ES24.15E3)'
         1         0       180  5.000000000000000E-001  2.500000000000000E-001  2.500000000000000E-001  2.500000000000000E-001
         2         0       186  1.500000000000000E+000  2.500000000000000E-001  2.500000000000000E-001  7.500000000000000E-001
# Number of G-vectors. Format: '(i10)'
      1459
# Number of PW-vectors. Format: '(i10)'
       266
# Number of min - max values of fft grid in x, y and z axis. Format: '(6i10)'
        -7         7        -7         7        -7         7
# Cell (a.u.). Format: '(a5, 3ES24.15E3)'
# a1  -5.100000000000000E+000  0.000000000000000E+000  5.100000000000000E+000
# a2   0.000000000000000E+000  5.100000000000000E+000  5.100000000000000E+000
# a3  -5.100000000000000E+000  5.100000000000000E+000  0.000000000000000E+000
# Reciprocal cell (a.u.). Format: '(a5, 3ES24.15E3)'
# b1  -6.159985595274104E-001 -6.159985595274104E-001  6.159985595274104E-001
# b2   6.159985595274104E-001  6.159985595274104E-001  6.159985595274104E-001
# b3  -6.159985595274104E-001  6.159985595274104E-001 -6.159985595274104E-001
# Number of Atoms. Format: '(i10)'
         2
# Number of Types. Format: '(i10)'
         1
# Atoms type, position(1:3) (a.u.). Format: '(i10,3ES24.15E3)'
         1  0.000000000000000E+000  0.000000000000000E+000  0.000000000000000E+000
         1  2.550000000000000E+000  2.550000000000000E+000  2.550000000000000E+000
# Number of Bands. Format: '(i10)'
         4
# Spin. Format: '(i10)'
         1
# Element
 Si
# Number of Atoms of this type. Format: '(i10)'
         2
# Number of projectors. Format: '(i10)'
         6
# Angular momentum, index of the projectors. Format: '(2i10)'
         0         1
         0         2
         1         3
         1         4
         2         5
         2         6
# Number of channels. Format: '(i10)'
        18
# Number of radial mesh points. Format: '(2i10)'
      1141       837
# Radial grid, Integratable grid. Format: '(2ES24.15E3)'
  6.513442611103688E-005  8.141803263879609E-007
  6.595371633350159E-005  8.244214541687697E-007
  6.678331195832729E-005  8.347913994790911E-007
  6.762334261151813E-005  8.452917826439768E-007
  6.847393954957285E-005  8.559242443696606E-007
  ...
# AE, PS radial wfc for each beta function. Format: '(2ES24.15E3)'
  5.130578464215461E-004  3.370121859129863E-005
  5.194773426646861E-004  3.412512773645108E-005
  5.259771611021033E-004  3.455436900225269E-005
  5.325583067433258E-004  3.498900945852498E-005
  5.392217087420185E-004  3.542911701872407E-005
  ...
# Fermi Energy (Hartree). Format: '(ES24.15E3)'
  2.090654946391985E-001

```

Currently, the program does not read in cell volume, `groundState`, number of G vectors, the fft grid min and max values, real space or reciprocal space lattice vectors, number of bands, or spin format.
