The first line of the file contains the number of q points (`nOfqPoints`) and atoms (`nAtoms`).
```
<number of q points>  <number of atoms>
```

The file then has a blank line before and array that is 66x4. The first 3 columns belong to the `atomD` array
and the last column belongs to the `atomM` vector.
```

<atomD array> <atomM vector>

```

After that, the file begins looping. For each q point, there will first be a line like
```
     q = (    <phonF column> ) 
```
where the vales in the parentheses go into a column in the `phonF` array. After a blank line, the file has another inner loop
over the modes. For each mode, there will be

```

   <frequency in THz> THz   <?> 2PiTHz <?> cm-1   <?> meV
             X         Y         Z           dx          dy          dz
```
where the first line has the frequency in THz that goes into the `freqInTHz` variable and several other numbers that aren't 
read by the program. The second line contains column labels for the following data. The file then has another inner loop over 
the number of atoms. For each atom, there will be

```
      < X > < Y >  < Z >    < dx >    < dy >   < dz >  
```
where `X`, `Y`, and `Z` are ignored by the program and `dx`,`dy`, and `dz` go into the `phonD` array that holds the phonon
displacements. In each block, there will be `nAtoms` lines holding the positions and displacements of all atoms. There will be
`nModes` blocks where `nModes = 3*nAtoms - 3`. The first line with `q = ( <phonF column> )` and the `nModes` blocks will repeat
`nOfqPoints` times.
