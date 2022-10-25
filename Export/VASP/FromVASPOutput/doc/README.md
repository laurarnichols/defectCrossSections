# Export From VASP Output

This program exports the data from VASP output files into a form that can be processed by the matrix element code (`TME`). The program needs to know where the output files are and where to put the exported data.

The input file should look like
```
&inputParams
  VASPDir = 'path-to-VASP-output'
  exportDir = 'path-to-put-exported-files'
/
```

The default values are `VASPDir = './'` and `exportDir = './Export'`. If the directory used for the exported files does not exist, it will be created.

Currently, the I/O is only done by the root node (`ionode`), but the processing of the data is split across all processes.
