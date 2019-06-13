``` 
pull in MjModule

start a timer

call readInputs()

call computeGeneralizedDisplacements()

call computeVariables()

call displaceAtoms()

if ( readQEInput ) then
  call exportQEInput()
else
  call writeNewAtomicPositions()
endif
```
