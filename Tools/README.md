# Tools

* `kpointBinningByEnergy` -- Take in eigenvalue files for each k point and energy bins from the VfisVsE file, bin each band for each kpoint by energy, and output the 
results
* `PARCHG_sp_splitter.py` -- Split spin channels in partial charge density files (from Andy)
* `PARCHG_sp_fixer.py`
  * From Andy
  * Most of the time, running the splitter is sufficient
  * But you may get an error message that it ran into a bunch of \*\*\*\*\*'s which happens due to a formatting bug in VASP
  * Then you can run the fixer script which basically replaces the \*\*\*\*\*'s by the average of the grid points around it 
  * Obviously this is not 100% accurate but it allows the splitter to run since everything is a number again
  * 1 random data point being slightly off won't affect the visualization that we're after
