# Export Program

This program is largely copied from the QE file `pw_export.f90` found in the `PP/src` folder. This program is obsolete in QE, but you can still see a description of its inputs [here](https://www.quantum-espresso.org/Doc/INPUT_pw_export.html). As described on that page, the purpose of the program is 

<blockquote>
   Writes PWSCF data for postprocessing purposes in XML format using IOTK lib. <br/>
   Wave-functions are collected and written using IO_BASE module.
</blockquote>

The QE version heavily uses the `iotk` module ([learn more](http://web.mit.edu/espresso_v6.1/i386_linux26/qe-6.1/iotk/doc/manpages)).

There are some changes made to the QE version to fit our needs, but it is not completely clear what the purpose is. I would eventually like to walk through this program like all of the others, but it includes functions from QE that are incredibly dense. The main gist is that the program pulls the needed information from the QE output xml file and puts it in files in the `exportDir` folder.

# Notes

* Uses output directory instead of just output file
* Uses manual output more than `iotk` module
* Doesn't output some variables that base program does
* _What does the program pull from the QE output files?_
* Some functions are included directly from QE using library archive that is included in the `Makefile`
