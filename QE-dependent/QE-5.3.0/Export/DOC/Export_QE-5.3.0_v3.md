# Export Program

This program is largely copied from the QE file `pw_export.f90` found in the `PP/src` folder. This program is obsolete in QE, but you can still see a description of its inputs [here](https://www.quantum-espresso.org/Doc/INPUT_pw_export.html). As described on that page, the purpose of the program is 

<blockquote>
   Writes PWSCF data for postprocessing purposes in XML format using IOTK lib. <br/>
   Wave-functions are collected and written using IO_BASE module.
</blockquote>

The QE version heavily uses the `iotk` module ([learn more](http://web.mit.edu/espresso_v6.1/i386_linux26/qe-6.1/iotk/doc/manpages)).

There are a few changes to the QE version that tailor the program more to our specific needs. The changes are detailed below:
<table>
   <tr>
      <th>Change</th>
      <th>Explanation</th>
   </tr>
   <tr>
      <td>Add <code>exportDir</code> variable</td>
      <td></td>
   </tr>
   <tr>
      <td>Add <code>ig</code> variable</td>
      <td></td>
   </tr>
   <tr>
      <td></td>
      <td></td>
   </tr>
</table>
