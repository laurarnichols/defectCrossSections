# Version History

I was not the best in the past at having specific demarcations of working versions before and after different features, so the version history isn't as straightforward as I would like. However, I want to include this information here in case it is helpful for future testing.

* VASP Export after `wfc` migration, original TME -- [29483c2](https://github.com/laurarnichols/defectCrossSections/tree/29483c231924d3829a2447eba9e9627b2a9f30f1) (10/28/2022)
  * Still separate codes
  * Projectors and projections written from `ExportFromSrc`, wave functions and others from `ExportFromOutput`
  * Must have a copy of VASP source to use for `ExportFromSrc` compilation
  * Make sure to re-run `make initialize` in the main directory to get the `make.sys` file set up with the needed paths
  * Copy over original VASP source to VASP directory used for `ExportFromSrc` compilation
  * Exported files have old format, so will need to update if want to use with newer version of `TME`
* Original VASP Export and TME-- [0ba39ed](https://github.com/laurarnichols/defectCrossSections/tree/0ba39edb8a74b204f7d25e2a194288f46837bb99) (10/25/2022)
  * Still separate codes
  * Wave functions, projectors, and projections written from `ExportFromSrc`, others from `ExportFromOutput`
  * Must have a copy of VASP source to use for `ExportFromSrc` compilation
  * Make sure to re-run `make initialize` in the main directory to get the `make.sys` file set up with the needed paths
  * Copy over original VASP source to VASP directory used for `ExportFromSrc` compilation
  * Exported files have old format, so will need to update if want to use with newer version of `TME`