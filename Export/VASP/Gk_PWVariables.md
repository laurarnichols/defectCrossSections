# $G+k$ and Planewave Variables

This document details the different variables in the `ExportFromVASPOutput` code -- specifically, those related to the G-vectors/plane waves.

## `calculateGvecs`

* `gVecMillerIndicesGlobal_tmp(ix,ig)` -- every possible G-vector in Miller indices (unsorted)
* `gVecMillerIndicesGlobal(ix,ig)` -- every possible G-vector in Miller indices (sorted by magnitude of G-vector)
* `gVecInCart(ix,ig)` -- G-vectors in Cartesian coordinates for only G-vectors on this process

## `distributeGvecsOverProcessors`

* `gIndexLocalToGlobal(ig)` -- index of this local G-vector in the global G-vector array `gVecMillerIndicesGlobal`
* `mill_local(ix,ig)` -- G-vectors in Miller indices for only this process

## `reconstructFFTGrid`
* `xkCart(ix)` -- k-point position in Cartesian coordinates
* `q` -- $|G+k|$ in Cartesian coordinates
* `gkMod(ik,ig)` -- $|G+k|$ for G-vectors on this process where $|G+k|$ is less than the cutoff 
* `gToGkIndexMap(ik,ig)` -- index of G-vector that satisfies $|G+k| <$ cutoff within the local array `gVecInCart`
* `maxNumPWsLocal` -- maximum number of $G+k$ vectors less than the cutoff for only k-points on this process
* `nGkLessECutLocal(ik)` -- number of $G+k$ vectors less than the cutoff for this k-point (only holds k-points within pool)
* `nGkLessECutGlobal(ik)` -- number of $G+k$ vectors less than the cutoff for this k-point (holds all k-points)
* `nPWs1kGlobal(ik)` -- number of $G+k$ vectors less than the cutoff for this k-point, as read from the WAVECAR file (holds all k-points)
* `gKIndexLocalToGlobal(ig,ik)`
  * Index of this local G-vector in the global G-vector array `gVecMillerIndicesGlobal`
  * G-vectors only stored if $G+k$ for this local k-point is less than the cutoff
  * In order of increasing $|G+k|$
* `maxGIndexLocal` -- max index of local G-vectors in global G-vector array `gVecMillerIndicesGlobal` s.t. $|G+k| < $ cutoff
* `maxGIndexGlobal` -- max index of all G-vectors in global G-vector array `gVecMillerIndicesGlobal` s.t. $|G+k| < $ cutoff
* `maxNumPWsGlobal` -- max number of $G+k$ vectors less than the cutoff across all k-points

## `getGlobalGkIndices`
* `itmp1` -- holds the index in spots where $|G+k| <$ cutoff and zero elsewhere
* `gkIndexGlobal` -- all indices of G-vectors in `gVecMillerIndicesGlobal` s.t. $|G+k| <$ cutoff (`itmp1` with zeros removed)
