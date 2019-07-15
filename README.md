![](https://github.com/felixbartel/fcv/raw/master/fcv.png "fcv logo")
# Fast Cross-validation in Harmonic Approximation

## About

The `fcv` package is a fast implementation for evaluating the ordinary cross-validation score and the generalized cross-validation score in various settings.
 * The `demos/demo_*` scripts create the figures from the paper below
 * The `demos/demo_ad_*` scripts use the combination of the fast evaluation with the a minimization technique to provide a atomatic denoising scheme


## Requirements

 * the equispaced and rank-1 lattice examples on the torus and chebyshev example on the unit interval work out of the box
 * for the other scripts we need the [`nfft`-library](https://www-user.tu-chemnitz.de/~potts/nfft)
   * for Matlab we need the `--with-matlab=PATH_TO_MATLAB` flag
   * for Octave we need the `--with-octave` flag
   * for the scripts on the unit interval we need additionally the `--enable-ndct` flag
   * for the scripts on the two-dimensional sphere we need the `--enable-nfsft` flag and the [`mtex` toolbox](https://github.com/mtex-toolbox) for the quadrature grid, Voronoi weights, and the plotting
   * for the scripts on the rotation group we need the `--enable-nfsoft` flag and the [`mtex` toolbox](https://github.com/mtex-toolbox) for the demos

## Installation

* after you have all libraries just run `fcv_install('t')`, `fcv_install('i')`, `fcv_install('s2')`, or `fcv_install('so3')`
* if you didn't install the libraries in `~/repo/` you additionally have to specify theirs paths with `fcv_install('s2',nfft_path,mtex_path)`

## Citing

If you use `fcv` in your work, please cite the following:

```tex
@article{,
author = {Felix Bartel and Ralf Hielscher and Daniel Potts},
title = {Fast Cross-validation in Harmonic Approximation},
year = {2019},
eprint = {arXiv:1903.10206},
}
```
