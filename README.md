# Fast Cross-validation in Harmonic Approximation

## Todo

[ ] update README
[ ] automatic denoising examples for the unit interval
[ ] proper install/uninstall scripts


## Overview

The `fct` package is a fast implementation for evaluating the ordinary cross-validation score and the generalized cross-validation score.
 * In the `matlab/examples` folder you can find the scripts for creating the figures from the paper below
 * In the `matlab/smoothn` folder you can find functions which combine the fast evaluation with the simple bisection method for automated soothing.


## Requirements

 * `equispaced_1d.m`, `equispaced_2d.m` and `r1l.m` on the torus and `equispaced.m` on the unit interval work out of the box
 * for the other scripts we need the [`nfft`-library](https://www-user.tu-chemnitz.de/~potts/nfft)
   * for Matlab we need the `--with-matlab=PATH_TO_MATLAB` flag
   * for Octave we need the `--with-octave` flag
   * for julia we need the `--with-julia` flag
   * for the scripts on the unit interval we need additionally the `--enable-ndct` flag
   * and for the scripts on the two-dimensional sphere we need the `--enable-nfsft` flag and the [`mtex` toolbox](https://github.com/mtex-toolbox) for the quadrature grid and the plotting


## Citing

If you use `fcv` in your work, please cite the following.

```tex
@article{,
author  = {},
title   = {},
journal = {},
year    = {},
volume  = {},
number  = {},
pages   = {},
doi     = {}
}
```
