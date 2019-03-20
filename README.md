# Fast Cross-validation in Harmonic Approximation


## Overview

The `fct` package is a fast implementation for evaluating the ordinary cross-validation score and the generalized cross-validation score.
 * The `matlab/*/examples/plot_*` scripts create the figures from the paper below
 * The `matlab/*/examples/ad_*` scripts use the combination of the fast evaluation with the a minimization technique to provide a atomatic denoising scheme


## Requirements

 * the equispaced and rank-1 lattice examples on the torus and chebyshev example on the unit interval work out of the box
 * for the other scripts we need the [`nfft`-library](https://www-user.tu-chemnitz.de/~potts/nfft)
   * for Matlab we need the `--with-matlab=PATH_TO_MATLAB` flag
   * for Octave we need the `--with-octave` flag
   * for julia we need the `--with-julia` flag
   * for the scripts on the unit interval we need additionally the `--enable-ndct` flag
   * and for the scripts on the two-dimensional sphere we need the `--enable-nfsft` flag and the [`mtex` toolbox](https://github.com/mtex-toolbox) for the quadrature grid, Voronoi weights, and the plotting


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
