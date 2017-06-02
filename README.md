# gcvspline

Python wrapper of the gcv-spl Fortran library gcvspl.f, created by H.J. Woltring.

Reference: Woltring, 1986, A FORTRAN package for generalized, cross-validatory spline smoothing and differentiation. Adv. Eng. Softw. 8:104-113. 

## Contributors:

Charles Le Losq, the Australian National University (RSES), Canberra.

Yu Feng, University of California, Berkeley.

contact: charles.lelosq@anu.edu.au or yfeng1@berkeley.edu

## Licence information

The gcvspline wrapper is provided under a GNU GPL3 license. The license of the fortran code GCVSPL.f is different, retraining commercial use. See the LICENSE file for further details.

## Disclaimer

gcvspline is provided as is, use at your own risks.

## Requirements

numpy >= 1.12.1 (due to a f2py bug in 1.12.0)

## Installation

pip install . or to get the latest tagged version:

	pip install gcvspline

or

	conda install gcvspline

## Documentation

Documentation is provided at [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://charlesll.github.io/gcvspline/)
