# gcvspline

https://travis-ci.org/charlesll/gcvspline.svg?branch=master

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

Installation through pip is recommended:

	pip install gcvspline

Pip wheels for Python 2.7, 3.4 and 3.6 are provided for Windows users (only for Python 3.6 for 32 bit systems). Mac OS and Linux version are built from source and requires gfortran.

Another option for Linux 64 bits is to us Conda: 

	conda install -c charlesll gcvspline=0.4

If the installation fails and this seems related to a problem with FORTRAN compilation, please check the status of your FORTRAN compiler.

The fastest way will be to upload any fortran code and try building it. 

OSX Sierra and High Sierra may run into problems with the assembler in some case, fixed by adding the line

export PATH="/usr/bin/$PATH"

in your .bash_profile file.

## Documentation

Documentation is provided at [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://charlesll.github.io/gcvspline/)
