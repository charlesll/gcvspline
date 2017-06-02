# Welcome to gcvspline's documentation!

## Introduction

gcvspline is a small Python package wrapping the gcvspl.f FORTRAN library, created by H.J. Woltring and available on the netlib.

Code available at https://github.com/charlesll/gcvspline

## Installation

pip install . or to get the latest tagged version:

	pip install gcvspline

or

	conda install -c charlesll gcvspline=0.4

Note that the conda install only works for linux 64 at this time

## Licence information

The gcvspline python wrapper is provided under a GNU GPL3 licence, however please take into account the Licence header in gcvspl.f for full licence information:

From the header of gcvspl.f:

"(C) COPYRIGHT 1985, 1986: H.J. Woltring
This software is copyrighted, and may be  copied  for  exercise,
study  and  use  without authorization from the copyright owner(s), in
compliance with paragraph 16b of  the  Dutch  Copyright  Act  of  1912
("Auteurswet  1912").  Within the constraints of this legislation, all
forms of academic and research-oriented excercise, study, and use  are
allowed,  including  any  necessary modifications.  Copying and use as
object for commercial exploitation are not allowed without  permission
of  the  copyright owners, including those upon whose work the package
is based."

## Disclaimer

gcvspline is provided as is, use at your own risks.

## Available functions

Functions are documented and a call with help() of pydoc will provide all the relevant informations.

- gcvspline and splderivative that are direct wrappers of the FORTRAN functions. See gcvspline_fortran_interface.py for an example.

- Scipy-like interface of several functions: GCVSmoothedNSpline, MSESmoothedNSpline, DOFSmoothedNSpline, SmoothedNSpline. See ./examples/gcvspline_scipylike_interface.ipynb for examples.


## Examples

Examples are available in the example folder of gcvspline.

A short example is the fitting of a spline to noisy data:


	from gcvspline import GCVSmoothedNSpline

	x = np.linspace(-3, 3, 50)
	y0 = np.exp(-x**2)
	np.random.seed(1234)

	n = 0.1 * np.random.normal(size=50)
	w = 1.0 / (np.ones_like(n) * std(n))
	y = y0 + n

	xs = np.linspace(-3, 3, 1000)

	GCV_auto = GCVSmoothedNSpline(x, y, w=w) # gcvspline fit

	y_smoothed = GCV_auto(xs) # retrieving the smoothed values

## Contributors:

Charles Le Losq, the Australian National University (RSES), Canberra.

Yu Feng, University of California, Berkeley.

contact: charles.lelosq@anu.edu.au or yfeng1@berkeley.edu

## References

See:

Woltring, 1986, A FORTRAN package for generalized, cross-validatory spline smoothing and differentiation. Adv. Eng. Softw. 8:104-113. 

and references cited therein.

gcvspl.f code available on www.netlib.org, other versions are available on https://isbweb.org/software/sigproc.html
The gcvspl.f code shipped in this package is a Fortran 77 version downloaded on https://isbweb.org/software/sigproc.html

