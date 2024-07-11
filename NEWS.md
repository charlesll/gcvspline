# gcvspline: summary of the versions

The changes performed between the different versions are references briefly below.

# v0.5

- Add a minimal pyproject.toml file

# v0.4

- Fix of a documentation issue under Python 2.7.

# v0.3

- major change in API 

- gcvspline and splderivative that are direct wrappers of the FORTRAN functions. See gcvspline_fortran_interface.py for an example.

- Scipy-like interface of several functions: GCVSmoothedNSpline, MSESmoothedNSpline, DOFSmoothedNSpline, SmoothedNSpline. See ./examples/gcvspline_scipylike_interface.ipynb for examples.

# v0.2

- gcvspline and splderivative are direct wrappers of the FORTRAN functions. See example.py for an example.

- Wrap of the GCVSPL.f code with f2py, broken with new Python versions.

# v0.1

- Registration of gcvspline in PyPi. 