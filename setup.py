from setuptools import setup

setup(name='gcvspline',
      version='0.2',
      description='A wrapper for gcvspl.f, a FORTRAN package for generalized, cross-validatory spline smoothing and differentiation.',
      url='http://github.com/storborg/funniest',
      author='Charles Le Losq',
      author_email='charles.lelosq@anu.edu.au',
      license='GNU-GPLv3',
      packages=['gcvspline'],
      zip_safe=False)