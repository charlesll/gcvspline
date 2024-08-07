from numpy.distutils.core import setup
from numpy.distutils.extension import Extension

extension = Extension("gcvspline._gcvspl", ["gcvspline/gcvspl.f"])

setup(name='gcvspline',
      version='0.5',
      description='A Python wrapper for gcvspl.f, a FORTRAN package for generalized, cross-validatory spline smoothing and differentiation.',
      url='https://github.com/charlesll/gcvspline',
      author='Charles Le Losq and Yu Feng',
      author_email='lelosq@ipgp.fr',
      license='GNU-GPLv3',
      packages=['gcvspline', 'gcvspline.tests'],
      install_requires=['numpy>=1.12'],
      python_requires='>=3.6, <3.12', 
      ext_modules=[extension],
      zip_safe=False)
