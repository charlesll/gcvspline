from gcvspline.gcvspline2 import SmoothedNSpline, GCVSmoothedNSpline, MSESmoothedNSpline, DOFSmoothedNSpline

import numpy as np
from numpy.testing import assert_allclose, assert_raises

def test_spline_cv(): 
    def f(x):
        return 4 * x ** 3 + 10 * x ** 2 + 4 * x + 3
    x = np.linspace(-10, 10, 40)
    y = f(x)

    spl = GCVSmoothedNSpline(x, y, kind='cubic')
    y1 = spl(x)

    # since the function is cubic, cubic spl shall work perfectly.
    assert_allclose(y1, y)

def test_spline_fixed(): 
    def f(x):
        return 4 * x ** 3 + 10 * x ** 2 + 4 * x + 3
    x = np.linspace(-10, 10, 40)
    y = f(x)

    spl = SmoothedNSpline(x, y, kind='cubic', p=0.001)
    y1 = spl(x)

    # since the function is cubic, cubic spl shall work perfectly.
    assert_allclose(y1, y, rtol=1e-3)

def test_spline_mse(): 
    def f(x):
        return 4 * x ** 3 + 10 * x ** 2 + 4 * x + 3
    x = np.linspace(-10, 10, 40)
    y = f(x)

    spl = MSESmoothedNSpline(x, y, kind='cubic', variance=0.001)
    y1 = spl(x)

    # since the function is cubic, cubic spl shall work perfectly.
    assert_allclose(y1, y, rtol=1e-3)

def test_spline_dof(): 
    def f(x):
        return 4 * x ** 3 + 10 * x ** 2 + 4 * x + 3
    x = np.linspace(-10, 10, 40)
    y = f(x)

    spl = DOFSmoothedNSpline(x, y, kind='cubic', dof=3)
    y1 = spl(x)

    # since the function is cubic, cubic spl shall work perfectly.
    assert_allclose(y1, y, rtol=1e-3)

def test_spline_raise(): 
    def f(x):
        return 4 * x ** 3 + 10 * x ** 2 + 4 * x + 3
    x = np.linspace(-10, 10, 40)
    y = f(x)

    spl = GCVSmoothedNSpline(x, y, kind='cubic',  bbox=(-1, 1), ext=2)

    assert_raises(ValueError, spl, x)

def test_spline_zeros(): 
    def f(x):
        return 4 * x ** 3 + 10 * x ** 2 + 4 * x + 3
    x = np.linspace(-10, 10, 40)
    y = f(x)

    spl = GCVSmoothedNSpline(x, y, kind='cubic',  bbox=(-1, 1), ext=1)

    y1 = spl(x)
    assert_allclose(y1[x < -1], 0)
    assert_allclose(y1[x > 1], 0)

def test_spline_consts():
    def f(x):
        return 4 * x ** 3 + 10 * x ** 2 + 4 * x + 3
    x = np.linspace(-10, 10, 40)
    y = f(x)

    spl = GCVSmoothedNSpline(x, y, kind='cubic',  bbox=(-1, 1), ext=3)

    y1 = spl(x)
    assert_allclose(y1[x < -1], spl(-1))
    assert_allclose(y1[x > 1], spl(1))
