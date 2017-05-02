from gcvspline.gcvspline2 import GCVSpline

import numpy as np
from numpy.testing import assert_allclose, assert_raises

def test_spline_cv(): 
    def f(x):
        return 4 * x ** 3 + 10 * x ** 2 + 4 * x + 3
    x = np.linspace(-10, 10, 40)
    y = f(x)

    spl = GCVSpline(x, y, kind='cubic', mode=2)
    y1 = spl(x)

    # since the function is cubic, cubic spl shall work perfectly.
    assert_allclose(y1, y)

def test_spline_mode1(): 
    def f(x):
        return 4 * x ** 3 + 10 * x ** 2 + 4 * x + 3
    x = np.linspace(-10, 10, 40)
    y = f(x)

    spl = GCVSpline(x, y, kind='cubic', mode=1, prior=0.001)
    y1 = spl(x)

    # since the function is cubic, cubic spl shall work perfectly.
    assert_allclose(y1, y, rtol=1e-3)

def test_spline_mode3(): 
    def f(x):
        return 4 * x ** 3 + 10 * x ** 2 + 4 * x + 3
    x = np.linspace(-10, 10, 40)
    y = f(x)

    spl = GCVSpline(x, y, kind='cubic', mode=3, prior=0.001)
    y1 = spl(x)

    # since the function is cubic, cubic spl shall work perfectly.
    assert_allclose(y1, y, rtol=1e-3)

def test_spline_raise(): 
    def f(x):
        return 4 * x ** 3 + 10 * x ** 2 + 4 * x + 3
    x = np.linspace(-10, 10, 40)
    y = f(x)

    spl = GCVSpline(x, y, kind='cubic', mode=2, bbox=(-1, 1), ext=2)

    assert_raises(ValueError, spl, x)

def test_spline_zeros(): 
    def f(x):
        return 4 * x ** 3 + 10 * x ** 2 + 4 * x + 3
    x = np.linspace(-10, 10, 40)
    y = f(x)

    spl = GCVSpline(x, y, kind='cubic', mode=2, bbox=(0, 0), ext=1)

    y1 = spl(x)
    assert_allclose(y1, 0)

def test_spline_consts():
    def f(x):
        return 4 * x ** 3 + 10 * x ** 2 + 4 * x + 3
    x = np.linspace(-10, 10, 40)
    y = f(x)

    spl = GCVSpline(x, y, kind='cubic', mode=2, bbox=(0, 1), ext=3)

    y1 = spl(x)
    assert_allclose(y1[x < 0], spl(0))
    assert_allclose(y1[x > 1], spl(1))
