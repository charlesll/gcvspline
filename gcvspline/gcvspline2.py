# -*- coding: utf-8 -*-
"""
scipy.interpolate style wrapper for gcvspl.

"""

__author__ = "Yu Feng"
__email__ = 'rainwoodman@gmail.com'

from . import _gcvspl as gcvspl
import numpy as np

PARAMETERS= \
"""
            w: array_like
                inverted standard derivation of observations, e.g. 1.0 / std(y, axis=1)
            w1: array_like
                inverted standard derivation of datasets, e.g. 1.0 / std(y, axis=0)
            kind : int, or str
                kind / order of spline, NOTE: different from the internal gcvspline parameter.
            nc  : int
                number of control points.
            bbox : tuple, or None
                bounding box of valid region for evaluation. default is None, (x[0], x[-1])
            ext : int or str
                extrapolation mode.
"""
class GCVSplineBase(object):
    SPLORDER = {'linear' : 1, 'cubic' : 2, 'quintic' : 3, 'heptic': 4, 1 : 1, 3: 2, 5 :3, 7:4}
    EXT = {'extrapolate': 0, 'zeros' : 0, 'raise' : 2, 'const' : 3}
    def __init__(self, x, y, w=1, w1=1, MD=2, VAL=1, nc=None, kind='cubic', bbox=None, ext=0):
        """ Parameters
            ----------
            %s
            MD: int
                   |MD| = 1: Prior given value for p in VAL
                             (VAL.ge.ZERO). This is the fastest
                             use of GCVSPL, since no iteration
                             is performed in p.
                   |MD| = 2: Generalized cross validation.
                   |MD| = 3: True predicted mean-squared error,
                             with prior given variance in VAL.
                   |MD| = 4: Prior given number of degrees of
                             freedom in VAL (ZERO.le.VAL.le.N-M).
                    MD  < 0: It is assumed that the contents of
                            X, W, M, N, and WK have not been
                             modified since the previous invoca-
                             tion of GCVSPL. If MD < -1, WK(4)
                             is used as an initial estimate for
                             the smoothing parameter p.  At the
                             first call to GCVSPL, MD must be > 0.
                   Other values for |MD|, and inappropriate values
                   for VAL will result in an error condition, or
                   cause a default value for VAL to be selected.
                   After return from MD.ne.1, the same number of
                   degrees of freedom can be obtained, for identical
                   weight factors and knot positions, by selecting
                   |MD|=1, and by copying the value of p from WK(4)
                   into VAL. In this way, no iterative optimization
                   is required when processing other data in Y.
            VAL: float
                parameter controls the behavior of the prior, see MD.

        """

        if not (x[:-1] < x[1:]).all():
            raise ValueError("x must be strictly increasing")

        if len(y.shape) > 2:
            raise ValueError("y must be at most two dimensional (len(x), ndataset)")

        if len(y.shape) == 1:
            y = y.reshape(-1, 1)

        if bbox is None:
            bbox = (x[0], x[-1])

        wx = np.ones(y.shape[0])
        wx[...] = w ** 2
        wy = np.ones(y.shape[1])
        wy[...] = w1 ** 2

        if nc is None: nc = len(x)

        splorder = self.SPLORDER.get(kind)

        if splorder not in (1, 2, 3, 4):
            raise ValueError("splorder must be 1, 2, 3, 4.")

        if MD not in (1, 2, 3, 4):
            raise ValueError("MD must be 1, 2, 3, 4")

        ext = self.EXT.get(ext, ext)
        if ext not in (0, 1, 2, 3):
            raise ValueError("exts must be 1, 2, 3, 4")

        c, wk, ier = gcvspl.gcvspl(x, y, wx, wy, splorder, MD, VAL, nc)
        if ier == 1:
            raise RuntimeError(" M <= 0 or N < 2 * M ")
        elif ier == 2:
            raise RuntimeError("Knot sequence is not strictly increasing"
                    + " or some weight factor is not positive.")
        elif ier == 3:
            raise RuntimeError("Wrong mode parameter or value")

        self.x = x
        self.c = c
        self.wk = wk
        self.splorder = splorder
        self.L = 1
        self.bbox = np.array(bbox, dtype='f8')
        self.ext = ext

    @property
    def gcv_value(self):
        """ Generalized Cross Validation value """
        return self.wk[0]

    @property
    def msr(self):
        """ Mean Squared Residual """
        return self.wk[1]

    @property
    def dof(self):
        """ Estimated DOF """
        return self.wk[2]

    @property
    def p(self):
        """ Smoothing Parameter p """
        return self.wk[3]

    @property
    def mse(self):
        """ Estimated MSE. Different formula for MD == 3 """
        return self.wk[4]

    @property
    def error_variance(self):
        """ Gauss-Markov error variance """
        return self.wk[5]

    @property
    def variance(self):
        return self.wk[6:] / self.wk[5]

    def __call__(self, x, nu=0, ext=None):
        x = np.array(x)
        oldshape = x.shape
        x = np.ravel(x)
        if ext is None:
            ext = self.ext
        if ext != 0:
            bmask1 = (x < self.bbox[0])
            bmask2 = (x > self.bbox[1])
        if ext == 3:
            x = np.copy(x)
            x[bmask1] = self.bbox[0]
            x[bmask2] = self.bbox[1]
        elif ext == 2:
            if bmask1.any() or bmask2.any():
                raise ValueError("some values are out of bound")

        y = gcvspl.splderv(nu, self.splorder, x, self.x, self.c, self.L)

        if ext == 1:
            y[bmask1 | bmask2] = 0

        return y.reshape(oldshape)

class SmoothedNSpline(GCVSplineBase):
    def __init__(self, x, y, p, w=1, w1=1, nc=None, kind='cubic', bbox=None, ext=0):
        """ Natural Spline fitting with fixed smoothing

            Parameters
            ----------
            %s
            p : float
                smoothing parameter
        """
        GCVSplineBase.__init__(self, x=x, y=y, w=w, w1=w1, nc=nc, kind=kind, bbox=bbox, ext=ext, VAL=p, MD=1)

class GCVSmoothedNSpline(GCVSplineBase):
    def __init__(self, x, y, w=1, w1=1, nc=None, kind='cubic', bbox=None, ext=0):
        """ Natural Spline fitting with smoothing selected by cross validation.

            Parameters
            ----------
            %s
        """
        GCVSplineBase.__init__(self, x=x, y=y, w=w, w1=w1, nc=nc, kind=kind, bbox=bbox, ext=ext, VAL=0, MD=2)

class MSESmoothedNSpline(GCVSplineBase):
    def __init__(self, x, y, variance, w=1, w1=1, nc=None, kind='cubic', bbox=None, ext=0):
        """ Natural Spline fitting with known error variance.

            Parameters
            ----------
            %s
            variance : float
                smoothing terms of expected variance

        """
        GCVSplineBase.__init__(self, x=x, y=y, w=w, w1=w1, nc=nc, kind=kind, bbox=bbox, ext=ext, VAL=variance, MD=3)

class DOFSmoothedNSpline(GCVSplineBase):
    def __init__(self, x, y, dof, w=1, w1=1, nc=None, kind='cubic', bbox=None, ext=0):
        """ Nautural Spline fitting with known degrees of freedom.

            Parameters
            ----------
            %s
            dof : float
                desired effective degrees of freedom.
        """
        GCVSplineBase.__init__(self, x=x, y=y, w=w, w1=w1, nc=nc, kind=kind, bbox=bbox, ext=ext, VAL=dof, MD=3)

def adddocstring(klass):
    klass.__init__.__doc__ = klass.__init__.__doc__ % PARAMETERS.strip()

adddocstring(GCVSplineBase)
adddocstring(SmoothedNSpline)
adddocstring(GCVSmoothedNSpline)
adddocstring(MSESmoothedNSpline)
adddocstring(DOFSmoothedNSpline)

