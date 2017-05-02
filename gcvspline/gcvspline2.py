# -*- coding: utf-8 -*-
"""
scipy.interpolate style wrapper for gcvspl.

"""

__author__ = "Yu Feng"
__email__ = 'rainwoodman@gmail.com'

from . import _gcvspl as gcvspl
import numpy as np

class GCVSpline(object):
    SPLORDER = {'linear' : 1, 'cubic' : 2, 'quintic' : 3, 'heptic': 4, 1 : 1, 3: 2, 5 :3, 7:4}
    EXT = {'extrapolate': 0, 'zeros' : 0, 'raise' : 2, 'const' : 3}
    def __init__(self, x, y, w=1, w1=1, mode=2, prior=1, nc=None, kind='cubic', bbox=None, ext=0):
        """ Parameters
            ----------
            w: array_like
                inverted variance of observations, e.g. 1.0 / var(y, axis=1)
            w1: array_like
                inverted variance of datasets, e.g. 1.0 / var(y, axis=0)
            kind : int, or str
                kind / order of spline, NOTE: different from the internal gcvspline parameter.
            nc  : int
                number of control points.
            mode : int
                2 for cross-validation. not well understood
            prior : float
                parameter controls the behavior of the prior. not well understood.

        """

        splmode = mode # FIXME: after understanding what these modes are
                       # give them proper string names.
        if not (x[:-1] < x[1:]).all():
            raise ValueError("x must be strictly increasing")

        if len(y.shape) > 2:
            raise ValueError("y must be at most two dimensional (len(x), ndataset)")

        if len(y.shape) == 1:
            y = y.reshape(-1, 1)

        if bbox is None:
            bbox = (x[0], x[-1])

        wx = np.ones(y.shape[0])
        wx[...] = w
        wy = np.ones(y.shape[1])
        wy[...] = w1

        if nc is None: nc = len(x)

        splorder = self.SPLORDER.get(kind)

        if splorder not in (1, 2, 3, 4):
            raise ValueError("splorder must be 1, 2, 3, 4.")

        if splmode not in (1, 2, 3, 4):
            raise ValueError("splmode must be 1, 2, 3, 4")

        ext = self.EXT.get(ext, ext)
        if ext not in (0, 1, 2, 3):
            raise ValueError("exts must be 1 2, 3, 4")

        c, wk, ier = gcvspl.gcvspl(x, y, wx, wy, splorder, splmode, prior, nc)
        if ier != 0:
            raise RuntimeError(" gcvspl returned error code %d" % ierr)
        self.x = x
        self.c = c
        self.wk = wk
        self.splorder = splorder
        self.L = 1
        self.bbox = np.array(bbox, dtype='f8')
        self.ext = ext

    def __call__(self, x, nu=0, ext=None):
        x = np.array(x)
        oldshape = x.shape
        x = np.ravel(x)
        if self.ext != 0:
            bmask1 = (x < self.bbox[0])
            bmask2 = (x > self.bbox[1])
        if self.ext == 3:
            x = np.copy(x)
            x[bmask1] = self.bbox[0]
            x[bmask2] = self.bbox[1]
        elif self.ext == 2:
            if bmask1.any() or bmask2.any():
                raise ValueError("some values are out of bound")

        y = gcvspl.splderv(nu, self.splorder, x, self.x, self.c, self.L)

        if self.ext == 1:
            y[bmask1 | bmask2] = 0

        return y.reshape(oldshape)

