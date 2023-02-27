# -*- coding: utf-8 -*-
import numpy as np
from objects.CSCG._2d.mesh.domain.regions.region.interpolations.base import InterpolationBase


class Crazy(InterpolationBase):
    """ The Transfinite interpolation in 2D."""
    def __init__(self, region):
        super().__init__(region)
        assert region._domain_input_.domain_name in ('Crazy', 'CrazyPeriodic')
        self._c_ = region._domain_input_.c
        self._bounds_ = region._domain_input_.bounds
        self._freeze_self_()

    @property
    def c(self):
        return self._c_

    @property
    def bounds(self):
        return self._bounds_

    def mapping(self, r, s):
        """ r, s, t be in [0, 1]. """
        r, s = self.___check_rs___(r, s)
        a, b = self.bounds[0]
        c, d = self.bounds[1]
        x = (b - a) * (r + 0.5 * self.c * np.sin(2 * np.pi * r) * np.sin(2 * np.pi * s)) + a
        y = (d - c) * (s + 0.5 * self.c * np.sin(2 * np.pi * r) * np.sin(2 * np.pi * s)) + c
        return x, y

    def mapping_X(self, r, s):
        r, s = self.___check_rs___(r, s)
        a, b = self.bounds[0]
        x = (b - a) * (r + 0.5 * self.c * np.sin(2 * np.pi * r) * np.sin(2 * np.pi * s)) + a
        return x

    def mapping_Y(self, r, s):
        r, s = self.___check_rs___(r, s)
        c, d = self.bounds[1]
        y = (d - c) * (s + 0.5 * self.c * np.sin(2 * np.pi * r) * np.sin(2 * np.pi * s)) + c
        return y


    def Jacobian_matrix(self, r, s):
        """ r, s, t be in [0, 1]. """
        r, s = self.___check_rs___(r, s)
        a, b = self.bounds[0]
        c, d = self.bounds[1]
        xr = (b - a) + (b - a) * 2 * np.pi * 0.5 * self.c * np.cos(2 * np.pi * r) * np.sin(2 * np.pi * s)
        xs = (b - a) * 2 * np.pi * 0.5 * self.c * np.sin(2 * np.pi * r) * np.cos(2 * np.pi * s)
        yr = (d - c) * 2 * np.pi * 0.5 * self.c * np.cos(2 * np.pi * r) * np.sin(2 * np.pi * s)
        ys = (d - c) + (d - c) * 2 * np.pi * 0.5 * self.c * np.sin(2 * np.pi * r) * np.cos(2 * np.pi * s)
        return ((xr, xs),
                (yr, ys))

    def Jacobian_Xr(self, r, s):
        r, s = self.___check_rs___(r, s)
        a, b = self.bounds[0]
        xr = (b - a) + (b - a) * 2 * np.pi * 0.5 * self.c * np.cos(2 * np.pi * r) * np.sin(2 * np.pi * s)
        return xr

    def Jacobian_Xs(self, r, s):
        r, s = self.___check_rs___(r, s)
        a, b = self.bounds[0]
        xs = (b - a) * 2 * np.pi * 0.5 * self.c * np.sin(2 * np.pi * r) * np.cos(2 * np.pi * s)
        return xs

    def Jacobian_Yr(self, r, s):
        r, s = self.___check_rs___(r, s)
        c, d = self.bounds[1]
        yr = (d - c) * 2 * np.pi * 0.5 * self.c * np.cos(2 * np.pi * r) * np.sin(2 * np.pi * s)
        return yr

    def Jacobian_Ys(self, r, s):
        r, s = self.___check_rs___(r, s)
        c, d = self.bounds[1]
        ys = (d - c) + (d - c) * 2 * np.pi * 0.5 * self.c * np.sin(2 * np.pi * r) * np.cos(2 * np.pi * s)
        return ys