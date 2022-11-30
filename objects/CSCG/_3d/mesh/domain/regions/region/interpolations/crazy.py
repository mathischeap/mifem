# -*- coding: utf-8 -*-
import numpy as np
from objects.CSCG._3d.mesh.domain.regions.region.interpolations.base import InterpolationBase

class Crazy(InterpolationBase):
    """ The Crazy interpolation in 3D. """

    def __init__(self, region):
        """
        To initialize a crazy interpolation, we take a regions as input.

        Parameters
        ----------
        region : Region

        """
        super().__init__(region)
        self._c_ = region._domain_input_.c
        self._bounds_ = region._domain_input_.bounds
        self._freeze_self_()

    @property
    def c(self):
        return self._c_

    @property
    def bounds(self):
        return self._bounds_

    def mapping(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)

        a, b = self.bounds[0]
        c, d = self.bounds[1]
        e, f = self.bounds[2]
        if self.c == 0:
            x = (b - a) * r + a
            y = (d - c) * s + c
            z = (f - e) * t + e
        else:
            x = (b - a) * (r + 0.5 * self.c * np.sin(2 * np.pi * r) * np.sin(2 * np.pi * s) * np.sin(2 * np.pi * t)) + a
            y = (d - c) * (s + 0.5 * self.c * np.sin(2 * np.pi * r) * np.sin(2 * np.pi * s) * np.sin(2 * np.pi * t)) + c
            z = (f - e) * (t + 0.5 * self.c * np.sin(2 * np.pi * r) * np.sin(2 * np.pi * s) * np.sin(2 * np.pi * t)) + e

        return x, y, z

    def mapping_X(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)
        a, b = self.bounds[0]
        if self.c == 0:
            x = (b - a) * r + a
        else:
            x = (b - a) * (r + 0.5 * self.c * np.sin(2 * np.pi * r) * np.sin(2 * np.pi * s) * np.sin(2 * np.pi * t)) + a
        return x

    def mapping_Y(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)
        c, d = self.bounds[1]
        if self.c == 0:
            y = (d - c) * s + c
        else:
            y = (d - c) * (s + 0.5 * self.c * np.sin(2 * np.pi * r) * np.sin(2 * np.pi * s) * np.sin(2 * np.pi * t)) + c
        return y

    def mapping_Z(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)
        e, f = self.bounds[2]
        if self.c == 0:
            z = (f - e) * t + e
        else:
            z = (f - e) * (t + 0.5 * self.c * np.sin(2 * np.pi * r) * np.sin(2 * np.pi * s) * np.sin(2 * np.pi * t)) + e
        return z

    def Jacobian_matrix(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)

        a, b = self.bounds[0]
        c, d = self.bounds[1]
        e, f = self.bounds[2]

        if self.c == 0:
            xr = (b - a) * np.ones_like(r)  # have to do this to make it an array.
            xs = 0  # np.zeros_like(r)
            xt = 0

            yr = 0
            ys = (d - c) * np.ones_like(r)
            yt = 0

            zr = 0
            zs = 0
            zt = (f - e) * np.ones_like(r)
        else:
            xr = (b - a) + (b - a) * 2 * np.pi * 0.5 * self.c * np.cos(2 * np.pi * r) * np.sin(2 * np.pi * s) * np.sin(
                2 * np.pi * t)
            xs = (b - a) * 2 * np.pi * 0.5 * self.c * np.sin(2 * np.pi * r) * np.cos(2 * np.pi * s) * np.sin(
                2 * np.pi * t)
            xt = (b - a) * 2 * np.pi * 0.5 * self.c * np.sin(2 * np.pi * r) * np.sin(2 * np.pi * s) * np.cos(
                2 * np.pi * t)

            yr = (d - c) * 2 * np.pi * 0.5 * self.c * np.cos(2 * np.pi * r) * np.sin(2 * np.pi * s) * np.sin(
                2 * np.pi * t)
            ys = (d - c) + (d - c) * 2 * np.pi * 0.5 * self.c * np.sin(2 * np.pi * r) * np.cos(2 * np.pi * s) * np.sin(
                2 * np.pi * t)
            yt = (d - c) * 2 * np.pi * 0.5 * self.c * np.sin(2 * np.pi * r) * np.sin(2 * np.pi * s) * np.cos(
                2 * np.pi * t)

            zr = (f - e) * 2 * np.pi * 0.5 * self.c * np.cos(2 * np.pi * r) * np.sin(2 * np.pi * s) * np.sin(
                2 * np.pi * t)
            zs = (f - e) * 2 * np.pi * 0.5 * self.c * np.sin(2 * np.pi * r) * np.cos(2 * np.pi * s) * np.sin(
                2 * np.pi * t)
            zt = (f - e) + (f - e) * 2 * np.pi * 0.5 * self.c * np.sin(2 * np.pi * r) * np.sin(2 * np.pi * s) * np.cos(
                2 * np.pi * t)

        return [(xr, xs, xt),
                (yr, ys, yt),
                (zr, zs, zt)]

    def Jacobian_Xr(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)
        a, b = self.bounds[0]
        if self.c == 0:
            xr = (b - a) * np.ones_like(r)  # have to do this to make it an array.
        else:
            xr = (b - a) + (b - a) * 2 * np.pi * 0.5 * self.c * np.cos(2 * np.pi * r) * np.sin(2 * np.pi * s) * np.sin(
                2 * np.pi * t)
        return xr

    def Jacobian_Xs(self, r, s, t):
        r, s, t = self.___check_rst___(r, s, t)
        a, b = self.bounds[0]
        if self.c == 0:
            xs = np.zeros_like(r)
        else:
            xs = (b - a) * 2 * np.pi * 0.5 * self.c * np.sin(2 * np.pi * r) * np.cos(2 * np.pi * s) * np.sin(
                2 * np.pi * t)
        return xs

    def Jacobian_Xt(self, r, s, t):
        r, s, t = self.___check_rst___(r, s, t)
        a, b = self.bounds[0]
        if self.c == 0:
            xt = np.zeros_like(r)
        else:
            xt = (b - a) * 2 * np.pi * 0.5 * self.c * np.sin(2 * np.pi * r) * np.sin(2 * np.pi * s) * np.cos(
                2 * np.pi * t)
        return xt

    def Jacobian_Yr(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)
        c, d = self.bounds[1]
        if self.c == 0:
            yr = np.zeros_like(r)
        else:
            yr = (d - c) * 2 * np.pi * 0.5 * self.c * np.cos(2 * np.pi * r) * np.sin(2 * np.pi * s) * np.sin(
                2 * np.pi * t)
        return yr

    def Jacobian_Ys(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)
        c, d = self.bounds[1]
        if self.c == 0:
            ys = (d - c) * np.ones_like(r)
        else:
            ys = (d - c) + (d - c) * 2 * np.pi * 0.5 * self.c * np.sin(2 * np.pi * r) * np.cos(2 * np.pi * s) * np.sin(
                2 * np.pi * t)
        return ys

    def Jacobian_Yt(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)
        c, d = self.bounds[1]
        if self.c == 0:
            yt = np.zeros_like(r)
        else:
            yt = (d - c) * 2 * np.pi * 0.5 * self.c * np.sin(2 * np.pi * r) * np.sin(2 * np.pi * s) * np.cos(
                2 * np.pi * t)
        return yt

    def Jacobian_Zr(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)
        e, f = self.bounds[2]
        if self.c == 0:
            zr = np.zeros_like(r)
        else:
            zr = (f - e) * 2 * np.pi * 0.5 * self.c * np.cos(2 * np.pi * r) * np.sin(2 * np.pi * s) * np.sin(
                2 * np.pi * t)
        return zr

    def Jacobian_Zs(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)
        e, f = self.bounds[2]
        if self.c == 0:
            zs = np.zeros_like(r)
        else:
            zs = (f - e) * 2 * np.pi * 0.5 * self.c * np.sin(2 * np.pi * r) * np.cos(2 * np.pi * s) * np.sin(
                2 * np.pi * t)
        return zs

    def Jacobian_Zt(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)
        e, f = self.bounds[2]
        if self.c == 0:
            zt = (f - e) * np.ones_like(r)
        else:
            zt = (f - e) + (f - e) * 2 * np.pi * 0.5 * self.c * np.sin(2 * np.pi * r) * np.sin(2 * np.pi * s) * np.cos(
                2 * np.pi * t)
        return zt

    def Jacobian_X_(self, r, s, t):
        return self.Jacobian_Xr(r, s, t), self.Jacobian_Xs(r, s, t), self.Jacobian_Xt(r, s, t)

    def Jacobian_Y_(self, r, s, t):
        return self.Jacobian_Yr(r, s, t), self.Jacobian_Ys(r, s, t), self.Jacobian_Yt(r, s, t)

    def Jacobian_Z_(self, r, s, t):
        return self.Jacobian_Zr(r, s, t), self.Jacobian_Zs(r, s, t), self.Jacobian_Zt(r, s, t)
