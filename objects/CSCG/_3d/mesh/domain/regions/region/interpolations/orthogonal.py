
import numpy as np
from objects.CSCG._3d.mesh.domain.regions.region.interpolations.base import InterpolationBase



class Orthogonal(InterpolationBase):
    """ The orthogonal interpolation in 3D. """

    def __init__(self, region):
        """
        To initialize an orthogonal interpolation, we take a region as input.

        Parameters
        ----------
        region : Region

        """
        super().__init__(region)
        cc = region.corner_coordinates
        self._bounds_ = ([cc[0][0], cc[1][0]],
                         [cc[0][1], cc[2][1]],
                         [cc[0][2], cc[4][2]])
        self._freeze_self_()

    @property
    def bounds(self):
        return self._bounds_

    def mapping(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)

        a, b = self.bounds[0]
        c, d = self.bounds[1]
        e, f = self.bounds[2]
        x = (b - a) * r + a
        y = (d - c) * s + c
        z = (f - e) * t + e
        return x, y, z

    def mapping_X(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)
        a, b = self.bounds[0]
        x = (b - a) * r + a
        return x

    def mapping_Y(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)
        c, d = self.bounds[1]
        y = (d - c) * s + c
        return y

    def mapping_Z(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)
        e, f = self.bounds[2]
        z = (f - e) * t + e
        return z

    def Jacobian_matrix(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)

        a, b = self.bounds[0]
        c, d = self.bounds[1]
        e, f = self.bounds[2]

        xr = (b - a) * np.ones_like(r)  # have to do this to make it an array.
        xs = 0  # np.zeros_like(r)
        xt = 0

        yr = 0
        ys = (d - c) * np.ones_like(r)
        yt = 0

        zr = 0
        zs = 0
        zt = (f - e) * np.ones_like(r)

        return [(xr, xs, xt),
                (yr, ys, yt),
                (zr, zs, zt)]

    def Jacobian_Xr(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)
        a, b = self.bounds[0]
        xr = (b - a) * np.ones_like(r)  # have to do this to make it an array.
        return xr

    def Jacobian_Xs(self, r, s, t):
        r, s, t = self.___check_rst___(r, s, t)
        xs = np.zeros_like(r)
        return xs

    def Jacobian_Xt(self, r, s, t):
        r, s, t = self.___check_rst___(r, s, t)
        xt = np.zeros_like(r)
        return xt

    def Jacobian_Yr(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)
        yr = np.zeros_like(r)
        return yr

    def Jacobian_Ys(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)
        c, d = self.bounds[1]
        ys = (d - c) * np.ones_like(r)
        return ys

    def Jacobian_Yt(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)
        yt = np.zeros_like(r)
        return yt

    def Jacobian_Zr(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)
        zr = np.zeros_like(r)
        return zr

    def Jacobian_Zs(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)
        zs = np.zeros_like(r)
        return zs

    def Jacobian_Zt(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)
        e, f = self.bounds[2]
        zt = (f - e) * np.ones_like(r)
        return zt

    def Jacobian_X_(self, r, s, t):
        return self.Jacobian_Xr(r, s, t), self.Jacobian_Xs(r, s, t), self.Jacobian_Xt(r, s, t)

    def Jacobian_Y_(self, r, s, t):
        return self.Jacobian_Yr(r, s, t), self.Jacobian_Ys(r, s, t), self.Jacobian_Yt(r, s, t)

    def Jacobian_Z_(self, r, s, t):
        return self.Jacobian_Zr(r, s, t), self.Jacobian_Zs(r, s, t), self.Jacobian_Zt(r, s, t)