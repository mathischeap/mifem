



from _3dCSCG.mesh.regions.region.interpolations.base import InterpolationBase

from SCREWS.functions._2d_transfinite import StraightLine
from SCREWS.functions._2d_transfinite import ArcClockWise

from SCREWS.functions._2d_transfinite import TransfiniteMapping

import numpy as np


# noinspection PyAbstractClass
class BridgeArchCracked(InterpolationBase):
    """ The BridgeArchCracked interpolation in 3D. """

    def __init__(self, region):
        """
        To initialize a BridgeArchCracked interpolation, we take a regions as input.

        Parameters
        ----------
        region : Region

        """
        super().__init__(region)
        l = region._domain_input_.l
        self._w_ = region._domain_input_.w
        h = region._domain_input_.h
        c = region._domain_input_.center
        beta = self._region_._domain_input_.beta
        d = self._region_._domain_input_.d
        self._rdr_ = (beta - d) / beta
        self._dr_ = d / beta
        gamma0, dgamma0 = StraightLine((0, 0), (h, 0))()
        gamma1, dgamma1 = ArcClockWise(c, (h, 0), (beta, l / 2))()
        gamma2, dgamma2 = StraightLine((0, l / 2), (beta, l / 2))()
        gamma3, dgamma3 = StraightLine((0, 0), (0, l / 2))()
        self._TFM2D_L_ = TransfiniteMapping((gamma0, gamma1, gamma2, gamma3),
                                            dgamma=(dgamma0, dgamma1, dgamma2, dgamma3))
        gamma0, dgamma0 = StraightLine((0, l / 2), (beta, l / 2))()
        gamma1, dgamma1 = ArcClockWise(c, (beta, l / 2), (h, l))()
        gamma2, dgamma2 = StraightLine((0, l), (h, l))()
        gamma3, dgamma3 = StraightLine((0, l / 2), (0, l))()
        self._TFM2D_R_ = TransfiniteMapping((gamma0, gamma1, gamma2, gamma3),
                                            dgamma=(dgamma0, dgamma1, dgamma2, dgamma3))
        self._freeze_self_()

    def mapping(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)
        if self._region_.name == 'R:R_left_up':
            x, y = self._TFM2D_L_.mapping(self._rdr_ * r, s)
        elif self._region_.name == 'R:R_left_down':
            x, y = self._TFM2D_L_.mapping(self._rdr_ + self._dr_ * r, s)
        elif self._region_.name == 'R:R_right_up':
            x, y = self._TFM2D_R_.mapping(self._rdr_ * r, s)
        elif self._region_.name == 'R:R_right_down':
            x, y = self._TFM2D_R_.mapping(self._rdr_ + self._dr_ * r, s)
        else:
            raise Exception()
        z = self._w_ * t
        return x, y, z

    def Jacobian_matrix(self, r, s, t):
        """ r, s, t be in [0, 1]. """
        r, s, t = self.___check_rst___(r, s, t)
        zeros_shape = np.shape(r)
        if self._region_.name == 'R:R_left_up':
            xr = self._rdr_ * self._TFM2D_L_.dx_dr(self._rdr_ * r, s)
            xs = self._TFM2D_L_.dx_ds(self._rdr_ * r, s)
            xt = np.zeros(zeros_shape)
            yr = self._rdr_ * self._TFM2D_L_.dy_dr(self._rdr_ * r, s)
            ys = self._TFM2D_L_.dy_ds(self._rdr_ * r, s)
            yt = np.zeros(zeros_shape)
        elif self._region_.name == 'R:R_left_down':
            xr = self._dr_ * self._TFM2D_L_.dx_dr(self._rdr_ + self._dr_ * r, s)
            xs = self._TFM2D_L_.dx_ds(self._rdr_ + self._dr_ * r, s)
            xt = np.zeros(zeros_shape)
            yr = self._dr_ * self._TFM2D_L_.dy_dr(self._rdr_ + self._dr_ * r, s)
            ys = self._TFM2D_L_.dy_ds(self._rdr_ + self._dr_ * r, s)
            yt = np.zeros(zeros_shape)
        elif self._region_.name == 'R:R_right_up':
            xr = self._rdr_ * self._TFM2D_R_.dx_dr(self._rdr_ * r, s)
            xs = self._TFM2D_R_.dx_ds(self._rdr_ * r, s)
            xt = np.zeros(zeros_shape)
            yr = self._rdr_ * self._TFM2D_R_.dy_dr(self._rdr_ * r, s)
            ys = self._TFM2D_R_.dy_ds(self._rdr_ * r, s)
            yt = np.zeros(zeros_shape)
        elif self._region_.name == 'R:R_right_down':
            xr = self._dr_ * self._TFM2D_R_.dx_dr(self._rdr_ + self._dr_ * r, s)
            xs = self._TFM2D_R_.dx_ds(self._rdr_ + self._dr_ * r, s)
            xt = np.zeros(zeros_shape)
            yr = self._dr_ * self._TFM2D_R_.dy_dr(self._rdr_ + self._dr_ * r, s)
            ys = self._TFM2D_R_.dy_ds(self._rdr_ + self._dr_ * r, s)
            yt = np.zeros(zeros_shape)
        else:
            raise Exception()
        zr = np.zeros(zeros_shape)
        zs = np.zeros(zeros_shape)
        zt = np.ones(zeros_shape) * self._w_
        return ((xr, xs, xt),
                (yr, ys, yt),
                (zr, zs, zt))
