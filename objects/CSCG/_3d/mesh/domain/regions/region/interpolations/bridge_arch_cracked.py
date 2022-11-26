



from objects.CSCG._3d.mesh.domain.regions.region.interpolations.base import InterpolationBase

from components.functions._2dSpace.geometrical_functions.straight_line import StraightLine
from components.functions._2dSpace.geometrical_functions.clockwise_arc import ArcClockWise

from components.functions._2dSpace.transfinite import TransfiniteMapping

import numpy as np


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


    def mapping_X(self, r, s, t):
        """"""
        r, s, t = self.___check_rst___(r, s, t)
        if self._region_.name == 'R:R_left_up':
            x = self._TFM2D_L_.x(self._rdr_ * r, s)
        elif self._region_.name == 'R:R_left_down':
            x = self._TFM2D_L_.x(self._rdr_ + self._dr_ * r, s)
        elif self._region_.name == 'R:R_right_up':
            x = self._TFM2D_R_.x(self._rdr_ * r, s)
        elif self._region_.name == 'R:R_right_down':
            x = self._TFM2D_R_.x(self._rdr_ + self._dr_ * r, s)
        else:
            raise Exception()
        return x


    def mapping_Y(self, r, s, t):
        """"""
        r, s, t = self.___check_rst___(r, s, t)
        if self._region_.name == 'R:R_left_up':
            y = self._TFM2D_L_.y(self._rdr_ * r, s)
        elif self._region_.name == 'R:R_left_down':
            y = self._TFM2D_L_.y(self._rdr_ + self._dr_ * r, s)
        elif self._region_.name == 'R:R_right_up':
            y = self._TFM2D_R_.y(self._rdr_ * r, s)
        elif self._region_.name == 'R:R_right_down':
            y = self._TFM2D_R_.y(self._rdr_ + self._dr_ * r, s)
        else:
            raise Exception()
        return y


    def mapping_Z(self, r, s, t):
        """"""
        r, s, t = self.___check_rst___(r, s, t)
        z = self._w_ * t
        return z


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

    def Jacobian_Xr(self, r, s, t):
        r, s, t = self.___check_rst___(r, s, t)
        if self._region_.name == 'R:R_left_up':
            xr = self._rdr_ * self._TFM2D_L_.dx_dr(self._rdr_ * r, s)
        elif self._region_.name == 'R:R_left_down':
            xr = self._dr_ * self._TFM2D_L_.dx_dr(self._rdr_ + self._dr_ * r, s)
        elif self._region_.name == 'R:R_right_up':
            xr = self._rdr_ * self._TFM2D_R_.dx_dr(self._rdr_ * r, s)
        elif self._region_.name == 'R:R_right_down':
            xr = self._dr_ * self._TFM2D_R_.dx_dr(self._rdr_ + self._dr_ * r, s)
        else:
            raise Exception()
        return xr

    def Jacobian_Xs(self, r, s, t):
        r, s, t = self.___check_rst___(r, s, t)
        if self._region_.name == 'R:R_left_up':
            xs = self._TFM2D_L_.dx_ds(self._rdr_ * r, s)
        elif self._region_.name == 'R:R_left_down':
            xs = self._TFM2D_L_.dx_ds(self._rdr_ + self._dr_ * r, s)
        elif self._region_.name == 'R:R_right_up':
            xs = self._TFM2D_R_.dx_ds(self._rdr_ * r, s)
        elif self._region_.name == 'R:R_right_down':
            xs = self._TFM2D_R_.dx_ds(self._rdr_ + self._dr_ * r, s)
        else:
            raise Exception()
        return xs

    def Jacobian_Xt(self, r, s, t):
        r, s, t = self.___check_rst___(r, s, t)
        zeros_shape = np.shape(r)
        if self._region_.name == 'R:R_left_up':
            xt = np.zeros(zeros_shape)
        elif self._region_.name == 'R:R_left_down':
            xt = np.zeros(zeros_shape)
        elif self._region_.name == 'R:R_right_up':
            xt = np.zeros(zeros_shape)
        elif self._region_.name == 'R:R_right_down':
            xt = np.zeros(zeros_shape)
        else:
            raise Exception()
        return xt

    def Jacobian_Yr(self, r, s, t):
        r, s, t = self.___check_rst___(r, s, t)
        if self._region_.name == 'R:R_left_up':
            yr = self._rdr_ * self._TFM2D_L_.dy_dr(self._rdr_ * r, s)
        elif self._region_.name == 'R:R_left_down':
            yr = self._dr_ * self._TFM2D_L_.dy_dr(self._rdr_ + self._dr_ * r, s)
        elif self._region_.name == 'R:R_right_up':
            yr = self._rdr_ * self._TFM2D_R_.dy_dr(self._rdr_ * r, s)
        elif self._region_.name == 'R:R_right_down':
            yr = self._dr_ * self._TFM2D_R_.dy_dr(self._rdr_ + self._dr_ * r, s)
        else:
            raise Exception()
        return yr
    def Jacobian_Ys(self, r, s, t):
        r, s, t = self.___check_rst___(r, s, t)
        if self._region_.name == 'R:R_left_up':
            ys = self._TFM2D_L_.dy_ds(self._rdr_ * r, s)
        elif self._region_.name == 'R:R_left_down':
            ys = self._TFM2D_L_.dy_ds(self._rdr_ + self._dr_ * r, s)
        elif self._region_.name == 'R:R_right_up':
            ys = self._TFM2D_R_.dy_ds(self._rdr_ * r, s)
        elif self._region_.name == 'R:R_right_down':
            ys = self._TFM2D_R_.dy_ds(self._rdr_ + self._dr_ * r, s)
        else:
            raise Exception()
        return ys
    def Jacobian_Yt(self, r, s, t):
        r, s, t = self.___check_rst___(r, s, t)
        zeros_shape = np.shape(r)
        if self._region_.name == 'R:R_left_up':
            yt = np.zeros(zeros_shape)
        elif self._region_.name == 'R:R_left_down':
            yt = np.zeros(zeros_shape)
        elif self._region_.name == 'R:R_right_up':
            yt = np.zeros(zeros_shape)
        elif self._region_.name == 'R:R_right_down':
            yt = np.zeros(zeros_shape)
        else:
            raise Exception()
        return yt

    def Jacobian_Zr(self, r, s, t):
        r, s, t = self.___check_rst___(r, s, t)
        zeros_shape = np.shape(r)
        zr = np.zeros(zeros_shape)
        return zr
    def Jacobian_Zs(self, r, s, t):
        r, s, t = self.___check_rst___(r, s, t)
        zeros_shape = np.shape(r)
        zs = np.zeros(zeros_shape)
        return zs
    def Jacobian_Zt(self, r, s, t):
        r, s, t = self.___check_rst___(r, s, t)
        zeros_shape = np.shape(r)
        zt = np.ones(zeros_shape) * self._w_
        return zt

    def Jacobian_X_(self, r, s, t):
        return self.Jacobian_Xr(r, s, t), self.Jacobian_Xs(r, s, t), self.Jacobian_Xt(r, s, t)

    def Jacobian_Y_(self, r, s, t):
        return self.Jacobian_Yr(r, s, t), self.Jacobian_Ys(r, s, t), self.Jacobian_Yt(r, s, t)

    def Jacobian_Z_(self, r, s, t):
        return self.Jacobian_Zr(r, s, t), self.Jacobian_Zs(r, s, t), self.Jacobian_Zt(r, s, t)