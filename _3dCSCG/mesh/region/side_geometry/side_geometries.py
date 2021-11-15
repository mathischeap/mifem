# -*- coding: utf-8 -*-
"""
With SideGeometry, we store the mapping form (p, q) to the side geometry, while
(p, q)= ((0,1), (0,1)).

Yi Zhang (C)
Created on Tue Sep  4 21:30:37 2018
Aerodynamics, AE
TU Delft
"""
import numpy as np
from SCREWS.frozen import FrozenOnly
from sympy import Point3D
from sympy import Plane as sympyPlane
from SCREWS.decorators import memoize1

class SideGeometry(FrozenOnly):
    def __init__(self, corner_coordinates, side_type):
        self._corner_coordinates_ = corner_coordinates
        self._side_type_ = side_type
        assert np.shape(self.corner_coordinates) == (4, 3), \
            " <SideGeometryPlane> : corner_coordinates={} wrong.".format(self.corner_coordinates)
        self._freeze_self_()
    
    @property
    def corner_coordinates(self):
        return self._corner_coordinates_
    
    @property
    def side_type(self):
        return self._side_type_



# --- particular side geometries below ------------------------------------------------------------------

class Customized(SideGeometry):
    """ """

    def __init__(self, cc, st):
        """ """
        super().__init__(cc, st)
        assert self.side_type[0] == 'customized', \
            " <SideGeometryPlane> : side_type[0]={} wrong.".format(self.side_type)
        self._melt_self_()
        self._mapping_ = self.side_type[1]
        self._Jacobian_ = self.side_type[2]
        self._freeze_self_()

    # X __
    def X(self, p, q):
        return self._mapping_[0](p, q)

    def Xp(self, p, q):
        return self._Jacobian_[0][0](p, q)

    def Xq(self, p, q):
        return self._Jacobian_[0][1](p, q)

    # Y __
    def Y(self, p, q):
        return self._mapping_[1](p, q)

    def Yp(self, p, q):
        return self._Jacobian_[1][0](p, q)

    def Yq(self, p, q):
        return self._Jacobian_[1][1](p, q)

    # Z __
    def Z(self, p, q):
        return self._mapping_[2](p, q)

    def Zp(self, p, q):
        return self._Jacobian_[2][0](p, q)

    def Zq(self, p, q):
        return self._Jacobian_[2][1](p, q)


class Free(SideGeometry):
    """
    A free side geometry is a surface we do not really care about its
    classification. A region of such side geometry(s) usually will call a
    specific interpolator to generate the mapping and therefore Jacobian and
    so on, like the crazy mapping. The crazy mapping is a analytical mapping
    that we only need to know the bounds and c which is stored in the
    `domain_input`.

    While for the transfinite mapping, the side geometries is essential. We can
    not set it to be free.

    """

    def __init__(self, cc, st):
        """ """
        super().__init__(cc, st)
        assert self.side_type == ('free',), \
            " <SideGeometry> <Free> : side_type={} wrong.".format(self.side_type)


class Plane(SideGeometry):
    """ """

    def __init__(self, cc, st):
        """ """
        super().__init__(cc, st)
        assert self.side_type == ('plane',), \
            " <SideGeometryPlane> : side_type={} wrong.".format(self.side_type)
        self._check_if_four_corners_form_a_plane_()

    def _check_if_four_corners_form_a_plane_(self):
        data = self.corner_coordinates
        distance = sympyPlane(Point3D(data[0]), Point3D(data[1]), Point3D(data[2])).distance(Point3D(data[3]))
        assert distance <= 10 ** -12, " <SideGeometryPlane> : four corners do not fit a plane."

    @memoize1
    def _fit_abcd_(self, _data_):
        if _data_ == 'x':
            _data_ = self.corner_coordinates[:, 0]
        elif _data_ == 'y':
            _data_ = self.corner_coordinates[:, 1]
        elif _data_ == 'z':
            _data_ = self.corner_coordinates[:, 2]
        else:
            pass
        d = _data_[0]
        b = _data_[1] - d
        c = _data_[2] - d
        a = _data_[3] - b - c - d
        return a, b, c, d

    # X________________________________________________________________________
    def X(self, p, q):
        a, b, c, d = self._fit_abcd_('x')
        return a * p * q + b * p + c * q + d

    def Xp(self, p, q):
        a, b, _, _ = self._fit_abcd_('x')
        return a * q + b + 0 * p

    def Xq(self, p, q):
        a, _, c, _ = self._fit_abcd_('x')
        return a * p + c + 0 * q

    # Y________________________________________________________________________
    def Y(self, p, q):
        a, b, c, d = self._fit_abcd_('y')
        return a * p * q + b * p + c * q + d

    def Yp(self, p, q):
        a, b, _, _ = self._fit_abcd_('y')
        return a * q + b + 0 * p

    def Yq(self, p, q):
        a, _, c, _ = self._fit_abcd_('y')
        return a * p + c + 0 * q

    # Z________________________________________________________________________
    def Z(self, p, q):
        a, b, c, d = self._fit_abcd_('z')
        return a * p * q + b * p + c * q + d

    def Zp(self, p, q):
        a, b, _, _ = self._fit_abcd_('z')
        return a * q + b + 0 * p

    def Zq(self, p, q):
        a, _, c, _ = self._fit_abcd_('z')
        return a * p + c + 0 * q