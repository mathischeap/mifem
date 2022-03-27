

from sympy import Point3D
from sympy import Plane as sympyPlane
from screws.decorators.accepts import memoize1
from objects.CSCG._3d.mesh.domain.regions.region.side_geometries.base import SideGeometryBase

class Plane(SideGeometryBase):
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