# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/26/2022 1:57 AM
"""

from components.freeze.base import FrozenOnly


class _3dCSCG_T_perp_BW(FrozenOnly):
    """"""

    def __init__(self, mesh, bn, bn_func):
        """"""
        self.mesh = mesh
        self.bn = bn
        self.w = bn_func
        self.n = mesh.boundaries[bn].coordinate_transformation.constant.outward_unit_normal_vector
        self._freeze_self_()

    def fx(self, t, x, y, z):
        w1, w2, w3 = self.w
        n1, n2, n3 = self.n

        w2 = w2(t, x, y, z)
        w3 = w3(t, x, y, z)

        return w2 * n3 - w3 * n2

    def fy(self, t, x, y, z):
        w1, w2, w3 = self.w
        n1, n2, n3 = self.n

        w1 = w1(t, x, y, z)
        w3 = w3(t, x, y, z)

        return w3 * n1 - w1 * n3

    def fz(self, t, x, y, z):
        w1, w2, w3 = self.w
        n1, n2, n3 = self.n

        w1 = w1(t, x, y, z)
        w2 = w2(t, x, y, z)

        return w1 * n2 - w2 * n1
