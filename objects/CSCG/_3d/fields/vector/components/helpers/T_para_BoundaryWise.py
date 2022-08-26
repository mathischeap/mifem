# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/26/2022 12:18 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class _3dCSCG_T_para_BW(FrozenOnly):
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

        w1 = w1(t, x, y, z)
        w2 = w2(t, x, y, z)
        w3 = w3(t, x, y, z)

        wXn2 = w3 * n1 - w1 * n3
        wXn3 = w1 * n2 - w2 * n1

        return n2 * wXn3 - n3 * wXn2

    def fy(self, t, x, y, z):
        w1, w2, w3 = self.w
        n1, n2, n3 = self.n

        w1 = w1(t, x, y, z)
        w2 = w2(t, x, y, z)
        w3 = w3(t, x, y, z)

        wXn1 = w2 * n3 - w3 * n2
        wXn3 = w1 * n2 - w2 * n1

        return n3 * wXn1 - n1 * wXn3

    def fz(self, t, x, y, z):
        w1, w2, w3 = self.w
        n1, n2, n3 = self.n

        w1 = w1(t, x, y, z)
        w2 = w2(t, x, y, z)
        w3 = w3(t, x, y, z)

        wXn1 = w2 * n3 - w3 * n2
        wXn2 = w3 * n1 - w1 * n3

        return n1 * wXn2 - n2 * wXn1





if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
