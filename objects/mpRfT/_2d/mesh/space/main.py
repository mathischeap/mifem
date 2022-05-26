# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/19 1:44 PM
"""
import sys

if './' not in sys.path: sys.path.append('../')

from screws.freeze.base import FrozenOnly
from objects.CSCG._2d.spaces.polynomials import _2dCSCG_PolynomialSpace
from objects.mpRfT._2d.mesh.space.do import mpRfT2_Mesh_Space_do


class mpRfT2_Mesh_Space(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._pool_ = dict()
        self._do_ = mpRfT2_Mesh_Space_do(self)
        self._freeze_self_()

    def __getitem__(self, N):
        if N in self._pool_:
            pass
        else:
            self._pool_[N] = _2dCSCG_PolynomialSpace(('Lobatto', N), ndim=2)
        return self._pool_[N]

    @property
    def do(self):
        return self._do_




if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
