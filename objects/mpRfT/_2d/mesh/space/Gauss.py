# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/13 12:01 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from screws.quadrature import Quadrature


class mpRfT2_Mesh_Space_Gauss(FrozenOnly):
    """"""

    def __init__(self, space):
        """"""
        self._space_ = space
        self._GN_ = dict()
        self._freeze_self_()

    def __call__(self, N):

        if N not in self._GN_:
            self._GN_[N] =  Quadrature(N, category='Gauss').quad
        return self._GN_[N]




if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/mesh/space/Gauss.py
    pass
