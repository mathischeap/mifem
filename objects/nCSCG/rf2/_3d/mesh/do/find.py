# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly

class _3nCSCG_MeshDoFind(FrozenOnly):
    """"""
    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._freeze_self_()

    @staticmethod
    def local_indices_of_sub_cell(i):
        """"""
        if i == 0:
            return 0, 0, 0
        elif i == 1:
            return 1, 0, 0
        elif i == 2:
            return 0, 1, 0
        elif i == 3:
            return 1, 1, 0
        elif i == 4:
            return 0, 0, 1
        elif i == 5:
            return 1, 0, 1
        elif i == 6:
            return 0, 1, 1
        elif i == 7:
            return 1, 1, 1
        else:
            raise Exception()

    def origin_and_delta(self, indices):
        """consider the level-0 cell (cscg mesh-element) as a [-1,1]^3 reference domain.
        """
        delta = 2 * 0.5 ** (len(indices) - 1)
        if len(indices) == 1:  # a level-0 cell, namely, a mesh-element of the cscg mesh.
            origin = (-1, -1, -1)
        else:
            o = self.origin_and_delta(indices[:-1])[0]
            i, j, k = self.local_indices_of_sub_cell(indices[-1])
            origin = o[0] + i * delta, o[1] + j * delta, o[2] + k * delta
        return origin, delta




if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
