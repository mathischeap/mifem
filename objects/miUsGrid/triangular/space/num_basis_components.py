# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 7:22 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class miUsGrid_TriangularFunctionSpace_NumBasisComponents(FrozenOnly):
    """"""

    def __init__(self, space):
        """"""
        self._space_ = space
        self._freeze_self_()


    @property
    def miUsTriangular_S1F_Outer(self):
        p = self._space_.p
        return p * (p+1), p ** 2

    @property
    def miUsTriangular_S1F_Inner(self):
        p = self._space_.p
        return p ** 2, p * (p+1)

if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
