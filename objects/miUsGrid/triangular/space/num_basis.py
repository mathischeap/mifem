# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/20 2:31 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from components.freeze.base import FrozenOnly


class miUsGrid_TriangularFunctionSpace_NumBasis(FrozenOnly):
    """"""

    def __init__(self, space):
        """"""
        self._space_ = space
        self._freeze_self_()


    @property
    def miUsTriangular_S0F_Outer(self):
        p = self._space_.p
        return p * (p+1) + 1

    @property
    def miUsTriangular_S0F_Inner(self):
        p = self._space_.p
        return p * (p+1) + 1


    @property
    def miUsTriangular_S1F_Outer(self):
        p = self._space_.p
        return p * p + p * (p+1)

    @property
    def miUsTriangular_S1F_Inner(self):
        p = self._space_.p
        return p * p + p * (p+1)


    @property
    def miUsTriangular_S2F_Outer(self):
        p = self._space_.p
        return p * p

    @property
    def miUsTriangular_S2F_Inner(self):
        p = self._space_.p
        return p * p



if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
