# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/13/2022 5:49 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from components.freeze.main import FrozenOnly


class t3d_ScalarMultiply(FrozenOnly):
    """"""

    def __init__(self, v0, v1):
        """"""
        self._v0_ = v0
        self._v1_ = v1
        self._freeze_self_()

    def __call__(self, t, x, y, z):
        return self._v0_(t, x, y, z) * self._v1_(t, x, y, z)



if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
