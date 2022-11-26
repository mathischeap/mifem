# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/13/2022 4:57 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from components.freeze.main import FrozenOnly


class t2d_ScalarSub(FrozenOnly):
    """"""

    def __init__(self, s0, s1):
        """"""
        self._s0_ = s0
        self._s1_ = s1
        self._freeze_self_()

    def __call__(self, t, x, y):
        return self._s0_(t, x, y) - self._s1_(t, x, y)


if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
