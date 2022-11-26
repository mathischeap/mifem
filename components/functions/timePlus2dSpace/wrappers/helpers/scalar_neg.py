# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/13/2022 7:13 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from components.freeze.main import FrozenOnly


class t2d_ScalarNeg(FrozenOnly):
    """"""

    def __init__(self, s):
        """"""
        self._s_ = s
        self._freeze_self_()

    def __call__(self, t, x, y):
        return - self._s_(t, x, y)

if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass