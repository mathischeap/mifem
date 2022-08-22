# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/20/2022 4:38 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from screws.freeze.main import FrozenOnly


class _3dCSCG_VecMulHelper(FrozenOnly):
    def __init__(self, func, number):
        """"""
        self._func_ = func
        self._number_ = number
        self._freeze_self_()

    def __call__(self, t, x, y, z):
        return self._number_ * self._func_(t, x, y, z)


if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
