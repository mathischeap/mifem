# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/20/2022 5:51 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from screws.freeze.main import FrozenOnly


class _3dCSCG_ScaMulHelper(FrozenOnly):
    def __init__(self, func, number):
        """"""
        self._func_ = func
        self._number_ = number
        self._freeze_self_()

    def __call__(self, t, x,  y, z):
        return self._func_(t, x, y, z) * self._number_


class _3dCSCG_ScaMulHelper1(FrozenOnly):
    def __init__(self, sfc, vfc):
        self._sfc_ = sfc
        self._vfc_ = vfc
        self._freeze_self_()

    def __call__(self, t, x, y, z):
        return self._sfc_(t, x, y, z) * self._vfc_(t, x, y, z)

if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
