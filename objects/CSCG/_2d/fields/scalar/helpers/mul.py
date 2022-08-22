# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/21/2022 4:37 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from screws.freeze.main import FrozenOnly


class _2dCSCG_ScaMulHelper(FrozenOnly):
    def __init__(self, func, number):
        """"""
        self._f_ = func
        self._n_ = number
        self._freeze_self_()

    def __call__(self, t, x, y):
        return self._f_(t, x, y) * self._n_

class _2dCSCG_ScaMulHelper1(FrozenOnly):
    def __init__(self, sf, vf):
        """"""
        self._sf_ = sf
        self._vf_ = vf
        self._freeze_self_()

    def __call__(self, t, x, y):
        return self._sf_(t, x, y) * self._vf_(t, x, y)

if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
