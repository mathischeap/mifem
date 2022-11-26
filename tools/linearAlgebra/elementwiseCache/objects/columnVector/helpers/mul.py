# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/12/2022 12:34 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from components.freeze.main import FrozenOnly


class ColVec_MUL(FrozenOnly):
    def __init__(self, vec, f):
        """"""
        self._vec_ = vec
        self._f_ = f
        self._freeze_self_()

    def __call__(self, basic_unit):
        """"""
        return self._f_ * self._vec_[basic_unit]


if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
