# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/8/10 20:47
"""
import sys

if './' not in sys.path: sys.path.append('./')
from screws.freeze.main import FrozenOnly


class MDM_MUL(FrozenOnly):
    def __init__(self, MDM, other):
        """"""
        assert isinstance(other, (int, float))
        self._mdm_ = MDM
        self._other_ = other
        self._freeze_self_()

    def __call__(self, basic_unit):
        return self._other_ * self._mdm_[basic_unit]


if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
