# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/8/8 22:55
"""
import sys

if './' not in sys.path: sys.path.append('/')
from components.freeze.main import FrozenOnly


class ColVec_PRIVATE_sum(FrozenOnly):
    def __init__(self, ColVec):
        """"""
        self._cv_ = ColVec
        self._freeze_self_()

    def __call__(self, basic_unit):
        """"""
        SUM = list()
        for v in self._cv_:
            SUM.append(v[basic_unit])
        SUM = sum(SUM)
        return SUM

    def __KG_call__(self, basic_unit):
        """"""
        KEY = list()
        for v in self._cv_:
            KEY.append(v._KG_(basic_unit))
        return ''.join(KEY)


if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
