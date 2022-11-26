# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/8/8 23:22
"""
import sys

if './' not in sys.path: sys.path.append('./')
from components.freeze.main import FrozenOnly


class SpaMat_PRIVATE_sum(FrozenOnly):
    def __init__(self, MATs):
        """"""
        self._Ms_ = MATs
        self._freeze_self_()

    def __call__(self, basic_unit):
        """"""
        SUM = list()
        for m in self._Ms_:
            SUM.append(m[basic_unit])
        return sum(SUM)

    def __KG_call__(self, basic_unit):
        """"""
        KEY = list()
        for v in self._Ms_:
            KEY.append(v._KG_(basic_unit))
        return ''.join(KEY)



if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
