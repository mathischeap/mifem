# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 6/2/2022 11:34 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly


class mpRfT2_SF_N(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._freeze_self_()

    def __getitem__(self, rp):
        return self._f_.mesh[rp].N




if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
