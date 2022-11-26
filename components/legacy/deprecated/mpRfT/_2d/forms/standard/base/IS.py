# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/19 6:11 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly


class mpRfT2_SF_IS(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._freeze_self_()

    @property
    def hybrid(self):
        return self._f_._hybrid_

    @property
    def inner(self):
        return True if self._f_.orientation == 'inner' else False



if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
