# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 5/17/2022 10:49 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class mpRfT2_Mesh_Cell_IS(FrozenOnly):
    """"""

    def __init__(self, cell):
        """"""
        self._cell_ = cell
        self._freeze_self_()

    @property
    def root(self):
        return self._cell_.___isroot___




if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
