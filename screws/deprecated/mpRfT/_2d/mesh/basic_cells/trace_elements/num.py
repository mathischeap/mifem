# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/31 4:43 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class mpRfT2_Mesh_BasicCells_Num(FrozenOnly):
    """"""

    def __init__(self, elements):
        """"""
        self._elements_ = elements
        self._freeze_self_()


    @property
    def local_segments(self):
        """{int} : the amount of local segments on the trace-elements."""
        raise NotImplementedError()

if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
