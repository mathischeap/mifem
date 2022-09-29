# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/19 9:56 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class miUsGrid_TriangularMesh_Elements_Num(FrozenOnly):
    """"""

    def __init__(self, elements):
        """"""
        self._elements_ = elements
        self._freeze_self_()

    @property
    def GLOBAL_cells(self):
        """{int} The number of global cells (elements)."""
        return self._elements_.distributions[-1].stop

    @property
    def cells(self):
        """{int} amount of local cells (elements)."""
        return len(self._elements_.range)

    @property
    def GLOBAL_elements(self):
        """{int} The number of global cells (elements)."""
        return self.GLOBAL_cells

    @property
    def elements(self):
        """{int} amount of local cells (elements)."""
        return self.cells

    @property
    def GLOBAL_points(self):
        return self._elements_.__num_GLOBAL_points__


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
