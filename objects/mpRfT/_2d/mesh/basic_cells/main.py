# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/17 7:54 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.mpRfT._2d.mesh.cell.main import mpRfT2_Mesh_Cell


class mpRfT2_Mesh_BasicCells(FrozenOnly):
    """"""
    def __init__(self, mesh):
        self._mesh_ = mesh
        self._cells_ = dict()
        for i in self:
            self._cells_[i] = mpRfT2_Mesh_Cell(self._mesh_, 0, (i,))
        self._freeze_self_()

    def __getitem__(self, i):
        return self._cells_[i]

    def __iter__(self):
        """go through all local 0-level cells (cscg mesh elements.)"""
        for i in self._mesh_.cscg.elements:
            yield i

    def __contains__(self, item):
        return item in self._mesh_.cscg.elements





if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/mesh/basic_cells/main.py
    pass
