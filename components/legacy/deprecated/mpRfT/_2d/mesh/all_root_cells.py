# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 5/23/2022 1:07 AM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly
from root.config.main import COMM


class mpRfT2_Mesh_AllRootCells(FrozenOnly):
    """root-cells across all cores."""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._freeze_self_()

    def __iter__(self):
        """"""
        local_root_cells = list()
        for i in self._mesh_:
            local_root_cells.append(self._mesh_[i].__repr__())

        local_root_cells = COMM.allgather(local_root_cells)
        ___ = list()
        for _ in local_root_cells:
            ___.extend(_)
        local_root_cells = ___
        for rp in local_root_cells:
            yield rp





if __name__ == '__main__':
    # mpiexec -n 4 python objects/mpRfT/_2d/mesh/all_root_cells.py
    pass
