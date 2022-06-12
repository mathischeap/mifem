# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 5/17/2022 10:47 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.mpRfT._2d.mesh.cell.sub_cells import mpRfT2_Mesh_Cell_SubCells


class mpRfT2_Mesh_Cell_Do(FrozenOnly):
    """"""

    def __init__(self, cell):
        """"""
        self._cell_ = cell
        self._freeze_self_()

    def ___Pr_h_refine___(self):
        """Active the subgrid of this cell. Should only be called during initialization of the mesh."""
        assert self._cell_.___sub_cells___ is None, f"sub-cells of cell {self._cell_.indices} are already on."
        self._cell_.___sub_cells___ = mpRfT2_Mesh_Cell_SubCells(self._cell_)
        self._cell_.___isroot___ = False



if __name__ == '__main__':
    # mpiexec -n 4 python objects/mpRfT/_2d/mesh/cell/do.py
    pass
