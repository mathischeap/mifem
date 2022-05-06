# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.nCSCG.rf2._2d.mesh.cell.sub_cells.main import _2nCSCG_SubCells

class _2nCSCG_CellDo(FrozenOnly):
    """"""
    def __init__(self, cell):
        """"""
        self._cell_ = cell
        self._freeze_self_()

    def refine(self):
        """Active the subgrid of this cell."""
        assert self._cell_.___sub_cells___ is None, f"sub-grid of cell {self._cell_.indices} is already on."
        self._cell_.___sub_cells___ = _2nCSCG_SubCells(self._cell_)

    def dilute(self):
        """kill the subgrid of this cell."""
        assert self._cell_.___sub_cells___ is not None, f"sub-grid of cell {self._cell_.indices} is already off."
        self._cell_.___sub_cells___ = None


if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
