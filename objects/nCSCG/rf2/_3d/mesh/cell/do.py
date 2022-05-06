# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.nCSCG.rf2._3d.mesh.cell.sub_cells.main import _3nCSCG_SubCells

class _3nCSCG_CellDo(FrozenOnly):
    """"""
    def __init__(self, cell):
        """"""
        self._cell_ = cell
        self._freeze_self_()

    def refine(self):
        """Active the subgrid of this cell."""
        assert self._cell_.sub_cells is None, f"sub-grid of cell {self._cell_.indices} is already on."
        self._cell_.___sub_cells___ = _3nCSCG_SubCells(self._cell_)

    def dilute(self):
        """kill the subgrid of this cell."""
        assert self._cell_.sub_cells is not None, f"sub-grid of cell {self._cell_.indices} is already off."
        self._cell_.___sub_cells___ = None





if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
