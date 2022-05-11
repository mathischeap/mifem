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
from objects.nCSCG.rf2._2d.mesh.cell.do.find import _2nCSCG_CellDoFind

class _2nCSCG_CellDo(FrozenOnly):
    """"""
    def __init__(self, cell):
        """"""
        self._cell_ = cell
        self._find_ = None
        self._freeze_self_()

    def refine(self):
        """Active the subgrid of this cell."""
        assert not self._cell_.mesh._locker_, f"mesh is locked, unlock it first!"
        assert self._cell_.___sub_cells___ is None, f"sub-cells of cell {self._cell_.indices} are already on."
        self._cell_.___sub_cells___ = _2nCSCG_SubCells(self._cell_)
        self._cell_.___isroot___ = False
        self._cell_.mesh.base_cells.internal_segments._2bud_Sgs_[self._cell_.indices[0]] = True


    def dilute(self):
        """kill the subgrid of this cell."""
        assert not self._cell_.mesh._locker_, f"mesh is locked, unlock it first!"
        assert self._cell_.___sub_cells___ is not None, f"sub-cells of cell {self._cell_.indices} are already off."
        self._cell_.___sub_cells___ = None
        self._cell_.___isroot___ = True
        self._cell_.mesh.base_cells.internal_segments._2bud_Sgs_[self._cell_.indices[0]] = True

    @property
    def find(self):
        if self._find_ is None:
            self._find_ = _2nCSCG_CellDoFind(self._cell_)
        return self._find_




if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
