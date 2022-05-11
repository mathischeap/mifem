# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.nCSCG.rf2._2d.mesh.cell.main import _2nCSCG_RF2_MeshCell
from objects.nCSCG.rf2._2d.mesh.base_cells.internal_segments.main import BaseCellInternalSegments

class _2nCSCG_RF2_BaseMeshCells(FrozenOnly):
    """"""
    def __init__(self, mesh):
        self._mesh_ = mesh
        self._cells_ = dict()
        self._InternalSegments_ = None
        self._freeze_self_()

    def __getitem__(self, i):
        if i in self._cells_:
            pass
        elif i in self._mesh_.cscg.elements:
            self._cells_[i] = _2nCSCG_RF2_MeshCell(self._mesh_, 0, (i,))
        else:
            raise Exception()
        return self._cells_[i]

    def __iter__(self):
        """go through all local 0-level cells (cscg mesh elements.)"""
        for i in self._mesh_.cscg.elements:
            yield i

    def __contains__(self, item):
        return item in self._mesh_.cscg.elements

    @property
    def internal_segments(self):
        if self._InternalSegments_ is None:
            self._InternalSegments_ = BaseCellInternalSegments(self._mesh_)
        return self._InternalSegments_



if __name__ == '__main__':
    # mpiexec -n 4 python objects/nCSCG/rf2/_2d/mesh/base_cells.py
    pass
