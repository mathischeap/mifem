# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/17 7:54 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly
from objects.mpRfT._2d.mesh.cell.main import mpRfT2_Mesh_Cell
from objects.mpRfT._2d.mesh.basic_cells.trace_elements.main import mpRfT2_Mesh_BasicCells_TraceElements
from objects.mpRfT._2d.mesh.basic_cells.internal_segments import mpRfT2_Mesh_BasicCells_InternalSegments


class mpRfT2_Mesh_BasicCells(FrozenOnly):
    """"""
    def __init__(self, mesh):
        self._mesh_ = mesh
        self._cells_ = dict()
        for i in self._mesh_.cscg.elements:
            self._cells_[i] = mpRfT2_Mesh_Cell(self._mesh_, 0, (i,))
        self._trace_elements_ = None
        self._internal_segments_ = None
        self._trace_segments_ = None
        self._freeze_self_()

    def __getitem__(self, i):
        return self._cells_[i]

    def __iter__(self):
        """go through all local 0-level cells (cscg mesh elements.)"""

        for i in self._mesh_.cscg.elements:
            yield i

    def __contains__(self, item):
        return item in self._mesh_.cscg.elements

    @property
    def trace_elements(self):
        if self._trace_elements_ is None:
            self._trace_elements_ = mpRfT2_Mesh_BasicCells_TraceElements(self._mesh_)
        return self._trace_elements_

    @property
    def internal_segments(self):
        if self._internal_segments_ is None:
            self._internal_segments_ = mpRfT2_Mesh_BasicCells_InternalSegments(self._mesh_)
        return self._internal_segments_

    @property
    def trace_segments(self):
        return self.trace_elements.segments




if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/mesh/basic_cells/main.py
    pass
