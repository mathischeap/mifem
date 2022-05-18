# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/08 5:17 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.nCSCG.rf2._2d.mesh.segments.do import _2nCSCG_SegmentsDo
from objects.nCSCG.rf2._2d.mesh.segments.BCW import _2nCSCG_SegmentsBCW
from objects.nCSCG.rf2._2d.mesh.segments.visualize import _2nCSCG_SegmentsVisualize


class _2nCSCG_Segments(FrozenOnly):
    """All the local segments; segments of all local (sub-)cells."""
    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._do_ = None
        self._BCW_ = None
        self._visualize_ = None
        self._freeze_self_()

    @property
    def do(self):
        if self._do_ is None:
            self._do_ = _2nCSCG_SegmentsDo(self)
        return self._do_

    def __iter__(self):
        """Go through all local segments regardless of the base cells and lv0-trace-elements.

        Returns
        -------

        """
        #---- first go through all local segment on lv0-trace-elements
        for i in self._mesh_.Lv0trace.elements:
            segments = self._mesh_.Lv0trace.elements.segments[i]
            for seg in segments:
                yield seg

        #--- then go through all internal segments.
        for i in self._mesh_.base_cells.internal_segments:
            segments = self._mesh_.base_cells.internal_segments[i]
            for seg in segments:
                yield seg

    @property
    def BCW(self):
        """Base-Cell-Wise."""
        if self._BCW_ is None:
            self._BCW_ = _2nCSCG_SegmentsBCW(self)
        return self._BCW_

    @property
    def visualize(self):
        if self._visualize_ is None:
            self._visualize_ = _2nCSCG_SegmentsVisualize(self)
        return self._visualize_









if __name__ == "__main__":
    # mpiexec -n 4 python objects/nCSCG/rfT2/_2d/mesh/segments/main.py
    from objects.nCSCG.rf2._2d.__tests__.Random.mesh import random_mesh_of_elements_around as rm2
    mesh = rm2(100, refinement_intensity=0.5)
    # from root.read.main import read
    # mesh = read('test_mesh.mi')



