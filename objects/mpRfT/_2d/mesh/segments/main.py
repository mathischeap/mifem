# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 5/22/2022 10:04 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.mpRfT._2d.mesh.segments.bcW import mpRfT2_Mesh_Segments_bcW
from objects.mpRfT._2d.mesh.segments.visualize import mpRfT2_Mesh_Segments_Visualize


class mpRfT2_Mesh_Segments(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._bcW_ = None
        self._visualize_ = None
        self._freeze_self_()

    def __iter__(self):
        """Go through all local segments

        Returns
        -------

        """
        #---- first go through all local segment on lv0-trace-elements
        for i in self._mesh_.basic_cells.trace_elements:
            segments = self._mesh_.basic_cells.trace_segments[i]
            for seg in segments:
                yield seg

        #--- then go through all internal segments.
        for i in self._mesh_.basic_cells.internal_segments:
            segments = self._mesh_.basic_cells.internal_segments[i]
            for seg in segments:
                yield seg

    @property
    def bcW(self):
        """Basic-cell-wise segments."""
        if self._bcW_ is None:
            self._bcW_ = mpRfT2_Mesh_Segments_bcW(self)
        return self._bcW_

    @property
    def visualize(self):
        """Basic-cell-wise segments."""
        if self._visualize_ is None:
            self._visualize_ = mpRfT2_Mesh_Segments_Visualize(self)
        return self._visualize_



if __name__ == '__main__':
    # mpiexec -n 4 python objects/mpRfT/_2d/mesh/segments/main.py

    # from objects.mpRfT._2d.master import MeshGenerator
    # mesh = MeshGenerator('rectangle')([3,3], 2, show_info=True)

    from __init__ import rfT2
    mesh = rfT2.rm(50)
    # # for seg in mesh.segments:
    # #     print(seg)
    # mesh.segments.visualize()
    bcW = mesh.segments.bcW

    for i in bcW:
        print(i, bcW[i])
