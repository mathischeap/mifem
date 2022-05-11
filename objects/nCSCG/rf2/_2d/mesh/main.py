# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/04 2:56 PM
"""
import sys
if './' not in sys.path: sys.path.append('./')

from objects.nCSCG.rf2.base.mesh.main import nCSCG_RF2_MeshBase
from objects.nCSCG.rf2._2d.mesh.base_cells.main import _2nCSCG_RF2_BaseMeshCells
from objects.nCSCG.rf2._2d.mesh.do.main import _2nCSCG_MeshDo
from objects.nCSCG.rf2._2d.mesh.visualize.main import _2nCSCG_MeshVisualize
from objects.nCSCG.rf2._2d.mesh.Lv0_trace.main import _2nCSCG_MeshLv0Trace
from objects.nCSCG.rf2._2d.mesh.segments.main import _2nCSCG_Segments
from objects.nCSCG.rf2._2d.mesh.boundaries.main import _2nCSCG_RF2_MeshBoundaries
from objects.nCSCG.rf2._2d.mesh.IS import _2nCSCG_Mesh_RF2_IS

class _2nCSCG_RF2_Mesh(nCSCG_RF2_MeshBase):
    """"""
    def __init__(self, cscg):
        """"""
        super(_2nCSCG_RF2_Mesh, self).__init__(cscg)
        self._do_ = None
        self._visualize_ = None
        self._Lv0trace_ = None
        self._segments_ = None
        self._boundaries_ = None
        self._IS_ = None
        self._freeze_self_()
        self.___base_mesh_cells___ = _2nCSCG_RF2_BaseMeshCells(self)

    @property
    def do(self):
        if self._do_ is None:
            self._do_ = _2nCSCG_MeshDo(self)
        return self._do_

    @property
    def visualize(self):
        if self._visualize_ is None:
            self._visualize_ = _2nCSCG_MeshVisualize(self)
        return self._visualize_

    @property
    def Lv0trace(self):
        """Trace"""
        if self._Lv0trace_ is None:
            self._Lv0trace_ = _2nCSCG_MeshLv0Trace(self)
        return self._Lv0trace_

    @property
    def segments(self):
        if self._segments_ is None:
            self._segments_ = _2nCSCG_Segments(self)
        return self._segments_

    @property
    def boundaries(self):
        if self._boundaries_ is None:
            self._boundaries_ = _2nCSCG_RF2_MeshBoundaries(self)
        return self._boundaries_

    @property
    def IS(self):
        if self._IS_ is None:
            self._IS_ = _2nCSCG_Mesh_RF2_IS(self)
        return self._IS_







if __name__ == '__main__':
    # mpiexec -n 4 python objects/nCSCG/rf2/_2d/mesh/main.py
    # from objects.nCSCG.rf2._2d.__tests__.Random.mesh import random_mesh_of_elements_around as rm2
    # mesh = rm2(100, refinement_intensity=0.5)

    from objects.nCSCG.rf2._2d.master import MeshGenerator

    mesh = MeshGenerator('rectangle', region_layout=[2,2])([3, 3], EDM=None, show_info=True)

    mesh.do.unlock()

    if 6 in mesh.cscg.elements:
        c0 = mesh(6)
        c0.do.refine()

    mesh.do.lock()
    mesh.visualize()