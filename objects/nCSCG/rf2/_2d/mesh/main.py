# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/04 2:56 PM
"""
import sys
if './' not in sys.path: sys.path.append('./')

from objects.nCSCG.rf2.base.mesh.base import nCSCG_RF2_MeshBase
from objects.nCSCG.rf2._2d.mesh.base_cells import _2nCSCG_RF2_BaseMeshCells
from objects.nCSCG.rf2._2d.mesh.do.main import _2nCSCG_MeshDo
from objects.nCSCG.rf2._2d.mesh.visualize.main import _2nCSCG_MeshVisualize


class _2nCSCG_RF2_Mesh(nCSCG_RF2_MeshBase):
    """"""
    def __init__(self, cscg):
        """"""
        super(_2nCSCG_RF2_Mesh, self).__init__(cscg)
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






if __name__ == '__main__':
    # mpiexec -n 4 python objects/nCSCG/rf2/_2d/mesh/main.py
    from objects.nCSCG.rf2._2d.master import MeshGenerator

    mesh = MeshGenerator('crazy')([3, 3], EDM='chaotic', show_info=True)
    if 0 in mesh.cscg.elements:
        c0 = mesh(0)
        c0.do.refine()

    mesh.visualize()