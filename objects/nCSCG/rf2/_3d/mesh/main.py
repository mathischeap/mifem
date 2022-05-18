# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/04 2:56 PM
"""
import sys
if './' not in sys.path: sys.path.append('./')

from objects.nCSCG.rf2.base.mesh.main import nCSCG_RF2_MeshBase
from objects.nCSCG.rf2._3d.mesh.base_cells import _3nCSCG_RF2_BaseMeshCells
from objects.nCSCG.rf2._3d.mesh.do.main import _3nCSCG_MeshDo
from objects.nCSCG.rf2._3d.mesh.visualize.main import _3nCSCG_MeshVisualize
from objects.nCSCG.rf2._3d.mesh.facets.main import _3nCSCG_MeshFacets
from objects.nCSCG.rf2._3d.mesh.IS import _3nCSCG_Mesh_RF2_IS
from objects.nCSCG.rf2._3d.mesh.space.allocator import _3nCSCG_SpaceAllocator


class _3nCSCG_RF2_Mesh(nCSCG_RF2_MeshBase):
    """"""
    def __init__(self, cscg, dN, space_type='polynomials', **space_kwargs):
        """"""
        super(_3nCSCG_RF2_Mesh, self).__init__(cscg)
        self._do_ = None
        self._visualize_ = None
        self._facets_ = None
        self._IS_ = None
        self._freeze_self_()
        self.___base_mesh_cells___ = _3nCSCG_RF2_BaseMeshCells(self)
        self._space_ = _3nCSCG_SpaceAllocator(space_type)(dN, **space_kwargs)
        self._space_._mesh_ = self

    @property
    def do(self):
        if self._do_ is None:
            self._do_ = _3nCSCG_MeshDo(self)
        return self._do_

    @property
    def visualize(self):
        if self._visualize_ is None:
            self._visualize_ = _3nCSCG_MeshVisualize(self)
        return self._visualize_

    @property
    def facets(self):
        if self._facets_ is None:
            self._facets_ = _3nCSCG_MeshFacets(self)
        return self._facets_

    @property
    def IS(self):
        if self._IS_ is None:
            self._IS_ = _3nCSCG_Mesh_RF2_IS(self)
        return self._IS_



if __name__ == '__main__':
    # mpiexec -n 4 python objects/nCSCG/rfT2/_3d/mesh/main.py
    from objects.nCSCG.rf2._3d.master import MeshGenerator

    mesh = MeshGenerator('crazy')([3, 3, 3], 2, EDM='chaotic', show_info=True)

    print(mesh.space.mesh is mesh, mesh.space._PRM)