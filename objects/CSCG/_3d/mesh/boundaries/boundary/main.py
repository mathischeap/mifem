# -*- coding: utf-8 -*-
"""
"""
import sys
if './' not in sys.path:
    sys.path.append('./')

from components.freeze.main import FrozenOnly
from objects.CSCG._3d.mesh.boundaries.boundary.visualize import _3dCSCG_Mesh_Boundary_VIS
from objects.CSCG._3d.mesh.boundaries.boundary.IS import _3dCSCG_MeshBoundaryIs
from objects.CSCG._3d.mesh.boundaries.boundary.coordinate_transformation.main import _3dCSCG_MeshBoundaryCT


class _3dCSCG_Mesh_Boundary(FrozenOnly):
    def __init__(self, bdrs, name):
        self._bdrs_ = bdrs
        self._name_ = name
        self._visualize_ = None
        self._ct_ = _3dCSCG_MeshBoundaryCT(self)
        self._IS_ = _3dCSCG_MeshBoundaryIs(self)
        self._freeze_self_()

    @property
    def coordinate_transformation(self):
        return self._ct_

    @property
    def IS(self):
        return self._IS_

    @property
    def name(self):
        return self._name_

    @property
    def mesh(self):
        return self._bdrs_._mesh_

    @property
    def element_sides(self):
        """This mesh boundary covers these local mesh element sides."""
        return self._bdrs_.range_of_element_sides[self._name_]

    @property
    def trace_elements(self):
        """This mesh boundary covers these local trace elements."""
        return self._bdrs_.range_of_trace_elements[self._name_]

    @property
    def region_sides(self):
        """This mesh boundary globally covers these region sides."""
        return self._bdrs_.range_of_region_sides[self._name_]

    @property
    def visualize(self):
        if self._visualize_ is None:
            self._visualize_ = _3dCSCG_Mesh_Boundary_VIS(self)
        return self._visualize_


if __name__ == '__main__':
    # mpiexec -n 8 python _3dCSCG\mesh\boundaries\boundary\main.py
    from objects.CSCG._3d.master import MeshGenerator
    elements = [5, 5, 5]
    # mesh = MeshGenerator('crazy', c=0.0, bounds=([0,1], [0,1], [0,1]))(elements)
    mesh = MeshGenerator('bridge_arch_cracked')(elements)
    boundaries = mesh.boundaries
    boundary = boundaries['Bottom']

    boundary.visualize()
