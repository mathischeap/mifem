

"""
Sub-geometries of a mesh. Like we can pick a point, a slice or a volume from a mesh. Such a sub-geometry should
consist of sub-geometries of one or several elements.

"""

from screws.freeze.main import FrozenOnly
from _3dCSCG.mesh.sub_geometry.mesh_perpendicular_slice import _3dCSCG_MeshPerpendicularSlice



class _3dCSCG_Mesh_SubGeometry(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._MeshPerpendicularSlice_ = _3dCSCG_MeshPerpendicularSlice
        self._freeze_self_()



    def make_a_perpendicular_slice_object_on(self, *args, **kwargs):
        return self._MeshPerpendicularSlice_(self._mesh_, *args, **kwargs)
