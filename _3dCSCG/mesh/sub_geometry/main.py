

"""
Sub-geometries of a mesh. Like we can pick a point, a slice or a volume from a mesh. Such a sub-geometry should be
consist of sub-geometries of one or several elements.

"""


from SCREWS.frozen import FrozenOnly
from _3dCSCG.mesh.sub_geometry.perpendicular_slice import _3dCSCG_MeshPerpendicularSlice


class _3dCSCG_Mesh_SubGeometry(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh


        self._freeze_self_()





    def GENERATE_perpendicular_slice_object(self, *args, **kwargs):
        return _3dCSCG_MeshPerpendicularSlice(self._mesh_, *args, **kwargs)