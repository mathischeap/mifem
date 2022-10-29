# -*- coding: utf-8 -*-
"""
Sub-geometries of a mesh. Like we can pick a point, a slice or a volume from a mesh. Such a sub-geometry should
consist of sub-geometries of one or several elements.

"""

import sys
if './' not in sys.path: sys.path.append('./')

from screws.freeze.main import FrozenOnly
from objects.CSCG._3d.mesh.subGeometry.mesh_perpendicular_slice import _3dCSCG_MeshPerpendicularSlice



class _3dCSCG_Mesh_SubGeometry(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._MeshPerpendicularSlice_ = _3dCSCG_MeshPerpendicularSlice
        self._freeze_self_()



    def make_a_perpendicular_slice_object_on(self, *args, **kwargs):
        return self._MeshPerpendicularSlice_(self._mesh_, *args, **kwargs)




if __name__ == '__main__':
    # mpiexec -n 8 python objects/CSCG/_3d/mesh/sub_geometry/main.py

    from objects.CSCG._3d.master import MeshGenerator

    mesh = MeshGenerator('cuboid', region_layout=(2,2,2))([3, 3, 3], EDM='chaotic', show_info=True)

    MSG = mesh.sub_geometry

    MPS = MSG.make_a_perpendicular_slice_object_on(y=0.25, x=None, z=None)


    for i in MPS:
        print(i, MPS[i], MPS.perpendicular_to_axis)
