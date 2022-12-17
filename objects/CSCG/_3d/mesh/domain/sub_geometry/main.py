# -*- coding: utf-8 -*-
import sys
if './' not in sys.path:
    sys.path.append('./')

from components.freeze.main import FrozenOnly
from objects.CSCG._3d.mesh.domain.sub_geometry.perpendicular_slice import _3dCSCG_DomainPerpendicularSlice


class _3dCSCG_DomainSubGeometry(FrozenOnly):
    """"""

    def __init__(self, domain):
        """"""
        self._domain_ = domain
        self._3dCSCG_DomainPerpendicularSlice_ = _3dCSCG_DomainPerpendicularSlice

        self._freeze_self_()

    def make_a_perpendicular_slice_object_on(self, x=None, y=None, z=None):
        """"""
        num_None = 0
        for _ in (x, y, z):
            if _ is None:
                num_None += 1
        assert num_None == 2, f"x={x}, y={y}, z={z} wrong, only one of them must not be None."

        return self._3dCSCG_DomainPerpendicularSlice_(self._domain_, x, y, z)


if __name__ == '__main__':
    # mpiexec -n 8 python objects/CSCG/_3d/mesh/domain/sub_geometry/main.py

    from objects.CSCG._3d.master import MeshGenerator

    mesh = MeshGenerator('cuboid', region_layout=(2, 2, 2))([3, 3, 3], EDM='chaotic', show_info=True)

    DSG = mesh.domain.sub_geometry

    DSG.make_a_perpendicular_slice_object_on(x=0.25)
