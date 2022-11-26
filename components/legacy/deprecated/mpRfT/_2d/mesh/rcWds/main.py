# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/18 5:10 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly
from objects.mpRfT._2d.mesh.rcWds.helpers.scalar.main import mpRfT2_Mesh_rcWds_Scalar
from objects.mpRfT._2d.mesh.rcWds.helpers.vector.main import mpRfT2_Mesh_rcWds_Vector


class mpRfT2_Mesh_RootCellWiseDataStructure(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._freeze_self_()

    @property
    def scalar(self):
        return mpRfT2_Mesh_rcWds_Scalar(self._mesh_)

    @property
    def vector(self):
        return mpRfT2_Mesh_rcWds_Vector(self._mesh_)


if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/mesh/rcWds/main.py
    pass
