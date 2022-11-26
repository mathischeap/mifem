# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/08/26 2:31 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from components.freeze.base import FrozenOnly
from objects.CSCG._3d.mesh.boundaries.boundary.coordinate_transformation.constant import _3dCSCG_MeshBoundaryCT_constant


class _3dCSCG_MeshBoundaryCT(FrozenOnly):
    """"""

    def __init__(self, boundary):
        """"""
        self._boundary_ = boundary
        self._constant_ = _3dCSCG_MeshBoundaryCT_constant(boundary)
        self._freeze_self_()


    @property
    def constant(self):
        return self._constant_

if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
