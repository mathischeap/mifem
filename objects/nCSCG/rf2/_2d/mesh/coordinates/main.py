# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/13 7:06 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.nCSCG.rf2._2d.mesh.coordinates.distributions.homogeneous import Homogeneous




class _2nCSCG_MeshRF2_Coordinates(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._freeze_self_()

    @property
    def homogeneous(self):
        return Homogeneous(self._mesh_)





if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
