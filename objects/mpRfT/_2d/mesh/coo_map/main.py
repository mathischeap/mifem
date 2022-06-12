# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/18 4:32 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.mpRfT._2d.mesh.coo_map.helpers.uniform import mpRfT2_Mesh_UniformCooMap
from objects.mpRfT._2d.mesh.coo_map.helpers.Gauss import mpRfT2_Mesh_GaussCooMap
from objects.mpRfT._2d.mesh.coo_map.helpers.Lobatto import mpRfT2_Mesh_LobattoCooMap




class mpRfT2_Mesh_CooMap(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._freeze_self_()

    @property
    def uniform(self):
        return mpRfT2_Mesh_UniformCooMap(self._mesh_)

    @property
    def Gauss(self):
        return mpRfT2_Mesh_GaussCooMap(self._mesh_)

    @property
    def Lobatto(self):
        return mpRfT2_Mesh_LobattoCooMap(self._mesh_)






if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/mesh/coo_map/main.py
    pass
