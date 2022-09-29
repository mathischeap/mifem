# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/22 6:11 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.mpRfT._2d.mesh.coo_map.helpers.base import mpRfT2_CooMapBase
import numpy as np



class mpRfT2_Mesh_CustomizedCooMap(mpRfT2_CooMapBase):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._coo_ = dict()
        self._freeze_self_()







if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
