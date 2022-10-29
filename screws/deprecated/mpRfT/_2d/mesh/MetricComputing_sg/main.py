# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/01 3:12 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.mpRfT._2d.mesh.MetricComputing_sg.helpers.Jacobian import mpRfT2_Mesh_sgMC_Jacobian


class mpRfT2_Mesh_sgMC(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._freeze_self_()

    @property
    def Jacobian(self):
        return mpRfT2_Mesh_sgMC_Jacobian(self._mesh_)




if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/mesh/sgMetricComputing/main.py
    pass
