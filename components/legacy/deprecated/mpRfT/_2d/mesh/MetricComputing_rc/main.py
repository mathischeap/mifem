# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 5/24/2022 9:00 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly
from objects.mpRfT._2d.mesh.MetricComputing_rc.helpers.inverse_Jacobian_matrix import mpRfT2_Mesh_rcMC_iJM
from objects.mpRfT._2d.mesh.MetricComputing_rc.helpers.inverse_Jacobian import mpRfT2_Mesh_rcMC_inverseJacobian
from objects.mpRfT._2d.mesh.MetricComputing_rc.helpers.Jacobian import mpRfT2_Mesh_rcMC_Jacobian
from objects.mpRfT._2d.mesh.MetricComputing_rc.helpers.inverse_metric_matrix import mpRfT2_Mesh_rcMC_iMM


class mpRfT2_Mesh_rcMC(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._freeze_self_()

    @property
    def Jacobian(self):
        return mpRfT2_Mesh_rcMC_Jacobian(self._mesh_)

    @property
    def inverse_Jacobian_matrix(self):
        return mpRfT2_Mesh_rcMC_iJM(self._mesh_)

    @property
    def inverse_Jacobian(self):
        return mpRfT2_Mesh_rcMC_inverseJacobian(self._mesh_)

    @property
    def inverse_metric_matrix(self):
        return mpRfT2_Mesh_rcMC_iMM(self._mesh_)





if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
