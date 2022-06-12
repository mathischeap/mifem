# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/19 8:38 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from  objects.mpRfT._2d.mesh.space.helpers.s0f import mpRfT2_Mesh_Space_S0F_Basis
from  objects.mpRfT._2d.mesh.space.helpers.so1f import mpRfT2_Mesh_Space_So1F_Basis
from  objects.mpRfT._2d.mesh.space.helpers.s2f import mpRfT2_Mesh_Space_S2F_Basis
from  objects.mpRfT._2d.mesh.space.helpers.nsg import mpRfT2_Mesh_Space_NSgF_Basis


class mpRfT2_Mesh_Space_do(FrozenOnly):
    """"""

    def __init__(self, space):
        """"""
        self._space_ = space
        self._freeze_self_()


    def evaluate_basis(self, f, coo_map):
        """"""
        assert coo_map.___Pr_is_mpRfT2_mesh_coo_map___, f"I need a coo distribution object."

        if f.__class__.__name__ in ('mpRfT2_Si0F', 'mpRfT2_So0F'):
            basis = mpRfT2_Mesh_Space_S0F_Basis(self._space_._mesh_, coo_map)
        elif f.__class__.__name__ == 'mpRfT2_So1F':
            basis = mpRfT2_Mesh_Space_So1F_Basis(self._space_._mesh_, coo_map)
        elif f.__class__.__name__ in ('mpRfT2_Si2F', 'mpRfT2_So2F'):
            basis = mpRfT2_Mesh_Space_S2F_Basis(self._space_._mesh_, coo_map)
        elif f.__class__.__name__ == 'mpRfT2_NSgF':
            basis = mpRfT2_Mesh_Space_NSgF_Basis(f, self._space_._mesh_, coo_map)
        else:
            raise NotImplementedError(f"no basis for {f.__class__.__name__}")

        return basis












if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
