# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/13 2:00 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.nCSCG.rf2._2d.mesh.refinement.p_.allocator import _2nCSCG_RF2_Mesh_p_Refinement_Allocator



class _2nCSCG_Refinement(FrozenOnly):
    """"""
    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._freeze_self_()

    @property
    def p(self):
        return _2nCSCG_RF2_Mesh_p_Refinement_Allocator(self._mesh_)




if __name__ == "__main__":
    # mpiexec -n 4 python objects/nCSCG/rfT2/_2d/mesh/refinement/main.py
    from objects.nCSCG.rf2._2d.master import MeshGenerator

    mesh = MeshGenerator('crazy')([3, 3], 2, EDM='chaotic', show_info=True)
    pr = mesh.refinement.p('simple')

