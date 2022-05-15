# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/13 2:00 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.nCSCG.rf2._2d.mesh.refinement.h_.main import _2nCSCG_RF2_h_Refinement
from objects.nCSCG.rf2._2d.mesh.refinement.p_.main import _2nCSCG_RF2_p_Refinement



class _2nCSCG_Refinement(FrozenOnly):
    """"""
    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._freeze_self_()

    @property
    def h(self):
        return _2nCSCG_RF2_h_Refinement(self._mesh_)

    @property
    def p(self):
        return _2nCSCG_RF2_p_Refinement(self._mesh_)





if __name__ == "__main__":
    # mpiexec -n 4 python objects/nCSCG/rf2/_2d/mesh/refinement/main.py
    pass
