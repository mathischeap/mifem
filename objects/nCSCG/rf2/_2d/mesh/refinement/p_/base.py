# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/12 10:58 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class _2nCSCG_RF2_p_Refinement_Base(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._signature_ = mesh.signature
        self._freeze_self_()

    @property
    def mesh(self):
        return self._mesh_

    @property
    def signature(self):
        return self._signature_

    def ___Pr_apply___(self):
        """"""
        raise NotImplementedError()







if __name__ == "__main__":
    # mpiexec -n 4 python objects/nCSCG/rfT2/_2d/mesh/refinement/p_/base/main.py
    pass
