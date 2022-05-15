# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/12 10:58 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class _2nCSCG_RF2_p_Refinement(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._freeze_self_()

    @property
    def mesh(self):
        return self._mesh_






if __name__ == "__main__":
    # mpiexec -n 4 python objects/nCSCG/rf2/_2d/refinement/p_/main.py
    pass
