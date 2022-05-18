# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 12:27 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class _2nCSCG_RF2_h_Refinement(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._freeze_self_()

    @property
    def mesh(self):
        return self._mesh_






if __name__ == '__main__':
    # mpiexec -n 4 python objects/nCSCG/rfT2/_2d/refinement/h_/main.py
    pass
