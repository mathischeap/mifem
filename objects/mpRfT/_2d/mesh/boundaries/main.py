# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/30 5:54 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class mpRfT2_Mesh_Boundaries(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._names_ = mesh.cscg.boundaries.names
        self._freeze_self_()

    @property
    def names(self):
        return self._names_


if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/mesh/boundaries/main.py
    pass
