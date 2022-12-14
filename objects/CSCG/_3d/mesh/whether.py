# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 12/12/2022 10:50 PM
"""
from components.freeze.main import FrozenOnly


class _3dCSCG_Mesh_Whether(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._freeze_self_()

    @property
    def orthogonal(self):
        return self._mesh_.elements.whether.all_orthogonal
