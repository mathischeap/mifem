# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 4:22 PM
"""
from components.freeze.base import FrozenOnly


class miUs_Triangular_SF_Num(FrozenOnly):
    """"""

    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self._freeze_self_()

    @property
    def basis(self):
        return getattr(self._sf_.space.num_basis, self._sf_.__class__.__name__)

    @property
    def GLOBAL_dofs(self):
        return self._sf_.numbering.gathering.global_num_dofs