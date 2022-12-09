# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 10/2/2022 10:30 PM
"""
from components.freeze.main import FrozenOnly


class miUs_Triangular_SF_whether(FrozenOnly):
    """"""

    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self._freeze_self_()

    @property
    def hybrid(self):
        return self._sf_._hybrid_

    @property
    def volume_form(self):
        return self._sf_.k == self._sf_.ndim

    @property
    def standard_form(self):
        return True