# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/13/2022 4:54 PM
"""
from components.freeze.main import FrozenOnly


class t3d_3ScalarAdd(FrozenOnly):
    """"""

    def __init__(self, s0, s1, s2):
        """"""
        self._s0_ = s0
        self._s1_ = s1
        self._s2_ = s2
        self._freeze_self_()

    def __call__(self, t, x, y, z):
        return self._s0_(t, x, y, z) + self._s1_(t, x, y, z) + self._s2_(t, x, y, z)
