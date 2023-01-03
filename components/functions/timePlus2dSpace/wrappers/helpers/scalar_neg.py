# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/13/2022 7:13 PM
"""
from components.freeze.main import FrozenOnly


class t2d_ScalarNeg(FrozenOnly):
    """"""

    def __init__(self, s):
        """"""
        self._s_ = s
        self._freeze_self_()

    def __call__(self, t, x, y):
        return - self._s_(t, x, y)
