# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/22/2022 12:26 AM
"""

from components.freeze.base import FrozenOnly


class _2dCSCG_ScalarNeg(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._freeze_self_()

    def __call__(self, t, x, y):
        return - self._f_(t, x, y)
