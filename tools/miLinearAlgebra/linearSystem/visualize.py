# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 10/2/2022 7:14 PM
"""
from components.freeze.main import FrozenOnly


class LinearSystem_Visualize(FrozenOnly):
    """"""

    def __init__(self, ls):
        """"""
        self._ls_ = ls
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        return self.spy(*args, **kwargs)

    def spy(self, *args, **kwargs):
        """"""
        A = self._ls_.A.assembled
        return A.visualize.spy(*args, **kwargs)
