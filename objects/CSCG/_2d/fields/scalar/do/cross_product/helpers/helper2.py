# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/8/11 0:30
"""


class ___SF_CROSS_PRODUCT_HELPER_2___(object):
    """"""

    def __init__(self, f0, f1):
        self._f0_ = f0
        self._f1_ = f1

    def __call__(self, t, x, y):
        return self._f0_(t, x, y) * self._f1_(t, x, y)
