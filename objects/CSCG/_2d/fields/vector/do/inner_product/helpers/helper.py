# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/8/11 0:30
"""





class VF2_IP_H(object):
    def __init__(self, f0, f1, f2, f3):
        self._f0_ = f0
        self._f1_ = f1
        self._f2_ = f2
        self._f3_ = f3

    def __call__(self, t, x, y):
        return self._f0_(t, x, y) * self._f2_(t, x, y) + self._f1_(t, x, y) * self._f3_(t, x, y)

