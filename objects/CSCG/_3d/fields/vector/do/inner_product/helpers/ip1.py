# -*- coding: utf-8 -*-
"""
INTRO

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""




class ___VF_INNER_PRODUCT_HELPER_1___(object):
    def __init__(self, w0, w1, w2, u0, u1, u2):
        self._w0_ = w0
        self._w1_ = w1
        self._w2_ = w2
        self._u0_ = u0
        self._u1_ = u1
        self._u2_ = u2

    def __call__(self, t, x, y, z):
        return self._w0_(t, x, y, z) * self._u0_(t, x, y, z) + \
               self._w1_(t, x, y, z) * self._u1_(t, x, y, z) + \
               self._w2_(t, x, y, z) * self._u2_(t, x, y, z)