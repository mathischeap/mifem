# -*- coding: utf-8 -*-


class ___VECTOR_DIVERGENCE_HELPER___(object):
    def __init__(self, f0, f1, f2):
        self._f0_ = f0
        self._f1_ = f1
        self._f2_ = f2

    def __call__(self, t, x, y, z):
        return self._f0_(t, x, y, z) + self._f1_(t, x, y, z) + self._f2_(t, x, y, z)
