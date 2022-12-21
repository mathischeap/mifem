# -*- coding: utf-8 -*-
from components.freeze.main import FrozenOnly


class ___CV_SUB___(FrozenOnly):
    def __init__(self, v1, v2):
        self._v1_ = v1
        self._v2_ = v2
        self._freeze_self_()

    def __call__(self, i):
        """"""
        return self._v1_[i] - self._v2_[i]

    def __KG_call__(self, i):
        """"""
        return self._v1_._KG_(i) + self._v1_._KG_(i)
