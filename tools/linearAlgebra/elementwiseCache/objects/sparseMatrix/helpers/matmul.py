# -*- coding: utf-8 -*-


from components.freeze.base import FrozenOnly



class ___MATMUL___(FrozenOnly):
    def __init__(self, EWC1, EWC2):
        self._ewc1_ = EWC1
        self._ewc2_ = EWC2
        self._freeze_self_()

    def __DG_call__(self, item):
        return self._ewc1_[item] @ self._ewc2_[item]

    def __KG_call__(self, item):
        return self._ewc1_._KG_(item) + self._ewc2_._KG_(item)