# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly

class ___MUL___(FrozenOnly):
    def __init__(self, ewc, number):
        assert isinstance(number, (int, float))
        self._ewc_ = ewc
        self._number_ = number
        self._freeze_self_()

    def __call__(self, item):
        return self._ewc_[item] * self._number_