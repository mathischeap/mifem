# -*- coding: utf-8 -*-
from screws.freeze.main import FrozenOnly

class ___LinearSystem_Condition___(FrozenOnly):
    """Used to define customizations to A and b simultaneously."""
    def __init__(self, ls):
        self._LS_ = ls
        self._freeze_self_()

    @property
    def condition_number(self):
        """The condition number of the A matrix"""
        return self._LS_.A.assembled.condition.condition_number

    @property
    def shape(self):
        return self._LS_.A.assembled.shape