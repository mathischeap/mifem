# -*- coding: utf-8 -*-

from screws.freeze.main import FrozenOnly

class _3Form_Special(FrozenOnly):
    def __init__(self, _3sf):
        self._sf_ = _3sf
        self._freeze_self_()