# -*- coding: utf-8 -*-

from components.freeze.base import FrozenOnly


class _1Form_Outer_Special(FrozenOnly):
    def __init__(self, _1sf):
        self._sf_ = _1sf
        self._freeze_self_()
