# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly


class CSCG_standard_form_Whether(FrozenOnly):
    """"""
    def __init__(self, f):
        self._f_ = f
        self._freeze_self_()

    @property
    def hybrid(self):
        return self._f_._IS_hybrid_

    @property
    def inner_oriented(self):
        return True if self._f_._orientation_ == 'inner' else False

    @property
    def volume_form(self):
        return self._f_.ndim == self._f_.k