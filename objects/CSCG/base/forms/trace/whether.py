# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly


class CSCG_trace_form_Whether(FrozenOnly):
    """"""
    def __init__(self, tf):
        self._tf_ = tf
        self._freeze_self_()

    @property
    def hybrid(self):
        """A hybrid trace form or not."""
        return self._tf_._hybrid_

    @property
    def inner_oriented(self):
        return True if self._tf_._orientation_ == 'inner' else False
