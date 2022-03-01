# -*- coding: utf-8 -*-

from screws.frozen import FrozenOnly

from _2dCSCG.mesh.trace.visualize import _2dCSCG_Trace_Visualize
from _2dCSCG.mesh.trace.elements.main import _2dCSCG_Trace_Elements





class _2dCSCG_Trace(FrozenOnly):
    def __init__(self, mesh):
        self._mesh_ = mesh
        self._elements_ = _2dCSCG_Trace_Elements(self)
        self._visualize_ = _2dCSCG_Trace_Visualize(self)
        self.___PRIVATE_reset_cache___()
        self._freeze_self_()

    def ___PRIVATE_reset_cache___(self):
        self.elements.___PRIVATE_reset_cache___()

    @property
    def elements(self):
        return self._elements_

    @property
    def visualize(self):
        return self._visualize_

