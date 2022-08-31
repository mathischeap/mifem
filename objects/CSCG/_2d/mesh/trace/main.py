# -*- coding: utf-8 -*-

from screws.freeze.main import FrozenOnly
from objects.CSCG._2d.mesh.trace.visualize import _2dCSCG_Trace_Visualize
from objects.CSCG._2d.mesh.trace.elements.main import _2dCSCG_Trace_Elements



class _2dCSCG_Trace(FrozenOnly):
    def __init__(self, mesh):
        self._mesh_ = mesh
        self._elements_ = _2dCSCG_Trace_Elements(self)
        self._visualize_ = _2dCSCG_Trace_Visualize(self)
        self._freeze_self_()

    def RESET_cache(self):
        self.elements.RESET_cache()

    @property
    def elements(self):
        return self._elements_

    @property
    def visualize(self):
        return self._visualize_