# -*- coding: utf-8 -*-
from screws.freeze.base import FrozenOnly


class _3dCSCG_ADT_TF_Matrices(FrozenOnly):
    """"""
    def __init__(self, adt):
        self._adt_ = adt
        self._freeze_self_()


    @property
    def trace(self):
        """Return the trace matrix."""
        return self._adt_.coboundary.trace_matrix