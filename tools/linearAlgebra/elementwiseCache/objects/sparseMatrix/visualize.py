# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly


class EWC_SparseMatrix_Vis(FrozenOnly):
    """"""
    def __init__(self, MAT):
        """"""
        self._MAT_ = MAT
        self._freeze_self_()

    def spy(self, *args, **kwargs):
        """"""
        self._MAT_.assembled.visualize.spy(*args, **kwargs)