# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly


class EWC_SparseMatrix_Whether(FrozenOnly):
    """"""
    def __init__(self, MAT):
        """"""
        self._MAT_ = MAT
        self._freeze_self_()

    def assembled_matrix_locked(self):
        """If the assembled matrix locked?"""
        return self._MAT_.do.___locker___
