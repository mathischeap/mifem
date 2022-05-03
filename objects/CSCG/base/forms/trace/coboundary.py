# -*- coding: utf-8 -*-
from screws.freeze.main import FrozenOnly
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.main import EWC_SparseMatrix




class CSCG_Trace_Form_Coboundary_BASE(FrozenOnly):
    """"""
    def __init__(self, tf):
        self._tf_ = tf
        self._T_ = None
        self.___PRIVATE_reset_cache___()
        self._freeze_self_()

    def ___PRIVATE_reset_cache___(self):
        pass

    @property
    def trace_matrix(self):
        if self._T_ is None:
            formName = self._tf_.__class__.__name__
            T = getattr(self._tf_.space.trace_matrix, formName)[0]
            self._T_ = \
                EWC_SparseMatrix(self._tf_.mesh.elements, T, 'constant')
        return self._T_