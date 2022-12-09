# -*- coding: utf-8 -*-
from components.freeze.main import FrozenOnly
from tools.elementwiseCache.dataStructures.objects.sparseMatrix.main import EWC_SparseMatrix


class _3dCSCG_Algebra_DUAL_Trace_Form_Coboundary(FrozenOnly):
    """This is one of the key properties of a standard algebraic dual form. To perform the coboundary, now we need the
    help of another dual trace form.
    """
    def __init__(self, dt):
        self._dt_ = dt
        self._T_ = None
        self._freeze_self_()

    @property
    def trace_matrix(self):
        """"""
        if self._T_ is None:
            if self._dt_.k == 2:
                formName = '_3dCSCG_2Trace'
            elif self._dt_.k == 1:
                formName = '_3dCSCG_1Trace'
            elif self._dt_.k == 0:
                formName = '_3dCSCG_0Trace'
            else:
                raise Exception()
            T = getattr(self._dt_.space.trace_matrix, formName)[0][self._dt_.prime.whether.hybrid]
            self._T_ = \
                EWC_SparseMatrix(self._dt_.mesh.elements, T, 'constant')
        return self._T_