# -*- coding: utf-8 -*-

from importlib import import_module
from objects.CSCG.base.forms.standard.coboundary import CSCG_Standard_Form_Coboundary_BASE
from tools.elementwiseCache.dataStructures.objects.sparseMatrix.main import EWC_SparseMatrix


class _2dCSCG_Standard_Form_Coboundary(CSCG_Standard_Form_Coboundary_BASE):
    def __init__(self, sf):
        super().__init__(sf)

    @property
    def incidence_matrix(self):
        """(scipy.sparse.csc_matrix) Return ths incidence matrix of the standard form."""
        assert self._sf_.k < 2, "volume form has no incidence matrix."
        if self._incidenceMatrix_ is None:
            formName = self._sf_.__class__.__name__
            _incidenceMatrix_ = getattr(self._sf_.space.incidence_matrix, formName)
            self._incidenceMatrix_ = \
                EWC_SparseMatrix(self._sf_.mesh.elements, _incidenceMatrix_, 'constant')
        return self._incidenceMatrix_

    def ___PRIVATE_next_class___(self):
        assert self._sf_.k < 2, "volume form has no next prime space."
        k = self._sf_.k
        base_path = '.'.join(str(self.__class__).split(' ')[1][1:-2].split('.')[:-3]) + '.'
        nextPath = base_path + f'_{k+1}_form.' + self._sf_.orientation + '.main'
        nextName = f'_2dCSCG_{k+1}Form_'
        if self._sf_.orientation == 'inner':
            nextName += 'Inner'
        elif self._sf_.orientation == 'outer':
            nextName += 'Outer'
        else:
            raise Exception()
        return getattr(import_module(nextPath), nextName)
