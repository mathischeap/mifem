# -*- coding: utf-8 -*-

from screws.freeze.base import FrozenOnly
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.main import EWC_SparseMatrix




class _2dCSCG_Standard_Form_Matrices(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._freeze_self_()

    @property
    def mass(self):
        """(Dict[int, scipy.sparse.csr_matrix]) The mass matrix."""
        return self._sf_.operators.inner(self._sf_)

    @property
    def incidence(self):
        return self._sf_.coboundary.incidence_matrix

    @property
    def identity(self):
        """Return a identity matrix of local shape equal to the mass matrix; (#local dofs, l#ocal dofs)."""
        return EWC_SparseMatrix(self._sf_.mesh, ('identity', self._sf_.num.basis))