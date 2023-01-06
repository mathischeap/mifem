# -*- coding: utf-8 -*-

from components.freeze.base import FrozenOnly
from tools.elementwiseCache.dataStructures.objects.sparseMatrix.main import EWC_SparseMatrix
from objects.CSCG._2d.forms.standard.base.matrices.helpers.mass import MassMatrixHelper


class _2dCSCG_Standard_Form_Matrices(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._freeze_self_()

    @property
    def mass(self):
        """(Dict[int, scipy.sparse.csr_matrix]) The mass matrix."""
        data_generator = MassMatrixHelper(self._sf_)

        if self._sf_.mesh.elements.whether.homogeneous_according_to_types_wrt_metric:
            return EWC_SparseMatrix(
                self._sf_.mesh.elements,
                data_generator,
                'constant',
            )
        else:
            return EWC_SparseMatrix(
                self._sf_.mesh.elements,
                data_generator,
                self._sf_.mesh.elements.___PRIVATE_elementwise_cache_metric_key___,
            )

    @property
    def incidence(self):
        return self._sf_.coboundary.incidence_matrix

    @property
    def identity(self):
        """Return a identity matrix of local shape equal to the mass matrix; (#local dofs, l#ocal dofs)."""
        return EWC_SparseMatrix(self._sf_.mesh, ('identity', self._sf_.num.basis))
