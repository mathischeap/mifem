# -*- coding: utf-8 -*-
from components.freeze.main import FrozenOnly
from tools.elementwiseCache.dataStructures.objects.sparseMatrix.main import EWC_SparseMatrix
from objects.CSCG._3d.forms.standard.base.matrices.helpers.mass import MassMatrixHelper

class _3dCSCG_Standard_Form_Matrices(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._T_ = None
        self._S_ = None
        self._freeze_self_()

    @property
    def mass(self):
        """(Dict[int, scipy.sparse.csr_matrix]) The mass matrix."""
        data_generator = MassMatrixHelper(self._sf_)

        if self._sf_.mesh.elements.whether.homogeneous_according_to_types_wrt_metric:
            M = EWC_SparseMatrix(
                self._sf_.mesh.elements,
                data_generator,
                'constant',
            )

        else:
            # note that even all mesh elements are unique, we still cache the output because we may use it for multiple times.
            M = EWC_SparseMatrix(
                self._sf_.mesh.elements,
                data_generator,
                self._sf_.mesh.elements.___PRIVATE_elementwise_cache_metric_key___,
            )

        M.gathering_matrices = (self._sf_, self._sf_)

        return M

    @property
    def incidence(self):
        return self._sf_.coboundary.incidence_matrix

    @property
    def identity(self):
        """Return an identity matrix of local shape equal to the mass matrix; (#local dofs, l#ocal dofs)."""
        return EWC_SparseMatrix(self._sf_.mesh, ('identity', self._sf_.num.basis))

    @property
    def trace(self):
        """Return the trace matrix."""
        if self._T_ is None:
            k = self._sf_.k
            formName = f'_3dCSCG_{int(k)}Trace'
            T = getattr(self._sf_.space.trace_matrix, formName)[0]
            self._T_ = \
                EWC_SparseMatrix(self._sf_.mesh.elements, T, 'constant')
        return self._T_
