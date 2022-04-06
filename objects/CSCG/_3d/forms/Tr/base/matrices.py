from screws.freeze.main import FrozenOnly
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.main import EWC_SparseMatrix




class _3dCSCG_Tr_Matrices(FrozenOnly):
    def __init__(self, Tr):
        self._Tr_ = Tr
        self._T_ = None
        self._S_ = None
        self._mass_TEW_ = None
        self._freeze_self_()

    @property
    def trace(self):
        """EWC_SparseMatrix: Return the trace matrix.
        """
        if self._T_ is None:
            k = self._Tr_.k
            formName = f'_3dCSCG_{int(k)}Tr'
            T = getattr(self._Tr_.space.trace_matrix, formName)[0]
            self._T_ = \
                EWC_SparseMatrix(self._Tr_.mesh.elements, T, 'constant')
        return self._T_

    @property
    def selective(self):
        """EWC_SparseMatrix: Return the selective (mesh-element -> trace element) matrix.

        Like the trace matrix but without minus sign.
        """
        if self._S_ is None:
            k = self._Tr_.k
            formName = f'_3dCSCG_{int(k)}Tr'
            S = getattr(self._Tr_.space.selective_matrix, formName)[0]
            self._S_ = \
                EWC_SparseMatrix(self._Tr_.mesh.elements, S, 'constant')
        return self._S_

    @property
    def mass_TEW(self):
        """EWC_SparseMatrix: The mesh-element-wise mass matrix.
        """
        if self._mass_TEW_ is None:
            self._mass_TEW_ = self._Tr_.___PRIVATE_generate_TEW_mass_matrices___()
        return self._mass_TEW_