from screws.frozen import FrozenOnly
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.main import EWC_SparseMatrix




class _3dCSCG_Trace_Matrices(FrozenOnly):
    def __init__(self, tf):
        self._tf_ = tf
        self._T_ = None
        self._S_ = None
        self._mass_ = None
        self._freeze_self_()

    @property
    def trace(self):
        """Return the trace matrix."""
        if self._T_ is None:
            k = self._tf_.k
            formName = f'_3dCSCG_{int(k)}Trace'
            T = getattr(self._tf_.space.trace_matrix, formName)[0]
            self._T_ = \
                EWC_SparseMatrix(self._tf_.mesh.elements, T, 'constant')
        return self._T_

    @property
    def selective(self):
        """Return the selective (mesh-element -> trace element) matrix.

        Like the trace matrix but without minus sign.
        """
        if self._S_ is None:
            k = self._tf_.k
            formName = f'_3dCSCG_{int(k)}Trace'
            S = getattr(self._tf_.space.selective_matrix, formName)[0]
            self._S_ = \
                EWC_SparseMatrix(self._tf_.mesh.elements, S, 'constant')
        return self._S_

    @property
    def mass(self):
        """Return the mass matrix. It is a dict, keys are trace-element numbers, and values are the mass
        matrices on the trace elements."""
        if self._mass_ is None:
            self._mass_ = self._tf_.___PRIVATE_generate_TEW_mass_matrices___()
        return self._mass_
