
from screws.freeze.main import FrozenOnly
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.main import EWC_SparseMatrix





class _3dCSCG_Standard_Form_Matrices(FrozenOnly):
    def __init__(self, sf):
        self._sf_ = sf
        self._T_ = None
        self._S_ = None
        self._freeze_self_()

    @property
    def mass(self):
        """(Dict[int, scipy.sparse.csr_matrix]) The mass matrix."""
        M = self._sf_.operators.inner(self._sf_)
        # note that even all mesh elements are unique, we still cache the local mass matrices because we may use them for multiple times.
        M.gathering_matrices = (self._sf_, self._sf_)
        return M

    @property
    def incidence(self):
        return self._sf_.coboundary.incidence_matrix

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
