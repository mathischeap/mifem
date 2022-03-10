from screws.freeze.inheriting.frozen_only import FrozenOnly





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

