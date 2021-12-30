



from SCREWS.frozen import FrozenOnly
from root.config import rAnk
from _3dCSCG.form.trace.dofs.dof.basis_function import _3dCSCG_TF_DOF_BF

class _3dCSCG_Trace_forms_DOF(FrozenOnly):
    """A dof of a trace form."""
    def __init__(self, dofs, i):
        """"""
        # we first check if dof #i is a local dof, if not, raise Error.
        ELEMENTS, INDICES = dofs.___PRIVATE_FIND_local_mesh_elements_and_local_indices_of_dof___(i)
        assert len(ELEMENTS) > 0, f"dof #{i} is not a local dof in RANK {rAnk}."
        self._local_positions_ = list()
        for E, I in zip(ELEMENTS, INDICES):
            self._local_positions_.append((E, I))
        self._i_ = i # I am the #i dof.
        self._dofs_ = dofs
        self._tf_ = dofs._tf_
        self._bf_ = None
        self._freeze_self_()

    @property
    def basis_function(self):
        """The local basis function(s) of this dof."""
        if self._bf_ is None:
            self._bf_ = _3dCSCG_TF_DOF_BF(self)
        return self._bf_



