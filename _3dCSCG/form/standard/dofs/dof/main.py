

import sys
if './' not in sys.path: sys.path.append('./')
from root.config import rAnk

from SCREWS.frozen import FrozenOnly
from _3dCSCG.form.standard.dofs.dof.basis_function import _3dCSCG_SF_DOF_BF
from _3dCSCG.form.standard.dofs.dof.visualize.main import _3dCSCG_SF_DOF_VISUALIZE




class _3dCSCG_Standard_forms_DOF(FrozenOnly):
    """A dof of a standard form."""
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
        self._sf_ = dofs._sf_
        self._visualize_ = None
        self._bf_ = None
        self._freeze_self_()


    @property
    def positions(self):
        """Return a list of tuples which represent the "LOCAL" positions. For example,
            local_positions = [[(3, 4), (4, 3), (5, 2), (6, 1), (7, 0)]]
        Then we know, this dof is at GM[3][4], GM[4][3], GM[5][2], GM[6][1], GM[7][0].

        For each item of this list, for example, (5,2) means this dof is at mesh element #5, and the local
        numbering of this dof is 2, so the third local numbering. We then can identify where is it is
        according to the degree of the space and the type (k) of the form.

        """
        return self._local_positions_

    @property
    def GLOBAL_positions(self):
        """The "GLOBAL" positions of this dof. So if it is shared by multiple cores, we return all its
        positions. The positions are indicated in the same way as local positions, see `positions`."""
        raise NotImplementedError()


    @property
    def visualize(self):
        """To visualize this dof."""
        if self._visualize_ is None:
            self._visualize_ = _3dCSCG_SF_DOF_VISUALIZE(self)
        return self._visualize_

    @property
    def basis_function(self):
        """The local basis function(s) of this dof."""
        if self._bf_ is None:
            self._bf_ = _3dCSCG_SF_DOF_BF(self)
        return self._bf_






