


from root.config import *
from SCREWS.frozen import FrozenOnly
from _3dCSCG.form.trace.base.dofs.dof.basis_function import _3dCSCG_TF_DOF_BF

from _3dCSCG.form.trace.base.dofs.dof.visualize import _3dCSCG_Trace_forms_DOF_VISUALIZE


class _3dCSCG_Trace_forms_DOF(FrozenOnly):
    """A dof of a trace form."""
    def __init__(self, dofs, i):
        """"""
        # we first check if dof #i is a local dof, if not, raise Error.
        ELEMENTS, INDICES = dofs.___PRIVATE_FIND_local_mesh_elements_and_local_indices_of_dof___(i)
        self._local_positions_ = list()
        for E, I in zip(ELEMENTS, INDICES):
            self._local_positions_.append((E, I))
        self._i_ = i # I am the #i dof.
        self._dofs_ = dofs
        self._tf_ = dofs._tf_
        self._bf_ = None
        self._visualize_ = None
        self._GLOBAL_positions_ = None
        self._freeze_self_()

    @property
    def i(self):
        """I am the ith dof (I am numbered i in the gathering matrix)."""
        return self._i_

    @property
    def positions(self):
        """Return a list of tuples which represent the "LOCAL" positions. For example,
            positions = [[(3, 4), (4, 3), (5, 2), (6, 1), (7, 0)]]
        Then we know, this dof is at GM[3][4], GM[4][3], GM[5][2], GM[6][1], GM[7][0]. And mesh elements
        3, 4, 5, 6, 7 are all in this core.

        For each item of this list, for example, (5,2) means this dof is at mesh element #5, and the local
        numbering of this dof is 2, so the third local numbering. We then can identify where is it is
        according to the degree of the space and the type (k) of the form.

        If
            positions = []
        Then we know this dof has is not in this core.

        """
        return self._local_positions_

    @property
    def GLOBAL_positions(self):
        """The "GLOBAL" positions of this dof. So if it is shared by multiple cores, we return all its
        positions. The positions are indicated in the same way as the local positions, see `positions`."""
        if self._GLOBAL_positions_ is None:
            positions = self.positions
            positions = cOmm.gather(positions, root=mAster_rank)
            if rAnk == mAster_rank:
                GP = list()
                for PS in positions:
                    GP.extend(PS)
            else:
                GP = None
            self._GLOBAL_positions_ = cOmm.bcast(GP, root=mAster_rank)
        return self._GLOBAL_positions_

    @property
    def basis_function(self):
        """The local basis function(s) of this dof."""
        if self._bf_ is None:
            self._bf_ = _3dCSCG_TF_DOF_BF(self)
        return self._bf_


    @property
    def visualize(self):
        if self._visualize_ is None:
            self._visualize_ = _3dCSCG_Trace_forms_DOF_VISUALIZE(self)
        return self._visualize_
