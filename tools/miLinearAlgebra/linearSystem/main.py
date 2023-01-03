# -*- coding: utf-8 -*-
"""We make a class to represent a linear system.
"""
from components.freeze.main import FrozenClass
from tools.miLinearAlgebra.linearSystem.customize import ___LinearSystem_Customize___
from tools.miLinearAlgebra.linearSystem.condition import ___LinearSystem_Condition___
from tools.miLinearAlgebra.linearSystem.solve.main import ___LinearSystem_Solve___
from tools.miLinearAlgebra.linearSystem.visualize import LinearSystem_Visualize


class LinearSystem(FrozenClass):
    """"""
    def __init__(self, A, b):
        """Ax = b.

        `A` and `b` are given as (local) blocks, we will bmat and concatenate them into (local) EWC_SparseMatrix and
        EWC_ColumnVector for the local elements during initializing the linear system.

        :param A: A 2-d list or tuple of matrix blocks.
        :param b: A 1-d list or tuple of vector blocks.
        """
        assert A.__class__.__name__ == 'EWC_SparseMatrix', "A type wrong, must be an EWC_SparseMatrix."
        assert b.__class__.__name__ == 'EWC_ColumnVector', "b type wrong, must be an EWC_ColumnVector."

        assert A.bmat_shape is not False, f"A must be a bmat EWC to protect the blocks"
        assert b.con_shape is not False, f"b must be a concatenate EWC to protect the blocks"

        assert A.bmat_shape[0] == b.con_shape[0], f"A: {A.bmat_shape}, b: {b.con_shape} shapes do not match"

        GMr, GMc = A.gathering_matrices
        assert GMr is not None, f"A has no row gathering_matrix."
        assert GMc is not None, f"A has no col gathering_matrix."
        assert b.gathering_matrix is not None, f"b has no gathering_matrix."
        assert b.gathering_matrix == GMr, f"row gathering matrices of A and b do not match."

        self._GMr_ = GMr
        self._GMc_ = GMc

        self._A_ = A
        self._b_ = b

        self._customize_ = ___LinearSystem_Customize___(self)
        self._condition_ = ___LinearSystem_Condition___(self)
        self._solve_ = ___LinearSystem_Solve___(self)
        self._visualize_ = LinearSystem_Visualize(self)

        self._freeze_self_()

    # ------- functional ---------------------------------------------------------------------------
    @property
    def customize(self):
        """Adjust Ax=b by adjusting EWC matrix A and EWC vector b."""
        return self._customize_

    @property
    def condition(self):
        """Get access to the conditions of this Linear System."""
        return self._condition_

    @property
    def solve(self):
        return self._solve_

    @property
    def visualize(self):
        return self._visualize_

    # ------- properties ---------------------------------------------------------------------------
    @property
    def assembled(self):
        return self._A_.assembled, self._b_.assembled

    @property
    def A(self):
        """A of Ax=b."""
        return self._A_

    @property
    def b(self):
        """b of Ax=b."""
        return self._b_

    @property
    def GMr(self):
        return self._GMr_

    @property
    def GMc(self):
        return self._GMc_

    @property
    def block_shape(self):
        """Return (m,n). The A is m by n block-wise."""
        return self._A_.bmat_shape

    @property
    def GLOBAL_shape(self):
        """Return (P, Q): A.assembled will be a sparse matrix of shape = (P, Q)."""
        # noinspection PyUnresolvedReferences
        return self.GMr.global_num_dofs, self.GMc.global_num_dofs
