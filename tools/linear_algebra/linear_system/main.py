# -*- coding: utf-8 -*-
"""
We make a class to represent a linear system.

"""

import sys
if './' not in sys.path: sys.path.append('./')

from screws.freeze.main import FrozenClass
from tools.linear_algebra.elementwise_cache.objects.sparse_matrix.main import EWC_ColumnVector, EWC_SparseMatrix
from tools.linear_algebra.linear_system.customize import ___LinearSystem_Customize___
from tools.linear_algebra.linear_system.condition import ___LinearSystem_Condition___
from tools.linear_algebra.linear_system.solve.main import ___LinearSystem_Solve___



class LinearSystem(FrozenClass):
    """"""
    def __init__(self, A, b):
        """Ax = b.

        `A` and `b` are given as (local) blocks, we will bmat and concatenate them into (local) EWC_SparseMatrix and
        EWC_ColumnVector for the local elements during initializing the linear system.

        :param A: A 2-d list or tuple of matrix blocks.
        :param b: A 1-d list or tuple of vector blocks.
        """
        assert A.__class__ is EWC_SparseMatrix, "A type wrong, must be an EWC_SparseMatrix."
        assert b.__class__ is EWC_ColumnVector, "b type wrong, must be an EWC_ColumnVector."

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

        # noinspection PyUnresolvedReferences
        self._local_distribution_ = [ self._GMr_.local_dofs_distribution,
                                      self._GMc_.local_dofs_distribution ]

        assert A.bmat_shape[0] == len(self._local_distribution_[0]), "bmat_shape dis-match GMr shape."
        assert A.bmat_shape[1] == len(self._local_distribution_[1]), "bmat_shape dis-match GMc shape."

        self._customize_ = ___LinearSystem_Customize___(self)
        self._condition_ = ___LinearSystem_Condition___(self)
        self._solve_ = ___LinearSystem_Solve___(self)

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
    def local_distribution(self):
        """The local num of dofs of the variables compositing A.
        For example, if A = ([E, F], [C, D]), and
        self.local_distribution = [[29, 6], [29, 6]],
        then, locally, E is of shape (29, 29), F is of shape (19, 6) and so on for C and D.

        """
        return self._local_distribution_

    @property
    def block_shape(self):
        """Return (m,n). The A is m by n block-wise."""
        return self._A_.bmat_shape

    @property
    def local_shape(self):
        """Return (a, b): A[i] locally in each mesh element is a sparse matrix of shape = (a, b)."""
        return self._A_.shape[1:]

    @property
    def GLOBAL_shape(self):
        """Return (P, Q): A.assembled will be a sparse matrix of shape = (P, Q)."""
        # noinspection PyUnresolvedReferences
        return self.GMr.GLOBAL_num_dofs, self.GMc.GLOBAL_num_dofs











if __name__ == '__main__':
    # mpiexec -n 5 python TOOLS\linear_algebra\linear_system\main.py

    from objects.CSCG._3d.master import FormCaller
    from objects.CSCG._3d.__tests__.random_objects.form_caller import random_mesh_and_space_of_total_load_around
    from tools.linear_algebra.elementwise_cache.operators.concatenate.main import bmat, concatenate

    mesh, space = random_mesh_and_space_of_total_load_around(500, exclude_periodic=True, mesh_boundary_num='>=2')
    FC = FormCaller(mesh, space)

    f2 = FC('2-f', is_hybrid=True)
    f3 = FC('3-f', is_hybrid=True)
    t2 = FC('2-t')

    M2 = f2.matrices.mass
    E32 = f2.matrices.incidence
    T2 = t2.matrices.trace
    BMAT = [[M2, E32.T, T2.T],
            [E32, None, None],
            [T2, None, None]]
    A = bmat(BMAT)
    SHAPE = A.shape
    assert SHAPE[0] == len(mesh.elements)
    assert SHAPE[1] == f2.NUM_basis + f3.NUM_basis + t2.NUM_basis
    assert SHAPE[2] == f2.NUM_basis + f3.NUM_basis + t2.NUM_basis
    A.gathering_matrices = ([f2, f3, t2], [f2, f3, t2])

    b0 = EWC_ColumnVector(mesh, f2.NUM_basis)
    b1 = EWC_ColumnVector(mesh, f3.NUM_basis)
    b2 = EWC_ColumnVector(mesh, t2.NUM_basis)
    b = concatenate([b0, b1, b2])
    SHAPE = b.shape
    assert SHAPE[0] == len(mesh.elements)
    assert SHAPE[1] == f2.NUM_basis + f3.NUM_basis + t2.NUM_basis
    assert SHAPE[2] == 1
    b.gathering_matrix = [f2, f3, t2]

    Axb = LinearSystem(A, b)

    # print(Axb.GLOBAL_shape)