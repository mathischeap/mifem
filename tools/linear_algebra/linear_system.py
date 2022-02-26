"""
We make a class to represent a linear system.

"""

import sys
if './' not in sys.path: sys.path.append('./')


from screws.frozen import FrozenClass, FrozenOnly
from tools.linear_algebra.elementwise_cache import EWC_ColumnVector, EWC_SparseMatrix


class LinearSystem(FrozenClass):
    """"""
    def __init__(self, A, b):
        """Ax = b.

        `A` and `b` are given as (local) blocks, we will bmat and concatenate them into (local) EWC_SparseMatrix and
        EWC_ColumnVector for the local elements during initializing the linear system.

        :param A: A 2-d list or tuple of matrix blocks.
        :param b: A 1-d list or tuple of vector blocks.
        """
        assert A.__class__ is EWC_SparseMatrix, "A type wrong."
        assert b.__class__ is EWC_ColumnVector, "b type wrong."

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
        self._local_distribution_ = [self._GMr_.local_dofs_distribution,
                                     self._GMc_.local_dofs_distribution,]

        assert A.bmat_shape[0] == len(self._local_distribution_[0]), "bmat_shape dis-match GMr shape."
        assert A.bmat_shape[1] == len(self._local_distribution_[1]), "bmat_shape dis-match GMc shape."

        self._customize_ = ___LinearSystem_Customize___(self)

        self._freeze_self_()


    def assemble_and_solve(self, solver_name, **solver_kwargs):
        """
        Assemble and solve self.

        We first assemble self into a global system and then solve the global system.

        :param solver_name:
        :param solver_kwargs:
        :return: a tuple of 5 outputs:
            1. (DistributedVector or DistributedVector) results -- The result vector.
            2. (int) info -- The info which provides convergence information:

                * 0 : successful exit
                * >0 : convergence to tolerance not achieved, number of iterations
                * -1 : divergence

            3. (float) beta -- The residual.
            4. (int) ITER -- The number of used iterations.
            5. (str) message

        """

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
        then, locally, E is of shape (29, 29)

        """
        return self._local_distribution_

    @property
    def block_shape(self):
        """Return (m,n). The A is m by n block-wise."""
        return self._A_.bmat_shape

    @property
    def local_shape(self):
        """Return (a, b): A[i] is a sparse matrix of shape = (a, b)."""
        return self._A_.shape[1:]

    @property
    def GLOBAL_shape(self):
        """Return (P, Q): A.assembled is a sparse matrix of shape = (P, Q)."""
        # noinspection PyUnresolvedReferences
        return self.GMr.GLOBAL_num_dofs, self.GMc.GLOBAL_num_dofs

    @property
    def customize(self):
        """Adjust Ax=b by adjusting EWC matrix A and EWC vector b."""
        return self._customize_




class ___LinearSystem_Customize___(FrozenOnly):
    """Used to define customizations to A and b simultaneously."""
    def __init__(self, ls):
        self._LS_ = ls
        self._freeze_self_()


    def apply_strong_BC(self, i, j, pd, pc=None, interpreted_as='local_dofs'):
        """

        :param i: We apply it to block [i][j].
        :param j: We apply it to block [i][j].
        :param pd: partial degrees of freedom
        :param pc: partial cochain of freedom
        :param interpreted_as: how we interpret the `pd` and `pc`.
        :return:
        """
        if i == j:
            assert pc is None, f"when pc is None, we must have i==j, " \
                               f"now i={i}, j={j}."
        if pc is None:
            assert i == j, \
                f"when do not provide pc, we must set diagonal block, " \
                f"so i == j, now i={i}, j={j}."
            assert pd.__class__.__name__ == 'PartialCochain', \
                "I need a PartialCochain when pc is None."
            pc = pd

        assert pc.__class__.__name__ == 'PartialCochain', f"pc must be a PartialCochain."

        I, J = self._LS_.block_shape
        assert i % 1 == 0, f"i={i}({i.__class__.__name__}) cannot be an index!"
        assert j % 1 == 0, f"j={j}({j.__class__.__name__}) cannot be an index!"
        assert 0 <= i < I and 0 <= j < J, f"(i,j)= ({i},{j}) is out of range!"

        if i == j:
            self._LS_.A.customize.\
                identify_global_rows_according_to_CSCG_partial_dofs(
                i, pd, interpreted_as=interpreted_as)
            self._LS_.b.customize.\
                set_entries_according_to_two_CSCG_partial_cochains(
                i, pd, interpreted_as=interpreted_as)
        else:
            self._LS_.A.customize.\
                off_diagonally_identify_rows_according_to_two_CSCG_partial_dofs(
                i, j, pd, pc, interpreted_as=interpreted_as)
            self._LS_.b.customize.\
                set_entries_according_to_two_CSCG_partial_cochains(
                i, pd, pc=pc, interpreted_as=interpreted_as)





if __name__ == '__main__':
    # mpiexec -n 5 python TOOLS\linear_algebra\linear_system.py

    from _3dCSCG.main import FormCaller
    from _3dCSCG.tests.random_objects import random_3D_mesh_and_space_of_total_load_around
    from tools.linear_algebra.functions import bmat, concatenate

    mesh, space = random_3D_mesh_and_space_of_total_load_around(500, exclude_periodic=True, mesh_boundary_num='>=2')
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

    print(Axb.GLOBAL_shape)