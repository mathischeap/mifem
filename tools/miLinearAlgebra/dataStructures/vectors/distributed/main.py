# -*- coding: utf-8 -*-

from components.freeze.main import FrozenOnly
from scipy import sparse as spspa
from root.config.main import RANK, MASTER_RANK, COMM, np, SECRETARY_RANK
from tools.miLinearAlgebra.dataStructures.vectors.distributed.do import DistributedVectorDo
from tools.miLinearAlgebra.dataStructures.vectors.distributed.whether import DistributedVectorWhether


class DistributedVector(FrozenOnly):
    """
    An entry cannot be stored in multiple cores.

    This is convenient for, for example, storing the cochain of forms and the result vector of linear system.

    We will not adjust a DistributedVector (we will not do like clear values, change values and so on, but
    we may gather it to master and some other things which do not change the data). So we need it to be
    csc_matrix.
    """
    def __init__(self, V):
        """

        :param V:
            - csc_matrix or csr_matrix of shape (x, 1)
            - 1d array
            - 2d array of shape (x, 1)
            - None (cannot be the master core): we will make an empty sparse vector then.

        """
        # ------- parse input ---------------------------------------------------------------
        if V.__class__.__name__ == 'ndarray':
            if np.ndim(V) == 1:
                V = spspa.csr_matrix(V).T
            elif np.ndim(V) == 2:
                assert V.shape[1] == 1, f"Need a shape = (n,1) ndarray, now it is {V.shape}."
                V = spspa.csc_matrix(V)
            else:
                raise Exception(f"Only accept 1- or 2-d array, now it is {np.ndim(V)}.")
        elif V is None:
            assert RANK != MASTER_RANK, "in master core, can not give None."
        elif V.__class__.__name__ == 'csr_matrix':
            V = V.tocsc()
            assert V.shape[1] == 1
        elif V.__class__.__name__ == 'csc_matrix':
            assert V.shape[1] == 1
        else:
            raise Exception()

        # --------------- check input ------------------------------------------------------
        if RANK == MASTER_RANK:
            assert spspa.isspmatrix_csc(V), "I need a scipy csc sparse matrix"
            shape = V.shape
        else:
            shape = None

        shape = COMM.bcast(shape, root=MASTER_RANK)
        if RANK != MASTER_RANK:
            if V is None:
                V = spspa.csc_matrix(shape)
            else:
                pass

        assert spspa.isspmatrix_csc(V) and V.shape[1] == 1, "V must be a csc_matrix of shape (x, 1)."

        # --------------------------------------------------------------------------
        self._V_ = V
        SHAPE = COMM.gather(self.shape, root=SECRETARY_RANK)
        if RANK == SECRETARY_RANK:
            for i, sp in enumerate(SHAPE):
                assert sp == SHAPE[0], f"shape in core {i} is different from shape in core 0."

        # check if it is a distributed vector ------------- BELOW ----------------------------------
        indices = V.indices
        indices = COMM.gather(indices, root=MASTER_RANK)
        if RANK == MASTER_RANK:
            measure = np.zeros(self.shape[0], dtype=int)
            for ind in indices:
                measure[ind] += 1

            assert np.max(measure) <= 1, f"vector is not distributed. {measure}"
        # ================================================ ABOVE ===================================
        self._do_ = None
        self._whether_ = None
        self._freeze_self_()

    def __repr__(self):
        return f"DistributedVector{self.shape}:{id(self)}"

    @property
    def V(self):
        return self._V_

    @property
    def shape(self):
        return self.V.shape

    @property
    def indices(self):
        return self.V.indices

    def __len__(self):
        return self.shape[0]

    @property
    def nnz(self):
        return self.V.nnz

    @property
    def do(self):
        if self._do_ is None:
            self._do_ = DistributedVectorDo(self)
        return self._do_

    @property
    def whether(self):
        if self._whether_ is None:
            self._whether_ = DistributedVectorWhether(self)
        return self._whether_

    def ___PRIVATE_be_distributed_to___(self, *args, method='sequence'):
        """
        Consider this vector represents cochains of some forms in sequence, we can distribute this
        vector to the forms.

        :param args: the forms.
        :param method: It can be one of:

            1. ``sequence`` -- We distribute the values in sequence to *args.

        :return:
        """
        if method == 'sequence':
            indices = 0
            V = self.V
            for form in args:
                GLOBAL_num_dofs = form.num.global_dofs
                form.cochain.globe = DistributedVector(V[indices:indices+GLOBAL_num_dofs, 0])
                indices += GLOBAL_num_dofs
        else:
            raise NotImplementedError(f'distribution method: {method} not coded.')
