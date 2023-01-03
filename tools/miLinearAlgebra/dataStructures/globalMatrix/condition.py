# -*- coding: utf-8 -*-
from numpy import linalg as nplinalg
from components.freeze.main import FrozenOnly
from root.config.main import RANK, MASTER_RANK, COMM


class ___GM_CONDITION___(FrozenOnly):
    """Condition of course is global condition."""
    def __init__(self, gm):
        self._gm_ = gm
        self._freeze_self_()

    @property
    def eig(self):
        """
        The eigenvalues and eigenvectors of (global) A.

        :return: A tuple of two outputs (Only in master core):

            1. (numpy.ndarray) w -- eigen values.
            2. (numpy.ndarray) v -- eigen vectors: v[:,i] is the eigenvector corresponding to
                the eigenvalue w[i].
        """
        M = self._gm_.do.gather_M_to_core(core=MASTER_RANK)  # does not clear the local value.
        if RANK == MASTER_RANK:
            M = M.toarray()
            w, v = nplinalg.eig(M)
        else:
            w, v = None, None
        w, v = COMM.bcast([w, v], root=MASTER_RANK)
        return w, v

    @property
    def condition_number(self):
        M = self._gm_.do.gather_M_to_core(core=MASTER_RANK)  # does not clear the local value.
        if RANK == MASTER_RANK:
            M = M.toarray()
            cn = nplinalg.cond(M)
        else:
            cn = None
        cn = COMM.bcast(cn, root=MASTER_RANK)
        return cn

    @property
    def sparsity(self):
        """(float) The sparsity of (global) A."""
        M = self._gm_.do.gather_M_to_core(core=MASTER_RANK)
        if RANK == MASTER_RANK:
            nnz = M.nnz
            s = M.shape
            sparsity = nnz/(s[0]*s[1])
        else:
            sparsity = None
        return COMM.bcast(sparsity, root=MASTER_RANK)
