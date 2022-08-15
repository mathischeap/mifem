# -*- coding: utf-8 -*-


from numpy import linalg as nplinalg
from screws.freeze.main import FrozenOnly
from root.config.main import rAnk, mAster_rank, cOmm


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
        M = self._gm_.___PRIVATE_gather_M_to_core___(core=mAster_rank) # does not clear the local value.
        if rAnk == mAster_rank:
            M = M.toarray()
            w, v = nplinalg.eig(M)
        else:
            w, v = None, None
        w, v = cOmm.bcast([w, v], root=mAster_rank)
        return w, v

    @property
    def condition_number(self):
        M = self._gm_.___PRIVATE_gather_M_to_core___(core=mAster_rank) # does not clear the local value.
        if rAnk == mAster_rank:
            M = M.toarray()
            cn = nplinalg.cond(M)
        else:
            cn = None
        cn = cOmm.bcast(cn, root=mAster_rank)
        return cn

    @property
    def sparsity(self):
        """(float) The sparsity of (global) A."""
        M = self._gm_.___PRIVATE_gather_M_to_core___(core=mAster_rank)
        if rAnk == mAster_rank:
            nnz = M.nnz
            s = M.shape
            sparsity = nnz/(s[0]*s[1])
        else:
            sparsity = None
        return cOmm.bcast(sparsity, root=mAster_rank)