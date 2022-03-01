

from screws.frozen import FrozenOnly
from scipy import sparse as spspa
from root.config import rAnk, mAster_rank, cOmm, np, sEcretary_rank, MPI


class DistributedVector(FrozenOnly):
    """
    An entry cannot be stored in multiple cores.

    This is convenient for, for example, storing the cochain of forms and the result vector of linear system.

    We will not adjust a DistributedVector (we will not do like clear values, change values and so on, but
    we may gather it to master and some other things which do not change the data). So we need it to be
    csc_matrix.
    """
    def __init__(self, V):
        if V.__class__.__name__ == 'csr_matrix': V = V.tocsc()
        if V.__class__.__name__ == 'ndarray':
            if np.ndim(V) == 1:
                V = spspa.csr_matrix(V).T
            elif np.ndim(V) == 2:
                assert V.shape[1] == 1, f"Need a shape = (n,1) ndarray, now it is {V.shape}."
                V = spspa.csc_matrix(V)
            else:
                raise Exception(f"Only accept 1- or 2-d array, now it is {np.ndim(V)}.")
        elif V is None:
            assert rAnk != mAster_rank, "in master core, can not give None."
        else:
            pass

        if rAnk == mAster_rank:
            assert spspa.isspmatrix_csc(V), "I need a scipy csc sparse matrix"
            shape = V.shape
        else:
            shape = None
        shape = cOmm.bcast(shape, root=mAster_rank)
        if rAnk != mAster_rank:
            if V is None:
                V = spspa.csc_matrix(shape)
            else:
                pass

        assert spspa.isspmatrix_csc(V) and V.shape[1] == 1, "V must be a csc_matrix of shape (x, 1)."

        self._V_ = V
        SHAPE = cOmm.gather(self.shape, root=sEcretary_rank)
        if rAnk == sEcretary_rank:
            for i, sp in enumerate(SHAPE):
                assert sp == SHAPE[0], f"shape in core {i} is different from shape in core 0."

        # check if it is a distributed vector ------------- BELOW ----------------------------------
        indices = V.indices
        indices = cOmm.gather(indices, root=mAster_rank)
        if rAnk == mAster_rank:
            measure = np.zeros(self.shape[0], dtype=int)
            for ind in indices:
                measure[ind] += 1

            assert np.max(measure) <= 1, f"vector is not distributed. {measure}"
        #================================================= ABOVE ===================================

        self._freeze_self_()

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

    def DO_distribute_to(self, *args, method='sequence'):
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
                GLOBAL_num_dofs = form.GLOBAL_num_dofs
                form.cochain.globe = DistributedVector(V[indices:indices+GLOBAL_num_dofs, 0])
                indices += GLOBAL_num_dofs
        else:
            raise NotImplementedError(f'distribution method: {method} not coded.')

    def ___PRIVATE_gather_V_to_core___(self, core=None, clean_local=False):
        """
        Gather all data to one core.

        :param core:
        :param bool clean_local: If True, we clear the local V while gathering.
        :return: A 1-d numpy array.
        """
        if core is None: core = mAster_rank
        if self.IS_master_dominating:
            if rAnk == core:
                return self.V.toarray()[:,0]
            else:
                return None
        else:
            GV = cOmm.gather(self.V, root=core)
            if clean_local: self._V_ = None
            V = None
            if rAnk == core:
                for Vi in GV:
                    if V is None:
                        V = Vi.toarray()[:,0]
                    else:
                        Vi = Vi.tocsc() # it should already be csc, but still, we do this to make sure.
                        indices = Vi.indices
                        data = Vi.data
                        V[indices] = data
            return V

    @property
    def IS_globally_empty(self):
        local_judge = True if self.nnz == 0 else False
        return cOmm.allreduce(local_judge, op=MPI.LAND)

    @property
    def IS_master_dominating(self):
        """(bool) return True if all data are in master core, empty in slave cores."""
        if rAnk == mAster_rank:
            ToF = True
        else:
            if self.V is None:
                ToF = True
            else:
                nnz = self.nnz
                ToF = nnz == 0
        return cOmm.allreduce(ToF, op=MPI.LAND)


