


from screws.frozen import FrozenOnly
from scipy import sparse as spspa
from root.config import rAnk, mAster_rank, cOmm, np, sEcretary_rank, MPI
from tools.linear_algebra.data_structures.vectors.GLOBAL.adjust import ___GV_ADJUST___


class GlobalVector(FrozenOnly):
    """
    An entry can be split into parts and stored in multiple cores.

    This is convenient for, for example, the rhs of a linear system. To see the exact value of one entry,
    we must sum up that entry in all cores.

    GlobalVector may have to be adjusted, so we do not ask it to be csc_matrix.
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
            assert spspa.issparse(V), "I need a scipy sparse matrix"
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
        self._adjust_ = ___GV_ADJUST___(self)

        self._freeze_self_()

    @property
    def V(self):
        return self._V_

    @property
    def shape(self):
        return self.V.shape

    def __len__(self):
        return self.shape[0]

    @property
    def nnz(self):
        return self.V.nnz

    def __neg__(self):
        return GlobalVector(-self.V)

    def __sub__(self, other):
        """

        :param other:
        :return:
        """
        return GlobalVector(self.V - other.V)

    def __add__(self, other):
        """

        :param other:
        :return:
        """
        return GlobalVector(self.V + other.V)

    def __mul__(self, other):
        return GlobalVector(other*self.V)

    def __rmul__(self, other):
        return GlobalVector(other*self.V)

    def ___PRIVATE_gather_V_to_core___(self, core=None, clean_local=False):
        """
        Gather all vector to one core such that in all other core we have zero vector.

        :param core:
        :param bool clean_local: If True, we clear the local V while gathering.
        :return: A 1d ndarray that contains the vector in only one core.
        """
        if core is None: core = mAster_rank
        v = self.V
        v = cOmm.gather(v, root=core)
        if clean_local: self._V_ = None
        if rAnk == core:
            # noinspection PyUnresolvedReferences
            v = np.sum(v).toarray()[:, 0]
        return v

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

    @property
    def IS_globally_empty(self):
        local_judge = True if self.nnz == 0 else False
        return cOmm.allreduce(local_judge, op=MPI.LAND)

    def DO_resemble_row_distribution_of(self, GM):
        """
        We let self's distribution resemble that of a GM.

        :param GM:
        :return:
        """
        assert GM.shape[0] == self.shape[0], f"shape[0] does not match."
        already_match = set(self.V.indices) <= set(GM.nonempty_rows)
        already_match = cOmm.allreduce(already_match, op=MPI.LAND)
        if already_match: # already match, just stop here.
            return
        else:
            raise NotImplementedError(f"Not coded yet!")

    @property
    def adjust(self):
        return self._adjust_