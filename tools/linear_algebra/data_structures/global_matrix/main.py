# -*- coding: utf-8 -*-



import sys
if './' not in sys.path: sys.path.append('../')


from screws.freeze.main import FrozenOnly
from scipy import sparse as spspa
from root.config.main import *

from tools.linear_algebra.data_structures.global_matrix.do import ___GM_DO___
from tools.linear_algebra.data_structures.global_matrix.IS import ___GM_IS___
from tools.linear_algebra.data_structures.global_matrix.visualize import ___GM_VISUALIZE___
from tools.linear_algebra.data_structures.global_matrix.directed_graph import ___GM_Directed_Graph___
from tools.linear_algebra.data_structures.global_matrix.undirected_graph import ___GM_Undirected_Graph___
from tools.linear_algebra.data_structures.global_matrix.adjust import ___GM_ADJUST___
from tools.linear_algebra.data_structures.global_matrix.condition import ___GM_CONDITION___
from tools.linear_algebra.data_structures.global_matrix.helpers.split_A import ___split_A___

from tools.linear_algebra.data_structures.vectors.GLOBAL.main import GlobalVector
from tools.linear_algebra.data_structures.vectors.locally_full.main import LocallyFullVector
from tools.linear_algebra.data_structures.vectors.distributed.main import DistributedVector


class GlobalMatrix(FrozenOnly):
    """A wrapper of sparse matrix to adapt it to this library.

    :param M:
    :type M:
    """
    def __init__(self, M):
        """

        :param M: We will make the GlobalMatrix from this M.
            - scipy sparse matrix
            - a list or tuple of two positive integers: We will make an empty sparse lil_matrix
                from it.
        """
        #---------- parse input M -----------------------------------------------------------
        if isinstance(M, (tuple, list)): # a list or tuple of two positive integers
            assert len(M) == 2
            assert M[0]%1 == 0 and M[1]%1 == 0 and M[0] > 0 and M[1] > 0, \
                f"to initialize an empty GlobalMatrix of shape {M} is illegal."
            # we initialize an empty lil sparse matrix when provide tuple or list of 2 integers.
            M = spspa.lil_matrix(M)
        else:
            pass
        #------------- check M ---------------------------------------------------------------
        assert spspa.issparse(M), "I need a scipy sparse matrix"
        self._M_ = M
        # Default be False, may not the case. When matters, ``do.claim_distribution_pattern()`` first.
        SHAPE = cOmm.gather(self.shape, root=sEcretary_rank)
        if rAnk == sEcretary_rank:
            for i, sp in enumerate(SHAPE):
                assert sp == SHAPE[0], f"shape in core {i} is different from shape in core 0."
        #-----------------------------------------------------------------------------------------
        self._visualize_ = ___GM_VISUALIZE___(self)
        self._condition_ = ___GM_CONDITION___(self)
        self._undirected_graph_ = ___GM_Undirected_Graph___(self)
        self._directed_graph_ = ___GM_Directed_Graph___(self)
        self._adjust_ = ___GM_ADJUST___(self)
        self._DO_ = ___GM_DO___(self)
        self._IS_ = ___GM_IS___(self)
        self._freeze_self_()

    def __repr__(self):
        return f"GlobalMatrix{self.shape}:{id(self)}"

    @property
    def M(self):
        return self._M_

    @property
    def shape(self):
        return self.M.shape

    @property
    def nnz(self):
        """Local nnz."""
        return self.M.nnz

    @property
    def GLOBAL_approx_nnz(self):
        """Global nnz (approximately): how many nnz in all cores."""
        local_nnz = self.nnz
        return cOmm.allreduce(local_nnz, op=MPI.SUM)

    @property
    def mtype(self):
        """Matrix type."""
        return self.M.__class__.__name__[:3]

    @property
    def nonempty_columns(self):
        """The columns have non-zero."""
        if self.mtype == 'csc':
            _ = np.diff(self.M.indptr) != 0  # the columns have non-zero values.
            hnzc = np.argwhere(_ == True).ravel()
        elif self.mtype == 'csr':
            hnzc = self.M.indices
            hnzc = sorted(list(set(hnzc)))
        else:
            raise Exception(f"method 'nonempty_columns' does not work for mtype={self.mtype}")
        return hnzc

    @property
    def nonempty_rows(self):
        """The rows have non-zero."""
        if self.mtype == 'csr':
            _ = np.diff(self.M.indptr) != 0  # the rows have non-zero values.
            hnzr = np.argwhere(_ == True).ravel()
        elif self.mtype == 'csc':
            hnzr = self.M.indices
            hnzr = sorted(list(set(hnzr)))
        else:
            raise Exception(f"method 'nonempty_rows' does not work for mtype={self.mtype}")
        return hnzr

    @property
    def shared_rows(self):
        """
        Which rows are shared with other cores?

        While shared with which cores? This property has no clue.
        """
        _shared_rows_ = self.___PRIVATE_parse_distributions___(what='shared_rows')
        return _shared_rows_

    @property
    def T(self):
        """Transpose."""
        GMT = GlobalMatrix(self.M.T)
        if self.IS.regularly_distributed == 'row':
            GMT.IS.regularly_distributed = 'column'
        elif self.IS.regularly_distributed == 'column':
            GMT.IS.regularly_distributed = 'row'
        else:
            pass
        return GMT

    def __neg__(self):
        """ - self """
        GM = GlobalMatrix(-self.M)
        GM.IS.regularly_distributed = self.IS.regularly_distributed
        return GM

    def __mul__(self, other):
        """
        self * other

        :param other:
        :return:
        """
        assert isinstance(other, (int, float))
        GM = GlobalMatrix(other * self.M)
        GM.IS.regularly_distributed = self.IS.regularly_distributed
        return GM

    def __rmul__(self, other):
        """
        other * self

        :param other:
        :return:
        """
        assert isinstance(other, (int, float))
        GM = GlobalMatrix(other * self.M)
        GM.IS.regularly_distributed = self.IS.regularly_distributed
        return GM

    # noinspection PyMethodParameters
    def __matmul__(A, B):
        """
        @ operator; dot product: GlobalMatrix @ GlobalMatrix or GlobalMatrix @ GlobalVector.

        A @ B or A @ b.

        :param B:
        :return:

        """
        if B.__class__.__name__ == 'GlobalMatrix':
            if A.IS.master_dominating and B.IS.master_dominating: # this is important.
                if rAnk == mAster_rank:
                    M = A.M @ B.M
                else:
                    M = spspa.csc_matrix((A.shape[0], B.shape[1]))
            elif A.IS.regularly_distributed in (True, 'row', 'column') or \
                B.IS.regularly_distributed in (True, 'row', 'column'):
                M = A.___PRIVATE_MATMUL_rowMajor_dot_columnMajor___(A, B)
            else:
                raise Exception(f'Neither A, B is regular, dot them will be very expensive.')
            # ...
            GM = GlobalMatrix(M)
            GM.IS.regularly_distributed = True
            return GM
        elif B.__class__.__name__ == 'GlobalVector':
            # not fast but okay (only do it once), use the master to collect/distribute vector.
            v = B.___PRIVATE_gather_V_to_core___(mAster_rank)
            if A.IS.master_dominating:
                if rAnk == mAster_rank:
                    v = A.M @ v
                else:
                    v = spspa.csc_matrix((A.shape[0], 1))
            else:
                v = A.___PRIVATE_distribute_vector___(v)
                v = A.M @ v
            return GlobalVector(v)
        elif B.__class__.__name__ == 'DistributedVector':
            # not fast but okay (only do it once), use the master to collect/distribute vector.
            v = B.___PRIVATE_gather_V_to_core___(mAster_rank)
            if A.IS.master_dominating:
                if rAnk == mAster_rank:
                    v = A.M @ v
                else:
                    v = spspa.csc_matrix((A.shape[0], 1))
            else:
                v = A.___PRIVATE_distribute_vector___(v)
                v = A.M @ v
            return GlobalVector(v)
        elif B.__class__.__name__ == 'LocallyFullVector':
            v = A.M @ B.V
            return GlobalVector(v)
        else:
            raise NotImplementedError(f"Can not do GlobalMatrix @ {B.__class__.__name__}.")

    def __sub__(self, other):
        """
        GlobalMatrix - GlobalMatrix.

        :param other:
        :return:
        """
        assert other.__class__.__name__ == 'GlobalMatrix', \
            f'Can not do GlobalMatrix - {other.__class__.__name__}'
        assert self.shape == other.shape, f"self shape {self.shape} != other.shape {other.shape}"
        GM = GlobalMatrix(self.M - other.M)
        return GM

    def __add__(self, other):
        """
        GlobalMatrix + GlobalMatrix.

        :param other:
        :return:
        """
        assert other.__class__.__name__ == 'GlobalMatrix', \
            f'Can not do GlobalMatrix + {other.__class__.__name__}'
        assert self.shape == other.shape, f"self shape {self.shape} != other.shape {other.shape}"
        GM = GlobalMatrix(self.M + other.M)
        return GM


    @property
    def IS(self):
        return self._IS_

    @property
    def adjust(self):
        return self._adjust_

    @property
    def do(self):
        return self._DO_

    @property
    def undirected_graph(self):
        """For symmetric matrix."""
        return self._undirected_graph_

    @property
    def directed_graph(self):
        """For symmetric and un-symmetric matrix."""
        return self._directed_graph_

    @property
    def condition(self):
        return self._condition_

    @property
    def visualize(self):
        return self._visualize_

    @staticmethod
    def ___PRIVATE_MATMUL_rowMajor_dot_columnMajor___(A, B):
        """

        :param A:
        :param B:
        :return:
        """
        # C = A @ B writen as: C = A @ B
        C = None
        AM = A.M
        BM = B.M
        if A.IS.regularly_distributed: # True, row or column
            for i in range(sIze):
                N = cOmm.bcast(BM, root=i)
                if C is None:
                    C = AM @ N
                else:
                    C += AM @ N
        elif B.IS.regularly_distributed: # True, row or column
            for i in range(sIze):
                M = cOmm.bcast(AM, root=i)
                if C is None:
                    C = M @ BM
                else:
                    C += M @ BM
        else:
            raise Exception(f"At least one component must be regularly distributed")
        return C

    def ___PRIVATE_gather_M_to_core___(self, core=None,
        clean_local=False, splitting_factor=50000000):
        """
        We gather M to one core, leave M to be of all zero elements in all other cores.

        Unlike the gathering method for vector, this gathering method will damage the original data
        if ``clean_local`` is True.

        Do not easily use this, this can be very slow and can damage the data or cause
        failure for MPI.

        :param core: The core that recv all matrix elements.
        :param bool clean_local:
        :param int splitting_factor: When the matrix is too big, the MPI may fail (can not pass
            object larger than 2GB). So we need to split it. The criterion is: If the number of
            non-zero values are larger than ``splitting_factor`` we do the splitting. Remember,
            splitting will only be used when we do ``clean_local`` after gathering.
        :return:
        """
        if core is None: core = mAster_rank

        if core == mAster_rank and self.IS.master_dominating:
            # In this case, we directly return.
            return self.M

        splitting_factor = int(splitting_factor)
        assert splitting_factor > 0, f"splitting_factor={splitting_factor} wrong, must be > 0."

        if clean_local:
            if self.GLOBAL_approx_nnz < splitting_factor:
                M = cOmm.gather(self.M, root=core)
                if rAnk == core:
                    self._M_ = np.sum(M)
                else:
                    self._M_ = spspa.csr_matrix(self.shape)
            else:
                assert core == mAster_rank, "This routine only work for root=master yet!"
                tree = tRee(2)
                for Hi in tree:  # combine `A` to master core, `A` will be cleaned.
                    if Hi is None:
                        pass

                    elif Hi[0] == 'send':
                        if self.nnz <= splitting_factor:
                            ST, SL = 1, None
                        else:
                            if self.mtype == 'csr':
                                PS = self.shape[0]
                            elif self.mtype == 'csc':
                                PS = self.shape[1]
                            else:
                                raise Exception()
                            ST, SL = ___split_A___(self.M.indptr, splitting_factor, PS)

                        cOmm.send(ST, **Hi[1])
                        if ST == 1:
                            cOmm.send(self.M, **Hi[1])
                            self._M_ = spspa.csr_matrix(self.shape)
                        else:
                            for t in range(ST):
                                if self.mtype == 'csr':
                                    tbs = self.M[SL[t]:SL[t+1], :]
                                elif self.mtype == 'csc':
                                    tbs = self.M[:, SL[t]:SL[t+1]]
                                else:
                                    raise Exception()
                                cOmm.send(tbs, **Hi[1])
                            self._M_ = spspa.csr_matrix(self.shape)

                    elif Hi[0] == 'recv':
                        RT = cOmm.recv(**Hi[1])
                        if RT == 1:
                            self._M_ += cOmm.recv(**Hi[1])
                        else:
                            RECV = list()
                            for t in range(RT):
                                RECV.append(cOmm.recv(**Hi[1]))
                            if self.mtype == 'csr':
                                self._M_ += spspa.vstack(RECV, format='csr')
                            elif self.mtype == 'csc':
                                self._M_ += spspa.hstack(RECV, format='csc')
                            else:
                                raise Exception()

                    else:
                        raise Exception()

            if rAnk == mAster_rank: self._M_.sum_duplicates()

            return self._M_

        else:
            M = cOmm.gather(self.M, root=core)
            if rAnk == core:
                return np.sum(M)
            else:
                return spspa.csc_matrix(self.shape)

    def ___PRIVATE_distribute_vector___(self, v):
        """
        Distribute vector to all cores.

        Not very fast (But OK since we do not use it very common, maybe one or two times
        per iteration likely), maybe we can update this scheme later.

        :param v:
        :return:
        """
        nec = self.nonempty_columns
        NEC = cOmm.gather(nec, root=mAster_rank)
        if rAnk == mAster_rank:
            V = list()
            for i, nnc in enumerate(NEC):
                if v.__class__.__name__ == 'csc_matrix':
                    VI = v[nnc, 0].T.toarray()[0]
                else:
                    VI = v[nnc]
                V.append(VI)
        else:
            V = None
        V = cOmm.scatter(V, root=mAster_rank)
        V = spspa.csc_matrix((V, nec, [0, len(nec)]), shape=(self.shape[1], 1))
        return V

    def ___PRIVATE_parse_distributions___(self, what=None):
        """
        We use this private method to do some work of parse the structure of the matrix.

        :param what:
        :return:
        """
        all_can_know = ('shared_rows',)

        if what is None:
            what = all_can_know
        else:
            if isinstance(what, str):
                what = [what,]
            assert isinstance(what, (tuple, list)), "what input must be list or tuple"
            for w in what:
                assert w in all_can_know, f"Want to know{w}? Sorry, I cannot deal with"

        RETURN = tuple()

        for wt in what:
            if wt == 'shared_rows':
                # will update attribute ``self._shared_rows_``.
                hnzr = self.nonempty_rows
                fa = np.array([False for _ in range(self.shape[0])])
                fa[hnzr] = True
                FA = cOmm.gather(fa, root=sEcretary_rank)
                if rAnk == sEcretary_rank:
                    shared_rows = np.count_nonzero(FA, axis=0)
                    shared_rows = np.argwhere(shared_rows>=2).ravel()
                else:
                    shared_rows = None
                shared_rows = cOmm.bcast(shared_rows, root=sEcretary_rank)
                sr = np.array([False for _ in range(self.shape[0])])
                sr[shared_rows] = True
                RETURN += (np.all([fa, sr], axis=0),)
            else:
                raise NotImplementedError(f"what={wt} is not understandable.")

        if len(RETURN) == 1:
            return RETURN[0]
        else:
            return  RETURN

    def ___PRIVATE_check_row_major___(self):
        """Return True if I am row major."""
        if self.mtype != 'csr': return False
        indptr = self.M.indptr
        num_elements_per_row = np.diff(indptr)
        having_elements_rows = np.argwhere(num_elements_per_row>0)[:,0]
        having_elements_rows = cOmm.gather(having_elements_rows, root=mAster_rank)
        if rAnk == mAster_rank:
            POOL = set()
            L_her = 0
            for her in having_elements_rows:
                POOL.update(her)
                L_her += len(her)
            if len(POOL) == L_her:
                ToF = True
            else:
                ToF = False
        else:
            ToF = None
        ToF = cOmm.bcast(ToF, root=mAster_rank)
        return ToF

    def ___PRIVATE_check_col_major___(self):
        """Return True if I am column major."""
        if self.mtype != 'csc': return False
        indptr = self.M.indptr
        num_elements_per_col = np.diff(indptr)
        having_elements_cols = np.argwhere(num_elements_per_col>0)[:,0]
        having_elements_cols = cOmm.gather(having_elements_cols, root=mAster_rank)
        if rAnk == mAster_rank:
            POOL = set()
            L_her = 0
            for her in having_elements_cols:
                POOL.update(her)
                L_her += len(her)
            if len(POOL) == L_her:
                ToF = True
            else:
                ToF = False
        else:
            ToF = None
        ToF = cOmm.bcast(ToF, root=mAster_rank)
        return ToF,

    def ___PRIVATE_self_regularity_checker___(self):
        """"""
        if self.IS.regularly_distributed == 'row':
            assert self.___PRIVATE_check_row_major___()
        elif self.IS.regularly_distributed == 'column':
            assert self.___PRIVATE_check_col_major___()
        elif self.IS.regularly_distributed:
            assert self.___PRIVATE_check_row_major___() or self.___PRIVATE_check_col_major___()
        else:
            pass










if __name__ == '__main__':
    # mpiexec -n 6 python tools/linear_algebra/gathering/chain_matrix.py
    GV = DistributedVector
    GM = LocallyFullVector