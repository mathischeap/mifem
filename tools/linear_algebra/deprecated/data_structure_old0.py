# -*- coding: utf-8 -*-
from screws.frozen import FrozenOnly
from scipy import sparse as spspa
from root.config import *
import matplotlib.pyplot as plt
from numpy import linalg as nplinalg


class DistributedVector(FrozenOnly):
    """
    An entry can be stored in multiple cores, but they must be all the same.

    This is convenient for, for example, storing the cochain of forms and the result vector of linear system
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
            assert spspa.issparse(V), "I need a scipy csc sparse matrix"
            shape = V.shape
        else:
            shape = None
        shape = cOmm.bcast(shape, root=mAster_rank)
        if rAnk != mAster_rank:
            if V is None:
                V = spspa.csc_matrix(shape)
            else:
                pass

        assert spspa.issparse(V), "I need a scipy csc sparse matrix."
        assert V.shape[1] == 1
        self._V_ = V
        SHAPE = cOmm.gather(self.shape, root=sEcretary_rank)
        if rAnk == sEcretary_rank:
            for i, sp in enumerate(SHAPE):
                assert sp == SHAPE[0], f"shape in core {i} is different from shape in core 0."

        indices = V.indices
        indices = cOmm.gather(indices, root=mAster_rank)
        if rAnk == mAster_rank:
            measure = np.zeros(self.shape[0], dtype=int)
            for ind in indices:
                measure[ind] += 1

            assert np.max(measure) <= 1, f"vector is not distributed. {measure}"

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
        :return:
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

class GlobalVector(FrozenOnly):
    """
    An entry can be split into parts and stored in multiple cores.

    This is convenient for, for example, the rhs of a linear system. To see the exact value of one entry,
    we must sum up that entry in all cores.
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
            assert spspa.issparse(V), "I need a scipy csc sparse matrix"
            shape = V.shape
        else:
            shape = None
        shape = cOmm.bcast(shape, root=mAster_rank)
        if rAnk != mAster_rank:
            if V is None:
                V = spspa.csc_matrix(shape)
            else:
                pass

        assert spspa.issparse(V), "I need a scipy sparse matrix"
        assert V.shape[1] == 1
        self._V_ = V
        SHAPE = cOmm.gather(self.shape, root=sEcretary_rank)
        if rAnk == sEcretary_rank:
            for i, sp in enumerate(SHAPE):
                assert sp == SHAPE[0], f"shape in core {i} is different from shape in core 0."

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




def ___split_A___(indptr, spliting_factor, PS):
    """
    A private function to help spliting global matrix into parts.

    :param indptr:
    :param spliting_factor:
    :param PS:
    :return:
    """
    assert len(indptr) == PS+1
    ST = 0
    __ = spliting_factor
    SL = [0, ]
    for i, ind in enumerate(indptr):
        if ind >= __:
            ST += 1
            SL.append(i)
            __ += spliting_factor
    if SL[-1] != PS :
        ST += 1
        SL.append(PS)
    assert ST >= 1 and ST == len(SL)-1
    return ST, SL

# PS = 3
# indptr = [0,0,0,0]
# spliting_factor = 7
# print(___split_A___(indptr, spliting_factor, PS))





class GlobalMatrix(FrozenOnly):
    """
    A wrapper of sparse matrix to adapt it to this library.

    :param M:
    :type M:
    """
    def __init__(self, M):
        if isinstance(M, (tuple, list)):
            assert len(M) == 2
            assert M[0]%1 == 0 and M[1]%1 == 0 and M[0] > 0 and M[1] > 0, \
                f"to initialize an empty GlobalMatrix of shape {M} is illegal."
            # we initialize an empty lil sparse matrix when provide tuple or list of 2 integers.
            M = spspa.lil_matrix(M)
        assert spspa.issparse(M), "I need a scipy sparse matrix"
        self._M_ = M
        self._DO_ = ___GM_DO___(self)
        self._visualize_ = ___GM_VISUALIZE___(self)
        self._condition_ = ___GM_CONDITION___(self)
        self._undirected_graph_ = ___Undirected_Graph___(self)
        self._directed_graph_ = ___Directed_Graph___(self)
        self.IS_regularly_distributed = False
        # Default be False, may not the case. When matters, ``do.claim_distribution_pattern`` first.
        SHAPE = cOmm.gather(self.shape, root=sEcretary_rank)
        if rAnk == sEcretary_rank:
            for i, sp in enumerate(SHAPE):
                assert sp == SHAPE[0], f"shape in core {i} is different from shape in core 0."
        self._freeze_self_()

    @property
    def M(self):
        return self._M_

    @property
    def shape(self):
        return self.M.shape

    @property
    def nnz(self):
        return self.M.nnz

    @property
    def GLOBAL_approx_nnz(self):
        local_nnz = self.nnz
        return cOmm.allreduce(local_nnz, op=MPI.SUM)

    @property
    def mtype(self):
        return self.M.__class__.__name__[:3]

    @property
    def nonempty_columns(self):
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
    def IS_globally_empty(self):
        local_judge = True if self.nnz == 0 else False
        return cOmm.allreduce(local_judge, op=MPI.LAND)

    @property
    def IS_regularly_distributed(self):
        if self._IS_regularly_distributed_ is True:
            assert self.mtype in ('csr', 'csc'), \
                "M has to be csr or csc matrix when IS_regularly_distributed=True "
        elif self._IS_regularly_distributed_ == 'row':
            assert self.mtype == 'csr'
            #     pass
            # else:
            #     if saFe_mode:
            #         assert self.do.___PRIVATE_check_if_Iam_row_major___() == (True, 0), \
            #             "It is not a row major matrix."
        elif self._IS_regularly_distributed_ == 'column':
            assert self.mtype == 'csc'
            #     pass
            # else:
            #     if saFe_mode:
            #         assert self.do.___PRIVATE_check_if_Iam_column_major___() == (True, 0), \
            #             "It is not a column major matrix."
        else:
            pass
        return self._IS_regularly_distributed_

    @IS_regularly_distributed.setter
    def IS_regularly_distributed(self, IS_regularly_distributed):
        """
        Do not set this easily. If do so, make sure you have double checked it since we do not and the code
        does not neither (because the check is not fast at all).
        """
        assert IS_regularly_distributed in (True, 'row', 'column', False), \
            f"IS_regularly_distributed={IS_regularly_distributed} wrong, " \
            f"can only be one of (True, 'row', 'column', False)."
        self._IS_regularly_distributed_ = IS_regularly_distributed

    @property
    def IS_master_dominating(self):
        """(bool) return True if all data are in master core, empty in slave cores."""
        if rAnk == mAster_rank:
            ToF = True
        else:
            if self.M is None:
                ToF = True
            else:
                nnz = self.nnz
                ToF = nnz == 0
        return cOmm.allreduce(ToF, op=MPI.LAND)


    @property
    def IS_square(self):
        shape = self._M_.shape
        return shape[0] == shape[1]

    @property
    def T(self):
        GMT = GlobalMatrix(self.M.T)
        if self.IS_regularly_distributed == 'row':
            GMT.IS_regularly_distributed = 'column'
        elif self.IS_regularly_distributed == 'column':
            GMT.IS_regularly_distributed = 'row'
        else:
            pass
        return GMT

    def __neg__(self):
        GM = GlobalMatrix(-self.M)
        GM.IS_regularly_distributed = self.IS_regularly_distributed
        return GM

    def __mul__(self, other):
        """

        :param other:
        :return:
        """
        assert isinstance(other, (int, float))
        GM = GlobalMatrix(other * self.M)
        GM.IS_regularly_distributed = self.IS_regularly_distributed
        return GM

    def __rmul__(self, other):
        """

        :param other:
        :return:
        """
        assert isinstance(other, (int, float))
        GM = GlobalMatrix(other * self.M)
        GM.IS_regularly_distributed = self.IS_regularly_distributed
        return GM

    # noinspection PyMethodParameters
    def __matmul__(A, B):
        """
        @ operator; dot product: GlobalMatrix @ GlobalMatrix or GlobalMatrix @ GlobalVector.

        A @ B or A @ b.

        :param B:
        :return:

        .. seealso:: This function is tested with unittest named ``test_TOOLS_NO3_GlobalMatrix_matmul_test``
            in file ``TESTS.unittest_tools.py``.
        """
        if B.__class__.__name__ == 'GlobalMatrix':
            if A.IS_master_dominating and B.IS_master_dominating: # this is important.
                if rAnk == mAster_rank:
                    M = A.M @ B.M
                else:
                    M = spspa.csc_matrix((A.shape[0], B.shape[1]))
            elif A.IS_regularly_distributed in (True, 'row', 'column') or \
                B.IS_regularly_distributed in (True, 'row', 'column'):
                M = A.___PRIVATE_MATMUL_rowMajor_dot_columnMajor___(A, B)
            else:
                # PLEASE do NOT REACH HERE!
                s_m_r, s_m_c = A.DO.claim_distribution_pattern()
                if A.IS_regularly_distributed in ('row', 'column'):
                    pass
                else:
                    o_m_r, o_m_c = B.DO.claim_distribution_pattern()
                    if B.IS_regularly_distributed in ('row', 'column'):
                        pass
                    else:
                        SO_M_RC = [s_m_r, s_m_c, o_m_r, o_m_c]
                        MIN = min(SO_M_RC)
                        IND = SO_M_RC.index(MIN)
                        if IND == 0:
                            if MIN == 0:
                                A.DO.tocsr()
                                A.IS_regularly_distributed = 'row'
                            else:
                                A.DO.___PRIVATE_being_regularly_distributed___('row')
                        elif IND == 1:
                            if MIN == 0:
                                A.DO.tocsc()
                                A.IS_regularly_distributed = 'column'
                            else:
                                A.DO.___PRIVATE_being_regularly_distributed___('column')
                        elif IND == 2:
                            if MIN == 0:
                                B.DO.tocsr()
                                B.IS_regularly_distributed = 'row'
                            else:
                                B.DO.___PRIVATE_being_regularly_distributed___('row')
                        else:
                            if MIN == 0:
                                B.DO.tocsc()
                                B.IS_regularly_distributed = 'column'
                            else:
                                B.DO.___PRIVATE_being_regularly_distributed___('column')
                _ = A.IS_regularly_distributed # just to check self is OK.
                _ = B.IS_regularly_distributed # just to check other is OK.
                M = A.___PRIVATE_MATMUL_rowMajor_dot_columnMajor___(A, B)
            # ...
            GM = GlobalMatrix(M)
            GM.IS_regularly_distributed = True
            return GM
        elif B.__class__.__name__ == 'GlobalVector':
            # not fast but okay (only do it once), use the master to collect/distribute vector.
            v = B.___PRIVATE_gather_V_to_core___(mAster_rank)
            if A.IS_master_dominating:
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
            if A.IS_master_dominating:
                if rAnk == mAster_rank:
                    v = A.M @ v
                else:
                    v = spspa.csc_matrix((A.shape[0], 1))
            else:
                v = A.___PRIVATE_distribute_vector___(v)
                v = A.M @ v
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
    def undirected_graph(self):
        """For symmetric matrix."""
        return self._undirected_graph_

    @property
    def directed_graph(self):
        """For symmetric and unsymmetric matrix."""
        return self._directed_graph_

    @property
    def condition(self):
        return self._condition_

    @property
    def visualize(self):
        return self._visualize_




    @property
    def DO(self):
        return self._DO_

    def ___PRIVATE_row_major___(self):
        """
        Make the matrix row-major distributed: each entire row only exsits in one core; could not be
        shared.

        do NOT EASILY USE THIS! THIS IS VERY SLOW.

        :return:
        """
        self._M_ = self._M_.tocsr() # of course, first make it csr format.
        indptr = self._M_.indptr
        num_elements_per_row = np.diff(indptr)
        having_elements_rows = np.argwhere(num_elements_per_row>0)[:,0]
        elements_num_rows = num_elements_per_row[having_elements_rows]
        having_elements_rows = cOmm.gather(having_elements_rows, root=mAster_rank)
        elements_num_rows = cOmm.gather(elements_num_rows, root=mAster_rank)

        if rAnk == mAster_rank:
            data = list()
            PTR = [0,]
            for i, er in enumerate(having_elements_rows):
                # i: core number
                data.append(elements_num_rows[i])
                # noinspection PyUnresolvedReferences
                PTR.append(PTR[-1]+len(er))
            DT = spspa.csc_matrix((np.concatenate(data), np.concatenate(having_elements_rows), PTR),
                                  shape=(self.shape[0],sIze))
            del data, PTR, having_elements_rows, elements_num_rows
            in_core = spspa.csc_matrix.argmax(DT, axis=1)
            ROWS = list()
            for i in range(sIze):
                ROWS.append(np.argwhere(in_core==i)[:,0])
        else:
            ROWS = None

        ROWS = cOmm.bcast(ROWS, root=mAster_rank)

        self._M_ = self._M_.tolil()
        for i in range(sIze):
            if i == rAnk: # recv
                for j in range(sIze):
                    if j != i:
                        pass # recv
                        data_recv = cOmm.recv(source=j, tag=rAnk)
                        self._M_[ROWS[rAnk],:] += data_recv
            else: # send
                rows_2b_sent = ROWS[i]
                data_2b_sent = self._M_[rows_2b_sent, :]
                cOmm.send(data_2b_sent, dest=i, tag=i)
        for i, rows in enumerate(ROWS):
            if i != rAnk:
                self._M_[rows, :] = 0
        self._M_ = self._M_.tocsr()
        self.IS_regularly_distributed = 'row'

    def ___PRIVATE_column_major___(self):
        """
        Make the matrix column-major distributed: each entire column only exsits in one core; could not be
        shared.

        do NOT EASILY USE THIS! THIS IS VERY SLOW.

        :return:
        """
        self._M_ = self._M_.tocsc() # of course, first make it csc format.
        indptr = self._M_.indptr
        num_elements_per_column = np.diff(indptr)
        HER = np.argwhere(num_elements_per_column>0)[:,0]
        elements_num_columns = num_elements_per_column[HER]
        having_elements_columns = cOmm.gather(HER, root=mAster_rank)
        elements_num_columns = cOmm.gather(elements_num_columns, root=mAster_rank)

        if rAnk == mAster_rank:
            data = list()
            PTR = [0,]
            for i, er in enumerate(having_elements_columns):
                # i: core number
                data.append(elements_num_columns[i])
                # noinspection PyUnresolvedReferences
                PTR.append(PTR[-1]+len(er))
            DT = spspa.csr_matrix((np.concatenate(data), np.concatenate(having_elements_columns), PTR),
                                  shape=(sIze, self.shape[1]))
            del data, PTR, having_elements_columns, elements_num_columns
            in_core = spspa.csr_matrix.argmax(DT, axis=0)
            COLUMNS = list()
            for i in range(sIze):
                COLUMNS.append(np.argwhere(in_core==i)[:,1])
        else:
            COLUMNS = None

        COLUMNS = cOmm.bcast(COLUMNS, root=mAster_rank)

        self._M_ = self._M_.T.tolil()
        for i in range(sIze):
            if i == rAnk: # recv
                for j in range(sIze):
                    if j != i:
                        pass # recv
                        data_recv = cOmm.recv(source=j, tag=rAnk)
                        self._M_[COLUMNS[rAnk],:] += data_recv
            else: # send
                cols_2b_sent = COLUMNS[i]
                data_2b_sent = self._M_[cols_2b_sent,:]
                cOmm.send(data_2b_sent, dest=i, tag=i)

        for i, cols in enumerate(COLUMNS):
            if i != rAnk:
                self._M_[cols,:] = 0
        self._M_ = self._M_.tocsr().T
        self.IS_regularly_distributed = 'column'

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
        if A.IS_regularly_distributed: # True, row or cloumn
            for i in range(sIze):
                N = cOmm.bcast(BM, root=i)
                if C is None:
                    C = AM @ N
                else:
                    C += AM @ N
        elif B.IS_regularly_distributed: # True, row or cloumn
            for i in range(sIze):
                M = cOmm.bcast(AM, root=i)
                if C is None:
                    C = M @ BM
                else:
                    C += M @ BM
        else:
            raise Exception(f"At least one component must be regularly distributed")
        return C

    def ___PRIVATE_gather_M_to_core___(self, core=None, clean_local=False, splitting_factor=50000000):
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

        splitting_factor = int(splitting_factor)
        assert splitting_factor > 0, f"splitting_factor={splitting_factor} wrong, must be > 0."

        if clean_local:
            if self.GLOBAL_approx_nnz < splitting_factor:
                M = cOmm.gather(self.M, root=core)
                if rAnk == core:
                    self._M_ = np.sum(M)
                else:
                    self._M_ = spspa.csc_matrix(self.shape)
            else:
                assert core == mAster_rank, "This routine only work for root=master yet!"
                tree = tRee(2)
                for Hi in tree:  # combine A to master core, A will be cleaned.
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
                            self._M_ = spspa.csc_matrix(self.shape)
                        else:
                            for t in range(ST):
                                if self.mtype == 'csr':
                                    tbs = self.M[SL[t]:SL[t+1], :]
                                elif self.mtype == 'csc':
                                    tbs = self.M[:, SL[t]:SL[t+1]]
                                else:
                                    raise Exception()
                                cOmm.send(tbs, **Hi[1])
                            self._M_ = spspa.csc_matrix(self.shape)

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

        Not very fast (But OK since we do not use it very offen, maybe one or two times
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
        if what is None: what = all_can_know
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





class ___GM_DO___(FrozenOnly):
    def __init__(self, gm):
        self._gm_ = gm
        self._freeze_self_()

    # def clear_row(self, r):
    #     """
    #     Make row #i all zero.
    #
    #     This is done in all cores such that that row is for sure zero then.
    #
    #     :param r:
    #     :return: None
    #     """
    #     if self._gm_.mtype != 'lil':
    #         self._gm_._M_ = self._gm_._M_.tolil()
    #     self._gm_._M_[r, :] = 0
    #
    #
    # def clear_rows(self, rows):
    #     """
    #
    #     :param rows:
    #     :return:
    #     """
    #     if self._gm_.mtype != 'lil':
    #         self._gm_._M_ = self._gm_._M_.tolil()
    #     self._gm_._M_[rows, :] = 0
    #
    # def set_value(self, i, j, value):
    #     """
    #     Set M[i,j] = v only in master core, set M[i,j] = 0 in all other cores.
    #
    #     :param i:
    #     :param j:
    #     :param value:
    #     :return:
    #     """
    #     if rAnk == mAster_rank:
    #         self._gm_._M_[i, j] = value
    #     else:
    #         self._gm_._M_[i, j] = 0
    #
    #
    # def set_values(self, I, J, values):
    #     """
    #
    #     :param I:
    #     :param J:
    #     :param values:
    #     :return:
    #     """
    #     if self._gm_.mtype != 'lil':
    #         self._gm_._M_ = self._gm_._M_.tolil()
    #     if rAnk == mAster_rank:
    #         self._gm_._M_[I, J] = values
    #     else:
    #         self._gm_._M_[I, J] = 0
    #
    # def identify_rows(self, rows):
    #     """Make rows be diagonal 1.
    #
    #     :param rows:
    #     :return:
    #     """
    #     assert self._gm_.IS_square, "Can only do to square matrix."
    #
    #     if self._gm_.mtype != 'lil':
    #         self._gm_._M_ = self._gm_._M_.tolil()
    #
    #     self._gm_._M_[rows, :] = 0
    #
    #     if rAnk == mAster_rank:
    #         self._gm_._M_[rows, rows] = 1
    #     else:
    #         self._gm_._M_[rows, rows] = 0

    def tocsc(self):
        """"""
        if self._gm_.mtype == 'csc':
            pass
        else:
            self._gm_._M_ = self._gm_._M_.tocsc()

    def tocsr(self):
        """"""
        if self._gm_.mtype == 'csr':
            pass
        else:
            self._gm_._M_ = self._gm_._M_.tocsr()

    def claim_distribution_pattern(self):
        """
        We parse the structure of M and classify the global matirx into 'row', 'column' or False.

        :return:
        """
        IS_row_major, mixed_rows = self.___PRIVATE_check_if_Iam_row_major___()
        if IS_row_major:
            self._gm_.IS_regularly_distributed = 'row'
            return 'row', mixed_rows
        IS_column_major, mixed_cols = self.___PRIVATE_check_if_Iam_column_major___()
        if IS_column_major:
            self._gm_.IS_regularly_distributed = 'column'
            return 'column', mixed_cols
        self._gm_.IS_regularly_distributed = False
        return mixed_rows, mixed_cols

    def ___PRIVATE_being_regularly_distributed___(self, major):
        """
        do NOT EASILY USE THIS! THIS IS VERY SLOW.

        :param major: 'row' or 'column' Row major or column major. If row major: in each core,
            We store some sequential rows or the sparse matrix, otherwise, we store some
            sequential columns. No sharing anymore. So if Mij !=0, and Mij will only be non-zero
            in one core, and be zero in all other cores!
        :return:
        """
        # easy situations ...
        if major == 'row':
            if self._gm_.IS_regularly_distributed == 'row':
                return
            else:
                pass
        elif major == 'column':
            if self._gm_.IS_regularly_distributed == 'column':
                return
            else:
                pass
        else:
            raise Exception(f"sequential distribution major={major} wrong, "
                            f"should be 'row' or 'column'.")
        # difficulty situations ...
        if major == 'row':
            self._gm_.___PRIVATE_row_major___()
        elif major == 'column':
            self._gm_.___PRIVATE_column_major___()
        else:
            raise Exception(f"sequential distribution major={major} wrong, "
                            f"should be 'row' or 'column'.")

    def ___PRIVATE_check_if_Iam_row_major___(self):
        """
        Check if self._gm_ is row major distributed.

        :return: A tuple of two outputs:

            1. (bool) -- True of False
            2. (int) -- How many shared rows

        """
        if self._gm_.mtype != 'csr':
            self._gm_._M_ = self._gm_.M.tocsr()
        indptr = self._gm_.M.indptr
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
                mixed_rows = 0
            else:
                ToF = False
                mixed_rows = L_her - len(POOL)
        else:
            ToF, mixed_rows = None, None
        ToF = cOmm.bcast(ToF, root=mAster_rank)
        mixed_rows = cOmm.bcast(mixed_rows, root=mAster_rank)
        return ToF, mixed_rows

    def ___PRIVATE_check_if_Iam_column_major___(self):
        """
        Check if self._gm_ is column major distributed.

        :return: A tuple of two outputs:

            1. (bool) -- True of False
            2. (int) -- How many shared columns

        """
        if self._gm_.mtype != 'csc':
            self._gm_._M_ = self._gm_.M.tocsc()
        indptr = self._gm_.M.indptr
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
                mixed_columns = 0
            else:
                ToF = False
                mixed_columns = L_her - len(POOL)
        else:
            ToF, mixed_columns = None, None
        ToF = cOmm.bcast(ToF, root=mAster_rank)
        mixed_columns = cOmm.bcast(mixed_columns, root=mAster_rank)
        return ToF, mixed_columns


class ___GM_VISUALIZE___(FrozenOnly):
    def __init__(self, gm):
        self._gm_ = gm
        self._freeze_self_()

    def spy(self, markerfacecolor='k', markeredgecolor='g', markersize=6):
        """
        The spy plot of self.

        :param markerfacecolor:
        :param markeredgecolor:
        :param markersize:
        :return:
        """
        M = self._gm_.___PRIVATE_gather_M_to_core___(core=mAster_rank)
        if rAnk == mAster_rank:
            fig = plt.figure()
            plt.spy(M,
                    markerfacecolor=markerfacecolor, markeredgecolor=markeredgecolor,
                    markersize=markersize)
            plt.tick_params(axis='both', which='major', direction='out')
            plt.tick_params(which='both', top=True, right=True, labelbottom=True, labelright=True)
            plt.show()
            return fig


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
        M = self._gm_.___PRIVATE_gather_M_to_core___(core=mAster_rank)
        if rAnk == mAster_rank:
            M = M.toarray()
            w, v = nplinalg.eig(M)
        else:
            w, v = None, None
        return w, v

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







class ___Undirected_Graph___(FrozenOnly):
    """
    #TODO:
    """
    def __init__(self, gm):
        self._gm_ = gm
        self._freeze_self_()


class ___Directed_Graph___(FrozenOnly):
    """
    #TODO:
    """
    def __init__(self, gm):
        self._gm_ = gm
        self._freeze_self_()