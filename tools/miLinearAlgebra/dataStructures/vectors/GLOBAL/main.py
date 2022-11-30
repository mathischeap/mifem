# -*- coding: utf-8 -*-

from components.freeze.main import FrozenOnly
from scipy import sparse as spspa
from root.config.main import RANK, MASTER_RANK, COMM, np, SECRETARY_RANK
from tools.miLinearAlgebra.dataStructures.vectors.GLOBAL.adjust import ___GV_ADJUST___
from tools.miLinearAlgebra.dataStructures.vectors.GLOBAL.do import GlobalVectorDo
from tools.miLinearAlgebra.dataStructures.vectors.GLOBAL.whether import GlobalVectorWhether

class GlobalVector(FrozenOnly):
    """An entry can be split into parts and stored in multiple cores.

    This is convenient for, for example, the rhs of a linear system. To see the exact value of one entry,
    we must sum up that entry in all cores.

    GlobalVector may have to be adjusted, so we do not ask it to be csc_matrix.
    """
    def __init__(self, V):
        """

        :param V: it can be:
            - csc_matrix of shape (x, 1).
            - 1d array
            - 2d array of shape (x, 1)
            - int: we will make an empty sparse matrix.
            - None: (Cannot be in the master core) we will make an empty sparse matrix.

        """
        #------ parse input V -------------------------------------------
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
        elif V.__class__.__name__ in ('int', 'int32', 'int64'):
            V = spspa.csc_matrix((V, 1))
        else:
            pass
        #----------------- check V ---------------------------------------
        if V.__class__.__name__ == 'csr_matrix':
            V = V.tocsc()

        if RANK == MASTER_RANK:
            assert spspa.issparse(V), "I need a scipy sparse matrix"
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
        #--------------------------------------------------------------------

        self._V_ = V
        SHAPE = COMM.gather(self.shape, root=SECRETARY_RANK)
        if RANK == SECRETARY_RANK:
            for i, sp in enumerate(SHAPE):
                assert sp == SHAPE[0], f"shape in core {i} is different from shape in core 0."
        self._adjust_ = ___GV_ADJUST___(self)
        self._do_ = None
        self._whether_ = None
        self._freeze_self_()

    def __repr__(self):
        return f"GlobalVector:{id(self)}"

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

    @property
    def adjust(self):
        return self._adjust_

    @property
    def do(self):
        if self._do_ is None:
            self._do_ = GlobalVectorDo(self)
        return self._do_

    @property
    def whether(self):
        if self._whether_ is None:
            self._whether_ = GlobalVectorWhether(self)
        return self._whether_