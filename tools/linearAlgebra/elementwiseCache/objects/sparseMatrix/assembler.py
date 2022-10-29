# -*- coding: utf-8 -*-
"""The assembler for the EWC sparse matrix."""
import numpy as np
from scipy.sparse import csr_matrix, csc_matrix
from numpy import diff
from screws.freeze.base import FrozenOnly
from tools.linearAlgebra.dataStructures.global_matrix.main import GlobalMatrix
from root.config.main import ASSEMBLE_COST, COMM, MPI, RANK, MASTER_RANK


class EWC_SparseMatrix_Assembler(FrozenOnly):
    """"""
    def __init__(self, MAT):
        self._MAT_ = MAT
        self._chain_method_ = 'silly'
        self._format_ = 'csc' # by default, we use csc format.
        self.Reset_cache()
        self._freeze_self_()

    def Reset_cache(self):
        self._cache_ = None
        self.___AMC___ = None # assembled matrix cache

    @property
    def chain_method(self):
        return self._chain_method_

    @chain_method.setter
    def chain_method(self, chm):
        """"""
        self._chain_method_ = chm

    def __call__(self, routine=None, **kwargs):
        """Do the assembling."""

        if self.___AMC___ is not None:
            return self.___AMC___

        COMM.barrier()
        twS = MPI.Wtime()

        if routine is None:
            A = self.___default_routine___()
        else:
            raise Exception(f"Assembling routine = {routine} is wrong!")

        ___locker___ = self._MAT_.do.___locker___

        if ___locker___:
            self.___AMC___ = A
        else:
            pass

        COMM.barrier()
        twE = MPI.Wtime()

        if RANK == MASTER_RANK:
            ASSEMBLE_COST['recent'].append(twE-twS)

        return A

    @property
    def format(self):
        return self._format_

    @format.setter
    def format(self, fmt):
        assert fmt in ('csr', 'csc'), f"cannot set format to {fmt}."
        self._format_ = fmt

    def ___default_routine___(self):
        """The default routine of assembling a EWC sparse matrix.

        :return:
        :rtype GlobalMatrix:
        """
        GMs = self._MAT_.gathering_matrices

        assert GMs[0] is not None, "I have no gathering matrix"
        assert GMs[1] is not None, "I have no gathering matrix"
        GI = GMs[0]
        GJ = GMs[1]
        DEP = int(GI.GLOBAL_num_dofs)
        WID = int(GJ.GLOBAL_num_dofs)

        if self.format == 'csc':
            SPA_MATRIX = csc_matrix
        elif self.format == 'csr':
            SPA_MATRIX = csr_matrix
        else:
            raise Exception

        # --- NO Cache: no sparsity locker, sparsity may change, so do not cache!---------
        if self._MAT_.do.___sparsity_locker___ is False:
            ROW = list()
            COL = list()
            DAT = list()

            A = SPA_MATRIX((DEP, WID)) # initialize a sparse matrix

            for i in self._MAT_: # go through all local sparse matrices.
                Mi = self._MAT_[i] # get the local sparse matrix
                indices = Mi.indices
                indptr = Mi.indptr
                data = Mi.data
                nums = diff(indptr)
                if Mi.__class__.__name__ == 'csc_matrix':
                    for j, num in enumerate(nums):
                        idx = indices[indptr[j]:indptr[j+1]]
                        ROW.extend(GI[i][idx])
                        COL.extend([GJ[i][j],]*num)
                elif Mi.__class__.__name__ == 'csr_matrix':
                    for j, num in enumerate(nums):
                        idx = indices[indptr[j]:indptr[j+1]]
                        ROW.extend([GI[i][j],]*num)
                        COL.extend(GJ[i][idx])
                else:
                    raise Exception("I can not handle %r."%Mi)
                DAT.extend(data)

                if len(DAT) > 1e7: # every 10 million data, we make it into sparse matrix.
                    _ = SPA_MATRIX((DAT, (ROW, COL)), shape=(DEP, WID)) # we make it into sparse

                    del ROW, COL, DAT
                    A += _
                    del _
                    ROW = list()
                    COL = list()
                    DAT = list()

            _ = SPA_MATRIX((DAT, (ROW, COL)), shape=(DEP, WID))  # we make it into sparse

            del ROW, COL, DAT
            A += _

            return GlobalMatrix(A)

        # Cache on: sparsity don't change, entries could still change ------------------------------
        else:
            if self._cache_ is None or self._cache_[0] != 'default':
                ROW = list()
                COL = list()
                DAT = list()

                for i in self._MAT_: # go through all local sparse matrices.
                    Mi = self._MAT_[i] # get the local sparse matrix
                    indices = Mi.indices
                    indptr = Mi.indptr
                    data = Mi.data
                    nums = diff(indptr)
                    if Mi.__class__.__name__ == 'csc_matrix':
                        for j, num in enumerate(nums):
                            idx = indices[indptr[j]:indptr[j+1]]
                            ROW.extend(GI[i][idx])
                            COL.extend([GJ[i][j],]*num)
                    elif Mi.__class__.__name__ == 'csr_matrix':
                        for j, num in enumerate(nums):
                            idx = indices[indptr[j]:indptr[j+1]]
                            ROW.extend([GI[i][j],]*num)
                            COL.extend(GJ[i][idx])
                    else:
                        raise Exception("I can not handle %r."%Mi)
                    DAT.extend(data)

                self._cache_ = ('default', ROW, COL)

            else:
                ROW, COL = self._cache_[1], self._cache_[2]
                DAT = list()

                if len(self._MAT_) > 0:
                    for i in self._MAT_: # go through all local sparse matrices.
                        Mi = self._MAT_[i] # get the local sparse matrix
                        DAT.append(Mi.data)
                    DAT = np.concatenate(DAT)
                else:
                    pass

            DAT = SPA_MATRIX((DAT, (ROW, COL)), shape=(DEP, WID))

            return GlobalMatrix(DAT)