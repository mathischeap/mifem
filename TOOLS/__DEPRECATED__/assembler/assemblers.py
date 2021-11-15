# -*- coding: utf-8 -*-

from SCREWS.frozen import FrozenOnly
import numpy as np
from TOOLS.linear_algebra.data_structures import GlobalMatrix
from scipy import sparse as spspa


class AssemblerBase(FrozenOnly):
    def __init__(self, assembler):
        self._assembler_ = assembler
        self.___PRIVATE_check_G___()
        self._freeze_self_()

    def ___PRIVATE_check_M___(self, M):
        raise NotImplementedError('To be over-ridden in the children.')

    def ___PRIVATE_check_G___(self):
        raise NotImplementedError('To be over-ridden in the children.')







class _3dCSCG_EWC(AssemblerBase):
    def __init__(self, assembler):
        super(_3dCSCG_EWC, self).__init__(assembler)

    @property
    def default_method(self):
        return 'advanced_v0'

    def ___PRIVATE_check_G___(self):
        assert self._assembler_._Gi__.__class__.__name__ == 'Gathering_Matrix'
        assert self._assembler_._G_j_.__class__.__name__ == 'Gathering_Matrix'

    def ___PRIVATE_check_M___(self, M):
        assert self._assembler_._Gi__.GLOBAL_len == M.GLOBAL_len, "GLOBAL length do not match."
        assert self._assembler_._G_j_.GLOBAL_len == M.GLOBAL_len, "GLOBAL length do not match."

    def DO_assembling_van_advanced_v0(self, M):
        """"""
        GI = self._assembler_.Gi_
        GJ = self._assembler_.G_j
        DEP = GI.GLOBAL_num_dofs
        WID = GJ.GLOBAL_num_dofs
        ROW = list()
        COL = list()
        DAT = list()
        for i in M:
            Mi = M[i]
            indices = Mi.indices
            indptr = Mi.indptr
            data = Mi.data
            nums = np.diff(indptr)
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
        return GlobalMatrix(spspa.csc_matrix((DAT, (ROW, COL)), shape=(DEP, WID)))