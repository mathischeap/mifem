# -*- coding: utf-8 -*-
from screws.freeze.base import FrozenOnly
from root.config.main import rAnk, mAster_rank, MPI, cOmm

class DistributedVectorIS(FrozenOnly):
    """"""

    def __init__(self, GV):
        self._v_ = GV
        self._freeze_self_()

    @property
    def globally_empty(self):
        local_judge = True if self._v_.nnz == 0 else False
        return cOmm.allreduce(local_judge, op=MPI.LAND)

    @property
    def master_dominating(self):
        """(bool) return True if all data are in master core, empty in slave cores."""
        if rAnk == mAster_rank:
            ToF = True
        else:
            if self._v_.V is None:
                ToF = True
            else:
                nnz = self._v_.nnz
                ToF = nnz == 0
        return cOmm.allreduce(ToF, op=MPI.LAND)