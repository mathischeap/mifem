# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly
from root.config.main import RANK, MASTER_RANK, MPI, COMM

class DistributedVectorWhether(FrozenOnly):
    """"""

    def __init__(self, GV):
        self._v_ = GV
        self._freeze_self_()

    @property
    def globally_empty(self):
        local_judge = True if self._v_.nnz == 0 else False
        return COMM.allreduce(local_judge, op=MPI.LAND)

    @property
    def master_dominating(self):
        """(bool) return True if all data are in master core, empty in slave cores."""
        if RANK == MASTER_RANK:
            ToF = True
        else:
            if self._v_.V is None:
                ToF = True
            else:
                nnz = self._v_.nnz
                ToF = nnz == 0
        return COMM.allreduce(ToF, op=MPI.LAND)