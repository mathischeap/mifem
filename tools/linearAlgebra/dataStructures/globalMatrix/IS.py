# -*- coding: utf-8 -*-
from root.config.main import RANK, MASTER_RANK, COMM, MPI
from components.freeze.main import FrozenOnly



class ___GM_IS___(FrozenOnly):
    def __init__(self, gm):
        self._gm_ = gm
        self._regularly_distributed_ = False
        self._freeze_self_()

    @property
    def globally_empty(self):
        """I am an empty matrix in all cores."""
        local_judge = True if self._gm_.nnz == 0 else False
        return COMM.allreduce(local_judge, op=MPI.LAND)

    @property
    def regularly_distributed(self):
        if self._regularly_distributed_ is True:
            assert self._gm_.mtype in ('csr', 'csc'), \
                "M has to be csr or csc matrix when IS_regularly_distributed=True "
        elif self._regularly_distributed_ == 'row':
            assert self._gm_.mtype == 'csr'
            #     pass
            # else:
            #     if SAFE_MODE:
            #         assert self.do.___PRIVATE_check_if_Iam_row_major___() == (True, 0), \
            #             "It is not a row major matrix."
        elif self._regularly_distributed_ == 'column':
            assert self._gm_.mtype == 'csc'
            #     pass
            # else:
            #     if SAFE_MODE:
            #         assert self.do.___PRIVATE_check_if_Iam_column_major___() == (True, 0), \
            #             "It is not a column major matrix."
        else:
            pass
        return self._regularly_distributed_

    @regularly_distributed.setter
    def regularly_distributed(self, regularly_distributed):
        """
        Do not set this easily. If it does so, make sure you have checked it since we do not and the code
        does neither (because the check is not fast at all).
        """
        assert regularly_distributed in (True, 'row', 'column', False), \
            f"regularly_distributed={regularly_distributed} wrong, " \
            f"can only be one of (True, 'row', 'column', False)."
        self._regularly_distributed_ = regularly_distributed


    @property
    def master_dominating(self):
        """(bool) return True if all data are in master core, empty in slave cores."""
        if RANK == MASTER_RANK:
            ToF = True
        else:
            if self._gm_.M is None:
                ToF = True
            else:
                nnz = self._gm_.nnz
                ToF = nnz == 0
        return COMM.allreduce(ToF, op=MPI.LAND)

    @property
    def square(self):
        shape = self._gm_._M_.shape
        return shape[0] == shape[1]