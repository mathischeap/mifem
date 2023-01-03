# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly
from root.config.main import RANK, MASTER_RANK, COMM


class DistributedVectorDo(FrozenOnly):
    """"""

    def __init__(self, DV):
        self._v_ = DV
        self._freeze_self_()

    def distributed_to(self, *args, **kwargs):
        self._v_.___PRIVATE_be_distributed_to___(*args, **kwargs)

    def gather_V_to_core(self, core=None, clean_local=False):
        """
        Gather all data to one core.

        :param core:
        :param bool clean_local: If True, we clear the local V while gathering.
        :return: A 1-d numpy array.
        """
        if core is None:
            core = MASTER_RANK
        if self._v_.whether.master_dominating:
            if RANK == core:
                return self._v_.V.toarray()[:, 0]
            else:
                return None
        else:
            GV = COMM.gather(self._v_.V, root=core)
            if clean_local:
                self._v_._V_ = None
            V = None
            if RANK == core:
                for Vi in GV:
                    if V is None:
                        V = Vi.toarray()[:, 0]
                    else:
                        Vi = Vi.tocsc()  # it should already be csc, but still, we do this to make sure.
                        indices = Vi.indices
                        data = Vi.data
                        V[indices] = data
            return V
