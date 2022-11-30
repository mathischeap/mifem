# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly
from root.config.main import MASTER_RANK, RANK, COMM, np

class GlobalVectorDo(FrozenOnly):
    """"""

    def __init__(self, GV):
        self._v_ = GV
        self._freeze_self_()

    def gather_V_to_core(self, core=None, clean_local=False):
        """
        Gather all vector to one core such that in all other core we have zero vector.

        :param core:
        :param bool clean_local: If True, we clear the local V while gathering.
        :return: A 1d ndarray that contains the vector in only one core.
        """
        if core is None: core = MASTER_RANK
        v = self._v_.V
        v = COMM.gather(v, root=core)
        if clean_local: self._v_._V_ = None
        if RANK == core:
            # noinspection PyUnresolvedReferences
            v = np.sum(v).toarray()[:, 0]
        return v
