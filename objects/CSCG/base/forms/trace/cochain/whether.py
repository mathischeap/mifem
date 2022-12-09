# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/30/2022 10:22 PM
"""
from components.freeze.main import FrozenOnly
from root.config.main import COMM, MPI

class CSCG_cT_Cochain_Whether(FrozenOnly):
    """"""

    def __init__(self, tf):
        """"""
        self._tf_ = tf
        self._full_ = None
        self._freeze_self_()

    @property
    def full(self):
        if self._tf_.cochain._local_ is None:
            raise Exception(f'no local cochain, make no sense to tell full or not.')
        return COMM.allreduce(self._full_, op=MPI.LAND)
