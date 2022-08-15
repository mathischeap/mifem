# -*- coding: utf-8 -*-
from screws.freeze.base import FrozenOnly


class GlobalVectorDo(FrozenOnly):
    """"""

    def __init__(self, GV):
        self._v_ = GV
        self._freeze_self_()

    def gather_V_to_core(self, **kwargs):
        return self._v_.___PRIVATE_gather_V_to_core___(**kwargs)

    def resemble_row_distribution_of(self, GM):
        return self._v_.___PRIVATE_resemble_row_distribution_of___(GM)