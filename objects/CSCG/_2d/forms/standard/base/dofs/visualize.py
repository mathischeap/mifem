# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly


class _2dCSCG_SF_dofs_VIS(FrozenOnly):
    """"""
    def __init__(self, dofs):
        self._dofs_ = dofs
        self._freeze_self_()
