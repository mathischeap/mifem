# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly


class _2dCSCG_Mesh_Boundaries_Matplot(FrozenOnly):
    """"""
    def __init__(self, boundaries):
        """"""
        self._boundaries_ = boundaries
        self._freeze_self_()
