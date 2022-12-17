# -*- coding: utf-8 -*-

from components.freeze.base import FrozenOnly


class _3dCSCG_Regions_Topology(FrozenOnly):
    """"""
    def __init__(self, regions):
        self._regions_ = regions
        self._freeze_self_()
