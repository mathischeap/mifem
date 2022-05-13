# -*- coding: utf-8 -*-

from screws.freeze.base import FrozenOnly




class _2dCSCG_Regions_Topology(FrozenOnly):
    """"""
    def __init__(self, regions):
        self._regions_ = regions
        self._freeze_self_()