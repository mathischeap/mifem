# -*- coding: utf-8 -*-

from components.freeze.base import FrozenOnly


class _2dCSCG_Regions_Vis(FrozenOnly):
    """"""
    def __init__(self, regions):
        """"""
        self._regions_ = regions
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        """"""
