# -*- coding: utf-8 -*-
from screws.freeze.base import FrozenOnly

# import matplotlib.pyplot as plt


class _2dCSCG_ES_VIS_Matplot(FrozenOnly):
    """"""
    def __init__(self, es):
        self._es_ = es
        self._freeze_self_()