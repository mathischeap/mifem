# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly


class _3dCSCG_ADT_TF_IS(FrozenOnly):
    """"""
    def __init__(self, adt):
        self._adt_ = adt
        self._freeze_self_()

    @property
    def hybrid(self):
        return True

    @property
    def inner_oriented(self):
        return True if self._adt_._orientation_ == 'inner' else False

    @property
    def ADF(self):
        return self._adt_.___IS_ADF___