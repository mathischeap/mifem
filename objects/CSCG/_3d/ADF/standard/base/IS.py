# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly







class _3dCSCG_ADF_SF_IS(FrozenOnly):
    """"""
    def __init__(self, adf):
        self._adf_ = adf
        self._freeze_self_()

    @property
    def inner_oriented(self):
        return True if self._adf_._orientation_ == 'inner' else False

    @property
    def volume_form(self):
        """"""
        return self._adf_.ndim == self._adf_.k

    @property
    def hybrid(self):
        """An AD standard form must be a hybrid form."""
        return True

    @property
    def ADF(self):
        return self._adf_.___IS_ADF___