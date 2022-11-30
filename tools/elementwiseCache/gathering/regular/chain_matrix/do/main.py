# -*- coding: utf-8 -*-

from components.freeze.main import FrozenOnly

from tools.elementwiseCache.gathering.regular.chain_matrix.do.find import ___Chain_Gathering_Matrix_FIND___


class ___Chain_Gathering_Matrix_DO___(FrozenOnly):
    """"""
    def __init__(self, CGM):
        self._CGM_ = CGM
        self._find_ = None
        self._freeze_self_()

    @property
    def find(self):
        if self._find_ is None:
            self._find_ = ___Chain_Gathering_Matrix_FIND___(self._CGM_)
        return self._find_