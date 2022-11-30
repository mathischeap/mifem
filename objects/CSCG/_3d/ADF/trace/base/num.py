# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly



class _3dCSCG_ADF_T_NUM(FrozenOnly):
    """"""
    def __init__(self, adt):
        self._adt_ = adt
        self._freeze_self_()

    @property
    def global_dofs(self):
        return self._adt_.prime.num.global_dofs

    @property
    def basis(self):
        return self._adt_.prime.num.basis