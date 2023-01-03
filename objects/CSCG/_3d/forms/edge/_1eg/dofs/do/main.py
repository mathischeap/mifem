# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly
from objects.CSCG._3d.forms.edge._1eg.dofs.do.find import _3dCSCG_E1F_Dofs_DoFind


class _3dCSCG_E1F_Dofs_Do(FrozenOnly):
    """"""
    def __init__(self, dofs):
        """"""
        self._dofs_ = dofs
        self._find_ = None
        self._freeze_self_()

    @property
    def find(self):
        if self._find_ is None:
            self._find_ = _3dCSCG_E1F_Dofs_DoFind(self._dofs_)
        return self._find_
