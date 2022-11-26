# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly
from objects.CSCG._3d.forms.edge._1eg.dofs.do.main import _3dCSCG_E1F_Dofs_Do


class _3dCSCG_E1F_Dofs(FrozenOnly):
    """"""
    def __init__(self, ef):
        """"""
        self._ef_ = ef
        self._do_ = None
        self._freeze_self_()

    @property
    def do(self):
        if self._do_ is None:
            self._do_ = _3dCSCG_E1F_Dofs_Do(self)
        return self._do_

