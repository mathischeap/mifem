# -*- coding: utf-8 -*-
from screws.freeze.base import FrozenOnly


from objects.CSCG._3d.forms.edge.base.dofs.do.find import _3dCSCG_EdgeDofs_DoFIND


class _3dCSCG_EdgeDofs_Do(FrozenOnly):
    """"""
    def __init__(self, dofs):
        self._dofs_ = dofs
        self._find_ = None
        self._freeze_self_()

    @property
    def find(self):
        if self._find_ is None:
            self._find_ = _3dCSCG_EdgeDofs_DoFIND(self._dofs_)
        return self._find_

