# -*- coding: utf-8 -*-

from components.freeze.base import FrozenOnly
from objects.CSCG._2d.forms.standard.base.dofs.dof.main import _2dCSCG_SF_DOF
from objects.CSCG._2d.forms.standard.base.dofs.visualize import _2dCSCG_SF_dofs_VIS
from objects.CSCG._2d.forms.standard.base.dofs.do.main import _2dCSCG_SF_dofs_do


class _2dCSCG_SF_dofs(FrozenOnly):
    """"""
    def __init__(self, sf):
        self._sf_ = sf
        self._dofs_ = dict()
        self._visualize_ = None
        self._do_ = None
        self._freeze_self_()

    def __getitem__(self, i):
        """return and cache the global dof #i."""
        if i not in self._dofs_:
            assert i % 1 == 0 and 0 <= i < self.GLOBAL_num, \
                f"i={i} is out of range, we only have {self.GLOBAL_num} dofs."
            self._dofs_[i] = _2dCSCG_SF_DOF(self, i)
        return self._dofs_[i]

    def __iter__(self):
        for i in range(self.GLOBAL_num):
            yield i

    def __contains__(self, i):
        """If dof #`i` is contained in this core?"""
        if i % 1 != 0:
            return False
        else:
            if 0 <= i < self.GLOBAL_num:
                return True
            else:
                return False

    @property
    def do(self):
        if self._do_ is None:
            self._do_ = _2dCSCG_SF_dofs_do(self)
        return self._do_

    @property
    def visualize(self):
        if self._visualize_ is None:
            self._visualize_ = _2dCSCG_SF_dofs_VIS(self)
        return self._visualize_

    @property
    def GLOBAL_num(self):
        return self._sf_.num.global_dofs
