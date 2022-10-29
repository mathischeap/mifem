# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/13 9:32 PM
"""

import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.mpRfT.base.forms.base.partial.dofs.main import PartialDofs
from objects.mpRfT.base.forms.base.partial.cochain.main import PartialCochain


class mpRfT_Form_BoundaryCondition(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._analytic_expression_ = None
        self._valid_boundaries_ = None
        self._pd_ = None
        self._pc_ = None
        self._freeze_self_()

    @property
    def analytic_expression(self):
        return self._analytic_expression_

    @analytic_expression.setter
    def analytic_expression(self, analytic_expression):
        """"""
        self._f_.___Pr_check_BC_analytic_expression___(analytic_expression)
        self._analytic_expression_ = analytic_expression

    @property
    def valid_boundaries(self):
        return self._valid_boundaries_

    @valid_boundaries.setter
    def valid_boundaries(self, valid_boundaries):
        """"""
        if isinstance(valid_boundaries, str):
            valid_boundaries = [valid_boundaries,]
        bns = self._f_.mesh.boundaries.names
        for vb in valid_boundaries:
            assert vb in bns, f"valid boundary {vb} is not a boundary name."
        self._valid_boundaries_ = valid_boundaries
        self._pd_ = None
        self._pc_ = None

    @property
    def partial_dofs(self):
        assert self.valid_boundaries is not None, f"no valid boundaries"
        if self._pd_ is None:
            self._pd_ = PartialDofs(self._f_)
            self._pd_.include.boundaries(self.valid_boundaries)
        return self._pd_

    @property
    def partial_cochain(self):
        assert self.valid_boundaries is not None, f"no valid boundaries"
        if self._pc_ is None:
            self._pc_ = PartialCochain(self._f_)
            self._pc_.include.boundaries(self.valid_boundaries)
        return self._pc_




if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
