# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/10/20 2:20 PM
"""
from abc import ABC

from tools.miLinearAlgebra.preconditioners.base import Preconditioner


class spiLU(Preconditioner, ABC):
    """"""
    def __init__(self, A, drop_tol=1e-4, fill_factor=10, drop_rule='basic'):
        """To improve the better approximation to the inverse, you may need to increase
        fill_factor AND decrease drop_tol."""
        super(spiLU, self).__init__(A)
        self._A_ = None
        self._drop_tol_ = drop_tol
        self._fill_factor_ = fill_factor
        self._drop_rule_ = drop_rule
        self._freeze_self_()

    @property
    def drop_tol(self):
        """(0,1), the lower, the more accurate, the slower."""
        return self._drop_tol_

    @property
    def fill_factor(self):
        """[1,), the larger, the more accurate, the slower."""
        return self._fill_factor_

    @property
    def drop_rule(self):
        return self._drop_rule_
