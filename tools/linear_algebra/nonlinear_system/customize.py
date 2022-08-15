# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/8/9 14:28
"""
import sys

if './' not in sys.path: sys.path.append('./')
from screws.freeze.main import FrozenOnly


class nLS_Customize(FrozenOnly):
    def __init__(self, nLS):
        """"""
        self._nLS_ = nLS
        self.___customizations___ = list()
        self._freeze_self_()

    @property
    def customizations(self):
        """The current existing customizations."""
        return self.___customizations___

    def set_no_evaluation(self, r):
        """Let the nonlinear system do not affect the value of #r dof.

        So dof_i will always be equal to dof_0 (the initial value (or initial guess)).
        """
        self.___customizations___.append(('set_no_evaluation', r))


if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
