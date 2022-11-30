# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/8/4 21:12
"""
import sys

if './' not in sys.path: sys.path.append('./')
from components.freeze.main import FrozenOnly


class nLS_num(FrozenOnly):
    def __init__(self, nLS):
        """"""
        self._nLS_ = nLS
        self._freeze_self_()

    @property
    def equations(self):
        return self._nLS_._num_equations_

    @property
    def dofs(self):
        return self._nLS_._num_dofs_

if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
