# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/7/21 15:21
"""
import sys

if './' not in sys.path: sys.path.append('./')
from screws.freeze.main import FrozenOnly

from scipy.sparse import csr_matrix, csc_matrix

class ____PRIVATE_SparseMatrix_caller___(FrozenOnly):
    def __init__(self, mdm, transpose):
        """"""
        self._mdm_ = mdm
        self._T_ = transpose
        self._freeze_self_()

    def __call__(self, basic_unit):
        """"""
        if self._T_:
            return csc_matrix(self._mdm_[basic_unit]).T
        else:
            return csr_matrix(self._mdm_[basic_unit])

if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
