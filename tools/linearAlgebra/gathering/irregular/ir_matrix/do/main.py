# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 6/7/2022 3:59 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from tools.linearAlgebra.gathering.irregular.ir_matrix.do.find import iR_Gathering_Matrix_DoFind

class iR_Gathering_Matrix_DO(FrozenOnly):
    """"""

    def __init__(self, igm):
        """"""
        self._igm_ = igm
        self._find_ = iR_Gathering_Matrix_DoFind(igm)
        self._freeze_self_()

    @property
    def find(self):
        return self._find_




if __name__ == '__main__':
    # mpiexec -n 4 python tools/linear_algebra/gathering/irregular/ir_matrix/do/main.py
    pass
