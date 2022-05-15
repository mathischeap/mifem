# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/15 3:20 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class _2nCSCG_MRF2_CellSpaceDo(FrozenOnly):
    """"""

    def __init__(self, space):
        """"""
        self._space_ = space
        self._freeze_self_()






if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
