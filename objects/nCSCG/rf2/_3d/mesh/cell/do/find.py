# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/08 1:50 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class _3nCSCG_CellDoFind(FrozenOnly):
    """"""

    def __init__(self, cell):
        """"""
        self._cell_ = cell
        self._freeze_self_()


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
