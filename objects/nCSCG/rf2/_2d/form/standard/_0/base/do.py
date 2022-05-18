# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/17 3:57 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class _2nCSCG_RF2_S0F_Do(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._freeze_self_()

    def update(self):
        """"""
        if self._f_.cochain._local_ is not None:
            LCC = self._f_.cochain.local





if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
