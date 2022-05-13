# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/12 2:48 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class _3nCSCG_PolySpaceVisualize(FrozenOnly):
    """"""
    def __init__(self, space):
        """"""
        self._space_ = space
        self._freeze_self_()




if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
