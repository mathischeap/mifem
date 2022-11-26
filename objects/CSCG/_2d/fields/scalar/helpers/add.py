# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/22/2022 12:27 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from components.freeze.base import FrozenOnly


class _2dCSCG_ScalarFieldAddHelper(FrozenOnly):
    """"""

    def __init__(self, sf, of):
        """"""
        self._sf_ = sf
        self._of_ = of
        self._freeze_self_()

    def __call__(self, t, x, y):
        return self._of_(t, x, y) + self._sf_(t, x, y)




if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
