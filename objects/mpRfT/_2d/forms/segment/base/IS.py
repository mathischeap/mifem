# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 5/25/2022 9:20 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly



class mpRfT2_SgF_IS(FrozenOnly):
    """"""

    def __init__(self, t):
        """"""
        self._t_ = t
        self._freeze_self_()



if __name__ == '__main__':
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/segment/base/IS.py
    pass
