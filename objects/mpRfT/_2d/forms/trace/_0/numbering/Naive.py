# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 5/25/2022 9:19 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class Naive(FrozenOnly):
    """"""

    def __init__(self, t):
        """"""
        self._t_ = t
        self._freeze_self_()

if __name__ == '__main__':
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/trace/_0/numbering/Naive.py
    from __init__ import rfT2

    fc = rfT2.rf(100)

    t0 = fc('0-t')
