# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 5/25/2022 9:35 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

import numpy as np
from screws.freeze.base import FrozenOnly


class mpRfT2_T1F_Reconstruct(FrozenOnly):
    """"""

    def __init__(self, t):
        """"""
        self._t_ = t
        self._freeze_self_()

    def __call__(self, coo_map, ravel=False, rp=None, value_only=False):
        """

        Parameters
        ----------
        coo_map
        ravel
        rp :
            In which root-cells we are going to reconstruct the standard form.

        Returns
        -------

        """


if __name__ == '__main__':
    # mpiexec -n 4 python 
    from __init__ import rfT2

    fc = rfT2.rf(100)

    t1 = fc('1-t')
