# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/24 12:30 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly


class mpRfT2_Si1F_Discretize_Standard_Vector(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._freeze_self_()


    def __call__(self, target):
        """

        Parameters
        ----------
        target : str
            {'func',}

        Returns
        -------

        """

if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
