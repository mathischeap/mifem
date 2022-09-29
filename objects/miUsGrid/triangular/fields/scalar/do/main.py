# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 2:55 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class miUsGrid_Triangular_Scalar_Do(FrozenOnly):
    """A collection of other methods beside the standard ones like reconstruct, discretize, visualize, etc."""

    def __init__(self, scalar):
        """"""
        self._scalar_ = scalar
        self._freeze_self_()


    def evaluate_func_at_time(self, time=None):
        return self._scalar_.___DO_evaluate_func_at_time___(time=time)



if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
