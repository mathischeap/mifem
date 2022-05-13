# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/12 11:57 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from screws.freeze.base import FrozenOnly


class _2nCSCG_SpacePolyDo(FrozenOnly):
    """"""

    def __init__(self, space):
        """"""
        self._space_ = space
        self._freeze_self_()


    def add(self, N):
        """add the space of degree N to the pool."""
        pool = self._space_.___pool___
        assert isinstance(N, int) and N > 0, f"N={N} ({N.__class__.__name__}) is wrong, I need a positive integer."
        assert N not in pool, f"N={N} function space already exists."
        new_space = self._space_.___space_class___((self._space_.___node_type___, N), ndim=2)
        self._space_.___pool___[N] = new_space




if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
