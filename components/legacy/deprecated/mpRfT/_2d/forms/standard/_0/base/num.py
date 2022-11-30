# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/19 7:49 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly


class mpRfT2_S0F_Num(FrozenOnly):
    """"""

    def __init__(self, f):
        """"""
        self._f_ = f
        self._basis_ = ___mpRfT2_S0F_Num_Basis___(f)
        self._freeze_self_()

    @property
    def basis(self):
        return self._basis_

    @property
    def components(self):
        """0-form has 1 components."""
        return 1

    @property
    def local_dofs(self):
        if self._f_.numbering._num_local_dofs_ is None:
            _ = self._f_.numbering.gathering
        return self._f_.numbering._num_local_dofs_

    @property
    def GLOBAL_dofs(self):
        return self._f_.numbering.gathering.global_num_dofs



class ___mpRfT2_S0F_Num_Basis___(FrozenOnly):
    """"""
    def __init__(self, f):
        self._f_ = f
        self._cache_ = dict()
        self._freeze_self_()


    def __getitem__(self, rp):
        """"""
        N = self._f_.mesh[rp].N
        if N not in self._cache_:
            self._cache_[N] = (N+1)**2
        return self._cache_[N]




if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
