# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/12 1:20 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class mpRfT2_ESgF_Num(FrozenOnly):
    """"""

    def __init__(self, t):
        """"""
        self._t_ = t
        self._basis_ = ___Pr_Basis___(t)
        self._freeze_self_()

    @property
    def basis(self):
        return self._basis_

    @property
    def local_dofs(self):
        if self._t_.numbering._num_local_dofs_ is None:
            _ = self._t_.numbering.sgW_gathering
        return self._t_.numbering._num_local_dofs_


class ___Pr_Basis___(FrozenOnly):
    """"""
    def __init__(self, t):
        self._t_ = t
        self._freeze_self_()

    def __getitem__(self, seg):
        """"""
        assert seg.__class__.__name__ == 'mpRfT2_Segment', f"I need a mpRfT2_Segment"
        N = self._t_.N[seg]
        return N



if __name__ == "__main__":
    # mpiexec -n 4 python
    pass
