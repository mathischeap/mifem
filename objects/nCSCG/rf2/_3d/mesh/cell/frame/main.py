# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/07 3:44 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.nCSCG.rf2._3d.mesh.cell.frame.facets.main import FrameFacets


class _3nCSCG_CellFrame(FrozenOnly):
    """"""

    def __init__(self, cell):
        """"""
        assert cell.mesh._locker_, f"mesh should be locked before accessing frame. Update mesh first."
        self._cell_ = cell
        self._N = None
        self._S = None
        self._W = None
        self._E = None
        self._B = None
        self._F = None
        self._freeze_self_()

    @property
    def N(self):
        if self._N is None:
            self._N = FrameFacets(self._cell_, 'N')
        return self._N

    @property
    def S(self):
        if self._S is None:
            self._S = FrameFacets(self._cell_, 'S')
        return self._S

    @property
    def W(self):
        if self._W is None:
            self._W = FrameFacets(self._cell_, 'W')
        return self._W

    @property
    def E(self):
        if self._E is None:
            self._E = FrameFacets(self._cell_, 'E')
        return self._E

    @property
    def B(self):
        if self._B is None:
            self._B = FrameFacets(self._cell_, 'B')
        return self._B

    @property
    def F(self):
        if self._F is None:
            self._F = FrameFacets(self._cell_, 'F')
        return self._F






if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
