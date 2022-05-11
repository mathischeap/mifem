# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/08 4:05 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.nCSCG.rf2._2d.mesh.Lv0_trace.elements.do import _2dCSCG_MeshLv0TraceElementsDo
from objects.nCSCG.rf2._2d.mesh.Lv0_trace.elements.element.main import _2nCSCG_Lv0TraceElement


class _2nCSCG_Lv0TraceElements(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._do_ = None
        self._elements_ = dict()
        self._segments_ = None
        self._freeze_self_()

    @property
    def do(self):
        if self._do_ is None:
            self._do_ = _2dCSCG_MeshLv0TraceElementsDo(self)
        return self._do_

    def __iter__(self):
        """Go through all local level-0-trace-elements"""
        for i in self._mesh_.cscg.trace.elements:
            yield i

    def __getitem__(self, i):
        if i in self._elements_:
            pass
        else:
            self._elements_[i] = _2nCSCG_Lv0TraceElement(self, self._mesh_, i)

        return self._elements_[i]

    def __contains__(self, i):
        return i in self._mesh_.cscg.trace.elements

    @property
    def map(self):
        """dict : keys are level-0 cells, values are level-0-trace-elements around them."""
        return self._mesh_.cscg.trace.elements.map

    @property
    def segments(self):
        """"""
        assert self._segments_ is not None, f"segments is None, first do update."
        return self._segments_




if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
