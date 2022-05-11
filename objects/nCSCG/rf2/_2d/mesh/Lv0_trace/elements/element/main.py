# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/08 4:08 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class _2nCSCG_Lv0TraceElement(FrozenOnly):
    """"""

    def __init__(self, Lv0elements, mesh, i):
        """"""
        self._elements_ = Lv0elements
        self._te_ = mesh.cscg.trace.elements[i] # the corresponding trace-element in the cscg mesh.
        self._freeze_self_()

    @property
    def i(self):
        return self._te_.i

    @property
    def segments(self):
        """The segments on this lv0-trace-element."""
        return self._elements_.segments[self.i]


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
