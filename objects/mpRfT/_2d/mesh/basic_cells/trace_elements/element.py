# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 5/22/2022 10:37 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class mpRfT2_Mesh_BasicCells_TraceElement(FrozenOnly):
    """"""

    def __init__(self, elements, mesh, i):
        """"""
        self._elements_ = elements
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
