# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/08 3:09 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.nCSCG.rf2._2d.mesh.Lv0_trace.elements.main import _2nCSCG_Lv0TraceElements
from objects.nCSCG.rf2._2d.mesh.Lv0_trace.visualize import _2nCSCG_Lv0TraceVisualize



class _2nCSCG_MeshLv0Trace(FrozenOnly):
    """The trace mesh of the cscg mesh (level-0-cells)."""
    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._elements_ = _2nCSCG_Lv0TraceElements(mesh)
        self._visualize_ = None
        self._freeze_self_()

    @property
    def elements(self):
        return self._elements_

    @property
    def visualize(self):
        if self._visualize_ is None:
            self._visualize_ = _2nCSCG_Lv0TraceVisualize(self)
        return self._visualize_






if __name__ == "__main__":
    # mpiexec -n 4 python objects/nCSCG/rf2/_2d/mesh/Lv0_trace/main.py

    from root.read.main import read
    mesh = read('test_mesh.mi')

    trace = mesh.Lv0trace
    # print(trace.elements.segments)
    trace.visualize()