# -*- coding: utf-8 -*-
""" Internal segments of a base-cell plus the segments on the lv0-trace-element of the base cell
together are the all segments of a base cell.

@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/09 12:04 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.nCSCG.rf2._2d.mesh.base_cells.internal_segments.do import BaseCellsInternalSegmentsDo

class BaseCellInternalSegments(FrozenOnly):
    """"""
    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._do_ = BaseCellsInternalSegmentsDo(mesh)
        self._segments_ = dict() # should be updated before accessing
        self._2bud_Sgs_ = dict() # to be updated dict
        for i in self._mesh_.cscg.elements:
            self._2bud_Sgs_[i] = True
        self._freeze_self_()

    @property
    def do(self):
        """"""
        return self._do_

    def __getitem__(self, item):
        """"""
        return self._segments_[item]

    def __iter__(self):
        """go through all local 0-level cells (cscg mesh elements.)"""
        for i in self._mesh_.cscg.elements:
            yield i

    def __contains__(self, item):
        return item in self._mesh_.cscg.elements









if __name__ == "__main__":
    # mpiexec -n 4 python objects/nCSCG/rfT2/_2d/mesh/base_cells/internal_segments/main.py
    from objects.nCSCG.rf2._2d.__tests__.Random.mesh import random_mesh_of_elements_around as rm2
    mesh = rm2(100, refinement_intensity=0.5)

    # from root.read.main import read
    # mesh = read('test_mesh.mi')


