# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.nCSCG.rf2.base.mesh.do.main import nCSCG_MeshDoBase
from objects.nCSCG.rf2._2d.mesh.do.find import _2nCSCG_MeshDoFind
from time import time
from random import random

class _2nCSCG_MeshDo(nCSCG_MeshDoBase):
    """"""
    def __init__(self, mesh):
        """"""
        super(_2nCSCG_MeshDo, self).__init__(mesh)
        self._find_ = None
        self._digest_ = None
        self._freeze_self_()

    def update(self):
        """We update the mesh and clear all relevant cached data and properties."""
        self._mesh_.Lv0trace.elements.do._Pr_update()
        self._mesh_.base_cells.internal_segments.do._Pr_update()
        self._mesh_.segments.do._Pr_update()
        for ind in self._mesh_:
            self._mesh_(ind)._frame_ = None
        self._mesh_._signature_ = time() + random()
        for obj in self._mesh_.___ecosystem___:
            obj.do.update()
        self.lock() # lock self after updating!

    @property
    def find(self):
        if self._find_ is None:
            self._find_ = _2nCSCG_MeshDoFind(self._mesh_)
        return self._find_

    def refine(self, indices):
        """We will refine the cell(s) indicated by (group of) `indices`."""
        if isinstance(indices, int):
            self._mesh_(indices).do.refine()
        elif isinstance(indices, tuple) and isinstance(indices[0], int):
            # we get the indices of one cell
            self._mesh_(indices).do.refine()
        else:
            for ind in indices:
                self._mesh_(ind).do.refine()

    def dilute(self, indices):
        """We will dilute the cell(s) indicated by (group of) `indices`."""
        if isinstance(indices, int):
            self._mesh_(indices).do.dilute()
        elif isinstance(indices, tuple) and isinstance(indices[0], int):
            # we get the indices of one cell
            self._mesh_(indices).do.dilute()
        else:
            for ind in indices:
                self._mesh_(ind).do.dilute()

    @staticmethod
    def digest(refinement):
        """We digest a refinement object to refine the mesh accordingly."""
        refinement.___Pr_apply___() # all checks will be done in the particular refinement.









if __name__ == '__main__':
    # mpiexec -n 4 python objects/nCSCG/rfT2/_2d/mesh/do/main.py
    from objects.nCSCG.rf2._2d.master import MeshGenerator

    mesh = MeshGenerator('crazy')([3, 3], 2, EDM='chaotic', show_info=True)
    mesh.do.unlock()

    if 0 in mesh.cscg.elements:
        c0 = mesh(0)
        c0.do.refine()

    mesh.do.update()