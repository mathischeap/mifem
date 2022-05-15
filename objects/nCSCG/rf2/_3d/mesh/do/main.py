# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.nCSCG.rf2.base.mesh.do.main import nCSCG_MeshDoBase
from objects.nCSCG.rf2._3d.mesh.do.find import _3nCSCG_MeshDoFind
from time import time
from random import random


class _3nCSCG_MeshDo(nCSCG_MeshDoBase):
    """"""
    def __init__(self, mesh):
        """"""
        super(_3nCSCG_MeshDo, self).__init__(mesh)
        self._find_ = None
        self._freeze_self_()

    @property
    def find(self):
        if self._find_ is None:
            self._find_ = _3nCSCG_MeshDoFind(self._mesh_)
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

    def update(self):
        """We update the mesh and clear all relevant cached data and properties."""
        self._mesh_.facets.do._Pr_update()
        for ind in self._mesh_:
            cell = self._mesh_(ind)
            cell._frame_ = None
        self._mesh_._signature_ = time() + random()
        self.lock() # lock self after update!



if __name__ == '__main__':
    # mpiexec -n 4 python objects/nCSCG/rf2/_3d/mesh/do/main.py
    pass
