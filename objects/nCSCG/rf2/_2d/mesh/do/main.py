# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.nCSCG.rf2._2d.mesh.do.find import _2nCSCG_MeshDoFind



class _2nCSCG_MeshDo(FrozenOnly):
    """"""
    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._find_ = None
        self._freeze_self_()

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




if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
