# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.nCSCG.rf2._3d.mesh.cell.main import _3nCSCG_RF2_MeshCell


class _3nCSCG_RF2_BaseMeshCells(FrozenOnly):
    """"""
    def __init__(self, mesh):
        self._mesh_ = mesh
        self._cells_ = dict()
        self._freeze_self_()

    def __getitem__(self, i):
        if i in self._cells_:
            pass
        elif i in self._mesh_.cscg.elements:
            self._cells_[i] = _3nCSCG_RF2_MeshCell(self._mesh_, 0, (i,))
        else:
            raise Exception()
        return self._cells_[i]

    def __iter__(self):
        """go through all local 0-level cells (cscg mesh elements.)"""
        for i in self._mesh_.cscg.elements:
            yield i






if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
