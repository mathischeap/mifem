# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/10 6:39 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from root.config.main import cOmm, MPI


class _2nCSCG_Mesh_RF2_IS(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._freeze_self_()

    @property
    def locked(self):
        return self._mesh_._locker_

    @property
    def basic(self):
        """Return True if all base cells are root cells."""
        assert self.locked, f"lock the mesh first before checking if it is basic."
        basic = True
        for i in self._mesh_.base_cells:
            base_cell = self._mesh_.base_cells[i]
            if base_cell.IS.root:
                pass
            else:
                basic = False
                break
        basic = cOmm.allreduce(basic, op=MPI.LAND)
        return basic




if __name__ == "__main__":
    # mpiexec -n 4 python objects/nCSCG/rfT2/_2d/mesh/IS.py
    pass
