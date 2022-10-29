# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 10/6/2022 10:18 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from screws.freeze.main import FrozenOnly
from root.config.main import COMM, MPI


class miUsTriangle_Domain(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._area_ = None
        self._freeze_self_()

    @property
    def area(self):
        """The area of the domain."""
        if self._area_ is None:
            area = 0
            for e in self._mesh_.elements:
                element = self._mesh_.elements[e]
                area += element.area
            self._area_ = COMM.allreduce(area, MPI.SUM)
        return self._area_





if __name__ == '__main__':
    # mpiexec -n 4 python objects/miUsGrid/triangular/mesh/domain/main.py

    from __init__ import miTri
    fc = miTri.form('st2', 2)
    mesh = fc.mesh

    print(mesh.domain.area)
