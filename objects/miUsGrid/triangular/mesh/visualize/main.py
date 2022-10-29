# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/19 3:45 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.miUsGrid.triangular.mesh.visualize.matplot import miUsGrid_TriangularMesh_Matplot


class miUsGrid_TriangularMesh_Visualize(FrozenOnly):
    """"""

    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._matplot_ = None
        self._freeze_self_()

    @property
    def matplot(self):
        """"""
        if self._matplot_ is None:
            self._matplot_ = miUsGrid_TriangularMesh_Matplot(self._mesh_)
        return self._matplot_

    def __call__(self, *args, **kwargs):
        """"""
        return self.matplot(*args, **kwargs)






if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/triangular/mesh/visualize/main.py
    pass
