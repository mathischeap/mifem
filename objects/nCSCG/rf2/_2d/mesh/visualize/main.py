# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/05 12:40 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.nCSCG.rf2._2d.mesh.visualize.matplot import _2nCSCG_MeshVisualizeMatplot


class _2nCSCG_MeshVisualize(FrozenOnly):
    """"""
    def __init__(self, mesh):
        """"""
        self._mesh_ = mesh
        self._matplot_ = None
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        return self.matplot(*args, **kwargs)

    @property
    def matplot(self):
        if self._matplot_ is None:
            self._matplot_ = _2nCSCG_MeshVisualizeMatplot(self._mesh_)
        return self._matplot_






if __name__ == "__main__":
    # mpiexec -n 8 python objects/nCSCG/rfT2/_2d/mesh/visualize/main.py
    pass
