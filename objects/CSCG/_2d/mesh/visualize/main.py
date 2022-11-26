# -*- coding: utf-8 -*-

import sys
if './' not in sys.path: sys.path.append('./')

from components.freeze.main import FrozenOnly
from objects.CSCG._2d.mesh.visualize.matplot import _2dCSCG_Mesh_Visualize_Matplot



class _2dCSCG_Mesh_Visualize(FrozenOnly):
    def __init__(self, mesh):
        assert mesh.__class__.__name__ == '_2dCSCG_Mesh', " <MeshVisualize> "
        assert mesh.ndim == 2, " <MeshVisualize> "
        self._matplot_ = _2dCSCG_Mesh_Visualize_Matplot(mesh)
        self._freeze_self_()

    def __call__(self, **kwargs):
        return self.matplot(**kwargs)

    @property
    def matplot(self):
        return self._matplot_