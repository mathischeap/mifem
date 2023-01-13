# -*- coding: utf-8 -*-

from components.freeze.main import FrozenOnly
from objects.CSCG._2d.mesh.boundaries.visualize.matplot import _2dCSCG_Mesh_Boundaries_Matplot


class _2dCSCG_Mesh_Boundaries_Visualize(FrozenOnly):
    def __init__(self, bdrs):
        self._bdrs_ = bdrs
        self._matplot_ = None
        self._freeze_self_()

    def __call__(self, **kwargs):
        return self.matplot(**kwargs)

    @property
    def matplot(self):
        if self._matplot_ is None:
            self._matplot_ = _2dCSCG_Mesh_Boundaries_Matplot(self._bdrs_)
        return self._matplot_
