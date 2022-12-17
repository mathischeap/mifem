# -*- coding: utf-8 -*-
from components.freeze.main import FrozenOnly
from objects.CSCG._3d.fields.scalar.visualize.matplot import _3dCSCG_ScalarField_matplot_Visualize


class _3dCSCG_ScalarField_Visualize(FrozenOnly):
    """"""
    def __init__(self, f):
        """"""
        self._f_ = f
        self._matplot_ = _3dCSCG_ScalarField_matplot_Visualize(f)
        self._freeze_self_()


    def __call__(self, *args, **kwargs):
        """"""
        return self.matplot(*args, **kwargs)

    @property
    def matplot(self):
        return self._matplot_
