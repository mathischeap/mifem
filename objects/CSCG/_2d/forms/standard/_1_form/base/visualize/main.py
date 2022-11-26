# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly
from objects.CSCG._2d.forms.standard._1_form.base.visualize.matplot import _2dCSCG_S1F_VIS_Matplot


class _2dCSCG_S1F_VIS(FrozenOnly):
    """"""
    def __init__(self, sf):
        self._sf_ = sf
        self._matplot_ = _2dCSCG_S1F_VIS_Matplot(sf)
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        return self._matplot_(*args, **kwargs)

    @property
    def matplot(self):
        return self._matplot_