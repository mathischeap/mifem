# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/30/2022 11:50 PM
"""
from components.freeze.base import FrozenOnly
from objects.CSCG._2d.discreteFields.scalar.visualize.matplot import _2cCSCG_DS_VisualizeMatplot


class _2cCSCG_DS_Visualize(FrozenOnly):
    """"""

    def __init__(self, ds):
        """"""
        self._ds_ = ds
        self._matplot_ = _2cCSCG_DS_VisualizeMatplot(ds)
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        return self.matplot(*args, **kwargs)

    @property
    def matplot(self):
        return self._matplot_
