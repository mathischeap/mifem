# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/08/30 4:19 PM
"""
from components.freeze.base import FrozenOnly
from objects.CSCG._2d.discreteFields.vector.visualize.matplot import _2cCSCG_DV_VisualizeMatplot


class _2dCSCG_DV_Visualize(FrozenOnly):
    """"""

    def __init__(self, dv):
        """"""
        self._dv_ = dv
        self._matplot_ = _2cCSCG_DV_VisualizeMatplot(dv)
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        return self.matplot(*args, **kwargs)

    @property
    def matplot(self):
        return self._matplot_
