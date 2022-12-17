# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/27/2022 9:52 PM
"""
from components.freeze.main import FrozenOnly
from objects.CSCG._3d.forms.localTrace._2ltf.visualize.matplot import _3dCSCG_2LocalTrace_Visualize_Matplot


class _3dCSCG_2LocalTrace_Visualize(FrozenOnly):
    """"""

    def __init__(self, ltf):
        """"""
        self._ltf_ = ltf
        self._matplot_ = None
        self._freeze_self_()

    def __call__(self, *args, **kwargs):
        return self.matplot(*args, **kwargs)

    @property
    def matplot(self):
        if self._matplot_ is None:
            self._matplot_ = _3dCSCG_2LocalTrace_Visualize_Matplot(self._ltf_)
        return self._matplot_
