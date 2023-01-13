# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/27/2022 9:52 PM
"""
from components.freeze.main import FrozenOnly


class _3dCSCG_1LocalTrace_Visualize(FrozenOnly):
    """"""

    def __init__(self, ltf):
        """"""
        self._ltf_ = ltf
        self._freeze_self_()
