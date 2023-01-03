# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/08/29 9:17 PM
"""
from components.freeze.base import FrozenOnly


class _2dCSCG_Element_Constant(FrozenOnly):
    """"""

    def __init__(self, element):
        """"""
        self._element_ = element
        self._freeze_self_()
