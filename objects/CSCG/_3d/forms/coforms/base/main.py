# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/08/25 2:08 PM
"""

from components.freeze.base import FrozenOnly


class classname(FrozenOnly):
    """"""

    def __init__(self):
        """"""
        self._freeze_self_()
