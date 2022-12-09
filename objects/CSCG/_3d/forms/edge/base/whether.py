# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 12/3/2022 9:59 AM
"""
from components.freeze.main import FrozenOnly


class _3dCSCG_EdgeForm_Whether(FrozenOnly):
    """"""

    def __init__(self, ef):
        """"""
        self._ef_ = ef
        self._freeze_self_()

    @property
    def hybrid(self):
        return self._ef_._hybrid_
