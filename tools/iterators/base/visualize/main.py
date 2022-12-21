# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/16/2022 1:21 AM
"""
import sys

if './' not in sys.path:
    sys.path.append('./')

from components.freeze.base import FrozenOnly


class IteratorVisualize(FrozenOnly):
    """"""

    def __init__(self, iterator):
        """"""
        self._iterator_ = iterator
        self._freeze_self_()

    @property
    def iterator(self):
        return self._iterator_
