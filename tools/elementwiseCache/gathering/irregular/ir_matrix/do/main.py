# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 6/7/2022 3:59 PM
"""

from components.freeze.base import FrozenOnly
from tools.elementwiseCache.gathering.irregular.ir_matrix.do.find import iR_Gathering_Matrix_DoFind


class iR_Gathering_Matrix_DO(FrozenOnly):
    """"""

    def __init__(self, igm):
        """"""
        self._igm_ = igm
        self._find_ = iR_Gathering_Matrix_DoFind(igm)
        self._freeze_self_()

    @property
    def find(self):
        return self._find_
