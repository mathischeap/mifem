# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/12 5:02 PM
"""

from components.freeze.base import FrozenOnly
from tools.elementwiseCache.gathering.irregular.ir_chain_matrix.do.find import iR_CGM_DO_FIND


class iR_CGM_DO(FrozenOnly):
    """"""

    def __init__(self, ir_cgm):
        """"""
        self._CGM_ = ir_cgm
        self._find_ = iR_CGM_DO_FIND(ir_cgm)
        self._freeze_self_()

    @property
    def find(self):
        return self._find_
