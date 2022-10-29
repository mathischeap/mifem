# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/15 8:32 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.mpRfT._2d.forms.standard._0.base.numbering.do.find import mpRfT2_S0F_Numbering_DO_FIND


class mpRfT2_S0F_Numbering_DO(FrozenOnly):
    """"""

    def __init__(self, numbering):
        """"""
        self._numbering_ = numbering
        self._find_ = mpRfT2_S0F_Numbering_DO_FIND(numbering)
        self._freeze_self_()

    @property
    def find(self):
        return self._find_

if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
