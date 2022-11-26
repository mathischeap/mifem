# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/19 7:53 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from components.freeze.base import FrozenOnly
from objects.CSCG._2d.forms.trace.base.numbering.do.find import _2dCSCG_Trace_Numbering_DO_FIND


class _2dCSCG_Trace_Numbering_DO(FrozenOnly):
    """"""

    def __init__(self, numbering):
        """"""
        self._numbering_ = numbering
        self._find_ = _2dCSCG_Trace_Numbering_DO_FIND(numbering)
        self._freeze_self_()

    @property
    def find(self):
        return self._find_
if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
