# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 11/28/2022 9:27 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from components.freeze.main import FrozenOnly
from objects.CSCG._3d.forms.localTrace.base.numbering.do.find import _3dCSCG_LocalTrace_Numbering_DoFind

class _3dCSCG_LocalTrace_Numbering_Do(FrozenOnly):
    """"""

    def __init__(self, ltf):
        """"""
        self._ltf_ = ltf
        self._find_ = _3dCSCG_LocalTrace_Numbering_DoFind(ltf)
        self._freeze_self_()


    @property
    def find(self):
        return self._find_


if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
