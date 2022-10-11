# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 9/29/2022 10:42 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.miUsGrid.triangular.forms.standard.base.numbering.do.find import miUsTriangle_Numbering_DoFind


class miUsTriangle_Numbering_Do(FrozenOnly):
    """"""

    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self._find_ = None
        self._freeze_self_()

    @property
    def find(self):
        if self._find_ is None:
            self._find_ = miUsTriangle_Numbering_DoFind(self._sf_)
        return self._find_

if __name__ == '__main__':
    # mpiexec -n 4 python objects/miUsGrid/triangular/forms/standard/base/numbering/do/main.py
    pass
