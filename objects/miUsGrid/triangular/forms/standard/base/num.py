# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 4:22 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class miUs_Triangular_SF_Num(FrozenOnly):
    """"""

    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self._freeze_self_()

    @property
    def basis(self):
        return getattr(self._sf_.space.num_basis, self._sf_.__class__.__name__)


if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/triangular/forms/standard/base/num.py
    pass
