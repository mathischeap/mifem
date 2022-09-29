# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 4:34 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class miUs_Triangular_SF_Do(FrozenOnly):
    """"""

    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self._freeze_self_()

    def evaluate_basis_at_meshgrid(self, *args, **kwargs):
        """"""
        return getattr(self._sf_.space.evaluation, self._sf_.__class__.__name__)(*args, **kwargs)



if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
