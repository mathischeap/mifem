# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/30/2022 11:50 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.CSCG._2d.discrete_fields.scalar.visualize.matplot import _2cCSCG_DS_VisualizeMatplot

class _2cCSCG_DS_Visualize(FrozenOnly):
    """"""

    def __init__(self, ds):
        """"""
        self._ds_ = ds
        self._matplot_ = _2cCSCG_DS_VisualizeMatplot(ds)
        self._freeze_self_()

    @property
    def matplot(self):
        return self._matplot_


if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
