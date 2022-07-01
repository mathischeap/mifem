# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/12 6:30 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from objects.mpRfT._2d.forms.segment.node.dofs.visualization import mpRfT2_NSgF_Dofs_Visualization




class mpRfT2_NSgF_Dofs(FrozenOnly):
    """"""

    def __init__(self, t):
        """"""
        self._t_ = t
        self._visualization_ = mpRfT2_NSgF_Dofs_Visualization(self)
        self._freeze_self_()

    @property
    def visualization(self):
        return self._visualization_




if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/segment/node/dofs/main.py
    pass
