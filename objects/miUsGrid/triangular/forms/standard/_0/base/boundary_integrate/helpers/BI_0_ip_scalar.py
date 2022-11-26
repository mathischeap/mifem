# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 10/7/2022 7:17 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from components.freeze.main import FrozenOnly


class miUsTriangle_BI_S0F_ip_Scalar(FrozenOnly):
    """"""

    def __init__(self, s0f, scalar):
        """"""
        self._s0f_ = s0f
        self._scalar_ = scalar
        self._freeze_self_()


if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
