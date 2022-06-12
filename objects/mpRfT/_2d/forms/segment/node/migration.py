# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/12 2:58 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class mpRfT2_NSgF_Migration(FrozenOnly):
    """"""

    def __init__(self, nsg):
        """"""
        self._nsg_ = nsg
        self._freeze_self_()


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
