# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/15 9:32 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly


class mpRfT2_NSgF_Numbering_Local(FrozenOnly):
    """"""

    def __init__(self, numbering):
        """"""
        self._numbering_ = numbering
        self._cache_ = dict()
        self._freeze_self_()

    def __getitem__(self, rc_rp):
        """"""
        raise NotImplementedError()




if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
