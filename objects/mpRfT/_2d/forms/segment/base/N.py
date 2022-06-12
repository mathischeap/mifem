# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/02 8:46 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly



class mpRfT2_SgF_N(FrozenOnly):
    """"""

    def __init__(self, t):
        """"""
        self._t_ = t
        self._freeze_self_()


    def __getitem__(self, seg):
        """"""
        N = seg.N

        ndp = self._t_.ndp

        if isinstance(ndp, int):
            return N + ndp
        else:
            raise NotImplementedError(f"npp={ndp} not implemented.")





if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
