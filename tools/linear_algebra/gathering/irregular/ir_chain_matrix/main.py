# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/23 10:27 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class iR_Chain_Gathering_Matrix(FrozenOnly):
    """"""

    def __init__(self, GMs):
        """"""
        self._GMs_ = GMs
        self._freeze_self_()

    @property
    def ___Pr_IS_regular___(self):
        return False




if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
