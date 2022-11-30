# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/7/23 12:11
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly




class ___NEG___(FrozenOnly):
    def __init__(self, mdm):
        self._mdm_ = mdm
        self._freeze_self_()

    def __call__(self, item):
        return - self._mdm_[item]

if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
