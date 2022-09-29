# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 9/24/2022 10:17 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')
from screws.freeze.main import FrozenOnly
from importlib import import_module


class miUs_Triangular_SF_Numbering(FrozenOnly):
    """"""

    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self._routine_ = 'Naive'
        self._gathering_ = None
        self._freeze_self_()

    @property
    def gathering(self):
        if self._gathering_ is None:
            path = '.'.join(str(self.__class__).split(' ')[1][1:-2].split('.')[:-2]) + '.'
            path += self._routine_
            name = self._routine_
            number = getattr(import_module(path), name)(self._sf_)
            self._gathering_ = getattr(number, self._sf_.__class__.__name__)
        return self._gathering_







if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
