# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from importlib import import_module



class _3nCSCG_SpaceAllocator(FrozenOnly):
    """"""
    def __init__(self, ID):
        """"""
        assert ID in self.___space_name___(), f"space_type={ID} is wrong!"
        self._NAME_ = self.___space_name___()[ID]
        self._PATH_ = self.___space_path___()[ID]
        self._freeze_self_()

    def __call__(self, dN, **kwargs):
        """"""
        CLASS = getattr(import_module(self._PATH_), self._NAME_)
        return CLASS(dN, **kwargs)

    @classmethod
    def ___space_name___(cls):
        return {'polynomials': "_3nCSCG_PolynomialSpace"}


    @classmethod
    def ___space_path___(cls):
        base_path = '.'.join(str(cls).split(' ')[1][1:-2].split('.')[:-2]) + '.'
        return {'polynomials': base_path + "polynomials.main"}






if __name__ == '__main__':
    # mpiexec -n 4 python objects/nCSCG/rfT2/_2d/space/allocator.py
    pass
