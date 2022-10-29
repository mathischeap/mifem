# -*- coding: utf-8 -*-
from screws.freeze.base import FrozenOnly



class _2dCSCG_SpaceAllocator(FrozenOnly):
    """"""

    @classmethod
    def ___space_name___(cls):
        return {'polynomials': "_2dCSCG_PolynomialSpace"}


    @classmethod
    def ___space_path___(cls):
        base_path = '.'.join(str(cls).split(' ')[1][1:-2].split('.')[:-2]) + '.'
        return {'polynomials': base_path + "polynomials"}
