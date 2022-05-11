# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly



class _2nCSCG_SpaceAllocator(FrozenOnly):
    """"""



    @classmethod
    def ___space_name___(cls):
        return {'polynomials': "_2nCSCG_PolynomialSpace"}


    @classmethod
    def ___space_path___(cls):
        base_path = '.'.join(str(cls).split(' ')[1][1:-2].split('.')[:-2]) + '.'
        return {'polynomials': base_path + "polynomials"}



if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
