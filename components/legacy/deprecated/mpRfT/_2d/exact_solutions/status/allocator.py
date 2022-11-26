# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/17 3:46 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly



class mpRfT2_ExactSolutionAllocator(FrozenOnly):
    """"""

    @classmethod
    def ___exact_solution_name___(cls):
        return {'Poisson:sincos1': 'Poisson_SinCos1',
                }

    @classmethod
    def ___exact_solution_path___(cls):
        base_path = '.'.join(str(cls).split(' ')[1][1:-2].split('.')[:-2]) + '.'
        return {'Poisson:sincos1': base_path + 'Poisson.sin_cos',
                }


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
