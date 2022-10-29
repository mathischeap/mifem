# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/10/08 8:01 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly


class miUsGrid_ExactSolutionAllocator(FrozenOnly):
    """"""

    @classmethod
    def ___es_name___(cls):
        return {'diNS : dipole collision': "DipoleCollision",
                'Stokes : sin cos 1': "Stokes_SinCos1",
                }

    @classmethod
    def ___es_path___(cls):
        base_path = '.'.join(str(cls).split(' ')[1][1:-2].split('.')[:-2]) + '.'
        return {'diNS : dipole collision': base_path + "dimensionlessIncompressibleNavierStokes.dipole_collision",
                'Stokes : sin cos 1': base_path + "Stokes.sin_cos_1",
                }





if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/triangular/exact_solution/allocator.py
    pass
