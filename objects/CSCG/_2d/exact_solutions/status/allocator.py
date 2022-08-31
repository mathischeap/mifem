# -*- coding: utf-8 -*-
from screws.freeze.base import FrozenOnly


class _2dCSCG_ExactSolutionAllocator(FrozenOnly):
    """"""



    @classmethod
    def ___exact_solution_name___(cls):
        return {'sL:sincos1'              : 'SinCos1',
                'Euler:shear_layer_rollup': 'ShearLayerRollup',
                'icpsNS:TGV'              : 'TaylorGreenVortex',
                'icpsNS:dipole_collision' : 'DipoleCollision',
        }

    @classmethod
    def ___exact_solution_path___(cls):
        base_path = '.'.join(str(cls).split(' ')[1][1:-2].split('.')[:-2]) + '.'
        return {'sL:sincos1'              : base_path + 'scalar_Laplace.SinCos',
                'Euler:shear_layer_rollup': base_path + 'Euler.shear_layer_rollup',
                'icpsNS:TGV'              : base_path + 'incompressible_Navier_Stokes.Taylor_Green_vortex',
                'icpsNS:dipole_collision' : base_path + 'incompressible_Navier_Stokes.dipole_collision',
        }