# -*- coding: utf-8 -*-
from screws.freeze.base import FrozenOnly



class _3dCSCG_ExactSolution_Allocator(FrozenOnly):
    """"""


    @classmethod
    def ___exact_solution_name___(cls):
        """"""
        return {'icpsNS:TGV1'    : 'TGV1',
                'icpsNS:LDC'     : 'LidDrivenCavity',
                'icpsNS:sincosRC': 'SinCosRebholz_Conservation',
                'icpsNS:sincos_CCBF'  : 'SinCos_Conservation_Conservative_Body_Force',
                'icpsNS:sincos_CCBF1' : 'SinCos_Conservation_Conservative_Body_Force1',
                'icpsNS:sincos_CCBF_P': 'SinCos_Conservation_Conservative_Body_Force_POLYNOMIALS',
                'icpsNS:sincos_CCBF_C': 'SinCos_Conservation_Conservative_Body_Force_CONSTANT',
                'icpsNS:sincosRD': 'SinCosRebholz_Dissipation',
                'icpsNS:sincosMD': 'SinCos_Modified_Dissipation',
                'icpsNS:CUCD1': 'Closed_Unit_Cube_Disspation1',
                'icpsNS:CBFx' : 'Constant_X_direction_body_force',
                'icpsNS:still': 'Still',

                'Poisson:sincos1': "Poisson_SinCos1",

                'Stokes:sincos1': "Stokes_SinCos1",

                'TISE:sincos1': "TISE_SinCos1",

                'MHD:sincos1': "MHD_SinCos1",
                'MHD:as1': "AS1",

                }

    @classmethod
    def ___exact_solution_path___(cls):
        """"""
        base_path = '.'.join(str(cls).split(' ')[1][1:-2].split('.')[:-2]) + '.'
        _path_icpsNS_ =  base_path + 'incompressible_Navier_Stokes.'
        return {'icpsNS:TGV1'    : _path_icpsNS_ + 'Taylor_Green_vortex',
                'icpsNS:LDC'     : _path_icpsNS_ + 'lid_driven_cavity',
                'icpsNS:sincosRC': _path_icpsNS_ + 'Sin_Cos',
                'icpsNS:sincos_CCBF'  :  _path_icpsNS_ + 'Sin_Cos',
                'icpsNS:sincos_CCBF1' :  _path_icpsNS_ + 'Sin_Cos',
                'icpsNS:sincos_CCBF_P':  _path_icpsNS_ + 'Sin_Cos',
                'icpsNS:sincos_CCBF_C': _path_icpsNS_ + 'Sin_Cos',
                'icpsNS:sincosRD':  _path_icpsNS_ + 'Sin_Cos',
                'icpsNS:sincosMD':  _path_icpsNS_ + 'Sin_Cos',
                'icpsNS:CUCD1':  _path_icpsNS_ + 'others',
                'icpsNS:CBFx' :  _path_icpsNS_ + 'others',
                'icpsNS:still':  _path_icpsNS_ + 'others',

                'Poisson:sincos1': base_path + "Poisson.Sin_Cos",

                'Stokes:sincos1': base_path + "Stokes.sin_cos_1",

                'TISE:sincos1': base_path + "time_independent_Schrodinger_equation.sin_cos_1",

                'MHD:sincos1': base_path + "incompressible_MHD.sin_cos",
                'MHD:as1': base_path + "incompressible_MHD.analytic_solution_1",
                }