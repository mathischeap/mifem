# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly
from root.config.main import RANK, MASTER_RANK
from importlib import import_module


class _2dCSCG_ExactSolutionAllocator(FrozenOnly):
    """"""


    @classmethod
    def listing(cls, printing=True, returning=True):
        """"""
        if RANK != MASTER_RANK: return

        listing = ''
        names = cls.___exact_solution_name___()
        paths = cls.___exact_solution_path___()
        for ID in names:
            CLASS = getattr(import_module(paths[ID]), names[ID])
            doc = CLASS.__doc__
            if doc is not None and doc != '':
                doc = doc.split('\n')[0]
                listing += '>>> ' + ID + ' ~ ' + doc + '\n\n'
            else:
                listing += '>>> ' + ID + '\n\n'

        if printing:
            print(listing)
        else:
            pass

        if returning:
            return listing
        else:
            pass

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
        return {'sL:sincos1'              : base_path + 'scalarLaplace.SinCos',
                'Euler:shear_layer_rollup': base_path + 'Euler.shear_layer_rollup',
                'icpsNS:TGV'              : base_path + 'incompressibleNavierStokes.Taylor_Green_vortex',
                'icpsNS:dipole_collision' : base_path + 'incompressibleNavierStokes.dipole_collision',
        }