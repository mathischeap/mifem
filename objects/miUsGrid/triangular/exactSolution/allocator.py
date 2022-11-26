# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/10/08 8:01 PM
"""

from components.freeze.base import FrozenOnly
from root.config.main import RANK, MASTER_RANK
from importlib import import_module

class miUsGrid_ExactSolutionAllocator(FrozenOnly):
    """"""
    @classmethod
    def listing(cls, printing=True, returning=True):
        """"""
        if RANK != MASTER_RANK: return

        listing = ''
        names = cls.___es_name___()
        paths = cls.___es_path___()
        for ID in names:
            CLASS = getattr(import_module(paths[ID]), names[ID])
            doc = CLASS.__doc__
            if doc is not None and doc != '':
                doc = doc.split('\n')[0]
                listing += ">>> " + ID + ' ~ ' + doc + '\n\n'
            else:
                listing += ">>> " + ID + '\n\n'

        if printing:
            print(listing)
        else:
            pass

        if returning:
            return listing
        else:
            pass

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



