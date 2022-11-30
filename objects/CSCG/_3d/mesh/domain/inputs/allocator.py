# -*- coding: utf-8 -*-
"""
INTRO

Yi Zhang (C)
Created on Fri Dec 14 20:40:12 2018
Aerodynamics, AE
TU Delft
"""

from components.freeze.main import FrozenOnly
from importlib import import_module
from root.config.main import RANK, MASTER_RANK


class DomainInputAllocator(FrozenOnly):
    """ We use this finder to get a `DomainInput`."""
    def __init__(self, ID):
        assert ID in self.___defined_DI___(), f" <DomainInputFinder> : mesh ID = {ID} is wrong."
        cls_name = self.___defined_DI___()[ID]
        cls_path = self.___DI_path___()[ID]
        self._DomainInput_ = getattr(import_module(cls_path), cls_name)
        self._freeze_self_()
    
    def __call__(self, *args, **kwargs):
        """"""
        return self._DomainInput_(*args, **kwargs)

    @classmethod
    def listing(cls, printing=True, returning=True):
        """"""
        if RANK != MASTER_RANK: return

        listing = ''
        names = cls.___defined_DI___()
        paths = cls.___DI_path___()
        for ID in names:
            CLASS = getattr(import_module(paths[ID]), names[ID])
            doc = CLASS.__doc__
            if doc is not None and doc != '':
                doc = doc.split('\n')[0]
                listing += ">>> " + ID + ' ~ ' + CLASS.__name__ + ' -> ' + doc + '\n\n'
            else:
                listing += ">>> " + ID + ' ~ ' + CLASS.__name__ + '\n\n'

        if printing:
            print(listing)
        else:
            pass

        if returning:
            return listing
        else:
            pass

    @classmethod
    def ___defined_DI___(cls):
        """Here we store all defined meshComponents. Whenever we define a new meshComponents (actually, a new
        domain_input), we add a nickname for it here.
        
        """
        _dict_ = {
            'crazy': "Crazy",
            'crazy_periodic': "CrazyPeriodic",
            'bridge_arch_cracked': "BridgeArchCracked",
            'psc': "Periodic_Square_Channel",
            'pwc': "Parallel_Wall_Channel",
            'LDC': "Lid_Driven_Cavity",
            'cuboid': "Cuboid",
            'cuboid_periodic': "CuboidPeriodic",
            'ct' : 'CurvilinearTestMesh',
        }
        return _dict_

    @classmethod
    def ___DI_path___(cls):
        """ """
        base_path = '.'.join(str(cls).split(' ')[1][1:-2].split('.')[:-2]) + '.'
        return {
            'crazy': base_path + "crazy",
            'crazy_periodic': base_path + "crazy_periodic",
            'bridge_arch_cracked': base_path + "bridge_arch_cracked",
            'psc': base_path + "psc",
            'pwc': base_path + "pwc",
            'LDC': base_path + "LDC",
            'cuboid': base_path + "cuboid",
            'cuboid_periodic': base_path + "cuboid_periodic",
            'ct': base_path + "ct"
        }