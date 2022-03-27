# -*- coding: utf-8 -*-
"""
The interpolation maps [0, 1]^3 into the regions.

Yi Zhang (C)
Created on Fri Dec 14 20:13:27 2018
Aerodynamics, AE
TU Delft
"""
from screws.freeze.main import FrozenOnly
from importlib import import_module


class InterpolationAllocator(FrozenOnly):
    """ """
    def __init__(self, ID):
        """ """
        assert ID in self.___coded_interpolator___(), \
            " <InterpolationSearcher> : interpolation named '{}' is not coded.".format(ID)
        self._ID_ = ID
        cls_name = self.___coded_interpolator___()[ID]
        cls_path = self.___interpolator_path___()[ID]
        self._itp_ = getattr(import_module(cls_path), cls_name)
        self._freeze_self_()

    def __call__(self, region):
        """ """
        return self._itp_(region)
        
    @classmethod
    def ___coded_interpolator___(cls):
        """ Update this whenever we code a new Interpolator. """
        return {'crazy': 'Crazy',
                'transfinite': 'Transfinite',
                'bridge_arch_cracked': 'BridgeArchCracked',}
    
    @classmethod
    def ___interpolator_path___(cls):
        """ """
        base_path = '.'.join(str(cls).split(' ')[1][1:-2].split('.')[:-2]) + '.'
        return {'crazy': base_path + 'crazy',
                'transfinite': base_path + 'transfinite',
                'bridge_arch_cracked': base_path + 'bridge_arch_cracked',}