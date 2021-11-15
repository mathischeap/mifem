# -*- coding: utf-8 -*-
"""
The interpolation maps [0, 1]^3 into the region.

Yi Zhang (C)
Created on Fri Dec 14 20:13:27 2018
Aerodynamics, AE
TU Delft
"""
from SCREWS.frozen import FrozenOnly
from importlib import import_module


class InterpolationSearcher(FrozenOnly):
    """ """
    def __init__(self, ID):
        """ """
        assert ID in self.___coded_interpolator___(), \
            " <InterpolationSearcher> : interpolation named '{}' is not coded.".format(ID)
        self._ID_ = ID
        cls_name = self.___coded_interpolator___()[ID]
        cls_path = self.___interpolator_path___()
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
        return "_3dCSCG.mesh.region.interpolation.interpolations"