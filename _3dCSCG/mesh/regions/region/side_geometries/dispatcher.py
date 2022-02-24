# -*- coding: utf-8 -*-
"""
With SideGeometry, we store the mapping form (p, q) to the side geometry, while (p, q)=
((0,1), (0,1)).

Yi Zhang (C)
Created on Fri Dec 14 20:28:15 2018
Aerodynamics, AE
TU Delft
"""
from SCREWS.frozen import FrozenOnly
from importlib import import_module


class SideGeometryDispatcher(FrozenOnly):
    """ """
    def __init__(self, st):
        """ """
        assert st[0] in self.___coded_side_geometries___(), \
            " <Region> : SideGeometry: {} is not coded. ".format(st[0])
        self._st_ = st
        cls_name = self.___coded_side_geometries___()[st[0]]
        cls_path = self.___side_geometries_path___() + '.' + st[0]
        self._cls_ = getattr(import_module(cls_path), cls_name)
        self._freeze_self_()
        
    def __call__(self, cc):
        """ """
        return self._cls_(cc, self._st_)

    @classmethod
    def ___coded_side_geometries___(cls):
        """ update this whenever we code a new SideGeometry. """
        return {'free': 'Free',
                'plane': 'Plane',
                'customized' : 'Customized',}
    
    @classmethod
    def ___side_geometries_path___(cls):
        """ """
        return '_3dCSCG.mesh.regions.region.side_geometries'