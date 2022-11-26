# -*- coding: utf-8 -*-
"""
With SideGeometryBase, we store the mapping form (p, q) to the side geometry, while
(p, q)= ((0,1), (0,1)).

Yi Zhang (C)
Created on Tue Sep  4 21:30:37 2018
Aerodynamics, AE
TU Delft
"""
import numpy as np
from components.freeze.main import FrozenOnly

class SideGeometryBase(FrozenOnly):
    def __init__(self, corner_coordinates, side_type):
        self._corner_coordinates_ = corner_coordinates
        self._side_type_ = side_type
        assert np.shape(self.corner_coordinates) == (4, 3), \
            " <SideGeometryPlane> : corner_coordinates={} wrong.".format(self.corner_coordinates)
        self._freeze_self_()
    
    @property
    def corner_coordinates(self):
        return self._corner_coordinates_
    
    @property
    def side_type(self):
        return self._side_type_
