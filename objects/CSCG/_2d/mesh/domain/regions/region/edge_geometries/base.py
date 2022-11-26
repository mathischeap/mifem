# -*- coding: utf-8 -*-
"""
INTRO

@author: Yi Zhang. Created on Tue May 21 18:07:50 2019
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft,
         Delft, the Netherlands

"""
import numpy as np
from components.freeze.main import FrozenOnly


class EdgeGeometryBase(FrozenOnly):
    def __init__(self, corner_coordinates, edge_type):
        self._corner_coordinates_ = corner_coordinates
        self._edge_type_ = edge_type
        assert np.shape(self.corner_coordinates) == (2, 2), \
            " <SideGeometryPlane> : corner_coordinates={} wrong.".format(self.corner_coordinates)
        self._freeze_self_()

    @property
    def corner_coordinates(self):
        return self._corner_coordinates_

    @property
    def edge_type(self):
        return self._edge_type_


    def X(self, o):
        raise NotImplementedError()

    def Xo(self, o):
        raise NotImplementedError()

    def Y(self, o):
        raise NotImplementedError()

    def Yo(self, o):
        raise NotImplementedError()

    def XY(self, o):
        return self.X(o), self.Y(o)

    def XoYo(self, o):
        return self.Xo(o), self.Yo(o)

