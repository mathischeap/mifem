# -*- coding: utf-8 -*-
import numpy as np
from objects.CSCG._2d.mesh.domain.regions.region.edge_geometries.base import EdgeGeometryBase


class Straight(EdgeGeometryBase):
    def __init__(self, cc, et):
        super().__init__(cc, et)
        assert self.edge_type == ('straight',), \
            " <SideGeometryPlane> : side_type={} wrong.".format(self.edge_type)
        self._melt_self_()
        self._x1_, self._y1_ = self.corner_coordinates[0]
        self._x2_, self._y2_ = self.corner_coordinates[1]
        self._freeze_self_()

    def X(self, o):
        return self._x1_ + o * (self._x2_ - self._x1_)

    def Xo(self, o):
        return (self._x2_ - self._x1_) * np.ones(np.shape(o))

    def Y(self, o):
        return self._y1_ + o * (self._y2_ - self._y1_)

    def Yo(self, o):
        return (self._y2_ - self._y1_) * np.ones(np.shape(o))