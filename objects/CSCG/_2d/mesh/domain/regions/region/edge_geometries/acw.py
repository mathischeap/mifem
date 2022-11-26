# -*- coding: utf-8 -*-


from components.functions._2dSpace.angle import angle
import numpy as np
from objects.CSCG._2d.mesh.domain.regions.region.edge_geometries.base import EdgeGeometryBase

class ACW(EdgeGeometryBase):
    """ Arc ClockWise."""
    def __init__(self, cc, et):
        super().__init__(cc, et)
        assert self.edge_type[0] == 'acw', \
            " <SideGeometryPlane> : edge_type={} wrong.".format(self.edge_type[0])
        self._melt_self_()
        self.x1, self.y1 = self.corner_coordinates[0]
        self.x2, self.y2 = self.corner_coordinates[1]
        self.xc, self.yc = self.edge_type[1]
        self.r = np.sqrt((self.x1 - self.xc) ** 2 + (self.y1 - self.yc) ** 2)
        assert np.abs(self.r - (np.sqrt((self.x2 - self.xc) ** 2 +
                                        (self.y2 - self.yc) ** 2))) < 10e-13, 'center is not at proper place.'
        self.start_theta = angle((self.xc, self.yc), (self.x1, self.y1))
        self.end_theta = angle((self.xc, self.yc), (self.x2, self.y2))
        if self.end_theta > self.start_theta: self.end_theta -= 2 * np.pi
        self._freeze_self_()

    def X(self, o):
        theta = o * (self.end_theta - self.start_theta) + self.start_theta
        return self.xc + self.r * np.cos(theta)

    def Xo(self, o):
        theta = o * (self.end_theta - self.start_theta) + self.start_theta
        return -self.r * np.sin(theta) * (self.end_theta - self.start_theta)

    def Y(self, o):
        theta = o * (self.end_theta - self.start_theta) + self.start_theta
        return self.yc + self.r * np.sin(theta)

    def Yo(self, o):
        theta = o * (self.end_theta - self.start_theta) + self.start_theta
        return self.r * np.cos(theta) * (self.end_theta - self.start_theta)

