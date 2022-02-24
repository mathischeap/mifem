# -*- coding: utf-8 -*-
"""
INTRO

@author: Yi Zhang. Created on Tue May 21 18:07:50 2019
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft
         Delft, Netherlands

"""
import numpy as np
from SCREWS.frozen import FrozenOnly
from SCREWS.functions._2d import angle


class EdgeGeometry(FrozenOnly):
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





class AACW(EdgeGeometry):
    """ Arc AntiClockWise."""
    def __init__(self, cc, et):
        super().__init__(cc, et)
        assert self.edge_type[0] == 'aacw', \
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
        if self.end_theta < self.start_theta: self.end_theta += 2 * np.pi
        self._freeze_self_()

    def X(self, o):
        """ """
        theta = o * (self.end_theta - self.start_theta) + self.start_theta
        return self.xc + self.r * np.cos(theta)

    def Xo(self, o):
        theta = o * (self.end_theta - self.start_theta) + self.start_theta
        return -self.r * np.sin(theta) * (self.end_theta - self.start_theta)

    def Y(self, o):
        """ """
        theta = o * (self.end_theta - self.start_theta) + self.start_theta
        return self.yc + self.r * np.sin(theta)

    def Yo(self, o):
        theta = o * (self.end_theta - self.start_theta) + self.start_theta
        return self.r * np.cos(theta) * (self.end_theta - self.start_theta)


class ACW(EdgeGeometry):
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




class Customized(EdgeGeometry):
    """ """
    def __init__(self, cc, et):
        """ """
        super().__init__(cc, et)
        assert self.edge_type[0] == 'customized', \
            " <SideGeometryPlane> : edge_type={} wrong.".format(self.edge_type[0])
        self._melt_self_()
        self._mapping_ = self.edge_type[1]
        self._Jacobian_ = self.edge_type[2]
        self._freeze_self_()

    def X(self, o):
        return self._mapping_[0](o)

    def Xo(self, o):
        return self._Jacobian_[0](o)

    def Y(self, o):
        return self._mapping_[1](o)

    def Yo(self, o):
        return self._Jacobian_[1](o)




class Free(EdgeGeometry):
    """A free edge geometry is a line we do not really care about its classfication.

    A regions of such edge geometry(s) usually will call a specific interpolator to
    generate the mapping and therefore Jacobian and so on, like the crazy mapping. The
    crazy mapping is a analytical mapping that we only need to know the bounds and c
    which is stored in the `domain_input`.

    While for the transfinite mapping, the edge geometries is essential. We can not set
    it to be free.

    BUT, even for free EdgeGeometry (or SideGeometry in 3D), the corner_coordinates are
    still necessary because they are used to study the topology of the regions.

    Free edge does not mean it is straight! it can be everything (we just donot care).

    """
    def __init__(self, cc, st):
        """ """
        super().__init__(cc, st)
        assert self.edge_type == ('free',), \
            " <EdgeGeometry> <Free> : edge_type={} wrong.".format(self.edge_type)

    def X(self, o):
        raise NotImplementedError()

    def Xo(self, o):
        raise NotImplementedError()

    def Y(self, o):
        raise NotImplementedError()

    def Yo(self, o):
        raise NotImplementedError()



class Straight(EdgeGeometry):
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