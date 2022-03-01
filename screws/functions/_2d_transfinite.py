# -*- coding: utf-8 -*-
"""
These two-dimensional tools are currently working for three-dimensional cases. For
example, for a three-dimensional meshComponents, if on the third-dimension, it is just a
transformation mapping, then we can use a two-dimensional transfinite mapping to get 
the mapping.

They will be deprecated when the 2D transfinite interpolation and the edge functions
are officially coded. But it looks like this will never be happening.

Now used:
    1). 3D `bridge_arch_cracked` interpolation.
    2). 3D `bridge_arch` interpolation.
    3). 2D `transfinite` interpolation.

<unittests> <2D> <unittests_Interpolations> <test_No1_2D_TransfiniteMapping>.
    
Yi Zhang (C)
Created on Mon Feb 25 15:56:49 2019
Aerodynamics, AE
TU Delft
"""
import numpy as np
from screws.decorators import accepts
from screws.frozen import FrozenOnly
from screws.numerical._2d import NumericalJacobian_xy_t_21




class TransfiniteMapping(FrozenOnly):
    """"""
    def __init__(self, gamma, dgamma=None):
        """
         y          - 2 +
         ^     _______________
         |    |       R       |
         |    |               | 
         |  + |               | +
         |  3 | U           D | 1
         |  - |               | -
         |    |       L       | 
         |    |_______________|
         |          - 0 +
         |_______________________> x
        
        The indices of `gamma` and `dgamma` are as above. And the directs of the
        mappings are indicated as well.
        
        Parameters
        ----------
        gamma : 
            A tuple of the four boundary functions
        dgamma : optional
            A tuple of first derivative of gamma.
            
        """
        t = np.linspace(0,1,12)[1:-1]
        _dict_ = {0:'L', 1:'D', 2:'R', 3:'U'}
        for i in range(4):
            XY = gamma[i]
            XtYt = dgamma[i]
            NJ21 =NumericalJacobian_xy_t_21(XY)
            assert all(NJ21.check_Jacobian(XtYt, t)), \
            " <TransfiniteMapping> :  '{}' edge mapping or Jacobian wrong.".format(_dict_[i])
            
        self.gamma = gamma
        self.dgamma = dgamma
        self.gamma1_x0, self.gamma1_y0 = self.gamma[0](0.0)
        self.gamma1_x1, self.gamma1_y1 = self.gamma[0](1.0)
        self.gamma3_x0, self.gamma3_y0 = self.gamma[2](0.0)
        self.gamma3_x1, self.gamma3_y1 = self.gamma[2](1.0)
    
    def mapping(self, r, s):
        """
        mapping (r, s) = (0, 1)^2 into (x, y) using the transfinite mapping.
        
        """
        gamma1_xs, gamma1_ys = self.gamma[0](r)
        gamma2_xt, gamma2_yt = self.gamma[1](s)
        gamma3_xs, gamma3_ys = self.gamma[2](r)
        gamma4_xt, gamma4_yt = self.gamma[3](s)
        x = (1-r)*gamma4_xt + r*gamma2_xt + (1-s)*gamma1_xs + s*gamma3_xs - \
            (1-r)*((1-s)*self.gamma1_x0 + s*self.gamma3_x0) - r*((1-s)*self.gamma1_x1 + s*self.gamma3_x1)
        y = (1-r)*gamma4_yt + r*gamma2_yt + (1-s)*gamma1_ys + s*gamma3_ys - \
            (1-r)*((1-s)*self.gamma1_y0 + s*self.gamma3_y0) - r*((1-s)*self.gamma1_y1 + s*self.gamma3_y1)
        return x, y
        
    def dx_dr(self, r, s):
        """ """
        gamma2_xt, gamma2_yt = self.gamma[1](s)
        gamma4_xt, gamma_4yt = self.gamma[3](s)
        dgamma1_xds, dgamma1_yds = self.dgamma[0](r)
        dgamma3_xds, dgamma3_yds = self.dgamma[2](r)
        dx_dxi_result = (-gamma4_xt + gamma2_xt + (1-s)*dgamma1_xds + s*dgamma3_xds + 
                ((1-s)*self.gamma1_x0 + s*self.gamma3_x0) - ((1-s)*self.gamma1_x1 + s*self.gamma3_x1) )
        return dx_dxi_result
    
    def dx_ds(self, r, s):
        """ """
        gamma1_xs, gamma1_ys = self.gamma[0](r)
        gamma3_xs, gamma3_ys = self.gamma[2](r)
        dgamma2_xdt, dgamma2_ydt = self.dgamma[1](s)
        dgamma4_xdt, dgamma4_ydt = self.dgamma[3](s)
        dx_deta_result = ((1-r)*dgamma4_xdt + r*dgamma2_xdt - gamma1_xs + gamma3_xs - 
                (1-r)*(-self.gamma1_x0 + self.gamma3_x0) - r*(-self.gamma1_x1 + self.gamma3_x1))
        return dx_deta_result
    
    def dy_dr(self, r, s):
        """ """
        gamma2_xt, gamma2_yt = self.gamma[1](s)
        gamma4_xt, gamma4_yt = self.gamma[3](s)
        dgamma1_xds, dgamma1_yds = self.dgamma[0](r)
        dgamma3_xds, dgamma3_yds = self.dgamma[2](r)
        dy_dxi_result = (-gamma4_yt + gamma2_yt + (1-s)*dgamma1_yds + s*dgamma3_yds + 
                ((1-s)*self.gamma1_y0 + s*self.gamma3_y0) - ((1-s)*self.gamma1_y1 + s*self.gamma3_y1)) 
        return dy_dxi_result
    
    def dy_ds(self, r, s):
        """ """
        gamma1_xs, gamma1_ys = self.gamma[0](r)
        gamma3_xs, gamma3_ys = self.gamma[2](r)
        dgamma2_xdt, dgamma2_ydt = self.dgamma[1](s)
        dgamma4_xdt, dgamma4_ydt = self.dgamma[3](s)
        dy_deta_result = ((1-r)*dgamma4_ydt + r*dgamma2_ydt - gamma1_ys + gamma3_ys - 
                (1-r)*(-self.gamma1_y0 + self.gamma3_y0) - r*(-self.gamma1_y1 + self.gamma3_y1))
        return dy_deta_result




# FIT INTO A STRAIGHT LINE FROM TWO POINTS
class StraightLine(object):
    @accepts('self', (tuple, list), (tuple, list))
    def __init__(self, start_point, end_point):
        assert np.shape(start_point) == np.shape(end_point) == (2,)
        self.x1, self.y1 = start_point
        self.x2, self.y2 = end_point
        
    # o in [0, 1]
    def _gamma(self, o): 
        return self.x1 + o*(self.x2-self.x1), self.y1 + o*(self.y2-self.y1)
    
    def _dgamma(self, o): 
        return (self.x2-self.x1) * np.ones(np.shape(o)), (self.y2-self.y1) * np.ones(np.shape(o))
    
    def __call__(self):
        return self._gamma, self._dgamma


# ANGLE BETWEEN TWO LINES
@accepts((tuple, list), (tuple, list))
def Angle(origin, pt):
    """
    angle between the vector from origin to pt and the x-direction vector.
    """
    x1, y1 = (1, 0)
    x2, y2 = (pt[0]-origin[0], pt[1]-origin[1])
    inner_product = x1*x2 + y1*y2
    len1 = np.hypot(x1, y1)
    len2 = np.hypot(x2, y2)
    if y2 < 0:
        return 2*np.pi - np.arccos(inner_product/(len1*len2))
    else:
        return np.arccos(inner_product/(len1*len2))




# FIT INTO anti-clock-wise ARC WITH:center, start_point, end_point
class ArcAntiClockWise(object):
    """
    fit two center and radius int arc (up half) the return arc is ALWAYS
    anti-clock-wise!!
    """
    @accepts('self', (tuple, list), (tuple, list), (tuple, list))
    def __init__(self, center, start_point, end_point):
        """ """
        self.x0, self.y0 = center
        x1, y1 = start_point
        x2, y2 = end_point
        self.r = np.sqrt((x1-self.x0)**2 + (y1-self.y0)**2)
        assert np.abs(self.r - (np.sqrt((x2-self.x0)**2 + (y2-self.y0)**2))) < 10e-13, \
            'center is not at proper place'
        self.start_theta = Angle(center, start_point)
        self.end_theta = Angle(center, end_point)
        if self.end_theta < self.start_theta: self.end_theta += 2*np.pi
        
    # o in [0, 1]
    def _gamma(self, o):
        theta = o*(self.end_theta-self.start_theta) + self.start_theta
        return self.x0 + self.r*np.cos(theta), self.y0 + self.r*np.sin(theta)
    
    def _dgamma(self, o):
        theta = o*(self.end_theta-self.start_theta) + self.start_theta
        return -self.r*np.sin(theta) * (self.end_theta-self.start_theta), \
                self.r*np.cos(theta) * (self.end_theta-self.start_theta)
                
    def __call__(self):
        return self._gamma, self._dgamma




# FIT INTO clock-wise ARC WITH:center, start_point, end_point
class ArcClockWise(object):
    """ 
    fit two center and radius int arc (up half). The return arc is ALWAYS 
    clock-wise!!
    """
    @accepts('self', (tuple, list), (tuple, list), (tuple, list))
    def __init__(self, center, start_point, end_point):
        self.x0, self.y0 = center
        x1, y1 = start_point
        x2, y2 = end_point
        self.r = np.sqrt((x1-self.x0)**2 + (y1-self.y0)**2)
        assert np.abs(self.r - (np.sqrt((x2-self.x0)**2 + (y2-self.y0)**2))) < 10e-13, \
            'center is not at proper place'
        self.start_theta = Angle(center, start_point)
        self.end_theta = Angle(center, end_point)
        # o in [0, 1]
        if self.end_theta > self.start_theta: 
            self.end_theta -= 2*np.pi
            
    def _gamma(self, o):
        theta = o*(self.end_theta-self.start_theta) + self.start_theta
        return self.x0 + self.r*np.cos(theta), self.y0 + self.r*np.sin(theta)
    
    def _dgamma(self, o):
        theta = o*(self.end_theta-self.start_theta) + self.start_theta
        return -self.r*np.sin(theta) * (self.end_theta-self.start_theta), \
                self.r*np.cos(theta) * (self.end_theta-self.start_theta)
    
    def __call__(self):
        return self._gamma, self._dgamma