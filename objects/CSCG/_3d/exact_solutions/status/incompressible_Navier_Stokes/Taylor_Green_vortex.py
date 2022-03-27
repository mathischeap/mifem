# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from numpy import sin, cos
from objects.CSCG._3d.exact_solutions.status.incompressible_Navier_Stokes.base import incompressible_NavierStokes_Base


class TGV1(incompressible_NavierStokes_Base):
    """
    This is a 3D Taylor-Green Vortex initial condition. And for `t != 0`, the
    solution is not correct. So we can only use it as our initial condition.

    See the case in  [Inviscid and Viscous Simulations of the Taylor-Green Vortex Flow Using a
    Modal Discontinuous Galerkin Approach] with L=1, V0=1.

    Reynolds number
    """
    def __init__(self,es, nu=0, L=1, V0=1):
        self._L_ = L
        self._V0_ = V0
        super(TGV1, self).__init__(es, nu)

    @property
    def valid_time(self):
        return 'valid_only_at_its_first_instant'

    @property
    def L(self):
        return self._L_
    @property
    def V0(self):
        return self._V0_

    def u(self, t, x, y, z):
        return self.V0 * sin(x / self.L) * cos(y / self.L) * cos(z / self.L)
    def u_x(self, t, x, y, z):
        return (self.V0 / self.L) * cos(x / self.L) * cos(y / self.L) * cos(z / self.L)
    def u_y(self, t, x, y, z):
        return - (self.V0 / self.L) * sin(x / self.L) * sin(y / self.L) * cos(z / self.L)
    def u_z(self, t, x, y, z):
        return - (self.V0 / self.L) * sin(x / self.L) * cos(y / self.L) * sin(z / self.L)

    def v(self, t, x, y, z):
        return - self.V0 * cos(x / self.L) * sin(y / self.L) * cos(z / self.L)
    def v_x(self, t, x, y, z):
        return (self.V0 / self.L) * sin(x / self.L) * sin(y / self.L) * cos(z / self.L)
    def v_y(self, t, x, y, z):
        return - (self.V0 / self.L) * cos(x / self.L) * cos(y / self.L) * cos(z / self.L)
    def v_z(self, t, x, y, z):
        return (self.V0 / self.L) * cos(x / self.L) * sin(y / self.L) * sin(z / self.L)

    def w(self, t, x, y, z): return 0 * x
    def w_x(self, t, x, y, z): return 0 * x
    def w_y(self, t, x, y, z): return 0 * x
    def w_z(self, t, x, y, z): return 0 * x

    def fx(self, t, x, y, z): return 0 * x
    def fy(self, t, x, y, z): return 0 * x
    def fz(self, t, x, y, z): return 0 * x