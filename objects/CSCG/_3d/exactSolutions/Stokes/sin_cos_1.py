# -*- coding: utf-8 -*-
from numpy import sin, pi, cos
from objects.CSCG._3d.exactSolutions.Stokes.base import Stokes_Base


# noinspection PyAbstractClass
class Stokes_SinCos1(Stokes_Base):
    """The sin cos test case 1.
    """
    def __init__(self, mesh):
        super(Stokes_SinCos1, self).__init__(mesh)

    @property
    def valid_time(self):
        """Return None because this exact solution is valid at any time."""
        return None

    def p(self, t, x, y, z):
        return sin(2*pi*x) * sin(2*pi*y) * sin(2*pi*z)


    def u(self, t, x, y, z):
        return cos(2*pi*x) * sin(2*pi*y) * sin(2*pi*z)

    def v(self, t, x, y, z):
        return sin(2*pi*x) * cos(2*pi*y) * sin(2*pi*z)

    def w(self, t, x, y, z):
        return - 2 * sin(2*pi*x) * sin(2*pi*y) * cos(2*pi*z)
