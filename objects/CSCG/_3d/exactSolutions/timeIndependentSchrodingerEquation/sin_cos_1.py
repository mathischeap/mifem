# -*- coding: utf-8 -*-
from numpy import sin, pi
from objects.CSCG._3d.exactSolutions.timeIndependentSchrodingerEquation.base import \
    TimeIndependentSchrodingerEquationBase



class TISE_SinCos1(TimeIndependentSchrodingerEquationBase):
    """The sin cos test case 1.
    """
    def __init__(self, mesh, m, E):
        super(TISE_SinCos1, self).__init__(mesh, m, E)


    def psi(self, t, x, y, z):
        return sin(2*pi*x) * sin(2*pi*y) * sin(2*pi*z)

    def V(self, t, x, y, z):
        """"""
        return self.E - 12 * pi**2 * self._alpha_