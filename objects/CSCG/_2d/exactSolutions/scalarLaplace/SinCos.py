# -*- coding: utf-8 -*-
from numpy import sin, pi
from objects.CSCG._2d.exactSolutions.scalarLaplace.base import scalar_Laplace_Base



class SinCos1(scalar_Laplace_Base):
    """A sin-cos test case for the scalar Laplacian."""
    def __init__(self, mesh):
        super(SinCos1, self).__init__(mesh)

    def p(self, t, x, y): return sin(2*pi*x) * sin(2*pi*y)

