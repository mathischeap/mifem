# -*- coding: utf-8 -*-
from numpy import sin, pi
from objects.CSCG._2d.exact_solutions.status.scalar_Laplace.base import scalar_Laplace_Base



class SinCos1(scalar_Laplace_Base):
    """"""
    def __init__(self, es):
        super(SinCos1, self).__init__(es)

    def p(self, t, x, y): return sin(2*pi*x) * sin(2*pi*y)

