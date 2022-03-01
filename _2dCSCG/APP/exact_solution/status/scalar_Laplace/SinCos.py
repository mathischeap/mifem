
from numpy import sin, cos, pi
from _2dCSCG.APP.exact_solution.status.scalar_Laplace.base import scalar_Laplace_Base



class SinCos1(scalar_Laplace_Base):
    """"""
    def __init__(self, es):
        super(SinCos1, self).__init__(es)

    def p(self, t, x, y): return sin(2*pi*x) * sin(2*pi*y)
    def p_x(self, t, x, y): return 2*pi * cos(2*pi*x) * sin(2*pi*y)
    def p_xx(self, t, x, y): return - 4*pi**2 * sin(2*pi*x) * sin(2*pi*y)
    def p_y(self, t, x, y): return 2*pi * sin(2*pi*x) * cos(2*pi*y)
    def p_yy(self, t, x, y): return - 4*pi**2 * sin(2*pi*x) * sin(2*pi*y)

