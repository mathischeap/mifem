# -*- coding: utf-8 -*-
"""
Yi Zhang
zhangyi_aero@hotmail.com
created at: 2/10/2023 11:08 AM
"""
import sys

if './' not in sys.path:
    sys.path.append('./')
from objects.CSCG._3d.exactSolutions.incompressibleMHD.base import incompressible_MHD_Base

from numpy import sin, cos, pi, exp


class AS2(incompressible_MHD_Base):
    """The second analytical solution."""

    def __init__(self, mesh, Rf=10, Rm=100, c=0.01):
        """"""
        super(AS2, self).__init__(mesh, Rf, Rm, c)
        self._freeze_self_()


    def u(self, t, x, y, z):
        return sin(x) * cos(y) * cos(z) * exp(t/pi)

    def v(self, t, x, y, z):
        return cos(x) * sin(y) * cos(z) * exp(t/pi)

    def w(self, t, x, y, z):
        return -2 * cos(x) * cos(y) * sin(z) * exp(t/pi)


    def p(self, t, x, y, z):
        return sin(x) * sin(y) * sin(z) * exp(t/pi)


    def Bx(self, t, x, y, z):
        return cos(x) * cos(y) * cos(z) * exp(t/pi)

    def By(self, t, x, y, z):
        return sin(x) * sin(y) * cos(z) * exp(t/pi)

    def Bz(self, t, x, y, z):
        return cos(x) * sin(y) * exp(t/pi)


if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_3d/exactSolutions/incompressibleMHD/analytic_solution_2.py

    from objects.CSCG._3d.master import MeshGenerator, ExactSolutionSelector
    mesh = MeshGenerator('crazy', bounds=([0, 2*pi], [0, 2*pi], [0, 2*pi]), c=0.0)([5, 5, 5])
    es = ExactSolutionSelector(mesh)("MHD:as2", show_info=True)

    r = es.electric_source_term
    r.current_time = 1
    r.visualize(density=100000)
