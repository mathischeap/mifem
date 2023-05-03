# -*- coding: utf-8 -*-
"""
Yi Zhang
zhangyi_aero@hotmail.com
created at: 4/13/2023 5:22 PM
"""
import sys

if './' not in sys.path:
    sys.path.append('./')
from objects.CSCG._3d.exactSolutions.incompressibleMHD.base import incompressible_MHD_Base

from numpy import sin, cos


class AS3(incompressible_MHD_Base):
    """The second analytical solution."""

    def __init__(self, mesh, Rf=10, Rm=100, c=0.01):
        """"""
        super(AS3, self).__init__(mesh, Rf, Rm, c)
        self._freeze_self_()

    @property
    def valid_time(self):
        return 'valid_only_at_its_first_instant'

    def u(self, t, x, y, z):
        return sin(x) * cos(y) * cos(z)

    def v(self, t, x, y, z):
        return cos(x) * sin(y) * cos(z)

    def w(self, t, x, y, z):
        return -2 * cos(x) * cos(y) * sin(z)

    def Bx(self, t, x, y, z):
        return cos(x) * cos(y) * cos(z)

    def By(self, t, x, y, z):
        return sin(x) * sin(y) * cos(z)

    def Bz(self, t, x, y, z):
        return cos(x) * sin(y)

    def fx(self, t, x, y, z):
        """"""
        return 0 * x

    def fy(self, t, x, y, z):
        """"""
        return 0 * x

    def fz(self, t, x, y, z):
        """"""
        return 0 * x


if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_3d/exactSolutions/incompressibleMHD/analytic_solution_3.py

    from objects.CSCG._3d.master import MeshGenerator, ExactSolutionSelector
    mesh = MeshGenerator('crazy', bounds=([0, 2], [0, 2], [0, 2]), c=0.0)([5, 5, 5])
    es = ExactSolutionSelector(mesh)("MHD:as3", show_info=True)

    r = es.electric_source_term
    r.current_time = 1
    r.visualize(density=100000)
