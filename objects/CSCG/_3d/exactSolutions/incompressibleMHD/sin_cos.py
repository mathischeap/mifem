# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/20/2022 10:33 PM
"""
import sys

if './' not in sys.path:
    sys.path.append('./')
from objects.CSCG._3d.exactSolutions.incompressibleMHD.base import incompressible_MHD_Base

from numpy import sin, cos, pi


class MHD_SinCos1(incompressible_MHD_Base):
    """The first manufactured solution with sin cos functions."""
    def __init__(self, mesh, Rf=1, Rm=1, c=1):
        """"""
        super(MHD_SinCos1, self).__init__(mesh, Rf, Rm, c)
        self._freeze_self_()

    def u(self, t, x, y, z):
        """"""
        return cos(2 * pi * z)

    def v(self, t, x, y, z):
        """"""
        return sin(2 * pi * z)

    def w(self, t, x, y, z):
        """"""
        return sin(2 * pi * x)

    def p(self, t, x, y, z):
        """"""
        return sin(2 * pi * (x + y + z)) + t

    def Bx(self, t, x, y, z):
        """"""
        return cos(2 * pi * z)

    def By(self, t, x, y, z):
        """"""
        return sin(2 * pi * z)

    def Bz(self, t, x, y, z):
        """"""
        return sin(2 * pi * x)


if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_3d/exact_solutions/status/incompressible_MHD/sin_cos.py
    from objects.CSCG._3d.master import MeshGenerator, ExactSolutionSelector
    mesh = MeshGenerator('crazy', c=0.0)([5, 5, 5])
    es = ExactSolutionSelector(mesh)("MHD:sincos1", show_info=True)

    r = es.status.velocity
    r.current_time = 1
    r.visualize()
