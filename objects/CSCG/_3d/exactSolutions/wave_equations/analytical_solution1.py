# -*- coding: utf-8 -*-
"""
Yi Zhang
zhangyi_aero@hotmail.com
created at: 3/5/2023 2:45 PM
"""
import sys

if './' not in sys.path:
    sys.path.append('./')
from objects.CSCG._3d.exactSolutions.wave_equations.base import Wave_Base
from numpy import pi, sin, exp


class AS1(Wave_Base):
    """The first analytical solution."""
    def __init__(self, mesh):
        """"""
        super(AS1, self).__init__(mesh)
        self._freeze_self_()

    def phi(self, t, x, y, z):
        return sin(2*pi*x) * sin(2*pi*y) * sin(2*pi*z) * exp(t/pi)


if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_3d/exactSolutions/wave_equations/analytical_solution1.py

    from objects.CSCG._3d.master import MeshGenerator, ExactSolutionSelector
    mesh = MeshGenerator('crazy', bounds=([0, 1], [0, 1], [0, 1]), c=0.0)([5, 5, 5])
    es = ExactSolutionSelector(mesh)("WE:as1", show_info=True)

    r = es.velocity
    r.current_time = 0
    r.visualize(density=100000)
