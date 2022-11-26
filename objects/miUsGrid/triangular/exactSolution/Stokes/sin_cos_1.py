# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 10/10/2022 11:01 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from numpy import sin, pi, cos
from objects.miUsGrid.triangular.exactSolution.Stokes.base import Stokes



# noinspection PyAbstractClass
class Stokes_SinCos1(Stokes):
    """The sin-cos test case 1.
    """
    def __init__(self, mesh):
        super(Stokes_SinCos1, self).__init__(mesh)



    @property
    def valid_time(self):
        """Return None because this exact solution is valid at any time."""
        return None

    def p(self, t, x, y):
        return sin(2*pi*x) * sin(2*pi*y)


    def u(self, t, x, y):
        return cos(2*pi*x) * sin(2*pi*y)
    def v(self, t, x, y):
        return - sin(2*pi*x) * cos(2*pi*y)



if __name__ == '__main__':
    # mpiexec -n 4 python objects/miUsGrid/triangular/exact_solution/Stokes/sin_cos_1.py
    from __init__ import miTri

    fc = miTri.call('rand0', 2)

    es = fc('Stokes : sin cos 1')

