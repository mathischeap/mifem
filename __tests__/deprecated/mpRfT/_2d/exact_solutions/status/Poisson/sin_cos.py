# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/17 3:51 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')


from numpy import sin, pi
from objects.mpRfT._2d.exact_solutions.status.Poisson.base import Poisson_Base



class Poisson_SinCos1(Poisson_Base):
    """
    The sin cos test case 1.
    """
    def __init__(self, es):
        super(Poisson_SinCos1, self).__init__(es)
        self._es_.standard_properties.name = 'Poisson-sin-cos-1'


    @property
    def valid_time(self):
        """Return None because this exact solution is valid at any time."""
        return None

    def phi(self, t, x, y):
        return sin(pi*x) * sin(pi*y) + t


if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/exact_solutions/status/Poisson/sin_cos.py
    pass
