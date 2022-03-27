


from numpy import sin, pi
from _3dCSCG.APP.exact_solution.status.Poisson.base import Poisson_Base



# noinspection PyAbstractClass
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

    def phi(self, t, x, y, z):
        return sin(2*pi*x) * sin(2*pi*y) * sin(2*pi*z) + t