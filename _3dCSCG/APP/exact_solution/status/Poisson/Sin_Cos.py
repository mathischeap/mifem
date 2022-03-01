


from numpy import sin, cos, pi
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
        """Return None cause this exact solution is valid at any time."""
        return None

    def phi(self, t, x, y, z):
        return sin(2*pi*x) * sin(2*pi*y) * sin(2*pi*z) + t

    def u(self, t, x, y, z):
        """phi_x"""
        return 2*pi * cos(2*pi*x) * sin(2*pi*y) * sin(2*pi*z)
    def u_x(self, t, x, y, z):
        """phi_x"""
        return -4*pi**2 * sin(2*pi*x) * sin(2*pi*y) * sin(2*pi*z)

    def v(self, t, x, y, z):
        """phi_y"""
        return 2*pi * sin(2*pi*x) * cos(2*pi*y) * sin(2*pi*z)
    def v_y(self, t, x, y, z):
        return -4*pi**2 * sin(2*pi*x) * sin(2*pi*y) * sin(2*pi*z)

    def w(self, t, x, y, z):
        """phi_z"""
        return 2*pi * sin(2*pi*x) * sin(2*pi*y) * cos(2*pi*z)
    def w_z(self, t, x, y, z):
        return -4*pi**2 * sin(2*pi*x) * sin(2*pi*y) * sin(2*pi*z)

    def f(self, t, x, y, z):
        return 12*pi**2 * sin(2*pi*x) * sin(2*pi*y) * sin(2*pi*z)