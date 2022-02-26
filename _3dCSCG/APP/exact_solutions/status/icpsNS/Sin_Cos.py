# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from numpy import sin, cos, pi
from _3dCSCG.APP.exact_solutions.status.icpsNS.base import icpsNS_Base
from _3dCSCG.fields.vector import _3dCSCG_VectorField



# noinspection PyAbstractClass
class SinCosRebholz_Conservation(icpsNS_Base):
    """
    The sin cos test case for the conservation, see Section 5.2 of paper:
        [An Energy- and helicity-conserving finite element scheme for the Navier-Stokes
         equations, Leo G. Rebholz, 2007]
    """
    def __init__(self, es):
        super(SinCosRebholz_Conservation, self).__init__(es, 0, 1)

    @property
    def valid_time(self):
        return 'valid_only_at_its_first_instant'

    def u(self, t, x, y, z): return cos(2 * pi * z)

    def u_x(self, t, x, y, z): return 0 * x

    def u_y(self, t, x, y, z): return 0 * x

    def u_z(self, t, x, y, z): return -2 * pi * sin(2 * pi * z)

    def v(self, t, x, y, z): return sin(2 * pi * z)

    def v_x(self, t, x, y, z): return 0 * x

    def v_y(self, t, x, y, z): return 0 * x

    def v_z(self, t, x, y, z): return 2 * pi * cos(2 * pi * z)

    def w(self, t, x, y, z): return sin(2 * pi * x)

    def w_x(self, t, x, y, z): return 2 * pi * cos(2 * pi * x)

    def w_y(self, t, x, y, z): return 0 * x

    def w_z(self, t, x, y, z): return 0 * x


    def fx(self, t, x, y, z): return 0 * x # can not name it by _fx_

    def fy(self, t, x, y, z): return 0 * x # can not name it by _fy_

    def fz(self, t, x, y, z): return 0 * x # can not name it by _fz_


    @property
    def body_force(self):
        """This makes body force valid at all time instants."""
        if self._bodyForce_ is None:
            self._bodyForce_ = _3dCSCG_VectorField(self.mesh, (self.fx, self.fy, self.fz))
        return self._bodyForce_




class SinCosRebholz_Dissipation(icpsNS_Base):
    """
    The sin cos test case for the conservation, see Section 5.3 of paper:
        [An Energy- and helicity-conserving finite element scheme for the Navier-Stokes
         equations, Leo G. Rebholz, 2007]
    """
    def __init__(self, es, nu=1, rho=1):
        super(SinCosRebholz_Dissipation, self).__init__(es, nu, rho)

    def u(self, t, x, y, z): return (2 - t) * cos(2 * pi * z)

    def u_x(self, t, x, y, z): return 0 * x

    def u_y(self, t, x, y, z): return 0 * x

    def u_z(self, t, x, y, z): return - 2 * pi * (2 - t) * sin(2 * pi * z)

    def u_t(self, t, x, y, z): return - cos(2 * pi * z)

    def u_xx(self, t, x, y, z): return 0 * x

    def u_yy(self, t, x, y, z): return 0 * y

    def u_zz(self, t, x, y, z): return -4 * pi ** 2 * (2 - t) * cos(2 * pi * z)

    def v(self, t, x, y, z): return (1 + t) * sin(2 * pi * z)

    def v_x(self, t, x, y, z): return 0 * x

    def v_y(self, t, x, y, z): return 0 * x

    def v_z(self, t, x, y, z): return 2 * pi * (1 + t) * cos(2 * pi * z)

    def v_t(self, t, x, y, z): return sin(2 * pi * z)

    def v_xx(self, t, x, y, z): return 0 * x

    def v_yy(self, t, x, y, z): return 0 * x

    def v_zz(self, t, x, y, z): return - 4 * pi ** 2 * (1 + t) * sin(2 * pi * z)

    def w(self, t, x, y, z): return (1 - t) * sin(2 * pi * x)

    def w_x(self, t, x, y, z): return 2 * pi * (1 - t) * cos(2 * pi * x)

    def w_y(self, t, x, y, z): return 0 * x

    def w_z(self, t, x, y, z): return 0 * x

    def w_t(self, t, x, y, z): return - sin(2 * pi * x)

    def w_xx(self, t, x, y, z): return - 4 * pi ** 2 * (1 - t) * sin(2 * pi * x)

    def w_yy(self, t, x, y, z): return 0 * x

    def w_zz(self, t, x, y, z): return 0 * x

    def p(self, t, x, y, z): return sin(2 * pi * (x + y + z + t))

    def p_x(self, t, x, y, z): return 2 * pi * cos(2 * pi * (x + y + z + t))

    def p_y(self, t, x, y, z): return 2 * pi * cos(2 * pi * (x + y + z + t))

    def p_z(self, t, x, y, z): return 2 * pi * cos(2 * pi * (x + y + z + t))




class SinCos_Modified_Dissipation(icpsNS_Base):
    """A modified case that the solution along t is not linear."""
    def __init__(self, es, nu=1, rho=1):
        super(SinCos_Modified_Dissipation, self).__init__(es, nu, rho)

    def u(self, t, x, y, z): return (1 - sin(2*pi*t)) * cos(2 * pi * z)

    def u_x(self, t, x, y, z): return 0 * x

    def u_y(self, t, x, y, z): return 0 * x

    def u_z(self, t, x, y, z): return - 2 * pi * (1 - sin(2*pi*t)) * sin(2 * pi * z)

    def u_t(self, t, x, y, z): return - 2*pi*cos(2*pi*t) * cos(2 * pi * z)

    def u_xx(self, t, x, y, z): return 0 * x

    def u_yy(self, t, x, y, z): return 0 * y

    def u_zz(self, t, x, y, z): return -4 * pi ** 2 * (1 - sin(2*pi*t)) * cos(2 * pi * z)

    def v(self, t, x, y, z): return (1 + cos(2*pi*t)) * sin(2 * pi * z)

    def v_x(self, t, x, y, z): return 0 * x

    def v_y(self, t, x, y, z): return 0 * x

    def v_z(self, t, x, y, z): return 2 * pi * (1 + cos(2*pi*t)) * cos(2 * pi * z)

    def v_t(self, t, x, y, z): return -2*pi*sin(2*pi*t) * sin(2 * pi * z)

    def v_xx(self, t, x, y, z): return 0 * x

    def v_yy(self, t, x, y, z): return 0 * x

    def v_zz(self, t, x, y, z): return - 4 * pi ** 2 * (1 + cos(2*pi*t)) * sin(2 * pi * z)

    def w(self, t, x, y, z): return (1 - sin(2*pi*t)) * sin(2 * pi * x)

    def w_x(self, t, x, y, z): return 2 * pi * (1 - sin(2*pi*t)) * cos(2 * pi * x)

    def w_y(self, t, x, y, z): return 0 * x

    def w_z(self, t, x, y, z): return 0 * x

    def w_t(self, t, x, y, z): return - 2*pi*cos(2*pi*t) * sin(2 * pi * x)

    def w_xx(self, t, x, y, z): return - 4 * pi ** 2 * (1 - sin(2*pi*t)) * sin(2 * pi * x)

    def w_yy(self, t, x, y, z): return 0 * x

    def w_zz(self, t, x, y, z): return 0 * x

    def p(self, t, x, y, z): return sin(2 * pi * (x + y + z + t))

    def p_x(self, t, x, y, z): return 2 * pi * cos(2 * pi * (x + y + z + t))

    def p_y(self, t, x, y, z): return 2 * pi * cos(2 * pi * (x + y + z + t))

    def p_z(self, t, x, y, z): return 2 * pi * cos(2 * pi * (x + y + z + t))



# noinspection PyAbstractClass
class SinCos_Conservation_Conservative_Body_Force(icpsNS_Base):
    """
    The sin cos test case for the conservation, see Section 5.2 of paper:
        [An Energy- and helicity-conserving finite element scheme for the Navier-Stokes
         equations, Leo G. Rebholz, 2007]
    """
    def __init__(self, es):
        super(SinCos_Conservation_Conservative_Body_Force, self).__init__(es, 0, 1)

    @property
    def valid_time(self):
        return 'valid_only_at_its_first_instant'

    def u(self, t, x, y, z): return cos(2 * pi * z)

    def u_x(self, t, x, y, z): return 0 * x

    def u_y(self, t, x, y, z): return 0 * x

    def u_z(self, t, x, y, z): return -2 * pi * sin(2 * pi * z)

    def v(self, t, x, y, z): return sin(2 * pi * z)

    def v_x(self, t, x, y, z): return 0 * x

    def v_y(self, t, x, y, z): return 0 * x

    def v_z(self, t, x, y, z): return 2 * pi * cos(2 * pi * z)

    def w(self, t, x, y, z): return sin(2 * pi * x)

    def w_x(self, t, x, y, z): return 2 * pi * cos(2 * pi * x)

    def w_y(self, t, x, y, z): return 0 * x

    def w_z(self, t, x, y, z): return 0 * x


    # varphi(t,x,y,z) = t * sin(2 * pi * x) * sin(2 * pi * y) * sin(2 * pi * z)

    def fx(self, t, x, y, z): return 2 * pi * t * cos(2 * pi * x) * sin(2 * pi * y) * sin(2 * pi * z)

    def fy(self, t, x, y, z): return 2 * pi * t * sin(2 * pi * x) * cos(2 * pi * y) * sin(2 * pi * z)

    def fz(self, t, x, y, z): return 2 * pi * t * sin(2 * pi * x) * sin(2 * pi * y) * cos(2 * pi * z)


    @property
    def body_force(self):
        """This makes body force valid at all time instants."""
        if self._bodyForce_ is None:
            self._bodyForce_ = _3dCSCG_VectorField(self.mesh, (self.fx, self.fy, self.fz))
        return self._bodyForce_

# noinspection PyAbstractClass
class SinCos_Conservation_Conservative_Body_Force1(icpsNS_Base):
    """
    The sin cos test case for the conservation, see Section 5.2 of paper:
        [An Energy- and helicity-conserving finite element scheme for the Navier-Stokes
         equations, Leo G. Rebholz, 2007]
    """
    def __init__(self, es):
        super(SinCos_Conservation_Conservative_Body_Force1, self).__init__(es, 0, 1)

    @property
    def valid_time(self):
        return 'valid_only_at_its_first_instant'

    def u(self, t, x, y, z): return cos(2 * pi * z)

    def u_x(self, t, x, y, z): return 0 * x

    def u_y(self, t, x, y, z): return 0 * x

    def u_z(self, t, x, y, z): return -2 * pi * sin(2 * pi * z)

    def v(self, t, x, y, z): return sin(2 * pi * z)

    def v_x(self, t, x, y, z): return 0 * x

    def v_y(self, t, x, y, z): return 0 * x

    def v_z(self, t, x, y, z): return 2 * pi * cos(2 * pi * z)

    def w(self, t, x, y, z): return sin(2 * pi * x)

    def w_x(self, t, x, y, z): return 2 * pi * cos(2 * pi * x)

    def w_y(self, t, x, y, z): return 0 * x

    def w_z(self, t, x, y, z): return 0 * x


    # varphi(t,x,y,z) = sin(2 * pi * x) * sin(2 * pi * y) * sin(2 * pi * z)

    def fx(self, t, x, y, z): return 2 * pi * cos(2 * pi * x) * sin(2 * pi * y) * sin(2 * pi * z)

    def fy(self, t, x, y, z): return 2 * pi * sin(2 * pi * x) * cos(2 * pi * y) * sin(2 * pi * z)

    def fz(self, t, x, y, z): return 2 * pi * sin(2 * pi * x) * sin(2 * pi * y) * cos(2 * pi * z)


    @property
    def body_force(self):
        """This makes body force valid at all time instants."""
        if self._bodyForce_ is None:
            self._bodyForce_ = _3dCSCG_VectorField(self.mesh, (self.fx, self.fy, self.fz))
        return self._bodyForce_

# noinspection PyAbstractClass
class SinCos_Conservation_Conservative_Body_Force_POLYNOMIALS(icpsNS_Base):
    """
    The sin cos test case for the conservation, see Section 5.2 of paper:
        [An Energy- and helicity-conserving finite element scheme for the Navier-Stokes
         equations, Leo G. Rebholz, 2007]
    """
    def __init__(self, es):
        super(SinCos_Conservation_Conservative_Body_Force_POLYNOMIALS, self).__init__(es, 0, 1)

    @property
    def valid_time(self):
        return 'valid_only_at_its_first_instant'

    def u(self, t, x, y, z): return cos(2 * pi * z)

    def u_x(self, t, x, y, z): return 0 * x

    def u_y(self, t, x, y, z): return 0 * x

    def u_z(self, t, x, y, z): return -2 * pi * sin(2 * pi * z)

    def v(self, t, x, y, z): return sin(2 * pi * z)

    def v_x(self, t, x, y, z): return 0 * x

    def v_y(self, t, x, y, z): return 0 * x

    def v_z(self, t, x, y, z): return 2 * pi * cos(2 * pi * z)

    def w(self, t, x, y, z): return sin(2 * pi * x)

    def w_x(self, t, x, y, z): return 2 * pi * cos(2 * pi * x)

    def w_y(self, t, x, y, z): return 0 * x

    def w_z(self, t, x, y, z): return 0 * x


    # phi(t,x,y,z) = t * (x**3/3 - x**2/2 + y**3/3 - y**2/2 + z**3/3 - z**2/2)

    def fx(self, t, x, y, z): return t * x * (x-1)

    def fy(self, t, x, y, z): return t * y * (y-1)

    def fz(self, t, x, y, z): return t * z * (z-1)


    @property
    def body_force(self):
        """This makes body force valid at all time instants."""
        if self._bodyForce_ is None:
            self._bodyForce_ = _3dCSCG_VectorField(self.mesh, (self.fx, self.fy, self.fz))
        return self._bodyForce_




# noinspection PyAbstractClass
class SinCos_Conservation_Conservative_Body_Force_CONSTANT(icpsNS_Base):
    """
    The sin cos test case for the conservation, see Section 5.2 of paper:
        [An Energy- and helicity-conserving finite element scheme for the Navier-Stokes
         equations, Leo G. Rebholz, 2007]
    """
    def __init__(self, es):
        super(SinCos_Conservation_Conservative_Body_Force_CONSTANT, self).__init__(es, 0, 1)

    @property
    def valid_time(self):
        return 'valid_only_at_its_first_instant'

    def u(self, t, x, y, z): return cos(2 * pi * z)

    def u_x(self, t, x, y, z): return 0 * x

    def u_y(self, t, x, y, z): return 0 * x

    def u_z(self, t, x, y, z): return -2 * pi * sin(2 * pi * z)

    def v(self, t, x, y, z): return sin(2 * pi * z)

    def v_x(self, t, x, y, z): return 0 * x

    def v_y(self, t, x, y, z): return 0 * x

    def v_z(self, t, x, y, z): return 2 * pi * cos(2 * pi * z)

    def w(self, t, x, y, z): return sin(2 * pi * x)

    def w_x(self, t, x, y, z): return 2 * pi * cos(2 * pi * x)

    def w_y(self, t, x, y, z): return 0 * x

    def w_z(self, t, x, y, z): return 0 * x


    # phi(t,x,y,z) = x

    def fx(self, t, x, y, z): return 1 + 0 * x * y * z

    def fy(self, t, x, y, z): return 0 + 0 * x * y * z

    def fz(self, t, x, y, z): return 0 + 0 * x * y * z


    @property
    def body_force(self):
        """This makes body force valid at all time instants."""
        if self._bodyForce_ is None:
            self._bodyForce_ = _3dCSCG_VectorField(self.mesh, (self.fx, self.fy, self.fz))
        return self._bodyForce_
