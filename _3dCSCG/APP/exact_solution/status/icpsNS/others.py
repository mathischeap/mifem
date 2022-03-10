

from numpy import sin, pi
from _3dCSCG.fields.vector.main import _3dCSCG_VectorField
from _3dCSCG.APP.exact_solution.status.icpsNS.base import icpsNS_Base


# noinspection PyAbstractClass
class Closed_Unit_Cube_Disspation1(icpsNS_Base):
    """A modified case that the solution along t is not linear."""
    def __init__(self, es, nu=1, rho=1):
        super(Closed_Unit_Cube_Disspation1, self).__init__(es, nu, rho)
        name = es.mesh.domain.domain_input.__class__.__name__
        bounds = es.mesh.domain.domain_input.bounds
        c = es.mesh.domain.domain_input.c
        assert name == 'Crazy' and \
            bounds[0][0] == 0 and bounds[0][1] == 1 and \
            bounds[1][0] == 0 and bounds[1][1] == 1 and \
            bounds[2][0] == 0 and bounds[2][1] == 1 and \
            c == 0, \
            f"I only work in Crazy domain of [0,1]^3 and c=0. Now it is name={name}, " \
            f"bounds={bounds}, c={c}."

    @property
    def valid_time(self):
        return 'valid_only_at_its_first_instant'

    def u(self, t, x, y, z): return 0 * x

    def u_x(self, t, x, y, z): return 0 * x

    def u_y(self, t, x, y, z): return 0 * x

    def u_z(self, t, x, y, z): return 0 * x

    def v(self, t, x, y, z): return 0 * x

    def v_x(self, t, x, y, z): return 0 * x

    def v_y(self, t, x, y, z): return 0 * x

    def v_z(self, t, x, y, z): return 0 * x

    def w(self, t, x, y, z): return 0 * x

    def w_x(self, t, x, y, z): return 0 * x

    def w_y(self, t, x, y, z): return 0 * x

    def w_z(self, t, x, y, z): return 0 * x


    def fx(self, t, x, y, z): return sin(2 * pi * x)

    def fy(self, t, x, y, z): return sin(2 * pi * y)

    def fz(self, t, x, y, z): return sin(2 * pi * z)

    @property
    def body_force(self):
        """This makes body force valid at all time instants."""
        if self._bodyForce_ is None:
            self._bodyForce_ = _3dCSCG_VectorField(self.mesh, (self.fx, self.fy, self.fz))
        return self._bodyForce_




# noinspection PyAbstractClass
class Constant_X_direction_body_force(icpsNS_Base):
    """"""
    def __init__(self, es, nu=1, rho=1, f=1):
        super(Constant_X_direction_body_force, self).__init__(es, nu, rho)
        self._melt_self_()
        assert isinstance(f, (int, float)), f"f={f} wrong, need to be int or float."
        self._f_ = f
        self._freeze_self_()


    @property
    def valid_time(self):
        return 'valid_only_at_its_first_instant'

    def u(self, t, x, y, z): return 0 * x

    def u_x(self, t, x, y, z): return 0 * x

    def u_y(self, t, x, y, z): return 0 * x

    def u_z(self, t, x, y, z): return 0 * x

    def v(self, t, x, y, z): return 0 * x

    def v_x(self, t, x, y, z): return 0 * x

    def v_y(self, t, x, y, z): return 0 * x

    def v_z(self, t, x, y, z): return 0 * x

    def w(self, t, x, y, z): return 0 * x

    def w_x(self, t, x, y, z): return 0 * x

    def w_y(self, t, x, y, z): return 0 * x

    def w_z(self, t, x, y, z): return 0 * x


    def fx(self, t, x, y, z): return 0 * x + self._f_

    def fy(self, t, x, y, z): return 0 * x

    def fz(self, t, x, y, z): return 0 * x

    @property
    def body_force(self):
        """This makes body force valid at all time instants."""
        if self._bodyForce_ is None:
            self._bodyForce_ = _3dCSCG_VectorField(self.mesh, (self.fx, self.fy, self.fz),
                                                   ftype='constant', valid_time=None)
        return self._bodyForce_



# noinspection PyAbstractClass
class Still(icpsNS_Base):
    """"""
    def __init__(self, es):
        super(Still, self).__init__(es, 1, 1)
        self._melt_self_()
        self._freeze_self_()


    @property
    def valid_time(self):
        return 'valid_only_at_its_first_instant'

    def u(self, t, x, y, z): return 0 * x

    def u_x(self, t, x, y, z): return 0 * x

    def u_y(self, t, x, y, z): return 0 * x

    def u_z(self, t, x, y, z): return 0 * x

    def v(self, t, x, y, z): return 0 * x

    def v_x(self, t, x, y, z): return 0 * x

    def v_y(self, t, x, y, z): return 0 * x

    def v_z(self, t, x, y, z): return 0 * x

    def w(self, t, x, y, z): return 0 * x

    def w_x(self, t, x, y, z): return 0 * x

    def w_y(self, t, x, y, z): return 0 * x

    def w_z(self, t, x, y, z): return 0 * x


    def fx(self, t, x, y, z): return 0 * x

    def fy(self, t, x, y, z): return 0 * x

    def fz(self, t, x, y, z): return 0 * x

    @property
    def body_force(self):
        """This makes body force valid at all time instants."""
        if self._bodyForce_ is None:
            self._bodyForce_ = _3dCSCG_VectorField(self.mesh, (self.fx, self.fy, self.fz),
                                                   ftype='constant', valid_time=None)
        return self._bodyForce_