"""

The 2d Euler flows (incompressible, constant density, inviscid Navier-Stokes flows).

"""


from _2dCSCG.APP.exact_solution.status.base import Base
from screws.numerical.time_plus_2d_space.partial_derivative_as_functions import NumericalPartialDerivative_txy_Functions

from _2dCSCG.fields.vector.main import _2dCSCG_VectorField
from _2dCSCG.fields.scalar.main import _2dCSCG_ScalarField

class EulerBase(Base):
    """ Given a scalar field: potential,  source = laplace potential.

    or in mixed formulation:

        velocity = gradient potential
        source = divergence velocity

    """
    def __init__(self, es):
        super(EulerBase, self).__init__(es)

        self._velocity_ = None
        self._vorticity_ = None
        self._body_force_ = None
        self._pressure_ = None
        self._gradient_of_pressure_ = None
        self._total_pressure_ = None
        self._kinetic_energy_distribution_ = None

        self._NPDf_u_ = None
        self._NPDf_v_ = None
        self._NPDf_p_ = None
        self._NPDf_tp_ = None

        self._freeze_self_()

    # to be overridden (must) ...

    def u(self, t, x, y): # u-component of velocity
        raise NotImplementedError()

    def u_t(self, t, x, y):
        if self._NPDf_u_ is None:
            self._NPDf_u_ = NumericalPartialDerivative_txy_Functions(self.u)
        return self._NPDf_u_('t')(t, x, y)
    def u_x(self, t, x, y):
        if self._NPDf_u_ is None:
            self._NPDf_u_ = NumericalPartialDerivative_txy_Functions(self.u)
        return self._NPDf_u_('x')(t, x, y)
    def u_y(self, t, x, y):
        if self._NPDf_u_ is None:
            self._NPDf_u_ = NumericalPartialDerivative_txy_Functions(self.u)
        return self._NPDf_u_('y')(t, x, y)

    def v(self, t, x, y): # v-component of velocity
        raise NotImplementedError()
    def v_t(self, t, x, y):
        if self._NPDf_v_ is None:
            self._NPDf_v_ = NumericalPartialDerivative_txy_Functions(self.v)
        return self._NPDf_v_('t')(t, x, y)
    def v_x(self, t, x, y):
        if self._NPDf_v_ is None:
            self._NPDf_v_ = NumericalPartialDerivative_txy_Functions(self.v)
        return self._NPDf_v_('x')(t, x, y)
    def v_y(self, t, x, y):
        if self._NPDf_v_ is None:
            self._NPDf_v_ = NumericalPartialDerivative_txy_Functions(self.v)
        return self._NPDf_v_('y')(t, x, y)

    @property
    def velocity(self):
        """A scalar field of the kinetic energy distribution."""
        if self._velocity_ is None:
            self._velocity_ =_2dCSCG_VectorField(
                self.mesh, [self.u, self.v], valid_time=self.valid_time, name='velocity')
        return self._velocity_

    def omega(self, t, x, y):
        return self.v_x(t, x, y) - self.u_y(t, x, y)

    @property
    def vorticity(self):
        """A scalar field of the kinetic energy distribution."""
        if self._vorticity_ is None:
            self._vorticity_ =_2dCSCG_ScalarField(
                self.mesh, self.omega, valid_time=self.valid_time, name='vorticity')
        return self._vorticity_

    def p(self, t, x, y): # static pressure
        raise NotImplementedError()
    @property
    def pressure(self):
        """A scalar field of the kinetic energy distribution."""
        if self._pressure_ is None:
            self._pressure_ =_2dCSCG_ScalarField(
                self.mesh, self.p, valid_time=self.valid_time, name='pressure')
        return self._pressure_
    def p_t(self, t, x, y):
        if self._NPDf_p_ is None:
            self._NPDf_p_ = NumericalPartialDerivative_txy_Functions(self.p)
        return self._NPDf_p_('t')(t, x, y)
    def p_x(self, t, x, y):
        if self._NPDf_p_ is None:
            self._NPDf_p_ = NumericalPartialDerivative_txy_Functions(self.p)
        return self._NPDf_p_('x')(t, x, y)
    def p_y(self, t, x, y):
        if self._NPDf_p_ is None:
            self._NPDf_p_ = NumericalPartialDerivative_txy_Functions(self.p)
        return self._NPDf_p_('y')(t, x, y)


    def _tp_(self, t, x, y): #total pressure
        return self.p(t, x, y) + self.___kinetic_energy_distribution___(t, x, y)
    def _tp_x_(self, t, x, y):
        if self._NPDf_tp_ is None:
            self._NPDf_tp_ = NumericalPartialDerivative_txy_Functions(self._tp_)
        return self._NPDf_tp_('x')(t, x, y)
    def _tp_y_(self, t, x, y):
        if self._NPDf_tp_ is None:
            self._NPDf_tp_ = NumericalPartialDerivative_txy_Functions(self._tp_)
        return self._NPDf_tp_('y')(t, x, y)


    def ___kinetic_energy_distribution___(self, t, x, y):
        return 0.5 * (self.u(t, x, y)**2 + self.v(t, x, y)**2)

    @property
    def kinetic_energy_distribution(self):
        """A scalar field of the kinetic energy distribution."""
        if self._kinetic_energy_distribution_ is None:
            self._kinetic_energy_distribution_ =_2dCSCG_ScalarField(
                self.mesh, self.___kinetic_energy_distribution___,
                valid_time=self.valid_time,
                name='kinetic_energy')
        return self._kinetic_energy_distribution_


    def fx(self, t, x, y):
        return self.u_t(t,x,y) - self.omega(t, x, y) * self.v(t, x, y)  + self._tp_x_(t, x, y)
    def fy(self, t, x, y):
        return self.v_t(t,x,y) + self.omega(t, x, y) * self.u(t, x, y)  + self._tp_y_(t, x, y)


    @property
    def body_force(self):
        """A scalar field of the kinetic energy distribution."""
        if self._body_force_ is None:
            self._body_force_ =_2dCSCG_VectorField(
                self.mesh, [self.fx, self.fy], valid_time=self.valid_time, name='body_force')
        return self._body_force_



    def ___PreFrozenChecker___(self):
        """Will be called before freezing self."""
