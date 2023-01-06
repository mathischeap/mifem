# -*- coding: utf-8 -*-
from objects.CSCG._2d.exactSolutions.base import Base
from components.numerical.timePlus2dSpace.partial_derivative_as_functions import \
    NumericalPartialDerivative_txy_Functions
from components.numerical.timePlus2dSpace.partial_derivative import NumericalPartialDerivative_txy

from objects.CSCG._2d.fields.vector.main import _2dCSCG_VectorField
from objects.CSCG._2d.fields.scalar.main import _2dCSCG_ScalarField

import numpy as np


class incompressibleNavierStokesBase(Base):
    """ """
    def __init__(self, mesh, nu):
        """

        :param mesh:
        :param nu: The kinematic viscosity nu = mu / rho
        """
        super(incompressibleNavierStokesBase, self).__init__(mesh)
        assert isinstance(nu, (int, float)) and nu >= 0, f'kinematic viscosity  = {nu} is wrong.'
        self._nu_ = nu

        self._velocity_ = None
        self._vorticity_ = None
        self._body_force_ = None
        self._pressure_ = None
        self._gradient_of_pressure_ = None
        self._total_pressure_ = None
        self._kinetic_energy_distribution_ = None
        self._curl_of_omega_ = None
        self._divergence_of_velocity_ = None

        self._NPDf_u_ = None
        self._NPDf_v_ = None
        self._NPDf_p_ = None
        self._NPDf_tp_ = None
        self._NPDf_omega_ = None
        self._freeze_self_()

    @property
    def nu(self):
        """The kinematic viscosity."""
        return self._nu_

    @property
    def Re(self):
        """Reynolds number. It will also depend on the computational domain and velocity."""
        raise NotImplementedError()

    def u(self, t, x, y):  # u, x-component of velocity
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

    def v(self, t, x, y):  # v, y-component of velocity
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
            self._velocity_ = _2dCSCG_VectorField(
                self.mesh, [self.u, self.v], valid_time=self.valid_time, name='velocity')
        return self._velocity_

    @property
    def divergence_of_velocity(self):
        if self._divergence_of_velocity_ is None:
            # this condition must be valid for all time.
            self._divergence_of_velocity_ = _2dCSCG_ScalarField(self.mesh,
                                                                self.___div_of_velocity___,
                                                                valid_time=None)
        return self._divergence_of_velocity_

    def ___div_of_velocity___(self, t, x, y):
        # must be zero
        return self.u_x(t, x, y) + self.v_y(t, x, y)


    def omega(self, t, x, y):
        return self.v_x(t, x, y) - self.u_y(t, x, y)

    @property
    def vorticity(self):
        """A scalar field of the kinetic energy distribution."""
        if self._vorticity_ is None:
            self._vorticity_ = _2dCSCG_ScalarField(
                self.mesh, self.omega, valid_time=self.valid_time, name='vorticity')
        return self._vorticity_

    def omega_t(self, t, x, y):
        if self._NPDf_omega_ is None:
            self._NPDf_omega_ = NumericalPartialDerivative_txy_Functions(self.omega)
        return self._NPDf_omega_('t')(t, x, y)

    def omega_x(self, t, x, y):
        if self._NPDf_omega_ is None:
            self._NPDf_omega_ = NumericalPartialDerivative_txy_Functions(self.omega)
        return self._NPDf_omega_('x')(t, x, y)

    def omega_y(self, t, x, y):
        if self._NPDf_omega_ is None:
            self._NPDf_omega_ = NumericalPartialDerivative_txy_Functions(self.omega)
        return self._NPDf_omega_('y')(t, x, y)

    def ___minus_omega_x___(self, t, x, y):
        return - self.omega_x(t, x, y)

    @property
    def curl_of_omega(self):
        if self._curl_of_omega_ is None:
            self._curl_of_omega_ = _2dCSCG_VectorField(
                self.mesh,
                [self.omega_y, self.___minus_omega_x___],
                valid_time=self.valid_time,
                name='curl_of_omega')
        return self._curl_of_omega_

    @property
    def curl_of_vorticity(self):
        return self.curl_of_omega

    def p(self, t, x, y):  # static pressure
        raise NotImplementedError()

    @property
    def pressure(self):
        """A scalar field of the kinetic energy distribution."""
        if self._pressure_ is None:
            self._pressure_ = _2dCSCG_ScalarField(
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

    @property
    def total_pressure(self):
        """A scalar field of the kinetic energy distribution."""
        if self._total_pressure_ is None:
            self._total_pressure_ = _2dCSCG_ScalarField(
                self.mesh, self._tp_, valid_time=self.valid_time, name='total pressure')
        return self._total_pressure_

    def _tp_(self, t, x, y):  # total pressure
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
            self._kinetic_energy_distribution_ = _2dCSCG_ScalarField(
                self.mesh, self.___kinetic_energy_distribution___,
                valid_time=self.valid_time,
                name='kinetic_energy')
        return self._kinetic_energy_distribution_

    def fx(self, t, x, y):
        return self.u_t(t, x, y) - self.omega(t, x, y) * self.v(t, x, y)  \
               + self.nu * self.omega_y(t, x, y) + self._tp_x_(t, x, y)

    def fy(self, t, x, y):
        return self.v_t(t, x, y) + self.omega(t, x, y) * self.u(t, x, y)  \
               - self.nu * self.omega_x(t, x, y) + self._tp_y_(t, x, y)

    @property
    def body_force(self):
        """A scalar field of the kinetic energy distribution."""
        if self._body_force_ is None:
            self._body_force_ = _2dCSCG_VectorField(
                self.mesh, [self.fx, self.fy], valid_time=self.valid_time, name='body_force')
        return self._body_force_

    # noinspection PyUnusedLocal
    @staticmethod
    def source_term(t, x, y):
        """By default, we have divergence free condition; the source term is zero."""
        return 0 * x


    def ___PreFrozenChecker___(self):
        """Will be called before freezing self."""
        TS = self.___PRIVATE_generate_random_valid_time_instances___()
        x, y = self._mesh_.do.generate_random_coordinates()

        if len(x) == 0:
            return

        for t in TS:
            try:
                fx = self.u_t(t, x, y) - self.omega(t, x, y) * self.v(t, x, y) \
                    + self.nu * self.omega_y(t, x, y) + self._tp_x_(t, x, y)
                fy = self.v_t(t, x, y) + self.omega(t, x, y) * self.u(t, x, y) \
                    - self.nu * self.omega_x(t, x, y) + self._tp_y_(t, x, y)
                FX = self.fx(t, x, y)
                FY = self.fy(t, x, y)
                np.testing.assert_array_almost_equal(FX - fx, 0, decimal=5)
                np.testing.assert_array_almost_equal(FY - fy, 0, decimal=5)
            except NotImplementedError:
                pass

            try:
                Pu = NumericalPartialDerivative_txy(self.u, t, x, y)
                assert all(Pu.check_total(self.u_t, self.u_x, self.u_y))
            except NotImplementedError:
                pass

            try:
                Pu = NumericalPartialDerivative_txy(self.v, t, x, y)
                assert all(Pu.check_total(self.v_t, self.v_x, self.v_y))
            except NotImplementedError:
                pass

            try:
                Pu = NumericalPartialDerivative_txy(self.p, t, x, y)
                assert all(Pu.check_total(self.p_t, self.p_x, self.p_y))
            except NotImplementedError:
                pass

            try:
                np.testing.assert_array_almost_equal(
                    self.u_x(t, x, y) + self.v_y(t, x, y),
                    self.source_term(t, x, y),
                    decimal=5)
            except NotImplementedError:
                pass
