# -*- coding: utf-8 -*-
"""The 3d Poisson problem, k = 1.

The Poisson equation: (k=1)
- div (k grad phi) = f

And in the mixed formulation:

u = k grad phi
f = - div u

And we set k equal to 1; k = 1.

We will interpolate the Poisson problem as a flow problem (potential flow). So phi is the potential,
(u,v,w) is the velocity, and f is the source term. Note that they are just names, they can also be
for example (temperature, heat flux and heat source) and so on.

"""

import numpy as np
from functools import lru_cache
from objects.CSCG._3d.exact_solutions.status.base import Base
from objects.CSCG._3d.fields.vector.main import _3dCSCG_VectorField
from objects.CSCG._3d.fields.scalar.main import _3dCSCG_ScalarField

from screws.numerical.time_plus_3d_space.partial_derivative import NumericalPartialDerivative_txyz
from screws.numerical.time_plus_3d_space.partial_derivative_as_functions import \
    NumericalPartialDerivative_txyz_Functions


class Poisson_Base(Base):
    def __init__(self, es):
        super(Poisson_Base, self).__init__(es)
        self._potential_ = None
        self._velocity_ = None
        self._source_term_ = None
        self._kineticEnergyDistribution_ = None

        self._NPDf_p_ = None
        self._NPDf_px_ = None
        self._NPDf_py_ = None
        self._NPDf_pz_ = None

        self._freeze_self_()


    # to be overridden (must) ...

    def phi(self, t, x, y, z):
        raise NotImplementedError()

    def u(self, t, x, y, z):
        """phi_x"""
        if self._NPDf_p_ is None:
            self._NPDf_p_ = NumericalPartialDerivative_txyz_Functions(self.phi)
        return self._NPDf_p_('x')(t, x, y, z)
    def v(self, t, x, y, z):
        """phi_y"""
        if self._NPDf_p_ is None:
            self._NPDf_p_ = NumericalPartialDerivative_txyz_Functions(self.phi)
        return self._NPDf_p_('y')(t, x, y, z)
    def w(self, t, x, y, z):
        """phi_z"""
        if self._NPDf_p_ is None:
            self._NPDf_p_ = NumericalPartialDerivative_txyz_Functions(self.phi)
        return self._NPDf_p_('z')(t, x, y, z)


    def u_x(self, t, x, y, z):
        if self._NPDf_px_ is None:
            self._NPDf_px_ = NumericalPartialDerivative_txyz_Functions(self.u)
        return self._NPDf_px_('x')(t, x, y, z)

    def v_y(self, t, x, y, z):
        if self._NPDf_py_ is None:
            self._NPDf_py_ = NumericalPartialDerivative_txyz_Functions(self.v)
        return self._NPDf_py_('y')(t, x, y, z)

    def w_z(self, t, x, y, z):
        if self._NPDf_pz_ is None:
            self._NPDf_pz_ = NumericalPartialDerivative_txyz_Functions(self.w)
        return self._NPDf_pz_('z')(t, x, y, z)


    def f(self, t, x, y, z):
        return - self.u_x(t, x, y, z) - self.v_y(t, x, y, z) - self.w_z(t, x, y, z)


    @property
    def potential(self):
        if self._potential_ is None:
            self._potential_ = _3dCSCG_ScalarField(
                self.mesh,
                self.phi,
                valid_time=self.valid_time,
                name='potential')
        return self._potential_


    @property
    def velocity(self):
        if self._velocity_ is None:
            self._velocity_ = _3dCSCG_VectorField(
                self.mesh,
                (self.u, self.v, self.w),
                valid_time=self.valid_time,
                name='velocity')
        return self._velocity_

    @property
    def source_term(self):
        if self._source_term_ is None:
            self._source_term_ = _3dCSCG_ScalarField(
                 self.mesh,
                 self.f,
                 valid_time=self.valid_time,
                 name='source_term')
        return self._source_term_



    @property
    def kinetic_energy_distribution(self):
        """A scalar field of the kinetic energy distribution."""
        if self._kineticEnergyDistribution_ is None:
            self._kineticEnergyDistribution_ =_3dCSCG_ScalarField(
                self.mesh,
                self.___kinetic_energy_distribution___,
                valid_time=self.valid_time,
                name='kinetic_energy_distribution')
        return self._kineticEnergyDistribution_
    def ___kinetic_energy_distribution___(self, t, x, y, z):
        return 0.5 * (self.u(t, x, y, z)**2 + self.v(t, x, y, z)**2 + self.w(t, x, y, z)**2)
    @lru_cache(maxsize=8)
    def kinetic_energy(self, t):
        """Kinetic energy at time `t`."""
        return self._es_.do.compute_Ln_norm_of('kinetic_energy_distribution', time=t, n=1)


    def ___PreFrozenChecker___(self):
        """
        We use this general method to do the check, in particular exact solution, we can define particular
        check method by override this method.
        """
        TS = self.___PRIVATE_generate_random_valid_time_instances___()
        x, y, z = self._mesh_.do.generate_random_coordinates()

        if len(x) == 0: return

        for t in TS:

            t = float(t)

            try:
                Pu = NumericalPartialDerivative_txyz(self.phi, t, x, y, z)
                assert Pu.check_partial_x(self.u)
                assert Pu.check_partial_y(self.v)
                assert Pu.check_partial_z(self.w)
            except NotImplementedError:
                pass

            try:
                Pu = NumericalPartialDerivative_txyz(self.u, t, x, y, z)
                assert Pu.check_partial_x(self.u_x)
            except NotImplementedError:
                pass

            try:
                Pu = NumericalPartialDerivative_txyz(self.v, t, x, y, z)
                assert Pu.check_partial_y(self.v_y)
            except NotImplementedError:
                pass

            try:
                Pu = NumericalPartialDerivative_txyz(self.w, t, x, y, z)
                assert Pu.check_partial_z(self.w_z)
            except NotImplementedError:
                pass

            try:
                f = - self.u_x(t, x, y, z) - self.v_y(t, x, y, z) - self.w_z(t, x, y, z)
                F = self.f(t, x, y, z)
                np.testing.assert_array_almost_equal(F - f, 0, decimal=5)
            except NotImplementedError:
                pass