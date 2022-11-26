# -*- coding: utf-8 -*-
"""

The Poisson equation: (k=1)

- div (k grad phi) = f

And in the mixed formulation:

u = k grad phi
f = - div u

And we set k equal to 1; k = 1.

We will interpolate the Poisson problem as a flow problem (potential flow). So phi is the potential,
(u,v,w) is the velocity, and f is the source term. Note that they are just names, they can also be
for example (temperature, heat flux and heat source) and so on.
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/06/17 3:51 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

import numpy as np
from objects.mpRfT._2d.exact_solutions.status.base import Base
from objects.mpRfT._2d.cf.vector.main import mpRfT2_Vector
from objects.mpRfT._2d.cf.scalar.main import mpRfT2_Scalar

from components.numerical.timePlus2dSpace.partial_derivative import NumericalPartialDerivative_txy
from components.numerical.timePlus2dSpace.partial_derivative_as_functions import \
    NumericalPartialDerivative_txy_Functions


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

    def phi(self, t, x, y):
        raise NotImplementedError()

    def u(self, t, x, y):
        """phi_x"""
        if self._NPDf_p_ is None:
            self._NPDf_p_ = NumericalPartialDerivative_txy_Functions(self.phi)
        return self._NPDf_p_('x')(t, x, y)
    def v(self, t, x, y):
        """phi_y"""
        if self._NPDf_p_ is None:
            self._NPDf_p_ = NumericalPartialDerivative_txy_Functions(self.phi)
        return self._NPDf_p_('y')(t, x, y)



    def u_x(self, t, x, y):
        if self._NPDf_px_ is None:
            self._NPDf_px_ = NumericalPartialDerivative_txy_Functions(self.u)
        return self._NPDf_px_('x')(t, x, y)

    def v_y(self, t, x, y):
        if self._NPDf_py_ is None:
            self._NPDf_py_ = NumericalPartialDerivative_txy_Functions(self.v)
        return self._NPDf_py_('y')(t, x, y)


    def f(self, t, x, y):
        return - self.u_x(t, x, y) - self.v_y(t, x, y)


    @property
    def potential(self):
        if self._potential_ is None:
            self._potential_ = mpRfT2_Scalar(
                self.mesh,
                self.phi,
                valid_time=self.valid_time,
                name='potential')
        return self._potential_


    @property
    def velocity(self):
        if self._velocity_ is None:
            self._velocity_ = mpRfT2_Vector(
                self.mesh,
                (self.u, self.v),
                valid_time=self.valid_time,
                name='velocity')
        return self._velocity_

    @property
    def source_term(self):
        if self._source_term_ is None:
            self._source_term_ = mpRfT2_Scalar(
                 self.mesh,
                 self.f,
                 valid_time=self.valid_time,
                 name='source_term')
        return self._source_term_



    @property
    def kinetic_energy_distribution(self):
        """A scalar field of the kinetic energy distribution."""
        if self._kineticEnergyDistribution_ is None:
            self._kineticEnergyDistribution_ =mpRfT2_Scalar(
                self.mesh,
                self.___kinetic_energy_distribution___,
                valid_time=self.valid_time,
                name='kinetic_energy_distribution')
        return self._kineticEnergyDistribution_
    def ___kinetic_energy_distribution___(self, t, x, y):
        return 0.5 * (self.u(t, x, y)**2 + self.v(t, x, y)**2)



    def kinetic_energy(self, t):
        """Kinetic energy at time `t`."""
        return self._es_.do.compute_Ln_norm_of('kinetic_energy_distribution', time=t, n=1)



    def ___PreFrozenChecker___(self):
        """
        We use this general method to do the check, in particular exact solution, we can define particular
        check method by override this method.
        """
        TS = self.___Pr_generate_random_valid_time_instances___()
        x, y = self._mesh_.do.generate_random_coordinates()

        if len(x) == 0: return

        for t in TS:

            t = float(t)

            try:
                Pu = NumericalPartialDerivative_txy(self.phi, t, x, y)
                assert Pu.check_partial_x(self.u)
                assert Pu.check_partial_y(self.v)
            except NotImplementedError:
                pass

            try:
                Pu = NumericalPartialDerivative_txy(self.u, t, x, y)
                assert Pu.check_partial_x(self.u_x)
            except NotImplementedError:
                pass

            try:
                Pu = NumericalPartialDerivative_txy(self.v, t, x, y)
                assert Pu.check_partial_y(self.v_y)
            except NotImplementedError:
                pass


            try:
                f = - self.u_x(t, x, y) - self.v_y(t, x, y)
                F = self.f(t, x, y)
                np.testing.assert_array_almost_equal(F - f, 0, decimal=5)
            except NotImplementedError:
                pass



if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
