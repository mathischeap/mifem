"""The 3d Poisson problem, k = 1.

The Poisson equation:
- div (k grad phi) = f

And in the mixed formulation:

u = k grad phi
f = - div u

And we set k equal to 1; k = 1.

We will interpolate the Poisson problem as a flow problem (potential flow). So phi is the potential,
(u,v,w) is the velocity, and f is the source term. Note that they are just names, they can also be
for example (temperature, heat flux and heat source) and so on.

"""



import random
import numpy as np
from functools import partial, lru_cache
from _3dCSCG.APP.exact_solution.status.base import Base
from screws.numerical._3d_space.Jacobian_33 import NumericalPartialDerivative_xyz
from _3dCSCG.fields.vector.main import _3dCSCG_VectorField
from _3dCSCG.fields.scalar.main import _3dCSCG_ScalarField


class Poisson_Base(Base):
    def __init__(self, es):
        super(Poisson_Base, self).__init__(es)
        self._potential_ = None
        self._velocity_ = None
        self._source_term_ = None
        self._kineticEnergyDistribution_ = None
        self.___check_self___()
        self._freeze_self_()


    @property
    def ___check_domain___(self):
        """
        We use this general domain to do the check, in particular exact solution, we can define new domain
        by override this method.

        """
        times = [random.uniform(-5, 5) for _ in range(3)]
        r = np.linspace(random.uniform(-1, -0.1), random.uniform(0.1, 1), random.randint(3, 5))
        s = np.linspace(random.uniform(-1, -0.2), random.uniform(0.2, 1), random.randint(3, 5))
        t = np.linspace(random.uniform(-1, -0.3), random.uniform(0.3, 1), random.randint(3, 5))
        r, s, t = np.meshgrid(r, s, t, indexing='ij')
        dt = 0.000001
        return (r, s, t), times, dt


    def ___check_self___(self):
        """
        A default checker. If you do not want to run this check, override ___PRIVATE_check_self___ method in the
        particular class.

        :return:
        """
        rst, times, dt = self.___check_domain___
        for time in times:
            phi = partial(self.phi, time)
            u = partial(self.u, time)
            v = partial(self.v, time)
            w = partial(self.w, time)
            P_phi = NumericalPartialDerivative_xyz(phi, *rst)
            assert all(P_phi.check_total(u, v, w))

            ux = partial(self.u_x, time)
            vy = partial(self.v_y, time)
            wz = partial(self.w_z, time)

            P_u = NumericalPartialDerivative_xyz(u, *rst)
            P_v = NumericalPartialDerivative_xyz(v, *rst)
            P_w = NumericalPartialDerivative_xyz(w, *rst)

            assert P_u.check_partial_x(ux)
            assert P_v.check_partial_y(vy)
            assert P_w.check_partial_z(wz)

            f = partial(self.f, time)(*rst)
            np.testing.assert_array_almost_equal(f, -(ux(*rst) + vy(*rst) + wz(*rst)))


    # to be overridden (must) ...

    def phi(self, t, x, y, z):
        raise NotImplementedError()

    def u(self, t, x, y, z):
        """phi_x"""
        raise NotImplementedError()
    def u_x(self, t, x, y, z):
        raise NotImplementedError()

    def v(self, t, x, y, z):
        """phi_y"""
        raise NotImplementedError()
    def v_y(self, t, x, y, z):
        raise NotImplementedError()

    def w(self, t, x, y, z):
        """phi_z"""
        raise NotImplementedError()
    def w_z(self, t, x, y, z):
        raise NotImplementedError()


    def _f_(self, t, x, y, z):
        """To define a specific f, we override self.f, please do not self._f_."""
        return - self.u_x(t, x, y, z) - self.v_y(t, x, y, z) - self.w_z(t, x, y, z)

    def f(self, t, x, y, z):
        return self._f_(t, x, y, z)


    # BELOW: properties ....................................................................



    @property
    def pressure(self):
        if self._potential_ is None:
            self._potential_ = _3dCSCG_ScalarField(self.mesh, self.phi, valid_time=self.valid_time)
        return self._potential_

    @property
    def velocity(self):
        if self._velocity_ is None:
            self._velocity_ = _3dCSCG_VectorField(self.mesh, (self.u, self.v, self.w), valid_time=self.valid_time)
        return self._velocity_

    @property
    def source_term(self):
        if self._source_term_ is None:
            self._source_term_ = _3dCSCG_ScalarField(self.mesh, self.phi, valid_time=self.valid_time)
        return self._source_term_



    @property
    def kinetic_energy_distribution(self):
        """A scalar field of the kinetic energy distribution."""
        if self._kineticEnergyDistribution_ is None:
            self._kineticEnergyDistribution_ =_3dCSCG_ScalarField(
                self.mesh, self.___kinetic_energy_distribution___, valid_time=self.valid_time)
        return self._kineticEnergyDistribution_
    def ___kinetic_energy_distribution___(self, t, x, y, z):
        return 0.5 * (self.u(t, x, y, z)**2 + self.v(t, x, y, z)**2 + self.w(t, x, y, z)**2)
    @lru_cache(maxsize=8)
    def kinetic_energy(self, t):
        """Kinetic energy at time `t`."""
        return self._es_.do.compute_Ln_norm_of('kinetic_energy_distribution', time=t, n=1)

