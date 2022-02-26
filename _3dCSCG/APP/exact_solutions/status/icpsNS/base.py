# -*- coding: utf-8 -*-
"""
The parent of all exact solution for incompressible NS.

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import random
import numpy as np
from functools import partial, lru_cache
from _3dCSCG.APP.exact_solutions.status.base import Base
from screws.numerical._3d import NumericalPartialDerivative_xyz
from _3dCSCG.fields.vector import _3dCSCG_VectorField
from _3dCSCG.fields.scalar import _3dCSCG_ScalarField



class icpsNS_Base(Base):
    def __init__(self, es, nu, rho):
        super(icpsNS_Base, self).__init__(es)
        self._nu_ = nu
        self._rho_ = rho
        self._velocity_ = None
        self._vorticity_ = None
        self._bodyForce_ = None
        self._pressure_ = None
        self._gradientOfPressure_ = None
        self._totalPressure_ = None
        self._gradientOfTotalPressure_ = None
        self._kineticEnergyDistribution_ = None
        self._helicityDistribution_ = None
        self._enstrophyDistribution_ = None
        self._divergence_of_velocity_ = None
        self._divergence_of_vorticity_ = None
        self._curl_of_vorticity_ = None
        self.___check_self___()
        self._freeze_self_()

    @property
    def nu(self):
        return self._nu_

    @property
    def rho(self):
        return self._rho_

    # to be overridden (must) ...

    def u(self, t, x, y, z):
        raise NotImplementedError()
    def u_x(self, t, x, y, z):
        raise NotImplementedError()
    def u_y(self, t, x, y, z):
        raise NotImplementedError()
    def u_z(self, t, x, y, z):
        raise NotImplementedError()

    def v(self, t, x, y, z):
        raise NotImplementedError()
    def v_x(self, t, x, y, z):
        raise NotImplementedError()
    def v_y(self, t, x, y, z):
        raise NotImplementedError()
    def v_z(self, t, x, y, z):
        raise NotImplementedError()

    def w(self, t, x, y, z):
        raise NotImplementedError()
    def w_x(self, t, x, y, z):
        raise NotImplementedError()
    def w_y(self, t, x, y, z):
        raise NotImplementedError()
    def w_z(self, t, x, y, z):
        raise NotImplementedError()

    # optional ...

    def p(self, t, x, y, z): raise NotImplementedError()
    def p_x(self, t, x, y, z): raise NotImplementedError()
    def p_y(self, t, x, y, z): raise NotImplementedError()
    def p_z(self, t, x, y, z): raise NotImplementedError()

    # if p, p_x, p_y, p_z are provided, we need also provided following, such that we can compute the exact body force.

    def u_t(self, t, x, y, z):
        raise NotImplementedError()
    def u_xx(self, t, x, y, z):
        raise NotImplementedError()
    def u_yy(self, t, x, y, z):
        raise NotImplementedError()
    def u_zz(self, t, x, y, z):
        raise NotImplementedError()

    def v_t(self, t, x, y, z):
        raise NotImplementedError()
    def v_xx(self, t, x, y, z):
        raise NotImplementedError()
    def v_yy(self, t, x, y, z):
        raise NotImplementedError()
    def v_zz(self, t, x, y, z):
        raise NotImplementedError()

    def w_t(self, t, x, y, z):
        raise NotImplementedError()
    def w_xx(self, t, x, y, z):
        raise NotImplementedError()
    def w_yy(self, t, x, y, z):
        raise NotImplementedError()
    def w_zz(self, t, x, y, z):
        raise NotImplementedError()

    # .............................................................................

    @property
    def velocity(self):
        if self._velocity_ is None:
            self._velocity_ =  _3dCSCG_VectorField(self.mesh, (self.u, self.v, self.w), valid_time=self.valid_time)
        return self._velocity_

    @property
    def curl_of_velocity(self):
        return self.vorticity

    @property
    def divergence_of_velocity(self):
        if self._divergence_of_velocity_ is None:
            # this condition must be valid for all time.
            self._divergence_of_velocity_ = _3dCSCG_ScalarField(self.mesh, self.___div_of_velocity___,
                                                                valid_time=None)
        return self._divergence_of_velocity_

    def ___div_of_velocity___(self, t, x, y, z):
        # must be zero
        return self.u_x(t, x, y, z) + self.v_y(t, x, y, z) + self.w_z(t, x, y, z)

    @staticmethod
    def source_term(t, x, y, z):
        """By default, we have divergence free condition; the source term is zero."""
        return 0 * (x + y + z + t)


    @property
    def vorticity(self):
        if self._vorticity_ is None:
            self._vorticity_ = _3dCSCG_VectorField(self.mesh, (self.wx, self.wy, self.wz), valid_time=self.valid_time)
        return self._vorticity_
    def wx(self, t, x, y, z):
        return self.w_y(t, x, y, z) - self.v_z(t, x, y, z)
    def wy(self, t, x, y, z):
        return self.u_z(t, x, y, z) - self.w_x(t, x, y, z)
    def wz(self, t, x, y, z):
        return self.v_x(t, x, y, z) - self.u_y(t, x, y, z)



    @property
    def curl_of_vorticity(self):
        """ curl_of_vorticity = curl of curl u = grad(div u) - Vector_Laplace u.

        Since div u must be zero (conservation of mass), we have

        curl_of_vorticity = - Vector_Laplace u = - ( Laplace u_x,  Laplace u_y,  Laplace u_z) ^ T.

        """
        if self._curl_of_vorticity_ is None:
            self._curl_of_vorticity_ = _3dCSCG_VectorField(
                self.mesh, (self.___m_laplace_u___, self.___m_laplace_v___, self.___m_laplace_w___),
                valid_time=self.valid_time)
        return self._curl_of_vorticity_

    def ___m_laplace_u___(self, t, x, y, z):
        return - ( self.u_xx(t, x, y, z) + self.u_yy(t, x, y, z) + self.u_zz(t, x, y, z) )
    def ___m_laplace_v___(self, t, x, y, z):
        return - ( self.v_xx(t, x, y, z) + self.v_yy(t, x, y, z) + self.v_zz(t, x, y, z) )
    def ___m_laplace_w___(self, t, x, y, z):
        return - ( self.w_xx(t, x, y, z) + self.w_yy(t, x, y, z) + self.w_zz(t, x, y, z) )


    @property
    def divergence_of_vorticity(self):
        if self._divergence_of_vorticity_ is None:
            # this condition must be valid for all time.
            self._divergence_of_vorticity_ = _3dCSCG_ScalarField(self.mesh, self.___div_of_vorticity___,
                                                                valid_time=None)
        return self._divergence_of_vorticity_

    @staticmethod
    def ___div_of_vorticity___(t, x, y, z):
        # must be zero
        return 0 * t * x * y * z

    # Want to use specific body_force? then override fx, fy, fz in particular exact
    # solution class. Notice that, to have a exact solution for body force, we need to
    # know the exact solution for velocity and pressure.
    def _fx_(self, t, x, y, z):
        return self.u_t(t,x,y,z) + \
               self.wy(t,x,y,z)*self.w(t,x,y,z) - self.wz(t,x,y,z)*self.v(t,x,y,z) + \
               self._tp_x_(t, x, y, z) - \
               self.nu * (self.u_xx(t,x,y,z) + self.u_yy(t,x,y,z) + self.u_zz(t,x,y,z))
    def _fy_(self, t, x, y, z):
        return self.v_t(t,x,y,z) + \
               self.wz(t,x,y,z)*self.u(t,x,y,z) - self.wx(t,x,y,z)*self.w(t,x,y,z) + \
               self._tp_y_(t, x, y, z) - \
               self.nu * (self.v_xx(t,x,y,z) + self.v_yy(t,x,y,z) + self.v_zz(t,x,y,z))
    def _fz_(self, t, x, y, z):
        return self.w_t(t,x,y,z) + \
               self.wx(t,x,y,z)*self.v(t,x,y,z) - self.wy(t,x,y,z)*self.u(t,x,y,z) + \
               self._tp_z_(t, x, y, z) - \
               self.nu * (self.w_xx(t,x,y,z) + self.w_yy(t,x,y,z) + self.w_zz(t,x,y,z))

    @property
    def body_force(self):
        if self._bodyForce_ is None:
            self._bodyForce_ = _3dCSCG_VectorField(self.mesh, (self.fx, self.fy, self.fz),
                                                   valid_time=self.valid_time)
        return self._bodyForce_
    # @property
    # def fx(self):
    #     return self._fx_
    # @property
    # def fy(self):
    #     return self._fy_
    # @property
    # def fz(self):
    #     return self._fz_

    def fx(self, t, x, y, z):
        return self._fx_(t, x, y, z)

    def fy(self, t, x, y, z):
        return self._fy_(t, x, y, z)

    def fz(self, t, x, y, z):
        return self._fz_(t, x, y, z)

    @property
    def pressure(self):
        if self._pressure_ is None:
            self._pressure_ = _3dCSCG_ScalarField(self.mesh, self.p, valid_time=self.valid_time)
        return self._pressure_

    @property
    def gradient_of_pressure(self):
        if self._gradientOfPressure_ is None:
            self._gradientOfPressure_ =  _3dCSCG_VectorField(self.mesh, (self.p_x, self.p_y, self.p_z),
                                                             valid_time=self.valid_time)
        return self._gradientOfPressure_

    @property
    def total_pressure(self):
        if self._totalPressure_ is None:
            self._totalPressure_ = _3dCSCG_ScalarField(self.mesh, self.___total_pressure___,
                                                       valid_time=self.valid_time)
        return self._totalPressure_
    def ___total_pressure___(self, t, x, y, z):
        return self.p(t, x, y, z) / self.rho + self.___kinetic_energy_distribution___(t, x, y, z)
    def _tp_x_(self, t, x, y, z):
        """ """
        return self.p_x(t, x, y, z)/self.rho + self._ke_x_(t, x, y, z)
    def _tp_y_(self, t, x, y, z):
        """ """
        return self.p_y(t, x, y, z)/self.rho + self._ke_y_(t, x, y, z)
    def _tp_z_(self, t, x, y, z):
        """ """
        return self.p_z(t, x, y, z)/self.rho + self._ke_z_(t, x, y, z)

    @property
    def gradient_of_total_pressure(self):
        """A vector field of the total pressure distribution."""
        if self._gradientOfTotalPressure_ is None:
            self._gradientOfTotalPressure_ = _3dCSCG_VectorField(
                self.mesh, (self._tp_x_, self._tp_y_, self._tp_z_), valid_time=self.valid_time
            )
        return self._gradientOfTotalPressure_

    @property
    def kinetic_energy_distribution(self):
        """A scalar field of the kinetic energy distribution."""
        if self._kineticEnergyDistribution_ is None:
            self._kineticEnergyDistribution_ =_3dCSCG_ScalarField(
                self.mesh, self.___kinetic_energy_distribution___, valid_time=self.valid_time)
        return self._kineticEnergyDistribution_
    def ___kinetic_energy_distribution___(self, t, x, y, z):
        return 0.5 * (self.u(t, x, y, z)**2 + self.v(t, x, y, z)**2 + self.w(t, x, y, z)**2)
    def _ke_x_(self, t, x, y, z):
        """ d(kinetic_energy)/dx """
        return self.u(t, x, y, z)*self.u_x(t, x, y, z) + \
               self.v(t, x, y, z)*self.v_x(t, x, y, z) + \
               self.w(t, x, y, z)*self.w_x(t, x, y, z)
    def _ke_y_(self, t, x, y, z):
        """ d(kinetic_energy)/dy """
        return self.u(t, x, y, z)*self.u_y(t, x, y, z) + \
               self.v(t, x, y, z)*self.v_y(t, x, y, z) + \
               self.w(t, x, y, z)*self.w_y(t, x, y, z)
    def _ke_z_(self, t, x, y, z):
        """ d(kinetic_energy)/dz """
        return self.u(t, x, y, z)*self.u_z(t, x, y, z) + \
               self.v(t, x, y, z)*self.v_z(t, x, y, z) + \
               self.w(t, x, y, z)*self.w_z(t, x, y, z)

    @property
    def helicity_distribution(self):
        """A scalar field of the helicity distribution."""
        if self._helicityDistribution_ is None:
            self._helicityDistribution_ =_3dCSCG_ScalarField(
                self.mesh, self.___helicity_distribution___, valid_time=self.valid_time)
        return self._helicityDistribution_
    def ___helicity_distribution___(self, t, x, y, z):
        return self.u(t,x,y,z)*self.wx(t,x,y,z) + \
               self.v(t,x,y,z)*self.wy(t,x,y,z) + \
               self.w(t,x,y,z)*self.wz(t,x,y,z)

    @property
    def enstrophy_distribution(self):
        """A scalar field of the enstrophy distribution."""
        if self._enstrophyDistribution_ is None:
            self._enstrophyDistribution_ =_3dCSCG_ScalarField(
                self.mesh, self.___enstrophy_distribution___, valid_time=self.valid_time)
        return self._enstrophyDistribution_
    def ___enstrophy_distribution___(self, t, x, y, z):
        return 0.5 * (self.wx(t, x, y, z)**2 + self.wy(t, x, y, z)**2 + self.wz(t, x, y, z)**2)



    @lru_cache(maxsize=8)
    def helicity(self, t):
        """Helicity at time `t`."""
        return self.___PRIVATE_compute_L_norm_of___('helicity_distribution', time=t, n=1)

    @lru_cache(maxsize=8)
    def kinetic_energy(self, t):
        """Kinetic energy at time `t`."""
        return self.___PRIVATE_compute_L_norm_of___('kinetic_energy_distribution', time=t, n=1)

    @lru_cache(maxsize=8)
    def enstrophy(self, t):
        """Enstrophy at time `t`."""
        return self.___PRIVATE_compute_L_norm_of___('enstrophy_distribution', time=t, n=1)



    @property
    def ___check_domain___(self):
        """
        We use this general domain to do the check, in particular exact solution, we can define new domain
        by override this method.

        """
        time0, time1, time2 = random.uniform(-1, 1), random.uniform(-2, 2), random.uniform(-3, 3)
        r = np.linspace(random.uniform(-1, -0.9), random.uniform(0.82, 1), random.randint(1, 3))
        s = np.linspace(random.uniform(-0.99, -0.95), random.uniform(0.88, 0.95), random.randint(2, 3))
        t = np.linspace(-1, 1, random.randint(1, 2))
        r, s, t = np.meshgrid(r, s, t, indexing='ij')
        dt = 0.000001
        return (r, s, t), (time0, time1, time2), dt

    def ___check_self___(self):
        """
        We use this general method to do the check, in particular exact solution, we can define particular
        check method by override this method.
        """
        rst, times, dt = self. ___check_domain___
        for time in times:
            u = partial(self.u, time)
            u_x = partial(self.u_x, time)
            u_y = partial(self.u_y, time)
            u_z = partial(self.u_z, time)
            Pu = NumericalPartialDerivative_xyz(u, *rst)
            assert all(Pu.check_total(u_x, u_y, u_z))
            v = partial(self.v, time)
            v_x = partial(self.v_x, time)
            v_y = partial(self.v_y, time)
            v_z = partial(self.v_z, time)
            Pv = NumericalPartialDerivative_xyz(v, *rst)
            assert all(Pv.check_total(v_x, v_y, v_z))
            w = partial(self.w, time)
            w_x = partial(self.w_x, time)
            w_y = partial(self.w_y, time)
            w_z = partial(self.w_z, time)
            Pw = NumericalPartialDerivative_xyz(w, *rst)
            assert all(Pw.check_total(w_x, w_y, w_z))
            assert np.all(np.isclose(u_x(*rst) + v_y(*rst) + w_z(*rst), self.source_term(time, *rst))), \
                " <exact solution> <icpNS> : velocity divergence condition."
            assert np.all(np.isclose(self.___div_of_velocity___(time, *rst), self.source_term(time, *rst))), \
                " <exact solution> <icpNS> : velocity divergence free condition."

            wx = partial(self.wx, time)
            wy = partial(self.wy, time)
            wz = partial(self.wz, time)
            Pwx = NumericalPartialDerivative_xyz(wx, *rst)
            Pwy = NumericalPartialDerivative_xyz(wy, *rst)
            Pwz = NumericalPartialDerivative_xyz(wz, *rst)
            Pwx_x = Pwx.scipy_partial('x')
            Pwy_y = Pwy.scipy_partial('y')
            Pwz_z = Pwz.scipy_partial('z')
            A = np.sum(np.abs(Pwx_x + Pwy_y + Pwz_z))
            np.testing.assert_almost_equal(A, 0)

        try:
            _ = self.p(0,0,0,0)
        except NotImplementedError:
            # then the exact solution is not FULL!
            # The idea is: when we do not give pressure, then we have no exact solutions
            # for all variables, so we will not do future check.
            pass
        else:
            # else, we have pressure. Then we can derive the body force. We need to firstly
            # check the second-order partial differential of the velocity
            for time in times:
                p = partial(self.p, time)
                p_x = partial(self.p_x, time)
                p_y = partial(self.p_y, time)
                p_z = partial(self.p_z, time)
                Pp = NumericalPartialDerivative_xyz(p, *rst)
                assert all(Pp.check_total(p_x, p_y, p_z))
                u_x = partial(self.u_x, time)
                u_xx = partial(self.u_xx, time)
                u_y = partial(self.u_y, time)
                u_yy = partial(self.u_yy, time)
                u_z = partial(self.u_z, time)
                u_zz = partial(self.u_zz, time)
                v_x = partial(self.v_x, time)
                v_xx = partial(self.v_xx, time)
                v_y = partial(self.v_y, time)
                v_yy = partial(self.v_yy, time)
                v_z = partial(self.v_z, time)
                v_zz = partial(self.v_zz, time)
                w_x = partial(self.w_x, time)
                w_xx = partial(self.w_xx, time)
                w_y = partial(self.w_y, time)
                w_yy = partial(self.w_yy, time)
                w_z = partial(self.w_z, time)
                w_zz = partial(self.w_zz, time)
                Pu_x = NumericalPartialDerivative_xyz(u_x, *rst)
                assert Pu_x.check_partial_x(u_xx)
                Pu_y = NumericalPartialDerivative_xyz(u_y, *rst)
                assert Pu_y.check_partial_y(u_yy)
                Pu_z = NumericalPartialDerivative_xyz(u_z, *rst)
                assert Pu_z.check_partial_z(u_zz)
                Pv_x = NumericalPartialDerivative_xyz(v_x, *rst)
                assert Pv_x.check_partial_x(v_xx)
                Pv_y = NumericalPartialDerivative_xyz(v_y, *rst)
                assert Pv_y.check_partial_y(v_yy)
                Pv_z = NumericalPartialDerivative_xyz(v_z, *rst)
                assert Pv_z.check_partial_z(v_zz)
                Pw_x = NumericalPartialDerivative_xyz(w_x, *rst)
                assert Pw_x.check_partial_x(w_xx)
                Pw_y = NumericalPartialDerivative_xyz(w_y, *rst)
                assert Pw_y.check_partial_y(w_yy)
                Pw_z = NumericalPartialDerivative_xyz(w_z, *rst)
                assert Pw_z.check_partial_z(w_zz)
                R, S, T = rst
                R = R.ravel('F')
                S = S.ravel('F')
                T = T.ravel('F')
                for i, r in enumerate(R):
                    s = S[i]
                    t = T[i]
                    du_dt = (self.u(time+dt, r, s, t) - self.u(time-dt, r, s, t))/(2*dt)
                    assert np.isclose(du_dt, self.u_t(time,r,s,t)), " <> : u_t wrong."
                    dv_dt = (self.v(time+dt, r, s, t) - self.v(time-dt, r, s, t))/(2*dt)
                    assert np.isclose(dv_dt, self.v_t(time,r,s,t)), " <> : v_t wrong."
                    dw_dt = (self.w(time+dt, r, s, t) - self.w(time-dt, r, s, t))/(2*dt)
                    assert np.isclose(dw_dt, self.w_t(time,r,s,t)), " <> : w_t wrong."
                _fx_, fx = self._fx_(time, *rst), self.fx(time, *rst)
                _fy_, fy = self._fy_(time, *rst), self.fy(time, *rst)
                _fz_, fz = self._fz_(time, *rst), self.fz(time, *rst)
                np.testing.assert_almost_equal(fx, _fx_)
                np.testing.assert_almost_equal(fy, _fy_)
                np.testing.assert_almost_equal(fz, _fz_)