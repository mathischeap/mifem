# -*- coding: utf-8 -*-
"""
The 3d stokes flow.

w = curl u

curl w + grad p = f

div u = 0


"""



from functools import lru_cache, partial
import numpy as np

from objects.CSCG._3d.exactSolutions.base import Base
from objects.CSCG._3d.fields.vector.main import _3dCSCG_VectorField
from objects.CSCG._3d.fields.scalar.main import _3dCSCG_ScalarField
from screws.numerical._3dSpace.partial_derivative import NumericalPartialDerivative_xyz
from screws.numerical.timePlus3dSpace.partial_derivative_as_functions import \
    NumericalPartialDerivative_txyz_Functions




class Stokes_Base(Base):
    """"""
    def __init__(self, mesh):
        """"""
        super(Stokes_Base, self).__init__(mesh)
        self._pressure_ = None
        self._velocity_ = None
        self._vorticity_ = None
        self._divergence_of_velocity_ = None
        self._curl_of_vorticity_ = None
        self._divergence_of_vorticity_ = None
        self._kineticEnergyDistribution_ = None
        self._bodyForce_ = None
        self._gradientOfPressure_ = None

        self._NPDf_u_ = None
        self._NPDf_v_ = None
        self._NPDf_w_ = None

        self._NPDf_p_ = None

        self._NPDf_ux_ = None
        self._NPDf_uy_ = None
        self._NPDf_uz_ = None

        self._NPDf_vx_ = None
        self._NPDf_vy_ = None
        self._NPDf_vz_ = None

        self._NPDf_wx_ = None
        self._NPDf_wy_ = None
        self._NPDf_wz_ = None

        self._freeze_self_()



    def u(self, t, x, y, z):
        raise NotImplementedError()
    def u_x(self, t, x, y, z):
        if self._NPDf_u_ is None:
            self._NPDf_u_ = NumericalPartialDerivative_txyz_Functions(self.u)
        return self._NPDf_u_('x')(t, x, y, z)
    def u_y(self, t, x, y, z):
        if self._NPDf_u_ is None:
            self._NPDf_u_ = NumericalPartialDerivative_txyz_Functions(self.u)
        return self._NPDf_u_('y')(t, x, y, z)
    def u_z(self, t, x, y, z):
        if self._NPDf_u_ is None:
            self._NPDf_u_ = NumericalPartialDerivative_txyz_Functions(self.u)
        return self._NPDf_u_('z')(t, x, y, z)


    def v(self, t, x, y, z):
        raise NotImplementedError()
    def v_x(self, t, x, y, z):
        if self._NPDf_v_ is None:
            self._NPDf_v_ = NumericalPartialDerivative_txyz_Functions(self.v)
        return self._NPDf_v_('x')(t, x, y, z)
    def v_y(self, t, x, y, z):
        if self._NPDf_v_ is None:
            self._NPDf_v_ = NumericalPartialDerivative_txyz_Functions(self.v)
        return self._NPDf_v_('y')(t, x, y, z)
    def v_z(self, t, x, y, z):
        if self._NPDf_v_ is None:
            self._NPDf_v_ = NumericalPartialDerivative_txyz_Functions(self.v)
        return self._NPDf_v_('z')(t, x, y, z)


    def w(self, t, x, y, z):
        raise NotImplementedError()
    def w_x(self, t, x, y, z):
        if self._NPDf_w_ is None:
            self._NPDf_w_ = NumericalPartialDerivative_txyz_Functions(self.w)
        return self._NPDf_w_('x')(t, x, y, z)
    def w_y(self, t, x, y, z):
        if self._NPDf_w_ is None:
            self._NPDf_w_ = NumericalPartialDerivative_txyz_Functions(self.w)
        return self._NPDf_w_('y')(t, x, y, z)
    def w_z(self, t, x, y, z):
        if self._NPDf_w_ is None:
            self._NPDf_w_ = NumericalPartialDerivative_txyz_Functions(self.w)
        return self._NPDf_w_('z')(t, x, y, z)


    def p(self, t, x, y, z): raise NotImplementedError()
    def p_x(self, t, x, y, z):
        if self._NPDf_p_ is None:
            self._NPDf_p_ = NumericalPartialDerivative_txyz_Functions(self.p)
        return self._NPDf_p_('x')(t, x, y, z)
    def p_y(self, t, x, y, z):
        if self._NPDf_p_ is None:
            self._NPDf_p_ = NumericalPartialDerivative_txyz_Functions(self.p)
        return self._NPDf_p_('y')(t, x, y, z)
    def p_z(self, t, x, y, z):
        if self._NPDf_p_ is None:
            self._NPDf_p_ = NumericalPartialDerivative_txyz_Functions(self.p)
        return self._NPDf_p_('z')(t, x, y, z)


    def u_xx(self, t, x, y, z):
        if self._NPDf_ux_ is None:
            self._NPDf_ux_ = NumericalPartialDerivative_txyz_Functions(self.u_x)
        return self._NPDf_ux_('x')(t, x, y, z)
    def u_yy(self, t, x, y, z):
        if self._NPDf_uy_ is None:
            self._NPDf_uy_ = NumericalPartialDerivative_txyz_Functions(self.u_y)
        return self._NPDf_uy_('y')(t, x, y, z)
    def u_zz(self, t, x, y, z):
        if self._NPDf_uz_ is None:
            self._NPDf_uz_ = NumericalPartialDerivative_txyz_Functions(self.u_z)
        return self._NPDf_uz_('z')(t, x, y, z)

    def v_xx(self, t, x, y, z):
        if self._NPDf_vx_ is None:
            self._NPDf_vx_ = NumericalPartialDerivative_txyz_Functions(self.v_x)
        return self._NPDf_vx_('x')(t, x, y, z)
    def v_yy(self, t, x, y, z):
        if self._NPDf_vy_ is None:
            self._NPDf_vy_ = NumericalPartialDerivative_txyz_Functions(self.v_y)
        return self._NPDf_vy_('y')(t, x, y, z)
    def v_zz(self, t, x, y, z):
        if self._NPDf_vz_ is None:
            self._NPDf_vz_ = NumericalPartialDerivative_txyz_Functions(self.v_z)
        return self._NPDf_vz_('z')(t, x, y, z)

    def w_xx(self, t, x, y, z):
        if self._NPDf_wx_ is None:
            self._NPDf_wx_ = NumericalPartialDerivative_txyz_Functions(self.w_x)
        return self._NPDf_wx_('x')(t, x, y, z)
    def w_yy(self, t, x, y, z):
        if self._NPDf_wy_ is None:
            self._NPDf_wy_ = NumericalPartialDerivative_txyz_Functions(self.w_y)
        return self._NPDf_wy_('y')(t, x, y, z)
    def w_zz(self, t, x, y, z):
        if self._NPDf_wz_ is None:
            self._NPDf_wz_ = NumericalPartialDerivative_txyz_Functions(self.w_z)
        return self._NPDf_wz_('z')(t, x, y, z)

    # ---------------------------------------------------------------

    @property
    def pressure(self):
        if self._pressure_ is None:
            self._pressure_ = _3dCSCG_ScalarField(
                self.mesh, self.p, valid_time=self.valid_time)
        return self._pressure_

    @property
    def gradient_of_pressure(self):
        if self._gradientOfPressure_ is None:
            self._gradientOfPressure_ =  _3dCSCG_VectorField(
                self.mesh, (self.p_x, self.p_y, self.p_z), valid_time=self.valid_time)
        return self._gradientOfPressure_

    @property
    def velocity(self):
        if self._velocity_ is None:
            self._velocity_ =  _3dCSCG_VectorField(
                self.mesh, (self.u, self.v, self.w), valid_time=self.valid_time)
        return self._velocity_

    def omega_x(self, t, x, y, z):
        return self.w_y(t, x, y, z) - self.v_z(t, x, y, z)
    def omega_y(self, t, x, y, z):
        return self.u_z(t, x, y, z) - self.w_x(t, x, y, z)
    def omega_z(self, t, x, y, z):
        return self.v_x(t, x, y, z) - self.u_y(t, x, y, z)
    @property
    def vorticity(self):
        if self._vorticity_ is None:
            self._vorticity_ = _3dCSCG_VectorField(
                self.mesh, (self.omega_x, self.omega_y, self.omega_z), valid_time=self.valid_time)
        return self._vorticity_

    @property
    def curl_of_velocity(self):
        return self.vorticity

    @property
    def divergence_of_velocity(self):
        if self._divergence_of_velocity_ is None:
            # this condition must be valid for all time.
            self._divergence_of_velocity_ = _3dCSCG_ScalarField(
                self.mesh, self.___div_of_velocity___, valid_time=None)
        return self._divergence_of_velocity_

    def ___div_of_velocity___(self, t, x, y, z):
        # must be zero
        return self.u_x(t, x, y, z) + self.v_y(t, x, y, z) + self.w_z(t, x, y, z)

    @staticmethod
    def source_term(t, x, y, z):
        """By default, we have divergence free condition; the source term is zero."""
        return 0 * (x + y + z + t)

    @property
    def curl_of_vorticity(self):
        """ curl_of_vorticity = curl of curl u = grad(div u) - Vector_Laplace u.

        Since div u must be zero (conservation of mass), we have

        curl_of_vorticity = - Vector_Laplace u = - ( Laplace u_x,  Laplace u_y,  Laplace u_z) ^ T.

        """
        if self._curl_of_vorticity_ is None:
            self._curl_of_vorticity_ = _3dCSCG_VectorField(
                self.mesh,
                (self.___m_laplace_u___, self.___m_laplace_v___, self.___m_laplace_w___),
                valid_time=self.valid_time
            )
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
            self._divergence_of_vorticity_ = _3dCSCG_ScalarField(
                self.mesh, self.___div_of_vorticity___, valid_time=None)
        return self._divergence_of_vorticity_

    @staticmethod
    def ___div_of_vorticity___(t, x, y, z):
        # must be zero
        return 0 * t * x * y * z



    def fx(self, t, x, y, z):
        return self.___m_laplace_u___(t,x,y,z) + self.p_x(t,x,y,z)
    def fy(self, t, x, y, z):
        return self.___m_laplace_v___(t,x,y,z) + self.p_y(t,x,y,z)
    def fz(self, t, x, y, z):
        return self.___m_laplace_w___(t,x,y,z) + self.p_z(t,x,y,z)

    @property
    def body_force(self):
        if self._bodyForce_ is None:
            self._bodyForce_ = _3dCSCG_VectorField(self.mesh,
                                                   (self.fx, self.fy, self.fz),
                                                   valid_time=self.valid_time)
        return self._bodyForce_

    @property
    def kinetic_energy_distribution(self):
        """A scalar field of the kinetic energy distribution."""
        if self._kineticEnergyDistribution_ is None:
            self._kineticEnergyDistribution_ =_3dCSCG_ScalarField(
                self.mesh,
                self.___kinetic_energy_distribution___,
                valid_time=self.valid_time)
        return self._kineticEnergyDistribution_
    def ___kinetic_energy_distribution___(self, t, x, y, z):
        return 0.5 * (self.u(t, x, y, z)**2 + self.v(t, x, y, z)**2 + self.w(t, x, y, z)**2)


    @lru_cache(maxsize=8)
    def kinetic_energy(self, t):
        """Kinetic energy at time `t`."""
        return self.___Pr_compute_Ln_norm_of___('kinetic_energy_distribution', time=t, n=1)


    def ___PreFrozenChecker___(self):
        """Will be called before freezing self."""
        TS = self.___PRIVATE_generate_random_valid_time_instances___()
        x, y, z = self._mesh_.do.generate_random_coordinates()

        rst = (x, y, z)

        if len(x) == 0: return

        for t in TS:
            fx = self.___m_laplace_u___(t,x,y,z) + self.p_x(t,x,y,z)

            fy = self.___m_laplace_v___(t,x,y,z) + self.p_y(t,x,y,z)

            fz = self.___m_laplace_w___(t,x,y,z) + self.p_z(t,x,y,z)

            FX = self.fx(t, x, y, z)
            FY = self.fy(t, x, y, z)
            FZ = self.fz(t, x, y, z)
            np.testing.assert_array_almost_equal(FX - fx, 0, decimal=5)
            np.testing.assert_array_almost_equal(FY - fy, 0, decimal=5)
            np.testing.assert_array_almost_equal(FZ - fz, 0, decimal=5)

            time = t

            assert np.all(np.isclose(
                self.___div_of_velocity___(time, *rst),
                self.source_term(time, *rst))), \
                " <exact solution> <icpNS> : velocity divergence free condition."

            wx = partial(self.omega_x, time)
            wy = partial(self.omega_y, time)
            wz = partial(self.omega_z, time)
            Pwx = NumericalPartialDerivative_xyz(wx, *rst)
            Pwy = NumericalPartialDerivative_xyz(wy, *rst)
            Pwz = NumericalPartialDerivative_xyz(wz, *rst)
            Pwx_x = Pwx.scipy_partial('x')
            Pwy_y = Pwy.scipy_partial('y')
            Pwz_z = Pwz.scipy_partial('z')
            A = np.sum(np.abs(Pwx_x + Pwy_y + Pwz_z))
            np.testing.assert_almost_equal(A, 0)

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

