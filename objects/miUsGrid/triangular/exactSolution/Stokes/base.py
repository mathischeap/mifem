# -*- coding: utf-8 -*-
"""

w = curl u

curl w + grad p = f

div u = 0


@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 10/10/2022 1:27 AM
"""
import sys

if './' not in sys.path: sys.path.append('./')

import numpy as np
from objects.miUsGrid.triangular.exactSolution.base import miUsTriangle_ExactSolutionBase

from functools import partial

from components.numerical.timePlus2dSpace.partial_derivative_as_functions import \
    NumericalPartialDerivative_txy_Functions

from objects.miUsGrid.triangular.fields.vector.main import miUsGrid_Triangular_Vector
from objects.miUsGrid.triangular.fields.scalar.main import miUsGrid_Triangular_Scalar

from components.numerical._2dSpace.partial_derivative import NumericalPartialDerivative_xy

class Stokes(miUsTriangle_ExactSolutionBase):
    """"""

    def __init__(self, mesh):
        """"""
        super(Stokes, self).__init__(mesh)

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
        self._NPDf_p_ = None

        self._NPDf_ux_ = None
        self._NPDf_uy_ = None

        self._NPDf_vx_ = None
        self._NPDf_vy_ = None

        self._freeze_self_()



    def u(self, t, x, y):
        raise NotImplementedError()
    def u_x(self, t, x, y):
        if self._NPDf_u_ is None:
            self._NPDf_u_ = NumericalPartialDerivative_txy_Functions(self.u)
        return self._NPDf_u_('x')(t, x, y)
    def u_y(self, t, x, y):
        if self._NPDf_u_ is None:
            self._NPDf_u_ = NumericalPartialDerivative_txy_Functions(self.u)
        return self._NPDf_u_('y')(t, x, y)



    def v(self, t, x, y):
        raise NotImplementedError()
    def v_x(self, t, x, y):
        if self._NPDf_v_ is None:
            self._NPDf_v_ = NumericalPartialDerivative_txy_Functions(self.v)
        return self._NPDf_v_('x')(t, x, y)
    def v_y(self, t, x, y):
        if self._NPDf_v_ is None:
            self._NPDf_v_ = NumericalPartialDerivative_txy_Functions(self.v)
        return self._NPDf_v_('y')(t, x, y)




    def p(self, t, x, y): raise NotImplementedError()
    def p_x(self, t, x, y):
        if self._NPDf_p_ is None:
            self._NPDf_p_ = NumericalPartialDerivative_txy_Functions(self.p)
        return self._NPDf_p_('x')(t, x, y)
    def p_y(self, t, x, y):
        if self._NPDf_p_ is None:
            self._NPDf_p_ = NumericalPartialDerivative_txy_Functions(self.p)
        return self._NPDf_p_('y')(t, x, y)


    def u_xx(self, t, x, y):
        if self._NPDf_ux_ is None:
            self._NPDf_ux_ = NumericalPartialDerivative_txy_Functions(self.u_x)
        return self._NPDf_ux_('x')(t, x, y)
    def u_yy(self, t, x, y):
        if self._NPDf_uy_ is None:
            self._NPDf_uy_ = NumericalPartialDerivative_txy_Functions(self.u_y)
        return self._NPDf_uy_('y')(t, x, y)

    def v_xx(self, t, x, y):
        if self._NPDf_vx_ is None:
            self._NPDf_vx_ = NumericalPartialDerivative_txy_Functions(self.v_x)
        return self._NPDf_vx_('x')(t, x, y)
    def v_yy(self, t, x, y):
        if self._NPDf_vy_ is None:
            self._NPDf_vy_ = NumericalPartialDerivative_txy_Functions(self.v_y)
        return self._NPDf_vy_('y')(t, x, y)

    # ---------------------------------------------------------------

    @property
    def pressure(self):
        if self._pressure_ is None:
            self._pressure_ = miUsGrid_Triangular_Scalar(
                self.mesh, self.p, valid_time=self.valid_time)
        return self._pressure_

    @property
    def gradient_of_pressure(self):
        if self._gradientOfPressure_ is None:
            self._gradientOfPressure_ =  miUsGrid_Triangular_Vector(
                self.mesh, (self.p_x, self.p_y), valid_time=self.valid_time)
        return self._gradientOfPressure_

    @property
    def velocity(self):
        if self._velocity_ is None:
            self._velocity_ =  miUsGrid_Triangular_Vector(
                self.mesh, (self.u, self.v), valid_time=self.valid_time)
        return self._velocity_

    def omega_z(self, t, x, y):
        return self.v_x(t, x, y) - self.u_y(t, x, y)
    @property
    def vorticity(self):
        if self._vorticity_ is None:
            self._vorticity_ = miUsGrid_Triangular_Scalar(
                self.mesh, self.omega_z, valid_time=self.valid_time)
        return self._vorticity_

    @property
    def curl_of_velocity(self):
        return self.vorticity

    @property
    def divergence_of_velocity(self):
        if self._divergence_of_velocity_ is None:
            # this condition must be valid for all time.
            self._divergence_of_velocity_ = miUsGrid_Triangular_Scalar(
                self.mesh, self.___div_of_velocity___, valid_time=None)
        return self._divergence_of_velocity_

    def ___div_of_velocity___(self, t, x, y):
        # must be zero
        return self.u_x(t, x, y) + self.v_y(t, x, y)

    @staticmethod
    def source_term(t, x, y):
        """By default, we have divergence free condition; the source term is zero."""
        return 0 * (x + y + t)

    @property
    def curl_of_vorticity(self):
        """ curl_of_vorticity = curl of curl u = grad(div u) - Vector_Laplace u.

        Since div u must be zero (conservation of mass), we have

        curl_of_vorticity = - Vector_Laplace u = - ( Laplace u_x,  Laplace u_y) ^ T.

        """
        if self._curl_of_vorticity_ is None:
            self._curl_of_vorticity_ = miUsGrid_Triangular_Vector(
                self.mesh,
                (self.___m_laplace_u___, self.___m_laplace_v___),
                valid_time=self.valid_time
            )
        return self._curl_of_vorticity_

    def ___m_laplace_u___(self, t, x, y):
        return - ( self.u_xx(t, x, y) + self.u_yy(t, x, y) )
    def ___m_laplace_v___(self, t, x, y):
        return - ( self.v_xx(t, x, y) + self.v_yy(t, x, y) )


    @property
    def divergence_of_vorticity(self):
        if self._divergence_of_vorticity_ is None:
            # this condition must be valid for all time.
            self._divergence_of_vorticity_ = miUsGrid_Triangular_Scalar(
                self.mesh, self.___div_of_vorticity___, valid_time=None)
        return self._divergence_of_vorticity_

    @staticmethod
    def ___div_of_vorticity___(t, x, y):
        # must be zero
        return 0 * t * x * y



    def fx(self, t, x, y):
        return self.___m_laplace_u___(t,x,y) + self.p_x(t,x,y)
    def fy(self, t, x, y):
        return self.___m_laplace_v___(t,x,y) + self.p_y(t,x,y)

    @property
    def body_force(self):
        if self._bodyForce_ is None:
            self._bodyForce_ = miUsGrid_Triangular_Vector(self.mesh,
                                                   (self.fx, self.fy),
                                                   valid_time=self.valid_time)
        return self._bodyForce_

    @property
    def kinetic_energy_distribution(self):
        """A scalar field of the kinetic energy distribution."""
        if self._kineticEnergyDistribution_ is None:
            self._kineticEnergyDistribution_ =miUsGrid_Triangular_Scalar(
                self.mesh,
                self.___kinetic_energy_distribution___,
                valid_time=self.valid_time)
        return self._kineticEnergyDistribution_
    def ___kinetic_energy_distribution___(self, t, x, y):
        return 0.5 * (self.u(t, x, y)**2 + self.v(t, x, y)**2 )




    def ___PreFrozenChecker___(self):
        """Will be called before freezing self."""
        TS = self.___PRIVATE_generate_random_valid_time_instances___()
        x, y = self._mesh_.do.generate_random_coordinates()

        rst = (x, y)

        if len(x) == 0: return

        for t in TS:
            fx = self.___m_laplace_u___(t,x,y) + self.p_x(t,x,y)

            fy = self.___m_laplace_v___(t,x,y) + self.p_y(t,x,y)

            FX = self.fx(t, x, y)
            FY = self.fy(t, x, y)
            np.testing.assert_array_almost_equal(FX - fx, 0, decimal=5)
            np.testing.assert_array_almost_equal(FY - fy, 0, decimal=5)

            time = t

            assert np.all(np.isclose(
                self.___div_of_velocity___(time, *rst),
                self.source_term(time, *rst))), \
                " <exact solution> <icpNS> : velocity divergence free condition."

            p = partial(self.p, time)
            p_x = partial(self.p_x, time)
            p_y = partial(self.p_y, time)
            Pp = NumericalPartialDerivative_xy(p, *rst)
            assert all(Pp.check_total(p_x, p_y))

            u_x = partial(self.u_x, time)
            u_xx = partial(self.u_xx, time)
            u_y = partial(self.u_y, time)
            u_yy = partial(self.u_yy, time)
            v_x = partial(self.v_x, time)
            v_xx = partial(self.v_xx, time)
            v_y = partial(self.v_y, time)
            v_yy = partial(self.v_yy, time)

            Pu_x = NumericalPartialDerivative_xy(u_x, *rst)
            assert Pu_x.check_partial_x(u_xx)
            Pu_y = NumericalPartialDerivative_xy(u_y, *rst)
            assert Pu_y.check_partial_y(u_yy)
            Pv_x = NumericalPartialDerivative_xy(v_x, *rst)
            assert Pv_x.check_partial_x(v_xx)
            Pv_y = NumericalPartialDerivative_xy(v_y, *rst)
            assert Pv_y.check_partial_y(v_yy)




if __name__ == '__main__':
    # mpiexec -n 4 python 
    pass
