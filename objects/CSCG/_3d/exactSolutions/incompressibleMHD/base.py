# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 8/20/2022 3:09 PM

incompressible MHD reads (dimensionless):

Pu/Pt + w X u + (1/Re)*curl(w) - c * (j X B) + grad(P) = f
w = curl(u)
div(u) = 0
PB/Pt + curl(E) = 0
(1/Rm) * j - E - u X B = 0
j = curl(B)

"""
from functools import lru_cache, partial
from objects.CSCG._3d.exactSolutions.base import Base
from objects.CSCG._3d.fields.vector.main import _3dCSCG_VectorField
from objects.CSCG._3d.fields.scalar.main import _3dCSCG_ScalarField
from components.numerical._3dSpace.partial_derivative import NumericalPartialDerivative_xyz
from components.numerical.timePlus3dSpace.partial_derivative import NumericalPartialDerivative_txyz
from components.numerical.timePlus3dSpace.partial_derivative_as_functions import \
    NumericalPartialDerivative_txyz_Functions
import numpy as np


class incompressible_MHD_Base(Base):
    def __init__(self, mesh, Rf, Rm, c):
        """

        Parameters
        ----------
        mesh :
            The status instance.
        Rf :
            The fluid Reynolds number.
        Rm :
            The magnetic Reynolds number.
        c :
            c:=V^2_A/V^2 is the coupling number, where V_A and V are the scales of the
            Alfvén speed and of the flow, respectively.

            c = Al^{-2} where Al is the Alfvén number.

            Al = sqrt(1/c)

        """
        super(incompressible_MHD_Base, self).__init__(mesh)

        self._Rf_ = Rf
        self._Rm_ = Rm

        if isinstance(Rf, (int, float)) and Rf > 0:
            self._1_over_Rf_ = 1 / Rf

        elif Rf == 'infty':
            self._1_over_Rf_ = 0
        else:
            raise Exception(f"Rf={Rf} wrong, be a positive number or 'infty'.")

        if isinstance(Rm, (int, float)) and Rm > 0:
            self._1_over_Rm_ = 1 / Rm
        elif Rm == 'infty':
            self._1_over_Rm_ = 0
        else:
            raise Exception(f"Rm={Rm} wrong, be a positive number or 'infty'.")

        self._c_ = c
        self._Alfven_number_ = np.sqrt(1 / c)

        self._velocity_ = None              # u
        self._vorticity_ = None             # w
        self._pressure_ = None              # p
        self._body_force_ = None            # f

        self._volume_current_density_ = None  # j
        self._magnetic_field_ = None          # B
        self._electric_field_ = None          # E

        self._divergence_of_magnetic_field_ = None

        self._curl_of_magnetic_field_ = None

        self._NPDf_u_ = None
        self._NPDf_v_ = None
        self._NPDf_w_ = None
        self._NPDf_Bx_ = None
        self._NPDf_By_ = None
        self._NPDf_Bz_ = None

        self._NPDf_ux_ = None
        self._NPDf_uy_ = None
        self._NPDf_uz_ = None

        self._NPDf_vx_ = None
        self._NPDf_vy_ = None
        self._NPDf_vz_ = None

        self._NPDf_wx_ = None
        self._NPDf_wy_ = None
        self._NPDf_wz_ = None

        self._NPDf_p_ = None

        self._NPDf_Ex_ = None
        self._NPDf_Ey_ = None
        self._NPDf_Ez_ = None

        self._NPDf_Omx_ = None
        self._NPDf_Omy_ = None
        self._NPDf_Omz_ = None

        self._gradientOfPressure_ = None
        self._totalPressure_ = None              # p + 0.5 * u * u
        self._gradientOfTotalPressure_ = None
        self._kineticEnergyDistribution_ = None
        self._helicityDistribution_ = None
        self._enstrophyDistribution_ = None
        self._divergence_of_velocity_ = None
        self._divergence_of_vorticity_ = None
        self._curl_of_vorticity_ = None
        self._Laplace_of_velocity_ = None

        self._cross_helicity_distribution_ = None
        self._magnetic_helicity_distribution_ = None

        self._mass_source_term_ = None      # s, div(u) = s, should be 0;
        self._magnetic_source_term_ = None  # m, div(B) = m, should be 0;
        self._electric_source_term_ = None  # r, P_t(B) + curl(E) = r, should be 0;

        self._freeze_self_()

    @property
    def Rf(self):
        return self._Rf_

    @property
    def Rm(self):
        return self._Rm_

    @property
    def c(self):
        return self._c_

    @property
    def Alfven_number(self):
        return self._Alfven_number_

    def u(self, t, x, y, z):
        raise NotImplementedError()
    
    def u_t(self, t, x, y, z):
        if self._NPDf_u_ is None:
            self._NPDf_u_ = NumericalPartialDerivative_txyz_Functions(self.u)
        return self._NPDf_u_('t')(t, x, y, z)
    
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
    
    def v_t(self, t, x, y, z):
        if self._NPDf_v_ is None:
            self._NPDf_v_ = NumericalPartialDerivative_txyz_Functions(self.v)
        return self._NPDf_v_('t')(t, x, y, z)
    
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
    
    def w_t(self, t, x, y, z):
        if self._NPDf_w_ is None:
            self._NPDf_w_ = NumericalPartialDerivative_txyz_Functions(self.w)
        return self._NPDf_w_('t')(t, x, y, z)
    
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

    def Bx(self, t, x, y, z):
        raise NotImplementedError()
    
    def Bx_t(self, t, x, y, z):
        if self._NPDf_Bx_ is None:
            self._NPDf_Bx_ = NumericalPartialDerivative_txyz_Functions(self.Bx)
        return self._NPDf_Bx_('t')(t, x, y, z)
    
    def Bx_x(self, t, x, y, z):
        if self._NPDf_Bx_ is None:
            self._NPDf_Bx_ = NumericalPartialDerivative_txyz_Functions(self.Bx)
        return self._NPDf_Bx_('x')(t, x, y, z)
    
    def Bx_y(self, t, x, y, z):
        if self._NPDf_Bx_ is None:
            self._NPDf_Bx_ = NumericalPartialDerivative_txyz_Functions(self.Bx)
        return self._NPDf_Bx_('y')(t, x, y, z)
    
    def Bx_z(self, t, x, y, z):
        if self._NPDf_Bx_ is None:
            self._NPDf_Bx_ = NumericalPartialDerivative_txyz_Functions(self.Bx)
        return self._NPDf_Bx_('z')(t, x, y, z)

    def By(self, t, x, y, z):
        raise NotImplementedError()
    
    def By_t(self, t, x, y, z):
        if self._NPDf_By_ is None:
            self._NPDf_By_ = NumericalPartialDerivative_txyz_Functions(self.By)
        return self._NPDf_By_('t')(t, x, y, z)
    
    def By_x(self, t, x, y, z):
        if self._NPDf_By_ is None:
            self._NPDf_By_ = NumericalPartialDerivative_txyz_Functions(self.By)
        return self._NPDf_By_('x')(t, x, y, z)
    
    def By_y(self, t, x, y, z):
        if self._NPDf_By_ is None:
            self._NPDf_By_ = NumericalPartialDerivative_txyz_Functions(self.By)
        return self._NPDf_By_('y')(t, x, y, z)
    
    def By_z(self, t, x, y, z):
        if self._NPDf_By_ is None:
            self._NPDf_By_ = NumericalPartialDerivative_txyz_Functions(self.By)
        return self._NPDf_By_('z')(t, x, y, z)

    def Bz(self, t, x, y, z):
        raise NotImplementedError()
    
    def Bz_t(self, t, x, y, z):
        if self._NPDf_Bz_ is None:
            self._NPDf_Bz_ = NumericalPartialDerivative_txyz_Functions(self.Bz)
        return self._NPDf_Bz_('t')(t, x, y, z)
    
    def Bz_x(self, t, x, y, z):
        if self._NPDf_Bz_ is None:
            self._NPDf_Bz_ = NumericalPartialDerivative_txyz_Functions(self.Bz)
        return self._NPDf_Bz_('x')(t, x, y, z)
    
    def Bz_y(self, t, x, y, z):
        if self._NPDf_Bz_ is None:
            self._NPDf_Bz_ = NumericalPartialDerivative_txyz_Functions(self.Bz)
        return self._NPDf_Bz_('y')(t, x, y, z)
    
    def Bz_z(self, t, x, y, z):
        if self._NPDf_Bz_ is None:
            self._NPDf_Bz_ = NumericalPartialDerivative_txyz_Functions(self.Bz)
        return self._NPDf_Bz_('z')(t, x, y, z)

    def Ex_x(self, t, x, y, z):
        if self._NPDf_Ex_ is None:
            self._NPDf_Ex_ = NumericalPartialDerivative_txyz_Functions(self.Ex)
        return self._NPDf_Ex_('x')(t, x, y, z)
    
    def Ex_y(self, t, x, y, z):
        if self._NPDf_Ex_ is None:
            self._NPDf_Ex_ = NumericalPartialDerivative_txyz_Functions(self.Ex)
        return self._NPDf_Ex_('y')(t, x, y, z)
    
    def Ex_z(self, t, x, y, z):
        if self._NPDf_Ex_ is None:
            self._NPDf_Ex_ = NumericalPartialDerivative_txyz_Functions(self.Ex)
        return self._NPDf_Ex_('z')(t, x, y, z)

    def Ey_x(self, t, x, y, z):
        if self._NPDf_Ey_ is None:
            self._NPDf_Ey_ = NumericalPartialDerivative_txyz_Functions(self.Ey)
        return self._NPDf_Ey_('x')(t, x, y, z)
    
    def Ey_y(self, t, x, y, z):
        if self._NPDf_Ey_ is None:
            self._NPDf_Ey_ = NumericalPartialDerivative_txyz_Functions(self.Ey)
        return self._NPDf_Ey_('y')(t, x, y, z)
    
    def Ey_z(self, t, x, y, z):
        if self._NPDf_Ey_ is None:
            self._NPDf_Ey_ = NumericalPartialDerivative_txyz_Functions(self.Ey)
        return self._NPDf_Ey_('z')(t, x, y, z)

    def Ez_x(self, t, x, y, z):
        if self._NPDf_Ez_ is None:
            self._NPDf_Ez_ = NumericalPartialDerivative_txyz_Functions(self.Ez)
        return self._NPDf_Ez_('x')(t, x, y, z)
    
    def Ez_y(self, t, x, y, z):
        if self._NPDf_Ez_ is None:
            self._NPDf_Ez_ = NumericalPartialDerivative_txyz_Functions(self.Ez)
        return self._NPDf_Ez_('y')(t, x, y, z)
    
    def Ez_z(self, t, x, y, z):
        if self._NPDf_Ez_ is None:
            self._NPDf_Ez_ = NumericalPartialDerivative_txyz_Functions(self.Ez)
        return self._NPDf_Ez_('z')(t, x, y, z)

    def omega_x_x(self, t, x, y, z):
        if self._NPDf_Omx_ is None:
            self._NPDf_Omx_ = NumericalPartialDerivative_txyz_Functions(self.omega_x)
        return self._NPDf_Omx_('x')(t, x, y, z)

    def omega_x_y(self, t, x, y, z):
        if self._NPDf_Omx_ is None:
            self._NPDf_Omx_ = NumericalPartialDerivative_txyz_Functions(self.omega_x)
        return self._NPDf_Omx_('y')(t, x, y, z)

    def omega_x_z(self, t, x, y, z):
        if self._NPDf_Omx_ is None:
            self._NPDf_Omx_ = NumericalPartialDerivative_txyz_Functions(self.omega_x)
        return self._NPDf_Omx_('z')(t, x, y, z)

    def omega_y_x(self, t, x, y, z):
        if self._NPDf_Omy_ is None:
            self._NPDf_Omy_ = NumericalPartialDerivative_txyz_Functions(self.omega_y)
        return self._NPDf_Omy_('x')(t, x, y, z)

    def omega_y_y(self, t, x, y, z):
        if self._NPDf_Omy_ is None:
            self._NPDf_Omy_ = NumericalPartialDerivative_txyz_Functions(self.omega_y)
        return self._NPDf_Omy_('y')(t, x, y, z)

    def omega_y_z(self, t, x, y, z):
        if self._NPDf_Omy_ is None:
            self._NPDf_Omy_ = NumericalPartialDerivative_txyz_Functions(self.omega_y)
        return self._NPDf_Omy_('z')(t, x, y, z)

    def omega_z_x(self, t, x, y, z):
        if self._NPDf_Omz_ is None:
            self._NPDf_Omz_ = NumericalPartialDerivative_txyz_Functions(self.omega_z)
        return self._NPDf_Omz_('x')(t, x, y, z)

    def omega_z_y(self, t, x, y, z):
        if self._NPDf_Omz_ is None:
            self._NPDf_Omz_ = NumericalPartialDerivative_txyz_Functions(self.omega_z)
        return self._NPDf_Omz_('y')(t, x, y, z)

    def omega_z_z(self, t, x, y, z):
        if self._NPDf_Omz_ is None:
            self._NPDf_Omz_ = NumericalPartialDerivative_txyz_Functions(self.omega_z)
        return self._NPDf_Omz_('z')(t, x, y, z)

    # .............................................................................
    @property
    def magnetic_field(self):
        """B"""
        if self._magnetic_field_ is None:
            self._magnetic_field_ = _3dCSCG_VectorField(
                self.mesh,
                (self.Bx, self.By, self.Bz),
                valid_time=self.valid_time,
            )
        return self._magnetic_field_

    @property
    def divergence_of_magnetic_field(self):
        """div(B)"""
        if self._divergence_of_magnetic_field_ is None:
            # this condition must be valid for all time.
            self._divergence_of_magnetic_field_ = _3dCSCG_ScalarField(
                self.mesh,
                self.___div_of_B___,
                valid_time=self.valid_time,
            )
        return self._divergence_of_magnetic_field_

    def ___div_of_B___(self, t, x, y, z):
        # must be zero
        return self.Bx_x(t, x, y, z) + self.By_y(t, x, y, z) + self.Bz_z(t, x, y, z)
    
    @property
    def volume_current_density(self):
        """j"""
        if self._volume_current_density_ is None:
            self._volume_current_density_ = _3dCSCG_VectorField(
                self.mesh,
                (self.j_x, self.j_y, self.j_z),
                valid_time=self.valid_time,
            )
        return self._volume_current_density_

    def j_x(self, t, x, y, z):
        """ j = curl(B) """
        return self.Bz_y(t, x, y, z) - self.By_z(t, x, y, z)
    
    def j_y(self, t, x, y, z):
        return self.Bx_z(t, x, y, z) - self.Bz_x(t, x, y, z)
    
    def j_z(self, t, x, y, z):
        return self.By_x(t, x, y, z) - self.Bx_y(t, x, y, z)

    @property
    def electric_field(self):
        """E

        E = (1/Rm) * j - u X B
        """
        if self._electric_field_ is None:
            self._electric_field_ = _3dCSCG_VectorField(
                self.mesh,
                (self.Ex, self.Ey, self.Ez),
                valid_time=self.valid_time,
            )
        return self._electric_field_

    def _uXB_x_(self, t, x, y, z):
        """x-component of u X B"""
        return self.v(t, x, y, z) * self.Bz(t, x, y, z) - self.w(t, x, y, z) * self.By(t, x, y, z)

    def _uXB_y_(self, t, x, y, z):
        """y-component of u X B"""
        return self.w(t, x, y, z) * self.Bx(t, x, y, z) - self.u(t, x, y, z) * self.Bz(t, x, y, z)

    def _uXB_z_(self, t, x, y, z):
        """z-component of u X B"""
        return self.u(t, x, y, z) * self.By(t, x, y, z) - self.v(t, x, y, z) * self.Bx(t, x, y, z)

    def Ex(self, t, x, y, z):
        """"""
        return self._1_over_Rm_ * self.j_x(t, x, y, z) - self._uXB_x_(t, x, y, z)
    
    def Ey(self, t, x, y, z):
        """"""
        return self._1_over_Rm_ * self.j_y(t, x, y, z) - self._uXB_y_(t, x, y, z)
    
    def Ez(self, t, x, y, z):
        """"""
        return self._1_over_Rm_ * self.j_z(t, x, y, z) - self._uXB_z_(t, x, y, z)

    @property
    def velocity(self):
        if self._velocity_ is None:
            self._velocity_ = _3dCSCG_VectorField(
                self.mesh,
                (self.u, self.v, self.w),
                valid_time=self.valid_time,
            )
        return self._velocity_

    @property
    def curl_of_velocity(self):
        return self.vorticity

    @property
    def divergence_of_velocity(self):
        if self._divergence_of_velocity_ is None:
            # this condition must be valid for all time.
            self._divergence_of_velocity_ = _3dCSCG_ScalarField(
                self.mesh,
                self.___div_of_velocity___,
                valid_time=None,
            )
        return self._divergence_of_velocity_

    def ___div_of_velocity___(self, t, x, y, z):
        # must be zero
        return self.u_x(t, x, y, z) + self.v_y(t, x, y, z) + self.w_z(t, x, y, z)

    # noinspection PyUnusedLocal
    @staticmethod
    def s(t, x, y, z):
        """The mass source term.

        By default, we have divergence free condition for velocity; the source term is zero.
        But we make define a non-zero source term through method `s` for some particular
        manufactured solutions.
        """
        return 0 * x

    @property
    def mass_source_term(self):
        if self._mass_source_term_ is None:
            self._mass_source_term_ = _3dCSCG_ScalarField(
                self.mesh,
                self.s,
                valid_time=self.valid_time,
                name='mass-source-term')
        return self._mass_source_term_

    # noinspection PyUnusedLocal
    @staticmethod
    def m(t, x, y, z):
        """the magnetic source term. So, div(B) = m"""
        return 0 * x

    @property
    def magnetic_source_term(self):
        if self._magnetic_source_term_ is None:
            self._magnetic_source_term_ = _3dCSCG_ScalarField(
                self.mesh,
                self.m,
                valid_time=self.valid_time,
                name='magnetic-source-term',
            )
        return self._magnetic_source_term_

    @property
    def vorticity(self):
        if self._vorticity_ is None:
            self._vorticity_ = _3dCSCG_VectorField(
                self.mesh,
                (self.omega_x, self.omega_y, self.omega_z),
                valid_time=self.valid_time,
            )
        return self._vorticity_

    def omega_x(self, t, x, y, z):
        return self.w_y(t, x, y, z) - self.v_z(t, x, y, z)
    
    def omega_y(self, t, x, y, z):
        return self.u_z(t, x, y, z) - self.w_x(t, x, y, z)
    
    def omega_z(self, t, x, y, z):
        return self.v_x(t, x, y, z) - self.u_y(t, x, y, z)


    def ___curl_omega_x___(self, t, x, y, z):
        return self.omega_z_y(t, x, y, z) - self.omega_y_z(t, x, y, z)

    def ___curl_omega_y___(self, t, x, y, z):
        return self.omega_x_z(t, x, y, z) - self.omega_z_x(t, x, y, z)

    def ___curl_omega_z___(self, t, x, y, z):
        return self.omega_y_x(t, x, y, z) - self.omega_x_y(t, x, y, z)


    def _CurlBx_(self, t, x, y, z):
        return self.Bz_y(t, x, y, z) - self.By_z(t, x, y, z)
    
    def _CurlBy_(self, t, x, y, z):
        return self.Bx_z(t, x, y, z) - self.Bz_x(t, x, y, z)
    
    def _CurlBz_(self, t, x, y, z):
        return self.By_x(t, x, y, z) - self.Bx_y(t, x, y, z)

    def _CurlEx_(self, t, x, y, z):
        return self.Ez_y(t, x, y, z) - self.Ey_z(t, x, y, z)
    
    def _CurlEy_(self, t, x, y, z):
        return self.Ex_z(t, x, y, z) - self.Ez_x(t, x, y, z)
    
    def _CurlEz_(self, t, x, y, z):
        return self.Ey_x(t, x, y, z) - self.Ex_y(t, x, y, z)

    @property
    def curl_of_magnetic_field(self):
        if self._curl_of_magnetic_field_ is None:
            self._curl_of_magnetic_field_ = _3dCSCG_VectorField(
                self.mesh,
                (self._CurlBx_, self._CurlBy_, self._CurlBz_),
                valid_time=self.valid_time)
        return self._curl_of_magnetic_field_

    @property
    def curl_of_vorticity(self):
        """curl_of_vorticity = curl of curl u = grad(div u) - Vector_Laplace u.

        If div(u) = zero (conservation of mass), we have

        curl_of_vorticity = - Vector_Laplace u = - ( Laplace u_x,  Laplace u_y,  Laplace u_z) ^ T.

        """
        if self._curl_of_vorticity_ is None:
            self._curl_of_vorticity_ = _3dCSCG_VectorField(
                self.mesh,
                (self.___curl_omega_x___, self.___curl_omega_y___, self.___curl_omega_z___),
                valid_time=self.valid_time,
            )
        return self._curl_of_vorticity_

    @property
    def Laplace_of_velocity(self):
        if self._Laplace_of_velocity_ is None:
            self._Laplace_of_velocity_ = _3dCSCG_VectorField(
                self.mesh,
                (self.___laplace_u___, self.___laplace_v___, self.___laplace_w___),
                valid_time=self.valid_time,
            )
        return self._Laplace_of_velocity_

    def ___laplace_u___(self, t, x, y, z):
        return self.u_xx(t, x, y, z) + self.u_yy(t, x, y, z) + self.u_zz(t, x, y, z)
    
    def ___laplace_v___(self, t, x, y, z):
        return self.v_xx(t, x, y, z) + self.v_yy(t, x, y, z) + self.v_zz(t, x, y, z)
    
    def ___laplace_w___(self, t, x, y, z):
        return self.w_xx(t, x, y, z) + self.w_yy(t, x, y, z) + self.w_zz(t, x, y, z)

    @property
    def divergence_of_vorticity(self):
        if self._divergence_of_vorticity_ is None:
            # this condition must be valid for all time.
            self._divergence_of_vorticity_ = _3dCSCG_ScalarField(
                self.mesh,
                self.___div_of_vorticity___,
                valid_time=None,
            )
        return self._divergence_of_vorticity_

    def ___div_of_vorticity___(self, t, x, y, z):
        # must be zero
        return self.omega_x_x(t, x, y, z) + self.omega_y_y(t, x, y, z) + self.omega_z_z(t, x, y, z)

    def fx(self, t, x, y, z):
        return self.u_t(t, x, y, z) \
                    + self.omega_y(t, x, y, z) * self.w(t, x, y, z) - self.omega_z(t, x, y, z) * self.v(t, x, y, z) \
                    + self._tp_x_(t, x, y, z) \
                    - self._1_over_Rf_ * self.___laplace_u___(t, x, y, z) \
                    - self.c * (self.j_y(t, x, y, z) * self.Bz(t, x, y, z) - self.j_z(t, x, y, z) * self.By(t, x, y, z))
    
    def fy(self, t, x, y, z):
        return self.v_t(t, x, y, z) \
                    + self.omega_z(t, x, y, z) * self.u(t, x, y, z) - self.omega_x(t, x, y, z) * self.w(t, x, y, z) \
                    + self._tp_y_(t, x, y, z) \
                    - self._1_over_Rf_ * self.___laplace_v___(t, x, y, z) \
                    - self.c * (self.j_z(t, x, y, z) * self.Bx(t, x, y, z) - self.j_x(t, x, y, z) * self.Bz(t, x, y, z))
    
    def fz(self, t, x, y, z):
        return self.w_t(t, x, y, z) \
                    + self.omega_x(t, x, y, z) * self.v(t, x, y, z) - self.omega_y(t, x, y, z) * self.u(t, x, y, z) \
                    + self._tp_z_(t, x, y, z) \
                    - self._1_over_Rf_ * self.___laplace_w___(t, x, y, z) \
                    - self.c * (self.j_x(t, x, y, z) * self.By(t, x, y, z) - self.j_y(t, x, y, z) * self.Bx(t, x, y, z))

    @property
    def body_force(self):
        """In fact, it is sort of momentum source term."""
        if self._body_force_ is None:
            self._body_force_ = _3dCSCG_VectorField(
                self.mesh,
                (self.fx, self.fy, self.fz),
                valid_time=self.valid_time)
        return self._body_force_

    def rx(self, t, x, y, z):
        return self.Bx_t(t, x, y, z) + self._CurlEx_(t, x, y, z)

    def ry(self, t, x, y, z):
        return self.By_t(t, x, y, z) + self._CurlEy_(t, x, y, z)

    def rz(self, t, x, y, z):
        return self.Bz_t(t, x, y, z) + self._CurlEz_(t, x, y, z)

    @property
    def electric_source_term(self):
        """r; partial_t(B) + curl(E) = r"""
        if self._electric_source_term_ is None:

            self._electric_source_term_ = _3dCSCG_VectorField(
                self.mesh,
                (self.rx, self.ry, self.rz),
                valid_time=self.valid_time,
            )
        return self._electric_source_term_

    @property
    def pressure(self):
        if self._pressure_ is None:
            self._pressure_ = _3dCSCG_ScalarField(
                self.mesh,
                self.p,
                valid_time=self.valid_time,
            )
        return self._pressure_

    @property
    def gradient_of_pressure(self):
        if self._gradientOfPressure_ is None:
            self._gradientOfPressure_ = _3dCSCG_VectorField(
                self.mesh,
                (self.p_x, self.p_y, self.p_z),
                valid_time=self.valid_time,
            )
        return self._gradientOfPressure_

    @property
    def total_pressure(self):
        if self._totalPressure_ is None:
            self._totalPressure_ = _3dCSCG_ScalarField(
                self.mesh,
                self.___total_pressure___,
                valid_time=self.valid_time,
            )
        return self._totalPressure_
    
    def ___total_pressure___(self, t, x, y, z):
        return self.p(t, x, y, z) + self.___kinetic_energy_distribution___(t, x, y, z)
    
    def _tp_x_(self, t, x, y, z):
        """ """
        return self.p_x(t, x, y, z) + self._ke_x_(t, x, y, z)
    
    def _tp_y_(self, t, x, y, z):
        """ """
        return self.p_y(t, x, y, z) + self._ke_y_(t, x, y, z)
    
    def _tp_z_(self, t, x, y, z):
        """ """
        return self.p_z(t, x, y, z) + self._ke_z_(t, x, y, z)

    @property
    def gradient_of_total_pressure(self):
        """A vector field of the total pressure distribution."""
        if self._gradientOfTotalPressure_ is None:
            self._gradientOfTotalPressure_ = _3dCSCG_VectorField(
                self.mesh,
                (self._tp_x_, self._tp_y_, self._tp_z_),
                valid_time=self.valid_time
            )
        return self._gradientOfTotalPressure_

    @property
    def kinetic_energy_distribution(self):
        """A scalar field of the kinetic energy distribution."""
        if self._kineticEnergyDistribution_ is None:
            self._kineticEnergyDistribution_ = _3dCSCG_ScalarField(
                self.mesh,
                self.___kinetic_energy_distribution___,
                valid_time=self.valid_time)
        return self._kineticEnergyDistribution_
    
    def ___kinetic_energy_distribution___(self, t, x, y, z):
        return 0.5 * (self.u(t, x, y, z)**2 + self.v(t, x, y, z)**2 + self.w(t, x, y, z)**2)
    
    def _ke_x_(self, t, x, y, z):
        """ d(kinetic_energy)/dx """
        return self.u(t, x, y, z) * self.u_x(t, x, y, z) + \
            self.v(t, x, y, z) * self.v_x(t, x, y, z) + \
            self.w(t, x, y, z) * self.w_x(t, x, y, z)
    
    def _ke_y_(self, t, x, y, z):
        """ d(kinetic_energy)/dy """
        return self.u(t, x, y, z) * self.u_y(t, x, y, z) + \
            self.v(t, x, y, z) * self.v_y(t, x, y, z) + \
            self.w(t, x, y, z) * self.w_y(t, x, y, z)
    
    def _ke_z_(self, t, x, y, z):
        """ d(kinetic_energy)/dz """
        return self.u(t, x, y, z) * self.u_z(t, x, y, z) + \
            self.v(t, x, y, z) * self.v_z(t, x, y, z) + \
            self.w(t, x, y, z) * self.w_z(t, x, y, z)

    @property
    def fluid_helicity_distribution(self):
        """A scalar field of the helicity distribution."""
        if self._helicityDistribution_ is None:
            self._helicityDistribution_ = _3dCSCG_ScalarField(
                self.mesh,
                self.___fluid_helicity_distribution___,
                valid_time=self.valid_time)
        return self._helicityDistribution_

    def ___fluid_helicity_distribution___(self, t, x, y, z):
        return self.u(t, x, y, z) * self.omega_x(t, x, y, z) + \
               self.v(t, x, y, z) * self.omega_y(t, x, y, z) + \
               self.w(t, x, y, z) * self.omega_z(t, x, y, z)

    @property
    def enstrophy_distribution(self):
        """A scalar field of the enstrophy distribution."""
        if self._enstrophyDistribution_ is None:
            self._enstrophyDistribution_ = _3dCSCG_ScalarField(
                self.mesh,
                self.___enstrophy_distribution___,
                valid_time=self.valid_time)
        return self._enstrophyDistribution_

    def ___enstrophy_distribution___(self, t, x, y, z):
        return 0.5 * (self.omega_x(t, x, y, z) ** 2 +
                      self.omega_y(t, x, y, z) ** 2 +
                      self.omega_z(t, x, y, z) ** 2)

    @property
    def cross_helicity_distribution(self):
        if self._cross_helicity_distribution_ is None:
            self._cross_helicity_distribution_ = \
                self.velocity.do.inner_product(self.magnetic_field)
        return self._cross_helicity_distribution_

    @lru_cache(maxsize=8)
    def fluid_helicity(self, t):
        """Helicity at time `t`."""
        return self.___Pr_compute_Ln_norm_of___('fluid_helicity_distribution', time=t, n=1)

    @lru_cache(maxsize=8)
    def kinetic_energy(self, t):
        """Kinetic energy at time `t`."""
        return self.___Pr_compute_Ln_norm_of___('kinetic_energy_distribution', time=t, n=1)

    @lru_cache(maxsize=8)
    def enstrophy(self, t):
        """Enstrophy at time `t`."""
        return self.___Pr_compute_Ln_norm_of___('enstrophy_distribution', time=t, n=1)

    def ___PreFrozenChecker___(self):
        """Will be called before freezing self."""
        TS = self.___PRIVATE_generate_random_valid_time_instances___()
        x, y, z = self._mesh_.do.generate_random_coordinates()

        rst = (x, y, z)

        if len(x) == 0:
            return

        for t in TS:

            assert np.all(np.isclose(
                self.___div_of_velocity___(t, *rst),
                self.s(t, *rst))), \
                " <exact solution> <MHD> : velocity divergence free condition."

            try:
                fx = self.u_t(t, x, y, z) + \
                    self.omega_y(t, x, y, z) * self.w(t, x, y, z) - self.omega_z(t, x, y, z) * self.v(t, x, y, z) + \
                    self._tp_x_(t, x, y, z) - \
                    self._1_over_Rf_ * (self.u_xx(t, x, y, z) + self.u_yy(t, x, y, z) + self.u_zz(t, x, y, z)) - \
                    self.c * (self.j_y(t, x, y, z) * self.Bz(t, x, y, z) - self.j_z(t, x, y, z) * self.By(t, x, y, z))

                fy = self.v_t(t, x, y, z) + \
                    self.omega_z(t, x, y, z) * self.u(t, x, y, z) - self.omega_x(t, x, y, z) * self.w(t, x, y, z) + \
                    self._tp_y_(t, x, y, z) - \
                    self._1_over_Rf_ * (self.v_xx(t, x, y, z) + self.v_yy(t, x, y, z) + self.v_zz(t, x, y, z)) - \
                    self.c * (self.j_z(t, x, y, z) * self.Bx(t, x, y, z) - self.j_x(t, x, y, z) * self.Bz(t, x, y, z))

                fz = self.w_t(t, x, y, z) + \
                    self.omega_x(t, x, y, z) * self.v(t, x, y, z) - self.omega_y(t, x, y, z) * self.u(t, x, y, z) + \
                    self._tp_z_(t, x, y, z) - \
                    self._1_over_Rf_ * (self.w_xx(t, x, y, z) + self.w_yy(t, x, y, z) + self.w_zz(t, x, y, z)) - \
                    self.c * (self.j_x(t, x, y, z) * self.By(t, x, y, z) - self.j_y(t, x, y, z) * self.Bx(t, x, y, z))

                FX = self.fx(t, x, y, z)
                FY = self.fy(t, x, y, z)
                FZ = self.fz(t, x, y, z)

                np.testing.assert_array_almost_equal(FX - fx, 0, decimal=5)
                np.testing.assert_array_almost_equal(FY - fy, 0, decimal=5)
                np.testing.assert_array_almost_equal(FZ - fz, 0, decimal=5)

                rx = self.Bx_t(t, x, y, z) + self._CurlEx_(t, x, y, z)
                ry = self.By_t(t, x, y, z) + self._CurlEy_(t, x, y, z)
                rz = self.Bz_t(t, x, y, z) + self._CurlEz_(t, x, y, z)

                RX = self.rx(t, x, y, z)
                RY = self.ry(t, x, y, z)
                RZ = self.rz(t, x, y, z)

                np.testing.assert_array_almost_equal(RX, rx, decimal=4)
                np.testing.assert_array_almost_equal(RY, ry, decimal=4)
                np.testing.assert_array_almost_equal(RZ, rz, decimal=4)

            except NotImplementedError:
                pass

            time = t
            try:
                Pu = NumericalPartialDerivative_txyz(self.u, t, x, y, z)
                assert all(Pu.check_total(self.u_t, self.u_x, self.u_y, self.u_z))
            except NotImplementedError:
                pass

            try:
                Pu = NumericalPartialDerivative_txyz(self.v, t, x, y, z)
                assert all(Pu.check_total(self.v_t, self.v_x, self.v_y, self.v_z))
            except NotImplementedError:
                pass
            try:
                Pu = NumericalPartialDerivative_txyz(self.w, t, x, y, z)
                assert all(Pu.check_total(self.w_t, self.w_x, self.w_y, self.w_z))
            except NotImplementedError:
                pass

            try:
                PB = NumericalPartialDerivative_txyz(self.Bx, t, x, y, z)
                assert all(PB.check_total(self.Bx_t, self.Bx_x, self.Bx_y, self.Bx_z))
            except NotImplementedError:
                pass
            try:
                PB = NumericalPartialDerivative_txyz(self.By, t, x, y, z)
                assert all(PB.check_total(self.By_t, self.By_x, self.By_y, self.By_z))
            except NotImplementedError:
                pass
            try:
                PB = NumericalPartialDerivative_txyz(self.Bz, t, x, y, z)
                assert all(PB.check_total(self.Bz_t, self.Bz_x, self.Bz_y, self.Bz_z))
            except NotImplementedError:
                pass

            assert np.all(np.isclose(
                self.u_x(t, *rst) + self.v_y(t, *rst) + self.w_z(t, *rst),
                self.s(time, *rst))), \
                " <exact solution> <MHD> : velocity divergence condition."

            assert np.all(np.isclose(
                self.___div_of_B___(time, *rst),
                self.m(time, *rst))), \
                " <exact solution> <MHD> : magnetic field divergence free condition."

            assert np.all(np.isclose(
                self.___div_of_vorticity___(time, *rst),
                0)), \
                " <exact solution> <MHD> : magnetic field divergence free condition."

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

            try:
                p = partial(self.p, time)
                p_x = partial(self.p_x, time)
                p_y = partial(self.p_y, time)
                p_z = partial(self.p_z, time)
                Pp = NumericalPartialDerivative_xyz(p, *rst)
                assert all(Pp.check_total(p_x, p_y, p_z))
            except NotImplementedError:
                pass

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
