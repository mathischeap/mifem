# -*- coding: utf-8 -*-
"""
The time independent Schrödinger equation.

    - (bh / 2m) nabla^2 psi(x,y,z) + V(x,y,z) psi(x,y,z) = E psi(x,y,z)

psi: the time-independent wave equation.
bh : reduced Planck constant:  bh = h / (2pi)
m : mass
V(x, y, z): potential field
E : energy constant


If we introduce u = gradient psi, we will have

u = grad psi
- (bh / 2m) div u + V(x,y,z) psi(x,y,z) - E psi(x,y,z) = 0

"""
import numpy as np
from scipy.constants import Planck, hbar
from objects.CSCG._3d.exactSolutions.base import Base
from objects.CSCG._3d.fields.vector.main import _3dCSCG_VectorField
from objects.CSCG._3d.fields.scalar.main import _3dCSCG_ScalarField

from components.numerical.timePlus3dSpace.partial_derivative import NumericalPartialDerivative_txyz
from components.numerical.timePlus3dSpace.partial_derivative_as_functions import \
    NumericalPartialDerivative_txyz_Functions

class TimeIndependentSchrodingerEquationBase(Base):
    def __init__(self, mesh, m, E):
        super(TimeIndependentSchrodingerEquationBase, self).__init__(mesh)

        self._m_ = m # mass
        self._E_ = E # energy constant

        self._h_ = Planck # Planck constant.
        self._h_bar_ = hbar # reduced Planck constant.

        self._alpha_ = hbar ** 2 / (2 * m)

        self._V_ = None
        self._wave_function_ = None
        self._flux_ = None
        self._source_term_ = None

        self._NPDf_p_ = None
        self._NPDf_px_ = None
        self._NPDf_py_ = None
        self._NPDf_pz_ = None
        self._freeze_self_()

    @property
    def h(self):
        """Planck constant"""
        return self._h_

    @property
    def hbar(self):
        """reduced Planck constant."""
        return self._h_bar_

    @property
    def m(self):
        """Mass"""
        return self._m_

    @property
    def E(self):
        """Energy constant."""
        return self._E_

    def V(self, t, x, y, z):
        """Potential field"""
        raise NotImplementedError()

    def psi(self, t, x, y, z):
        raise NotImplementedError()

    def u(self, t, x, y, z):
        """phi_x"""
        if self._NPDf_p_ is None:
            self._NPDf_p_ = NumericalPartialDerivative_txyz_Functions(self.psi)
        return self._NPDf_p_('x')(t, x, y, z)
    def v(self, t, x, y, z):
        """phi_y"""
        if self._NPDf_p_ is None:
            self._NPDf_p_ = NumericalPartialDerivative_txyz_Functions(self.psi)
        return self._NPDf_p_('y')(t, x, y, z)
    def w(self, t, x, y, z):
        """phi_z"""
        if self._NPDf_p_ is None:
            self._NPDf_p_ = NumericalPartialDerivative_txyz_Functions(self.psi)
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


    @property
    def potential(self):
        if self._V_ is None:
            self._V_ = _3dCSCG_ScalarField(
                self.mesh,
                self.V,
                valid_time=self.valid_time,
                name='potential')
        return self._V_


    @property
    def wave_function(self):
        if self._wave_function_ is None:
            self._wave_function_ = _3dCSCG_ScalarField(
                self.mesh,
                self.psi,
                valid_time=self.valid_time,
                name='wave function')
        return self._wave_function_

    @property
    def flux(self):
        if self._flux_ is None:
            self._flux_ = _3dCSCG_VectorField(
                self.mesh,
                (self.u, self.v, self.w),
                valid_time=self.valid_time,
                name='flux')
        return self._flux_



    def ___source___(self, t, x, y, z):
        """"""
        P = self.psi(t, x, y, z)
        EP = self.E * P
        VP = self.V(t, x, y, z) * P
        return (VP - EP) / self._alpha_


    @property
    def source_term(self):
        if self._source_term_ is None:
            self._source_term_ = _3dCSCG_ScalarField(
                self.mesh,
                self.___source___,
                valid_time=self.valid_time,
                name='wave function')
        return self._source_term_


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
                Pu = NumericalPartialDerivative_txyz(self.psi, t, x, y, z)
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
                laplacian_psi = self.u_x(t, x, y, z) + self.v_y(t, x, y, z) + self.w_z(t, x, y, z)

                left = - self._alpha_ * laplacian_psi + self.V(t, x, y, z) * self.psi(t, x, y, z)

                right = self.E * self.psi(t, x, y, z)

                np.testing.assert_array_almost_equal(left, right, decimal=3)
            except NotImplementedError:
                pass