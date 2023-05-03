# -*- coding: utf-8 -*-
"""The 3d wave problem

The wave equation: \partial^2 / \partial t^2 - Laplace \phi = f

And in the mixed formulation:

\partial \phi / \partial t = \psi
\partial \psi / \partial t  - div u = f
u = grad \phi

"""

import numpy as np
from objects.CSCG._3d.exactSolutions.base import Base
from objects.CSCG._3d.fields.vector.main import _3dCSCG_VectorField
from objects.CSCG._3d.fields.scalar.main import _3dCSCG_ScalarField

from components.numerical.timePlus3dSpace.partial_derivative import NumericalPartialDerivative_txyz
from components.numerical.timePlus3dSpace.partial_derivative_as_functions import \
    NumericalPartialDerivative_txyz_Functions


class Wave_Base(Base):
    def __init__(self, mesh):
        super(Wave_Base, self).__init__(mesh)
        self._wave_function_ = None
        self._velocity_ = None
        self._source_term_ = None
        self._potential_ = None

        self._NPDf_p_ = None
        self._NPDf_px_ = None
        self._NPDf_py_ = None
        self._NPDf_pz_ = None

        self._NPDf_psi_ = None

        self._freeze_self_()


    def phi(self, t, x, y, z):
        raise NotImplementedError()

    def psi(self, t, x, y, z):
        """phi_t"""
        if self._NPDf_p_ is None:
            self._NPDf_p_ = NumericalPartialDerivative_txyz_Functions(self.phi)
        return self._NPDf_p_('t')(t, x, y, z)

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

    def psi_t(self, t, x, y, z):
        """psi_t"""
        if self._NPDf_psi_ is None:
            self._NPDf_psi_ = NumericalPartialDerivative_txyz_Functions(self.psi)
        return self._NPDf_psi_('t')(t, x, y, z)

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
        return self.psi_t(t, x, y, z) - self.u_x(t, x, y, z) - self.v_y(t, x, y, z) - self.w_z(t, x, y, z)


    @property
    def wave_function(self):
        if self._wave_function_ is None:
            self._wave_function_ = _3dCSCG_ScalarField(
                self.mesh,
                self.phi,
                valid_time=self.valid_time,
                name='wave equation')
        return self._wave_function_

    @property
    def potential(self):
        if self._potential_ is None:
            self._potential_ = _3dCSCG_ScalarField(
                self.mesh,
                self.psi,
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

    def ___PreFrozenChecker___(self):
        """
        We use this general method to do the check, in particular exact solution, we can define particular
        check method by override this method.
        """
        TS = self.___PRIVATE_generate_random_valid_time_instances___()
        x, y, z = self._mesh_.do.generate_random_coordinates()

        if len(x) == 0:
            return

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
                f = self.psi_t(t, x, y, z) - self.u_x(t, x, y, z) - self.v_y(t, x, y, z) - self.w_z(t, x, y, z)
                F = self.f(t, x, y, z)
                np.testing.assert_array_almost_equal(F - f, 0, decimal=5)
            except NotImplementedError:
                pass
