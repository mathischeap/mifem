# -*- coding: utf-8 -*-
"""
The parent of all exact solutions of scalar laplacian problem. Note that this is not a Poisson problem

source = Laplacian of (potential)

f = div grad (p)

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import numpy as np
from _2dCSCG.APP.exact_solution.status.base import Base
from _2dCSCG.fields.scalar.main import _2dCSCG_ScalarField
from _2dCSCG.fields.vector.main import _2dCSCG_VectorField

from screws.numerical.time_plus_2d_space.partial_derivative import NumericalPartialDerivative_txy
from screws.numerical.time_plus_2d_space.partial_derivative_as_functions import \
    NumericalPartialDerivative_txy_Functions


class scalar_Laplace_Base(Base):
    """ Given a scalar field: potential,  source = laplace potential.

    or in mixed formulation:

        velocity = gradient potential
        source = divergence velocity

    """
    def __init__(self, es):
        super(scalar_Laplace_Base, self).__init__(es)
        self._potential_ = None
        self._velocity_ = None
        self._source_ = None

        self._NPDf_p_ = None
        self._NPDf_px_ = None
        self._NPDf_py_ = None
        self._freeze_self_()


    def p(self, t, x, y): raise NotImplementedError()
    def p_x(self, t, x, y):
        if self._NPDf_p_ is None:
            self._NPDf_p_ = NumericalPartialDerivative_txy_Functions(self.p)
        return self._NPDf_p_('x')(t, x, y)
    def p_y(self, t, x, y):
        if self._NPDf_p_ is None:
            self._NPDf_p_ = NumericalPartialDerivative_txy_Functions(self.p)
        return self._NPDf_p_('y')(t, x, y)

    def p_xx(self, t, x, y):
        if self._NPDf_px_ is None:
            self._NPDf_px_ = NumericalPartialDerivative_txy_Functions(self.p_x)
        return self._NPDf_px_('x')(t, x, y)
    def p_yy(self, t, x, y):
        if self._NPDf_py_ is None:
            self._NPDf_py_ = NumericalPartialDerivative_txy_Functions(self.p_y)
        return self._NPDf_py_('y')(t, x, y)

    # .............................................................................

    @property
    def potential(self):
        if self._potential_ is None:
            self._potential_ = _2dCSCG_ScalarField(self.mesh, self.p, valid_time=self.valid_time)
        return self._potential_

    @property
    def velocity(self):
        """"""
        if self._velocity_ is None:
            self._velocity_ = _2dCSCG_VectorField(self.mesh, (self.p_x, self.p_y), valid_time=self.valid_time)
        return self._velocity_

    def f(self, t, x, y):
        """"""
        return self.p_xx(t, x, y) + self.p_yy(t, x, y)

    @property
    def source_term(self):
        if self._source_ is None:
            self._source_ = _2dCSCG_ScalarField(self.mesh,
                                                self.f,
                                                valid_time=self.valid_time,
                                                name='source_term')
        return self._source_



    def ___PreFrozenChecker___(self):
        """
        We use this general method to do the check, in particular exact solution, we can define particular
        check method by override this method.
        """
        TS = self.___PRIVATE_generate_random_valid_time_instances___()
        x, y = self._mesh_.do.generate_random_coordinates()

        if len(x) == 0: return

        for t in TS:

            t = float(t)

            try:
                Pu = NumericalPartialDerivative_txy(self.p, t, x, y)
                assert Pu.check_partial_x(self.p_x)
                assert Pu.check_partial_y(self.p_y)
            except NotImplementedError:
                pass

            try:
                Pu = NumericalPartialDerivative_txy(self.p_x, t, x, y)
                assert Pu.check_partial_x(self.p_xx)
            except NotImplementedError:
                pass

            try:
                Pu = NumericalPartialDerivative_txy(self.p_y, t, x, y)
                assert Pu.check_partial_y(self.p_yy)
            except NotImplementedError:
                pass

            try:
                f = self.p_xx(t, x, y) + self.p_yy(t, x, y)
                F = self.f(t, x, y)
                np.testing.assert_array_almost_equal(F - f, 0, decimal=5)
            except NotImplementedError:
                pass