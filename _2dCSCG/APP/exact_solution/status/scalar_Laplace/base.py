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
from functools import partial
from _2dCSCG.APP.exact_solution.status.base import Base
from screws.numerical._2d_space.Jacobian_21 import NumericalPartialDerivative_xy
from _2dCSCG.fields.scalar.main import _2dCSCG_ScalarField
from _2dCSCG.fields.vector.main import _2dCSCG_VectorField



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
        self.___PRIVATE_check_self___()
        self._freeze_self_()

    # to be overridden (must) ...

    def p(self, t, x, y): raise NotImplementedError()
    def p_x(self, t, x, y): raise NotImplementedError()
    def p_xx(self, t, x, y): raise NotImplementedError()
    def p_y(self, t, x, y): raise NotImplementedError()
    def p_yy(self, t, x, y): raise NotImplementedError()

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

    def ___source___(self, t, x, y):
        """By default, we have divergence free condition; the source term is zero."""
        return self.p_xx(t, x, y) + self.p_yy(t, x, y)
    @property
    def source(self):
        if self._source_ is None:
            self._source_ = _2dCSCG_ScalarField(self.mesh, self.___source___, valid_time=self.valid_time)
        return self._source_



    @property
    def ___check_domain___(self):
        """
        We use this general domain to do the check, in particular exact solution, we can define new domain
        by override this method.

        """
        r = np.linspace(random.uniform(-1, -0.9), random.uniform(0.82, 1), random.randint(4, 6))
        s = np.linspace(random.uniform(-0.95, -0.88), random.uniform(0.75, 1), random.randint(3, 5))
        r, s = np.meshgrid(r, s, indexing='ij')
        return r, s

    def ___PRIVATE_check_self___(self):
        """
        We use this general method to do the check, in particular exact solution, we can define particular
        check method by override this method.
        """
        time = 0
        rs = self. ___check_domain___
        p = partial(self.p, time)
        p_x = partial(self.p_x, time)
        p_y = partial(self.p_y, time)
        Pp = NumericalPartialDerivative_xy(p, *rs)
        assert all(Pp.check_total(p_x, p_y))
        p_xx = partial(self.p_xx, time)
        p_yy = partial(self.p_yy, time)
        Ppx = NumericalPartialDerivative_xy(p_x, *rs)
        Ppy = NumericalPartialDerivative_xy(p_y, *rs)
        assert Ppx.check_partial_x(p_xx)
        assert Ppy.check_partial_y(p_yy)