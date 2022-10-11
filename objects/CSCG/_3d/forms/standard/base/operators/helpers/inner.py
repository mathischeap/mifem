# -*- coding: utf-8 -*-

import numpy as np
from screws.freeze.base import FrozenOnly



class ___Operators_3dCSCG_sf_Inner___(FrozenOnly):
    """The class for the inner product matrix."""
    def __init__(self, sf, of, quad_degree=None, quad_type='Lobatto'):
        assert sf.ndim == of.ndim and sf.k == of.k, " <___STORAGE_OPERATORS_INNER___> "
        assert sf.mesh == of.mesh, "Meshes do not match."
        self._mesh_ = sf.mesh
        self._sf_ = sf
        self._of_ = of
        if quad_degree is None:
            quad_degree = [int(np.max([sf.dqp[i], of.dqp[i]])) for i in range(3)]
        quad_nodes, _, quad_weights = sf.space.___PRIVATE_do_evaluate_quadrature___(
            quad_degree, quad_type=quad_type)
        xietasigma, bfSelf = sf.do.evaluate_basis_at_meshgrid(*quad_nodes)
        if of is sf:
            bfOther = bfSelf
        else:
            xietasigma, bfOther = of.do.evaluate_basis_at_meshgrid(*quad_nodes)
        self._xietasigma_ = xietasigma
        self._quad_weights_ = quad_weights
        self._bfSelf_ = bfSelf
        self._bfOther_ = bfOther
        self._freeze_self_()

    def __call__(self, i):
        Mi = self._sf_.___PRIVATE_operator_inner___(
            self._of_, i, self._xietasigma_, self._quad_weights_, self._bfSelf_, self._bfOther_
        )
        return Mi
