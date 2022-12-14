# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly

class MassMatrixHelper(FrozenOnly):
    """The class for the inner product matrix."""
    def __init__(self, sf):
        assert sf.ndim == 2, " <___STORAGE_OPERATORS_INNER___> "
        self._mesh_ = sf.mesh
        self._sf_ = sf

        quad_nodes, _, quad_weights = sf.space.___PRIVATE_do_evaluate_quadrature___(
            sf.dqp,
            quad_type='Gauss',
        )
        xietasigma, bfSelf = sf.do.evaluate_basis_at_meshgrid(*quad_nodes)
        self._O_xietasigma_ = xietasigma
        self._O_quad_weights_ = quad_weights
        self._O_bfSelf_ = bfSelf

        if not self._mesh_.whether.orthogonal:
            quad_nodes, _, quad_weights = sf.space.___PRIVATE_do_evaluate_quadrature___(
                [sf.dqp[i] + 1 for i in range(2)],
                quad_type='Gauss'
            )
            xietasigma, bfSelf = sf.do.evaluate_basis_at_meshgrid(*quad_nodes)
            self._C_xietasigma_ = xietasigma
            self._C_quad_weights_ = quad_weights
            self._C_bfSelf_ = bfSelf
        else:
            pass

        self._freeze_self_()

    def __call__(self, i):
        element = self._mesh_.elements[i]
        mark = element.type_wrt_metric.mark

        if isinstance(mark, str) and mark[:4] == 'Orth':

            Mi = self._sf_.___PRIVATE_operator_inner___(
                self._sf_,
                i,
                self._O_xietasigma_,
                self._O_quad_weights_,
                self._O_bfSelf_,
                self._O_bfSelf_,
            )

        else:

            Mi = self._sf_.___PRIVATE_operator_inner___(
                self._sf_,
                i,
                self._C_xietasigma_,
                self._C_quad_weights_,
                self._C_bfSelf_,
                self._C_bfSelf_,
            )

        return Mi