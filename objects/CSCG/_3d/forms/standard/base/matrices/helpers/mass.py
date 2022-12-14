# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 12/12/2022 11:14 PM
"""
from components.freeze.main import FrozenOnly


class MassMatrixHelper(FrozenOnly):
    """The class for the inner product matrix."""
    def __init__(self, sf):
        assert sf.ndim == 3, " <MassMatrixHelper> "
        self._mesh_ = sf.mesh
        self._sf_ = sf

        quad_nodes, _, quad_weights = sf.space.___PRIVATE_do_evaluate_quadrature___(
            sf.dqp,
            quad_type='Gauss',
        )
        xietasigma, bfSelf = sf.do.evaluate_basis_at_meshgrid(*quad_nodes)
        self._Orth_xietasigma_ = xietasigma
        self._Orth_quad_weights_ = quad_weights
        self._Orth_bfSelf_ = bfSelf

        if not self._mesh_.whether.orthogonal:
            quad_nodes, _, quad_weights = sf.space.___PRIVATE_do_evaluate_quadrature___(
                [sf.dqp[i] + 1 for i in range(3)],
                quad_type='Gauss',
            )
            xietasigma, bfSelf = sf.do.evaluate_basis_at_meshgrid(*quad_nodes)
            self._Cur_xietasigma_ = xietasigma
            self._Cur_quad_weights_ = quad_weights
            self._Cur_bfSelf_ = bfSelf
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
                self._Orth_xietasigma_,
                self._Orth_quad_weights_,
                self._Orth_bfSelf_,
                self._Orth_bfSelf_,
            )

        else:
            Mi = self._sf_.___PRIVATE_operator_inner___(
                self._sf_,
                i,
                self._Cur_xietasigma_,
                self._Cur_quad_weights_,
                self._Cur_bfSelf_,
                self._Cur_bfSelf_,
            )
        return Mi