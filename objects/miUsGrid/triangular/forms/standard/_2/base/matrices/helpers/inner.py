# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/29 8:51 PM
"""
from screws.freeze.base import FrozenOnly
import numpy as np
from scipy import sparse as spspa

class ___MassMatrix_Inner___(FrozenOnly):
    """The class for the inner product matrix."""
    def __init__(self, sf, xi_eta, quad_weights, bf):
        self._sf_ = sf
        self._xi_eta_ = xi_eta
        self._quad_weights_ = quad_weights
        self._bf_ = bf
        self._freeze_self_()

    def __call__(self, i):
        Mi = self.___Pr_operator_inner___(
            i, self._xi_eta_, self._quad_weights_, self._bf_
        )
        return Mi

    def ___Pr_operator_inner___(self, i, xietasigma, quad_weights, bf):
        """Note that here we only return a local matrix."""
        element = self._sf_.mesh.elements[i]
        detJ = element.coordinate_transformation.Jacobian(*xietasigma)

        Mi = np.einsum('im, jm, m -> ij',
            bf[0], bf[0], np.reciprocal(detJ)*quad_weights,
            optimize='greedy'
        )
        Mi = spspa.csr_matrix(Mi)
        return Mi