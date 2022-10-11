# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/29 8:51 PM
"""
from screws.freeze.base import FrozenOnly
import numpy as np
from scipy import sparse as spspa



class ___Operators_Inner___(FrozenOnly):
    """The class for the inner product matrix."""
    def __init__(self, sf, xi_eta, quad_weights, bf, of):
        self._sf_ = sf
        self._xi_eta_ = xi_eta
        self._quad_weights_ = quad_weights
        self._bf_ = bf
        self._of_ = of
        self._freeze_self_()

    def __call__(self, i):
        Mi = self.___Pr_operator_inner___(
            i, self._xi_eta_, self._quad_weights_, self._bf_, self._of_
        )
        return Mi

    def ___Pr_operator_inner___(self, i, xietasigma, quad_weights, bf, of):
        """Note that here we only return a local matrix."""
        element = self._sf_.mesh.elements[i]
        detJ = element.coordinate_transformation.Jacobian(*xietasigma)

        Mi = np.einsum('im, jm, m -> ij',
            bf[0], of[0], detJ*quad_weights,
            optimize='greedy'
        )
        Mi = spspa.csr_matrix(Mi)
        return Mi