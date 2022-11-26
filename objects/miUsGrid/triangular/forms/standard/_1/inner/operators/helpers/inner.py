# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/29 8:51 PM
"""
from components.freeze.base import FrozenOnly
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
        J = element.coordinate_transformation.Jacobian_matrix(*xietasigma)
        sqrt_g = element.coordinate_transformation.Jacobian(*xietasigma, itmD=J)
        iJ = element.coordinate_transformation.inverse_Jacobian_matrix(*xietasigma, itmD=J)
        g = element.coordinate_transformation.inverse_metric_matrix(*xietasigma, itmD=iJ)

        M00 = self.___PRIVATE_inner_Helper1___(quad_weights, sqrt_g*g[0][0], bf[i][0], of[i][0])
        M11 = self.___PRIVATE_inner_Helper1___(quad_weights, sqrt_g*g[1][1], bf[i][1], of[i][1])
        M01 = self.___PRIVATE_inner_Helper1___(quad_weights, sqrt_g*g[0][1], bf[i][0], of[i][1])
        M10 = self.___PRIVATE_inner_Helper1___(quad_weights, sqrt_g*g[1][0], bf[i][1], of[i][0])
        Mi = spspa.bmat([(M00, M01),
                         (M10, M11)], format='csr')
        return Mi


    @staticmethod
    def ___PRIVATE_inner_Helper1___(quad_weights, sqrt_g_g, bfS, bfO):
        M = np.einsum('m, im, jm -> ij', quad_weights * sqrt_g_g, bfS, bfO, optimize='greedy')
        return spspa.csc_matrix(M)