# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/05 9:20 PM
"""
import numpy as np
from scipy.sparse import csr_matrix
from screws.freeze.base import FrozenOnly


class ___3dCSCG_2Form_CrossProduct_2__ip_2_2M0___(FrozenOnly):
    """
    The class for the inner wedge matrix; representing :math:`(a \\times b, e)`.

    a, b, c are all 2-forms.

    """
    def __init__(self, w2, u2, e2, quad_degree=None, cache=None):
        """
        (w2 X u2, e2). The cochain of u2 must be known.

        :param w2: 2-standard-form
        :param u2: 2-standard-form
        :param e2: 2-standard-form
        :param quad_degree:
        :param cache:
            A dict to pass cache data to me.
        """

        assert u2.ndim == w2.ndim == e2.ndim, " <___3dCSCG_2Form_CrossProduct_2__ip_2___> "
        assert u2.k == e2.k == 2, " <___3dCSCG_2Form_CrossProduct_2__ip_2___> "
        assert u2.mesh == w2.mesh, "___3dCSCG_2Form_CrossProduct_2__ip_2___: Meshes do not match."
        assert u2.mesh == e2.mesh, "___3dCSCG_2Form_CrossProduct_2__ip_2___: Meshes do not match."

        if cache is not None and '3D_matrix' in cache:
            self._CP_IP_3dM_ = cache['3D_matrix']

        else:
            if quad_degree is None:
                quad_degree = [int(np.max([u2.dqp[i], w2.dqp[i], e2.dqp[i]]) * 1.5)  for i in range(3)]

            quad_nodes, _, quad_weights_1d = \
                u2.space.___PRIVATE_do_evaluate_quadrature___(quad_degree)

            RMw = w2.do.make_reconstruction_matrix_on_grid(*quad_nodes)
            RMu = u2.do.make_reconstruction_matrix_on_grid(*quad_nodes)
            RMe = e2.do.make_reconstruction_matrix_on_grid(*quad_nodes)

            xi, et, sg = np.meshgrid(*quad_nodes, indexing='ij')
            xi = xi.ravel()
            et = et.ravel()
            sg = sg.ravel()
            detJ = w2.mesh.elements.coordinate_transformation.Jacobian(xi, et, sg)

            CP_IP_3dM = dict()
            type_cache = dict()
            for i in RMw: # go through all local mesh-elements
                typeWr2Metric = w2.mesh.elements[i].type_wrt_metric.mark
                if isinstance(typeWr2Metric, str):
                    if typeWr2Metric in type_cache:
                        CP_IP_3dM[i] = type_cache[typeWr2Metric]
                    else:
                        wx, wy, wz = RMw[i]
                        u, v, w = RMu[i]
                        a, b, c = RMe[i]
                        dJi = detJ[i]
                        # w2 = [wx wy, wz]^T    u2= [u v w]^T   e2= [a b c]^T
                        # WXU = w1 X u2 = [wy*w - wz*v,   wz*u - wx*w,   wx*v - wy*u]^T = [A B C]^T
                        # WXU dot e2 = Aa + Bb + Cc
                        Aa = + np.einsum('li, lj, lk, l -> ijk', wy, w, a, quad_weights_1d * dJi, optimize='greedy')\
                             - np.einsum('li, lj, lk, l -> ijk', wz, v, a, quad_weights_1d * dJi, optimize='greedy')
                        Bb = + np.einsum('li, lj, lk, l -> ijk', wz, u, b, quad_weights_1d * dJi, optimize='greedy')\
                             - np.einsum('li, lj, lk, l -> ijk', wx, w, b, quad_weights_1d * dJi, optimize='greedy')
                        Cc = + np.einsum('li, lj, lk, l -> ijk', wx, v, c, quad_weights_1d * dJi, optimize='greedy')\
                             - np.einsum('li, lj, lk, l -> ijk', wy, u, c, quad_weights_1d * dJi, optimize='greedy')
                        CP_IP_3dM_i_ = Aa + Bb + Cc
                        CP_IP_3dM[i] = CP_IP_3dM_i_
                        type_cache[typeWr2Metric] = CP_IP_3dM_i_

                else:
                    wx, wy, wz = RMw[i]
                    u, v, w = RMu[i]
                    a, b, c = RMe[i]
                    dJi = detJ[i]
                    # w2 = [wx wy, wz]^T    u2= [u v w]^T   e2= [a b c]^T
                    # WXU = w1 X u2 = [wy*w - wz*v,   wz*u - wx*w,   wx*v - wy*u]^T = [A B C]^T
                    # WXU dot e2 = Aa + Bb + Cc
                    Aa = + np.einsum('li, lj, lk, l -> ijk', wy, w, a, quad_weights_1d * dJi, optimize='greedy')\
                         - np.einsum('li, lj, lk, l -> ijk', wz, v, a, quad_weights_1d * dJi, optimize='greedy')
                    Bb = + np.einsum('li, lj, lk, l -> ijk', wz, u, b, quad_weights_1d * dJi, optimize='greedy')\
                         - np.einsum('li, lj, lk, l -> ijk', wx, w, b, quad_weights_1d * dJi, optimize='greedy')
                    Cc = + np.einsum('li, lj, lk, l -> ijk', wx, v, c, quad_weights_1d * dJi, optimize='greedy')\
                         - np.einsum('li, lj, lk, l -> ijk', wy, u, c, quad_weights_1d * dJi, optimize='greedy')

                    CP_IP_3dM[i] = Aa + Bb + Cc

            self._CP_IP_3dM_ = CP_IP_3dM
        self._u2_ = u2
        self._freeze_self_()

    def __call__(self, e):
        """return 2d matrix of output = '2-M-0' type for mesh-element #e."""
        M = np.einsum('ijk, j -> ki', self._CP_IP_3dM_[e], self._u2_.cochain.local[e], optimize='greedy')
        return csr_matrix(M)