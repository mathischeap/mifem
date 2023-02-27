# -*- coding: utf-8 -*-
import numpy as np
from scipy.sparse import csr_matrix
from components.freeze.base import FrozenOnly
from tools.elementwiseCache.dataStructures.objects.multiDimMatrix.main import MultiDimMatrix


class ___3dCSCG_2Form_CrossProduct_1__ip_2___(FrozenOnly):
    """
    The class for the inner wedge matrix; representing :math:`(w \\times u, e)`.

    w, e are 2-forms. u is 1-form

    """
    def __init__(self, w2, u1, e2, quad_degree=None):
        """
        (w2 X u1, e2).

        :param w2: 2-standard-form
        :param u1: 1-standard-form
        :param e2: 2-standard-form
        :param quad_degree:
        """

        assert u1.ndim == w2.ndim == e2.ndim, " <___3dCSCG_2Form_CrossProduct_1__ip_2___> "
        assert u1.k + 1 == e2.k == 2, " <___3dCSCG_2Form_CrossProduct_1__ip_2___> "
        assert u1.mesh == w2.mesh, "___3dCSCG_2Form_CrossProduct_1__ip_2___: Meshes do not match."
        assert u1.mesh == e2.mesh, "___3dCSCG_2Form_CrossProduct_1__ip_2___: Meshes do not match."

        if quad_degree is None:
            quad_degree = [int(np.max([u1.dqp[i], w2.dqp[i], e2.dqp[i]]) * 1.5) for i in range(3)]

        quad_nodes, _, quad_weights_1d = \
            u1.space.___PRIVATE_do_evaluate_quadrature___(quad_degree)

        RMw = w2.do.make_reconstruction_matrix_on_grid(*quad_nodes)
        RMu = u1.do.make_reconstruction_matrix_on_grid(*quad_nodes)
        RMe = e2.do.make_reconstruction_matrix_on_grid(*quad_nodes)

        xi, et, sg = np.meshgrid(*quad_nodes, indexing='ij')
        xi = xi.ravel()
        et = et.ravel()
        sg = sg.ravel()
        detJ = w2.mesh.elements.coordinate_transformation.Jacobian(xi, et, sg)

        CP_IP_3dM = dict()
        type_cache = dict()
        for i in RMw:  # go through all local mesh-elements
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
                        - np.einsum('li, lj, lk, l -> ijk',  wz, v, a, quad_weights_1d * dJi, optimize='greedy')
                    Bb = + np.einsum('li, lj, lk, l -> ijk', wz, u, b, quad_weights_1d * dJi, optimize='greedy')\
                        - np.einsum('li, lj, lk, l -> ijk',  wx, w, b, quad_weights_1d * dJi, optimize='greedy')
                    Cc = + np.einsum('li, lj, lk, l -> ijk', wx, v, c, quad_weights_1d * dJi, optimize='greedy')\
                        - np.einsum('li, lj, lk, l -> ijk',  wy, u, c, quad_weights_1d * dJi, optimize='greedy')
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
                    - np.einsum('li, lj, lk, l -> ijk',  wz, v, a, quad_weights_1d * dJi, optimize='greedy')
                Bb = + np.einsum('li, lj, lk, l -> ijk', wz, u, b, quad_weights_1d * dJi, optimize='greedy')\
                    - np.einsum('li, lj, lk, l -> ijk',  wx, w, b, quad_weights_1d * dJi, optimize='greedy')
                Cc = + np.einsum('li, lj, lk, l -> ijk', wx, v, c, quad_weights_1d * dJi, optimize='greedy')\
                    - np.einsum('li, lj, lk, l -> ijk',  wy, u, c, quad_weights_1d * dJi, optimize='greedy')

                CP_IP_3dM[i] = Aa + Bb + Cc

        self._CP_IP_3dM_ = CP_IP_3dM
        self._w2_ = w2
        self._u1_ = u1
        self._e2_ = e2
        self._freeze_self_()

    def __call__(self, i):
        """return 2d matrix of output = '2-M-1' type for mesh-element #i."""
        M = np.einsum('ijk, i -> kj', self._CP_IP_3dM_[i], self._w2_.cochain.local[i], optimize='greedy')
        return csr_matrix(M)

    @property
    def MDM(self):
        """Return a multi-dimension matrix representing this triple-operator."""
        return MultiDimMatrix(self._w2_.mesh.elements, self._CP_IP_3dM_, [self._w2_, self._u1_, self._e2_], 'no_cache')

    def _2_M_0_(self, i):
        """return 2d matrix of output = '2-M-0' type for mesh-element #i."""
        M = np.einsum('ijk, j -> ki', self._CP_IP_3dM_[i], self._u1_.cochain.local[i], optimize='greedy')
        return csr_matrix(M)
