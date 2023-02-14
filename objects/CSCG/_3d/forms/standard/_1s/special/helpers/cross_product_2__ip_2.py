# -*- coding: utf-8 -*-

import numpy as np
from scipy.sparse import csc_matrix
from components.freeze.main import FrozenOnly
from tools.elementwiseCache.dataStructures.objects.multiDimMatrix.main import MultiDimMatrix


class ___3dCSCG_1Form_CrossProduct_2__ip_2___(FrozenOnly):
    """
    The class for the inner wedge matrix; representing :math:`(\\omega \\times u, e)`.

    To do:

        (w X u, e) where w is 1 form u, e are 2-forms.

    We will have to call function from w, the 1-form, for example,

        CP = w.special.cross_product(u, e)

    This will give an instance `CP `of class `___3dCSCG_1Form_CrossProduct_2__ip_2___`. And when we
    do
        CP[i]

    it will give a matrix for mesh-element #i, and the columns represent the local
    dofs of e and the rows represent the local dofs of u.

    """
    def __init__(self, w1, u2, e2, quad_degree=None):
        """
        (w1 X u2, e2).

        :param w1: 1-standard-form
        :param u2: 2-standard-form
        :param e2: 2-standard-form
        :param quad_degree:
        """

        assert u2.ndim == w1.ndim == e2.ndim, " <___3dCSCG_1Form_CrossProduct_2__ip_2___> "
        assert u2.k == e2.k == 2, " <___3dCSCG_1Form_CrossProduct_2__ip_2___> "
        assert u2.mesh == w1.mesh, "___3dCSCG_1Form_CrossProduct_2__ip_2___: Meshes do not match."
        assert u2.mesh == e2.mesh, "___3dCSCG_1Form_CrossProduct_2__ip_2___: Meshes do not match."

        if quad_degree is None:
            quad_degree = [int(np.max([u2.dqp[i], w1.dqp[i], e2.dqp[i]]) * 1.5) for i in range(3)]

        quad_nodes, _, quad_weights_1d = \
            u2.space.___PRIVATE_do_evaluate_quadrature___(quad_degree)

        RMw = w1.do.make_reconstruction_matrix_on_grid(*quad_nodes)
        RMu = u2.do.make_reconstruction_matrix_on_grid(*quad_nodes)
        RMe = e2.do.make_reconstruction_matrix_on_grid(*quad_nodes)

        xi, et, sg = np.meshgrid(*quad_nodes, indexing='ij')
        xi = xi.ravel()
        et = et.ravel()
        sg = sg.ravel()
        detJ = w1.mesh.elements.coordinate_transformation.Jacobian(xi, et, sg)

        CP_IP_3dM = dict()
        type_cache = dict()
        for i in RMw:  # go through all local mesh-elements
            typeWr2Metric = w1.mesh.elements[i].type_wrt_metric.mark
            if isinstance(typeWr2Metric, str):
                if typeWr2Metric in type_cache:
                    CP_IP_3dM[i] = type_cache[typeWr2Metric]
                else:
                    wx, wy, wz = RMw[i]
                    u, v, w = RMu[i]
                    a, b, c = RMe[i]
                    dJi = detJ[i]
                    # w1 = [wx wy, wz]^T    u2= [u v w]^T   e2= [a b c]^T
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
                # w1 = [wx wy, wz]^T    u2= [u v w]^T   e2= [a b c]^T
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
        self._w1_ = w1
        self._u2_ = u2
        self._e2_ = e2
        self._freeze_self_()

    def __call__(self, i):
        """return 2d matrix of output = '2-M-1' type for mesh-element #i."""
        M = np.einsum('ijk, i -> kj', self._CP_IP_3dM_[i], self._w1_.cochain.local[i], optimize='greedy')
        return csc_matrix(M)


    @property
    def MDM(self):
        """Return a multi-dimension matrix representing this triple-operator."""
        return MultiDimMatrix(self._w1_.mesh.elements, self._CP_IP_3dM_, [self._w1_, self._u2_, self._e2_], 'no_cache')
