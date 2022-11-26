# -*- coding: utf-8 -*-
import numpy as np
from scipy.sparse import csr_matrix
from components.freeze.base import FrozenOnly
from tools.linearAlgebra.elementwiseCache.objects.multiDimMatrix.main import MultiDimMatrix


class ___3dCSCG_2Form_CrossProduct_2__ip_2___(FrozenOnly):
    """
    The class for the inner wedge matrix; representing :math:`(a \\times b, e)`.

    a, b, c are all 2-forms.

    """
    def __init__(self, w2, u2, e2, quad_degree=None):
        """
        (w2 X u2, e2).

        :param w2: 2-standard-form
        :param u2: 2-standard-form
        :param e2: 2-standard-form
        :param quad_degree:
        """

        assert u2.ndim == w2.ndim == e2.ndim, " <___3dCSCG_2Form_CrossProduct_2__ip_2___> "
        assert u2.k == e2.k == 2, " <___3dCSCG_2Form_CrossProduct_2__ip_2___> "
        assert u2.mesh == w2.mesh, "___3dCSCG_2Form_CrossProduct_2__ip_2___: Meshes do not match."
        assert u2.mesh == e2.mesh, "___3dCSCG_2Form_CrossProduct_2__ip_2___: Meshes do not match."

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
        self._w2_ = w2
        self._u2_ = u2
        self._e2_ = e2
        self._freeze_self_()

    def __call__(self, i):
        """return 2d matrix of output = '2-M-1' type for mesh-element #i."""
        M = np.einsum('ijk, i -> kj', self._CP_IP_3dM_[i], self._w2_.cochain.local[i], optimize='greedy')
        return csr_matrix(M)

    @property
    def MDM(self):
        """Return a multi-dimension matrix representing this triple-operator."""
        return MultiDimMatrix(self._w2_.mesh.elements, self._CP_IP_3dM_, [self._w2_, self._u2_, self._e2_], 'no_cache')

    #     assert a.ndim == b.ndim == e.ndim, " <___3dCSCG_2Form_CrossProduct_2__ip_2___> "
    #     assert a.k == b.k == e.k == 2, " <___3dCSCG_2Form_CrossProduct_2__ip_2___> "
    #     assert a.mesh == b.mesh, "___3dCSCG_2Form_CrossProduct_2__ip_2___, Meshes do not match."
    #     assert a.mesh == e.mesh, "___3dCSCG_2Form_CrossProduct_2__ip_2___, Meshes do not match."
    #     self._mesh_ = a.mesh
    #     self._a_ = a
    #     self._b_ = b
    #     self._e_ = e
    #
    #     if quad_degree is None:
    #         quad_degree = [int(np.max([a.dqp[i], b.dqp[i], e.dqp[i]])) for i in range(3)]
    #     quad_nodes, _, quad_weights = a.space.___PRIVATE_do_evaluate_quadrature___(quad_degree)
    #     xietasigma, abf = a.do.evaluate_basis_at_meshgrid(*quad_nodes, compute_xietasigma=True)
    #     _         , bbf = b.do.evaluate_basis_at_meshgrid(*quad_nodes, compute_xietasigma=False)
    #     _         , ebf = e.do.evaluate_basis_at_meshgrid(*quad_nodes, compute_xietasigma=False)
    #     self._qw_ = quad_weights
    #     self._abf_ = abf
    #     self._bbf_ = bbf
    #     self._ebf_ = ebf
    #
    #     self._JM_    = self._mesh_.elements.coordinate_transformation.Jacobian_matrix(*xietasigma)
    #     self._sqrtg_ = self._mesh_.elements.coordinate_transformation.Jacobian(*xietasigma, J=self._JM_)
    #     self._iJ_    = self._mesh_.elements.coordinate_transformation.inverse_Jacobian_matrix(*xietasigma, J=self._JM_)
    #     self._g_     = self._mesh_.elements.coordinate_transformation.inverse_metric_matrix(*xietasigma, iJ=self._iJ_)
    #
    #     self.DO_reset_cache()
    #     self._freeze_self_()
    #
    # def DO_reset_cache(self):
    #     self._J_cache_ = dict()
    #     self._G_cache_ = dict()
    #
    # def _J_(self, i):
    #     element = self._mesh_.elements[i]
    #     typeWr2Metric = element.type_wrt_metric.mark
    #     if typeWr2Metric in self._J_cache_:
    #         return self._J_cache_[typeWr2Metric]
    #     else:
    #         iJ = self._iJ_[i]
    #         if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
    #             J00 = iJ[1][1] * iJ[2][2]
    #             J01 = 0
    #             J02 = 0
    #             J10 = 0
    #             J11 = iJ[2][2] * iJ[0][0]
    #             J12 = 0
    #             J20 = 0
    #             J21 = 0
    #             J22 = iJ[0][0] * iJ[1][1]
    #         else:
    #             J00 = iJ[1][1] * iJ[2][2] - iJ[1][2] * iJ[2][1]
    #             J01 = iJ[2][1] * iJ[0][2] - iJ[2][2] * iJ[0][1]
    #             J02 = iJ[0][1] * iJ[1][2] - iJ[0][2] * iJ[1][1]
    #             J10 = iJ[1][2] * iJ[2][0] - iJ[1][0] * iJ[2][2]
    #             J11 = iJ[2][2] * iJ[0][0] - iJ[2][0] * iJ[0][2]
    #             J12 = iJ[0][2] * iJ[1][0] - iJ[0][0] * iJ[1][2]
    #             J20 = iJ[1][0] * iJ[2][1] - iJ[1][1] * iJ[2][0]
    #             J21 = iJ[2][0] * iJ[0][1] - iJ[2][1] * iJ[0][0]
    #             J22 = iJ[0][0] * iJ[1][1] - iJ[0][1] * iJ[1][0]
    #         J = (J00, J01, J02, J10, J11, J12, J20, J21, J22)
    #         # cache it even for unique mesh cells (because we may use them multiple times when do temporal iterations.)
    #         self._J_cache_[typeWr2Metric] = J
    #         return J
    #
    # def _G_(self, i):
    #     element = self._mesh_.elements[i]
    #     typeWr2Metric = element.type_wrt_metric.mark
    #     if typeWr2Metric in self._G_cache_:
    #         return self._G_cache_[typeWr2Metric]
    #     else:
    #         sqrtg = self._sqrtg_[i]
    #         g = self._g_[i]
    #         JM = self._JM_[i]
    #         if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
    #             G00 = sqrtg * g[1][1] * g[2][2] * JM[1][1] * JM[2][2]
    #             G01 = 0
    #             G02 = 0
    #             G10 = 0
    #             G11 = sqrtg * g[2][2] * g[0][0] * JM[2][2] * JM[0][0]
    #             G12 = 0
    #             G20 = 0
    #             G21 = 0
    #             G22 = sqrtg * g[0][0] * g[1][1] * JM[0][0] * JM[1][1]
    #         else:
    #             G00 = sqrtg * (g[1][1] * g[2][2] - g[1][2] * g[2][1])
    #             G01 = sqrtg * (g[1][2] * g[2][0] - g[1][0] * g[2][2])
    #             G02 = sqrtg * (g[1][0] * g[2][1] - g[1][1] * g[2][0])
    #             G10 = G01
    #             G11 = sqrtg * (g[2][2] * g[0][0] - g[2][0] * g[0][2])
    #             G12 = sqrtg * (g[2][0] * g[0][1] - g[2][1] * g[0][0])
    #             G20 = G02
    #             G21 = G12
    #             G22 = sqrtg * (g[0][0] * g[1][1] - g[0][1] * g[1][0])
    #         G = (G00, G01, G02, G10, G11, G12, G20, G21, G22)
    #         # cache it even for unique mesh cells (because we may use them multiple times when do temporal iterations.)
    #         self._G_cache_[typeWr2Metric] = G
    #         return G
    #
    # def __call__(self, i):
    #     typeWr2Metric = self._mesh_.elements[i].type_wrt_metric.mark
    #
    #     a0, a1, a2 =  self._abf_ # a; given
    #     b0, b1, b2 =  self._bbf_ # b
    #     e0, e1, e2 =  self._ebf_ # epsilon
    #
    #     a0p = np.einsum('ij, i -> j', a0, self._a_.cochain.___PRIVATE_local_on_axis___('x', i), optimize='greedy')
    #     a1p = np.einsum('ij, i -> j', a1, self._a_.cochain.___PRIVATE_local_on_axis___('y', i), optimize='greedy')
    #     a2p = np.einsum('ij, i -> j', a2, self._a_.cochain.___PRIVATE_local_on_axis___('z', i), optimize='greedy')
    #
    #     J00, J01, J02, J10, J11, J12, J20, J21, J22 = self._J_(i)
    #     G00, G01, G02, G10, G11, G12, G20, G21, G22 = self._G_(i)
    #
    #     if isinstance(typeWr2Metric, str) and typeWr2Metric[:4] == 'Orth':
    #         a0 = a0p * J00
    #         a1 = a1p * J11
    #         a2 = a2p * J22
    #
    #         A01, A02 =  J00*a2, -J00*a1
    #         A10, A12 = -J11*a2,  J11*a0
    #         A20, A21 =  J22*a1, -J22*a0
    #
    #         m01 = A10*G00
    #         m02 = A20*G00
    #         m10 = A01*G11
    #         m12 = A21*G11
    #         m20 = A02*G22
    #         m21 = A12*G22
    #
    #         M01 = np.einsum('iw, jw, w -> ij', e0, b1, m01*self._qw_, optimize='greedy')
    #         M02 = np.einsum('iw, jw, w -> ij', e0, b2, m02*self._qw_, optimize='greedy')
    #         M10 = np.einsum('iw, jw, w -> ij', e1, b0, m10*self._qw_, optimize='greedy')
    #         M12 = np.einsum('iw, jw, w -> ij', e1, b2, m12*self._qw_, optimize='greedy')
    #         M20 = np.einsum('iw, jw, w -> ij', e2, b0, m20*self._qw_, optimize='greedy')
    #         M21 = np.einsum('iw, jw, w -> ij', e2, b1, m21*self._qw_, optimize='greedy')
    #
    #         M = ([None                 , spspa.csc_matrix(M01), spspa.csc_matrix(M02)],
    #              [spspa.csc_matrix(M10), None                 , spspa.csc_matrix(M12)],
    #              [spspa.csc_matrix(M20), spspa.csc_matrix(M21), None                 ])
    #
    #     else:
    #         #TODO: correct the following ...
    #         raise NotImplementedError()
    #
    #         # a0 = a0p * J00 + a1p * J01 + a2p * J02
    #         # a1 = a0p * J10 + a1p * J11 + a2p * J12
    #         # a2 = a0p * J20 + a1p * J21 + a2p * J22
    #         #
    #         # A00, A01, A02 = (J20*a1 - J10*a2), (J00*a2 - J20*a0), (J10*a0 - J00*a1)
    #         # A10, A11, A12 = (J21*a1 - J11*a2), (J01*a2 - J21*a0), (J11*a0 - J01*a1)
    #         # A20, A21, A22 = (J22*a1 - J12*a2), (J02*a2 - J22*a0), (J12*a0 - J02*a1)
    #         #
    #         # m00 = A00*G00 + A01*G01 + A02*G02
    #         # m01 = A10*G00 + A11*G01 + A12*G02
    #         # m02 = A20*G00 + A21*G01 + A22*G02
    #         # m10 = A00*G10 + A01*G11 + A02*G12
    #         # m11 = A10*G10 + A11*G11 + A12*G12
    #         # m12 = A20*G10 + A21*G11 + A22*G12
    #         # m20 = A00*G20 + A01*G21 + A02*G22
    #         # m21 = A10*G20 + A11*G21 + A12*G22
    #         # m22 = A20*G20 + A21*G21 + A22*G22
    #         #
    #         # M00 = np.einsum('iw, jw, w -> ij', e0, b0, m00*self._qw_, optimize='greedy')
    #         # M01 = np.einsum('iw, jw, w -> ij', e0, b1, m01*self._qw_, optimize='greedy')
    #         # M02 = np.einsum('iw, jw, w -> ij', e0, b2, m02*self._qw_, optimize='greedy')
    #         # M10 = np.einsum('iw, jw, w -> ij', e1, b0, m10*self._qw_, optimize='greedy')
    #         # M11 = np.einsum('iw, jw, w -> ij', e1, b1, m11*self._qw_, optimize='greedy')
    #         # M12 = np.einsum('iw, jw, w -> ij', e1, b2, m12*self._qw_, optimize='greedy')
    #         # M20 = np.einsum('iw, jw, w -> ij', e2, b0, m20*self._qw_, optimize='greedy')
    #         # M21 = np.einsum('iw, jw, w -> ij', e2, b1, m21*self._qw_, optimize='greedy')
    #         # M22 = np.einsum('iw, jw, w -> ij', e2, b2, m22*self._qw_, optimize='greedy')
    #         # M = ([spspa.csc_matrix(M00), spspa.csc_matrix(M01), spspa.csc_matrix(M02)],
    #         #      [spspa.csc_matrix(M10), spspa.csc_matrix(M11), spspa.csc_matrix(M12)],
    #         #      [spspa.csc_matrix(M20), spspa.csc_matrix(M21), spspa.csc_matrix(M22)])
    #
    #     MW = spspa.bmat(M, format='csc')
    #     return MW
