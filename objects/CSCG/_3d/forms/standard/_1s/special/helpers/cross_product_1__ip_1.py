# -*- coding: utf-8 -*-
import numpy as np
from scipy.sparse import csr_matrix
from components.freeze.main import FrozenOnly


class ___3dCSCG_1Form_CrossProduct_1__ip_1___(FrozenOnly):
    """
    The class for the inner wedge matrix; representing :math:`(\\omega \\times u, e)`.

    To do:

        (w X u, e) where w, u, e all are 1-forms.

    We will have to call function from w, the 1-form. For example

        CP = w.special.cross_product(u, e)

    This will give an instance `CP `of class `___3dCSCG_1Form_CrossProduct___`. And when we
    do
        CP[i]

    it will give a matrix for mesh-element #i, and the columns represent the local
    dofs of e and the rows represent the local dofs of u.

    :param w1:
    :param u1:
    :param e1:
    :param quad_degree:
    """
    def __init__(self, w1, u1, e1, quad_degree=None):
        """
        (w1 X u1, e1).

        :param w1: 1-standard-form
        :param u1: 1-standard-form
        :param e1: 1-standard-form
        :param quad_degree:
        """

        assert u1.ndim == w1.ndim == e1.ndim, " <___3dCSCG_1Form_CrossProduct_1__ip_1___> "
        assert u1.k == e1.k == 1, " <___3dCSCG_1Form_CrossProduct_1__ip_1___> "
        assert u1.mesh == w1.mesh, "___3dCSCG_1Form_CrossProduct_1__ip_1___: Meshes do not match."
        assert u1.mesh == e1.mesh, "___3dCSCG_1Form_CrossProduct_1__ip_1___: Meshes do not match."

        if quad_degree is None:
            quad_degree = [int(np.max([u1.dqp[i], w1.dqp[i], e1.dqp[i]]) * 1.5) for i in range(3)]

        quad_nodes, _, quad_weights_1d = \
            u1.space.___PRIVATE_do_evaluate_quadrature___(quad_degree)

        RMw = w1.do.make_reconstruction_matrix_on_grid(*quad_nodes)
        RMu = u1.do.make_reconstruction_matrix_on_grid(*quad_nodes)
        RMe = e1.do.make_reconstruction_matrix_on_grid(*quad_nodes)

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
        self._u1_ = u1
        self._e1_ = e1
        self._freeze_self_()

    def __call__(self, i):
        """return 2d matrix of output = '2-M-1' type for mesh-element #i."""
        M = np.einsum('ijk, i -> kj', self._CP_IP_3dM_[i], self._w1_.cochain.local[i], optimize='greedy')
        return csr_matrix(M)

    def _2_M_0_(self, i):
        """return 2d matrix of output = '2-M-1' type for mesh-element #i."""
        M = np.einsum('ijk, i -> ki', self._CP_IP_3dM_[i], self._u1_.cochain.local[i], optimize='greedy')
        return csr_matrix(M)

    #     """
    #     (w1 X u1, e1).
    #     """
    #     assert u1.ndim == w1.ndim == e1.ndim, " <_3dCSCG_1Form_CrossProduct> "
    #     assert u1.k == w1.k == e1.k == 1, " <_3dCSCG_1Form_CrossProduct> "
    #     assert u1.mesh == w1.mesh, "Meshes do not match."
    #     assert u1.mesh == e1.mesh, "Meshes do not match."
    #     self._mesh_ = u1.mesh
    #     self._u_ = u1
    #     self._w_ = w1
    #     self._e_ = e1
    #
    #     if quad_degree is None:
    #         quad_degree = [int(np.max([u1.dqp[i], w1.dqp[i], e1.dqp[i]])) for i in range(3)]
    #     quad_nodes, _, quad_weights = u1.space.___PRIVATE_do_evaluate_quadrature___(quad_degree)
    #     xietasigma, wbf = w1.do.evaluate_basis_at_meshgrid(*quad_nodes, compute_xietasigma=True)
    #     _, ubf = u1.do.evaluate_basis_at_meshgrid(*quad_nodes, compute_xietasigma=False)
    #     if e1 is u1:
    #         ebf = ubf
    #     else:
    #         _, ebf = e1.do.evaluate_basis_at_meshgrid(*quad_nodes, compute_xietasigma=False)
    #     self._xietasigma_ = xietasigma
    #     self._qw_ = quad_weights
    #     self._wbf_ = wbf
    #     self._ubf_ = ubf
    #     self._ebf_ = ebf
    #     self._JM_ = self._mesh_.elements.coordinate_transformation.Jacobian_matrix(*xietasigma)
    #     self._iJ_ = self._mesh_.elements.coordinate_transformation.inverse_Jacobian_matrix(*xietasigma, J=self._JM_)
    #     self.___PRIVATE_reset_cache___()
    #     self._freeze_self_()
    #
    # def ___PRIVATE_reset_cache___(self):
    #     self._J_cache_ = dict()
    #
    # def _J_(self, i):
    #     element = self._mesh_.elements[i]
    #     typeWr2Metric = element.type_wrt_metric.mark
    #     if typeWr2Metric in self._J_cache_:
    #         return self._J_cache_[typeWr2Metric]
    #     else:
    #         JM = self._JM_[i]
    #         if isinstance(typeWr2Metric, str) and  typeWr2Metric[:4] == 'Orth':
    #             J00 = JM[1][1] * JM[2][2]
    #             J01 = 0
    #             J02 = 0
    #             J10 = 0
    #             J11 = JM[2][2] * JM[0][0]
    #             J12 = 0
    #             J20 = 0
    #             J21 = 0
    #             J22 = JM[0][0] * JM[1][1]
    #         else:
    #             J00 = JM[1][1] * JM[2][2] - JM[1][2] * JM[2][1]
    #             J01 = JM[2][1] * JM[0][2] - JM[2][2] * JM[0][1]
    #             J02 = JM[0][1] * JM[1][2] - JM[0][2] * JM[1][1]
    #             J10 = JM[1][2] * JM[2][0] - JM[1][0] * JM[2][2]
    #             J11 = JM[2][2] * JM[0][0] - JM[2][0] * JM[0][2]
    #             J12 = JM[0][2] * JM[1][0] - JM[0][0] * JM[1][2]
    #             J20 = JM[1][0] * JM[2][1] - JM[1][1] * JM[2][0]
    #             J21 = JM[2][0] * JM[0][1] - JM[2][1] * JM[0][0]
    #             J22 = JM[0][0] * JM[1][1] - JM[0][1] * JM[1][0]
    #         J = (J00, J01, J02, J10, J11, J12, J20, J21, J22)
    #         # cache it even for unique mesh cells (because we may use them multiple times when do temporal iterations.)
    #         self._J_cache_[typeWr2Metric] = J
    #         return J
    #
    # def __call__(self, i):
    #     typeWr2Metric = self._mesh_.elements[i].type_wrt_metric.mark
    #
    #     u0, u1, u2 =  self._ubf_ # a
    #     w0, w1, w2 =  self._wbf_ # b; given
    #     e0, e1, e2 =  self._ebf_ # epsilon
    #
    #     b0p = np.einsum('ij, i -> j', w0, self._w_.cochain.___PRIVATE_local_on_axis___('x', i), optimize='greedy')
    #     b1p = np.einsum('ij, i -> j', w1, self._w_.cochain.___PRIVATE_local_on_axis___('y', i), optimize='greedy')
    #     b2p = np.einsum('ij, i -> j', w2, self._w_.cochain.___PRIVATE_local_on_axis___('z', i), optimize='greedy')
    #
    #     iJ = self._iJ_[i]
    #     J00, J01, J02, J10, J11, J12, J20, J21, J22 = self._J_(i)
    #
    #     if isinstance(typeWr2Metric, str) and  typeWr2Metric[:4] == 'Orth':
    #         b0 = b0p * iJ[0][0]
    #         b1 = b1p * iJ[1][1]
    #         b2 = b2p * iJ[2][2]
    #
    #         B01, B02 = -iJ[0][0]*b2, iJ[0][0]*b1
    #         B10, B12 = iJ[1][1]*b2, -iJ[1][1]*b0
    #         B20, B21 = -iJ[2][2]*b1, iJ[2][2]*b0
    #
    #         m01 = B10*J00
    #         m02 = B20*J00
    #         m10 = B01*J11
    #         m12 = B21*J11
    #         m20 = B02*J22
    #         m21 = B12*J22
    #
    #         # put `-` because a x b = - b x a
    #         M01 = - np.einsum('iw, jw, w -> ij', e0, u1, m01*self._qw_, optimize='greedy')
    #         M02 = - np.einsum('iw, jw, w -> ij', e0, u2, m02*self._qw_, optimize='greedy')
    #
    #         M10 = - np.einsum('iw, jw, w -> ij', e1, u0, m10*self._qw_, optimize='greedy')
    #         M12 = - np.einsum('iw, jw, w -> ij', e1, u2, m12*self._qw_, optimize='greedy')
    #
    #         M20 = - np.einsum('iw, jw, w -> ij', e2, u0, m20*self._qw_, optimize='greedy')
    #         M21 = - np.einsum('iw, jw, w -> ij', e2, u1, m21*self._qw_, optimize='greedy')
    #
    #         M = ([None, spspa.csc_matrix(M01), spspa.csc_matrix(M02)],
    #              [spspa.csc_matrix(M10), None, spspa.csc_matrix(M12)],
    #              [spspa.csc_matrix(M20), spspa.csc_matrix(M21), None])
    #
    #     else:
    #         b0 = b0p * iJ[0][0] + b1p * iJ[1][0] + b2p * iJ[2][0]
    #         b1 = b0p * iJ[0][1] + b1p * iJ[1][1] + b2p * iJ[2][1]
    #         b2 = b0p * iJ[0][2] + b1p * iJ[1][2] + b2p * iJ[2][2]
    #
    #         B00, B01, B02 = (iJ[0][1]*b2-iJ[0][2]*b1), (iJ[0][2]*b0-iJ[0][0]*b2), (iJ[0][0]*b1-iJ[0][1]*b0)
    #         B10, B11, B12 = (iJ[1][1]*b2-iJ[1][2]*b1), (iJ[1][2]*b0-iJ[1][0]*b2), (iJ[1][0]*b1-iJ[1][1]*b0)
    #         B20, B21, B22 = (iJ[2][1]*b2-iJ[2][2]*b1), (iJ[2][2]*b0-iJ[2][0]*b2), (iJ[2][0]*b1-iJ[2][1]*b0)
    #
    #         m00 = B00*J00 + B01*J01 + B02*J02
    #         m01 = B10*J00 + B11*J01 + B12*J02
    #         m02 = B20*J00 + B21*J01 + B22*J02
    #         m10 = B00*J10 + B01*J11 + B02*J12
    #         m11 = B10*J10 + B11*J11 + B12*J12
    #         m12 = B20*J10 + B21*J11 + B22*J12
    #         m20 = B00*J20 + B01*J21 + B02*J22
    #         m21 = B10*J20 + B11*J21 + B12*J22
    #         m22 = B20*J20 + B21*J21 + B22*J22
    #
    #         # put `-` because a x b = - b x a
    #         M00 = - np.einsum('iw, jw, w -> ij', e0, u0, m00*self._qw_, optimize='greedy')
    #         M01 = - np.einsum('iw, jw, w -> ij', e0, u1, m01*self._qw_, optimize='greedy')
    #         M02 = - np.einsum('iw, jw, w -> ij', e0, u2, m02*self._qw_, optimize='greedy')
    #
    #         M10 = - np.einsum('iw, jw, w -> ij', e1, u0, m10*self._qw_, optimize='greedy')
    #         M11 = - np.einsum('iw, jw, w -> ij', e1, u1, m11*self._qw_, optimize='greedy')
    #         M12 = - np.einsum('iw, jw, w -> ij', e1, u2, m12*self._qw_, optimize='greedy')
    #
    #         M20 = - np.einsum('iw, jw, w -> ij', e2, u0, m20*self._qw_, optimize='greedy')
    #         M21 = - np.einsum('iw, jw, w -> ij', e2, u1, m21*self._qw_, optimize='greedy')
    #         M22 = - np.einsum('iw, jw, w -> ij', e2, u2, m22*self._qw_, optimize='greedy')
    #
    #         M = ([spspa.csc_matrix(M00), spspa.csc_matrix(M01), spspa.csc_matrix(M02)],
    #              [spspa.csc_matrix(M10), spspa.csc_matrix(M11), spspa.csc_matrix(M12)],
    #              [spspa.csc_matrix(M20), spspa.csc_matrix(M21), spspa.csc_matrix(M22)])
    #     MW = spspa.bmat(M, format='csc')
    #
    #     return MW