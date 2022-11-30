# -*- coding: utf-8 -*-

from components.freeze.base import FrozenOnly
import numpy as np
from scipy.sparse import csr_matrix

from tools.elementwiseCache.dataStructures.objects.multiDimMatrix.main import MultiDimMatrix


class ___0_x_1__ip__1___(FrozenOnly):
    """To compute (w0 X u1, e1). Cochain of w0 must be known, return mesh-element-wise
    matrices whose columns refer to local cochain of u1, rows' refer to local cochain
    of e1.
    """
    def __init__(self, w0, u1, e1, quad_degree=None):
        """
        We do not need to care about the orientations of the forms in this
        class.

        Compute (w0 X u1, e1) over the whole domain.

        We will call this from w0, the 0-form. u1 and e1 are two 1-form,
        we should do:

            CP = w0.special.cross_product(u1, e1)

        This will give an instance `CP `of class `___2dCSCG_0Form_CrossProduct_0_X_1__inner__1___`.
        And when we do

            CP[i]

        it will give a matrix for mesh-element #i, and the rows represent the local
        dofs of e1 and the cols represent the local dofs of u1.

        :param w0: will be the 0-form itself. We call this function form it.
        :param u1:
        :param e1:
        :param quad_degree:
        """

        assert u1.ndim == w0.ndim == e1.ndim, "___0_x_1__ip__1___"
        assert u1.k == e1.k == 1, "___0_x_1__ip__1___ "
        assert u1.mesh == w0.mesh, "___0_x_1__ip__1___: Meshes do not match."
        assert u1.mesh == e1.mesh, "___0_x_1__ip__1___: Meshes do not match."

        if quad_degree is None:
            _ = int(max([u1.p, w0.p, e1.p]) * 2)
            quad_degree = [_, _]

        quad_nodes, _, quad_weights_1d = u1.space.evaluation.quadrature(quad_degree)

        self._qn_ = quad_nodes
        self._qw_ = quad_weights_1d

        self._w0_ = w0
        self._u1_ = u1
        self._e1_ = e1

        self._3dM_cache_ = dict()
        self._freeze_self_()

    @property
    def MDM(self):
        """Return a multi-dimension matrix representing this triple-operator."""
        return MultiDimMatrix(self._w0_.mesh.elements,
                              self.___Pr_realtime_3dM___,
                              [self._w0_, self._u1_, self._e1_], 'no_cache')

    def ___Pr_realtime_3dM___(self, e):
        """Compute the 3d data for element #e."""

        cKey = self._u1_.___Pr_EWC_cache_key___

        if cKey != 'no_cache':
            cKey = cKey(e)
        else:
            pass

        if cKey != 'no_cache' and cKey in self._3dM_cache_:
            return self._3dM_cache_[cKey]

        RMw = self._w0_.do.make_reconstruction_matrix_on_grid(*self._qn_, element_range=e)[e]
        RMu = self._u1_.do.make_reconstruction_matrix_on_grid(*self._qn_, element_range=e)[e]
        RMe = self._e1_.do.make_reconstruction_matrix_on_grid(*self._qn_, element_range=e)[e]

        xi, et = np.meshgrid(*self._qn_, indexing='ij')
        xi = xi.ravel('F')
        et = et.ravel('F')
        detJ = self._w0_.mesh.elements[e].coordinate_transformation.Jacobian(xi, et)
        w = RMw
        u, v = RMu
        a, b = RMe
        # so, w0 = [0 0 w]^T, u1 = [u, v, 0]^T, e1 = [a b 0]^T, A = w0 X u1 = [-wv wu 0]^T, (A, e1) = -wva + wub
        CP_IP_3dM_i_ = - np.einsum('li, lj, lk, l -> ijk', w, v, a, self._qw_ * detJ, optimize='greedy')\
                       + np.einsum('li, lj, lk, l -> ijk', w, u, b, self._qw_ * detJ, optimize='greedy')

        if cKey != 'no_cache':
            self._3dM_cache_[cKey] = CP_IP_3dM_i_

        return CP_IP_3dM_i_

    def ___2_M_0___(self, e):
        """"""
        M = np.einsum('ijk, j -> ki', self.___Pr_realtime_3dM___(e), self._u1_.cochain.local[e], optimize='greedy')
        return csr_matrix(M)

    def ___2_M_1___(self, e):
        """"""
        M = np.einsum('ijk, i -> kj', self.___Pr_realtime_3dM___(e), self._w0_.cochain.local[e], optimize='greedy')
        return csr_matrix(M)

