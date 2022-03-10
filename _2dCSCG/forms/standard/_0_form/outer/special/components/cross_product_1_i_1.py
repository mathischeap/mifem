

from screws.freeze.inheriting.frozen_only import FrozenOnly
import numpy as np
from scipy.sparse import csc_matrix



class ___2dCSCG_0_o_Form_CrossProduct_0_X_1__ip_1___(FrozenOnly):
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

        it will give a matrix for mesh-element #i, and the columns represent the local
        dofs of e1 and the rows represent the local dofs of u1.

        :param w0: will be the 0-form itself. We call this function form it.
        :param u1:
        :param e1:
        :param quad_degree:
        """

        assert u1.ndim == w0.ndim == e1.ndim, " <___2dCSCG_0_o_Form_CrossProduct_0_X_1__ip_1___> "
        assert u1.k == e1.k == 1, " <___2dCSCG_0_o_Form_CrossProduct_0_X_1__ip_1___> "
        assert u1.mesh == w0.mesh, "___2dCSCG_0_o_Form_CrossProduct_0_X_1__ip_1___: Meshes do not match."
        assert u1.mesh == e1.mesh, "___2dCSCG_0_o_Form_CrossProduct_0_X_1__ip_1___: Meshes do not match."

        if quad_degree is None:
            quad_degree = [int(np.max([u1.dqp[i], w0.dqp[i], e1.dqp[i]])) * 2 for i in range(2)]

        quad_nodes, _, quad_weights_1d = \
            u1.space.___PRIVATE_do_evaluate_quadrature___(quad_degree)

        RMw = w0.do.make_reconstruction_matrix_on_grid(*quad_nodes)
        RMu = u1.do.make_reconstruction_matrix_on_grid(*quad_nodes)
        RMe = e1.do.make_reconstruction_matrix_on_grid(*quad_nodes)

        xi, et = np.meshgrid(*quad_nodes, indexing='ij')
        xi = xi.ravel()
        et = et.ravel()
        detJ = w0.mesh.elements.coordinate_transformation.Jacobian(xi, et)

        CP_IP_3dM = dict()
        type_cache = dict()
        for i in RMw: # go through all local mesh-elements
            typeWr2Metric = w0.mesh.elements[i].type_wrt_metric.mark
            if isinstance(typeWr2Metric, str):
                if typeWr2Metric in type_cache:
                    CP_IP_3dM[i] = type_cache[typeWr2Metric]
                else:
                    w = RMw[i]
                    u, v = RMu[i]
                    a, b = RMe[i]
                    dJi = detJ[i]
                    # so, w0 = [0 0 w]^T, u1 = [u, v, 0]^T, e1 = [a b 0]^T, A = w0 X u1 = [-wv wu 0]^T, (A, e1) = -wva + wub
                    CP_IP_3dM_i_ = - np.einsum('li, lj, lk, l -> ijk', w, v, a, quad_weights_1d * dJi, optimize='greedy')\
                                   + np.einsum('li, lj, lk, l -> ijk', w, u, b, quad_weights_1d * dJi, optimize='greedy')
                    CP_IP_3dM[i] = CP_IP_3dM_i_
                    type_cache[typeWr2Metric] = CP_IP_3dM_i_

            else:
                w = RMw[i]
                u, v = RMu[i]
                a, b = RMe[i]
                dJi = detJ[i]
                # so, w0 = [0 0 w]^T, u1 = [u, v, 0]^T, e1 = [a b 0]^T, A = w0 X u1 = [-wv wu 0]^T, (A, e1) = -wva + wub
                CP_IP_3dM[i] = - np.einsum('li, lj, lk, l -> ijk', w, v, a, quad_weights_1d * dJi, optimize='greedy')\
                               + np.einsum('li, lj, lk, l -> ijk', w, u, b, quad_weights_1d * dJi, optimize='greedy')

        self._CP_IP_3dM_ = CP_IP_3dM
        self._w0_ = w0
        self._freeze_self_()

    def __call__(self, i):
        """return 2d matrix of output = '1-M-2' type for mesh-element #i."""
        M = np.einsum('ijk, i -> kj', self._CP_IP_3dM_[i], self._w0_.cochain.local[i], optimize='greedy')
        return csc_matrix(M)