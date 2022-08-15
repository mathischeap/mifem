# -*- coding: utf-8 -*-
"""
INTRO

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
import sys
if './' not in sys.path: sys.path.append('./')
import numpy as np
from scipy.sparse import csr_matrix
from screws.freeze.base import FrozenOnly


class ___3dCSCG_2Form_CrossProduct_2__ip_1___(FrozenOnly):
    """
    The class for the inner wedge matrix; representing :math:`(a \\times b, e)`.

    u, B are both 2-forms, j is 1-form.

    :param u:
    :param B:
    :param j:
    :param quad_degree:
    """
    def __init__(self, u, B, j, quad_degree=None):
        assert u.ndim == B.ndim == j.ndim, " <___3dCSCG_2Form_CrossProduct_2__ip_1___> "
        assert u.k == B.k == j.k+1 == 2, " <___3dCSCG_2Form_CrossProduct_2__ip_1___> "
        assert u.mesh == B.mesh, "___3dCSCG_2Form_CrossProduct_2__ip_1___, Meshes do not match."
        assert j.mesh == B.mesh, "___3dCSCG_2Form_CrossProduct_2__ip_1___, Meshes do not match."

        if quad_degree is None:
            quad_degree = [int(np.max([u.dqp[i], B.dqp[i], j.dqp[i]]) * 1.5)  for i in range(3)]

        quad_nodes, _, quad_weights_1d = \
            u.space.___PRIVATE_do_evaluate_quadrature___(quad_degree)

        RMB = B.do.make_reconstruction_matrix_on_grid(*quad_nodes)
        RMu = u.do.make_reconstruction_matrix_on_grid(*quad_nodes)
        RMe = j.do.make_reconstruction_matrix_on_grid(*quad_nodes)

        xi, et, sg = np.meshgrid(*quad_nodes, indexing='ij')
        xi = xi.ravel()
        et = et.ravel()
        sg = sg.ravel()
        detJ = u.mesh.elements.coordinate_transformation.Jacobian(xi, et, sg)

        CP_IP_3dM = dict()
        type_cache = dict()
        for i in RMB: # go through all local mesh-elements
            typeWr2Metric = B.mesh.elements[i].type_wrt_metric.mark
            if isinstance(typeWr2Metric, str):
                if typeWr2Metric in type_cache:
                    CP_IP_3dM[i] = type_cache[typeWr2Metric]
                else:
                    wx, wy, wz = RMu[i]
                    U, v, w = RMB[i]
                    a, b, c = RMe[i]
                    dJi = detJ[i]
                    # w1 = [wx wy, wz]^T    u2= [u v w]^T   e2= [a b c]^T
                    # WXU = w1 X u2 = [wy*w - wz*v,   wz*u - wx*w,   wx*v - wy*u]^T = [A B C]^T
                    # WXU dot e2 = Aa + Bb + Cc
                    Aa = + np.einsum('li, lj, lk, l -> ijk', wy, w, a, quad_weights_1d * dJi, optimize='greedy')\
                         - np.einsum('li, lj, lk, l -> ijk', wz, v, a, quad_weights_1d * dJi, optimize='greedy')
                    Bb = + np.einsum('li, lj, lk, l -> ijk', wz, U, b, quad_weights_1d * dJi, optimize='greedy')\
                         - np.einsum('li, lj, lk, l -> ijk', wx, w, b, quad_weights_1d * dJi, optimize='greedy')
                    Cc = + np.einsum('li, lj, lk, l -> ijk', wx, v, c, quad_weights_1d * dJi, optimize='greedy')\
                         - np.einsum('li, lj, lk, l -> ijk', wy, U, c, quad_weights_1d * dJi, optimize='greedy')
                    CP_IP_3dM_i_ = Aa + Bb + Cc
                    CP_IP_3dM[i] = CP_IP_3dM_i_
                    type_cache[typeWr2Metric] = CP_IP_3dM_i_

            else:
                wx, wy, wz = RMu[i]
                U, v, w = RMB[i]
                a, b, c = RMe[i]
                dJi = detJ[i]
                # w1 = [wx wy, wz]^T    u2= [u v w]^T   e2= [a b c]^T
                # WXU = w1 X u2 = [wy*w - wz*v,   wz*u - wx*w,   wx*v - wy*u]^T = [A B C]^T
                # WXU dot e2 = Aa + Bb + Cc
                Aa = + np.einsum('li, lj, lk, l -> ijk', wy, w, a, quad_weights_1d * dJi, optimize='greedy')\
                     - np.einsum('li, lj, lk, l -> ijk', wz, v, a, quad_weights_1d * dJi, optimize='greedy')
                Bb = + np.einsum('li, lj, lk, l -> ijk', wz, U, b, quad_weights_1d * dJi, optimize='greedy')\
                     - np.einsum('li, lj, lk, l -> ijk', wx, w, b, quad_weights_1d * dJi, optimize='greedy')
                Cc = + np.einsum('li, lj, lk, l -> ijk', wx, v, c, quad_weights_1d * dJi, optimize='greedy')\
                     - np.einsum('li, lj, lk, l -> ijk', wy, U, c, quad_weights_1d * dJi, optimize='greedy')

                CP_IP_3dM[i] = Aa + Bb + Cc

        self._CP_IP_3dM_ = CP_IP_3dM
        self._u_ = u
        self._freeze_self_()

    def __call__(self, i):
        """return 2d matrix of output = '1-M-2' type for mesh-element #i."""
        M = np.einsum('ijk, i -> kj', self._CP_IP_3dM_[i], self._u_.cochain.local[i], optimize='greedy')
        return csr_matrix(M)


if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_3d/forms/standard/_2s/special/helpers/cross_product_2__ip_1.py

    from objects.CSCG._3d.master import MeshGenerator, SpaceInvoker, FormCaller#, ExactSolutionSelector

    mesh = MeshGenerator('cuboid', region_layout=[2,2,2])([3,3,3])
    space = SpaceInvoker('polynomials')((2,2,2))
    FC = FormCaller(mesh, space)

    def u(t,x,y,z): return np.sin(np.pi*x)*np.cos(2*np.pi*y)*np.cos(np.pi*z) + t
    def v(t,x,y,z): return np.cos(np.pi*x)*np.sin(np.pi*y)*np.cos(2*np.pi*z) + t
    def w(t,x,y,z): return np.cos(np.pi*x)*np.cos(np.pi*y)*np.sin(2*np.pi*z) + t

    velocity = FC('vector', (u,v,w))
    U = FC('scalar', u)
    V = FC('scalar', v)
    W = FC('scalar', w)

    u = FC('2-f', is_hybrid=False)
    B = FC('2-f', is_hybrid=False)
    j = FC('1-f', is_hybrid=False)


    u.TW.func.do.set_func_body_as(velocity)
    u.TW.current_time = 0
    u.TW.do.push_all_to_instant()
    u.discretize()

    CM = u.special.cross_product_2f__ip_1f(B, j)
