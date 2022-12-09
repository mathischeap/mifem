# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 7/5/2022 10:45 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from components.freeze.base import FrozenOnly
from objects.CSCG._3d.forms.standard._2s.main import _3dCSCG_2Form
import numpy as np
from scipy.sparse import csr_matrix


class ___3dCSCG_curl1_CrossProduct_1__ip_2___(FrozenOnly):
    """"""

    def __init__(self, u, e, quad_degree=None):
        """
        w = curl u  or w = d(u) is a 2-form, and e is a 2-form. The cochain of u must be given.

        We compute (w X u, e). As the cochain of u is given, we can obtain a vector whose local vectors
        refer to the local cochains of e.

        Parameters
        ----------
        u
        e
        quad_degree
        """
        assert u.ndim == e.ndim, " <___3dCSCG_curl1_CrossProduct_1__ip_2___> "
        assert u.k+1 == e.k == 2, " <___3dCSCG_curl1_CrossProduct_1__ip_2___> "
        assert u.mesh == e.mesh, "___3dCSCG_curl1_CrossProduct_1__ip_2___: Meshes do not match."

        w = _3dCSCG_2Form(u.mesh, u.space, hybrid=u.whether.hybrid)

        if quad_degree is None:
            quad_degree = [int(np.max([u.dqp[i], e.dqp[i]])) * 2 for i in range(3)]
        else:
            pass

        quad_nodes, _, quad_weights_1d = \
            u.space.___PRIVATE_do_evaluate_quadrature___(quad_degree)

        RMw = w.do.make_reconstruction_matrix_on_grid(*quad_nodes)
        RMu = u.do.make_reconstruction_matrix_on_grid(*quad_nodes)
        RMe = e.do.make_reconstruction_matrix_on_grid(*quad_nodes)

        xi, et, sg = np.meshgrid(*quad_nodes, indexing='ij')
        xi = xi.ravel('F')
        et = et.ravel('F')
        sg = sg.ravel('F')
        detJ = u.mesh.elements.coordinate_transformation.Jacobian(xi, et, sg)

        CP_IP_3dM = dict()
        type_cache = dict()

        for i in RMw: # go through all local mesh-elements
            typeWr2Metric = u.mesh.elements[i].type_wrt_metric.mark
            if isinstance(typeWr2Metric, str):
                if typeWr2Metric in type_cache:
                    CP_IP_3dM[i] = type_cache[typeWr2Metric]
                else:
                    # W = [wx wy, wz]^T    U= [U V W]^T  e2= [a b c]^T
                    # W X U = [wy*w - wz*v,   wz*u - wx*w,   wx*v - wy*u]^T = [A B C]^T
                    # W X U dot e2 = Aa + Bb + Cc
                    #              =  wy*w*a - wz*v*a + wz*u*b - wx*w*b + wx*v*c - wy*u*c
                    wx, wy, wz = RMw[i]
                    U, V, W = RMu[i]
                    a, b, c = RMe[i]
                    dJi = detJ[i]
                    CP_IP_3dM_i_ = + np.einsum('li, lj, lk, l -> ijk', wy, W, a, quad_weights_1d * dJi, optimize='optimal')\
                                   - np.einsum('li, lj, lk, l -> ijk', wz, V, a, quad_weights_1d * dJi, optimize='optimal')\
                                   + np.einsum('li, lj, lk, l -> ijk', wz, U, b, quad_weights_1d * dJi, optimize='optimal')\
                                   - np.einsum('li, lj, lk, l -> ijk', wx, W, b, quad_weights_1d * dJi, optimize='optimal')\
                                   + np.einsum('li, lj, lk, l -> ijk', wx, V, c, quad_weights_1d * dJi, optimize='optimal')\
                                   - np.einsum('li, lj, lk, l -> ijk', wy, U, c, quad_weights_1d * dJi, optimize='optimal')

                    CP_IP_3dM[i] = CP_IP_3dM_i_
                    type_cache[typeWr2Metric] = CP_IP_3dM_i_

            else:
                wx, wy, wz = RMw[i]
                U, V, W = RMu[i]
                a, b, c = RMe[i]
                dJi = detJ[i]
                CP_IP_3dM[i] = + np.einsum('li, lj, lk, l -> ijk', wy, W, a, quad_weights_1d * dJi, optimize='optimal')\
                               - np.einsum('li, lj, lk, l -> ijk', wz, V, a, quad_weights_1d * dJi, optimize='optimal')\
                               + np.einsum('li, lj, lk, l -> ijk', wz, U, b, quad_weights_1d * dJi, optimize='optimal')\
                               - np.einsum('li, lj, lk, l -> ijk', wx, W, b, quad_weights_1d * dJi, optimize='optimal')\
                               + np.einsum('li, lj, lk, l -> ijk', wx, V, c, quad_weights_1d * dJi, optimize='optimal')\
                               - np.einsum('li, lj, lk, l -> ijk', wy, U, c, quad_weights_1d * dJi, optimize='optimal')

        self._CP_IP_3dM_ = CP_IP_3dM
        self._u_ = u
        self._E21_ = u.matrices.incidence
        self._freeze_self_()

    def __call__(self, i):
        """return a vector of output = 'e' type for mesh-element #i."""
        uLCC = self._u_.cochain.local[i]
        V = np.einsum('ijk, i, j -> k', self._CP_IP_3dM_[i], self._E21_[i] @ uLCC, uLCC, optimize='optimal')
        return csr_matrix(V).T



if __name__ == '__main__':
    # mpiexec -n 4 python objects/CSCG/_3d/forms/standard/_1s/special/helpers/curl1_cross_product_1__ip_2.py
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

    f1 = FC('1-f', is_hybrid=False)
    u2 = FC('2-f', is_hybrid=False)


    f1.TW.func.do.set_func_body_as(velocity)
    f1.TW.current_time = 0
    f1.TW.do.push_all_to_instant()
    f1.discretize()

    CP = f1.special.curl_self_cross_product_self__ip_2f(u2)


    E21 = f1.matrices.incidence
    E12 = E21.T

    V1 = E12 @ CP

    for i in V1:
        print(V1[i].shape)
