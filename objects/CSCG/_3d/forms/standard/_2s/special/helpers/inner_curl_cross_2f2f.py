# -*- coding: utf-8 -*-
"""
Yi Zhang
zhangyi_aero@hotmail.com
created at: 2/22/2023 12:36 PM
"""
from components.freeze.base import FrozenOnly
import numpy as np

from tools.elementwiseCache.dataStructures.objects.multiDimMatrix.main import MultiDimMatrix
from importlib import import_module
from scipy.sparse import csc_matrix, bmat, csr_matrix

_global_cache = {
    'key': '',
    'data': dict(),
}


class _InnerCurlCrossProduct_2f2f(FrozenOnly):
    """"""

    def __init__(self, u2, B2, b2, quad_degree=None):
        """( curl(u2 X B2), b1 )."""
        assert u2.ndim == B2.ndim == b2.ndim, " <_InnerCurlCrossProduct_2f2f> "
        assert b2.k == B2.k == u2.k == 2, " <_InnerCurlCrossProduct_2f2f> "
        assert u2.mesh == B2.mesh, "_InnerCurlCrossProduct_2f2f: Meshes do not match."
        assert u2.mesh == b2.mesh, "_InnerCurlCrossProduct_2f2f: Meshes do not match."

        key = "{}{}{}{}".format(id(u2), id(B2), id(b2), quad_degree)
        if key == _global_cache['key']:
            _3dM = _global_cache['data']

        else:
            if quad_degree is None:
                quad_degree = [int(np.max([u2.dqp[i], B2.dqp[i], b2.dqp[i]]) * 3) for i in range(3)]
            else:
                pass

            base_path = '.'.join(str(self.__class__).split(' ')[1][1:-2].split('.')[:-5]) + '.'

            up = u2.space.p
            Bp = B2.space.p

            p = [up[0]+Bp[0], up[1]+Bp[1], up[2]+Bp[2]]
            new_space = u2.space.__class__(p, None)

            _1form_class = getattr(import_module(base_path + '_1s.main'), '_3dCSCG_1Form')

            f1 = _1form_class(u2.mesh, new_space, hybrid=False)
            f2 = u2.__class__(u2.mesh, new_space, hybrid=False)

            cpm = u2.special.cross_product_2f__ip_1f(B2, f1, quad_degree=quad_degree, output='MDM')

            invM = f1.matrices.mass.inv

            E21 = f1.matrices.incidence

            quad_nodes, _, quad_weights = f2.space.___PRIVATE_do_evaluate_quadrature___(
                f2.dqp,
                quad_type='Gauss',
            )
            xietasigma, bfSelf = f2.do.evaluate_basis_at_meshgrid(*quad_nodes)
            bfOther = b2.do.evaluate_basis_at_meshgrid(*quad_nodes)[1]

            mesh = u2.mesh
            elements = mesh.elements
            _3dM = dict()
            type_cache = dict()

            for i in elements:
                element = elements[i]
                typeWr2Metric = element.type_wrt_metric.mark
                if isinstance(typeWr2Metric, str):

                    if typeWr2Metric in type_cache:
                        _3dM[i] = type_cache[typeWr2Metric]

                    else:

                        Mi = self.___inner___(element, xietasigma, quad_weights, bfSelf, bfOther)  # (b, f2)

                        invMi = invM[i].toarray()   # (f1, f1)
                        cpm_i = cpm[i]  # (u, B, f1)

                        E21i = E21[i].toarray()   # (f2, f1)
                        _ = np.einsum('ij, mnj -> imn', invMi, cpm_i, optimize='greedy')  # f1 u B
                        _ = np.einsum('ji, imn -> jmn', E21i, _, optimize='greedy')   # f2 u B
                        _ = np.einsum('kj, jmn -> mnk', Mi, _, optimize='greedy')     # u B b

                        # print(_.__class__.__name__, _.shape)

                        _3dM[i] = _
                        type_cache[typeWr2Metric] = _

                else:
                    raise NotImplementedError()

            _global_cache['key'] = key
            _global_cache['data'] = _3dM

        self._3dM = _3dM
        self._u = u2
        self._B = B2
        self._b = b2
        self._freeze_self_()

    def ___inner___(self, element, xietasigma, quad_weights, bfSelf, bfOther):
        """"""
        typeWr2Metric = element.type_wrt_metric.mark
        J = element.coordinate_transformation.Jacobian_matrix(*xietasigma)
        sqrtg = element.coordinate_transformation.Jacobian(*xietasigma, J=J)
        iJ = element.coordinate_transformation.inverse_Jacobian_matrix(*xietasigma, J=J)
        g = element.coordinate_transformation.inverse_metric_matrix(*xietasigma, iJ=iJ)
        if typeWr2Metric[:4] == 'Orth':
            M00 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg, g[1][1]*g[2][2], bfOther[0], bfSelf[0])
            M11 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg, g[2][2]*g[0][0], bfOther[1], bfSelf[1])
            M22 = self.___PRIVATE_inner_H1___(quad_weights, sqrtg, g[0][0]*g[1][1], bfOther[2], bfSelf[2])
            M01 = None
            M02 = None
            M12 = None
            M10 = None
            M20 = None
            M21 = None
        else:
            raise NotImplementedError()

        Mi = bmat([(M00, M01, M02),
                   (M10, M11, M12),
                   (M20, M21, M22)], format='csc')

        return Mi.toarray()

    @staticmethod
    def ___PRIVATE_inner_H1___(quad_weights, sqrt_g, g, bfO, bfS):
        M = np.einsum('m, im, jm -> ij', quad_weights * sqrt_g * g, bfO, bfS, optimize='optimal')
        return csc_matrix(M)

    def _2_M_1_(self, i):
        """"""
        M = np.einsum('ijk, i -> kj', self._3dM[i], self._u.cochain.local[i], optimize='greedy')
        return csr_matrix(M)

    def _2_M_0_(self, i):
        """"""
        M = np.einsum('ijk, j -> ki', self._3dM[i], self._B.cochain.local[i], optimize='greedy')
        return csr_matrix(M)

    @property
    def MDM(self):
        return MultiDimMatrix(self._u.mesh.elements, self._3dM, [self._u, self._B, self._b], 'no_cache')
