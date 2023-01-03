# -*- coding: utf-8 -*-
from components.freeze.base import FrozenOnly
import numpy as np
from scipy import sparse as spspa


class _VF_InnerWith2Form(FrozenOnly):
    def __init__(self, vf, _2f, quad_degree):
        if quad_degree is None:
            quad_degree = [_2f.dqp[i]+1 for i in range(3)]
        quad_nodes, _, quad_weights = _2f.space.___PRIVATE_do_evaluate_quadrature___(quad_degree)
        _, bf2 = _2f.do.evaluate_basis_at_meshgrid(*quad_nodes, compute_xietasigma=False)
        self._g0_, self._g1_, self._g2_ = bf2
        self._JM_ = _2f.mesh.elements.coordinate_transformation.QUAD_1d.Jacobian_matrix(quad_degree, 'Gauss')
        self._mapping_ = _2f.mesh.elements.coordinate_transformation.QUAD_1d.mapping(quad_degree, 'Gauss')
        self._vf_ = vf
        self._qw_ = quad_weights
        self._mesh_ = _2f.mesh
        self._freeze_self_()

    def __call__(self, i):
        """
        :param i: # element.
        :return:
        """
        mark = self._mesh_.elements[i].type_wrt_metric.mark
        xyz = self._mapping_[i]
        _f0_, _f1_, _f2_ = self._vf_.do.evaluate_func_at_time()
        f0, f1, f2 = _f0_(*xyz), _f1_(*xyz), _f2_(*xyz)
        g0, g1, g2 = self._g0_, self._g1_, self._g2_
        JM = self._JM_[i]
        if isinstance(mark, str) and mark[:4] == 'Orth':
            v0 = np.einsum('w, iw -> i', f0 * JM[0][0] * self._qw_, g0, optimize='greedy')
            v1 = np.einsum('w, iw -> i', f1 * JM[1][1] * self._qw_, g1, optimize='greedy')
            v2 = np.einsum('w, iw -> i', f2 * JM[2][2] * self._qw_, g2, optimize='greedy')
        else:
            v0 = np.einsum('w, iw -> i', (f0*JM[0][0] + f1*JM[1][0] + f2*JM[2][0])*self._qw_, g0, optimize='greedy')
            v1 = np.einsum('w, iw -> i', (f0*JM[0][1] + f1*JM[1][1] + f2*JM[2][1])*self._qw_, g1, optimize='greedy')
            v2 = np.einsum('w, iw -> i', (f0*JM[0][2] + f1*JM[1][2] + f2*JM[2][2])*self._qw_, g2, optimize='greedy')
        RETURN = spspa.csr_matrix(np.concatenate([v0, v1, v2])).T

        return RETURN
