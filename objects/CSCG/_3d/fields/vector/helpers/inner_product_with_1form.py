
from screws.freeze.base import FrozenOnly
import numpy as np
from scipy import sparse as spspa




class _VF_InnerWith1Form(FrozenOnly):
    def __init__(self, vf, _1f, quad_degree):
        if quad_degree is None: quad_degree = [_1f.dqp[i]+1 for i in range(3)]
        quad_nodes, _, quad_weights = _1f.space.___PRIVATE_do_evaluate_quadrature___(quad_degree)
        _, bf1 = _1f.do.evaluate_basis_at_meshgrid(*quad_nodes, compute_xietasigma=False)
        self._g0_, self._g1_, self._g2_ = bf1
        self._JM_ = _1f.mesh.elements.coordinate_transformation.QUAD_1d.Jacobian_matrix(quad_degree, 'Gauss')
        self._mapping_ = _1f.mesh.elements.coordinate_transformation.QUAD_1d.mapping(quad_degree, 'Gauss')
        self._qw_ = quad_weights
        self._mesh_ = _1f.mesh
        self._vf_ = vf
        self.RESET_cache()
        self._freeze_self_()

    def RESET_cache(self):
        self._J_cache_ = dict()

    def ___PRIVATE_J___(self, i, mark):
        """"""
        if mark in self._J_cache_:
            return self._J_cache_[mark]
        else:
            JM = self._JM_[i]
            if isinstance(mark, str) and mark[:4] == 'Orth':
                J00 = JM[1][1] * JM[2][2]
                J01 = None
                J02 = None
                J10 = None
                J11 = JM[2][2] * JM[0][0]
                J12 = None
                J20 = None
                J21 = None
                J22 = JM[0][0] * JM[1][1]
            else:
                J00 = JM[1][1]*JM[2][2] - JM[1][2]*JM[2][1]
                J01 = JM[2][1]*JM[0][2] - JM[2][2]*JM[0][1]
                J02 = JM[0][1]*JM[1][2] - JM[0][2]*JM[1][1]
                J10 = JM[1][2]*JM[2][0] - JM[1][0]*JM[2][2]
                J11 = JM[2][2]*JM[0][0] - JM[2][0]*JM[0][2]
                J12 = JM[0][2]*JM[1][0] - JM[0][0]*JM[1][2]
                J20 = JM[1][0]*JM[2][1] - JM[1][1]*JM[2][0]
                J21 = JM[2][0]*JM[0][1] - JM[2][1]*JM[0][0]
                J22 = JM[0][0]*JM[1][1] - JM[0][1]*JM[1][0]
            J = (J00, J01, J02, J10, J11, J12, J20, J21, J22)
            self._J_cache_[mark] = J
            return J

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
        J00, J01, J02, J10, J11, J12, J20, J21, J22 = self.___PRIVATE_J___(i, mark)
        if isinstance(mark, str) and mark[:4] == 'Orth':
            v0 = np.einsum('w, iw -> i', f0 * J00 * self._qw_, g0, optimize='greedy')
            v1 = np.einsum('w, iw -> i', f1 * J11 * self._qw_, g1, optimize='greedy')
            v2 = np.einsum('w, iw -> i', f2 * J22 * self._qw_, g2, optimize='greedy')
        else:
            v0 = np.einsum('w, iw -> i', (f0*J00 + f1*J01 + f2*J02) * self._qw_, g0, optimize='greedy')
            v1 = np.einsum('w, iw -> i', (f0*J10 + f1*J11 + f2*J12) * self._qw_, g1, optimize='greedy')
            v2 = np.einsum('w, iw -> i', (f0*J20 + f1*J21 + f2*J22) * self._qw_, g2, optimize='greedy')
        RETURN = spspa.csr_matrix(np.concatenate([v0, v1, v2])).T

        return RETURN



