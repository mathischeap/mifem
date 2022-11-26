# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/25 5:41 PM
"""
import sys

if './' not in sys.path: sys.path.append('/')

from components.freeze.base import FrozenOnly
from numpy import einsum
from scipy.sparse import csr_matrix, bmat



class OuterMassMatrixGenerator(FrozenOnly):
    """"""
    def __init__(self, f):
        self._f_ = f
        mesh = self._f_.mesh
        self._mesh_ = mesh
        self._coo_ = mesh.coo_map.Gauss(2)
        self._basis_ = mesh.space.do.evaluate_basis(self._f_, self._coo_)
        self._sqrtg_ = mesh.rcMC.Jacobian(self._basis_)
        self._g_ = mesh.rcMC.inverse_metric_matrix(self._basis_)
        self._freeze_self_()

    def __call__(self, rp):
        """"""
        basis = self._basis_[rp][1]
        sqrtg = self._sqrtg_[rp]
        g = self._g_[rp]
        weights = self._coo_[rp][1][2]
        M00 = self.___Pr_mh1___(weights, sqrtg*g[1][1], basis[0], basis[0])
        M11 = self.___Pr_mh1___(weights, sqrtg*g[0][0], basis[1], basis[1])
        mark = self._mesh_[rp].type_wrt_metric.mark
        if isinstance(mark, str) and mark[:4] == 'Orth':
            M01 = None
            M10 = None
        else:
            M01 = self.___Pr_mh1___(weights, -sqrtg*g[1][0], basis[0], basis[1])
            M10 = M01.T
        Mi = bmat([(M00, M01),
                   (M10, M11)], format='csr')
        return Mi

    @staticmethod
    def ___Pr_mh1___(quad_weights, sqrtg_g, bfO, bfS):
        M = einsum('m, im, jm -> ij', quad_weights*sqrtg_g, bfO, bfS, optimize='greedy')
        return csr_matrix(M)




if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/standard/_1/base/matrices/helpers/outer_mass.py
    pass
