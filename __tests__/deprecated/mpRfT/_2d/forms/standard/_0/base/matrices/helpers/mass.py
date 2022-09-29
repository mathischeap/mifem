# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/05/25 5:39 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from numpy import einsum
from scipy.sparse import csr_matrix

class MassMatrixGenerator(FrozenOnly):
    """"""
    def __init__(self, f):
        self._f_ = f
        mesh = self._f_.mesh
        self._coo_ = mesh.coo_map.Gauss(2)
        self._basis_ = mesh.space.do.evaluate_basis(self._f_, self._coo_)
        self._detJ_ = mesh.rcMC.Jacobian(self._basis_)
        self._freeze_self_()

    def __call__(self, rp):
        """"""
        basis = self._basis_[rp][1][0]
        detJ = self._detJ_[rp]
        weights = self._coo_[rp][1][2]
        Mi = einsum('im, jm, m -> ij', basis, basis, detJ * weights, optimize='greedy')
        Mi = csr_matrix(Mi)
        return Mi

if __name__ == "__main__":
    # mpiexec -n 4 python objects/mpRfT/_2d/forms/standard/_0/base/matrices/helpers/mass.py
    pass
