# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/29 8:39 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from components.freeze.base import FrozenOnly
from components.quadrature import Quadrature
from objects.miUsGrid.triangular.forms.standard._0.base.matrices.helpers.inner import ___MassMatrix_Inner___
from tools.linearAlgebra.elementwiseCache.objects.sparseMatrix.main import EWC_SparseMatrix


class miUs_Triangular_S0F_Matrices(FrozenOnly):
    """"""

    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self._freeze_self_()

    @property
    def incidence(self):
        return self._sf_.coboundary.incidence_matrix

    @property
    def mass(self):
        """"""
        quad_degree = [self._sf_.space.p + 1 for _ in range(2)]
        Quad = Quadrature(quad_degree, category='Gauss')
        quad_nodes = Quad.quad[0]
        quad_weights_ravel = Quad.quad_ndim_ravel[-1]
        xi_eta, bf = self._sf_.do.evaluate_basis_at_meshgrid(*quad_nodes)
        DG = ___MassMatrix_Inner___(self._sf_, xi_eta, quad_weights_ravel, bf)
        return EWC_SparseMatrix(self._sf_.mesh.elements, DG,
                                self._sf_.mesh.elements.___Pr_EWC_cache_key___)




if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/triangular/forms/standard/_0/base/matrices/main.py
    pass
