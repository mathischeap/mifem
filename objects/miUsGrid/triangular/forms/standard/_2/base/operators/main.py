# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 10/2/2022 9:09 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
from screws.quadrature import Quadrature
from objects.miUsGrid.triangular.forms.standard._2.base.operators.helpers.inner import ___Operators_Inner___
from tools.linearAlgebra.elementwiseCache.objects.sparseMatrix.main import EWC_SparseMatrix


class miUs_Triangular_S2F_Operators(FrozenOnly):
    """"""

    def __init__(self, sf):
        """"""
        self._sf_ = sf
        self._freeze_self_()

    def inner(self, other):
        """"""
        quad_degree = [self._sf_.space.p + 1 for _ in range(2)]
        Quad = Quadrature(quad_degree, category='Gauss')
        quad_nodes = Quad.quad[0]
        quad_weights_ravel = Quad.quad_ndim_ravel[-1]
        xi_eta, bf = self._sf_.do.evaluate_basis_at_meshgrid(*quad_nodes)
        _     , of =     other.do.evaluate_basis_at_meshgrid(*quad_nodes)
        DG = ___Operators_Inner___(self._sf_, xi_eta, quad_weights_ravel, bf, of)
        return EWC_SparseMatrix(self._sf_.mesh.elements, DG, 'no_cache')


if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
