# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/29 8:41 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from components.freeze.base import FrozenOnly
from components.quadrature import Quadrature
from objects.miUsGrid.triangular.forms.standard._1.outer.operators.helpers.inner import ___Operators_Inner___
from tools.elementwiseCache.dataStructures.objects.sparseMatrix.main import EWC_SparseMatrix


class miUs_Triangular_oS1F_Operators(FrozenOnly):
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
        ______, of =     other.do.evaluate_basis_at_meshgrid(*quad_nodes)
        DG = ___Operators_Inner___(self._sf_, xi_eta, quad_weights_ravel, bf, of)
        return EWC_SparseMatrix(self._sf_.mesh.elements, DG, 'no_cache')


if __name__ == "__main__":
    # mpiexec -n 4 python
    pass

