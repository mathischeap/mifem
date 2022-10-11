# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/20 2:16 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from objects.CSCG.base.spaces._1d_basis.polynomials import _1dPolynomial

from screws.freeze.main import FrozenClass
from objects.miUsGrid.triangular.space.num_basis import miUsGrid_TriangularFunctionSpace_NumBasis
from objects.miUsGrid.triangular.space.num_basis_components import miUsGrid_TriangularFunctionSpace_NumBasisComponents
from objects.miUsGrid.triangular.space.incidence_matrix import miUsGrid_TriangularFunctionSpace_IncidenceMatrix
from objects.miUsGrid.triangular.space.local_numbering import miUsGrid_TriangularFunctionSpace_LocalNumbering
from objects.miUsGrid.triangular.space.evaluation import miUsGrid_TriangularFunctionSpace_Evaluation



class miUsGrid_TriangularFunctionSpace(FrozenClass):
    """"""

    def __init__(self, p):
        """

        Parameters
        ----------
        p : int
            The degree of the basis functions (polynomials).

        """
        assert isinstance(p, int) and p >= 1, f"p={p} invalid, need a positive integer."
        self._p_ = p
        self._1db_ = _1dPolynomial(f'Lobatto-{p}')

        self._num_basis_ = miUsGrid_TriangularFunctionSpace_NumBasis(self)
        self._incidence_matrix_ = miUsGrid_TriangularFunctionSpace_IncidenceMatrix(self)
        self._local_numbering_ = miUsGrid_TriangularFunctionSpace_LocalNumbering(self)
        self._evaluation_ = miUsGrid_TriangularFunctionSpace_Evaluation(self)
        self._num_basis_components_ = miUsGrid_TriangularFunctionSpace_NumBasisComponents(self)
        self._freeze_self_()

    def __repr__(self):
        """"""
        return f"miUsGrid_TriangularFunctionSpace{self.p}"

    def __eq__(self, other):
        return other.__class__.__name__ == 'miUsGrid_TriangularFunctionSpace' and other.p == self.p

    @property
    def ndim(self):
        return 2

    @property
    def nodes(self):
        return self._1db_.nodes

    @property
    def p(self):
        """the degree of the basis functions."""
        return self._p_

    @property
    def num_basis(self):
        return self._num_basis_

    @property
    def num_basis_components(self):
        return self._num_basis_components_

    @property
    def evaluation(self):
        return self._evaluation_

    @property
    def incidence_matrix(self):
        return self._incidence_matrix_

    @property
    def local_numbering(self):
        return self._local_numbering_


if __name__ == "__main__":
    # mpiexec -n 4 python objects/miUsGrid/triangular/space/main.py
    pass
