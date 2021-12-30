# -*- coding: utf-8 -*-
"""

@author: Yi Zhang.
         Department of Aerodynamics
         Faculty of Aerospace Engineering
         TU Delft, Delft, Netherlands

"""
from _3dCSCG.space.main import _3dCSCG_Space
from SCREWS.quadrature import Quadrature
from root.config import *

from BASE.CSCG.space._1d_basis.polynomials import _1dPolynomial





class _3dCSCG_PolynomialSpace(_3dCSCG_Space):
    """"""
    def __init__(self, inputs, ndim):
        """
        INPUTS for all spaces must be the same as shown here.

        If `ndim` is not None, we will repeat the inputs for `ndim` times.

        `ndim` to be None or 3!
        """
        self.___1D_basis___ = _1dPolynomial
        self._quadrature_cache_ = [-1, None, None, None, None]
        super().__init__(inputs, ndim)




    def DO_evaluate_quadrature(self, quad_degree, quad_type=None):
        """
        We only cache the results of the last call.

        :param quad_degree:
        :param quad_type:
        :return:
        """
        if quad_type is None: quad_type = 'Gauss'

        if [quad_degree, quad_type] == self._quadrature_cache_[:2]:
            pass
        else:
            assert np.shape(quad_degree) == (3,), " <Polynomials> "
            _Quadrature_ = Quadrature(quad_degree, category=quad_type)
            quad_nodes, quad_weights = _Quadrature_.quad
            quad_weights_ravel = _Quadrature_.quad_ndim_ravel[-1]
            # return quad_nodes, quad_weights, quad_weights_ravel
            self._quadrature_cache_ = [quad_degree, quad_type,
                                       quad_nodes,
                                       quad_weights,
                                       quad_weights_ravel]

        return self._quadrature_cache_[2:]







if __name__ == "__main__":
    space = _3dCSCG_PolynomialSpace([('Lobatto', 3), ('Lobatto', 4), ('Lobatto', 5)], None)
    print(space.p)
    space = _3dCSCG_PolynomialSpace([([-1,0,1],), ('Lobatto', 4), ('Lobatto', 5)], None)
    print(space.p)
    space = _3dCSCG_PolynomialSpace(([-1,0,1],), 3) # ndim = 3
    print(space.p)
    space = _3dCSCG_PolynomialSpace([([-1,0,0.5,1],), ([-1,-0.5,0,0.5,1],), ([-1,1],)], None) # ndim = 3
    print(space.p)
    space = _3dCSCG_PolynomialSpace((2,3,4), None) # ndim = 3
    print(space.p)
    # space = _3dCSCG_PolynomialSpace([3,3,3], None) # ndim = 3
    # print(space.p)