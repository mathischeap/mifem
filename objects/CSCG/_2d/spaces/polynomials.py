# -*- coding: utf-8 -*-


from objects.CSCG._2d.spaces.base.main import _2dCSCG_Space
from screws.quadrature import Quadrature
from root.config.main import *
from objects.CSCG.base.spaces._1d_basis.polynomials import _1dPolynomial





class _2dCSCG_PolynomialSpace(_2dCSCG_Space):
    """"""
    def __init__(self, inputs, ndim):
        self.___1D_basis___ = _1dPolynomial
        self._quadrature_cache_ = [-1, None, None, None, None]
        self._r = None
        super().__init__(inputs, ndim)


    def ___PRIVATE_do_evaluate_quadrature___(self, quad_degree, quad_type=None):
        """
        We only do cache the results for last call.

        :param quad_degree:
        :param quad_type:
        :return:
        """
        if quad_type is None: quad_type = 'Gauss'

        if [quad_degree, quad_type] == self._quadrature_cache_[:2]:
            pass
        else:
            assert np.shape(quad_degree) == (2,), " <Polynomials> "
            _Quadrature_ = Quadrature(quad_degree, category=quad_type)
            quad_nodes, quad_weights = _Quadrature_.quad
            quad_weights_ravel = _Quadrature_.quad_ndim_ravel[-1]
            # return quad_nodes, quad_weights, quad_weights_ravel
            self._quadrature_cache_ = [quad_degree, quad_type,
                                       quad_nodes, quad_weights, quad_weights_ravel]

        return self._quadrature_cache_[2:]





if __name__ == "__main__":
    pass

