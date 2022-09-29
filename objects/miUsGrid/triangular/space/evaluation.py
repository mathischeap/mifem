# -*- coding: utf-8 -*-
"""
@author: Yi Zhang
@contact: zhangyi_aero@hotmail.com
@time: 2022/09/21 4:40 PM
"""
import sys

if './' not in sys.path: sys.path.append('./')

from screws.freeze.base import FrozenOnly
import numpy as np
from screws.quadrature import Quadrature


class miUsGrid_TriangularFunctionSpace_Evaluation(FrozenOnly):
    """"""

    def __init__(self, space):
        """"""
        self._space_ = space
        self._ndim_ = space.ndim
        self._quadrature_cache_ = [-1, '', None, None, None]
        self._freeze_self_()

    def quadrature(self, quad_degree):
        """
        We only do cache the results for last call. We must use Gauss quadrature.

        :param quad_degree:
        :return:
        """

        if [quad_degree, 'Gauss'] == self._quadrature_cache_[:2]:
            pass
        else:
            assert np.shape(quad_degree) == (2,), " <quadrature> "
            _Quadrature_ = Quadrature(quad_degree, category='Gauss')
            quad_nodes, quad_weights = _Quadrature_.quad
            quad_weights_ravel = _Quadrature_.quad_ndim_ravel[-1]
            # return quad_nodes, quad_weights, quad_weights_ravel
            self._quadrature_cache_ = [quad_degree, 'Gauss',
                                       quad_nodes, quad_weights, quad_weights_ravel]

        return self._quadrature_cache_[2:]


    @staticmethod
    def ___Pr_check_domain___(*domain):

        assert len(domain) == 2, \
            "domain shape={} wrong.".format(np.shape(domain))
        for i in range(2):
            assert domain[i].__class__.__name__ in ('list', 'ndarray'), \
                " domain[{}].type={} is wrong.".format(i, domain[i].__class__.__name__)
            assert np.ndim(domain[i]) == 1, \
                " ndim(domain[{}])={} is wrong.".format(i, np.ndim(domain[i]))
            if np.size(domain[i]) > 1:
                assert np.all(np.diff(domain[i]) > 0) and \
                       np.max(domain[i]) <= 1 and \
                       np.min(domain[i]) >= -1, \
                    " <2dCSCG Space> : domain[i]={} wrong, need to be " \
                    "increasing and bounded in [-1, 1].".format( domain[i])
            else:
                pass

    @staticmethod
    def ___Pr_compute_xi_eta___(xi, et):
        xi, et = np.meshgrid(xi, et, indexing='ij')
        xi = xi.ravel('F')
        et = et.ravel('F')
        return xi, et






    def miUsTriangular_S0F_Outer(self, xi, et):
        """"""
        self.___Pr_check_domain___(xi, et)
        node_basis_xi = self._space_._1db_.lagrange_basis(x=xi)
        node_basis_et = self._space_._1db_.lagrange_basis(x=et)
        basis_singular = np.sum(np.kron(node_basis_et[0,:], node_basis_xi),axis=0)[np.newaxis, :]
        basis_regular = np.kron(node_basis_et[1:,:], node_basis_xi)
        basis = np.vstack((basis_singular, basis_regular))
        return self.___Pr_compute_xi_eta___(xi, et), [basis,]

    def  miUsTriangular_S0F_Inner(self, *args, **kwargs):
        return self.miUsTriangular_S0F_Outer(*args, **kwargs)





    def miUsTriangular_S1F_Outer(self, xi, et):
        """"""
        self.___Pr_check_domain___(xi, et)
        node_basis_xi = self._space_._1db_.lagrange_basis(x=xi)
        node_basis_et = self._space_._1db_.lagrange_basis(x=et)

        edge_basis_xi = self._space_._1db_.edge_basis(x=xi)
        edge_basis_et = self._space_._1db_.edge_basis(x=et)


        basis_dy = np.kron(edge_basis_et, node_basis_xi)
        basis_dx = np.kron(node_basis_et[1:,:], edge_basis_xi)

        basis = (basis_dy, basis_dx)

        return self.___Pr_compute_xi_eta___(xi, et), basis


    def miUsTriangular_S1F_Inner(self, xi, et):
        """"""
        self.___Pr_check_domain___(xi, et)
        node_basis_xi = self._space_._1db_.lagrange_basis(x=xi)
        node_basis_et = self._space_._1db_.lagrange_basis(x=et)

        edge_basis_xi = self._space_._1db_.edge_basis(x=xi)
        edge_basis_et = self._space_._1db_.edge_basis(x=et)

        basis_dy = np.kron(edge_basis_et, node_basis_xi)
        basis_dx = np.kron(node_basis_et[1:,:], edge_basis_xi)

        basis = (basis_dx, basis_dy)

        return self.___Pr_compute_xi_eta___(xi, et), basis








    def miUsTriangular_S2F_Outer(self, xi, et):
        """"""
        self.___Pr_check_domain___(xi, et)
        edge_basis_xi = self._space_._1db_.edge_basis(x=xi)
        edge_basis_et = self._space_._1db_.edge_basis(x=et)
        basis = np.kron(edge_basis_et, edge_basis_xi)
        return self.___Pr_compute_xi_eta___(xi, et), [basis,]

    def  miUsTriangular_S2F_Inner(self, *args, **kwargs):
        return self.miUsTriangular_S2F_Outer(*args, **kwargs)



if __name__ == "__main__":
    # mpiexec -n 4 python 
    pass
