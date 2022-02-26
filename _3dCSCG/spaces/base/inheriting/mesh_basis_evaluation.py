# -*- coding: utf-8 -*-
"""
INTRO

Yi Zhang (C)
Created on Mon May  6 19:37:11 2019
Aerodynamics, AE
TU Delft
"""
import numpy as np
from screws.exceptions import DimensionError
from screws.decorators import memoize1

# noinspection PyUnresolvedReferences
class EvaluatingMeshBasis:
    def DO_evaluate_form_basis_at_meshgrid(self, k, *domain, compute_xietasigma=True):
        """
        Parameters
        ---------
        k : int
        domain : tuple
            The domain we are going to evaluate basis. Notice that here we only
            accept tuple. So even for 1-D Polynomials, we have to put xi in a
            tuple: (xi,).
        compute_xietasigma :
        """
        assert 0 <= k <= self.ndim, " <Polynomials> : k={} is wrong.".format(k)
        assert len(domain) == self.ndim, \
            " <Polynomials> : domain shape={} wrong.".format(np.shape(domain))
        for i in range(self.ndim):
            assert domain[i].__class__.__name__ in ('list', 'ndarray'), \
                " <Polynomials> : domain[{}].type={} is wrong.".format(
                        i, domain[i].__class__.__name__)
            assert np.ndim(domain[i]) == 1, \
                " <Polynomials> : ndim(domain[{}])={} is wrong.".format(
                        i, np.ndim(domain[i]))
            if np.size(domain[i]) > 1:
                assert np.all(np.diff(domain[i]) > 0) and np.max(domain[i]) <= 1 and np.min(domain[i]) >= -1, \
                    " <Polynomials> : domain[i]={} wrong, need to be increasing and bounded in [-1, 1].".format(
                            domain[i])
            else:
                pass

        if self.ndim == 3:
            if compute_xietasigma:
                _xietasigma_ = self.___PRIVATE_evaluate_basis_functions_meshgrid_3d___(domain)
            if k == 0:
                bf_xi = self.basises[0].node_basis(x=domain[0])
                bf_et = self.basises[1].node_basis(x=domain[1])
                bf_si = self.basises[2].node_basis(x=domain[2])
                bf = np.kron(np.kron(bf_si, bf_et), bf_xi)

                _basis_ = (bf,)

            elif k == 1:
                ed_xi = self.basises[0].edge_basis(x=domain[0])
                lb_et = self.basises[1].node_basis(x=domain[1])
                lb_si = self.basises[2].node_basis(x=domain[2])
                bf_edge_dxi = np.kron(np.kron(lb_si, lb_et), ed_xi)

                lb_xi = self.basises[0].node_basis(x=domain[0])
                ed_et = self.basises[1].edge_basis(x=domain[1])
                lb_si = self.basises[2].node_basis(x=domain[2])
                bf_edge_det = np.kron(np.kron(lb_si, ed_et), lb_xi)

                lb_xi = self.basises[0].node_basis(x=domain[0])
                lb_et = self.basises[1].node_basis(x=domain[1])
                ed_si = self.basises[2].edge_basis(x=domain[2])
                bf_edge_dsi = np.kron(np.kron(ed_si, lb_et), lb_xi)

                _basis_ = (bf_edge_dxi, bf_edge_det, bf_edge_dsi)

            elif k == 2:
                lb_xi = self.basises[0].node_basis(x=domain[0])
                ed_et = self.basises[1].edge_basis(x=domain[1])
                ed_si = self.basises[2].edge_basis(x=domain[2])
                bf_face_det_dsi = np.kron(np.kron(ed_si, ed_et), lb_xi)

                ed_xi = self.basises[0].edge_basis(x=domain[0])
                lb_et = self.basises[1].node_basis(x=domain[1])
                ed_si = self.basises[2].edge_basis(x=domain[2])
                bf_face_dsi_dxi = np.kron(np.kron(ed_si, lb_et), ed_xi)

                ed_xi = self.basises[0].edge_basis(x=domain[0])
                ed_et = self.basises[1].edge_basis(x=domain[1])
                lb_si = self.basises[2].node_basis(x=domain[2])
                bf_face_dxi_det = np.kron(np.kron(lb_si, ed_et), ed_xi)

                _basis_ = (bf_face_det_dsi, bf_face_dsi_dxi, bf_face_dxi_det)

            elif k == 3:
                bf_xi = self.basises[0].edge_basis(x=domain[0])
                bf_et = self.basises[1].edge_basis(x=domain[1])
                bf_si = self.basises[2].edge_basis(x=domain[2])
                bf = np.kron(np.kron(bf_si, bf_et), bf_xi)

                _basis_ = (bf,)
            else:
                raise DimensionError()
        else:
            raise DimensionError()
        
        if compute_xietasigma:
            # noinspection PyUnboundLocalVariable
            return _xietasigma_, _basis_
        else:
            return None, _basis_

    @staticmethod
    def ___PRIVATE_evaluate_basis_functions_meshgrid_3d___(domain):
        xi, eta, sigma = np.meshgrid(*domain, indexing='ij')
        xi = xi.ravel('F')
        eta = eta.ravel('F')
        sigma = sigma.ravel('F')
        return xi, eta, sigma

    @memoize1
    def DO_evaluate_form_basis_at_quadrature(
        self, k, quad_degree, quad_type=None, compute_xietasigma=True):
        """"""
        quad_nodes, quad_weights, quad_weights_ravel = \
            self.DO_evaluate_quadrature(quad_degree, quad_type=quad_type)
        _xietasigma_, _basis_ = \
            self.DO_evaluate_form_basis_at_meshgrid(
                k, *quad_nodes, compute_xietasigma=compute_xietasigma)
        return _xietasigma_, _basis_, quad_weights, quad_weights_ravel