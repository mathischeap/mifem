# -*- coding: utf-8 -*-
"""
INTRO

Yi Zhang (C)
Created on Mon May  6 19:37:11 2019
Aerodynamics, AE
TU Delft
"""
import numpy as np
from SCREWS.exceptions import DimensionError
from SCREWS.decorators import accepts

# noinspection PyUnresolvedReferences
class EvaluatingTraceBasis:
    """"""
    @accepts('self', int)
    def DO_evaluate_trace_basis_at_meshgrid(self, k, *domain, compute_xietasigma=True):
        """ """
        assert 0 <= k <= self.ndim-1, " <Polynomials> : k={} is wrong.".format(k)
        assert len(domain) == self.ndim, \
            " <Polynomials> : domain shape={} wrong.".format(np.shape(domain))
        for i in range(self.ndim):
            assert domain[i].__class__.__name__ in ('list', 'ndarray'), \
                " <Polynomials> : domain[{}].type={} is wrong.".format(i, domain[i].__class__.__name__)
            assert np.ndim(domain[i]) == 1, \
                " <Polynomials> : ndim(domain[{}])={} is wrong.".format(i, np.ndim(domain[i]))
            if np.size(domain[i]) > 1:
                assert np.all(np.diff(domain[i]) > 0) and np.max(domain[i]) <= 1 and np.min(domain[i]) >= -1, \
                    " <Polynomials> : domain[i]={} wrong, need to be increasing and bounded in [-1, 1].".format(
                            domain[i])
            else:
                pass
        if self.ndim == 3:
            if compute_xietasigma:
                _meshgrid_ravel_ = \
                    {'N': [item.ravel('F') for item in np.meshgrid([-1], domain[1], domain[2], indexing='ij')],
                     'S': [item.ravel('F') for item in np.meshgrid([+1], domain[1], domain[2], indexing='ij')],
                     'W': [item.ravel('F') for item in np.meshgrid(domain[0], [-1], domain[2], indexing='ij')],
                     'E': [item.ravel('F') for item in np.meshgrid(domain[0], [+1], domain[2], indexing='ij')],
                     'B': [item.ravel('F') for item in np.meshgrid(domain[0], domain[1], [-1], indexing='ij')],
                     'F': [item.ravel('F') for item in np.meshgrid(domain[0], domain[1], [+1], indexing='ij')]
                     }
            else:
                _meshgrid_ravel_ = None
            _basis_ = {}
            if k == 0:
                bf_xi = self.basises[0].node_basis(x=domain[0])
                bf_et = self.basises[1].node_basis(x=domain[1])
                bf_si = self.basises[2].node_basis(x=domain[2])
                _basis_['N'] = (np.kron(bf_si, bf_et),)
                _basis_['W'] = (np.kron(bf_si, bf_xi),)
                _basis_['B'] = (np.kron(bf_et, bf_xi),)
            elif k == 1:
                bfn_xi = self.basises[0].node_basis(x=domain[0])
                bfn_et = self.basises[1].node_basis(x=domain[1])
                bfn_si = self.basises[2].node_basis(x=domain[2])
                bfe_xi = self.basises[0].edge_basis(x=domain[0])
                bfe_et = self.basises[1].edge_basis(x=domain[1])
                bfe_si = self.basises[2].edge_basis(x=domain[2])
                _basis_['N'] = (np.kron(bfn_si, bfe_et), np.kron(bfe_si, bfn_et))
                _basis_['W'] = (np.kron(bfn_si, bfe_xi), np.kron(bfe_si, bfn_xi))
                _basis_['B'] = (np.kron(bfn_et, bfe_xi), np.kron(bfe_et, bfn_xi))
            elif k == 2:
                bf_xi = self.basises[0].edge_basis(x=domain[0])
                bf_et = self.basises[1].edge_basis(x=domain[1])
                bf_si = self.basises[2].edge_basis(x=domain[2])
                _basis_['N'] = (np.kron(bf_si, bf_et),)
                _basis_['W'] = (np.kron(bf_si, bf_xi),)
                _basis_['B'] = (np.kron(bf_et, bf_xi),)
            else:
                raise DimensionError()
        else:
            raise DimensionError(
                " <Polynomials> : evaluate_basis_functions for ndim={} not coded".format(self.ndim))
        _basis_['S'] = _basis_['N']
        _basis_['E'] = _basis_['W']
        _basis_['F'] = _basis_['B']
        return _meshgrid_ravel_, _basis_