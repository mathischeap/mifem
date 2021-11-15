

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
        assert len(domain) == self.ndim, " <Polynomials> : domain shape={} wrong.".format(np.shape(domain))
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

        assert self.ndim == 2

        if compute_xietasigma:
            _meshgrid_ravel_ = \
                {'U': [item.ravel('F') for item in np.meshgrid([-1], domain[1], indexing='ij')],
                 'D': [item.ravel('F') for item in np.meshgrid([+1], domain[1], indexing='ij')],
                 'L': [item.ravel('F') for item in np.meshgrid(domain[0], [-1], indexing='ij')],
                 'R': [item.ravel('F') for item in np.meshgrid(domain[0], [+1], indexing='ij')]}
        else:
            _meshgrid_ravel_ = None

        _basis_ = dict()
        if k == 0:
            bf_xi = self.basises[0].node_basis(x=domain[0])
            bf_et = self.basises[1].node_basis(x=domain[1])
            _basis_['U'] = (bf_et,)
            _basis_['L'] = (bf_xi,)
        elif k == 1:
            bf_xi = self.basises[0].edge_basis(x=domain[0])
            bf_et = self.basises[1].edge_basis(x=domain[1])
            _basis_['U'] = (bf_et,)
            _basis_['L'] = (bf_xi,)
        else:
            raise DimensionError()

        _basis_['D'] = _basis_['U']
        _basis_['R'] = _basis_['L']
        return _meshgrid_ravel_, _basis_