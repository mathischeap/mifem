
"""

"""


import numpy as np

# noinspection PyUnresolvedReferences
class EvaluatingEdgeBasis:
    """"""
    def DO_evaluate_edge_basis_at_meshgrid(self, k, *domain):
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
            if k == 0:
                bf_NS = self.basises[0].node_basis(x=domain[0])
                bf_WE = self.basises[1].node_basis(x=domain[1])
                bf_BF = self.basises[2].node_basis(x=domain[2])
            elif k == 1:
                bf_NS = self.basises[0].edge_basis(x=domain[0])
                bf_WE = self.basises[1].edge_basis(x=domain[1])
                bf_BF = self.basises[2].edge_basis(x=domain[2])
            else:
                raise Exception()
        else:
            raise Exception()

        bNS = (bf_NS,)
        bWE = (bf_WE,)
        bBF = (bf_BF,)

        _basis_ = {'WB': bNS,
                   'EB': bNS,
                   'WF': bNS,
                   'EF': bNS,
                   'NB': bWE,
                   'SB': bWE,
                   'NF': bWE,
                   'SF': bWE,
                   'NW': bBF,
                   'SW': bBF,
                   'NE': bBF,
                   'SE': bBF}

        return _basis_
